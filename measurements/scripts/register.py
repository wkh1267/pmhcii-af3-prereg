import json, gzip, collections, numpy as np, sys
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.join(_HERE, "..")


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║ SUP_RMSD_MAX 是 §02A 的 S1 參數，**已停用**。                            ║
# ║ 論文的現行引用值一律取自 registers_nothresh.json（SUP_RMSD_MAX=inf）。   ║
# ║ 預設仍為 3.0，只為了讓 data/registers.json 與 expected/ 的黃金檔逐位元   ║
# ║ 重現 —— 可重現性、正確性、納入範圍是三個獨立的主張（§01）。              ║
# ╚══════════════════════════════════════════════════════════════════════════╝
#
# [2026-10-07] 疊合門檻參數化。定義八已撤銷該門檻：3WEX 以 3.47 Å 被排除，
# 而其 P1 指認正確。設 SUP_RMSD_MAX=inf 取得全集。
#
# [2026-10-09] 偏誤已量化，並查明機制。**它量的不是結構品質，是序列離 1DLH
# （HLA-DR1）參考框架多遠。**
#
#     locus    實例   被砍   淘汰率    疊合錨點數中位
#     DR        229     0     0.0%        168
#     DQ         92     5     5.4%        104
#     DP         31     4    12.9%        110      ← 砍掉的四筆全是 3WEX 的複本
#     H2-IA      60    18    30.0%        106
#     H2-IE      33    12    36.4%        132
#     BL          7     6    85.7%        105
#
# 機制：下面的 nw() 只回傳**殘基相同**的對位（`if a[i-1] == b[j-1]`），錯配被
# 靜默丟棄。序列越偏離 DR1 → 錨點越少 → Kabsch 越不穩 → rms 越高 → 被此門檻砍掉。
# 因此門檻對非 DR 的 locus 有系統性偏誤，而本論文的焦點正是非 DR。
#
# ⚠ 砍掉 3WEX 的後果特別嚴重：那是焦點 allele（DPA1*02:02/DPB1*05:01）在全 PDB
# 的唯一結構，套用門檻後焦點 allele 的 ground-truth register 數 = 0。
#
# 另：撤銷門檻後，有兩個數字往**削弱本提案自身主張**的方向動（變異比 1.69→1.61×、
# |orient| 97.1→95.8%）。那是採全集的理由，不是保留門檻的理由。
#
# 完整分析：reports/DP_STRUCTURE_AUDIT_20261007.md §5
# 替代判準（p1_dist）待 M1 以外部驗證集決定，見 §02A S1 與定義八。
SUP_RMSD_MAX = float(os.environ.get("SUP_RMSD_MAX", "3.0"))

W  = os.path.join(_ROOT, "data") + os.sep
D  = os.path.join(_ROOT, "cif") + os.sep
rows = {r["pdb"]: r for r in json.load(open(W + "rows7.json"))}

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
       'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
       'THR':'T','TRP':'W','TYR':'Y','VAL':'V','MSE':'M','SEC':'U','PYL':'O'}

def parse_ca(pid):
    """(entity_id, label_seq_id, asym_id, comp_id, x,y,z) for CA of model 1."""
    out = []
    with gzip.open(D + pid + ".cif.gz", "rt", errors="ignore") as fh:
        cols, hdr = [], False
        for line in fh:
            if line.startswith("_atom_site."):
                cols.append(line.strip().split(".")[1]); hdr = True; continue
            if hdr and (line.startswith("ATOM") or line.startswith("HETATM")):
                p = line.split()
                if len(p) < len(cols): continue
                d = dict(zip(cols, p))
                if d.get("label_atom_id") != "CA": continue
                if d.get("pdbx_PDB_model_num", "1") not in ("1", "."): continue
                alt = d.get("label_alt_id", ".")
                if alt not in (".", "A"): continue
                try: seq = int(d["label_seq_id"])
                except: continue
                out.append((d["label_entity_id"], seq, d["label_asym_id"], d["label_comp_id"],
                            float(d["Cartn_x"]), float(d["Cartn_y"]), float(d["Cartn_z"])))
            elif hdr and line.startswith("#") and out:
                break
    return out

def nw(a, b, match=2, mis=-1, gap=-3):
    """global alignment -> list of (i,j) 0-based matched index pairs (no gaps)."""
    n, m = len(a), len(b)
    S = np.zeros((n+1, m+1), dtype=np.int32)
    S[:, 0] = np.arange(n+1) * gap
    S[0, :] = np.arange(m+1) * gap
    P = np.zeros((n+1, m+1), dtype=np.int8)   # 0 diag 1 up 2 left
    P[:, 0] = 1; P[0, :] = 2
    av = np.frombuffer(a.encode(), dtype=np.uint8)
    bv = np.frombuffer(b.encode(), dtype=np.uint8)
    for i in range(1, n+1):
        sc = np.where(av[i-1] == bv, match, mis)
        row_prev = S[i-1]
        diag = row_prev[:-1] + sc
        up   = row_prev[1:] + gap
        cur  = S[i]
        for j in range(1, m+1):
            best = diag[j-1]; pt = 0
            if up[j-1] > best: best = up[j-1]; pt = 1
            left = cur[j-1] + gap
            if left > best: best = left; pt = 2
            cur[j] = best; P[i, j] = pt
    pairs = []
    i, j = n, m
    while i > 0 and j > 0:
        pt = P[i, j]
        if pt == 0:
            if a[i-1] == b[j-1]: pairs.append((i-1, j-1))
            i -= 1; j -= 1
        elif pt == 1: i -= 1
        else: j -= 1
    return pairs[::-1]

def kabsch(P, Q):
    """rotation+translation mapping P onto Q."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    Dm = np.diag([1, 1, d])
    R = Vt.T @ Dm @ U.T
    return R, qc - R @ pc

def chain_map(ca):
    """(entity,asym) -> {label_seq_id: (comp, xyz)}"""
    g = collections.defaultdict(dict)
    for e, s, a, c, x, y, z in ca:
        g[(e, a)][s] = (c, np.array([x, y, z]))
    return g

def centroid(d): return np.mean([v[1] for v in d.values()], axis=0)

# ---------------- reference: 1DLH, HA306-318 PKYVKQNTLKLAT, P1 = Y (index 2) ------------
REF = "1DLH"; REF_P1_IDX = 2
ref_row = rows[REF]
ref_ca  = parse_ca(REF)
ref_g   = chain_map(ref_ca)
ref_a_seq = ref_row["alpha"][0][4]
ref_b_seq = ref_row["beta"][0][4]
ref_a_eid = ref_row["alpha"][0][0].split("_")[1]
ref_b_eid = ref_row["beta"][0][0].split("_")[1]
ref_p_eid = ref_row["sep"][0][0].split("_")[1]
ref_a = next(d for (e, a), d in ref_g.items() if e == ref_a_eid)
ref_b = next(d for (e, a), d in ref_g.items() if e == ref_b_eid)
ref_p = next(d for (e, a), d in ref_g.items() if e == ref_p_eid)
ref_p_ids = sorted(ref_p)
REF_P1 = ref_p[ref_p_ids[REF_P1_IDX]][1]
REF_P9 = ref_p[ref_p_ids[REF_P1_IDX + 8]][1]
REF_DIR = REF_P9 - REF_P1; REF_DIR /= np.linalg.norm(REF_DIR)
A1_MAX, B1_MAX = 84, 95          # groove-forming domains only

def seq_from(d):
    return "".join(AA3.get(d[k][0], "X") for k in sorted(d))

def assign(pid):
    r = rows[pid]
    if not r["alpha"] or not r["beta"]: return []
    ca = parse_ca(pid)
    if not ca: return []
    g = chain_map(ca)
    a_eid = r["alpha"][0][0].split("_")[1]; a_seq = r["alpha"][0][4]
    b_eid = r["beta"][0][0].split("_")[1];  b_seq = r["beta"][0][4]
    pa = {i: j for i, j in nw(a_seq, ref_a_seq) if j < A1_MAX}
    pb = {i: j for i, j in nw(b_seq, ref_b_seq) if j < B1_MAX}
    if len(pa) < 40 or len(pb) < 40: return []
    peps = []
    for eid, L, seq, desc in r["sep"] + r["ii"]:
        peps.append((eid.split("_")[1], None, seq, desc, "sep"))
    for f in r["fused"]:
        if f["kind"] == "fused_peptide":
            peps.append((f["eid"].split("_")[1], f.get("plen"), f.get("pep", ""), f.get("desc", ""), "fused"))
    a_inst = [(k, v) for k, v in g.items() if k[0] == a_eid]
    b_inst = [(k, v) for k, v in g.items() if k[0] == b_eid]
    if not a_inst or not b_inst: return []
    out = []
    for pe, plen, pseq, desc, kind in peps:
        for (e, asym), pd in g.items():
            if e != pe: continue
            ids = sorted(pd)
            if plen: ids = [i for i in ids if i <= plen]
            if len(ids) < 9: continue
            pc = np.mean([pd[i][1] for i in ids], axis=0)
            ak, ad = min(a_inst, key=lambda kv: np.linalg.norm(centroid(kv[1]) - pc))
            ac = centroid(ad)
            bk, bd = min(b_inst, key=lambda kv: np.linalg.norm(centroid(kv[1]) - ac))
            if np.linalg.norm(ac - pc) > 40: continue
            P, Q = [], []
            for i, j in pa.items():
                if (i + 1) in ad: P.append(ad[i + 1][1]); Q.append(ref_a[sorted(ref_a)[j]][1])
            for i, j in pb.items():
                if (i + 1) in bd: P.append(bd[i + 1][1]); Q.append(ref_b[sorted(ref_b)[j]][1])
            if len(P) < 60: continue
            P = np.array(P); Q = np.array(Q)
            R, t = kabsch(P, Q)
            rms = float(np.sqrt(((P @ R.T + t - Q) ** 2).sum(1).mean()))
            if rms > SUP_RMSD_MAX: continue
            xyz = {i: pd[i][1] @ R.T + t for i in ids}
            dists = {i: float(np.linalg.norm(xyz[i] - REF_P1)) for i in ids}
            p1_id = min(dists, key=dists.get)
            k = ids.index(p1_id)
            v = xyz[ids[-1]] - xyz[ids[0]]
            nv = np.linalg.norm(v)
            if nv < 1e-6: continue
            orient = float(np.dot(v / nv, REF_DIR))
            rev = orient < 0
            start = (k - 8) if rev else k
            seq_obs = "".join(AA3.get(pd[i][0], "X") for i in ids)
            ok = 0 <= start <= len(ids) - 9
            out.append(dict(pdb=pid, asym=asym, kind=kind, L=len(ids), pep=seq_obs,
                            p1_index=k, start=start if ok else None,
                            core=(seq_obs[start:start+9] if ok else None),
                            reverse=bool(rev), orient=round(orient, 3),
                            p1_dist=round(dists[p1_id], 2), sup_rmsd=round(rms, 2),
                            loci=r.get("loci"), rel=r.get("rel"), res=r.get("res"),
                            has_tcr=bool(r.get("tcr") or r.get("tcrfuse"))))
    return out

res, bad = [], []
for pid in sorted(rows):
    try:
        res += assign(pid)
    except Exception as ex:
        bad.append((pid, repr(ex)[:80]))
_OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_ROOT, "data", "registers.json")
json.dump(res, open(_OUT, "w"), indent=0)
print("assigned instances:", len(res), " structures:", len({r['pdb'] for r in res}), " errors:", len(bad))
for b in bad[:10]: print("  ERR", b)
