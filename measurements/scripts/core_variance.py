"""How much of the pMHC-II interface is register-INDEPENDENT?

Class II anchors peptides mainly through conserved H-bonds to the peptide MAIN
CHAIN. Those contacts are made whatever residue sits there, so they are identical
under a register shift. Only SIDE-CHAIN contacts (chiefly P1/P4/P6/P9) can tell
registers apart. Quantify the split on real structures.
"""
import json, gzip, collections, numpy as np
import os, sys
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.join(_HERE, "..")


W = os.path.join(_ROOT, "data") + os.sep
D = os.path.join(_ROOT, "cif") + os.sep
rows = {r["pdb"]: r for r in json.load(open(W + "rows7.json"))}
# [改 2026-10-11] 預設改讀 registers_nothresh.json：初稿 3.2.1 節明文不採用疊合
# 品質門檻，而本檔原本讀套了門檻的檔，跑出來是 1.69 倍 / n = 360，與初稿的
# 1.61 倍 / n = 404 對不上 —— 等於把已撤銷的門檻又裝回去。--thresh 可切回舊基準。
_SRC = "registers.json" if "--thresh" in sys.argv else "registers_nothresh.json"
REG  = json.load(open(os.path.join(_ROOT, "data", _SRC)))

# [FIX 2026-10-03] 共價融合構築的胜肽長度。rows7.json 的 fused_peptide 條目帶有 plen，
# register.py:138 已用它截斷殘基清單；contacts.py 原本沒有，導致整條融合鏈（胜肽 + linker
# + 完整 beta 鏈）被當成胜肽，使 beta-alpha 的蛋白內部接觸被計入界面統計。
# 鍵為 (pdb, entity_id)，與 register.py:129 的 f["eid"].split("_")[1] 同一個來源。
PLEN = {}
for _r in rows.values():
    for _f in _r.get("fused", []):
        if _f.get("kind") == "fused_peptide" and _f.get("plen"):
            PLEN[(_r["pdb"], _f["eid"].split("_")[1])] = int(_f["plen"])

BB = {"N", "CA", "C", "O", "OXT"}

def parse_atoms(pid):
    """asym -> (seq_ids, atom_names, elements, coords) for heavy atoms, model 1."""
    g = collections.defaultdict(lambda: ([], [], [], []))
    with gzip.open(D + pid + ".cif.gz", "rt", errors="ignore") as fh:
        cols, hdr = [], False
        for line in fh:
            if line.startswith("_atom_site."):
                cols.append(line.strip().split(".")[1]); hdr = True; continue
            if hdr and (line.startswith("ATOM") or line.startswith("HETATM")):
                p = line.split()
                if len(p) < len(cols): continue
                d = dict(zip(cols, p))
                if d.get("type_symbol") in ("H", "D"): continue
                if d.get("pdbx_PDB_model_num", "1") not in ("1", "."): continue
                if d.get("label_alt_id", ".") not in (".", "A"): continue
                try: s = int(d["label_seq_id"])
                except: continue
                a = g[d["label_asym_id"]]
                a[0].append(s); a[1].append(d["label_atom_id"]); a[2].append(d["label_entity_id"])
                a[3].append((float(d["Cartn_x"]), float(d["Cartn_y"]), float(d["Cartn_z"])))
            elif hdr and line.startswith("#") and g:
                break
    return {k: (np.array(v[0]), np.array(v[1]), np.array(v[2]), np.array(v[3])) for k, v in g.items()}

CUT = 4.0
out = []
for r in REG:
    if r["start"] is None or "X" in r["pep"]: continue
    pid = r["pdb"]; row = rows[pid]
    try: A = parse_atoms(pid)
    except Exception: continue
    if r["asym"] not in A: continue
    a_eid = row["alpha"][0][0].split("_")[1]; b_eid = row["beta"][0][0].split("_")[1]
    ps, pn, pe, pc = A[r["asym"]]
    if len(pc) == 0: continue
    # [FIX 2026-10-03] 融合構築只取前 plen 個殘基作為胜肽。
    # plen 必須與 register.py:138 用的值完全相同，否則下方 ids[start:start+9] 的
    # P1-P9 指派會整體位移。
    pl = PLEN.get((pid, pe[0])) if r.get("kind") == "fused" else None
    if pl is not None:
        keep = ps <= pl
        ps, pn, pe, pc = ps[keep], pn[keep], pe[keep], pc[keep]
        if len(pc) == 0: continue
    pcen = pc.mean(0)
    # MHC atoms = nearest alpha instance + nearest beta instance
    mhc = []
    for eid in (a_eid, b_eid):
        best, bd = None, 1e9
        for asym, (s, n, e, c) in A.items():
            if e[0] != eid or len(c) == 0: continue
            # [FIX 2026-10-03] 融合構築中，beta 的 entity 就是融合鏈本身。
            # 取 MHC 原子時要反向排除胜肽段（只留 seq_id > plen），否則胜肽會與自己比距離。
            if pl is not None and asym == r["asym"]:
                mk = s > pl
                if not mk.any(): continue
                c = c[mk]
            d = np.linalg.norm(c.mean(0) - pcen)
            if d < bd: bd, best = d, c
        if best is None: break
        mhc.append(best)
    if len(mhc) != 2: continue
    M = np.vstack(mhc)
    d2 = ((pc[:, None, :] - M[None, :, :]) ** 2).sum(-1)
    near = (d2 <= CUT * CUT)
    ncon = near.sum(1)                      # contacts per peptide atom
    isbb = np.array([n in BB for n in pn])
    tot = int(ncon.sum())
    if tot == 0: continue
    bbc = int(ncon[isbb].sum()); scc = tot - bbc
    # per-core-position side-chain contacts
    ids = sorted(set(ps.tolist()))
    idx = {s: k for k, s in enumerate(ids)}
    start = r["start"]
    core_ids = ids[start:start + 9]
    percore, percore_bb = [], []
    for k, sid in enumerate(core_ids):
        msc = (ps == sid) & (~isbb)
        mbb = (ps == sid) & isbb
        percore.append(int(ncon[msc].sum()))
        percore_bb.append(int(ncon[mbb].sum()))
    if r["reverse"]:
        percore = percore[::-1]; percore_bb = percore_bb[::-1]   # P1..P9 order
    out.append(dict(pdb=pid, asym=r["asym"], L=r["L"], start=start, reverse=r["reverse"],
                    loci=r["loci"], tot=tot, bb=bbc, sc=scc,
                    percore=percore, percore_bb=percore_bb))

# ⚠ 輸出檔名必須跟著輸入基準走。2026-10-12 改預設為 nothresh 時沒改這一行，
#   結果把有門檻的黃金檔 core_variance.json（360 筆）覆蓋成 404 筆。
_OUT = "core_variance_nothresh.json" if "nothresh" in _SRC else "core_variance.json"
json.dump(out, open(os.path.join(_ROOT, "data", _OUT), "w"), indent=0)

# ─────────────────────────────────────────────────────────────────────
# core 內的主鏈 vs 側鏈變異分析（只用 P1–P9 九個殘基，排除 flanking 干擾）
# ─────────────────────────────────────────────────────────────────────
import statistics as _st  # 標準差一律取樣本標準差（ddof=1），與初稿同慣例
def _cv(v): return _st.stdev(v) / _st.mean(v) if _st.mean(v) else float("nan")

core_bb  = [sum(o["percore_bb"]) for o in out]          # core 九個殘基的主鏈接觸
core_sc  = [sum(o["percore"])    for o in out]          # core 九個殘基的側鏈接觸
core_tot = [a + b for a, b in zip(core_bb, core_sc)]
anc_sc   = [sum(o["percore"][i] for i in (0, 3, 5, 8)) for o in out]

print(f"instances analysed: {len(out)}")
print("CORE-ONLY analysis (P1-P9, 9 residues, no flanking)")
print()
print(f"{'quantity':34s}{'mean':>9s}{'sd':>9s}{'CV':>9s}")
print("-" * 61)
for lab, v in [("core backbone contacts", core_bb),
               ("core side-chain contacts", core_sc),
               ("core anchor side-chain (P1/4/6/9)", anc_sc),
               ("core total contacts", core_tot)]:
    print(f"{lab:34s}{_st.mean(v):9.1f}{_st.stdev(v):9.1f}{_cv(v):9.3f}")
print()
print(f"variance ratio  anchor / backbone = {_cv(anc_sc)/_cv(core_bb):.2f}x")
print()
fb = [b / t for b, t in zip(core_bb, core_tot)]
print(f"backbone share of CORE contacts: mean {_st.mean(fb)*100:.1f}%  "
      f"median {_st.median(fb)*100:.1f}%  CV {_cv(fb):.3f}")
print()
print("per-position CV (side-chain):")
for i in range(9):
    v = [o["percore"][i] for o in out]
    tag = " *anchor" if i in (0, 3, 5, 8) else ""
    print(f"  P{i+1}: mean {_st.mean(v):5.1f}  CV {_cv(v):.3f}{tag}")
print()
print("per-position CV (backbone):")
for i in range(9):
    v = [o["percore_bb"][i] for o in out]
    print(f"  P{i+1}: mean {_st.mean(v):5.1f}  CV {_cv(v):.3f}")
