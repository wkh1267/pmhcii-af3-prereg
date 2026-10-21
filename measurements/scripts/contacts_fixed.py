"""How much of the pMHC-II interface is register-INDEPENDENT?

Class II anchors peptides mainly through conserved H-bonds to the peptide MAIN
CHAIN. Those contacts are made whatever residue sits there, so they are identical
under a register shift. Only SIDE-CHAIN contacts (chiefly P1/P4/P6/P9) can tell
registers apart. Quantify the split on real structures.
"""
import json, gzip, collections, numpy as np
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.join(_HERE, "..")


W = os.path.join(_ROOT, "data") + os.sep
D = os.path.join(_ROOT, "cif") + os.sep
rows = {r["pdb"]: r for r in json.load(open(W + "rows7.json"))}
# [2026-10-07] 輸入／輸出可覆寫，預設與原行為相同。
_REG_IN = os.environ.get("REGISTERS") or os.path.join(_ROOT, "data", "registers.json")
REG  = json.load(open(_REG_IN))

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

# 截斷距離。預設 4.0 Å 為重原子接觸的通用值；以 CUTOFF 覆寫可做敏感度，
# 見 scripts/cutoff_sensitivity.py。不設時行為與原版完全相同。
CUT = float(os.environ.get("CUTOFF", "4.0"))
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
    percore = []
    for k, sid in enumerate(core_ids):
        m = (ps == sid) & (~isbb)
        percore.append(int(ncon[m].sum()))
    if r["reverse"]: percore = percore[::-1]     # express in P1..P9 order
    out.append(dict(pdb=pid, asym=r["asym"], L=r["L"], start=start, reverse=r["reverse"],
                    loci=r["loci"], tot=tot, bb=bbc, sc=scc, percore=percore))

# ⚠ 只設 REGISTERS 而忘了 CONTACTS_OUT，會把有門檻的黃金檔覆蓋成 nothresh 的結果。
_DEF = "contacts_nothresh.json" if "nothresh" in _REG_IN else "contacts_fixed.json"
_OUT = os.environ.get("CONTACTS_OUT") or os.path.join(_ROOT, "data", _DEF)
json.dump(out, open(_OUT, "w"), indent=0)
n = len(out)
fb = np.array([o["bb"] / o["tot"] for o in out])
print(f"instances analysed: {n}")
print(f"peptide->MHC heavy-atom contacts <= {CUT} A")
print(f"  MAIN-CHAIN mediated: mean {fb.mean()*100:.1f}%  median {np.median(fb)*100:.1f}%  IQR {np.percentile(fb,25)*100:.1f}-{np.percentile(fb,75)*100:.1f}%")
print(f"  side-chain mediated: mean {(1-fb.mean())*100:.1f}%")
P = np.array([o["percore"] for o in out], dtype=float)
tot_sc = np.array([o["sc"] for o in out], dtype=float)
frac = P / np.maximum(tot_sc[:, None], 1)
print("\nside-chain contacts by core position (share of ALL peptide side-chain contacts):")
for i in range(9):
    print(f"  P{i+1}: {frac[:,i].mean()*100:5.1f}%   (mean {P[:,i].mean():5.1f} contacts)")
print(f"\n  P1+P4+P6+P9 together: {frac[:,[0,3,5,8]].sum(1).mean()*100:.1f}% of side-chain contacts,"
      f" {(P[:,[0,3,5,8]].sum(1)/np.array([o['tot'] for o in out])).mean()*100:.1f}% of ALL contacts")
