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
REG  = json.load(open(os.path.join(_ROOT, "data", "registers.json")))

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
    pcen = pc.mean(0)
    # MHC atoms = nearest alpha instance + nearest beta instance
    mhc = []
    for eid in (a_eid, b_eid):
        best, bd = None, 1e9
        for asym, (s, n, e, c) in A.items():
            if e[0] != eid or len(c) == 0: continue
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

json.dump(out, open(os.path.join(_ROOT, "data", "contacts.json"), "w"), indent=0)
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
