"""Frame-sensitivity test [2026-10-06].

Question: register.py defines the groove reference frame with 1DLH, an HLA-DR1.
Is the P1 assignment of HLA-DP structures biased by that choice?

Method: rebuild the frame from 7ZAK (HLA-DP, separate peptide chain, forward,
p1_dist 0.94 A) and re-assign every instance.  Compare P1 index and direction.
The 3.0 A superposition cutoff of register.py is deliberately NOT applied here
(revoked, see proposal_v6 definition 8).

Runtime ~90 s.  Writes frame_sensitivity.json next to data/.
"""
import sys, os, json, importlib.util
PKG = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PKG, "scripts"))

# load register.py but stop it from executing the main loop: easiest is to copy the
# functions we need by exec'ing the file up to the "# ---------------- reference" marker
src = open(os.path.join(PKG, "scripts", "register.py")).read()
head = src.split("# ---------------- reference")[0]
g = {"__file__": os.path.join(PKG, "scripts", "register.py")}
exec(compile(head, "register_head", "exec"), g)

import numpy as np
rows, parse_ca, chain_map, centroid, nw, kabsch, AA3 = (
    g["rows"], g["parse_ca"], g["chain_map"], g["centroid"], g["nw"], g["kabsch"], g["AA3"])
A1_MAX, B1_MAX = 84, 95

def build_ref(ref_pdb, p1_idx):
    r = rows[ref_pdb]
    ca = parse_ca(ref_pdb); gg = chain_map(ca)
    a_seq, b_seq = r["alpha"][0][4], r["beta"][0][4]
    a_eid = r["alpha"][0][0].split("_")[1]
    b_eid = r["beta"][0][0].split("_")[1]
    p_eid = r["sep"][0][0].split("_")[1]
    ra = next(d for (e, a), d in gg.items() if e == a_eid)
    rb = next(d for (e, a), d in gg.items() if e == b_eid)
    rp = next(d for (e, a), d in gg.items() if e == p_eid)
    ids = sorted(rp)
    P1 = rp[ids[p1_idx]][1]; P9 = rp[ids[p1_idx + 8]][1]
    D = P9 - P1; D /= np.linalg.norm(D)
    seq = "".join(AA3.get(rp[k][0], "X") for k in ids)
    return dict(a=ra, b=rb, a_seq=a_seq, b_seq=b_seq, P1=P1, DIR=D,
                core=seq[p1_idx:p1_idx + 9])

def assign(pid, R):
    r = rows[pid]
    if not r["alpha"] or not r["beta"]: return []
    ca = parse_ca(pid)
    if not ca: return []
    gg = chain_map(ca)
    a_eid = r["alpha"][0][0].split("_")[1]; a_seq = r["alpha"][0][4]
    b_eid = r["beta"][0][0].split("_")[1];  b_seq = r["beta"][0][4]
    pa = {i: j for i, j in nw(a_seq, R["a_seq"]) if j < A1_MAX}
    pb = {i: j for i, j in nw(b_seq, R["b_seq"]) if j < B1_MAX}
    if len(pa) < 40 or len(pb) < 40: return []
    peps = [(e.split("_")[1], None) for e, L, s, d in r["sep"] + r["ii"]]
    peps += [(f["eid"].split("_")[1], f.get("plen")) for f in r["fused"]
             if f["kind"] == "fused_peptide"]
    a_inst = [(k, v) for k, v in gg.items() if k[0] == a_eid]
    b_inst = [(k, v) for k, v in gg.items() if k[0] == b_eid]
    if not a_inst or not b_inst: return []
    ra_ids, rb_ids = sorted(R["a"]), sorted(R["b"])
    out = []
    for pe, plen in peps:
        for (e, asym), pd in gg.items():
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
                if (i + 1) in ad: P.append(ad[i + 1][1]); Q.append(R["a"][ra_ids[j]][1])
            for i, j in pb.items():
                if (i + 1) in bd: P.append(bd[i + 1][1]); Q.append(R["b"][rb_ids[j]][1])
            if len(P) < 60: continue
            P = np.array(P); Q = np.array(Q)
            rot, t = kabsch(P, Q)
            rms = float(np.sqrt(((P @ rot.T + t - Q) ** 2).sum(1).mean()))
            # NOTE: the 3.0 A cutoff of register.py is deliberately NOT applied here
            xyz = {i: pd[i][1] @ rot.T + t for i in ids}
            dists = {i: float(np.linalg.norm(xyz[i] - R["P1"])) for i in ids}
            p1 = min(dists, key=dists.get); k = ids.index(p1)
            v = xyz[ids[-1]] - xyz[ids[0]]; nv = np.linalg.norm(v)
            if nv < 1e-6: continue
            orient = float(np.dot(v / nv, R["DIR"]))
            start = (k - 8) if orient < 0 else k
            seq = "".join(AA3.get(pd[i][0], "X") for i in ids)
            ok = 0 <= start <= len(ids) - 9
            out.append(dict(pdb=pid, asym=asym, p1_index=k,
                            core=(seq[start:start + 9] if ok else None),
                            reverse=orient < 0, orient=round(orient, 3),
                            p1_dist=round(dists[p1], 2), sup_rmsd=round(rms, 2),
                            loci=tuple(sorted(set(r.get("loci") or [])))))
    return out

REFS = {"1DLH(DR1)": ("1DLH", 2), "7ZAK(DP4)": ("7ZAK", 4)}
res = {}
for name, (pdb, idx) in REFS.items():
    R = build_ref(pdb, idx)
    print(f"reference {name}: core={R['core']}")
    acc = []
    for pid in sorted(rows):
        try: acc += assign(pid, R)
        except Exception as ex: pass
    res[name] = {(r["pdb"], r["asym"]): r for r in acc}
    print(f"  assigned instances: {len(acc)}  structures: {len({r['pdb'] for r in acc})}")

json.dump({k: list(v.values()) for k, v in res.items()},
          open(os.path.join(PKG, "data", "frame_sensitivity.json"), "w"), indent=0)

a, b = res["1DLH(DR1)"], res["7ZAK(DP4)"]
keys = sorted(set(a) & set(b))
print(f"\ncomparable instances: {len(keys)}  (DR1frame-only {len(set(a)-set(b))}, DPframe-only {len(set(b)-set(a))})")
import collections
agree = collections.Counter(); disagree = []
for k in keys:
    loc = a[k]["loci"]
    same = (a[k]["p1_index"] == b[k]["p1_index"]) and (a[k]["reverse"] == b[k]["reverse"])
    agree[(loc, same)] += 1
    if not same: disagree.append((k, loc, a[k]["p1_index"], b[k]["p1_index"],
                                  a[k]["core"], b[k]["core"], a[k]["reverse"], b[k]["reverse"]))
loci = sorted({l for l, _ in agree})
print(f"\n{'locus':10s} {'n':>4s} {'same P1':>8s} {'differ':>7s} {'rate':>7s}")
for l in loci:
    s = agree[(l, True)]; d = agree[(l, False)]
    print(f"{str(l):10s} {s+d:4d} {s:8d} {d:7d} {d/(s+d):7.1%}")
print("\ndisagreements (pdb, asym, locus, P1_DR1, P1_DP5, core_DR1, core_DP5, rev_DR1, rev_DP5):")
for d in disagree[:40]: print("  ", d)
