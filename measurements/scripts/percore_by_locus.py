#!/usr/bin/env python3
"""逐位置側鏈接觸佔比，按 locus 分層。

contacts_fixed.py 的 stdout 只報全體平均。全體平均由 DR 主導（185/360 條），
因此看不出 DP 的分布與 DR 是兩個形狀。本腳本把同一份 contacts_fixed.json
按 locus 切開重算，聚合方式與 §05 相同（逐結構比例再平均，mean-of-ratios）。

這支腳本回答的問題是「本提案取 P1/P4/P6/P9 四個口袋，在 DP 上是否仍然合理」。
結論見 README.md「逐 locus 分層」一節：接觸數不等於專一性，四口袋的依據
應取自 DP 的 motif 與結構文獻，不是接觸數排序。

用法:  python3 percore_by_locus.py
輸入:  ../data/contacts_fixed.json
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from locus_fix import locus_of, note as locus_note
import json, os, statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

# [2026-10-07] 輸入可覆寫，預設與原行為相同。
# ⚠ 預設改讀 **nothresh**（2026-10-09）。contacts_fixed.json 套了 §02A 的 S1 疊合門檻，
# 而該參數已停用 —— 它在 DP 上唯一砍掉的是焦點 allele 的唯一結構 3WEX。
# 論文的現行引用值一律取自全集。設 CONTACTS=../data/contacts_fixed.json 可取回舊行為。
_IN = os.environ.get("CONTACTS") or os.path.join(DATA, "contacts_nothresh.json")
rows = json.load(open(_IN, encoding="utf-8"))

ANCHOR4 = (0, 3, 5, 8)   # P1 P4 P6 P9
P7 = 6
LABELS = [f"P{i + 1}" for i in range(9)]


def share(rs, i):
    """每條鏈的 percore[i]/sc 再平均。sc=0 的鏈不計（無側鏈接觸）。"""
    v = [r["percore"][i] / r["sc"] for r in rs if r["sc"] > 0]
    return 100 * st.mean(v) if v else 0.0


def block(label, rs):
    if not rs:
        print(f"{label}: n=0")
        return
    per = [share(rs, i) for i in range(9)]
    a4 = sum(per[i] for i in ANCHOR4)
    order = sorted(range(9), key=lambda i: (-per[i], i))
    print(f"{label}  (n={len(rs)} chains, {len({r['pdb'] for r in rs})} entries)")
    print("  " + "  ".join(f"{LABELS[i]} {per[i]:5.1f}%" for i in range(9)))
    print(f"  P1+P4+P6+P9: {a4:5.1f}%   P7: {per[P7]:5.1f}% (rank {order.index(P7) + 1})"
          f"   with P7: {a4 + per[P7]:5.1f}%")


print("side-chain contact share by core position, stratified by locus")
print("aggregation: mean-of-ratios over peptide chains (same as section 05)")
print(locus_note())
print()
# ⚠ locus 一律經 locus_of()，否則 3 筆蝙蝠會算進 DR（見 scripts/locus_fix.py）
# ALL **包含**蝙蝠 —— §01 本來就是跨物種量測（已含小鼠、雞），排掉蝙蝠反而不一致。
# 要排掉的只有 DR 那一列：蝙蝠不是人類 DR。
block("ALL", rows)
print()
for loc in ("DR", "DQ", "DP"):
    block(loc, [r for r in rows if locus_of(r) == loc])
    print()
bat = [r for r in rows if locus_of(r) == "BAT"]
if bat:
    block("BAT（蝙蝠，先前誤計入 DR）", bat)
    print()

dp = [r for r in rows if locus_of(r) == "DP" and r["sc"] > 0]
print("DP robustness")
print("  per-position spread (share of DP side-chain contacts)")
for i in range(9):
    v = sorted(100 * r["percore"][i] / r["sc"] for r in dp)
    tag = " <- in four-pocket set" if i in ANCHOR4 else (" <- P7" if i == P7 else "")
    print(f"    {LABELS[i]}  mean {st.mean(v):5.1f}%  median {st.median(v):5.1f}%"
          f"  min {v[0]:5.1f}%  max {v[-1]:5.1f}%{tag}")

f = lambda r, i: r["percore"][i] / r["sc"]
print(f"  chain-level: P7 > P4 in {sum(1 for r in dp if f(r, 6) > f(r, 3))}/{len(dp)},"
      f"  P7 > P9 in {sum(1 for r in dp if f(r, 6) > f(r, 8))}/{len(dp)},"
      f"  P7 > P1 in {sum(1 for r in dp if f(r, 6) > f(r, 0))}/{len(dp)}")

by_entry = {}
for r in dp:
    by_entry.setdefault(r["pdb"], []).append(r)
per = [100 * st.mean(st.mean(f(r, i) for r in rs) for rs in by_entry.values())
       for i in range(9)]
order = sorted(range(9), key=lambda i: (-per[i], i))
print(f"  entry-level (n={len(by_entry)}, each entry averaged first):")
print("    " + "  ".join(f"{LABELS[i]} {per[i]:5.1f}%" for i in range(9)))
print(f"    P1+P4+P6+P9: {sum(per[i] for i in ANCHOR4):5.1f}%"
      f"   P7: {per[P7]:5.1f}% (rank {order.index(P7) + 1})")
