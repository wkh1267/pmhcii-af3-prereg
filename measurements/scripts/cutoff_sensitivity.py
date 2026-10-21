#!/usr/bin/env python3
"""界面接觸截斷距離的敏感度。[新增 2026-10-13]

  python3 scripts/cutoff_sensitivity.py            重算四個截斷並寫出 data 下的摘要
  python3 scripts/cutoff_sensitivity.py --check    只比對既有檔與提案所載的值

為什麼要有這一支：4 Å 決定了表 2、圖 2、圖 3 與「HLA-DP 的 P6 為 HLA-DR 的 2.6 倍」
的每一個數字，而該比值隨截斷單調上升（3.5 Å 1.2 倍、4.0 Å 2.6 倍、5.0 Å 3.1 倍）。
跳幅集中在 3.5 → 4.0：3.5 Å 大致為氫鍵距離，所收者以極性接觸為主而非側鏈堆積，
與「口袋容不容得下大側鏈」所問者不同，故 4.0 至 5.0 Å 才是本宣稱的敏感度範圍。

⚠ 每個截斷都要重讀 291 個 mmCIF，四個截斷約數分鐘。
"""
import collections, json, os, statistics as st, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
OUT  = os.path.join(ROOT, "data", "contact_cutoff_sensitivity.json")
CUTS = ("3.5", "4.0", "4.5", "5.0")
sys.path.insert(0, HERE)
from locus_fix import locus_of

# 提案 4.2 節所載者，作為 assert 錨點
EXPECT = {"bb_4.0": 44.4, "ratio_4.0": 2.6, "ratio_3.5": 1.2, "ratio_5.0": 3.1}


def summarise(path, rows):
    cn = json.load(open(path, encoding="utf-8"))
    bb = sum(r["bb"] / r["tot"] for r in cn if r["tot"]) / sum(1 for r in cn if r["tot"])
    by = collections.defaultdict(list)
    for r in cn:
        if r.get("sc", 0) > 0:
            by[locus_of(r, rows)].append(r)
    pos = lambda rs, i: 100 * st.mean(x["percore"][i] / x["sc"] for x in rs)
    return {"backbone": round(bb * 100, 1),
            "P6_DP": round(pos(by["DP"], 5), 1), "P6_DR": round(pos(by["DR"], 5), 1),
            "ratio_P6": round(pos(by["DP"], 5) / pos(by["DR"], 5), 2)}


def main():
    rows = {r["pdb"]: r for r in json.load(open(os.path.join(ROOT, "data", "rows7.json"),
                                                encoding="utf-8"))}
    if "--check" in sys.argv:
        rec = json.load(open(OUT, encoding="utf-8"))["cutoffs"]
    else:
        rec = {}
        with tempfile.TemporaryDirectory() as td:
            for c in CUTS:
                dst = os.path.join(td, f"c{c}.json")
                env = dict(os.environ,
                           REGISTERS=os.path.join(ROOT, "data", "registers_nothresh.json"),
                           CONTACTS_OUT=dst, CUTOFF=c)
                subprocess.run([sys.executable, "-I", os.path.join(HERE, "contacts_fixed.py")],
                               env=env, check=True, capture_output=True)
                rec[c] = summarise(dst, rows)
                print(f"  ↻ {c} Å")
        json.dump({"unit": "百分比為逐胜肽鏈的比例再平均；ratio_P6 為 P6_DP ÷ P6_DR",
                   "note": "3.5 Å 大致為氫鍵距離，所收者以極性接觸為主，不列入本宣稱的敏感度範圍",
                   "cutoffs": rec},
                  open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
        print(f"  寫出 {os.path.relpath(OUT)}")
    bad = []
    for k, want in EXPECT.items():
        q, c = k.split("_")
        got = rec[c]["backbone" if q == "bb" else "ratio_P6"]
        ok = abs(got - want) < 0.06
        bad += [] if ok else [k]
        print(f"  {'✓' if ok else '✗'} {c} Å {q:5s} {got}"
              f"{'' if ok else f'　提案所載 {want}'}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
