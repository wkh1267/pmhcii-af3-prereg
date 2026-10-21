#!/usr/bin/env python3
"""由 IPD-IMGT/HLA 取得各 class II 基因的相異蛋白質數。[新增 2026-10-13]

  python3 scripts/fetch_imgt.py            重新取得並寫出 data/imgt_hla_<release>.json
  python3 scripts/fetch_imgt.py --check    只比對既有檔與提案所載的六個數，不連網

為什麼要有這一支：這六個數原本**硬編在 proposal/draft/figures.py 裡**，沒有來源檔、
沒有取得日期、沒有計數慣例，而提案的每結構 α/β 組合（HLA-DP 56,283、HLA-DQ 16,712、
HLA-DR 358）與 1.2.2 的「HLA-DR 的 157 倍」全由它們推導。2026-10-13 追回來源：

    IPD-IMGT/HLA release 3.65.0（git tag 3650）的 Allelelist.txt，
    逐基因數相異的二欄位（蛋白質層級）名稱，排除空效等位基因（後綴 N）。

六個數在此慣例下精確命中。⚠ 不排除 N 會得到 DRB1 2756；把 Q/L/S/C/A 也一併排除會得到
2600 —— 兩者都不是提案所載的值，故慣例必須照寫，不可憑印象重算。
"""
import collections, json, os, re, sys, urllib.request

HERE    = os.path.dirname(os.path.abspath(__file__))
DATA    = os.path.join(HERE, "..", "data")
RELEASE = "3.65.0"
TAG     = "3650"
URL     = f"https://raw.githubusercontent.com/ANHIG/IMGTHLA/{TAG}/Allelelist.txt"
GENES   = ("DRA", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1")
# 提案 4.5 節與圖 5 所載者，作為 assert 錨點
EXPECT  = {"DRA": 20, "DRB1": 2610, "DQA1": 535, "DQB1": 1843, "DPA1": 438, "DPB1": 1799}
OUT     = os.path.join(DATA, f"imgt_hla_{RELEASE}.json")


def tally(text):
    """逐基因數相異的二欄位名稱，排除後綴為 N（空效）者。"""
    prot = collections.defaultdict(set)
    for line in text.split("\n")[1:]:
        if "," not in line:
            continue
        gene, _, rest = line.split(",", 1)[1].partition("*")
        if gene not in GENES:
            continue
        m = re.search(r"([A-Z])$", rest)
        if m and m.group(1) == "N":          # 空效等位基因不計
            continue
        f = rest.split(":")
        if len(f) >= 2:
            prot[gene].add(f[0] + ":" + f[1])
    return {g: len(prot[g]) for g in GENES}


def main():
    if "--check" in sys.argv:
        rec = json.load(open(OUT, encoding="utf-8"))
        got = rec["proteins"]
    else:
        print(f"  取得 {URL}")
        got = tally(urllib.request.urlopen(URL, timeout=120).read().decode())
        rec = {
            "source":   "IPD-IMGT/HLA",
            "release":  RELEASE,
            "url":      URL,
            "accessed": os.environ.get("ACCESSED", "2026-10-13"),
            "unit":     "相異的二欄位（蛋白質層級）等位基因名稱數",
            "excluded": "空效等位基因（名稱後綴 N）",
            "proteins": got,
        }
        json.dump(rec, open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
        print(f"  寫出 {os.path.relpath(OUT)}")
    bad = [g for g in GENES if got[g] != EXPECT[g]]
    for g in GENES:
        print(f"  {'✓' if got[g] == EXPECT[g] else '✗'} {g:5s} {got[g]:5d}"
              f"{'' if got[g] == EXPECT[g] else f'  提案所載 {EXPECT[g]}'}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
