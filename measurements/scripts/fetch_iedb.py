#!/usr/bin/env python3
"""由 IEDB query API 取得 MEDi 的逐等位基因陽／陰筆數。[新增 2026-10-13]

  python3 scripts/fetch_iedb.py            連網取得並寫出 data/medi_iedb_1040686.json
  python3 scripts/fetch_iedb.py --check    只比對既有檔與提案所載的六個數，不連網

為什麼要有這一支：SOURCES.md 在 2026-10-09 記為「工作目錄沒有任何 IEDB / MEDi 的
資料檔，§08 引用的六個數字目前一個都不能重跑」。2026-10-13 以本檔的查詢全部重現，
該狀態因此解除。

查詢：mhc_search 端點，以 reference_id = 1040686（Obermair 2022, Sci Adv，MEDi）過濾，
再以 mhc_allele_name 與 qualitative_measure 分組計數。筆數取自 PostgREST 的
Content-Range 標頭（需 Prefer: count=exact），不必把三萬七千列全部拉下來。

⚠ 兩個查核時發現的陷阱：
  · MEDi 只有九個等位基因，其中 DQ 只有一個，即 DQA1*03:01/DQB1*03:02。初稿原本把它
    寫成 DQA1*05:01/DQB1*02:01（DQ2.5），該等位基因在這批資料裡根本不存在；
    189 / 1,069 屬於前者。2026-10-13 已更正 3.2.2 節。
  · SOURCES.md 另記，IEDB 取用不可用 mhc_search.quantitative_measure 作為數值端點
    （該欄混入 X 光解析度）。本檔只用 qualitative_measure 的陽／陰標籤，不碰該欄。
"""
import json, os, sys, urllib.parse, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
OUT  = os.path.join(HERE, "..", "data", "medi_iedb_1040686.json")
API  = "https://query-api.iedb.org/mhc_search"
REF  = 1040686
# 提案 3.2.2 節所載者，作為 assert 錨點
EXPECT = {
    "HLA-DPA1*02:02/DPB1*05:01": (259, 1304),
    "HLA-DQA1*03:01/DQB1*03:02": (189, 1069),
    "HLA-DRB1*15:01":            (1485, 8345),
}


def count(extra):
    req = urllib.request.Request(f"{API}?select=structure_id&limit=1&reference_id=eq.{REF}&{extra}")
    req.add_header("Prefer", "count=exact")
    with urllib.request.urlopen(req, timeout=60) as fh:
        return int(fh.headers["Content-Range"].split("/")[1])


def main():
    if "--check" in sys.argv:
        rec = json.load(open(OUT, encoding="utf-8"))
        got = {k: tuple(v) for k, v in rec["alleles"].items()}
    else:
        got = {}
        for allele in EXPECT:
            a = urllib.parse.quote(allele, safe="")
            tot = count(f"mhc_allele_name=eq.{a}")
            pos = count(f"mhc_allele_name=eq.{a}&qualitative_measure=eq.Positive")
            got[allele] = (pos, tot - pos)
            print(f"  ↓ {allele}")
        json.dump({"source": "IEDB query API", "endpoint": API, "reference_id": REF,
                   "citation": "Obermair FJ et al., Sci Adv 2022 (MEDi)",
                   "accessed": os.environ.get("ACCESSED", "2026-10-13"),
                   "unit": "逐等位基因的 (陽性, 陰性) 筆數，標籤取自 qualitative_measure",
                   "alleles": {k: list(v) for k, v in got.items()}},
                  open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
        print(f"  寫出 {os.path.relpath(OUT)}")
    bad = [a for a in EXPECT if got.get(a) != EXPECT[a]]
    for a in EXPECT:
        ok = got.get(a) == EXPECT[a]
        print(f"  {'✓' if ok else '✗'} {a:28s} 陽 {got.get(a, ('?', '?'))[0]:5} / "
              f"陰 {got.get(a, ('?', '?'))[1]:5}"
              f"{'' if ok else f'  提案所載 {EXPECT[a]}'}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
