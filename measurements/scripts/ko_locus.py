#!/usr/bin/env python3
"""Ko et al. 2026 benchmark 的 locus 分層 —— 原文無此分層，本提案推導。

來源：Ko et al. preprint（bioRxiv 2024.10.06.616783）的補充檔 media-1.xlsx。
⚠ 不是期刊版（BMC Med Genomics 2026;19(Suppl 1):120）的 Additional File 1——
  期刊版欄位不同（多 IEDB 負例與 Epitope ID、沒有 α/β 鏈欄），本支依賴 α/β 欄，
  餵期刊版會先在下方第 26–27 行的物種 assert 爆掉。詳 SOURCES.md 第 7 列。
      doi:10.1186/s12920-026-02432-4

用法：python3 -I ko_locus.py <path/to/media-1.xlsx>

論文自印的 51 human / 9 mouse / 12 chimeric = 72 以 assert 當作解析正確性的錨點。
"""
import collections
import re
import sys

import openpyxl

LOCUS = re.compile(r"HLA-D([PQR])")


def main(path):
    ws = openpyxl.load_workbook(path, read_only=True, data_only=True).worksheets[0]
    rows = [r for r in ws.iter_rows(values_only=True) if r[0]][1:]

    species = collections.Counter(r[3] for r in rows)
    assert len(rows) == 72, f"總數應為 72，得到 {len(rows)}"
    assert species["Homo sapiens"] == 51 and species["Mus musculus"] == 9 \
        and species["Chimeric"] == 12, f"物種組成與論文自印不符：{dict(species)}"

    # 一筆可能兩條鏈同 locus，取 set 避免重複計數
    loci = {r[0]: sorted({"D" + m for m in LOCUS.findall(f"{r[1] or ''} {r[2] or ''}")})
            for r in rows}

    human = collections.Counter()
    mixed = []          # 非 Homo sapiens 卻帶人類 HLA 鏈者
    for r in rows:
        if not loci[r[0]]:
            continue
        if r[3] == "Homo sapiens":
            human.update(loci[r[0]])
        else:
            mixed.append(r)

    print(f"全集 {len(rows)} 筆：", dict(species))
    print()
    print("只數 Homo sapiens（本提案採用）")
    for k in ("DR", "DQ", "DP"):
        print(f"  {k}  {human[k]:>2} / 72 = {human[k] / 72:.1%}"
              f"   （佔人類 51 筆的 {human[k] / 51:.1%}）")
    print(f"  合計 {sum(human.values())}，與人類 51 筆相符"
          if sum(human.values()) == 51 else "  ⚠ 合計對不上")
    print()
    print("⚠ 若把帶人類 HLA 鏈的 chimeric 一併計入，會多出這 "
          f"{len(mixed)} 筆（即 DR 35 / DQ 18 的來源）：")
    for r in mixed:
        print(f"    {r[0]}  α={r[1]}  β={r[2]}  [{r[3]}]  -> {','.join(loci[r[0]])}")
    print()
    dp = [r[0] for r in rows if "DP" in loci[r[0]]]
    print("DP 的那一筆：", ", ".join(dp))
    for r in rows:
        if r[0] in dp:
            print(f"    α={r[1]}  β={r[2]}")

    # 負例的構成（論文未在正文說明）
    ana = sum(1 for r in rows if sorted(r[4]) == sorted(r[5]))
    print()
    print(f"負例為對應正例的胺基酸重排（anagram）：{ana} / {len(rows)}")
    print(f"胜肽含 GS 接頭者：{sum(1 for r in rows if 'GGGS' in r[4] or 'GSGG' in r[4])} 筆")


if __name__ == "__main__":
    main(sys.argv[1])
