# 四個分析集合

`set_310_groundtruth.csv`（全集）、`set_187_precutoff.csv`、`set_81_postcutoff.csv`、
`set_49_main.csv`（主要分母）。由 `export_sets.py` 導出：

    python3 -I sets/export_sets.py measurements <輸出目錄>

⚠ **兩個參數不可省**——腳本的預設路徑指向本版本庫沒有的目錄。附 13 條 assert，
全部通過才寫出。

---

## 一、計數有兩個層級，兩個都對，混用才錯

| 層級 | 總數 | DP | DQ | DR | H2-IA | H2-IE | BL | BAT |
|---|---|---|---|---|---|---|---|---|
| **ligand-row**（CSV 的列） | 310 | 15 | 66 | 156 | 45 | 17 | 7 | 4 |
| **相異 PDB entry** | 278 | 14 | 56 | 140 | 44 | 16 | 5 | 3 |

引用任何比率時必須標明用的是哪一層。另有兩個在其他文件中出現的層級：
人類 `locus_of` 層 222，以及預先登記表 1 的過濾層 **268 組／240 結構**
（310 扣除 core 為空或含非標準殘基的 41 列、與 9A0U 的 2 列，兩者重疊 1 列，淨扣 42）。

⚠ **BL 是雞**（Gallus gallus 的 B-locus class II／BL-2），**BAT 是蝙蝠**。
本版本庫收錄的 291 個複合體中沒有任何牛（BoLA）的 class II 結構。
⚠ 這是**結構層**的陳述；`PREREGISTRATION.md:42` 提到的「牛」指 Racle 2023 在
**ligand/motif 層**涵蓋的 BoLA-DR 等位基因，兩者層級不同，並不衝突。

## 二、`reverse` 欄有三筆是偽陽性，不可直接統計

CSV 的 `reverse=True` 共 7 筆，其中**僅 4 筆可採信**：

| PDB | locus | 採信？ | 理由 |
|---|---|---|---|
| 3PGC | DR | ✓ | HLA-DR1＋CLIP106-120，X-ray 2.66 Å，RCSB 標題即「flipped peptide orientation」 |
| 4AEN | DR | ✓ | 同上，(G4S)₃ **共價接 α 鏈的人工構築**，2.20 Å。⚠ 預先登記 :112 自述其反向為設計所致，引用時須與自由胜肽的 3PGC 分開 |
| 7T6I | DP | ✓ | HLA-DP1＋CMV pp65，2.30 Å |
| 7ZFR | DP | ✓ | HLA-DP（DPA1\*02:01-DPB1\*01:01），2.90 Å |
| 9A0U | DR | ✗ | `method=Integrative`（對接模型，非實驗結構）；且被判為反向的鏈其序列為 `GILGFVFTL`，那是 MHC **class I** 表位，被誤派到 class II 框架。該結構真正的 class II 胜肽在另一條鏈且判為正向 |
| 9RUF | DQ | ✗ | 胜肽建模為 poly-UNK（序列未知），register 無定義 |
| 9RUH | DQ | ✗ | 同上 |

扣除後，HLA-DP 的結構層反向率為 **2/14 ＝ 14.3%**（相異 entry 層）。
⚠ HLA-DR 與 HLA-DQ 的對應數字**不在此並列**——理由見下方第三節。

⚠ **同批的 9RUG 條件完全相同（同日釋出、同樣 poly-UNK）卻被判為正向。** 三筆並列即可
看出：方向判定在退化輸入上的輸出是隨機的。使用本欄時必須先套上述排除規則，
或改用附有「不可判定」第三類輸出的版本。

## 三、DR 的兩筆反向不代表 HLA-DR 會反向呈現

3PGC 與 4AEN 是**同一條胜肽**（CLIP106-120）在**同一個等位基因**（DRB1\*01:01）上，
出自**同一研究線**（Günther et al. PNAS 2010；Schlundt et al. J Mol Biol 2012，
七位作者中六位重疊）。CLIP 是 invariant chain 的佔位胜肽且序列近似對稱，
原始文獻稱其為容許雙向結合的特殊情況。質譜洗脫配體層的分析（涵蓋 43 個 HLA-DR
等位基因、偵測下限約 5%）未測得 HLA-DR 的反向 motif。

→ 可據此陳述「class II 的開放式溝槽在原理上容許 C→N 方向結合」；
**不可**據此報出 HLA-DR 的反向率，或與 HLA-DP 的 2/14 並列比較。
