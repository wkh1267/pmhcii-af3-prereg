# pMHC-II 界面接觸統計 — 驗證包

打包日期 2026-10-03。對應碩論提案 §01「機制：為什麼『結構合理』在 class II 上幾乎是
免費的」所引用的兩個數字。

⚠ **2026-10-03 更新：獨立審查發現一個 bug，兩個引用值都已修正。**
本包現在保留**兩套基準**：原腳本（重現 v5 引用值）與修正後腳本（v6 採用值）。
三層驗證全部必須通過。詳見下方「兩套基準」一節與 `../../contacts_audit_20261003.md`。

## 要驗證什麼

提案主張：class II 以一條與序列無關的保守主鏈氫鍵梯固定 peptide，因此將近一半的
界面訊號在任何 peptide 上都會被滿足，「幾何看起來正確」在本系統中幾乎是免費的。

支撐這個主張的量測：

| 量 | v5 引用值（已作廢） | **v6 採用值** |
|---|---|---|
| 胜肽**主鏈**（N/CA/C/O/OXT）介導的接觸佔全部界面接觸 | 47.0%（中位 45.7%） | **44.3%**（中位 43.8%，IQR 39.9–48.0%） |
| **P1/P4/P6/P9 側鏈**佔全部界面接觸 | 23.7% | **29.7%** |
| P1/P4/P6/P9 佔側鏈接觸 | 42.4% | **53.4%** |

## 兩套基準

`contacts.py` 未對共價融合構築（胜肽–linker–β 鏈單鏈）套用胜肽長度截斷。360 條鏈中
有 82 條（23%）因此把整條 β 鏈當成胜肽，使 β–α 的蛋白質內部接觸被計入界面統計。
受影響條目的接觸總數中位為 20,310，分離胜肽鏈則為 189。污染在 locus 間不均勻，
**HLA-DP 受影響 13/27（48%）最嚴重**。

| 目錄／腳本 | 回答的問題 |
|---|---|
| `scripts/contacts.py` + `expected/` | 這個包能不能重現 v5 當初發表的數字？ |
| `scripts/contacts_fixed.py` + `expected_fixed/` | 修正後的數字是否正確？ |

**兩者都保留，因為可重現性與正確性是兩個獨立的主張，兩個都要能被檢驗。**
把黃金檔直接換掉，第一個問題就再也問不出來了。

`contacts_fixed.py` 相對 `contacts.py` 只有三個 hunk，全部標註 `[FIX 2026-10-03]`：
建立 `(pdb, entity_id) → plen` 對照表、胜肽側截斷、以及 MHC 側反向排除胜肽段
（融合構築的 β entity 就是融合鏈本身，不排除的話胜肽會與自己量距離）。
輸出寫到 `data/contacts_fixed.json`，不覆蓋原檔。

## 一個必須揭露的限制

`contacts.py` 對**每一條胜肽鏈**各輸出一筆，因此同一個 PDB entry 裡的晶體複本會被
重複計入。實際上是 **211 個 PDB entry 的 360 條鏈**，其中 109 個 entry 貢獻 ≥2 條
（最多一個貢獻 6 條）。這些複本不是獨立觀測。

`scripts/aggregate.py` 以三種彙總層級重算，確認此事不改變結論：

| 彙總層級 | n | 主鏈 mean | 主鏈 median | anchor mean |
|---|---|---|---|---|
| 每條胜肽鏈各算一筆（原始輸出） | 360 | 47.0% | 45.7% | 23.7% |
| 先彙總到 (結構, peptide) 再平均 | 237 | 46.5% | 45.2% | 25.0% |
| 先彙總到 PDB entry 層級再平均 | 211 | 46.8% | 45.7% | 24.8% |

三者差距 0.5 個百分點，遠小於 47% 與 24% 的對比。**但論文應改報
(結構, peptide) 層級的 46.5% / 25.0%，n=237** —— 與 register 分析的單位一致，
且避免「n=360 但只有 211 個結構」這個會被審查抓到的矛盾。

## 如何驗證

    sh verify.sh

只需 `python3` 與 `numpy`。分三層，**任一層不符即以非零碼結束**（約 50 秒）：

1. **不需結構檔，數秒。** 從 `data/contacts.json` 重算三種彙總，與
   `expected/aggregate.json` 逐欄比對。
2. **原腳本，從原始結構重新推導。** 以 `cif/` 下的 291 個 mmCIF 重跑
   `contacts.py`，與 `expected/contacts_stdout.txt` 比對 stdout，
   並確認重新產生的 `contacts.json` 與原檔逐筆相同。→ 應得 47.0% / 23.7%
3. **修正後腳本，同上。** 重跑 `contacts_fixed.py`，與
   `expected_fixed/contacts_fixed_stdout.txt` 比對。→ 應得 44.3% / 29.7%

## 方法

以 1DLH（HLA-DR1 + HA306–318）為參考，對每個結構做 α 鏈與 β 鏈的序列比對，
取 α1（前 84 殘基）與 β1（前 95 殘基）溝槽結構域的配對 Cα 做 Kabsch 疊合
（RMSD ≤ 3.0 Å 才採納），據此定出 binding register 與結合方向（`register.py`）。
接著對每條胜肽鏈取最近的 α/β 鏈實例，計算 4 Å 內的重原子接觸，依胜肽原子是否
屬於 N/CA/C/O 分為主鏈與側鏈，並將側鏈接觸按 core 位置 P1–P9 分解（`contacts.py`）。
反向結合者的 percore 已反轉為 P1→P9 方向。

## 內容

    verify.sh                 三層驗證的進入點
    scripts/register.py       register 與結合方向判定（contacts.py 的前置）
    scripts/contacts.py       界面接觸統計 —— 被驗證的對象（原版，未修正）
    scripts/contacts_fixed.py 同上，套用 fused 截斷修正（v6 採用）
    scripts/aggregate.py      彙總層級檢查 —— 不影響上游計算
    data/rows7.json           291 個結構的 entity 中介資料（plen 的來源）
    data/registers.json       每條鏈的 register、方向、疊合品質
    data/contacts.json        原版輸出
    data/contacts_fixed.json  修正後輸出
    cif/                      291 個 class II 結構的 mmCIF（gzip，約 86 MB）
    expected/                 原版的預期輸出、環境資訊、SHA256
    expected_fixed/           修正後的預期輸出

### 與原始腳本的唯一差異

`scripts/contacts.py` 與 `scripts/register.py` 只把寫死的絕對路徑改為相對於
腳本位置（在首個 import 後插入三行 `_HERE` / `_ROOT`，並改寫四處檔案路徑），
好讓本包解壓到任何位置都能執行。**計算邏輯一字未動。**

- `expected/SHA256SUMS.txt` — 本包所有檔案的雜湊，用來確認包的完整性：
  `cd <解壓目錄> && sha256sum -c expected/SHA256SUMS.txt`
- `expected/ORIGINAL_SCRIPT_HASHES.txt` — 改路徑**之前**原始腳本的雜湊，
  供與上游副本對照。兩者的 diff 應只包含上述路徑改動。

### 殘留限制（2026-10-03 新增）

1. **plen 切得不完全乾淨。**77 個 `fused_peptide` 條目中，有 11 個的 pep 字串在 plen
   窗口內仍含得出 linker motif（GGGGS / GGGS / GSGS / GGSG）。例如 1FNE plen=20、
   pep=GKKVITAFNDGLKGGGGSLV，linker 自第 14 位起。那些條目的分母仍含少量 linker 殘基。
   影響有限：修正後 fused 組主鏈 46.0% / anchor 27.5%，sep 組 43.8% / 30.4%，
   差 2.2 與 2.9 個百分點；且截斷後 fused 的觀測長度中位為 12、sep 為 13，
   代表多數 linker 殘基本來就沒被建模。**不另外處理**，因為收緊截斷會破壞與
   `register.py:138` 的一致性 —— 兩邊 plen 不同會讓 P1–P9 整體位移且不報錯。
2. **`expected/SHA256SUMS.txt` 於 2026-10-03 重算**，以涵蓋本次新增的檔案。
   原始分析腳本的雜湊記錄在 `expected/ORIGINAL_SCRIPT_HASHES.txt`，未更動。

## 附加驗證：參考框架的敏感度（2026-10-06）

```
python3 scripts/frame_sensitivity.py      # 約 90 秒
```

`register.py` 以 1DLH（HLA-DR1）定義溝槽參考框架。本腳本改以 7ZAK（HLA-DP，
獨立胜肽鏈、正向、p1_dist 0.94 Å）重建框架，重新指認全部結構的 P1 與結合方向，
逐一比對。`register.py` 的 3.0 Å 疊合門檻在此**刻意不套用**（已依 proposal_v6
定義八撤銷）。

結果：441 個可比對實例中 438 個（99.3%）完全相同，**HLA-DP 31/31 無一改變**。
三筆差異全部來自 8VSP（DQ）同一不對稱單位的三個複本，且為相鄰一格的位移
（MRMATPLLM 對 RMATPLLMQ）。

結論：DR1 框架對 HLA-DP 的 P1 指認沒有系統性偏誤。3WEX 被排除的原因是門檻，
不是框架。

## 逐 locus 分層（2026-10-07，verify.sh 第 4 層）

```
python3 scripts/percore_by_locus.py       # 數秒，不需結構檔
```

`contacts_fixed.py` 的 stdout 只報全體平均，而全體由 DR 主導（185/360 條），
因此看不出 DP 的分布與 DR 是兩個形狀。本腳本把同一份 `contacts_fixed.json`
按 locus 切開重算，聚合方式不變（mean-of-ratios）。

locus │  P1  │  P4  │  P6  │  P7  │  P9  │ 四口袋 │ n（鏈／entry）
──────┼──────┼──────┼──────┼──────┼──────┼────────┼───────────────
DR    │ 12.7 │ 15.2 │ 10.3 │  9.7 │ 14.0 │ 52.3%  │ 185 ／ 109
DQ    │ 10.5 │ 14.1 │ 15.9 │  9.3 │ 18.7 │ 59.2%  │  84 ／  50
DP    │ 12.2 │ 11.0 │ 26.6 │ 11.4 │  9.4 │ 59.1%  │  27 ／  13

DP 的 P6 是 DR 的 2.6 倍，P9 反而最低。P7 在 DP 排第三，高於被納入的 P4 與 P9
（逐鏈 P7 > P9 18/27、P7 > P4 14/27；entry 去重後仍排第三）。

**這不是把 P7 納入閘門的理由。** 接觸數量的是「有多少 MHC 重原子落在 4 Å 內」，
專一性量的是「換掉這個位置的殘基會不會改變結合」。兩者不等價，而 DP 的 motif
與結構文獻（Racle 2023、Sidney 2010、Klobuch 2022、Dai 2010）一致不把 P7 當
anchor。四口袋的依據取自那些文獻，不是本表的排序 —— 詳見 proposal_v6 §01
「四口袋在 DP 上的依據」。

本表的用途有二：(1) 佐證 pipeline —— 在未使用任何 motif 資料的情況下重現了已知的
DP anchor 階層（P1/P6 主導、P9 偏弱）；(2) 延伸主軸論證 —— 不只主鏈的 44.3%
是免費的，側鏈接觸裡也有一部分（DP 的 P7）不帶分辨力。

限制：DP 僅 13 個 unique entry；PDB 不是 DP ligand 的隨機樣本；DP 亦是當初 fused
截斷污染最嚴重的一組（13/27 條），本表為修正後資料，但 DP 對該修正最敏感。

## 第三套基準：撤銷疊合門檻（2026-10-07，verify.sh 第 5 層）

```
SUP_RMSD_MAX=inf python3 scripts/register.py ../data/registers_nothresh.json
REGISTERS=../data/registers_nothresh.json CONTACTS_OUT=../data/contacts_nothresh.json \
  python3 scripts/contacts_fixed.py
```

`register.py:154` 原為裸寫的 `if rms > 3.0: continue`。提案定義八已撤銷該門檻
（3WEX 以 3.47 Å 被排除，而其 P1 指認正確）。門檻現已參數化為環境變數
`SUP_RMSD_MAX`，**預設 3.0，輸出與 `data/registers.json` 位元相同**，
既有兩套基準不受影響。

### 為什麼撤銷：門檻沒有鑑別力

取未截斷的 452 筆判定比較兩組：

集合                  │  n  │ p1_dist 中位 │ 最差    │ >2 Å
──────────────────────┼─────┼──────────────┼─────────┼──────
通過舊門檻（≤3.0 Å）  │ 407 │ 0.53 Å       │ 44.28 Å │ 3 筆
被舊門檻排除（>3.0 Å）│  45 │ 0.77 Å       │  1.41 Å │ 0 筆

sup_rmsd 與 p1_dist 的 Pearson r = **+0.014**。舊門檻留下三筆 p1_dist > 2 Å 的
判定，卻排除了 45 筆最差僅 1.41 Å 的判定。若要設納入判準，變數應取 p1_dist。

### 重建後的數字

基準                     │ n   │ 主鏈         │ 四口袋佔全部 │ 佔側鏈
─────────────────────────┼─────┼──────────────┼──────────────┼────────
v5 引用值（原腳本）      │ 360 │ 47.0%        │ 23.7%        │ 42.4%
v6 採用值（門檻 3.0）    │ 360 │ 44.3%（43.8）│ 29.7%        │ 53.4%
撤銷門檻後的全集         │ 404 │ 44.4%（43.9）│ 29.8%        │ 53.6%

差 0.1 個百分點。納入 3WEX 與其餘 44 筆並未動搖 §01 的結論。

DP 格：n 由 27/13 增為 31/14（多的四筆即 3WEX 的四個複本），P6 26.6→26.7%、
P1 12.2→12.4%、P7 11.4→11.6%（仍排第三）、P9 9.4→11.4%、P4 11.0→10.2%。

三套基準全部保留：可重現性、正確性、納入範圍是三個獨立的主張，都要能被檢驗。
