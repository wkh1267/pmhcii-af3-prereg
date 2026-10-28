# pMHC-II register 與結合方向：提案階段的量測與預先登記

碩士論文「AlphaFold3 於 peptide–MHC class II 的結構準確度與結合判別，以 HLA-DP 為焦點」
提案階段（2026 年 10 月繳交）已執行量測的腳本、PDB 代碼清單與輸出表格，連同預先登記。

公開的目的是讓提案第四章引用的每一個數字都能被第三方重算，並為提案送出時已凍結的判準
留下時間戳記。

## 內容

| 目錄／檔案 | 內容 |
|---|---|
| [`PREREGISTRATION.md`](PREREGISTRATION.md) | 繳交版第三章「研究方法」全文，逐字收錄：操作型定義、計算流程與閘門規範、三個目標的設計、主要終點階層。送出後不再更動 |
| [`ERRATA.md`](ERRATA.md) | **兩條**勘誤：(1) 反向組的建模路徑在本研究的實測約束設定下未能指定 register（符合率 0/15）；⚠ 詳情與適用範圍見該檔，不宜據此斷言 Boltz-2 無法用於此目的，(2) 提案圖 5 所需的逐筆資料不存在於公開來源，已停止引用 |
| `measurements/` | 量測的腳本、資料、黃金檔與 291 個 mmCIF。內有自己的 `README.md` 與 `SOURCES.md` |
| [`sets/README.md`](sets/README.md) | 四個分析集合的使用須知：兩個計數層級、`reverse` 欄的三筆偽陽性與排除理由、HLA-DR 的兩筆反向為何不可外推 |
| `sets/` | 四個分析集合的逐筆 PDB 代碼清單，與產生它們的腳本 |

## 提案第四章的五項量測對應到哪些腳本

全部在 `measurements/scripts/`。

| 提案節次 | 量測 | 腳本 | 可重算？ |
|---|---|---|---|
| 4.1 | ground-truth 對照表與文獻一致 | `register.py`、`holdout_leakage.py`、`census.py` | ✅ |
| 4.2 | 界面接觸 44.3% 由胜肽主鏈介導 | `contacts.py`、`contacts_fixed.py`、`aggregate.py` | ✅ 見 `measurements/verify.sh` |
| 4.3 | register 平移一格的幾何代價 3.80 Å | `register_shift.py` | ✅ |
| 4.4 | 既有 class II benchmark 的分層 | `ko_locus.py` | ⚠ 只有 locus 分層可重算；圖 5 的十八個值已撤回，見 `ERRATA.md` |
| 4.5 | HLA-DP 的結構稀缺有結構性成因 | `dp_diversity.py`、`census.py`、`fetch_imgt.py` | ✅ |

## 怎麼重現

界面接觸那一項附有三層驗證，只需 `python3` 與 `numpy`：

```sh
cd measurements && ./verify.sh
```

第 1 層由輸出 JSON 重算彙總（不需結構檔）；第 2 層自 291 個 mmCIF 重新推導以重現提案
v5 的引用值；第 3 層套用融合構築的截斷修正以重現繳交版採用的數值。可重現性與正確性是
兩個獨立的主張，故兩層都保留。

### ⚠ 已知：第 4 層會失敗，腳本會在該處停止（退出碼 1）

**第 1 至 3 層全部通過**，它們認證的正是提案 4.2 節引用的兩個數字（主鏈介導 44.3%、
P1/P4/P6/P9 側鏈 29.7%）。

第 4 層與第 5 層的逐 locus 比對所用的黃金檔，凍結於一次 locus 更正之前：8JR4、8JRJ、
8JRK 三個條目為蝙蝠 MHC II，原先被誤推導為 HLA-DR。更正以具名清單實作於
`scripts/locus_fix.py`，**修在使用點而不動資料檔與黃金檔**，是當時刻意的決定——資料本身
沒錯（物種欄是對的），錯的是 locus 的推導，而重生黃金檔會抹掉那次更正的稽核痕跡。

因此這兩層所報的差異是已知且完整解釋的，逐項為：

- 多一行揭露「已套用 3 筆 locus 更正」
- DR 由 185 鏈 / 109 條目 降為 **181 鏈 / 106 條目**（移出 4 條鏈、3 個條目）
- 多一個 BAT 區塊，即移出去的那 3 個條目
- 其餘各列（ALL、DQ、DP）的差異皆為上述重新分組的連帶結果

**現行腳本的輸出是正確的，繳交稿引用的也是現行值**；舊的是黃金檔。想看第 5 層，
請逐層手動執行 `scripts/` 下的對應腳本。

其餘腳本各自以參數吃入 `measurements/data/` 下的檔案，用法見每支腳本的 docstring。
四個分析集合可由 `python3 -I sets/export_sets.py measurements <輸出目錄>` 重新導出，
該腳本附 13 條 assert，全部通過才寫出。
⚠ **兩個參數不可省**：腳本的預設路徑指向本版本庫沒有的 `proposal/`，無參數會失敗。

## 不在這個版本庫裡的東西

- **提案全文**。本版本庫只公開預先登記所需的三組內容，不含導論、文獻回顧、初步結果、
  預期結果與進度表
- **參考文獻的 PDF**。有版權，不可轉散布
- **AlphaFold3 與 Boltz-2 的輸出**。提案階段尚未執行；AlphaFold3 的輸出另受其
  `OUTPUT_TERMS_OF_USE.md` 約束，不能以 CC-BY 或 CC0 釋出

`measurements/` 底下的說明文字偶爾會提到 `proposal_v6`、`slides/` 這類名稱，那是未公開
的內部工作文件。這些檔案連同其雜湊值一併原樣保留，以免動到 `verify.sh` 所比對的黃金檔；
讀者可逕行忽略那些指涉。

## 來源與授權

各項資料的取得方式、取用日期與已知限制逐條記在
[`measurements/SOURCES.md`](measurements/SOURCES.md)。授權分三層，不以單一條款涵蓋全部：

| 範圍 | 授權 |
|---|---|
| 腳本（`measurements/scripts/`、`sets/export_sets.py`、`measurements/verify.sh`） | MIT，見 [`LICENSE`](LICENSE) |
| 量測輸出與清單（`measurements/data/`、`measurements/expected*/`、`sets/*.csv`）與本版本庫的文件 | CC BY 4.0 |
| `measurements/cif/` 的 291 個 mmCIF | 取自 RCSB PDB，公有領域（CC0），非本研究所有 |

分三層而不一刀切是刻意的：AlphaFold3 的 `OUTPUT_TERMS_OF_USE.md` 第 5.2 條禁止對其輸出
附加衝突條款，而 CC BY 與 CC0 皆允許商業使用，兩者直接衝突。本版本庫目前不含任何
AlphaFold3 輸出；日後若加入，該部分須另立聲明，不受上表涵蓋。
