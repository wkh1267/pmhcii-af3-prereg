"""資料夠不夠：逐 Aim 的結構需求盤點 [新增 2026-10-07]

起因：投影片前半花四頁說 DP 的結構有多缺，第十頁卻只列出來源清單，從未回答
「那到底夠不夠」。這支把答案算出來，並釘住每個分母的單位。

**「夠不夠」不是一個答案，是三個** —— 稀缺只卡住其中一件事：

    用途                        要的是什麼        稀缺的影響
    ──────────────────────────  ────────────────  ──────────────────────────────
    指派規範的 ground truth     涵蓋度            不受影響。每個結構各給一個
                                                  ground truth，不做跨組檢定
    AF3 結構準確度 benchmark    統計力            **卡住**。DP 臂 n=8
    反向結合的行為（Aim 3）     反向的 ground     **硬上限**。全 PDB 只有 4 個
                                truth             反向結構

單位必須分清楚，混用會產生假的不一致：
  · 「結構」= PDB entry
  · 「組」  = (結構, 胜肽序列)，同一個 entry 可含多條胜肽鏈
  例：DP 的 post-cutoff 是 **7 個結構 / 8 組**。提案 §03 寫 7、§04 寫 8，兩者都對。

    python3 scripts/data_sufficiency.py

⚠ 本腳本由 registers_nothresh.json（452 個實例）以 (pdb, pep) 去重重建，得 310 組。
  提案 Aim 1 記的是 236 組，其差異為 DR −39 / DQ −3，**DP 完全不受影響**（14 = 14，
  post-cutoff 8 = 8）。該缺口見 §08「236 仍未解釋」，M1 以記錄下來的查詢重新產生。
  故本表的 DP 欄可直接引用，DR/DQ 欄須標為工作值。
"""
import os, json, collections, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from locus_fix import locus_of, note as locus_note

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
CUT = "2021-09-30"          # AF3 訓練截止
# 方向判定已知的三個排除，見 §03 與 §16
EXCL_REV = {"9A0U": "class I 的 M1 胜肽誤納，orient −0.144（本工具唯一一次誤判，且自己標了出來）",
            "9RUF": "Fab–HLA cryo-EM，胜肽序列登記為 unknown",
            "9RUH": "同 9RUF"}

# ⚠ 讀 **nothresh**（未套疊合品質門檻）。registers.json 套了 SUP_RMSD_MAX = 3.0，
# 而 §02A 已把該參數（S1）標為「無依據、已證實有偏誤、**已撤銷**」。它在 DP 上唯一
# 砍掉的結構是 **3WEX —— 焦點 allele DPA1*02:02/DPB1*05:01**，套用後焦點 allele 的
# ground-truth register 數會變成 0。該門檻實際量的是「離 1DLH（HLA-DR1）參考框架多遠」：
# 淘汰率 DR 0.0% / DQ 5.4% / DP 12.9% / 小鼠 30–36% / 牛 85.7%。見 reports/
# DP_STRUCTURE_AUDIT_20261007.md §5.1 與 scripts/dp_diversity.py。
reg = json.load(open(os.path.join(ROOT, "data", "registers_nothresh.json")))
seen, rows = set(), []
for r in reg:                               # 去晶體複本：同結構同序列只算一次
    k = (r["pdb"], r["pep"])
    if k not in seen:
        seen.add(k); rows.append(r)

def locus(r):
    return locus_of(r)          # 套用已知的 locus 更正，見 scripts/locus_fix.py

def tally(sel):
    n = collections.Counter(locus(r) for r in sel)
    s = collections.defaultdict(set)
    for r in sel:
        s[locus(r)].add(r["pdb"])
    return n, s

post = [r for r in rows if r.get("rel") and r["rel"] > CUT]
nt, st = tally(rows)
np_, sp = tally(post)

print("結構資料的實際庫存　[本提案量測，可重現]")
print(f"  來源 data/registers_nothresh.json（{len(reg)} 實例 → 去 (pdb,胜肽序列) 重複後 "
      f"{len(rows)} 組 / {len({r['pdb'] for r in rows})} 結構）")
print(f"  holdout 切點 {CUT}（AF3 訓練截止）")
print(f"  {locus_note()}\n")
print(f"  {'locus':8s} {'總組數':>7s} {'總結構':>7s} {'post組':>7s} {'post結構':>9s}   反向結構")
for k in sorted(nt, key=lambda x: -nt[x]):
    rev = sorted({r["pdb"] for r in rows if locus(r) == k and r.get("reverse")})
    keep = [p for p in rev if p not in EXCL_REV]
    drop = [p for p in rev if p in EXCL_REV]
    tag = ", ".join(keep) if keep else "—"
    if drop:
        tag += f"　（排除 {', '.join(drop)}）"
    print(f"  {k:8s} {nt[k]:7d} {len(st[k]):7d} {np_.get(k,0):7d} "
          f"{len(sp.get(k,())):9d}   {tag}")

# ---- TCR 的有無 [2026-10-07 新增] ----------------------------------------
# Aim 1 的缺口陳述寫的是「**無 TCR 的** pMHC-II」—— 那正是與 Woods 2026 的 2.0 Å
# 區隔開來的那一點（該值量自三元複合體）。但 §04 記的 holdout「DP 8」含帶 TCR 者。
# 兩個數字必須分開報，否則主要賣點與主要數字互相矛盾。
#
# has_tcr 的判據（register.py:178）：該 PDB entry 的 `tcr`（獨立 TCR α/β 鏈）或
# `tcrfuse`（與 MHC 共價融合的 TCR 構築）任一非空。是 entry 層級屬性。
print("\nTCR 的有無　—— Aim 1 的缺口寫的是「無 TCR 的 pMHC-II」")
print(f"  {'locus':8s} {'post組':>7s} {'無TCR組':>8s} {'無TCR結構':>10s} {'其中正向':>9s}")
for k in sorted(np_, key=lambda x: -np_[x]):
    sel = [r for r in post if locus(r) == k]
    nt_ = [r for r in sel if not r.get("has_tcr")]
    fw = [r for r in nt_ if not r.get("reverse")]
    print(f"  {k:8s} {len(sel):7d} {len(nt_):8d} {len({r['pdb'] for r in nt_}):10d} "
          f"{len({r['pdb'] for r in fw}):9d}")
print("\n  DP 的 post-cutoff 逐筆：")
for r in sorted([r for r in post if locus(r) == "DP"], key=lambda x: x.get("rel") or ""):
    print(f"    {r['pdb']}  {r.get('rel')}  TCR={'有' if r.get('has_tcr') else '無'}"
          f"  {'反向' if r.get('reverse') else '正向'}  {r.get('core')}")
print("""
  → 照 Aim 1 自己寫的缺口（無 TCR）做，DP 臂是 **4 組 / 4 結構**，不是 8 / 7。
    再限定正向，只剩 **2 個結構**（7T2A、7ZAK）。
    §03 對總數誠實（「14 已含帶 TCR 者；無 TCR 僅 10」），但 holdout 的 8
    從未做過同樣的拆分。建議：TCR 有無改為**分層報告的共變項**，不是過濾器 ——
    與 §09 處理 assay type 的方式一致。""")

print("\n逐 Aim：要多少、有多少、夠不夠")
print("""
  Aim 1 主交付 · 指派規範的 ground truth
      要的是涵蓋度（每個結構各給一個 ground truth），不做跨組檢定 → **夠**
      可用全部 278 組，含小鼠 H2-IA/H2-IE 41 組 —— 規範本身不限人類

  Aim 1 benchmark · AF3 結構準確度，post-cutoff
      要的是統計力。DP 臂 n = 8 → **只能描述性**
      提案 §04 已預先宣告：locus 三分類的顯著性檢定只在 DR 與 DQ 之間做；
      DP 只報點估計與 CI。主要分析改用 FoldSeek qtmscore 作**連續**代表性
      變數，全部 holdout 進同一個模型 —— 這是繞開 n=8 的設計，不是補救

  Aim 3 · AF3 面對反向結合的行為
      要的是反向的 ground truth。全 PDB 只有 4 個（3PGC、4AEN、7T6I、7ZFR）
      → **硬上限**，提案記 n 上限 9（晶體 2 + Racle 已驗證 epitope 7）
      這個上限不是本提案的選擇，是世界上就只有這麼多

  Aim 2、Aim 4 · binder 判別與消融
      完全不用結構的 ground truth，用序列資料 → **不受結構稀缺影響**
""")

# n=8 的區間寬度（Wilson，α=0.05）—— 說明「描述性」有多描述性
import math
z = 1.959963985
def wilson(k, n):
    p, d = k / n, 1 + z*z/n
    c = (p + z*z/(2*n)) / d
    h = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / d
    return c-h, c+h
print("  DP 臂 n = 8 的 95% Wilson 區間寬度（提案 §04 記「最寬達 57 個百分點」）")
for k in range(9):
    lo, hi = wilson(k, 8)
    print(f"    {k}/8 = {k/8:4.2f}   [{lo:.2f}, {hi:.2f}]   寬 {100*(hi-lo):4.1f} 個百分點")
w = max(100*(hi-lo) for hi, lo in [(wilson(k,8)[1], wilson(k,8)[0]) for k in range(9)])
print(f"  → 最寬 {w:.0f} 個百分點，與 §04 的 57 相符")
