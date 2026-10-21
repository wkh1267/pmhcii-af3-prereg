# -*- coding: utf-8 -*-
"""locus 標註的已知更正 [新增 2026-10-09]

**問題**：`rows7.json` 的 locus 由 α 鏈對 DRA / DQA1 / DPA1 三條人類參考序列的 identity
指派，**沒有物種守門**。蝙蝠 MHC II 在三者中最接近 DRA，於是被標成 `DR`。

**影響**（已量化，見 reports/DP_STRUCTURE_AUDIT_20261007.md 與 §16 待辦第 1 項）：
    §01 的 DR 列       185 條 → 人類 181 條；四口袋 52.3% → 52.0%
    §01 全體           404 條 → 400 條；主鏈 44.40% → 44.46%
    post-cutoff DR     41 組 / 35 結構 → **37 組 / 32 結構**
    再套無 TCR         30 結構 → 27 結構
數值差 ≤0.3 pp，但那三筆**全是 post-cutoff**，因此進了 benchmark 的 DR holdout 臂；
而蝙蝠自己的形狀與人類 DR 不同（P4 24.6% 對 15.0%），不是可忽略的同類。
**DP 完全不受影響。**

**為什麼在使用點修而不改資料檔**：`rows7.json` 在 `expected/SHA256SUMS.txt` 的黃金檔內，
改它會讓四層驗證全部失效。資料本身沒錯（物種欄是對的），錯的是 locus 的推導。
故在此集中修正，所有報 locus 分層數字的腳本一律改用 `locus_of()`。

**為什麼是具名清單而不是通用的物種守門**：只有三筆已知誤標、成因已查明。寫一個通用
的物種→locus 規則需要對全部 291 個 entry 驗證，而本提案沒有做那件事。具名清單誠實、
可稽核，且明確記錄了改了什麼。新增誤標時往下面加，並在 §16 記一筆。

    from locus_fix import locus_of
    lab = locus_of(record, rows)       # record 需有 'pdb'；rows 為 rows7 的 dict
"""

# pdb → (錯誤標註, 更正後, 理由)
KNOWN = {
    "8JR4": ("DR", "BAT", "蝙蝠 MHC II，α 鏈最接近 DRA 而被誤標"),
    "8JRJ": ("DR", "BAT", "同上"),
    "8JRK": ("DR", "BAT", "同上"),
}


def locus_of(rec, rows=None):
    """取一筆紀錄的 locus，套用已知更正。rec 可為 registers/contacts 的列或 pdb 字串。"""
    pdb = rec if isinstance(rec, str) else rec.get("pdb")
    if pdb in KNOWN:
        return KNOWN[pdb][1]
    if isinstance(rec, str):
        if rows is None:
            raise ValueError("傳 pdb 字串時必須一併給 rows")
        rec = rows[pdb]
    loci = rec.get("loci")
    return (loci or ["?"])[0]


def note():
    """一行說明，供腳本在輸出中揭露。"""
    n = len(KNOWN)
    ids = "、".join(sorted(KNOWN))
    return (f"⚠ 已套用 {n} 筆 locus 更正（{ids}：蝙蝠 MHC II 被誤標為 DR）。"
            f"見 scripts/locus_fix.py。DP 不受影響。")
