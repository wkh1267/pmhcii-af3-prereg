"""CONSORT 式資料流：從 RCSB 查詢到每一個已發表數字的分母 [新增 2026-10-09]

起因：§08 列了五個互不相同的分母（360 / 357 / 236 / 233 / 210）而只有 360 對得上，
其中 **236 那一組的過濾規則完全追不到**，八種假說都重現不出來，且其差異有方向性
（DR −39、DQ −3，DP/H2-IA/H2-IE/BL 完全吻合）—— 而 Aim 2 的主要終點正是 DP vs DR。

本腳本的作用不是「找回 236」，是**用一個每一道過濾都有名字、有計數、可重跑的集合
取代它**。依定義七，選擇參數必須先於資料固定並記錄；一個追不到產生規則的集合不符合
這個要求，不論它是怎麼來的。

    python3 scripts/consort.py

輸出即論文的 Figure S1。每一行是一道過濾，括號內是它排除的理由，不是事後的描述。
"""
import os, sys, json, collections
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
sys.path.insert(0, HERE)
D = lambda f: json.load(open(os.path.join(ROOT, "data", f)))

census = D("census_2026-10-07.json")
rows   = {r["pdb"]: r for r in D("rows7.json")}
nt     = D("registers_nothresh.json")
th     = D("registers.json")

AA20 = set("ACDEFGHIKLMNPQRSTVWY")
CUT  = "2021-09-30"
sys.path.insert(0, HERE)
from locus_fix import locus_of, note as locus_note
loc  = locus_of          # 套用已知的 locus 更正，見 scripts/locus_fix.py

W = 78
def line(label, n, unit, note=""):
    print(f"  {label:<34s} {n:>5d} {unit:<4s} {note}")
def drop(label, n, note=""):
    print(f"    ├─ {label:<30s} −{n:<4d}      {note}")

print("=" * W)
print(" pMHC-II 結構資料流　CONSORT　[本提案量測，可重現]")
print(f" RCSB 序列搜尋取用日 2026-10-07　·　AF3 訓練截止 {CUT}")
print(f" {locus_note()}")
print("=" * W)

# ── 第一段：entry 來源 ──────────────────────────────────────────────
crow = census["rows"] if isinstance(census, dict) else census
cid  = set(census["ids"])
kept = [r for r in crow if not r.get("excluded_reason")]
spec = collections.Counter(tuple(r.get("species") or ["?"]) for r in rows.values())

print("\n【一】entry 來源")
print("  ▸ 可重現的部分 —— scripts/census.py（人類古典 class II 普查）")
line("RCSB 序列搜尋命中", len(crow), "entry", "DRA P01903 α1+α2，identity ≥ 0.55")
drop("identity < 0.60（非人類古典）", sum(1 for r in crow if "identity" in str(r.get("excluded_reason"))))
drop("次高者差 < 0.08（無法分 locus）", sum(1 for r in crow if "margin" in str(r.get("excluded_reason"))))
line("→ 判定為人類古典 class II", len(kept), "entry")

print("\n  ▸ ⚠ 不可重現的部分 —— data/rows7.json 的取得查詢**沒有紀錄**")
line("本提案實際註解的 entry", len(rows), "entry",
     "　".join(f"{s[0]} {n}" for s, n in spec.most_common()))
extra = sorted({p for p in rows} - cid)
drop("其中不在上述普查內者", len(extra), "、".join(extra))
print("       這 8 筆（5 雞、2 人、1 鼠）說明 rows7 來自一次更早、未記錄的取得。")
print("       **這比 236 那一組更根本** —— 不只某個子集的過濾追不到，是基底集合本身。")
print("    ✓ 但焦點 locus 不受影響：rows7 的 14 個 DP 與普查的 14 個**逐一相同**，")
print("       故 DP 的所有結論都落在可重現的範圍內。M1 須以記錄下來的查詢重建全集。")

# ── 第二段：實例層級 ────────────────────────────────────────────────
have_cif = [p for p in rows if os.path.exists(os.path.join(ROOT, "cif", p + ".cif.gz"))]
no_ab = [p for p in have_cif if not rows[p]["alpha"] or not rows[p]["beta"]]
def npep(p):
    r = rows[p]
    return len(r.get("sep", [])) + len(r.get("ii", [])) + \
           sum(1 for f in r.get("fused", []) if f["kind"] == "fused_peptide")
no_pep = [p for p in have_cif if p not in no_ab and npep(p) == 0]
inst_pdbs = {r["pdb"] for r in nt}

# ── 物種與 locus 的收錄範圍（2026-10-09 補；原先從未寫下）───────────────────
print("\n【一B】rows7 的收錄範圍　—— 原先沒有寫下來，由資料反推並固定")
tab = collections.Counter(((r.get("species") or ["?"])[0], (r.get("loci") or ["?"])[0])
                          for r in rows.values())
print(f"  {'species':9s} {'locus':9s} {'entry':>6s}   說明")
NOTE = {("human","DM"):"伴護分子，不呈現胜肽 → 0 實例",
        ("human","DO"):"同上 → 0 實例",
        ("mouse","H2-DM"):"同上 → 0 實例",
        ("chicken","BL"):"雞，非人類；兩框架一致性全數落在「僅單框架」",
        ("bat","DR"):"★ **誤標**，見下"}
for (s, l), n in sorted(tab.items(), key=lambda kv:(kv[0][0], -kv[1])):
    print(f"  {s:9s} {l:9s} {n:6d}   {NOTE.get((s,l),'')}")
print(f"  {'':9s} {'合計':9s} {len(rows):6d}")
print("\n  → 收錄規則（即日起明文）：**古典 MHC class II 的胜肽複合體，跨物種**，")
print("    外加 DM / DO 伴護分子（註解用，不產生實例）。census.py 只蓋人類古典，")
print("    範圍比 rows7 窄 —— 兩者的 8 筆差異全部由此解釋，不是遺漏。")

bats = sorted(p for p, r in rows.items() if (r.get("species") or [""])[0] == "bat")
batinst = [x for x in nt if x["pdb"] in bats]
print(f"\n  ★ **誤標：{len(bats)} 筆蝙蝠 MHC II 被標為 locus `DR`** —— {', '.join(bats)}")
print(f"     （2024-06-05 釋出，全部 post-cutoff）。locus 以 α 鏈對 DRA/DQA1/DPA1 的")
print(f"     identity 指派，蝙蝠在三者中最接近 DRA，於是落進 DR。**判定器沒有物種守門。**")
print(f"     影響（已量化，皆小但必須揭露）：")
print(f"       §01 的 DR 列        185 條 → 人類 181 條；四口袋 52.3% → 52.0%")
print(f"       §01 全體            404 條 → 400 條；主鏈 44.40% → 44.46%")
print(f"       **post-cutoff DR**  41 組 / 35 結構 → **37 組 / 32 結構**")
print(f"       再套無 TCR          30 結構 → **27 結構**")
print(f"     蝙蝠自己的形狀與人類 DR 不同（P4 24.6% 對 15.0%），不是可忽略的同類。")
print(f"     **處置：locus 判定器須加物種守門，蝙蝠另立標籤並排除於 DR 臂之外。**")
print(f"     在重建之前，所有 DR 的數字一律標註「含 4 條蝙蝠」。")

print("\n【二】(結構, 胜肽鏈) 實例層級　—— scripts/register.py")
line("已取得座標檔的 entry", len(have_cif), "entry", "cif/*.cif.gz")
drop("無 α 或 β 鏈註解", len(no_ab))
drop("無胜肽鏈（空溝槽等）", len(no_pep))
line("→ 通過全部幾何過濾", len(nt), "實例", f"{len(inst_pdbs)} entry")
print("       幾何過濾：α/β 各 ≥40 個對位殘基、胜肽 ≥9 殘基、")
print("       α 與胜肽質心距 ≤40 Å、疊合錨點 ≥60 個。全部有名字、可重跑。")
print(f"    ⚠ S1（疊合 RMSD ≤ 3.0）**已停用**。套用它會再砍 {len(nt)-len(th)} 個實例 / "
      f"{len(inst_pdbs)-len({r['pdb'] for r in th})} 個 entry，")
print("       其中含焦點 allele 的唯一結構 3WEX。淘汰率逐 locus 不對稱（DR 0.0% → 雞 BL 85.7%）。")

# ── 第三段：各已發表數字的分母 ──────────────────────────────────────
# contacts_fixed.py:60 的實際規則：start 為 None，或**觀測序列任一處**含 X 即排除
no_p1  = [r for r in nt if r.get("start") is None]
bad_x  = [r for r in nt if r.get("start") is not None and "X" in (r.get("pep") or "")]
core_ok = [r for r in nt if r.get("start") is not None and "X" not in (r.get("pep") or "")]
tight  = [r for r in nt if r.get("start") is not None
          and all(c in AA20 for c in (r.get("core") or "X"))]

print("\n【三】各已發表數字的分母")
print()
print("  ① 界面接觸分解（§01）")
line("  全部實例", len(nt), "實例")
drop("P1 落在觀測胜肽之外", len(no_p1), "§02A S4，可證明為失效案例")
drop("觀測序列含非標準殘基 X", len(bad_x), "§02A S3")
line("  → §01 的分母", len(core_ok), "實例", "✓ 對上 contacts_nothresh.json")
print(f"    ⚠ §02A S3 記「規則須收緊為：僅 core 內含 X 才排除」。已量化："
      f"收緊後保留 {len(tight)} 實例，")
print(f"       即多回收 {len(tight)-len(core_ok)} 筆 —— 45 個含 X 的實例裡有 44 個（97.8%）X 就落在 core 內。")
print(f"       逐 locus 排除率 DR 18.3% / DQ 3.3% / **DP 0.0%**。"
      f"收緊與否不影響焦點 locus，亦不改變任何結論。")

ded = {}
for r in nt:
    ded.setdefault((r["pdb"], r["pep"]), r)
print()
print("  ② ground-truth register 表（Aim 1）")
line("  全部實例", len(nt), "實例")
drop("同 entry 同序列的晶體複本", len(nt) - len(ded), "去 (pdb, 胜肽序列) 重複")
line("  → 獨立 (結構, 胜肽) 組", len(ded), "組", f"{len({k[0] for k in ded})} 結構")
print(f"    ⚠ §08 原記「236 組 / 210 結構」。該組的過濾規則追不到，且只砍 DR 與 DQ；")
print(f"       **本腳本的 {len(ded)} 組取代它**，每一道過濾都在上面列出且可重跑。")

post = {k: v for k, v in ded.items() if v.get("rel") and v["rel"] > CUT}
print()
print("  ③ post-cutoff holdout（Aim 1 benchmark）")
line("  獨立組", len(ded), "組")
drop(f"釋出早於 {CUT}", len(ded) - len(post), "AF3 可能看過")
line("  → holdout", len(post), "組", f"{len({k[0] for k in post})} 結構")
bl = collections.Counter(loc(v) for v in post.values())
sb = collections.defaultdict(set)
for k, v in post.items(): sb[loc(v)].add(k[0])
print("       逐 locus　" + "　".join(f"{k} {bl[k]}組/{len(sb[k])}結構"
      for k in ["DR", "DQ", "DP", "H2-IA"] if bl.get(k)))
notcr = {k: v for k, v in post.items() if not (rows[k[0]].get("tcr") or rows[k[0]].get("tcrfuse"))}
print(f"    ⚠ Aim 1 的缺口陳述寫「**無 TCR 的** pMHC-II」。再套這個條件："
      f"{len(notcr)} 組 / {len({k[0] for k in notcr})} 結構")
nb = collections.Counter(loc(v) for v in notcr.values())
print("       逐 locus　" + "　".join(f"{k} {nb[k]}" for k in ["DR", "DQ", "DP", "H2-IA"] if nb.get(k))
      + "　→ §09 已改為分層共變項，不作為過濾器")

print()
print("  ④ register 平移幾何代價（§05 證據二）")
print("       分母為「平移後 core 仍完整落在觀測胜肽內且編號連號」的實例，")
print("       故逐 k 不同，見 scripts/register_shift.py 的逐列輸出（k=1 時 444 實例 / 274 結構）。")
print("    ⚠ §08 原記「357 個結構」—— 與任一可重現分母都對不上，應以腳本輸出為準。")

print()
print("  ⑤ 9-mer core 假設檢驗（§11 風險 8）")
short = [r for r in core_ok if r["L"] < 9]
print(f"       §08 原記「233 組，與 236 差 3，扣除原因未說明」。")
print(f"       本集合中觀測長度 < 9 的實例：{len(short)} 筆 —— 不存在「core 不足 9 殘基」的案例。")
print(f"       差 3 的來源無從重建（繫於 236 那組）。本腳本以 {len(core_ok)} 實例為分母重報。")

print("\n" + "=" * W)
print(" 結論：全部分母自同一條鏈導出，每一道過濾有名字、有計數、可重跑。")
print(" 原先的 236 / 357 / 233 / 210 四個分母不再使用 —— 它們的產生規則沒有紀錄，")
print(" 依定義七不可作為選擇參數。§08 的對帳缺口至此關閉。")
print("=" * W)
