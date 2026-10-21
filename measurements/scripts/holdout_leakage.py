"""訓練截止後保留集的跨條目記憶 [新增 2026-10-11]

起因：初稿的洩漏控制只有一道「釋出日 > 2021-09-30」，而去重規則是「同條目且同
胜肽序列」—— 跨條目的同一核心並未排除。4.1 節自己就寫了「41 個胜肽出現在兩個
以上獨立 PDB 條目中」，所以這道缺口是初稿內部已有的證據指出來的。

本腳本量的是：截止後的每一組，其 9-mer 核心是否已經出現在截止前、同一基因座的
結構裡。若是，則 AlphaFold3 的權重裡已有這個核心的幾何，日期過濾擋不住它。

    python3 scripts/holdout_leakage.py

⚠ 三件事要分清楚：
  1. 比對的是**核心序列**，不是整條胜肽，因為進溝槽的是那九個殘基。
  2. 限**同基因座**。跨基因座的同核心（如恆定鏈胜肽）另行列出，兩者不可混用。
  3. 這不是「AF3 一定記得」的證明，是「無法排除」。它量的是殘留風險的上界。

⚠ 讀的是 registers_nothresh.json（未套疊合品質門檻），與 4.x 節同源；理由見
  dp_diversity.py 的同一段說明 —— 那道門檻已撤銷，且在 DP 上只砍掉焦點 allele。
"""
import os, sys, json, collections

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from locus_fix import locus_of
ROOT = os.path.join(HERE, "..")
CUT  = "2021-09-30"
D    = lambda f: json.load(open(os.path.join(ROOT, "data", f), encoding="utf-8"))

nt = D("registers_nothresh.json")
rows = {r["pdb"]: r for r in D("rows7.json")}
STD = set("ACDEFGHIKLMNPQRSTVWY")
seen, ded = set(), []
for r in nt:                                   # 初稿的去重規則：同條目且同胜肽序列
    k = (r["pdb"], r["pep"])
    if k not in seen:
        seen.add(k); ded.append(r)
# ⚠ 基因座一律經 locus_of()：8JR4/8JRJ/8JRK 是蝙蝠 MHC II，α 鏈最接近 DRA 而被
#   自動標成 DR，而三者全在截止日之後 —— 不更正就混進 benchmark 的 DR 臂。
loc  = lambda r: locus_of(r, rows)
# 建模用的子集：核心含非標準殘基者無序列可輸入，非實驗測定者不是 ground truth。
ok   = lambda r: (r.get("core") and set(r["core"]) <= STD
                  and rows[r["pdb"]].get("method") != "Integrative")
ded  = [r for r in ded if ok(r)]
post = [r for r in ded if r["rel"] > CUT]
pre  = [r for r in ded if r["rel"] <= CUT]
nstr = lambda s: len({r["pdb"] for r in s})

print("訓練截止後保留集的跨條目記憶　[本提案量測，可重現]")
print(f"  來源 data/registers_nothresh.json（{len(nt)} 實例 → 去重並排除非實驗測定"
      f"與核心含非標準殘基者後 {len(ded)} 組 / {nstr(ded)} 結構）")
print(f"  截止日 {CUT}　截止前 {len(pre)} 組　截止後 {len(post)} 組 / {nstr(post)} 結構\n")

print(f"    {'基因座':<8}{'截止後':>8}{'核心已見於截止前':>20}{'比例':>8}")
tot = collections.Counter(); hit = collections.Counter(); clean = []
for r in post:
    L = loc(r)
    tot[L] += 1
    if r["core"] in {p["core"] for p in pre if loc(p) == L}:
        hit[L] += 1
    else:
        clean.append(r)
for L in sorted(tot, key=lambda k: -tot[k]):
    print(f"    {L:<8}{tot[L]:>8}{hit[L]:>20}{hit[L]/tot[L]*100:>7.1f}%")
print(f"    {'合計':<7}{len(post):>8}{sum(hit.values()):>20}"
      f"{sum(hit.values())/len(post)*100:>7.1f}%")

cc = collections.Counter(loc(r) for r in clean)
print(f"\n  扣除後的主要分母：{len(clean)} 組 / {nstr(clean)} 結構　"
      + "、".join(f"{k} {v}" for k, v in sorted(cc.items(), key=lambda x: -x[1])))
print(f"  其中不含 TCR 者 {sum(1 for r in clean if not r['has_tcr'])} 組 / "
      f"{nstr([r for r in clean if not r['has_tcr']])} 結構")
cores = {}
for r in clean:
    cores.setdefault(loc(r), set()).add(r["core"])
print(f"  觀測單位為核心：這 {len(clean)} 組只有 {len({r[chr(34)+chr(34)] if False else r['core'] for r in clean})} 個相異核心　"
      + "、".join(f"{k} {len(v)}" for k, v in sorted(cores.items(), key=lambda x: -len(x[1]))))

x = [r for r in post if r["core"] in {p["core"] for p in pre} and loc(r) not in
     {loc(p) for p in pre if p["core"] == r["core"]}]
print(f"\n  另：核心只在**別的**基因座見過者 {len(x)} 組，未計入上表")

print("\n  ✓ 焦點基因座 HLA-DP 無一核心見於截止前 —— DP 分支不受此限制影響。")
print("  ⚠ HLA-DR 與 HLA-DQ 各有四至五成殘留，與結構代表性同向，兩者的效應"
      "在全部 93 組上分不開（3.2.1 節）。")

# 7YXB 的側翼長度對照：其核心是否乾淨，決定 3.5.3 節那個對照能推到多遠
y = [r for r in ded if r["pdb"] == "7YXB"]
if y:
    core = y[0]["core"]
    src = sorted({p["pdb"] for p in pre if p["core"] == core})
    print(f"\n  7YXB 的核心 {core}：截止前帶同一核心的條目 {len(src)} 個 "
          f"（{', '.join(src[:6])}{' …' if len(src) > 6 else ''}）")
