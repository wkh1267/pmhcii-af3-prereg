"""全 PDB 的 HLA-DP 結構：真實樣貌與多樣性剖析 [新增 2026-10-07]

起因：口試方向確立的第一個前提是「HLA-DP 有 14 個結構」。老師問的是下一個問題 ——
**那 14 個的多樣性夠不夠？胜肽或 allele 的重複性高不高？** 這直接決定夠不夠做驗證。

本腳本不談「有幾個」，談「那 14 個實際上是什麼」。三個層次：

  1. allele 多樣性   以**序列**分群（α1+β1 逐位元組比對），不用 UniProt accession ——
                     §08 陷阱 4：3WEX 與 7T6I/7ZAK/7ZFR 共用 TrEMBL S6B6U4 但 β1 實際不同
  2. 胜肽多樣性      獨特胜肽數、獨特 9-mer core 數、來源蛋白
  3. 實驗獨立性      同一篇論文同時沉積的多個 entry 不是獨立觀測

    python3 scripts/dp_diversity.py

⚠ 本腳本讀的是 registers_nothresh.json（未套用疊合品質門檻）。理由見下方「S1 門檻」段：
  那道門檻是 §02A 已標為「無依據、已證實有偏誤、**已撤銷**」的參數，而它在 DP 上
  唯一砍掉的就是 3WEX —— 焦點 allele。套用它會讓焦點 allele 的結構數變成 0。
"""
import os, json, collections, itertools

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
D = lambda f: json.load(open(os.path.join(ROOT, "data", f)))

rows = {r["pdb"]: r for r in D("rows7.json")}
nt   = D("registers_nothresh.json")
th   = D("registers.json")
DP   = sorted(p for p, r in rows.items() if r.get("loci") == ["DP"])

# 文獻給的 allele 名稱（來自 entry title 與原始論文；序列分群是客觀依據，名稱只是標籤）
# allele 名稱：獨立查核（2026-10-07，RCSB 序列搜尋 + 對 UniProt 參考逐殘基比對）
# 與本腳本的序列分群**獨立得到相同的分組**，兩條路徑交叉確認。
NAMED = {"3LQZ": "DPA1*01:03/DPB1*02:01 (DP2)", "4P4K": "DPA1*01:03/DPB1*02:01 (DP2)",
         "4P4R": "DPA1*01:03/DPB1*02:01 (DP2)", "4P5K": "DPA1*01:03/DPB1*02:01 (DP2)",
         "4P5M": "DPA1*01:03/DPB1*02:01 (DP2)", "4P57": "同上 + βE69K 工程突變",
         "7T2A": "DPA1*01:03/DPB1*04:01 (DP4)", "7T2B": "DPA1*01:03/DPB1*04:01 (DP4)",
         "7T2C": "DPA1*01:03/DPB1*04:01 (DP4)", "7T2D": "DPA1*01:03/DPB1*04:01 (DP4)",
         "7T6I": "DPA1*02:01/DPB1*01:01", "7ZAK": "DPA1*02:01/DPB1*01:01",
         "7ZFR": "DPA1*02:01/DPB1*01:01", "3WEX": "DPA1*02:02/DPB1*05:01 (DP5)"}
# ⚠ 名稱是文獻標籤，分群依據是序列。兩者有一處不一致，值得記下：
#   7T6I 題為「HLA-DP1」、7ZAK/7ZFR 題為「DPA1*02:01-DPB1*01:01」，
#   但三者的 α1 與 β1 **逐殘基完全相同**（0 處差異）—— 是同一個分子的兩種叫法。
#   4P57 題為「MHC TCR peptide complex」但**不含 TCR**（entity 只有 2 個），是沉積時的筆誤。
# 同一篇論文／同一次沉積（以 title 與 release date 判定，見下方列印）
PAPER = {"3LQZ": "Dai 2010 鈹", "4P4K": "Falta 2014 鈹", "4P4R": "Falta 2014 鈹",
         "4P57": "Falta 2014 鈹", "4P5K": "Falta 2014 鈹", "4P5M": "Falta 2014 鈹",
         "3WEX": "DP5 + Cry j 1", "7T2A": "DP4 + Ply", "7T2B": "DP4 + Ply",
         "7T2C": "DP4 + Ply", "7T2D": "DP4 + Ply",
         "7T6I": "Klobuch 2022 反向", "7ZAK": "Racle/DP 反向", "7ZFR": "Racle/DP 反向"}
CUT = "2021-09-30"

# 胜肽的來源蛋白與構築註記。獨立查核（2026-10-07）以 UniProt 全庫逐殘基比對 +
# 原始論文確認。⚠ 標 [改] 者為構築帶有非天然取代，不是野生型序列。
PEPSRC = {
 "RKFHYLPFLPST":    "HLA-DRA 自體肽 pDRα（P01903）",
 "KVTVAFNQF":       "Cry j 1 214–222（日本柳杉過敏原）",
 "QAFWIDLFETI":     "mimotope-2，合成 PS-library 模擬肽（非天然）",
 "NKFDTQLFHTIT":    "pRAS = Rab-9A 30–40（P51151）[改] 第 12 位 T，天然為 G",
 "QAYDGKDYIALK":    "pA28 = HLA class I 重鏈 α2 [改] N 端 Q 非天然側翼",
 "GATGLAWEWWRTVYE": "Pneumolysin 427–441（Q04IN8）[改] 天然為 ECTGLAW…，N 端兩殘基被取代",
 "LPVADAVIHASGKQMWQ": "HCMV pp65 142–158（P06725）",
 "KNLEKYKGKFVREID": "⚠ 來源不明。全人類 UniProt 無 ≥8 殘基完全相符者；論文付費牆",
 "IEFVFKNKAKEL":    "⚠ 來源不明。同上",
}

A1, B1 = 84, 95          # 溝槽結構域（與 register.py 的 A1_MAX / B1_MAX 一致）

# allele 分群不可用原始序列前綴 —— 構築的標籤與接頭長度不同會把同一個 allele 拆開
# （首次跑本腳本時 10 個「獨特配對」裡有一半是這個假象）。改用與 register.py 完全
# 相同的比對：把每條鏈 Needleman-Wunsch 對到 1DLH 參考，取參考座標 0..A1/0..B1 的
# 對位殘基串成「溝槽指紋」。構築差異因此被對掉，只剩真正的序列差異。
import importlib.util
_s = importlib.util.spec_from_file_location("_reg", os.path.join(HERE, "register.py"))
_src = open(os.path.join(HERE, "register.py")).read().split("\n")
_cut = next((i for i, l in enumerate(_src) if l.startswith("out = []") or l.startswith("res = [")), len(_src))
_g = {"__file__": os.path.join(HERE, "register.py"), "__name__": "_reg"}
exec(compile("\n".join(_src[:_cut]), "reg", "exec"), _g)
nw, REF_A_SEQ, REF_B_SEQ = _g["nw"], _g["ref_a_seq"], _g["ref_b_seq"]

# allele 分群踩過三個坑，每一個都會給出錯的多樣性數字：
#   1. 直接切序列前綴          構築的標籤／接頭長度不同 → 同一個 allele 被拆開（得 10）
#   2. 對到 1DLH（HLA-DR1）     參考是 DR，DP 自己的多型位置對不上而被丟掉（得 1）
#   3. 用 register.py 的 nw()   ★ 它只回傳**完全相同**的對位（`if a[i-1]==b[j-1]`），
#                              錯配被靜默丟棄，所以任何用它建的指紋必然塌成一個
# 正解：不用比對器。取 14 條序列的**最長共同子字串**當錨（那段必定在溝槽結構域內
# 且逐一出現），再以錨為原點切固定視窗 —— 構築切點被對齊，多型位置全部保留。

def _lcs_all(seqs):
    """所有序列共有的最長子字串。以最短者的子字串由長到短試。"""
    s0 = min(seqs, key=len)
    for L in range(min(60, len(s0)), 7, -1):
        for i in range(len(s0) - L + 1):
            sub = s0[i:i+L]
            if all(sub in s for s in seqs):
                return sub
    return None

def _dom(p, which, lim):
    e = rows[p][which][0]
    return e[4][max(0, e[3]):]

_aseqs = {p: _dom(p, "alpha", A1) for p in DP}
_bseqs = {p: _dom(p, "beta",  B1) for p in DP}
_aanc = _lcs_all(list(_aseqs.values()))
_banc = _lcs_all(list(_bseqs.values()))

def _window(seq, anchor, back, fwd):
    i = seq.index(anchor)
    lo, hi = i - back, i + len(anchor) + fwd
    if lo < 0 or hi > len(seq): return None
    return seq[lo:hi]

# 視窗必須**只蓋溝槽結構域**。第一版讓它往後延伸到 Ig 結構域與融合標籤，結果
# 7T6I 與 7ZAK 被一個溝槽外的差異拆成兩群（獨立查核把兩者都定為 DPA1*02:01/
# DPB1*01:01，是同一個分子）。改為：以錨回推一個對所有序列一致的起點，再取
# 恰好 A1 / B1 個殘基。
def _start(seqs, anchor):
    return min(s.index(anchor) for s in seqs.values())

_ab, _bb = _start(_aseqs, _aanc), _start(_bseqs, _banc)

def _grv(seq, anchor, back, lim):
    i = seq.index(anchor) - back
    return seq[i:i + lim]

def chains(p):
    return (_grv(_aseqs[p], _aanc, _ab, A1), _grv(_bseqs[p], _banc, _bb, B1))

def peptides(p):
    r, out = rows[p], []
    for e, L, seq, desc in r.get("sep", []) + r.get("ii", []):
        out.append((seq, desc or "", "分離鏈"))
    for f in r.get("fused", []):
        if f["kind"] == "fused_peptide":
            out.append((f.get("pep", ""), f.get("desc", "") or "", "融合"))
    return out

core_of = {}
for r in nt:
    if (r.get("loci") or [""])[0] == "DP" and r.get("core"):
        core_of.setdefault(r["pdb"], set()).add((r["core"], bool(r.get("reverse"))))

bar = "─" * 92
print(f"全 PDB 的 HLA-DP 胜肽複合體：{len(DP)} 個 entry　[本提案量測，可重現]")
print(f"  來源 data/rows7.json（RCSB 序列搜尋，取用 2026-10-07，scripts/census.py）")
print(f"  AF3 訓練截止 {CUT}\n{bar}")

# ---- 1. allele 多樣性：以 α1 / β1 序列分群 --------------------------------
ga, gb, gp = {}, {}, {}
for p in DP:
    a, b = chains(p)
    ka, kb = a, b
    ga.setdefault(ka, []).append(p); gb.setdefault(kb, []).append(p)
    gp.setdefault((ka, kb), []).append(p)

print(f"\n1 · allele 多樣性（以 α1 / β1 **序列**分群，不用 accession）")
print(f"   錨 α「{_aanc[:28]}…」({len(_aanc)}aa)、"
      f"β「{_banc[:28]}…」({len(_banc)}aa)；各取溝槽域 α1 {A1}aa、β1 {B1}aa")
print(f"   獨特 α1 序列 {len(ga)}　獨特 β1 序列 {len(gb)}　**獨特 α/β 配對 {len(gp)}**\n")
print(f"   {'#':>2s} {'結構':28s} {'文獻名稱':24s} {'n':>2s}")
for i, (k, ps) in enumerate(sorted(gp.items(), key=lambda kv: -len(kv[1])), 1):
    nm = {NAMED.get(p, "?") for p in ps}
    print(f"   {i:2d} {', '.join(sorted(ps)):28s} {'/'.join(sorted(nm)):24s} {len(ps):2d}")

reps = {}
for k, ps in sorted(gp.items(), key=lambda kv: -len(kv[1])):
    reps[sorted(ps)[0]] = ps
rk = list(reps)
print(f"\n   逐對差異（溝槽結構域，α1 / β1 各自的不同殘基數）")
print(f"   {'':14s}" + "".join(f"{r:>8s}" for r in rk))
for x in rk:
    line = f"   {x + '(' + str(len(reps[x])) + ')':14s}"
    for y in rk:
        ax, bx = chains(x); ay, by = chains(y)
        da = sum(1 for i, j in zip(ax, ay) if i != j)
        db = sum(1 for i, j in zip(bx, by) if i != j)
        line += f"{str(da) + '/' + str(db):>8s}"
    print(line)
print(f"   （格式 α差/β差。0/0 = 同一個分子）")

# ---- 2. 胜肽多樣性 ---------------------------------------------------------
pep2s = collections.defaultdict(list)
for p in DP:
    for seq, desc, kind in peptides(p):
        pep2s[seq].append(p)
cores = collections.defaultdict(list)
for p, cs in core_of.items():
    for c, rev in cs:
        cores[c].append((p, rev))

print(f"\n{bar}\n2 · 胜肽多樣性")
print(f"   獨特觀測胜肽序列 {len(pep2s)}　獨特 9-mer core {len(cores)}　（對 {len(DP)} 個結構）\n")
print(f"   {'胜肽序列':20s} {'core':11s} {'向':3s} {'n':>2s}  結構 / 來源")
for seq, ps in sorted(pep2s.items(), key=lambda kv: (-len(kv[1]), kv[0])):
    cs = set()
    for p in ps: cs |= core_of.get(p, set())
    c, rev = (sorted(cs)[0] if cs else ("—", False))
    print(f"   {seq:20s} {c:11s} {'反' if rev else '正':3s} {len(ps):2d}  {', '.join(sorted(ps))}")
    print(f"   {'':20s} └ {PEPSRC.get(seq, '（未查）')}")

cplx = {(chains(p), s) for p in DP for s, _, _ in peptides(p)}
print(f"\n   ★ **獨特 (allele 配對, 胜肽) 複合體 = {len(cplx)}** —— 這才是真正的獨立條件數")
nat = {(c, s) for c, s in cplx}
print(f"     14 個 entry 塌成 {len(cplx)} 個複合體；若把 4P57 的工程突變併回 DP2 則為 {len(cplx)-1}")

# ---- 3. 實驗獨立性 ---------------------------------------------------------
pg = collections.defaultdict(list)
for p in DP: pg[PAPER[p]].append(p)
print(f"\n{bar}\n3 · 實驗獨立性（同一次沉積的多個 entry 不是獨立觀測）")
print(f"   {len(DP)} 個 entry 來自 **{len(pg)} 組獨立沉積**\n")
for k, ps in sorted(pg.items(), key=lambda kv: -len(kv[1])):
    print(f"   {k:20s} {len(ps):2d}　{', '.join(sorted(ps))}")

# ---- 4. 逐層收斂 -----------------------------------------------------------
def subset(ps, **kw):
    out = []
    for p in ps:
        r = rows[p]
        if kw.get("post") and not (r.get("rel") and r["rel"] > CUT): continue
        if kw.get("notcr") and (r.get("tcr") or r.get("tcrfuse")): continue
        if kw.get("fwd"):
            cs = core_of.get(p, set())
            if not cs or all(rev for _, rev in cs): continue
        out.append(p)
    return out

def stat(ps):
    pr = {chains(p) for p in ps}
    pe = {s for p in ps for s, _, _ in peptides(p)}
    gpp = {PAPER[p] for p in ps}
    return len(ps), len(pr), len(pe), len(gpp)

print(f"\n{bar}\n4 · 逐層收斂：每加一個條件還剩多少\n")
print(f"   {'條件':34s} {'結構':>4s} {'allele配對':>10s} {'獨特胜肽':>9s} {'獨立沉積':>9s}")
for name, kw in [("全部", {}), ("＋ post-cutoff（AF3 未見）", dict(post=1)),
                 ("＋ 無 TCR", dict(post=1, notcr=1)),
                 ("＋ 僅正向結合", dict(post=1, notcr=1, fwd=1))]:
    s = subset(DP, **kw)
    n, a, pp, g = stat(s)
    print(f"   {name:34s} {n:4d} {a:10d} {pp:9d} {g:9d}   {', '.join(sorted(s))}")

# ---- 5. S1 門檻 ------------------------------------------------------------
kept = {r["pdb"] for r in th if (r.get("loci") or [""])[0] == "DP"}
lost = sorted(set(DP) - kept - {p for p in DP if p not in {r["pdb"] for r in nt}})
print(f"\n{bar}\n5 · ⚠ S1「疊合品質 RMSD ≤ 3.0」門檻的實際後果")
print(f"   §02A 已把 S1 標為「無依據、已證實有偏誤、**已撤銷**」，但 register.py 仍在套用。")
print(f"   它在 DP 上唯一砍掉的結構是 **{', '.join(lost) if lost else '（無）'}** —— 也就是焦點 allele HLA-DP5。")
print(f"   套用後焦點 allele 的 ground-truth register 數 = 0。")
rm = collections.defaultdict(list)
for r in nt:
    if r.get("sup_rmsd") is not None:
        rm[(r.get("loci") or ["?"])[0]].append(r["sup_rmsd"])
print(f"\n   該門檻按 locus 的淘汰率（它量的是「離 1DLH／HLA-DR1 參考框架多遠」）：")
print(f"   {'locus':8s} {'實例':>5s} {'被砍':>5s} {'淘汰率':>7s} {'RMSD 中位':>10s}")
for k in ["DR", "DQ", "DP", "H2-IA", "H2-IE", "BL"]:
    v = sorted(rm.get(k, []))
    if not v: continue
    bad = sum(1 for x in v if x > 3.0)
    print(f"   {k:8s} {len(v):5d} {bad:5d} {100*bad/len(v):6.1f}% {v[len(v)//2]:10.2f}")
print(f"\n   → 這不是品質過濾，是**與參考框架的相似度過濾**。DR 0%、DP 12.9%、雞 BL 85.7%。")
print(f"     它系統性地懲罰非 DR 的 locus，而本論文的焦點正是非 DR。")
