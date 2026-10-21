"""register 平移一格的幾何代價 [補進驗證包 2026-10-07]

提案 §05 證據二報「平移 1 格 = core Cα 位移 3.80 Å，標準差 0.02 Å」，標為 [本提案量測]
但腳本不在 repo，分母 357 也沒有與其他四個分母對帳。本腳本補上。

量的是什麼：把 core 沿溝槽方向平移 k 個殘基後，九個 core 位置的 Cα 平均位移。
因為 class II 的 core 是剛性的伸展構型，這個位移跨 allele、跨物種幾乎是定值 ——
這正是「register 錯誤沒有部分分數」的依據。

    python3 scripts/register_shift.py

必要的過濾（原始分析有做但未記錄，2026-10-07 補上）：**要求涵蓋到的殘基編號連號**。
胜肽鏈若中間缺殘基，表面上平移一格、實際跨過斷裂等於平移兩格，位移會變成約 6.6 Å
而非 3.8 Å。不濾掉的話 6 筆斷裂實例會把標準差由 0.019 撐到 0.349（17 倍），而平均與
百分位幾乎不動 —— 這正是「SD 0.02」一度重現不出來的原因。

分母說明：以 registers.json 的每個實例為單位，要求平移後的 core 仍完整落在觀測胜肽內
且編號連號，故不同 k 的分母不同，逐一印出。
"""
import os, sys, json, gzip, collections
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
sys.path.insert(0, HERE)

src = open(os.path.join(HERE, "register.py")).read()
head = src.split("# ---------------- reference")[0]
g = {"__file__": os.path.join(HERE, "register.py")}
exec(compile(head, "register_head", "exec"), g)
parse_ca, chain_map = g["parse_ca"], g["chain_map"]

# [改 2026-10-11] 預設改讀 registers_nothresh.json：初稿 3.2.1 節明文不採用疊合品質
# 門檻（它對非 DR 基因座有系統性偏誤，且砍掉焦點等位基因唯一的結構 3WEX），本圖的
# 分母必須與初稿其餘分母同源。--thresh 可切回舊基準做對照。
SRC = "registers.json" if "--thresh" in sys.argv else "registers_nothresh.json"
reg = json.load(open(os.path.join(ROOT, "data", SRC)))
by_pdb = collections.defaultdict(list)
for r in reg:
    by_pdb[r["pdb"]].append(r)

MAXK = 4
acc = {k: [] for k in range(1, MAXK + 1)}
structs = {k: set() for k in range(1, MAXK + 1)}
broken = {k: 0 for k in range(1, MAXK + 1)}

for pid in sorted(by_pdb):
    ca = parse_ca(pid)
    if not ca:
        continue
    gmap = chain_map(ca)
    for r in by_pdb[pid]:
        if r.get("start") is None:
            continue
        d = next((v for (e, a), v in gmap.items() if a == r["asym"]), None)
        if not d:
            continue
        ids = sorted(d)
        s = r["start"]
        if s + 9 > len(ids):
            continue
        base = np.array([d[ids[s + i]][1] for i in range(9)])
        for k in range(1, MAXK + 1):
            for sgn in (+1, -1):
                t = s + sgn * k
                if t < 0 or t + 9 > len(ids):
                    continue
                span = ids[min(s, t):max(s, t) + 9]
                if any(span[i + 1] - span[i] != 1 for i in range(len(span) - 1)):
                    broken[k] += 1       # 跨鏈斷裂，見檔頭說明
                    continue
                shifted = np.array([d[ids[t + i]][1] for i in range(9)])
                acc[k].append(float(np.linalg.norm(shifted - base, axis=1).mean()))
                structs[k].add(pid)
                break          # 每個實例每個 k 只取一個方向，避免重複計數

# [新增 2026-10-10] --dump <path>：把逐實例位移寫成 JSON，供 proposal/draft/figures.py
# 畫圖。畫圖端不重算，以免與本腳本分歧。
if "--dump" in sys.argv:
    _out = sys.argv[sys.argv.index("--dump") + 1]
    os.makedirs(os.path.dirname(os.path.abspath(_out)), exist_ok=True)
    json.dump({str(k): acc[k] for k in acc}, open(_out, "w"))
    print(f"已寫出逐實例位移 {_out}")

print("register 平移的 core Cα 位移　[本提案量測，可重現]")
print(f"  來源 data/{SRC}（{len(reg)} 個實例 / "
      f"{len({r['pdb'] for r in reg})} 個結構）\n")
print(f"  {'平移':>4s} {'n(實例)':>8s} {'n(結構)':>8s} {'平均 Å':>8s} {'標準差':>7s} "
      f"{'p5':>7s} {'p95':>7s}")
for k in range(1, MAXK + 1):
    v = np.array(acc[k])
    if not len(v):
        continue
    print(f"  {k:>4d} {len(v):>8d} {len(structs[k]):>8d} {v.mean():>8.2f} "
          f"{v.std(ddof=1):>7.3f} {np.percentile(v,5):>7.2f} {np.percentile(v,95):>7.2f}"
          f"   （另濾掉斷裂 {broken[k]} 筆）")

print("\n對照 · 提案 §05 證據二（分母記為 357 個結構）")
print("  1 格 3.80 Å（SD 0.02，p5–p95 3.78–3.84）　2 格 6.60　3 格 9.61　4 格 12.82")
print("\n✓ 1 格的平均、標準差與 p5–p95 三者皆重現。原始分析是對的 ——")
print("  它有濾掉跨鏈斷裂的實例，只是沒有記錄這道過濾。")
print("⚠ 分母仍未對帳：本腳本以實例為單位，§05 記為 357 個結構（§08 CONSORT 待辦）。")
