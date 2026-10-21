#!/usr/bin/env python3
"""把四個分析集合導成 CSV。[2026-10-15 隨交接包新增，不在 repo 內]

為什麼要有這一支：repo 裡四個集合（310 / 187 / 81 / 49 組）只有統計數字，沒有任何
PDB ID 清單檔。holdout_leakage.py 會把它們算出來，但只印逐基因座計數，一個 ID 都不印。
到了新機器要建 AlphaFold3 的輸入，第一件事就是拿到這些清單。

本檔不重新推導，只是把 holdout_leakage.py 既有的過濾鏈接上輸出，再用初稿已載的
計數當 assert。任一條對不上就整支失敗，不會寫出半套檔案。

    python3 -I tools/export_sets.py [驗證包路徑] [輸出目錄]

預設路徑：thesis/proposal/verify_contacts_20261003 與 tools/out/
"""
import os, sys, json, csv, collections

PKG = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "proposal", "verify_contacts_20261003")
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "out")
sys.path.insert(0, os.path.join(PKG, "scripts"))
from locus_fix import locus_of                                    # noqa: E402

CUT = "2021-09-30"
D = lambda f: json.load(open(os.path.join(PKG, "data", f), encoding="utf-8"))
STD = set("ACDEFGHIKLMNPQRSTVWY")

nt = D("registers_nothresh.json")
rows = {r["pdb"]: r for r in D("rows7.json")}

# ── 與 holdout_leakage.py 逐字相同的三道過濾 ───────────────────────────
seen, dedup = set(), []
for r in nt:                                   # 去重：同條目且同胜肽序列
    k = (r["pdb"], r["pep"])
    if k not in seen:
        seen.add(k); dedup.append(r)
loc = lambda r: locus_of(r, rows)              # 8JR4/8JRJ/8JRK 的蝙蝠更正在此
ok = lambda r: (r.get("core") and set(r["core"]) <= STD
                and rows[r["pdb"]].get("method") != "Integrative")
ded = [r for r in dedup if ok(r)]              # 可建模者
pre = [r for r in ded if r["rel"] <= CUT]
post = [r for r in ded if r["rel"] > CUT]
pre_cores_by_locus = {}
for p in pre:
    pre_cores_by_locus.setdefault(loc(p), set()).add(p["core"])
clean = [r for r in post
         if r["core"] not in pre_cores_by_locus.get(loc(r), set())]

nstr = lambda s: len({r["pdb"] for r in s})
bylocus = lambda s: dict(collections.Counter(loc(r) for r in s))

# ── assert：全部取自初稿，對不上就不要信這些 CSV ────────────────────────
CHECKS = [
    ("去重後 310 組",        len(dedup),                310),
    ("去重後 278 結構",      nstr(dedup),               278),
    ("可建模 268 組",        len(ded),                  268),
    ("可建模 240 結構",      nstr(ded),                 240),
    ("截止前 187 組",        len(pre),                  187),
    ("截止後 81 組",         len(post),                 81),
    ("截止後 68 結構",       nstr(post),                68),
    ("主要分母 49 組",       len(clean),                49),
    ("主要分母 42 結構",     nstr(clean),               42),
    ("相異 core 35 個",      len({r["core"] for r in clean}), 35),
    ("無 TCR 子集 35 組",    sum(1 for r in clean if not r["has_tcr"]), 35),
    ("81 組逐基因座",  bylocus(post),
     {"DR": 28, "DQ": 27, "H2-IA": 10, "DP": 8, "BL": 4, "BAT": 4}),
    ("49 組逐基因座",  bylocus(clean),
     {"DR": 15, "DQ": 11, "DP": 8, "H2-IA": 7, "BL": 4, "BAT": 4}),
]
bad = [(n, g, w) for n, g, w in CHECKS if g != w]
for n, g, w in CHECKS:
    print(f"  {'✓' if (n, g, w) not in bad else '✗'} {n:<18} 得 {g}"
          + ("" if (n, g, w) not in bad else f"　初稿為 {w}"))
if bad:
    sys.exit("\n✗ 與初稿不符，不寫出 CSV。先查過濾鏈再說。")

# ── 寫出 ──────────────────────────────────────────────────────────────
COLS = ["pdb", "asym", "kind", "locus", "L", "pep", "core", "p1_index",
        "start", "reverse", "rel", "res", "has_tcr", "method", "title"]
os.makedirs(OUT, exist_ok=True)
for name, rs in [("set_310_groundtruth", dedup), ("set_187_precutoff", pre),
                 ("set_81_postcutoff", post), ("set_49_main", clean)]:
    path = os.path.join(OUT, name + ".csv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS)
        w.writeheader()
        for r in sorted(rs, key=lambda r: (r["pdb"], r["asym"])):
            row = rows[r["pdb"]]
            w.writerow({**{c: r.get(c) for c in COLS},
                        "locus": loc(r), "method": row.get("method"),
                        "title": row.get("title")})
    print(f"  寫出 {os.path.relpath(path)}　{len(rs)} 列")

print("\n⚠ 這四個檔是 2026-10-15 由本腳本導出的，不是原本就在 repo 裡的資料。"
      "\n  過濾鏈抄自 scripts/holdout_leakage.py，計數與初稿 13 項逐一相符（上方）。"
      "\n  pep 欄是**觀測序列**（建入電子密度者），不是沉積 entity 序列；"
      "\n  kind 為 fused 者胜肽與 β 鏈共價相連，要套 contacts_fixed.py 的 plen 截斷規則。")
