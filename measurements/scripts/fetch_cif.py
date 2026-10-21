#!/usr/bin/env python3
"""取得或核對 cif/ 下的 291 個 mmCIF。[新增 2026-10-13]

  python3 scripts/fetch_cif.py --check     不連網，核對本機檔案與 rows7 的清單是否一致
  python3 scripts/fetch_cif.py --fetch     下載缺少者（RCSB，每檔之間間隔 0.2 秒）
  python3 scripts/fetch_cif.py --verify N  連網抽驗 N 個，比對本機與遠端的內容雜湊

為什麼要有這一支：cif/ 是全部分析的源頭（291 個條目、約 86 MB），但 2026-10-13 之前
repo 裡**沒有任何腳本說明它從哪裡來**，README 只寫「291 個 class II 結構的 mmCIF」。
清單的唯一真實來源是 data/rows7.json 的 pdb 欄；本檔以它為準重建或核對。

⚠ RCSB 的 mmCIF 會隨條目改版而變動（revision history 增列、座標修訂）。--verify 不符
  不必然代表本機檔有誤，可能是上游改版；遇到時記下日期與差異，不要逕自覆蓋黃金檔。
"""
import gzip, hashlib, json, os, random, sys, time, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
CIF  = os.path.join(ROOT, "cif")
URL  = "https://files.rcsb.org/download/{}.cif.gz"


def wanted():
    rows = json.load(open(os.path.join(ROOT, "data", "rows7.json"), encoding="utf-8"))
    return sorted({r["pdb"] for r in rows})


def local():
    return sorted(f[:-7] for f in os.listdir(CIF) if f.endswith(".cif.gz"))


def digest(path):
    """對解壓後的內容取雜湊：gzip 標頭含時間戳，壓縮檔本身的雜湊不穩定。"""
    with gzip.open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def main():
    need, have = wanted(), local()
    miss, extra = sorted(set(need) - set(have)), sorted(set(have) - set(need))
    print(f"  rows7 所列 {len(need)} 個條目；cif/ 現有 {len(have)} 個")
    if miss:  print(f"  ✗ 缺少 {len(miss)}：{'、'.join(miss[:8])}{' …' if len(miss) > 8 else ''}")
    if extra: print(f"  ⚠ 多出 {len(extra)}：{'、'.join(extra[:8])}{' …' if len(extra) > 8 else ''}")
    if not miss and not extra: print("  ✓ 清單與 rows7 逐一相同")

    if "--fetch" in sys.argv and miss:
        for i, pid in enumerate(miss, 1):
            dst = os.path.join(CIF, pid + ".cif.gz")
            with urllib.request.urlopen(URL.format(pid), timeout=120) as fh:
                open(dst, "wb").write(fh.read())
            print(f"  ↓ {i}/{len(miss)} {pid}")
            time.sleep(0.2)

    if "--verify" in sys.argv:
        n = int(sys.argv[sys.argv.index("--verify") + 1])
        # 固定亂數種子，抽樣可重現
        pick = random.Random(0).sample(have, min(n, len(have)))
        bad = []
        for pid in pick:
            with urllib.request.urlopen(URL.format(pid), timeout=120) as fh:
                open("/tmp/_fetch_cif.gz", "wb").write(fh.read())
            same = digest(os.path.join(CIF, pid + ".cif.gz")) == digest("/tmp/_fetch_cif.gz")
            print(f"  {'✓' if same else '✗'} {pid}")
            if not same: bad.append(pid)
            time.sleep(0.2)
        print(f"  抽驗 {len(pick)} 個，{len(pick) - len(bad)} 個與 RCSB 現行版本相同")
        return 1 if bad else 0
    return 1 if (miss or extra) else 0


if __name__ == "__main__":
    sys.exit(main())
