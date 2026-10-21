"""人類古典 class II 結構普查 —— 可重現版本 [建立 2026-10-07]

為什麼需要這支：提案 §03 的 locus 表（DR 145 / DQ 59 / DP 14 = 218）與 §04 的
「236 組」都沒有記錄 RCSB 查詢條件與取用日期，無法重現，而且與本驗證包自己的資料
對不上（rows7.json 給 DR 149）。本腳本把查詢、分類規則、取用日期全部寫死並存檔，
讓那張表之後一行指令就能重建。

    python3 scripts/census.py              # 需網路，約 2–3 分鐘
    python3 scripts/census.py --offline    # 重讀既有輸出，只重算統計

輸出 data/census_<取用日期>.json，內含查詢 JSON、取用日期、PDB ID 全清單、
逐筆 locus 判定與三個 identity 分數 —— 每一個數字都可以逐筆回查。

── 查詢定義 ──────────────────────────────────────────────────────────────
  以 UniProt DRA（P01903）的 α1+α2 做 RCSB **序列搜尋**，identity >= 0.55，
  experimental only。class II 的 α 鏈跨 locus 仍有約 55–60% identity，故單一條
  參考序列即可撈到 DR / DQ / DP 全部。

  ⚠ 不用 Pfam，也不用 UniProt accession，兩者都試過且都不可靠：
    · Pfam：3WEX 與 7T6I 的 entity **完全沒有 Pfam 註釋**，PF00993 查詢漏掉 33 筆，
      其中包含本研究焦點 allele 的唯一結構 3WEX。
    · accession：3WEX 與 7T6I/7ZAK/7ZFR 的 β 鏈共用 TrEMBL S6B6U4 但 β1 序列不同。
      以 accession 推 allele 會得到「DPB1*05:01 有 4 個結構」而非正確的 1 個。
    · taxonomy 9606：只要 entry 內**任一** entity 是人類就會命中，小鼠 MHC 配人類
      胜肽或 TCR 的結構會混進來（1JL4、1K2D、1MUJ 等）。

── locus 判定 ────────────────────────────────────────────────────────────
  取每條長度 150–300 的 polymer entity，與 UniProt 的 α 鏈正規序列比對，
  取 identity 最高的 locus。

      DRA   P01903      DQA1  P01909      DPA1  P20036

  不用 RCSB 的 entity 描述字串 —— 舊條目只寫「HLA CLASS II HISTOCOMPATIBILITY
  ANTIGEN」或「HLA-DR2」，不含 locus，201 筆中有 37 筆無法由描述判定。
  不用 accession，理由同上。

  比對為無 gap 的偏移掃描（class II α1+α2 無插入缺失，足夠）。判定同時輸出
  勝出者與次高者的差距；差距 < MARGIN 者標為低信心，列出供人工覆核。

── 已知界線 ──────────────────────────────────────────────────────────────
  · 非古典分子（DM、DO）會被 α 鏈比對歸到某個古典 locus。以 ENTITY_EXCLUDE
    的描述關鍵字排除，並單獨計數。
  · 本查詢回傳的是 **entry**，不是 (結構, peptide) 組，也不是胜肽鏈。
    提案的五個分母（360 / 357 / 236 / 233 / 210）單位各不相同，見 §08 對帳表。
"""
import json, os, sys, time, urllib.request, datetime, collections

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL = "https://data.rcsb.org/graphql"
UNIPROT = "https://rest.uniprot.org/uniprotkb/{}.fasta"

SEQ_IDENT = 0.55        # 序列搜尋的 identity 下限（跨 locus 的 α 鏈仍在此之上）
CALL_IDENT = 0.60       # identity 下限。不可設到 0.80 —— DQA1 比 DRA 多型，跨 allele
                        # 的 DQ 結構對 P01909 只有約 0.63–0.68（1JK8 = 0.683），0.80 會
                        # 把 26 筆 DQ 全部誤殺。小鼠 H-2 的判別靠 margin 不靠絕對值。
MARGIN_NOTE = ("小鼠 H-2 對三個人類參考序列幾乎等距（margin ~0.001），"
               "故 margin 是物種與 locus 的共同判準")


def build_query(seq):
    return {"query": {"type": "terminal", "service": "sequence", "parameters": {
                "evalue_cutoff": 1, "identity_cutoff": SEQ_IDENT,
                "sequence_type": "protein", "value": seq}},
            "return_type": "entry",
            "request_options": {"paginate": {"start": 0, "rows": 10000},
                                "results_content_type": ["experimental"],
                                "scoring_strategy": "sequence"}}

REF_ACC = {"DR": "P01903", "DQ": "P01909", "DP": "P20036"}
ENTITY_EXCLUDE = ("DM alpha", "DM beta", "DO alpha", "DO beta",
                  "DM ALPHA", "DM BETA", "DO ALPHA", "DO BETA")
MARGIN = 0.08          # 勝出者與次高者的 identity 差距下限
ALPHA_LEN = (150, 300)  # 候選 α 鏈的長度範圍
REF_WINDOW = (25, 205)  # UniProt 序列中 α1+α2 的大致範圍（跳過訊號肽）


def post(url, payload):
    req = urllib.request.Request(url, data=json.dumps(payload).encode(),
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as fh:
        return json.load(fh)


def ident(a, b):
    """無 gap 的偏移掃描 identity，取最佳偏移；分母為實際重疊長度。"""
    best = 0.0
    for off in range(-40, 41):
        m = n = 0
        for i, ch in enumerate(a):
            j = i + off
            if 0 <= j < len(b):
                n += 1
                m += (ch == b[j])
        if n >= 70:
            best = max(best, m / n)
    return best


def fetch_entities(ids):
    out, B = {}, 40
    q = """query($ids:[String!]!){ entries(entry_ids:$ids){ rcsb_id
             polymer_entities { entity_poly{pdbx_seq_one_letter_code_can}
               rcsb_polymer_entity{pdbx_description} } } }"""
    for i in range(0, len(ids), B):
        d = post(GRAPHQL, {"query": q, "variables": {"ids": ids[i:i + B]}})
        for e in d["data"]["entries"]:
            out[e["rcsb_id"]] = [
                (((pe.get("entity_poly") or {}).get("pdbx_seq_one_letter_code_can") or "").replace("\n", ""),
                 ((pe.get("rcsb_polymer_entity") or {}).get("pdbx_description") or ""))
                for pe in (e.get("polymer_entities") or [])]
        sys.stderr.write(f"\r  取序列 {min(i + B, len(ids))}/{len(ids)}")
        time.sleep(0.2)
    sys.stderr.write("\n")
    return out


def main():
    today = datetime.date.today().isoformat()
    path = os.path.join(DATA, f"census_{today}.json")

    if "--offline" in sys.argv or os.path.exists(path):
        if not os.path.exists(path):
            sys.exit(f"--offline 但找不到 {path}")
        rec = json.load(open(path))
        print(f"讀取既有輸出 {os.path.basename(path)}（取用日期 {rec['accessed']}）")
    else:
        print("取 UniProt 參考 α 鏈 …")
        refs = {}
        for loc, acc in REF_ACC.items():
            t = urllib.request.urlopen(UNIPROT.format(acc), timeout=60).read().decode()
            refs[loc] = "".join(t.split("\n")[1:])
            print(f"  {loc}  {acc}  {len(refs[loc])} aa")

        print("查詢 RCSB（序列搜尋，DRA α1+α2 參考）…")
        query = build_query(refs["DR"][REF_WINDOW[0]:REF_WINDOW[1]])
        res = post(SEARCH, query)
        ids = sorted(x["identifier"] for x in res.get("result_set", []))
        print(f"  total_count = {res.get('total_count')}，取回 {len(ids)}")
        ents = fetch_entities(ids)

        rows = []
        for pid in ids:
            descs = [d for _, d in ents.get(pid, [])]
            if any(k in d for d in descs for k in ENTITY_EXCLUDE):
                rows.append({"pdb": pid, "locus": "NONCLASSICAL",
                             "scores": {}, "margin": None})
                continue
            sc = {"DR": 0.0, "DQ": 0.0, "DP": 0.0}
            for s, _ in ents.get(pid, []):
                if not (ALPHA_LEN[0] <= len(s) <= ALPHA_LEN[1]):
                    continue
                for loc, r in refs.items():
                    sc[loc] = max(sc[loc], ident(s[:180], r[REF_WINDOW[0]:REF_WINDOW[1]]))
            order = sorted(sc, key=sc.get, reverse=True)
            top, second = sc[order[0]], sc[order[1]]
            if top < CALL_IDENT:
                call, reason = None, f"identity {top:.3f} < {CALL_IDENT}（非人類古典）"
            elif top - second < MARGIN:
                call, reason = None, f"margin {top - second:.3f} < {MARGIN}"
            else:
                call, reason = order[0], ""
            rows.append({"pdb": pid, "locus": call, "excluded_reason": reason,
                         "scores": {k: round(v, 3) for k, v in sc.items()},
                         "margin": round(top - second, 3)})
        rec = {"accessed": today, "search_url": SEARCH, "query": query,
               "ref_accessions": REF_ACC, "margin_threshold": MARGIN,
               "n_entries_returned": len(ids), "ids": ids, "rows": rows}
        json.dump(rec, open(path, "w"), indent=1)
        print(f"  已寫出 data/{os.path.basename(path)}")

    rows = rec["rows"]
    c = collections.Counter(r["locus"] for r in rows)
    tot = sum(c[k] for k in ("DR", "DQ", "DP"))
    print(f"\n人類古典 class II 結構普查")
    print(f"  查詢      RCSB 序列搜尋，DRA(P01903) α1+α2，identity >= 0.55，experimental")
    print(f"  取用日期  {rec['accessed']}")
    print(f"  判定門檻  identity >= {CALL_IDENT} 且 margin >= {MARGIN}")
    print(f"  回傳 entry {rec['n_entries_returned']}　→　排除非古典 DM/DO {c['NONCLASSICAL']} 筆"
          f"　→　人類古典 {tot} 筆\n")
    print(f"  {'locus':8s} {'entry':>6s} {'佔古典':>8s}")
    for k in ("DR", "DQ", "DP"):
        print(f"  HLA-{k:4s} {c[k]:6d} {c[k]/tot*100:7.1f}%")

    unk = [r for r in rows if r["locus"] is None and r.get("excluded_reason")]
    byr = collections.Counter(r["excluded_reason"].split("（")[0].split("<")[0].strip()
                              for r in unk)
    print(f"\n  排除 {len(unk)} 筆，原因分佈：")
    for k, v in byr.most_common():
        print(f"    {k:28s} {v}")

    print(f"\n對照")
    print(f"  提案 §03 locus 表    DR 145 / DQ  59 / DP 14 = 218")
    print(f"  rows7.json          DR 149 / DQ  59 / DP 14 = 222")
    print(f"  本次普查            DR {c['DR']:3d} / DQ {c['DQ']:3d} / DP {c['DP']:2d} = {tot}")


if __name__ == "__main__":
    main()
