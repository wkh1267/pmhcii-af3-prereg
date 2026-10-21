#!/bin/sh
# pMHC-II 界面接觸統計的重現驗證。只需 python3 與 numpy。
#
# 三層，全部必須通過：
#   第 1 層  從 contacts.json 重算彙總（不需結構檔）
#   第 2 層  從 291 個 mmCIF 重新推導 → 重現 v5 引用值 47.0% / 23.7%
#   第 3 層  同上但套用 fused 截斷修正 → 重現 v6 採用值 44.3% / 29.7%
#
# 第 2 與第 3 層測的是兩件不同的事：第 2 層問「這個包能不能重現當初發表的數字」，
# 第 3 層問「修正後的數字是否正確」。可重現性與正確性是兩個獨立的主張。
set -e
cd "$(dirname "$0")/scripts"

fail() { echo; echo "!! $1"; exit 1; }

echo "=============================================================="
echo " 第 1 層：從 contacts.json 重算彙總（不需結構檔，數秒）"
echo "=============================================================="
python3 aggregate.py /tmp/aggregate_rerun.json
echo
python3 - <<'PY'
import json, sys
a = json.load(open("../expected/aggregate.json", encoding="utf-8"))
b = json.load(open("/tmp/aggregate_rerun.json", encoding="utf-8"))
print("與 expected/aggregate.json 一致：", a == b)
sys.exit(0 if a == b else 1)
PY

echo
echo "=============================================================="
echo " 第 2 層：從 291 個 mmCIF 重新推導（原腳本，重現 v5 引用值）"
echo "=============================================================="
cp ../data/contacts.json /tmp/contacts_before.json
python3 contacts.py | tee /tmp/contacts_stdout.txt
echo
diff -q /tmp/contacts_stdout.txt ../expected/contacts_stdout.txt >/dev/null \
  || fail "第 2 層 stdout 與 expected/contacts_stdout.txt 不符"
echo "stdout 與 expected/contacts_stdout.txt 完全相同"
python3 - <<'PY'
import json, sys
a = json.load(open("/tmp/contacts_before.json", encoding="utf-8"))
b = json.load(open("../data/contacts.json", encoding="utf-8"))
print("重新推導的 contacts.json 與原檔一致：", a == b, f"({len(b)} 筆)")
sys.exit(0 if a == b else 1)
PY

echo
echo "=============================================================="
echo " 第 3 層：套用 fused 截斷修正（重現 v6 採用值）"
echo "=============================================================="
python3 contacts_fixed.py | tee /tmp/contacts_fixed_stdout.txt
echo
diff -q /tmp/contacts_fixed_stdout.txt ../expected_fixed/contacts_fixed_stdout.txt >/dev/null \
  || fail "第 3 層 stdout 與 expected_fixed/contacts_fixed_stdout.txt 不符"
echo "stdout 與 expected_fixed/contacts_fixed_stdout.txt 完全相同"
python3 - <<'PY'
import json, sys
a = json.load(open("../expected_fixed/contacts_fixed.json", encoding="utf-8"))
b = json.load(open("../data/contacts_fixed.json", encoding="utf-8"))
print("重新推導的 contacts_fixed.json 與基準一致：", a == b, f"({len(b)} 筆)")
sys.exit(0 if a == b else 1)
PY

echo
echo "=============================================================="
echo " 第 4 層：逐 locus 分層（不需結構檔，數秒）"
echo "=============================================================="
python3 percore_by_locus.py | tee /tmp/percore_by_locus_stdout.txt
echo
diff -q /tmp/percore_by_locus_stdout.txt ../expected_fixed/percore_by_locus_stdout.txt >/dev/null \
  || fail "第 4 層 stdout 與 expected_fixed/percore_by_locus_stdout.txt 不符"
echo "stdout 與 expected_fixed/percore_by_locus_stdout.txt 完全相同"

echo
echo "=============================================================="
echo " 第 5 層：撤銷疊合門檻後的全集（定義八；含焦點 allele 的 3WEX）"
echo "=============================================================="
SUP_RMSD_MAX=inf python3 register.py /tmp/registers_nothresh.json >/dev/null
python3 - <<'PY'
import json, sys
a = json.load(open("/tmp/registers_nothresh.json", encoding="utf-8"))
b = json.load(open("../data/registers_nothresh.json", encoding="utf-8"))
print("重新推導的 registers_nothresh.json 與基準一致：", a == b, f"({len(b)} 筆)")
sys.exit(0 if a == b else 1)
PY
REGISTERS=../data/registers_nothresh.json CONTACTS_OUT=/tmp/contacts_nothresh.json \
  python3 contacts_fixed.py | tee /tmp/contacts_nothresh_stdout.txt
echo
diff -q /tmp/contacts_nothresh_stdout.txt ../expected_nothresh/contacts_nothresh_stdout.txt >/dev/null \
  || fail "第 5 層 stdout 與 expected_nothresh/contacts_nothresh_stdout.txt 不符"
echo "stdout 與 expected_nothresh/contacts_nothresh_stdout.txt 完全相同"
CONTACTS=../data/contacts_nothresh.json python3 percore_by_locus.py \
  > /tmp/percore_nothresh_stdout.txt
diff -q /tmp/percore_nothresh_stdout.txt ../expected_nothresh/percore_by_locus_nothresh_stdout.txt >/dev/null \
  || fail "第 5 層逐 locus stdout 與基準不符"
echo "逐 locus stdout 與 expected_nothresh/percore_by_locus_nothresh_stdout.txt 完全相同"

echo
echo "五層全部通過。"
echo "  v5 引用值（原腳本）  ：主鏈 47.0% / anchor 23.7%   n=360"
echo "  v6 採用值（已修正）  ：主鏈 44.3% / anchor 29.7%   n=360"
echo "  撤銷疊合門檻後的全集 ：主鏈 44.4% / anchor 29.8%   n=404（含 3WEX）"
echo "  逐 locus（全集）：DP 的 P6 達 26.7%（DR 10.3%），P7 11.6% 排第三"
