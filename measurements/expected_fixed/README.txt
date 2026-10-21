v6 採用值的基準（bug 修正後）
================================

本目錄是 scripts/contacts_fixed.py 的預期輸出，對應 proposal_v6.txt §01 的引用值：
主鏈介導接觸 44.3%（中位 43.8%，IQR 39.9–48.0%）、P1/P4/P6/P9 佔全部接觸 29.7%。

與 ../expected/ 的差別
----------------------
../expected/        scripts/contacts.py 的輸出。重現 v5 引用值 47.0% / 23.7%。
                    保留不動，因為「這個包能不能重現 v5 發表的數字」本身是一個
                    需要能被檢驗的問題。
expected_fixed/     scripts/contacts_fixed.py 的輸出。v6 採用值。

兩者的差異來自一個 bug：contacts.py 未對共價融合構築套用胜肽長度截斷，
360 條鏈中有 82 條（23%）把整條 beta 鏈當成胜肽。詳見 ../../../contacts_audit_20261003.md。

可重現性與正確性是兩個不同的主張，本包同時保留兩者的基準，兩層都要通過。
