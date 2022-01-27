#!/bin/bash
xflgo 6
curl -v -X POST https://api.line.me/v2/bot/message/broadcast \
-H 'Content-Type: application/json' \
-H 'Authorization: Bearer SKtZA/K9ayC0aOlUffNJeLkY3hzE3zKLpZLo6fTo8NeyQLaFjycMF8O/GecWZWkTxJZAr3hGa6JcLDXKRmeLhl1NNaomjwySu8RDzK+8PkLczy9xej/01o7D6lSODWCheVb5q9/yTIcJzp3K9kwJMgdB04t89/1O/w1cDnyilFU=' \
-d '{
    "messages":[
        {
            "type":"text",
            "text":"計算終了通知\nフォルダ: ~/as\nファイル: fl6"
        }
    ]
}'
    