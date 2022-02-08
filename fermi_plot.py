from subprocess import call, STDOUT, Popen


def main():
    num = prinput("band番号: ")
    call(["./upload.sh", "wkfermi/wkplot/format.dat"])
    call(["./upload.sh", "wkfermi/wkplot/fort.1"])
    with open("notice.sh", "w") as f:
        f.write(line_notice_sh(num))
    call(["./upload.sh", "notice.sh"])
    Popen(
        f"ssh bern 'ssh bern6 \"cd ~/as && bash notice.sh\"'", shell=True)
    return


def prinput(s: str) -> str:
    print(s, end="")
    return input().strip()


def line_notice_sh(num: str) -> str:
    return f'''#!/bin/bash
cd ~/as/wkfermi/wkplot && ./frmpltmn && pig2ps
curl -v -X POST https://api.line.me/v2/bot/message/broadcast \\
-H 'Content-Type: application/json' \\
-H 'Authorization: Bearer SKtZA/K9ayC0aOlUffNJeLkY3hzE3zKLpZLo6fTo8NeyQLaFjycMF8O/GecWZWkTxJZAr3hGa6JcLDXKRmeLhl1NNaomjwySu8RDzK+8PkLczy9xej/01o7D6lSODWCheVb5q9/yTIcJzp3K9kwJMgdB04t89/1O/w1cDnyilFU=' \\
-d '{{
    "messages":[
        {{
            "type":"text",
            "text":"フェルミ面描画完了\\n以下を実行して確認\\n\\n./open.sh {num}"
        }}
    ]
}}'
'''


if __name__ == "__main__":
    main()
