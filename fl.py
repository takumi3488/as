from subprocess import call, STDOUT, check_output, Popen
from glob import glob
from re import search
from sys import argv
from os import getcwd


def main():
    # 計算結果をダウンロード,見つからなければ最終出力結果を取得
    dat_file = max_glob('fl*.dat')
    out_file = f"{dat_file[0:-4]}.out"
    try:
        check_output(["./download.sh", out_file], stderr=STDOUT)
    except:
        out_file = max_glob('fl*.out')
    # 計算を回す回数を取得
    flg = False
    iter_times = -1
    with open(dat_file) as f:
        for line in f:
            l = line.strip()
            if flg:
                if "-" in l:
                    break
                else:
                    iter_times += 1
            if "ITERATION DATA" in l:
                flg = True

    # 計算結果を出力
    ignore = True
    orbit_count = 0
    left = 0
    flg = False
    with open(out_file) as f:
        for line in f:
            l = line.strip()
            if search(rf"ITER=\s*{iter_times}", l):
                ignore = False
            if ignore:
                continue
            if "IND0=" in l and 0 <= orbit_count <= 3:
                orbit_count += 1
                if orbit_count <= 3:
                    print(f"\n{l}")
            if 0 < orbit_count <= 3 and l and l[0] in ["-", "0", "1"]:
                print(l)
            if left:
                if search(r"D-0[5-9]", l):
                    flg = True
                print(l)
                left -= 1
            if l == "== R.M.S. OF DIFFERNCE BETWEEN OLD AND NEW ==":
                left = 3
                print("")


    # 新しく設定ファイルを作成して開く(もしくは終了する)
    current_num = find_num(out_file)
    next_num = current_num + 1
    current = f"fl{current_num}.dat"
    _next = f"fl{next_num}.dat"
    print(f"\n次の設定ファイル({_next})を開きますか？(Yn): ", end="")
    if input().strip().lower() == "n":
        return
    call(["cp", current, _next])
    call(["code", _next])

    # 設定ファイルを送信
    print("入力が完了したらEnter: ", end="")
    input()
    call(["./upload.sh", _next])
    print("計算を開始")

    # 計算を開始,引数がbであれば計算完了まで待たずに終了
    dir = getcwd().split('/')[-1]
    print("\n1(default): このまま待つ\n2: 終了次第LINEに通知\n> ", end="")
    if input().strip().lower() == "2":
        with open("notice.sh", "w") as f:
            f.write(line_notice_sh(next_num))
        call(["./upload.sh", "notice.sh"])
        Popen(
            f"ssh bern 'ssh bern6 \"cd ~/{dir}/wk && mv ../notice.sh . && bash notice.sh\"'", shell=True)
        return
    call(f"ssh bern 'ssh bern6 \"cd ~/{dir}/wk && xflgo {next_num}\"'", shell=True)
    print("計算完了")

    # 再度最初から処理を始める
    print("\n続けますか？(Yn): ", end="")
    if input().strip().lower() == "n":
        return
    main()


def max_glob(s: str):
    l = glob(s)
    res = l[-1]
    for f in l:
        if find_num(res) < find_num(f):
            res = f
    return res


def find_num(txt: str) -> int:
    num_match = search(r"\d+", txt)
    if not num_match:
        raise ValueError("数値が見つかりませんでした。")
    num = int(num_match.group())
    return num


def line_notice_sh(num: int) -> str:
    dir = getcwd().split('/')[-1]
    return f'''#!/bin/bash
xflgo {num}
curl -v -X POST https://api.line.me/v2/bot/message/broadcast \\
-H 'Content-Type: application/json' \\
-H 'Authorization: Bearer SKtZA/K9ayC0aOlUffNJeLkY3hzE3zKLpZLo6fTo8NeyQLaFjycMF8O/GecWZWkTxJZAr3hGa6JcLDXKRmeLhl1NNaomjwySu8RDzK+8PkLczy9xej/01o7D6lSODWCheVb5q9/yTIcJzp3K9kwJMgdB04t89/1O/w1cDnyilFU=' \\
-d '{{
    "messages":[
        {{
            "type":"text",
            "text":"計算終了通知\\nフォルダ: ~/{dir}\\nファイル: fl{num}"
        }}
    ]
}}'
    '''


if __name__ == "__main__":
    main()
