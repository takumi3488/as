def main():
    print("File number: ", end="")
    num = input().strip()
    file_name = f"fl{num}esd.dat"
    out_text = ""
    with open(file_name) as f:
        for row in f:
            if row.startswith("***"):
                row = row.replace("***", "999")
            out_text += row
    with open(file_name, "w") as f:
        f.write(out_text)
    


if __name__ == "__main__":
    main()
