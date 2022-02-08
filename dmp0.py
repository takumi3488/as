def main():
    t = ""
    fn = "wkfermi/band29B.dat"
    with open(fn) as f:
        for row in f:
            if all([int(x) % 6 == 0 for x in row.split()[:3]]):
                t += row
    with open(fn, "w") as f:
        f.write(t)


if __name__ == "__main__":
    main()
