import os


input_file = "HW1.txt"
output_file = "HW1_ans.txt"


def lengthOfLine(lines, threshold):
    return sum(len(line) <= threshold for line in lines)


def lineTrim(lines):
    for line in lines:
        yield line.replace("!", "")\
                  .replace("?", "")\
                  .replace("+", "")\
                  .replace("-", "")


def readFile():
    if not os.path.exists(input_file):
        os.system("wget https://raw.githubusercontent.com"
                  "/starsong2266/2022summer-git-python_basics/main"
                  "/HW1/HW1.txt")
    return open(input_file).read().splitlines()  # remove \n or \r


def main():
    with open(output_file, "w") as f:
        total = len(readFile())
        f.write(f"There are {total} sentences in {input_file}. "
                f"{sum(a != b for a, b in zip(readFile(), lineTrim(readFile())))}"
                " of them include !?+-.\n"
                "length before trim after trim\n")
        for threshold in [40, 50, 60, 70, 80]:
            count      = lengthOfLine(         readFile() , threshold)
            trim_count = lengthOfLine(lineTrim(readFile()), threshold)
            f.write(f"  <={threshold} {count:11d} {trim_count:10d}\n")
            if threshold != 80:
                continue
            f.write(f"   >{threshold} {total - count:11d} {total - trim_count:10d}\n")

    os.system(f"cat {output_file}")


if __name__ == "__main__":
    main()
