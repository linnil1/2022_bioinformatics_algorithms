def getAllValues(filename):
    return map(int, open(filename))


def getN50(values):
    values = sorted(values, reverse=True)
    total = sum(values)
    print("Total contig length:", total)
    print("50% contig length:", total * .5)
    acc = 0
    for i in values:
        acc += i
        if acc >= total * .5:
            return i


if __name__:
    values = getAllValues("HW4.1.txt")
    N50 = getN50(values)
    print("N50:", N50)
