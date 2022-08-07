import os
import subprocess
from collections import Counter
import pysam
import pandas as pd


def runShell(cmd: str):
    print(cmd)
    proc = subprocess.run(cmd, shell=True)
    proc.check_returncode()


def download():
    if os.path.exists("chr17.fa"):
        return "chr17"
    runShell("wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr17.fa.gz")  # noqa
    runShell("gunzip chr17.fa.gz")
    return "chr17"


def runDocker(img: str, cmd: str):
    runShell(f"podman run -it --rm -v $PWD:/app -w /app {img} {cmd}")


def samtobam(name: str) -> str:
    if os.path.exists(f"{name}.bam"):
        return name
    runDocker("quay.io/biocontainers/samtools:1.15.1--h1170115_0",
              f"samtools sort -t 4 {name}.sam -o {name}.bam")
    runDocker("quay.io/biocontainers/samtools:1.15.1--h1170115_0",
              f"samtools index {name}.bam")
    return name


def pileUp(filename: str, chrom: str, pos: int) -> dict[str, str]:
    f = pysam.AlignmentFile(filename)
    filename_out = filename.rsplit(".")[0] + f".{chrom}.{pos}.strange_case.sam"
    f_out = pysam.AlignmentFile(filename_out, "w", template=f)

    reads = {}
    for column in f.pileup(chrom, pos - 1, pos, truncate=True):  # pileup is 0-based
        for read in column.pileups:
            id = read.alignment.query_name
            if read.indel != 0:
                # insertion
                if read.indel > 0:
                    reads[id] = read.alignment.query_sequence[
                                        read.query_position:
                                        read.query_position + read.indel + 1]
                # deletion
                else:
                    reads[id] = read.indel
                    f_out.write(read.alignment)
            elif read.is_del:
                reads[id] = "del"
                f_out.write(read.alignment)
            elif read.is_refskip:
                reads[id] = "N"
                f_out.write(read.alignment)
            elif read.query_position is None:
                reads[id] = None
                f_out.write(read.alignment)
            else:
                reads[id] = read.alignment.query_sequence[read.query_position]
        break  # only selected column
    return reads


def getTargetPosition(file_vcf: str) -> list[tuple[str, int]]:
    targets = []
    for i in open(file_vcf):
        if not i:
            continue
        chrom, pos, _, ref, alt = i.split("\t")[:5]  # not need to pos - 1
        targets.append((chrom, int(pos)))
    return targets


def calcPhaseOfVariants(sample_name: str, vcf_name: str):
    # read vcf's variants
    # and read the variant at the position via pileup
    var_stat = []
    for chrom, pos in getTargetPosition(vcf_name + ".vcf"):
        var_stat.append((chrom, pos, pileUp(sample_name + ".bam", chrom, pos)))
        print(f"{chrom}:{pos} -> {Counter(var_stat[-1][2].values())}")
        samtobam(f"{sample_name}.{chrom}.{pos}.strange_case")

    # intersect the read name
    read_set = set.intersection(*[set(stat[2].keys()) for stat in var_stat])
    print("phased_reads={len(read_set)}")
    for chrom, pos, stat in var_stat:
        print(f"{chrom}:{pos} reads={len(stat)} "
              f"unphased_reads={len(stat) - len(read_set)}")

    # merge variants to count
    read_phase_var = {}
    for read in read_set:
        read_phase_var[read] = tuple(stat[2][read] for stat in var_stat)

    # count
    count = Counter(read_phase_var.values())
    print(count)

    # to dataframe
    dfs = []
    for key, value in count.items():
        d = {f"{chrom}:{pos}": k for k, (chrom, pos, _) in zip(key, var_stat)}
        d['count'] = value
        dfs.append(d)
    df = pd.DataFrame(dfs)
    df = df.sort_values(by="count", ascending=False)
    print(df)
    df.to_csv(f"{sample_name}.{vcf_name}.count.tsv", sep="\t", index=False)


if __name__ == "__main__":
    # download()
    sample = "hw5"
    sample = samtobam(sample)
    # pileUp(sample + ".bam", "chr17", 42451720)
    # pileUp(sample + ".bam", "chr17", 42453079)
    calcPhaseOfVariants(sample, "ITGA2B")
