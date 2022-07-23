import os
import subprocess
import sys
import pandas as pd
import numpy as np
from Bio import SeqIO

sys.path.append("../week3")  # noqa
from week3_blast import runShell, runDocker, buildBlast

threads = 30


def download():
    name = "GCF_000001895.5_Rnor_6.0_genomic"
    if os.path.exists(name + ".fna"):
        return name
    runShell("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/895/GCF_000001895.5_Rnor_6.0/GCF_000001895.5_Rnor_6.0_genomic.fna.gz")  # noqa
    runShell(f"gunzip -k {name}.fna.gz")
    return name


def downloadRNA():
    name = "GCF_000001895.5_Rnor_6.0_rna"
    if os.path.exists(name + ".fna"):
        return name
    runShell("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/895/GCF_000001895.5_Rnor_6.0/GCF_000001895.5_Rnor_6.0_rna.fna.gz")  # noqa
    runShell(f"gunzip -k {name}.fna.gz")
    return name


def queryBlast(input_name, index):
    output_name = input_name + "." + index.replace("/", "_").replace('.', "_")
    if os.path.exists(f"{output_name}.csv"):
        return output_name
    runDocker("quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0", f""" \
              blastn -task blastn \
              -num_threads {threads} \
              -db {index} \
              -query {input_name}.fa \
              -outfmt 7 \
              -out {output_name}.csv
    """)
    return output_name


def statPair(df):
    if set(df['query_dir']) != set(["F", "R"]):
        return None
    # print(df)
    # I don't deal with multiple alignment on same query same subject
    pos = sorted(list(df["subject_end"]) + list(df["subject_start"]))
    assert len(pos) == 4
    df["insert_size"] = pos[2] - pos[1] - 1
    return df[["query_acc", "query_gene", "subject_acc",
               "subject_start", "subject_end",
               "insert_size"]].reset_index(drop=True)


def targetInPair(data_fr):
    if not len(data_fr):
        return pd.DataFrame()
    set_fr = data_fr.groupby(["query_gene", "query_dir"])["subject_acc"] \
                    .agg(lambda i: set(i)) \
                    .groupby("query_gene").agg(lambda i: set.intersection(*i))
    data_fr = data_fr[data_fr.apply(lambda i: i.subject_acc in set_fr[i.query_gene],
                                    axis=1)]
    data_fr_stat = data_fr.groupby(["query_gene", "subject_acc"]).apply(statPair)
    # print(data_fr_stat)
    return data_fr_stat


def stats(input_name):
    output_name = input_name + ".stat"
    if os.path.exists(f"{output_name}.csv"):
        return output_name

    data = pd.read_csv(input_name + ".csv", sep="\t", comment="#", names=[
        "query_acc", "subject_acc", "identity", "alignment_length", "mismatches",
        "gap_opens", "query_start", "query_end", "subject_start", "subject_end",
        "evalue", "bitscore"])
    data[["query_gene", "query_dir"]] = data['query_acc'].str.split('_', expand=True)

    # QA: 100% identity
    data = data[data['identity'] == 100]

    # QB: full alignment length
    query_length = {i.id: len(i.seq)
                    for i in SeqIO.parse(input_name.split('.')[0] + ".fa", "fasta")}
    data['real_length'] = list(map(query_length.get, data['query_acc']))
    data = data[data['alignment_length'] == data['real_length']]

    # QC: correct directions(Forward and Reverse)
    mask_same_dir = np.logical_or(
        np.logical_and(data['query_dir'] == "F",
                       data['subject_end'] - data['subject_start'] > 0),
        np.logical_and(data['query_dir'] == "R",
                       data['subject_end'] - data['subject_start'] < 0)
    )
    mask_diff_dir = np.logical_or(
        np.logical_and(data['query_dir'] == "R",
                       data['subject_end'] - data['subject_start'] > 0),
        np.logical_and(data['query_dir'] == "F",
                       data['subject_end'] - data['subject_start'] < 0)
    )

    # QD: in both F and R
    data_pair = pd.concat([
        targetInPair(data[mask_same_dir]),
        targetInPair(data[mask_diff_dir]),
    ], ignore_index=True)
    print(data_pair)
    data_pair = data_pair.groupby(["query_gene", "subject_acc"],
                                  as_index=False)["insert_size"].first()
    data_pair["ref"] = input_name
    # print(data_pair)
    data_pair.to_csv(output_name + ".csv", index=False)
    return output_name


def mergeResult(*names):
    input_name = "rat_sample.{}.stat"
    output_name = input_name.replace(".{}", "_merge")
    if os.path.exists(f"{output_name}.csv"):
        return output_name

    df = pd.concat([pd.read_csv(name + ".csv") for name in names])
    df_dna = df[df['ref'].str.contains("_genom")].rename(
            columns={'subject_acc': "chrom_ID"})
    df_rna = df[df['ref'].str.contains("_rna")].rename(
            columns={'subject_acc': "transcript_ID"})
    df = pd.concat([df_dna, df_rna]).drop(columns="ref")
    df.to_csv(output_name + ".csv", index=False)
    return output_name


if __name__ == "__main__":
    index = download()
    index = buildBlast(index)
    sample = "rat_sample"
    sample = queryBlast(sample, index)
    sample = stats(sample)
    sample_dna = sample

    index_rna = downloadRNA()
    index_rna = buildBlast(index_rna)
    sample = "rat_sample"
    sample = queryBlast(sample, index_rna)
    sample = stats(sample)
    sample_rna = sample

    # force to merge
    sample = mergeResult(sample_rna, sample_dna)
    print(pd.read_csv(sample + ".csv"))
