import os
import subprocess

threads = 4


def runShell(cmd):
    print(cmd)
    proc = subprocess.run(cmd, shell=True)
    proc.check_returncode()


def download():
    if os.path.exists("GRCh38_latest_rna.fna"):
        return "GRCh38_latest_rna"
    runShell("wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_rna.fna.gz")
    runShell("gunzip -k GRCh38_latest_rna.fna.gz")
    return "GRCh38_latest_rna"


def runDocker(img, cmd):
    runShell(f"podman run -it --rm -v $PWD:/app -w /app {img} {cmd}")


def buildBlast(input_name):
    output_name = input_name + ".blast"
    if os.path.exists(f"{output_name}.ndb"):
        return output_name
    runDocker("quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0", f""" \
              makeblastdb -dbtype nucl -parse_seqids \
              -in {input_name}.fna \
              -out {output_name}
    """)
    return output_name


def extractFasta(input_name, index):
    if os.path.exists(f"{input_name}.fa"):
        return input_name
    runDocker("quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0",
              f"blastdbcmd -db {index} -entry {input_name} > {input_name}.fa")
    return input_name


def queryBlast(input_name, index):
    output_name = input_name + "." + index.replace("/", "_")
    if os.path.exists(f"{output_name}.csv"):
        return output_name
    runDocker("quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0", f"""\
              blastn -num_threads {threads} \
              -task blastn \
              -query {input_name}.fa -db {index} \
              -evalue 0.05 \
              -out {output_name}.txt
    """)
    runDocker("quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0", f""" \
              blastn -task blastn \
              -db {index} \
              -query {input_name}.fa \
              -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" \
              -out {output_name}.csv
    """)
    return output_name


if __name__ == "__main__":
    index = download()
    index = buildBlast(index)
    sample = "NM_007297.4"
    sample = extractFasta(sample, index)
    sample = queryBlast(sample, index)
    """
    real    1m30.109s
    user    0m5.402s
    sys     0m2.703s
    """
