from week3_blast import runShell, threads, runDocker
import glob
import os
from namepipe import nt


def samToBam(name):
    runDocker("quay.io/biocontainers/samtools:1.15.1--h1170115_0",
              f"samtools sort -@{threads} {name}.sam -o {name}.bam")
    runDocker("quay.io/biocontainers/samtools:1.15.1--h1170115_0",
              f"samtools index            {name}.bam")


@nt
def downloadSamples(input_name):
    if os.path.exists("TruSeq3-PE.fa"):
        return "hcc1395.{}.{}"
    runShell("wget http://genomedata.org/rnaseq-tutorial/practical.tar")
    runShell("tar xf practical.tar")
    files = glob.glob("hcc1395*")
    for f in files:
        runShell(f"ln -s {f} {f.replace('_', '.')}")
    runShell("wget https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq3-PE.fa")
    return "hcc1395.{}.{}"


@nt
def downloadIndex(input_name):
    if os.path.exists("GRCh38_latest_genomic.fna.gz"):
        return "GRCh38_latest_genomic"
    runShell("wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz")
    runShell("gunzip -k GRCh38_latest_genomic.fna.gz")  # TODO
    return "GRCh38_latest_genomic"


@nt
def downloadRNAIndex(input_name):
    import week3_blast
    return week3_blast.download()


@nt
def fastQC(input_name):
    output_name = input_name + ".fastqc.{}"
    if os.path.exists(f"{output_name.format('r2')}.zip"):
        return output_name
    runDocker("quay.io/biocontainers/fastqc:0.11.9--0",
              f"fastqc -t {threads} {input_name}.r1.fastq.gz {input_name}.r2.fastq.gz")
    runShell(f"mv {input_name}.r1_fastqc.html {output_name.format('r1')}.html")
    runShell(f"mv {input_name}.r1_fastqc.zip  {output_name.format('r1')}.zip")
    runShell(f"mv {input_name}.r2_fastqc.html {output_name.format('r2')}.html")
    runShell(f"mv {input_name}.r2_fastqc.zip  {output_name.format('r2')}.zip")
    return output_name


@nt
def trimmomatic(input_name):
    output_name = input_name + ".trim"
    if os.path.exists(f"{output_name}.r2.fastq.gz"):
        return output_name
    pe = "TruSeq3-PE.fa"
    runDocker("quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2", f""" \
              trimmomatic PE -threads {threads} \
              -threads {threads} \
              {input_name}.r1.fastq.gz {input_name}.r2.fastq.gz \
              {output_name}.r1.fastq.gz {output_name}.unpair.r1.fastq.gz \
              {output_name}.r2.fastq.gz {output_name}.unpair.r2.fastq.gz \
              ILLUMINACLIP:{pe}:2:30:10:2:keepBothReads \
              LEADING:3 TRAILING:3 MINLEN:36
    """)
    return output_name


@nt
def bwaIndex(input_name):
    output_name = input_name + ".bwa"
    if os.path.exists(f"{output_name}.bwt"):
        return output_name
    runDocker("quay.io/biocontainers/bwa:0.7.17--hed695b0_7",
              f"bwa index {input_name}.fna.gz -p {output_name}")
    return output_name


@nt
def bwaMem(input_name, index):
    output_name = input_name + "." + index.replace("/", "_")
    if os.path.exists(f"{output_name}.bam"):
        return output_name
    runDocker("quay.io/biocontainers/bwa:0.7.17--hed695b0_7", f""" \
              bwa mem -t {threads} {index} \
              {input_name}.r1.fastq.gz {input_name}.r2.fastq.gz \
              -o {output_name}.sam
    """)
    samToBam(output_name)
    return output_name


@nt
def kallistoIndex(input_name):
    output_name = input_name + ".kallisto"
    if os.path.exists(output_name):
        return output_name
    runDocker("quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1", f""" \
        kallisto index {input_name}.fna.gz -i {output_name}
    """)
    return output_name


@nt
def kallistoMap(input_name, index):
    output_name = input_name + "." + index.replace("/", "_")
    if os.path.exists(f"{output_name}.abundance.tsv"):  # TODO
        return output_name
    dir_name = "tmp_" + input_name
    runDocker("quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1", f""" \
        kallisto quant -t {threads} \
        -i {index} \
        --pseudobam --bootstrap-samples 0 --seed 42 \
        -o {dir_name} \
        {input_name}.r1.fastq.gz {input_name}.r2.fastq.gz \
    """)
    runShell(f"mv {dir_name}/pseudoalignments.bam {output_name}.sam")
    samToBam(output_name)
    runShell(f"mv {dir_name}/run_info.json        {output_name}.info.json")
    runShell(f"mv {dir_name}/abundance.h5         {output_name}.abundance.h5")
    runShell(f"mv {dir_name}/abundance.tsv        {output_name}.abundance.tsv")
    runShell(f"rm -rf {dir_name}")
    return output_name


@nt
def hisatIndex(input_name):
    output_name = input_name + ".hisat"
    if os.path.exists(f"{output_name}.1.ht2"):
        return output_name
    runDocker("quay.io/biocontainers/hisat2:2.2.1--h87f3376_4", f""" \
        hisat2-build -p {threads} \
        {input_name}.fna {output_name}
    """)
    return output_name


@nt
def hisatMap(input_name, index):
    output_name = input_name + "." + index.replace("/", "_")
    if os.path.exists(f"{output_name}.bam"):
        return output_name
    runDocker("quay.io/biocontainers/hisat2:2.2.1--h87f3376_4", f""" \
        hisat2 -p {threads} \
        --dta -x {index} \
        -1 {input_name}.r1.fastq.gz -2 {input_name}.r2.fastq.gz \
        -S {output_name}.sam \
    """)
    samToBam(output_name)
    return output_name


if __name__ == "__main__":
    samples = None >> downloadSamples
    samples >> fastQC
    samples = samples >> trimmomatic
    samples >> fastQC
    index_bwa = None >> downloadRNAIndex >> bwaIndex
    index_kallisto = None >> downloadRNAIndex >> kallistoIndex
    index_hisat = None >> downloadIndex >> hisatIndex

    print(samples >> bwaMem.set_args(index=str(index_bwa)))
    print(samples >> kallistoMap.set_args(index=str(index_kallisto)))
    print(samples >> hisatMap.set_args(index=str(index_hisat)))
    """
    real    70m17.174s
    user    1m24.816s
    sys     0m44.045s
    """
