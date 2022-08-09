import os
import subprocess
from glob import glob

threads = 4


def runShell(cmd):
    print(cmd)
    proc = subprocess.run(cmd, shell=True, executable="/bin/bash")
    proc.check_returncode()


def runDocker(img, cmd, mount="-v $PWD:/app"):
    runShell(f"podman run -it --rm {mount} -w /app {img} sh -c '{cmd}'")


def download():
    if os.path.exists(f"data/P7741.r1.fastq.gz"):
        return "data/P7741"
    runShell("git clone https://github.com/vappiah/bacterial-genomics-tutorial")
    runShell("bash ./bacterial-genomics-tutorial/download_data.sh")
    runShell("ln -s P7741_R1.fastq.gz data/P7741.r1.fastq.gz")
    runShell("ln -s P7741_R2.fastq.gz data/P7741.r2.fastq.gz")
    return "data/P7741"


def fastQC(input_name):
    output_name = input_name + ".fastqc"
    if os.path.exists(f"{output_name}.r2.zip"):
        return output_name
    runDocker("quay.io/biocontainers/fastqc:0.11.9--0",
              f"fastqc -t {threads} {input_name}.r1.fastq.gz {input_name}.r2.fastq.gz")
    runShell(f"mv {input_name}.r1_fastqc.html {output_name}.r1.html")
    runShell(f"mv {input_name}.r1_fastqc.zip  {output_name}.r1.zip")
    runShell(f"mv {input_name}.r2_fastqc.html {output_name}.r2.html")
    runShell(f"mv {input_name}.r2_fastqc.zip  {output_name}.r2.zip")
    return output_name


def trim(input_name):
    output_name = input_name + ".trim"
    if os.path.exists(f"{output_name}.r2.fastq.gz"):
        return output_name

    runDocker("quay.io/biocontainers/sickle-trim:1.33--h7132678_7",
              f"sickle pe -g -f {input_name}.r1.fastq.gz -r {input_name}.r2.fastq.gz "
              f"-t sanger -q 20 -l 20 -s {input_name}.trim.s.fastq.gz "
              f"-o {input_name}.trim.r1.fastq.gz -p {input_name}.trim.r2.fastq.gz")
    return output_name


def assemble(input_name):
    output_name = input_name + ".spades"
    output_name1 = output_name + ".contig"
    if os.path.exists(f"{output_name1}.fasta"):
        return output_name1

    runDocker("quay.io/biocontainers/spades:3.15.5--h95f258a_0",
              "spades.py --careful "
              f"-1 {input_name}.r1.fastq.gz -2 {input_name}.r2.fastq.gz "
              f"-s {input_name}.s.fastq.gz  -o {output_name}")
    runShell(f"ln -s {output_name.rsplit('/', 1)[1]}/contigs.fasta {output_name1}.fasta")
    return output_name1


def reorder(input_name, genome):
    output_name = input_name + ".reorder"
    if os.path.exists(f"{output_name}.fasta"):
        return output_name
    runDocker("quay.io/biocontainers/ragtag:2.1.0--pyhb7b1952_0",
              f"ragtag.py scaffold {genome} {input_name}.fasta -o {output_name}")
    runShell("python bacterial-genomics-tutorial/extract_reordered.py "
             f"{output_name}/ragtag.scaffold.fasta NC_020133.1")
    runShell(f"mv P7741.reordered.fasta {output_name}.fasta")
    return output_name


def checkAMRGenes(input_name):
    output_name = input_name + ".amr"
    if os.path.exists(f"{output_name}.tab"):
        return output_name
    runDocker("quay.io/biocontainers/abricate:1.0.1--ha8f3691_1",
              f"abricate {input_name}.fasta > {output_name}.tab")
    runShell(f"cat {output_name}.tab")
    return output_name


def annotate(input_name):
    output_name = input_name + ".annotate"
    if os.path.exists(f"{output_name}.faa"):
        return output_name
    runDocker("quay.io/biocontainers/prokka:1.14.6--pl5321hdfd78af_4",
              f"prokka --kingdom Bacteria --locustag P7741 --prefix P7741 "
              f"--cpus {threads} --addgenes {input_name}.fasta --outdir {output_name}")
    runShell(f"perl bacterial-genomics-tutorial/get_pseudo.pl {output_name}/P7741.faa "
             f"> {output_name}/P7741.pseudo.txt")
    runShell(f"ln -s {output_name.rsplit('/', 1)[1]}/P7741.faa {output_name}.faa")
    return output_name


def showAnnotate(input_name):
    runShell(f"python bacterial-genomics-tutorial/get_annot_stats.py {input_name} P7741")


def plotDendogram(input_name):
    output_name = input_name + ".dendogram"
    if os.path.exists(f"{output_name}"):
        return output_name
    genomes = glob("bacterial-genomics-tutorial/genomes/*.fasta")
    runDocker("quay.io/biocontainers/drep:3.4.0--pyhdfd78af_0",
              f"dRep compare dendogram -p {threads} -g {input_name}.fasta "
              + " ".join(genomes))
    runShell(f"mv dendogram {output_name}")
    return output_name


if __name__ == "__main__":
    sample = download()
    fastQC(sample)
    sample = trim(sample)
    fastQC(sample)
    sample = assemble(sample)
    sample = reorder(sample, genome="bacterial-genomics-tutorial/genomes/Liflandii.fasta")
    checkAMRGenes(sample)
    sample_annot = annotate(sample)
    showAnnotate(sample_annot)
    plotDendogram(sample)
