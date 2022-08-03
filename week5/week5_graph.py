import os
import json
import subprocess
from pip._internal import main as pipmain


try:
    from pyhlamsa import HLAmsa, msaio
except ImportError:
    pipmain(["install", "git+https://github.com/linnil1/pyHLAMSA"])
    from pyhlamsa import HLAmsa, msaio


def runShell(cmd):
    print(cmd)
    proc = subprocess.run(cmd, shell=True, executable="/bin/bash")
    proc.check_returncode()


def runDocker(img, cmd, mount="-v $PWD:/app"):
    runShell(f"podman run -it --rm {mount} -w /app {img} sh -c '{cmd}'")


def prepareData(input_name):
    name = "c3"
    if os.path.exists(f"{name}.vcf.gz"):
        return name

    hla_c = HLAmsa("C", "gen")["C"]
    hla_c_e2e3 = hla_c.select_block(["exon2", "intron2", "exon3"])

    alleles = ["C*01:02:01:01", "C*02:03"]
    c3 = hla_c_e2e3.select_allele(alleles).shrink().reset_index()
    c3 = c3.append("C*consensus", c3.get_consensus())
    c3.set_reference("C*consensus")
    msaio.save_msa(c3, f"{name}.fa", f"{name}.json")
    msaio.to_fasta(c3, f"{name}.nogap.fa",                                 gap=False)
    msaio.to_fasta(c3.select_allele(["C*consensus"]), f"{name}.ref.fa",    gap=False)
    msaio.to_fasta(c3.select_allele(alleles),         f"{name}.sample.fa", gap=False)
    msaio.to_vcf  (c3, f"{name}.vcf.gz")
    return name


def simulateFastq(input_name):
    name = input_name
    if os.path.exists(f"{name}.fq"):
        return name
    runShell("wget https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier2016.06.05linux64.tgz")
    runShell("tar xvf artbinmountrainier2016.06.05linux64.tgz")
    runDocker("alpine",
              f"./art_bin_MountRainier/art_illumina -ss HS25 "
              f"-i {name}.fa  -l 150 -f 30 -sam -na -o {name}")
    return name


def buildHisat(input_name):
    output_name = input_name + ".hisat"
    if os.path.exists(f"{output_name}.8.ht2"):
        return output_name
    runShell(f"python hisat2_extract_snps_haplotypes_VCF.py "
             f"{input_name}.ref.fa {input_name}.vcf.gz {output_name}")
    runDocker("quay.io/biocontainers/hisat2:2.2.1--he1b5a44_2",
              f"hisat2-build {input_name}.ref.fa "
              f"--snp {output_name}.snp --haplotype {output_name}.haplotype "
              f"{output_name} ")
    return output_name


def hisatMap(input_name, index):
    output_name = input_name + "." + index.replace("/", "_").replace(".", "_")
    if os.path.exists(f"{output_name}.bam"):
        return output_name
    runDocker("quay.io/biocontainers/hisat2:2.2.1--he1b5a44_2",
              f"hisat2 -x {index} -U {input_name}.fq -S {output_name}.sam")
    runDocker("quay.io/biocontainers/samtools:1.15.1--h1170115_0",
              f"samtools sort {output_name}.sam -o      {output_name}.bam")
    runDocker("quay.io/biocontainers/samtools:1.15.1--h1170115_0",
              f"samtools index                          {output_name}.bam")
    return output_name


def buildVg(input_name):
    output_name = input_name + ".vg"
    if os.path.exists(f"{output_name}.vis.xg"):
        return output_name
    runDocker("quay.io/vgteam/vg:v1.40.0",
              f"vg autoindex --workflow giraffe "
              f"-r {input_name}.ref.fa -v {input_name}.vcf.gz -p {output_name}")
    runDocker("quay.io/vgteam/vg:v1.40.0",
              f"vg convert {output_name}.giraffe.gbz -x > {output_name}.vis.xg")
    return output_name


def vgMap(input_name, index):
    output_name = input_name + "." + index.replace("/", "_").replace(".", "_")
    output_name_sort = output_name + ".sort"
    if os.path.exists(f"{output_name_sort}.gam"):
        return output_name_sort

    # bam
    runDocker("quay.io/vgteam/vg:v1.40.0",
              f"vg giraffe -Z {index}.giraffe.gbz "
              f"-f {input_name}.fq -o BAM > {output_name}.bam")
    runDocker("quay.io/biocontainers/samtools:1.15.1--h1170115_0",
              f"samtools sort {output_name}.bam -o {output_name_sort}.bam")
    runDocker("quay.io/biocontainers/samtools:1.15.1--h1170115_0",
              f"samtools index                     {output_name_sort}.bam")
    # gam
    runDocker("quay.io/vgteam/vg:v1.40.0",
              f"vg giraffe -Z {index}.giraffe.gbz "
              f"-f {input_name}.fq -o gam > {output_name}.gam")
    runDocker("quay.io/vgteam/vg:v1.40.0",
              f"vg gamsort {output_name}.gam "
              f"-i {output_name_sort}.gam.gai > {output_name_sort}.gam")
    return output_name_sort


def setupTube(input_name, index):
    if os.path.exists(f"sequenceTubeMap/exampleData/{input_name}.gam"):
        return output_name_sort
    runShell("git clone https://github.com/vgteam/sequenceTubeMap.git")
    runShell("wget https://github.com/vgteam/vg/releases/download/v1.42.0/vg")
    runShell(f"cp {index}.vis.xg sequenceTubeMap/exampleData")
    runShell(f"cp {input_name}.gam sequenceTubeMap/exampleData")
    runShell(f"cp {input_name}.gam.gai sequenceTubeMap/exampleData")
    a = json.load(open("sequenceTubeMap/src/config.json"))
    a["DATA_SOURCES"].append({
        'name': input_name,
        'xgFile': index + ".vis.xg",
        'gamFile': input_name + ".gam",
        'dataPath': "mounted",
        'region': "C*consensus:1-100",
        'dataType': "built-in",
    })
    json.dump(a, open("sequenceTubeMap/src/config.json", "w"))
    runShell("chmod +x vg")
    runDocker("node:18-alpine", "yarn",         mount="-v $PWD/sequenceTubeMap:/app")
    runDocker("node:18-alpine", "yarn install", mount="-v $PWD/sequenceTubeMap:/app")
    runDocker("node:18-alpine", "yarn build",   mount="-v $PWD/sequenceTubeMap:/app")


if __name__ == "__main__":
    main = prepareData(None)
    sample = simulateFastq(main + '.sample')

    hisat_index = buildHisat(main)
    hisat_map = hisatMap(sample, index=hisat_index)
    # Skip IGV
    vg_index = buildVg(main)
    vg_map = vgMap(sample, index=vg_index)
    setupTube(vg_map, vg_index)
    # Start sequenceTubeMap service
    print("Run this")
    print("docker run -it --rm -w /app -v $PWD/sequenceTubeMap:/app -v $PWD/vg:/bin/vg -p 8000:3000 node:18-alpine yarn start")
