import os
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


table_name = "Standard"


def runShell(cmd):
    print(cmd)
    proc = subprocess.run(cmd, shell=True, executable="/bin/bash")
    proc.check_returncode()


def downloadProtein():
    name = "P08514"
    if os.path.exists(name + ".fasta"):
        return name
    runShell("wget https://alphafold.ebi.ac.uk/files/AF-P08514-F1-model_v3.pdb"
             f" -O {name}.alphafold.pdb")
    runShell(f"wget https://rest.uniprot.org/uniprotkb/{name}.fasta")
    return name


def downloadmRNA():
    name = "ITGA2B"
    if os.path.exists(name + ".fasta"):
        return name
    runShell("wget 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=1666305226'"
             f" -O {name}.gff")
    runShell("wget 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=Nucleotide&id=NM_000419&&rettype=fasta&retmode=text'"
             f" -O {name}.fasta")
    return name


def getCDSPos(gff_file) -> tuple[int, int]:
    gff = open(gff_file)
    for line in gff:
        if not line.strip() or line.startswith("#"):
            continue
        field = line.split("\t")
        if field[2] == "CDS":
            return (int(field[3]), int(field[4]))
    assert False


if __name__ == "__main__":
    name_nucl = downloadmRNA()
    name_prot = downloadProtein()
    output_name = "hw6_apply_variant"

    # read nuclotide
    # https://www.ncbi.nlm.nih.gov/nuccore/NM_000419
    record = list(SeqIO.parse(name_nucl + ".fasta", "fasta"))[0]
    cds_pos = getCDSPos(name_nucl + ".gff")
    seq = record.seq[cds_pos[0] - 1:cds_pos[1]]

    # transslate and double check
    # https://www.uniprot.org/uniprotkb/P08514/entry
    print("NM_000419.5 (CDS part)")
    seq_translate = seq.translate(table_name, to_stop=True)
    assert str(list(SeqIO.parse(name_prot + ".fasta", "fasta"))[0].seq) \
           == str(seq_translate)
    print(str(seq_translate))

    # case1
    print("NM_000419.5(ITGA2B):c.3060+2T>C = chr17:42451720A>G")
    # https://www.ncbi.nlm.nih.gov/clinvar/RCV001225236.1
    # assume it's deletion
    # WTF the position is not correct. It's related to CDS not NM_000419.5
    # cds_pos_del = (cds_pos[0], min(cds_pos[1], 3060 + cds_pos[0] - 1))
    # seq_del = record.seq[cds_pos_del[0] - 1:cds_pos_del[1]]
    seq_del = seq[:3060]
    seq_del_translate = seq_del.translate(table_name, to_stop=True)
    print(str(seq_del_translate))

    # case2
    print(">NM_000419.5_2605_dupT = chr17:42453079G>GT")
    # https://www.ncbi.nlm.nih.gov/snp/rs1251095196
    pos = 868 * 3 + 1  # I manually counted
    seq_ins = seq[:pos] + "A" + seq[pos:]
    print(seq[pos-5:pos+6].reverse_complement())
    print(seq_ins[pos-5:pos+7].reverse_complement())
    seq_ins_translate = seq_ins.translate(table_name, to_stop=True)
    print(str(seq_ins_translate))

    SeqIO.write([
        SeqRecord(seq,     id="NM_000419.5_cds",       description=""),
        SeqRecord(seq_del, id="NM_000419.5_3060_del",  description=""),
        SeqRecord(seq_ins, id="NM_000419.5_2605_dupT", description=""),
    ], output_name + ".nucl.fasta", "fasta")
    SeqIO.write([
        SeqRecord(seq_translate,     id="NM_000419.5_translate", description=""),
        SeqRecord(seq_del_translate, id="NM_000419.5_3060_del",  description=""),
        SeqRecord(seq_ins_translate, id="NM_000419.5_2605_dupT", description=""),
    ], output_name + ".prot.fasta", "fasta")
