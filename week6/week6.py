import copy
import subprocess
from Bio import SeqIO, Seq


def runShell(cmd):
    print(cmd)
    proc = subprocess.run(cmd, shell=True, executable="/bin/bash")
    proc.check_returncode()


if __name__ == "__main__":
    # download
    runShell("wget https://rest.uniprot.org/uniprotkb/Q96L58.fasta") 

    # read record
    new_record = None
    record = list(SeqIO.parse("Q96L58.fasta", "fasta"))
    assert len(record) == 1
    record = record[0]
    print(record.id.split("|")[1], str(record.seq))

    # R232C
    assert record.seq[231] == "R"
    new_record = copy.deepcopy(record)
    new_record.seq = Seq.Seq(record.seq[:231] + "C" + record.seq[232:])
    new_record.id = new_record.id.replace("Q96L58", "Q96L58_R232C")
    assert record.seq[231] == "C"
    print(new_record.id.split("|")[1], str(new_record.seq))

    # save
    SeqIO.write([new_record], "Q96L58_R232C.fasta", "fasta")
