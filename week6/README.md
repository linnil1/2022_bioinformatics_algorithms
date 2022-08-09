# Week6: Protein structure prediction

Run AlphaFold2 colab from https://github.com/sokrypton/ColabFold

Using two sequences (you can generate them from `python week6.py`):

* `Q96L58.fasta`: https://www.uniprot.org/uniprotkb/Q96L58/entry#sequences
* `Q96L58_R232C.fasta`: and change R in position 232 into C

The result of AlphaFold2 is in
* Q96L58.result.zip
* Q96L58_R232C.result.zip

Then run TM-align https://seq2fun.dcmb.med.umich.edu/TM-align/
to compare those two sequences
* Q96L58_40c08_unrelaxed_rank_1_model_3.pdb
* Q96L58_R232C_33493_unrelaxed_rank_1_model_3.pdb
