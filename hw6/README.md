# HW6:

Apply variant on sequences and compare the structural differences.

In the experiments, two variants are apply independently on the ITGA2B sequence:

1. `chr17:42451720A>G`
2. `chr17:42453079G>GT`

and we compare these to official released ITGA2B structure in uniprot.

Also, we change the parameters `num_recycles` and `num_models` to see the impact on AlphaFold2.

```
cd hw6
python hw6.py
```

input:  `ITGA2B.fasta`
Output: `hw6_apply_variant.nucl.fasta` and `hw6_apply_variant.prot.fasta`


## Questions

> What are the TM-score for both structures after insertion/deletion compared with the original structure?

GOTO https://zhanggroup.org/TM-score/

compare `p08514.alphafold.pdb` (treated it as answer)

to 

* `NM_0004195_3060_del_1b526.result.zip/NM_0004195_3060_del_1b526_unrelaxed_rank_1_model_3.pdb` (chr17:42451720A>G)
* `NM_0004195_2605_dupT_2bed2.result.zip/NM_0004195_2605_dupT_2bed2_unrelaxed_rank_1_model_3.pdb` (chr17:42453079G>GT)
* `NM_0004195_cds_f0957.result.zip/NM_0004195_cds_f0957_unrelaxed_rank_1_model_3.pdb` (num_recycles = 3, num_models = 3)
* `NM_0004195_cds_recycle_6_f0957.result.zip/NM_0004195_cds_recycle_6_f0957_unrelaxed_rank_1_model_3.pdb`  (num_recycles = 6)
* `NM_0004195_cds_model5_f0957.result.zip/NM_0004195_cds_model5_f0957_unrelaxed_rank_3_model_5.pdb`  (num_models = 5)

And the TM-score result saved in

* hw6_apply_variant.tmscore.2605_dupT.txt    (TM-socre: 0.7758, lots of difference occurs after the insertion)
* hw6_apply_variant.tmscore.3060_del.txt     (TM-score: 0.8832, most of difference located at last 20 bp)
* hw6_apply_variant.tmscore.cds.txt          (TM-score: 0.6054, I don't know why it's so bad)
* hw6_apply_variant.tmscore.cds_recycle6.txt (TM-score: 0.6088, I don't know why it's so bad)
* hw6_apply_variant.tmscore.cds_model5.txt   (TM-score: 0.8965, It's better)


> What are the amino acid sequences for both structures after insertion/deletion?

`hw6_apply_variant.nucl.fasta` and `hw6_apply_variant.prot.fasta`

> Which change affects the structure more? Explain.

The frameshift, 1-bp insertion, of course.

Also, the deletion occur at the end of structure, which may has minor effect.

> Does recycle process affect the performance? Share what you have done!

Recycle seems not work well in this case.

But increasing `num_models` can make the folding more close to the answer one.


> Three differences between RoseTTAFold and AlphaFold2. (It can be theoritically or differences on Colab)

Skip
