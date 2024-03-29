# 2022 Bioinformatics Algorithms

Setup (make sure you have python3.10)
```
pip install -r requirements.txt
```

## HW1: Python practice
```
cd hw1
python hw1.py
```

Input:  `hw1/HW1.txt`
Output: `hw1/HW1_ans.txt`


## HW2: Python practice II

Read Clinvar VCF

```
cd hw2
python hw2.py
```

Input:  `hw2/HW2_clinvar.txt`
Output: `hw2/hw2_clnsig.png` `hw2/hw2_rs.tsv`


## HW2_2: Bash practice
```
cd hw2_2
./hw2_2.sh
```

Input:  `/localpath/a_very_big_trash.tar.gz`
Output: `hw2_2/answer/home.txt`, `hw2_2/answer/location.txt`


## Week3: Blast and read mapping

see [week3/README.md](https://github.com/linnil1/2022_bioinformatics_algorithms/tree/main/week3)

## HW3: Read Blast Result

Some pair of sequences are mapped on the reference by blast,

and then calculate the insertion size of the pair.

Some pair should be removed if the mapping result doesn't meet some criteria.

```
cd hw3
python hw3.py
```

input:  `hw3/rat_sample.fa`
Output: `hw3/rat_sample_merge.stat.csv`


## Week4: Metagenomics

see [week4/README.md](https://github.com/linnil1/2022_bioinformatics_algorithms/tree/main/week4)


## HW4: N50 and Metagenomics Result Interpretation

```
cd hw4
# N50 practice
python hw4_1.py > hw4.1.txt
# Parse Result of https://github.com/metagenome-atlas/Tutorial.git and answer questions
python hw4_2.py > hw4.2.txt
```

input: `HW4.1.txt` `Tutorial/Example`
Output: `hw4.1.txt` `HW4.2.txt` `hw4.2.png`


## HW5: Read Sam/Bam format

Extract the phased genotype of specific position from sam file.

```
cd hw5
python hw5.py
```

input:  `hw5.sam` (Not provided), `ITGA2B.vcf`
Output: `hw5.ITGA2B.count.tsv` `hw5.chr17.*.strange_case.bam`


## Week5: Graph genome

see [week5/README.md](https://github.com/linnil1/2022_bioinformatics_algorithms/tree/main/week5)

```
cd week5
python week5_graph.py
```

input:  `NA`
Output: `week5/c3.sample.c3_vg.sort.gam.gai` `week5/c3.sample.c3_hisat.bam`


## Week6: Protein structure prediction

see [week6/README.md](https://github.com/linnil1/2022_bioinformatics_algorithms/tree/main/week6)

## HW6: Protein structure prediction

see [hw6/README.md](https://github.com/linnil1/2022_bioinformatics_algorithms/tree/main/hw6)

## Final: Final exam (coding)

see [final/README.md](https://github.com/linnil1/2022_bioinformatics_algorithms/tree/main/final)
