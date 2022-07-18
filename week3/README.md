# Week3

## Run it
```
cd week3
# blast takes 1.5 mins 
python week3_blast.py  
# fastqc + trimmomatic + bwa + hisat2 + kallisto takes 70 mins
python week3_mapping.py
```

## Input:

* `week3/GRCh38_latest_rna.fna.gz`
* `week3/GRCh38_latest_genomic.fna.gz`
* `http://genomedata.org/rnaseq-tutorial/practical.tar`

## Output:

```
.
├── GRCh38_latest_genomic.fna
├── GRCh38_latest_genomic.fna.gz
├── GRCh38_latest_genomic.hisat.1.ht2
├── GRCh38_latest_genomic.hisat.2.ht2
├── GRCh38_latest_genomic.hisat.3.ht2
├── GRCh38_latest_genomic.hisat.4.ht2
├── GRCh38_latest_genomic.hisat.5.ht2
├── GRCh38_latest_genomic.hisat.6.ht2
├── GRCh38_latest_genomic.hisat.7.ht2
├── GRCh38_latest_genomic.hisat.8.ht2
├── GRCh38_latest_rna.blast.ndb
├── GRCh38_latest_rna.blast.nhr
├── GRCh38_latest_rna.blast.nin
├── GRCh38_latest_rna.blast.nog
├── GRCh38_latest_rna.blast.nos
├── GRCh38_latest_rna.blast.not
├── GRCh38_latest_rna.blast.nsq
├── GRCh38_latest_rna.blast.ntf
├── GRCh38_latest_rna.blast.nto
├── GRCh38_latest_rna.bwa.amb
├── GRCh38_latest_rna.bwa.ann
├── GRCh38_latest_rna.bwa.bwt
├── GRCh38_latest_rna.bwa.pac
├── GRCh38_latest_rna.bwa.sa
├── GRCh38_latest_rna.fna
├── GRCh38_latest_rna.fna.gz
├── GRCh38_latest_rna.kallisto
├── hcc1395.normal.rep1.fastqc.r1.html
├── hcc1395.normal.rep1.fastqc.r1.zip
├── hcc1395.normal.rep1.fastqc.r2.html
├── hcc1395.normal.rep1.fastqc.r2.zip
├── hcc1395.normal.rep1.r1.fastq.gz -> hcc1395_normal_rep1_r1.fastq.gz
├── hcc1395_normal_rep1_r1.fastq.gz
├── hcc1395.normal.rep1.r2.fastq.gz -> hcc1395_normal_rep1_r2.fastq.gz
├── hcc1395_normal_rep1_r2.fastq.gz
├── hcc1395.normal.rep1.trim.fastqc.r1.html
├── hcc1395.normal.rep1.trim.fastqc.r1.zip
├── hcc1395.normal.rep1.trim.fastqc.r2.html
├── hcc1395.normal.rep1.trim.fastqc.r2.zip
├── hcc1395.normal.rep1.trim.GRCh38_latest_genomic.hisat.bam
├── hcc1395.normal.rep1.trim.GRCh38_latest_genomic.hisat.bam.bai
├── hcc1395.normal.rep1.trim.GRCh38_latest_rna.bwa.bam
├── hcc1395.normal.rep1.trim.GRCh38_latest_rna.bwa.bam.bai
├── hcc1395.normal.rep1.trim.GRCh38_latest_rna.kallisto.abundance.h5
├── hcc1395.normal.rep1.trim.GRCh38_latest_rna.kallisto.abundance.tsv
├── hcc1395.normal.rep1.trim.GRCh38_latest_rna.kallisto.bam
├── hcc1395.normal.rep1.trim.GRCh38_latest_rna.kallisto.bam.bai
├── hcc1395.normal.rep1.trim.GRCh38_latest_rna.kallisto.info.json
├── hcc1395.normal.rep1.trim.r1.fastq.gz
├── hcc1395.normal.rep1.trim.r2.fastq.gz
├── hcc1395.normal.rep1.trim.unpair.r1.fastq.gz
├── hcc1395.normal.rep1.trim.unpair.r2.fastq.gz
├── hcc1395.normal.rep2.fastqc.r1.html
├── hcc1395.normal.rep2.fastqc.r1.zip
├── hcc1395.normal.rep2.fastqc.r2.html
├── hcc1395.normal.rep2.fastqc.r2.zip
├── hcc1395.normal.rep2.r1.fastq.gz -> hcc1395_normal_rep2_r1.fastq.gz
├── hcc1395_normal_rep2_r1.fastq.gz
├── hcc1395.normal.rep2.r2.fastq.gz -> hcc1395_normal_rep2_r2.fastq.gz
├── hcc1395_normal_rep2_r2.fastq.gz
├── hcc1395.normal.rep2.trim.fastqc.r1.html
├── hcc1395.normal.rep2.trim.fastqc.r1.zip
├── hcc1395.normal.rep2.trim.fastqc.r2.html
├── hcc1395.normal.rep2.trim.fastqc.r2.zip
├── hcc1395.normal.rep2.trim.GRCh38_latest_genomic.hisat.bam
├── hcc1395.normal.rep2.trim.GRCh38_latest_genomic.hisat.bam.bai
├── hcc1395.normal.rep2.trim.GRCh38_latest_rna.bwa.bam
├── hcc1395.normal.rep2.trim.GRCh38_latest_rna.bwa.bam.bai
├── hcc1395.normal.rep2.trim.GRCh38_latest_rna.kallisto.abundance.h5
├── hcc1395.normal.rep2.trim.GRCh38_latest_rna.kallisto.abundance.tsv
├── hcc1395.normal.rep2.trim.GRCh38_latest_rna.kallisto.bam
├── hcc1395.normal.rep2.trim.GRCh38_latest_rna.kallisto.bam.bai
├── hcc1395.normal.rep2.trim.GRCh38_latest_rna.kallisto.info.json
├── hcc1395.normal.rep2.trim.r1.fastq.gz
├── hcc1395.normal.rep2.trim.r2.fastq.gz
├── hcc1395.normal.rep2.trim.unpair.r1.fastq.gz
├── hcc1395.normal.rep2.trim.unpair.r2.fastq.gz
├── hcc1395.normal.rep3.fastqc.r1.html
├── hcc1395.normal.rep3.fastqc.r1.zip
├── hcc1395.normal.rep3.fastqc.r2.html
├── hcc1395.normal.rep3.fastqc.r2.zip
├── hcc1395.normal.rep3.r1.fastq.gz -> hcc1395_normal_rep3_r1.fastq.gz
├── hcc1395_normal_rep3_r1.fastq.gz
├── hcc1395.normal.rep3.r2.fastq.gz -> hcc1395_normal_rep3_r2.fastq.gz
├── hcc1395_normal_rep3_r2.fastq.gz
├── hcc1395.normal.rep3.trim.fastqc.r1.html
├── hcc1395.normal.rep3.trim.fastqc.r1.zip
├── hcc1395.normal.rep3.trim.fastqc.r2.html
├── hcc1395.normal.rep3.trim.fastqc.r2.zip
├── hcc1395.normal.rep3.trim.GRCh38_latest_genomic.hisat.bam
├── hcc1395.normal.rep3.trim.GRCh38_latest_genomic.hisat.bam.bai
├── hcc1395.normal.rep3.trim.GRCh38_latest_rna.bwa.bam
├── hcc1395.normal.rep3.trim.GRCh38_latest_rna.bwa.bam.bai
├── hcc1395.normal.rep3.trim.GRCh38_latest_rna.kallisto.abundance.h5
├── hcc1395.normal.rep3.trim.GRCh38_latest_rna.kallisto.abundance.tsv
├── hcc1395.normal.rep3.trim.GRCh38_latest_rna.kallisto.bam
├── hcc1395.normal.rep3.trim.GRCh38_latest_rna.kallisto.bam.bai
├── hcc1395.normal.rep3.trim.GRCh38_latest_rna.kallisto.info.json
├── hcc1395.normal.rep3.trim.r1.fastq.gz
├── hcc1395.normal.rep3.trim.r2.fastq.gz
├── hcc1395.normal.rep3.trim.unpair.r1.fastq.gz
├── hcc1395.normal.rep3.trim.unpair.r2.fastq.gz
├── hcc1395.tumor.rep1.fastqc.r1.html
├── hcc1395.tumor.rep1.fastqc.r1.zip
├── hcc1395.tumor.rep1.fastqc.r2.html
├── hcc1395.tumor.rep1.fastqc.r2.zip
├── hcc1395.tumor.rep1.r1.fastq.gz -> hcc1395_tumor_rep1_r1.fastq.gz
├── hcc1395_tumor_rep1_r1.fastq.gz
├── hcc1395.tumor.rep1.r2.fastq.gz -> hcc1395_tumor_rep1_r2.fastq.gz
├── hcc1395_tumor_rep1_r2.fastq.gz
├── hcc1395.tumor.rep1.trim.fastqc.r1.html
├── hcc1395.tumor.rep1.trim.fastqc.r1.zip
├── hcc1395.tumor.rep1.trim.fastqc.r2.html
├── hcc1395.tumor.rep1.trim.fastqc.r2.zip
├── hcc1395.tumor.rep1.trim.GRCh38_latest_genomic.hisat.bam
├── hcc1395.tumor.rep1.trim.GRCh38_latest_genomic.hisat.bam.bai
├── hcc1395.tumor.rep1.trim.GRCh38_latest_rna.bwa.bam
├── hcc1395.tumor.rep1.trim.GRCh38_latest_rna.bwa.bam.bai
├── hcc1395.tumor.rep1.trim.GRCh38_latest_rna.kallisto.abundance.h5
├── hcc1395.tumor.rep1.trim.GRCh38_latest_rna.kallisto.abundance.tsv
├── hcc1395.tumor.rep1.trim.GRCh38_latest_rna.kallisto.bam
├── hcc1395.tumor.rep1.trim.GRCh38_latest_rna.kallisto.bam.bai
├── hcc1395.tumor.rep1.trim.GRCh38_latest_rna.kallisto.info.json
├── hcc1395.tumor.rep1.trim.r1.fastq.gz
├── hcc1395.tumor.rep1.trim.r2.fastq.gz
├── hcc1395.tumor.rep1.trim.unpair.r1.fastq.gz
├── hcc1395.tumor.rep1.trim.unpair.r2.fastq.gz
├── hcc1395.tumor.rep2.fastqc.r1.html
├── hcc1395.tumor.rep2.fastqc.r1.zip
├── hcc1395.tumor.rep2.fastqc.r2.html
├── hcc1395.tumor.rep2.fastqc.r2.zip
├── hcc1395.tumor.rep2.r1.fastq.gz -> hcc1395_tumor_rep2_r1.fastq.gz
├── hcc1395_tumor_rep2_r1.fastq.gz
├── hcc1395.tumor.rep2.r2.fastq.gz -> hcc1395_tumor_rep2_r2.fastq.gz
├── hcc1395_tumor_rep2_r2.fastq.gz
├── hcc1395.tumor.rep2.trim.fastqc.r1.html
├── hcc1395.tumor.rep2.trim.fastqc.r1.zip
├── hcc1395.tumor.rep2.trim.fastqc.r2.html
├── hcc1395.tumor.rep2.trim.fastqc.r2.zip
├── hcc1395.tumor.rep2.trim.GRCh38_latest_genomic.hisat.bam
├── hcc1395.tumor.rep2.trim.GRCh38_latest_genomic.hisat.bam.bai
├── hcc1395.tumor.rep2.trim.GRCh38_latest_rna.bwa.bam
├── hcc1395.tumor.rep2.trim.GRCh38_latest_rna.bwa.bam.bai
├── hcc1395.tumor.rep2.trim.GRCh38_latest_rna.kallisto.abundance.h5
├── hcc1395.tumor.rep2.trim.GRCh38_latest_rna.kallisto.abundance.tsv
├── hcc1395.tumor.rep2.trim.GRCh38_latest_rna.kallisto.bam
├── hcc1395.tumor.rep2.trim.GRCh38_latest_rna.kallisto.bam.bai
├── hcc1395.tumor.rep2.trim.GRCh38_latest_rna.kallisto.info.json
├── hcc1395.tumor.rep2.trim.r1.fastq.gz
├── hcc1395.tumor.rep2.trim.r2.fastq.gz
├── hcc1395.tumor.rep2.trim.unpair.r1.fastq.gz
├── hcc1395.tumor.rep2.trim.unpair.r2.fastq.gz
├── hcc1395.tumor.rep3.fastqc.r1.html
├── hcc1395.tumor.rep3.fastqc.r1.zip
├── hcc1395.tumor.rep3.fastqc.r2.html
├── hcc1395.tumor.rep3.fastqc.r2.zip
├── hcc1395.tumor.rep3.r1.fastq.gz -> hcc1395_tumor_rep3_r1.fastq.gz
├── hcc1395_tumor_rep3_r1.fastq.gz
├── hcc1395.tumor.rep3.r2.fastq.gz -> hcc1395_tumor_rep3_r2.fastq.gz
├── hcc1395_tumor_rep3_r2.fastq.gz
├── hcc1395.tumor.rep3.trim.fastqc.r1.html
├── hcc1395.tumor.rep3.trim.fastqc.r1.zip
├── hcc1395.tumor.rep3.trim.fastqc.r2.html
├── hcc1395.tumor.rep3.trim.fastqc.r2.zip
├── hcc1395.tumor.rep3.trim.GRCh38_latest_genomic.hisat.bam
├── hcc1395.tumor.rep3.trim.GRCh38_latest_genomic.hisat.bam.bai
├── hcc1395.tumor.rep3.trim.GRCh38_latest_rna.bwa.bam
├── hcc1395.tumor.rep3.trim.GRCh38_latest_rna.bwa.bam.bai
├── hcc1395.tumor.rep3.trim.GRCh38_latest_rna.kallisto.abundance.h5
├── hcc1395.tumor.rep3.trim.GRCh38_latest_rna.kallisto.abundance.tsv
├── hcc1395.tumor.rep3.trim.GRCh38_latest_rna.kallisto.bam
├── hcc1395.tumor.rep3.trim.GRCh38_latest_rna.kallisto.bam.bai
├── hcc1395.tumor.rep3.trim.GRCh38_latest_rna.kallisto.info.json
├── hcc1395.tumor.rep3.trim.r1.fastq.gz
├── hcc1395.tumor.rep3.trim.r2.fastq.gz
├── hcc1395.tumor.rep3.trim.unpair.r1.fastq.gz
├── hcc1395.tumor.rep3.trim.unpair.r2.fastq.gz
├── NM_007297.4.fa
├── NM_007297.4.GRCh38_latest_rna.blast.csv
├── NM_007297.4.GRCh38_latest_rna.blast.txt
├── practical.tar
├── TruSeq3-PE.fa
├── week3_blast.py
└── week3_mapping.py
```
