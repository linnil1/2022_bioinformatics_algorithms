# Week4 Metagenomics

Mostly follow https://github.com/vappiah/bacterial-genomics-tutorial

## Run it
```
cd week4
python week4.py
```

## Input:

* `week4/P7741.{r1,r2}.fasta.gz`


## Output
```
data
├── P7741.fastqc.r1.html
├── P7741.fastqc.r1.zip
├── P7741.fastqc.r2.html
├── P7741.fastqc.r2.zip
├── P7741.r1.fastq.gz -> P7741_R1.fastq.gz
├── P7741_R1.fastq.gz
├── P7741.r2.fastq.gz -> P7741_R2.fastq.gz
├── P7741_R2.fastq.gz
├── P7741.trim.fastqc.r1.html
├── P7741.trim.fastqc.r1.zip
├── P7741.trim.fastqc.r2.html
├── P7741.trim.fastqc.r2.zip
├── P7741.trim.r1.fastq.gz
├── P7741.trim.r2.fastq.gz
├── P7741.trim.s.fastq.gz
├── P7741.trim.spades
│   ├── assembly_graph_after_simplification.gfa
│   ├── assembly_graph.fastg
│   ├── assembly_graph_with_scaffolds.gfa
│   ├── before_rr.fasta
│   ├── contigs.fasta
│   ├── contigs.paths
│   ├── corrected
│   ├── dataset.info
│   ├── input_dataset.yaml
│   ├── K127
│   ├── K21
│   ├── K33
│   ├── K55
│   ├── K77
│   ├── K99
│   ├── misc
│   ├── mismatch_corrector
│   ├── params.txt
│   ├── pipeline_state
│   ├── run_spades.sh
│   ├── run_spades.yaml
│   ├── scaffolds.fasta
│   ├── scaffolds.paths
│   ├── spades.log
│   ├── tmp
│   └── warnings.log
├── P7741.trim.spades.contig.fasta -> P7741.trim.spades/contigs.fasta
├── P7741.trim.spades.contig.fasta.fai
├── P7741.trim.spades.contig.reorder
│   ├── ragtag.scaffold.agp
│   ├── ragtag.scaffold.asm.paf
│   ├── ragtag.scaffold.asm.paf.log
│   ├── ragtag.scaffold.confidence.txt
│   ├── ragtag.scaffold.err
│   ├── ragtag.scaffold.fasta
│   └── ragtag.scaffold.stats
├── P7741.trim.spades.contig.reorder.amr.tab
├── P7741.trim.spades.contig.reorder.annotate
│   ├── P7741.err
│   ├── P7741.faa
│   ├── P7741.ffn
│   ├── P7741.fna
│   ├── P7741.fsa
│   ├── P7741.gbk
│   ├── P7741.gff
│   ├── P7741.log
│   ├── P7741.pseudo.txt
│   ├── P7741.sqn
│   ├── P7741.tbl
│   ├── P7741.tsv
│   └── P7741.txt
├── P7741.trim.spades.contig.reorder.annotate.faa -> P7741.trim.spades.contig.reorder.annotate/P7741.faa
├── P7741.trim.spades.contig.reorder.dendogram
│   ├── data
│   ├── data_tables
│   ├── dereplicated_genomes
│   ├── figures
│   └── log
└── P7741.trim.spades.contig.reorder.fasta
```
