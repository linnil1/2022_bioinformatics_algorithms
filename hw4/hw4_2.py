import os
import sys
import subprocess
import numpy as np
import pandas as pd
import plotly.express as px
from sklearn.decomposition import PCA


def runShell(cmd):
    print(cmd)
    proc = subprocess.run(cmd, shell=True, executable="/bin/bash")
    proc.check_returncode()


def download():
    if os.path.exists("Tutorial"):
        return "Tutorial/Example"
    runShell("git clone https://github.com/metagenome-atlas/Tutorial.git")
    return "Tutorial/Example"


def getMappingRate():
    # read
    raw_counts = pd.read_csv("Tutorial/Example/genomes/counts/raw_counts_genomes.tsv",
                             index_col=0, sep='\t')
    read_stats = pd.read_csv("Tutorial/Example/stats/read_counts.tsv",
                             index_col=0, sep='\t')
    # filterd by QC
    # read_stats = read_stats[read_stats["Step"] == "QC"]

    # calculate mapping rate
    mapped = raw_counts.T.sum(axis=1)
    total = read_stats["Reads_pe"] * 2 + read_stats["Reads_se"]
    mapping_rate = mapped / total
    print("mapping rate", file=sys.stderr)
    print(mapping_rate, file=sys.stderr)
    return mapping_rate


def getAbundance():
    # read
    abund = pd.read_table("Tutorial/Example/genomes/counts/median_coverage_genomes.tsv",
                          index_col=0)
    taxonomy = pd.read_table("Tutorial/Example/genomes/taxonomy/gtdb_taxonomy.tsv",
                             index_col=0)
    # fill species == NA
    species_empty_index = taxonomy["species"].isnull()
    taxonomy.loc[species_empty_index, "species"] = taxonomy.index[species_empty_index]

    # MAG abundance
    abund_norm = (abund.T / abund.sum(axis=1)).T
    mag_abund = abund_norm.mean(axis=0).to_frame()
    mag_abund = mag_abund.merge(taxonomy, left_index=True, right_index=True)

    # group by species
    species_abund = mag_abund.groupby("species")[0].sum().sort_values(ascending=False)
    species_abund = species_abund.reset_index().rename(columns={0: "abundance"})
    print("abundacne", file=sys.stderr)
    print(species_abund, file=sys.stderr)
    return species_abund


def readFunctions():
    """ Read functions data and merge taxonomy """
    kegg_modules = pd.read_csv("Tutorial/Example/genomes/annotations/dram/kegg_modules.tsv",  # noqa
                               index_col=1, sep='\t')
    taxonomy = pd.read_table("Tutorial/Example/genomes/taxonomy/gtdb_taxonomy.tsv",
                             index_col=0)
    kegg_modules = kegg_modules.merge(taxonomy, how="left",
                                      left_index=True, right_index=True)
    return kegg_modules


def getMaxFunctions():
    # read
    kegg_modules = readFunctions()
    module_name = kegg_modules.groupby("module")["module_name"].first()

    # group by module
    module_avg_coverage = kegg_modules.groupby("module")["step_coverage"].mean()
    module_avg_coverage = pd.merge(module_name, module_avg_coverage,
                                   left_index=True, right_index=True)

    print(module_avg_coverage.sort_values("step_coverage", ascending=False),
          file=sys.stderr)
    max_value = module_avg_coverage["step_coverage"].max()
    max_module = module_avg_coverage[module_avg_coverage["step_coverage"] == max_value]
    return max_module


def getPhylumFunctions():
    # read and filter
    kegg_modules = readFunctions()[["module", "step_coverage", "phylum"]]
    kegg_modules_filtered = kegg_modules[kegg_modules["step_coverage"] > 0.8]

    # group by phylum
    phylum_count = kegg_modules.reset_index().groupby("index").first()\
                               .groupby("phylum").size()
    reserved_count = kegg_modules_filtered.groupby("phylum").size()
    # print(reserved_count, file=sys.stderr)
    # print(phylum_count, file=sys.stderr)
    return reserved_count / phylum_count


def functionsPCA():
    # read
    kegg_modules = readFunctions()
    phylums = ["Firmicutes_A", "Firmicutes"]
    kegg_modules = kegg_modules[kegg_modules["phylum"].isin(phylums)]

    # reshape to key=MAG value=all module's step_coverage
    kegg_modules = kegg_modules.reset_index().pivot(index=["index", "phylum"],
                                                    columns=["module"],
                                                    values=["step_coverage"])
    print(kegg_modules, file=sys.stderr) 

    # PCA
    pca = PCA(n_components=2)
    pca_data = pca.fit_transform(np.array(kegg_modules))
    genomes, phylum = list(zip(*kegg_modules.index))
    pca_kegg_modules = pd.DataFrame(zip(genomes, phylum, pca_data[:, 0], pca_data[:, 1]),
                                    columns=["genomes", "phylum", "x", "y"])
    print(pca_kegg_modules, file=sys.stderr)

    # plot
    fig = px.scatter(pca_kegg_modules, x="x", y="y", color="phylum", text="genomes",
                     title="PCA of step_coverage")
    fig.write_image("hw4.2.png")
    return "hw4.2.png"


if __name__ == "__main__":
    download()
    print("""
Hint
# start your jupyter notebook in docker
dk -p 8888:8888 quay.io/biocontainers/metagenome-atlas:2.11.0--pyhdfd78af_0  bash
mamba env create  -f Tutorial/Python/condaenv.yml -p conda_env
conda init bash
jupyter lab --allow-root --ip 0.0.0.0
    """)

    print(f"""
Q1
Quantification is based on mapping the reads to the genomes using bbmap.
Do you think the mapping rate is good? what could be done to improve it?
  Mean mapping rate: {np.mean(getMappingRate())}
  Altas map reads to contigs by bbmap written in `build_db_genomes` `align_reads_to_MAGs
  Maybe try BWA or Bowtie2 or HISAT2 for better performance (for short reads)
    """)

    print(f"""
Q2
We use the median coverage over the genomes to calculate the relative abundance.
What is the most abundant species in these metagenomes?
  Just group the relative abundance by species.
  And the species with maximum abundance is
{getAbundance().loc[0]}
    """)

    print(f"""
Q3
What are conserved functions among the genomes?
  I assume step_coverage higher -> more conserved
{getMaxFunctions()}
    """)

    print(f"""
Q4
Do you see trends (of conserved functions) by phylum?
  I assume conserved functions is the module that step_coverage > 0.8
  Then, I count the number of conserved functions group by phylum
{getPhylumFunctions()}
    """)

    print(f"""
Q5
Can you distinguish genomes from the Firmicutes_A and Firmicutes phyla?
  I run PCA on the step_coverage as MAG's features.
  See PCA result in {functionsPCA()}
  Most of Firmicutes can sepatated from Firmicutes_A, except MAG17.
    """)
