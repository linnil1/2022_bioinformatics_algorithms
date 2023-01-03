"""Final Exam"""
from pprint import pprint
from typing import Iterable
from pathlib import Path
from collections import defaultdict
import subprocess

from Bio import SeqIO, SeqRecord
from dash import Dash, dcc, html
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


def runShell(
    cmd: str, capture_output: bool = False, cwd: str | None = None
) -> subprocess.CompletedProcess[str]:
    """wrap os.system"""
    print(cmd)
    proc = subprocess.run(
        cmd,
        shell=True,
        capture_output=capture_output,
        cwd=cwd,
        check=True,
        universal_newlines=True,
    )
    return proc


def download(folder: str = "") -> str:
    """Download exam data"""
    if not folder:
        folder = "TWB2_and_GRCh37"
    if Path(folder).exists():
        return folder
    runShell("gdown --folder 'https://drive.google.com/drive/folders/1ljuYHB07vg3EUFbg_vmnblW6mVH5jSH5'")
    runShell("cd TWB2_and_GRCh37 && unzip ref_GRCh37.zip")
    if folder != "TWB2_and_GRCh37":
        runShell(f"mv TWB2_and_GRCh37 {folder}")
    return folder


def findRef(
    seqs: dict[str, SeqRecord], chrom: str, pos: int, length: int = 1, strand: str = ""
) -> str:
    """Query reference bases"""
    if strand == "-":
        seq = seqs[chrom].seq[pos - length : pos]
        return str(seq.reverse_complement())
    elif strand == "+":
        seq = seqs[chrom].seq[pos - 1 : pos - 1 + length]
        return str(seq)
    raise ValueError(f"invalid strand {strand}")


def countViolateSNV(seqs: dict[str, SeqRecord], twb2: pd.DataFrame) -> list[go.Figure]:
    """Q1"""
    count: dict[str, int] = defaultdict(int)
    for variant in twb2.itertuples():
        if variant.ref == "-":  # note1
            count["-"] += 1
            continue

        strand = variant.strand
        if variant.strand == "---":  # note2.b  # TODO: Answer: It's +
            strand = "+"

        if not all(i in "ATCG" for i in variant.ref):
            raise ValueError(f"ref contains strange alphabet: {variant.ref}")

        hg19_ref = findRef(
            seqs,
            variant.chrom,
            variant.pos,
            length=len(variant.ref),
            strand=strand,
        )
        if variant.ref != hg19_ref:
            # print(variant)
            # print(variant.ref, hg19_ref)
            count["violate"] += 1
        else:
            count["ok"] += 1
    pprint(count)
    return []


def countSNVTrend(twb2: pd.DataFrame) -> list[go.Figure]:
    """Q2"""
    nc_change = twb2.groupby(["ref", "alt"], as_index=False).size()
    nc_change = nc_change[~nc_change["ref"].str.contains("-")]  # type: ignore
    nc_change = nc_change[~nc_change["alt"].str.contains("-")]  # type: ignore
    nc_change = nc_change[ nc_change["ref"].str.len() == 1]  # type: ignore
    nc_change = nc_change[ nc_change["alt"].str.len() == 1]  # type: ignore
    print(nc_change)

    # TODO: strand? -> ANS: Ignore strand
    nc_change["variant"] = nc_change["ref"] + "->" + nc_change["alt"]  # type: ignore

    # plot
    fig1 = px.bar(nc_change, x="variant", y="size", text="size")
    fig1.update_layout(
        title="Q2-a",
        barmode="stack",
        yaxis_title="variant count",
        xaxis_categoryorder="total descending",
    )

    fig2 = px.bar(nc_change, x="variant", y="size", text="size")
    fig2.update_layout(
        title="Q2-b",
        yaxis_title="variant count",
        xaxis_categoryorder="category ascending",
    )
    return [fig1, fig2]


def countIndelLength(twb2: pd.DataFrame) -> list[go.Figure]:
    """Q3"""
    indel = twb2[
        np.logical_or.reduce(
            [
                twb2["ref"].str.contains("-"),
                twb2["alt"].str.contains("-"),
                # twb2["ref"].str.len() != 1,
                # twb2["alt"].str.len() != 1,
            ]
        )
    ]  # TODO: when ref and alt are not len=1  -> Answer: not indel, not snv
    indel["ref_len"] = indel["ref"].str.len()
    indel["alt_len"] = indel["alt"].str.len()
    indel["indel_len"] = indel[["ref_len", "alt_len"]].max(1)
    indel.loc[twb2["ref"].str.contains("-"), "type"] = "insertion"
    indel.loc[twb2["alt"].str.contains("-"), "type"] = "deletion"

    indel_size = indel.groupby(["indel_len", "type"], as_index=False).size()
    with pd.option_context("display.max_rows", None):
        print(indel_size)

    fig = px.bar(indel_size, x="indel_len", y="size", color="type", text="size")
    fig.update_layout(
        title="Q3",
        yaxis_title="variant count",
        xaxis_categoryorder="total descending",
        xaxis_range=(0, 50),  # most of them are in 0 to 50, so plot here
    )
    return [fig]


def showPlot(figs: Iterable[go.Figure]) -> None:
    """Show all the figure in website default localhost:8051"""
    app = Dash(__name__)
    app.layout = html.Div([dcc.Graph(figure=f) for f in figs])
    app.run_server(debug=True, port=8051)


def main(folder: str) -> None:
    """Run all 3 questions"""
    # read
    folder = download(folder)
    seqs = SeqIO.to_dict(
        SeqIO.parse(folder + "/ref_GRCh37/GRCh37.p13.genome.fa", "fasta")
    )
    twb2 = pd.read_csv(
        folder + "/TWB2_adjusted.csv",
        usecols=[
            "Ref Allele",
            "Alt Allele",
            "hg19 chromosome",
            "hg19 position",
            "hg19 strand",
        ],
    )
    twb2 = twb2.set_axis(["ref", "alt", "chrom", "pos", "strand"], axis=1)
    # print(seqs)
    # print(twb2)

    # main
    figs = []
    # 13:00 -
    countViolateSNV(seqs, twb2)
    # 14:00
    # 14:45 -
    figs.extend(countSNVTrend(twb2))
    figs.extend(countIndelLength(twb2))
    showPlot(figs)
    # 15:30


if __name__ == "__main__":
    # main("test")
    main("")
    # 16:27 Done
