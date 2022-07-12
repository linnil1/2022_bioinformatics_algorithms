import os
import gdown
import pandas as pd
import itertools
from collections import Counter
import plotly.express as px


data_source = "https://drive.google.com/file/d/1_p6SpoOPusGkKaG4sbyVoypvQFY0kDDU/view?usp=sharing"  # noqa
data_file = "HW2_clinvar.vcf.gz"
output_file1 = "hw2_clnsig.png"
output_file2 = "hw2_rs.tsv"


# main
if __name__ == "__main__":
    if not os.path.exists(data_file):
        print("Download vcf")
        gdown.download(data_source, data_file, fuzzy=True)

    print(f"read {data_file}")
    df = pd.read_csv(data_file, sep="\t", compression='gzip', comment="#", header=None)
    df.columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]

    print("INFO to table")
    # ALLELEID=1003021;CLNDISDB=MedGen:CN517202;CLNDN=not_provided
    # ->
    #   ALLELEID CLNDISDB        CLNDN  ...
    # 0 1003021  MedGen:CN517202 not_provided  ...
    # 1 1632777  MedGen:CN517202 not_provided  ...
    df_info = df["INFO"].apply(lambda i: dict(map(lambda j: j.split("="), i.split(";"))))
    # series of dict -> dataframe
    df_info = pd.DataFrame.from_records(df_info)

    print(f"Counting CLNSIG {output_file1}")
    # Conflicting_interpretations_of_pathogenicity|_other
    clnsig_list = df_info["CLNSIG"].dropna().str.split("|_", regex=False)
    # list of list of clnsig -> dict[type, count]
    count = Counter(itertools.chain.from_iterable(clnsig_list))
    df_count = pd.DataFrame(sorted(count.items()))
    df_count.columns = ["CLNSIG", "count"]
    print(df_count)

    print(f"Plot CLNSIG into {output_file1}")
    fig = px.bar(df_count, x="CLNSIG", y="count",
                 title=f"CLNSIG count in {output_file1}", text_auto="d")
    fig.write_image(output_file1)

    print("Merge RS to each record {output_file2}")
    out = df[["CHROM", "POS", "ID", "REF", "ALT"]]
    out["RS"] = df_info["RS"].fillna(".")
    out.to_csv(output_file2, index=False, sep="\t")
    print(out)
