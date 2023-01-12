#!/home/grahman/miniconda3/envs/evident-analyses/bin/python
#SBATCH --chdir=/home/grahman/projects/evident-analyses
#SBATCH --output=/home/grahman/projects/evident-analyses/log/supplemental/%x.log
#SBATCH --partition=short
#SBATCH --mem=8G
#SBATCH --time=6:00:00

import re

import biom
import numpy as np
import pandas as pd
import skbio

from src.helper import get_logger, time_function


@time_function
def main():
    logger = get_logger()

    tbl = biom.load_table("data/raw/supplemental/table.raw.ambig.biom")
    ids = [x for x in tbl.ids() if "BLANK" not in x]
    tbl.filter(ids)

    samp_rename_dict = {x: f"S{x}" for x in tbl.ids()}
    tbl.update_ids(samp_rename_dict)

    feat_rename_dict = {x: f"F{x}" for x in tbl.ids("observation")}
    tbl.update_ids(feat_rename_dict, "observation")

    na_vals = ["not provided", "not applicable"]
    orig_md_file = "data/raw/supplemental/metadata.raw.tsv"
    orig_md = pd.read_table(orig_md_file, sep="\t", index_col=0,
                            na_values=na_vals)
    orig_md = orig_md[~orig_md.index.str.contains("BLANK")]

    def determine_obese(row):
        if row["bmi_v2"] <= 25:
            return "No"
        elif row["bmi_v2"] >= 40:
            return "Yes"
        else:
            return np.nan
    orig_md["obese"] = orig_md.apply(determine_obese, axis=1)
    orig_md = orig_md.dropna()

    samp_reg = re.compile("(11666\.G\d+[A-Z])")
    def new_name(x):
        orig = samp_reg.search(x).groups()[0]
        return f"S{orig}"

    orig_md.index = [new_name(x) for x in orig_md.index]

    gg97_tax = pd.read_table("data/ref/97_otu_taxonomy.txt",
                             sep="\t", index_col=0, header=None)
    gg97_tax.columns = ["Taxon"]
    gg97_tax = gg97_tax["Taxon"].str.split( "; ", expand=True)
    gg97_tax.columns = list("kpcofgs")
    gg97_tax.index.name = "OTU_ID"
    gg97_tax.index = [f"F{x}" for x in gg97_tax.index]
    gg97_tax = gg97_tax.loc[tbl.ids("observation")]

    prevotella_feats = gg97_tax.query("g == 'g__Prevotella'").index
    bacteroides_feats = gg97_tax.query("g == 'g__Bacteroides'").index

    tbl_df = tbl.to_dataframe().T

    # Calculate log-ratio
    num_sum = tbl_df.loc[:, prevotella_feats].sum(axis=1)
    denom_sum = tbl_df.loc[:, bacteroides_feats].sum(axis=1)
    lr_df = pd.concat([num_sum, denom_sum], axis=1)
    lr_df.columns = ["num", "denom"]
    lr_df = lr_df.dropna(how="all")
    lr_df = lr_df + 1
    lr_df["log_ratio"] = np.log(lr_df["num"]/lr_df["denom"]).to_frame()
    lr_df = lr_df.join(orig_md, how="inner")

    outfile = "results/supplemental/prevotella_bacteroies_log_ratio.tsv"
    lr_df.to_csv(outfile, sep="\t", index=True)


if __name__ == "__main__":
    main()
