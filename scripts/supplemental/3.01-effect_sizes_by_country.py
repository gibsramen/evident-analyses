#!/home/grahman/miniconda3/envs/evident-analyses/bin/python
#SBATCH --chdir=/home/grahman/projects/evident-analyses
#SBATCH --output=/home/grahman/projects/evident-analyses/log/supplemental/%x.log
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16G
#SBATCH --time=6:00:00

import evident
from evident.effect_size import effect_size_by_category
import numpy as np
import pandas as pd
from skbio import DistanceMatrix

from src.helper import get_logger, time_function


@time_function
def main():
    logger = get_logger()

    md_file = "data/processed/metadata.disambig.filt.tsv"
    md = pd.read_table(md_file, sep="\t", index_col=0)
    md.index = md.index.astype(str)

    countries = ["United States", "United Kingdom"]
    md = md[md["country_of_birth"].isin(countries)]

    logger.info("Loading distance matrix...")
    dm = DistanceMatrix.read(
        "results/distance-matrix-u_unifrac.tsv"
    )
    logger.info("Finished loading distance matrix!")

    for country, country_md in md.groupby("country_of_birth"):
        logger.info(f"Country: {country}")
        bdh = evident.MultivariateDataHandler(dm, country_md)

        logger.info("Calculating effect sizes...")
        res = effect_size_by_category(
            bdh,
            bdh.metadata.columns,
            n_jobs=8,
            parallel_args={
                "backend": "multiprocessing",
                "verbose": 100
            }
        )
        res_df = res.to_dataframe()
        logger.info(f"\n{res_df}")
        country_str = country.replace(" ", "_")
        res_df.to_csv(
            f"results/supplemental/beta_effect_size_by_cat.{country_str}.tsv",
            sep="\t"
        )
        logger.info("Finished calculating effect sizes!")


if __name__ == "__main__":
    main()
