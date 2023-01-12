#!/home/grahman/miniconda3/envs/evident-analyses/bin/python
#SBATCH --chdir=/home/grahman/projects/evident-analyses
#SBATCH --output=/home/grahman/projects/evident-analyses/log/supplemental/%x.log
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16G
#SBATCH --time=6:00:00

import evident
import numpy as np
import pandas as pd

from src.helper import get_logger, time_function


@time_function
def main():
    logger = get_logger()

    df_file = "results/supplemental/prevotella_bacteroies_log_ratio.tsv"
    df = pd.read_table(df_file, sep="\t", index_col=0)

    udh = evident.UnivariateDataHandler(df["log_ratio"], df)

    difference = [0.5, 1.0, 1.5]
    total_observations = np.arange(100, 1001, step=100)

    results = udh.power_analysis(
        "obese",
        alpha=0.05,
        difference=difference,
        total_observations=total_observations
    ).to_dataframe()

    logger.info(f"\n{results.head()}")
    results.to_csv("results/supplemental/power_analysis.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
