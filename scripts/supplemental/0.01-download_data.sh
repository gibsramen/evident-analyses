#!/usr/bin/bash
#SBATCH --chdir=/home/grahman/projects/evident-analyses
#SBATCH --output=/home/grahman/projects/evident-analyses/log/supplemental/%x.log
#SBATCH --partition=long
#SBATCH --mem=16G
#SBATCH --time=12:00:00

set -e
date; pwd; hostname

echo "[ $(date) ] :: Start script"

source ~/miniconda3/bin/activate evident-analyses

mkdir -p "data/raw/supplemental"

CTX="Pick_closed-reference_OTUs-Greengenes-Illumina-16S-V4-90nt-44feac"
TBL_FILE="data/raw/supplemental/table.raw.ambig.biom"
MD_FILE="data/raw/supplemental/metadata.raw.tsv"
SAMPLE_FILE="data/raw/supplemental/samples.txt"
QUERY1="(qiita_study_id == 11666)"
QUERY2="(host_age != 'not_provided' and sample_type == 'feces')"

echo "[ $(date) ] :: Start redbiom fetch samples"
redbiom search metadata "where $QUERY1 and $QUERY2" > $SAMPLE_FILE

echo "[ $(date) ] :: Start redbiom fetch table"
redbiom fetch samples \
    --context $CTX \
    --from $SAMPLE_FILE \
    --resolve-ambiguities most-reads \
    --output $TBL_FILE

echo "[ $(date) ] :: Start redbiom fetch metadata"
redbiom fetch sample-metadata \
    --context $CTX \
    --from $SAMPLE_FILE \
    --force-category "abdominal_obesity_ncep_v2" \
    --force-category "abdominal_obesity_idf_v2" \
    --force-category "bmi_v2" \
    --output $MD_FILE

echo "[ $(date) ] :: Finished script!"
