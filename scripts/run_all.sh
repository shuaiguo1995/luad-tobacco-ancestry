#!/bin/sh
# Runs the analysis scripts in order and stops at the first failure.
#   sh run_all.sh
#   RSCRIPT=/path/to/Rscript sh run_all.sh
# 00_harmonise_cohorts.Rmd is not run here: it needs the raw downloads in
# ../raw/, which are not distributed, and its outputs are already in ../data/.
# The order matters: 03 reads a table written by 01, and 06 reads tables
# written by 01 and 05.
cd "$(dirname "$0")" || exit 1
RSCRIPT=${RSCRIPT:-Rscript}
for f in 01_four_groups.Rmd 02_within_cohort.Rmd 03_pooled_vs_stratified.Rmd \
         04_cohort_pairs.Rmd 05_sensitivity.Rmd 06_supplementary_tables.Rmd
do
  echo "########## $f"
  "$RSCRIPT" run.R "$f" || { echo "########## $f failed" >&2; exit 1; }
done
echo "########## all scripts finished"
