#!/bin/bash
#SBATCH -n 72
#SBATCH -N 1
#SBATCH --mem 299G
#SBATCH --time 48:0:0
#SBATCH --qos castles
#SBATCH --mail-type ALL
#SBATCH --account=croftap-stia-atac


export MRO_DISK_SPACE_CHECK=disable

#set -e
#module purge; module load bluebear


cd /rds/projects/c/croftap-celldive01/CosmX_analysis/Run2_boston/analysis_baysor/step1_remove_controls/S0

for f in *.csv; do /rds/projects/c/croftap-celldive01/CosmX_analysis/baysor0.6.0/baysor/bin/baysor run -s 50 -x x -y y -g target -m 10 --save-polygons=GeoJSON --prior-segmentation-confidence=0.7 --n-clusters=10 -o /rds/projects/c/croftap-celldive01/CosmX_analysis/Run2_boston/analysis_baysor/step2_runbaysor/S0/"$f" "$f" :CellId; done

cd /rds/projects/c/croftap-celldive01/CosmX_analysis/Run2_boston/analysis_baysor/step1_remove_controls/S1

for f in *.csv; do /rds/projects/c/croftap-celldive01/CosmX_analysis/baysor0.6.0/baysor/bin/baysor run -s 50 -x x -y y -g target -m 10 --save-polygons=GeoJSON --prior-segmentation-confidence=0.7 --n-clusters=10 -o /rds/projects/c/croftap-celldive01/CosmX_analysis/Run2_boston/analysis_baysor/step2_runbaysor/S1/"$f" "$f" :CellId; done
