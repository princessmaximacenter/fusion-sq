#!/usr/bin/env bash
root_dir=~/PycharmProjects/structuralvariation/
script_dir=${root_dir}fusion_sq/R/
utils_dir=${root_dir}utils/
resources_dir=${root_dir}resources/

## Make sure to set the paths in default.conf

Rscript ${script_dir}prepare_matching_intervals.R --patient_identifier test_patient1_PMABM000DKY
Rscript ${script_dir}run_match_wgs.R --patient_identifier test_patient1_PMABM000DKY --tool manta
Rscript ${script_dir}run_match_wgs.R --patient_identifier test_patient1_PMABM000DKY --tool delly 
Rscript ${script_dir}run_match_wgs.R --patient_identifier test_patient1_PMABM000DKY --tool gridss
Rscript ${script_dir}combine_wgs_support.R --patient_identifier test_patient1_PMABM000DKY
Rscript ${script_dir}collect_cohort.R

## retrieve resources first => see resources/get_resources.sh

Rscript ${utils_dir}pairwise_overlaps.population_svs.R

Rscript ${script_dir}annotate_supporting_svs.R
Rscript ${script_dir}annotate_cohort_report.R
Rscript ${script_dir}make_supplementary.R
