#!/usr/bin/env bash

test_identifier="fusion_sq"

wdir="/hpc/pmc_gen/ivanbelzen/github_fusion_sq/fusion-sq/"

cohort_wdir="${wdir}run/${test_identifier}/"
cohort_file="${cohort_wdir}cohort.${test_identifier}.tsv"
cohort_file="/hpc/pmc_gen/ivanbelzen/metadata/cohort.v9.metadata.20210726.tsv"

config_template="${wdir}patient.conf"

while read -r patient_id normal_id tumor_id rna_id rest; do
#while read -r patient_id tumor_id normal_id rna_id rest; do
      if [[ ${patient_id} == "patient_id" ]];  then continue; fi
      if [[ ${patient_id} == "PMC_ID" ]];  then continue; fi

      config_patient="${cohort_wdir}${test_identifier}.${patient_id}.conf"
      sed -e "s/\[patient_id]/${patient_id}/" -e "s/\[normal_id]/${normal_id}/"  -e "s/\[tumor_id]/${tumor_id}/" -e "s/\[rna_id]/${rna_id}/" ${config_template}  > ${config_patient}

done < $cohort_file
