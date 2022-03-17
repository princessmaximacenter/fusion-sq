#!/usr/bin/env bash
## Last update 2022-01-22
## Run fusion sq 

test_identifier="fusion_sq"
singularity_img="/hpc/pmc_gen/ivanbelzen/structuralvariation/structural_variation_latest.sif"

wdir="/hpc/pmc_gen/ivanbelzen/github_fusion_sq/fusion-sq/"

#script dir
fusionsq_dir="${wdir}R/"
config_default="${fusionsq_dir}default.conf"
config_docker="${fusionsq_dir}default.docker.conf"

cohort_dir="${wdir}run/${test_identifier}/"
cohort_file="${cohort_dir}cohort.${test_identifier}.tsv"
cohort_file="/hpc/pmc_gen/ivanbelzen/metadata/cohort.v9.metadata.20210726.tsv"

#run subset or full cohort
cohort_run_file="${cohort_file}"

#cohort_run_file="${cohort_dir}/run_cohort.finished_fusioncatcher.tsv"
#batch="/hpc/pmc_gen/ivanbelzen/fusion_catcher/finished_patients.lst"

#cohort_run_file="${cohort_dir}/run_cohort.last_two.tsv"
#batch="${wdir}last_two.lst"

#grep -f $batch $cohort_file > $cohort_run_file


config_general="${cohort_dir}${test_identifier}.conf"

output_singularity_dir="${cohort_dir}output/"

#retrieve vcfs first - not part of script now
input_vcf_data_dir="/hpc/pmc_gen/ivanbelzen/structuralvariation/tmp/"
host_server="gerrit"
source_dir="/data/groups/gen/ivanbelzen/"

TMPDIR="${input_vcf_data_dir}"

prepare_script="${fusionsq_dir}prepare_matching_intervals.docker.R"
match_wgs_script="${fusionsq_dir}run_match_wgs.docker.R"
combine_wgs_script="${fusionsq_dir}combine_wgs_support.docker.R"

analysis_type="starfusion"
runscript_prefix="run_fusion_sq.starfusion."

while read -r patient_id tumor_id normal_id rna_id; do
      #if [[ ${patient_id} == "PMCID827AAD" ]];  then continue; fi
      if [[ ${patient_id} == "patient_id" ]];  then continue; fi
      if [[ ${patient_id} == "PMC_ID" ]];  then continue; fi
      if [[ ${patient_id} == "" ]];  then continue; fi

      basename="${tumor_id}_${normal_id}_WGS"
      config_patient="${cohort_dir}${test_identifier}.${patient_id}.conf"
      runscript="${cohort_dir}${runscript_prefix}${analysis_type}.${patient_id}.sh"
      singularity_settings="source('${config_default}');source('${config_docker}');source('${config_general}');source('${config_patient}');"
      singularity_settings="${singularity_settings}analysis_type='${analysis_type}';"
    
    echo "#!/usr/bin/env bash" > ${runscript}
    echo "export TMPDIR=${TMPDIR}" >> ${runscript}

    echo "SINGULARITYENV_R_MAX_VSIZE=100Gb" >> ${runscript}
    echo "#singularity exec --bind /hpc/pmc_gen/ivanbelzen/ ${singularity_img} R -e \"${singularity_settings}source('${prepare_script}')\"" >> ${runscript}
    echo "singularity exec --bind /hpc/pmc_gen/ivanbelzen/ ${singularity_img} R -e \"${singularity_settings}run_tool='manta';source(paste0(script_dir,'functions.get_vcf_filepath.docker.R'));source('${match_wgs_script}')\"" >> ${runscript}
    echo "singularity exec --bind /hpc/pmc_gen/ivanbelzen/ ${singularity_img} R -e \"${singularity_settings}run_tool='delly';source(paste0(script_dir,'functions.get_vcf_filepath.docker.R'));source('${match_wgs_script}')\"" >> ${runscript}
    echo "singularity exec --bind /hpc/pmc_gen/ivanbelzen/ ${singularity_img} R -e \"${singularity_settings}run_tool='gridss';source(paste0(script_dir,'functions.get_vcf_filepath.docker.R'));source('${match_wgs_script}')\"" >> ${runscript}

    echo "singularity exec --bind /hpc/pmc_gen/ivanbelzen/ ${singularity_img} R -e \"${singularity_settings}source('${combine_wgs_script}')\"" >> ${runscript}
    
    sbatch --mem=120Gb --time=25:00:00 ${runscript}
  

done < $cohort_run_file

