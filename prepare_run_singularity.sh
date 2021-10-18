#!/usr/bin/env bash
## Last update 2021-10-11
## Run fusion sq

test_identifier=$1
singularity_img="/hpc/pmc_gen/ivanbelzen/structuralvariation/structural_variation_latest.sif"

fusionsq_dir="/hpc/pmc_gen/ivanbelzen/github_fusion_sq/fusion-sq/R/"
sv_wdir="/hpc/pmc_gen/ivanbelzen/structuralvariation/sv_functional_analysis/"

config_default="${fusionsq_dir}default.conf"
config_docker="${fusionsq_dir}hpc.default.conf"
cohort_dir="${sv_wdir}run/${test_identifier}/"
cohort_file="${cohort_dir}cohort.${test_identifier}.tsv"
cohort_run_file=${cohort_file}

#run subset or full cohort
cohort_run_file="${cohort_dir}run_cohort.wilms_v2_20211002.tsv"
config_general="${cohort_dir}${test_identifier}.conf"

output_singularity_dir="${cohort_dir}output/"

input_vcf_data_dir="/hpc/pmc_gen/ivanbelzen/structuralvariation/tmp/"
host_server="gerrit"
source_dir="/data/groups/gen/ivanbelzen/"
gridss_dir="${source_dir}data/gridss-2.7.2/blklst_dac_decoy/"
manta_dir="${source_dir}data/manta-1.6.0/"
delly_dir="${source_dir}data/delly-0.8.1/blklst_dac_decoy/"
archive_dir="${source_dir}structuralvariation/v1/${test_identifier}/"
#final_output_file_ext=".somatic.svs_merged.tsv"

TMPDIR="${input_vcf_data_dir}"

prepare_script="${fusionsq_dir}prepare_matching_intervals.docker.R"
match_wgs_script="${fusionsq_dir}run_match_wgs.docker.R"


while read -r patient_id tumor_id normal_id rna_id <&3; do
      if [[ ${patient_id} == "patient_id" ]];  then continue; fi
      if [[ ${patient_id} == "PMC_ID" ]];  then continue; fi
      if [[ ${patient_id} == "" ]];  then continue; fi

      basename="${tumor_id}_${normal_id}_WGS"
      config_patient="${cohort_dir}${test_identifier}.${patient_id}.conf"
      runscript="${cohort_dir}run_fusion_sq.${test_identifier}.${patient_id}.sh"
      singularity_settings="source('${config_default}');source('${config_docker}');source('${config_general}');source('${config_patient}');"

      #if  ls ${archive_dir}${tumor_id}_${normal_id}${final_output_file_ext} 1> /dev/null 2>&1; then
      #continue
      #fi

    echo "#!/usr/bin/env bash" > ${runscript}
    echo "export TMPDIR=${TMPDIR}" >> ${runscript}
    
    echo "SINGULARITYENV_R_MAX_VSIZE=100Gb" >> ${runscript}
    echo "singularity exec --bind /hpc/pmc_gen/ivanbelzen/ ${singularity_img} R -e \"${singularity_settings}source('${prepare_script}')\"" >> ${runscript}
    echo "singularity exec --bind /hpc/pmc_gen/ivanbelzen/ ${singularity_img} R -e \"${singularity_settings}run_tool='manta';source(paste0(script_dir,'functions.get_vcf_filepath.docker.R'));source('${match_wgs_script}')\"" >> ${runscript}
    echo "singularity exec --bind /hpc/pmc_gen/ivanbelzen/ ${singularity_img} R -e \"${singularity_settings}run_tool='delly';source(paste0(script_dir,'functions.get_vcf_filepath.docker.R'));source('${match_wgs_script}')\"" >> ${runscript}
    echo "singularity exec --bind /hpc/pmc_gen/ivanbelzen/ ${singularity_img} R -e \"${singularity_settings}run_tool='gridss';source(paste0(script_dir,'functions.get_vcf_filepath.docker.R'));source('${match_wgs_script}')\"" >> ${runscript}

    echo "#rsync -avhe ssh ${output_singularity_dir}${tumor_id}_${normal_id}.* ${host_server}:${archive_dir}" >> ${runscript}
    
    sbatch --mem=120Gb --time=5:00:00 ${runscript}

done 3<$cohort_run_file
