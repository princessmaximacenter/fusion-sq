#!/usr/bin/env bash
SINGULARITYENV_R_MAX_VSIZE=100Gb
singularity exec --bind /hpc/pmc_gen/ivanbelzen/ /hpc/pmc_gen/ivanbelzen/structuralvariation/structural_variation_latest.sif R -e "source('/hpc/pmc_gen/ivanbelzen/github_fusion_sq/fusion-sq/R/default.conf');source('/hpc/pmc_gen/ivanbelzen/github_fusion_sq/fusion-sq/R/hpc.default.conf');source(paste0(script_dir,'functions.get_vcf_filepath.docker.R'));source('/hpc/pmc_gen/ivanbelzen/structuralvariation/sv_functional_analysis/run/wilms_v2_20210923/wilms_v2_20210923.conf');source('/hpc/pmc_gen/ivanbelzen/structuralvariation/sv_functional_analysis/run/wilms_v2_20210923/wilms_v2_20210923.PMCID418AAA.conf');run_tool='gridss';source('/hpc/pmc_gen/ivanbelzen/github_fusion_sq/fusion-sq/R/run_match_wgs.docker.R');"

