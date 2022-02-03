## Fusion sq
## Combine RNA tools
## Last update: 2022-02-03
## 

if(FALSE) {
  
  #local
  
  source("~/fusion_sq/R/default.conf")
  source("~/fusion_sq/R/default.docker.local.conf")
  source("~/fusion_sq/run/fusion_sq/fusion_sq.conf")
  
  #cat /hpc/pmc_gen/ivanbelzen/github_fusion_sq/fusion-sq/run/fusion_sq/fusion_sq.PMCID203AAL.conf 
  patient = c()
  patient$patient_id = "PMCID203AAL"
  patient$tumor_id = "PMABM000HAS"
  patient$normal_id = "PMABM000HBB"
  patient$basename = paste0(patient$tumor_id,"_",patient$normal_id)
  patient$tumor_label = paste0(patient$tumor_id,"_WGS")
  patient$normal_label = paste0(patient$normal_id,"_WGS")
  patient$rna_id="PMABM000HAT"
  patient$patient_identifier= paste0(patient$patient_id,"_",patient$rna_id)
  
  analysis_type="fusioncatcher"
  analysis_type="starfusion"
  
}
if(FALSE){
  #set paths for hpc
  source("/hpc/pmc_gen/ivanbelzen/github_fusion_sq/fusion-sq/R/default.conf")
  ## HPC config overrides
  source("/hpc/pmc_gen/ivanbelzen/github_fusion_sq/fusion-sq/R/default.docker.conf")
  #HPC doesnt use argparser but patient specific config instead 
  source("/hpc/pmc_gen/ivanbelzen/github_fusion_sq/fusion-sq/run/fusion_sq/fusion_sq.conf")
  source("/hpc/pmc_gen/ivanbelzen/github_fusion_sq/fusion-sq/run/fusion_sq/fusion_sq.PMCID144AAK.conf")
}

suppressPackageStartupMessages({
  library(tidyverse, quietly=TRUE)
  library(stringr, quietly=TRUE)
  library(stringdist, quietly=TRUE)
  #library(argparser, quietly=TRUE)
  library(GenomicRanges, quietly=TRUE)
  library(dplyr)
  library(stringi)
  library(openssl)
  
})

source(paste0(script_dir,"functions.general.R"))
source(paste0(script_dir,"functions.svs.R"))


if(patient$patient_id =="") {
  print("patient$patient_id needs to be specified")
  #quit()
}


#order of arguments matters
map_template_vars=c('${input_dir}'=input_dir,'${output_dir}'=output_dir,'${cohort_identifier}'=cohort_identifier,'${cohort_wdir}'=cohort_wdir,'${patient_identifier}'=patient$patient_identifier)

reports_dir =  stri_replace_all_fixed(reports_dir_template,names(map_template_vars), map_template_vars,vectorize=F)

map_template_vars = c(map_template_vars,'${reports_dir}'=reports_dir)

#Paths 
#fc_fusion_level_results_path = "~/fusion_sq/run/fusion_sq/output/fusion_level_results.fusioncatcher.PMCID203AAL_PMABM000HAT.tsv"
#sf_fusion_level_results_path = "~/fusion_sq/run/fusion_sq/output/fusion_level_results.starfusion.PMCID203AAL_PMABM000HAT.tsv"
#merged_fusion_level_results_path = "~/fusion_sq/run/fusion_sq/output/fusion_level_results.merged.PMCID203AAL_PMABM000HAT.tsv"

map_template_vars_analysis = c(map_template_vars,'${analysis_type}'="starfusion.")
sf_fusion_level_results_path = stri_replace_all_fixed(fusion_level_results_path_template,names(map_template_vars_analysis), map_template_vars_analysis,vectorize=F)

map_template_vars_analysis = c(map_template_vars,'${analysis_type}'="fusioncatcher.")
fc_fusion_level_results_path = stri_replace_all_fixed(fusion_level_results_path_template,names(map_template_vars_analysis), map_template_vars_analysis,vectorize=F)

map_template_vars_analysis = c(map_template_vars,'${analysis_type}'="merged.")
merged_fusion_level_results_path = stri_replace_all_fixed(fusion_level_results_path_template,names(map_template_vars_analysis), map_template_vars_analysis,vectorize=F)


print("Combine fusion level files")
print(paste0("Running: patient: ",patient$patient_identifier))


## Integrate starfusion and fusioncatcher fusion level files ----

patient_starfusion_fusion_level = read.table(sf_fusion_level_results_path,sep="\t",header=T)
patient_starfusion_fusion_level$patient_id=patient$patient_id
patient_starfusion_fusion_level = patient_starfusion_fusion_level %>% 
  mutate(patient_fusion =  paste(patient_id,fusion_name,sep="_"),
         patient_fusion_sv =  paste(patient_id,gup_sv_merged,gdw_sv_merged,fusion_name,sep="_"),
         patient_sv_id =  paste(patient_id,gup_sv_merged,gdw_sv_merged,sep="_"))

patient_fusioncatcher_fusion_level = read.table(fc_fusion_level_results_path,sep="\t",header=T)
patient_fusioncatcher_fusion_level$patient_id=patient$patient_id
patient_fusioncatcher_fusion_level = patient_fusioncatcher_fusion_level %>% 
  mutate(patient_fusion =  paste(patient_id,fusion_name,sep="_"),
         patient_fusion_sv =  paste(patient_id,gup_sv_merged,gdw_sv_merged,fusion_name,sep="_"),
         patient_sv_id =  paste(patient_id,gup_sv_merged,gdw_sv_merged,sep="_"))

cols_1 = names(patient_starfusion_fusion_level)
cols_2 = names(patient_fusioncatcher_fusion_level)
shared_cols=cols_1[cols_1 %in% cols_2]

if(FALSE) {
  merge_cols=c("fusion_name","patient_fusion","patient_fusion_sv")
  gene_anno_cols =   c("gup_gene_id", "gdw_gene_id", "gup_gene_type", "gdw_gene_type", "gup_ensembl_id", "gdw_ensembl_id")
  sv_anno_cols = c("sv_names", "gup_coordinate", "gdw_coordinate", "tools", "tumor_af", "normal_af", "tumor_af_spread", "normal_af_spread",
                   "tumor_af_mean", "normal_af_mean","svtype", "svlen", "svlen_spread")
  
  gene_anno_cols = gene_anno_cols[gene_anno_cols %in% shared_cols]
  sv_anno_cols = sv_anno_cols[sv_anno_cols %in% shared_cols]
  
  merged_df = merge(patient_starfusion_fusion_level[,c(merge_cols,sv_anno_cols,gene_anno_cols)],
                    patient_fusioncatcher_fusion_level[,c(merge_cols,sv_anno_cols,gene_anno_cols)],
                    by=c(merge_cols,sv_anno_cols,gene_anno_cols))
  
  nrow(merged_df)
  
  missing_shared_cols = shared_cols[!shared_cols %in% names(merged_df)]
  missing_shared_cols
}
ignore_shared_cols=c("gup_sf_breakpoint","gdw_sf_breakpoint","fusion_predictions", "gup_location", "gdw_location","overlap_gup_gdw_genebody",
                     "gup_distance_mean","gdw_distance_mean")

merge_cols=shared_cols[!shared_cols %in% ignore_shared_cols]
merged_df = merge(patient_starfusion_fusion_level[,merge_cols],
                  patient_fusioncatcher_fusion_level[,merge_cols],
                  by=merge_cols)

#make union   
fusion_level_df = rbind(  merged_df,
                          patient_starfusion_fusion_level[!patient_starfusion_fusion_level$patient_fusion_sv %in% merged_df$patient_fusion_sv,merge_cols],
                          patient_fusioncatcher_fusion_level[!patient_fusioncatcher_fusion_level$patient_fusion_sv %in% merged_df$patient_fusion_sv,merge_cols] )

fusion_level_df = fusion_level_df %>% mutate(rna_tools= ifelse(patient_fusion_sv %in% merged_df$patient_fusion_sv,"fusioncatcher, starfusion",
                                                               ifelse(patient_fusion_sv %in% patient_starfusion_fusion_level$patient_fusion_sv, "starfusion",
                                                                      ifelse(patient_fusion_sv %in% patient_fusioncatcher_fusion_level$patient_fusion_sv, "fusioncatcher",NA))))


#add ignore shared cols with _starfusion and _fusioncatcher
#rest can be safely added
sf_metadata = patient_starfusion_fusion_level %>% 
  dplyr::rename_with(.cols=ignore_shared_cols,.fn=function(x){paste0(x,"_starfusion")})

fc_metadata = patient_fusioncatcher_fusion_level %>% 
  dplyr::rename_with(.cols=ignore_shared_cols,.fn=function(x){paste0(x,"_fusioncatcher")})

fusion_level_df = fusion_level_df %>% left_join(sf_metadata,by=merge_cols) %>% left_join(fc_metadata,by=merge_cols)

write.table(fusion_level_df,merged_fusion_level_results_path,col.names = T,row.names = F,sep="\t")

