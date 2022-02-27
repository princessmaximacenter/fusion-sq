## Fusion sq
## Combine RNA tools
## Last update: 2022-02-03
## 2022-02-16 removed tmp code

if(FALSE) {
  #local
  source("~/fusion_sq/R/default.conf")
  source("~/fusion_sq/R/default.docker.local.conf")
  source("~/fusion_sq/run/fusion_sq/fusion_sq.conf")
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
  library(GenomicRanges, quietly=TRUE)
  library(dplyr)
  library(stringi)
})

source(paste0(script_dir,"functions.general.R"))
source(paste0(script_dir,"functions.svs.R"))

#settings

prediction_ignore_shared_cols=c("annots","identifier","fusion_identifier","anno_healthy_chimera","fusion_bp_distance","gup_gene_type","gdw_gene_type")

#interaction between rna bp and sv bp
fusion_level_ignore_shared_cols=c("gup_sf_breakpoint","gdw_sf_breakpoint","fusion_predictions", "gup_location", "gdw_location","overlap_gup_gdw_genebody",
                                  "gup_distance_mean","gdw_distance_mean","gup_gene_type","gdw_gene_type", "tools_any_wgs",
                                  "location_precise", "precise_confident", "not_precise_confident","specific_sv")

cohort = read.table(patient_table_path,sep = "\t", header=T,stringsAsFactors = F)

for(pid in cohort$patient_id) {
  patient = cohort %>% filter(patient_id==pid)

  ## remove loop later!
  
if(patient$patient_id =="") {
  print("patient$patient_id needs to be specified")
  #quit()
}


#order of arguments matters
map_template_vars=c('${input_dir}'=input_dir,'${output_dir}'=output_dir,
                    '${cohort_identifier}'=cohort_identifier,'${cohort_wdir}'=cohort_wdir,
                    '${patient_identifier}'=patient$patient_identifier)

reports_dir =  stri_replace_all_fixed(reports_dir_template,names(map_template_vars), map_template_vars,vectorize=F)
base_dir =  stri_replace_all_fixed(base_dir_template,names(map_template_vars), map_template_vars,vectorize=F)

map_template_vars=c(map_template_vars,'${reports_dir}'=reports_dir,'${base_dir}'=base_dir)

#Paths 
## TODO: remove dots later
map_template_vars_analysis = c(map_template_vars,'${analysis_type}'="starfusion.")
sf_fusion_level_results_path = stri_replace_all_fixed(fusion_level_results_path_template,names(map_template_vars_analysis), map_template_vars_analysis,vectorize=F)
sf_fusion_anno_table_path = stri_replace_all_fixed(fusion_anno_table_path_template,names(map_template_vars_analysis), map_template_vars_analysis,vectorize=F)

map_template_vars_analysis = c(map_template_vars,'${analysis_type}'="fusioncatcher.")
fc_fusion_level_results_path = stri_replace_all_fixed(fusion_level_results_path_template,names(map_template_vars_analysis), map_template_vars_analysis,vectorize=F)
fc_fusion_anno_table_path = stri_replace_all_fixed(fusion_anno_table_path_template,names(map_template_vars_analysis), map_template_vars_analysis,vectorize=F)

map_template_vars_analysis = c(map_template_vars,'${analysis_type}'="merged.")
merged_fusion_level_results_path = stri_replace_all_fixed(fusion_level_results_path_template,names(map_template_vars_analysis), map_template_vars_analysis,vectorize=F)
merged_fusion_anno_table_path = stri_replace_all_fixed(fusion_anno_table_path_template,names(map_template_vars_analysis), map_template_vars_analysis,vectorize=F)


print(paste0("Running: patient: ",patient$patient_identifier))

## Integrate starfusion and fusioncatcher predictions / fusion annotation files ----

print("Combine predictions / fusion annotation files")


if(length(Sys.glob(sf_fusion_anno_table_path))==0 | length(Sys.glob(fc_fusion_anno_table_path))==0 ) {
  #quit()
  next()
}

patient_starfusion_predictions = read.table(sf_fusion_anno_table_path,sep="\t",header=T)
patient_starfusion_predictions$patient_id=patient$patient_id
patient_starfusion_predictions = patient_starfusion_predictions %>% 
  harmonize_gene_identifiers() %>%
  make_identifiers_cohort_analysis()
        

patient_fusioncatcher_predictions = read.table(fc_fusion_anno_table_path,sep="\t",header=T)
patient_fusioncatcher_predictions$patient_id=patient$patient_id
patient_fusioncatcher_predictions = patient_fusioncatcher_predictions %>% 
  harmonize_gene_identifiers() %>%
  make_identifiers_cohort_analysis()

cols_1 = names(patient_starfusion_predictions)
cols_2 = names(patient_fusioncatcher_predictions)
shared_cols=cols_1[cols_1 %in% cols_2]

merge_cols=shared_cols[!shared_cols %in% prediction_ignore_shared_cols]
merged_df = merge(patient_starfusion_predictions[,merge_cols],
                  patient_fusioncatcher_predictions[,merge_cols],
                  by=merge_cols)


#anti join to excl for merging => or just merge and use unique... => no cannot because of rna_tools
patient_starfusion_predictions_wo_merged = patient_starfusion_predictions[,merge_cols] %>% anti_join(merged_df[,merge_cols])
patient_fusioncatcher_predictions_wo_merged = patient_fusioncatcher_predictions[,merge_cols] %>% anti_join(merged_df[,merge_cols])

if(nrow(merged_df)>0)  merged_df$rna_tools = "fusioncatcher, starfusion"
patient_starfusion_predictions_wo_merged$rna_tools = "starfusion"
patient_fusioncatcher_predictions_wo_merged$rna_tools = "fusioncatcher"

#make union   
predictions_df = rbind(  merged_df,
                         patient_starfusion_predictions_wo_merged,
                         patient_fusioncatcher_predictions_wo_merged)

#add ignore shared cols with _starfusion and _fusioncatcher
#rest can be safely added
sf_metadata = patient_starfusion_predictions %>% 
  dplyr::rename_with(.cols=prediction_ignore_shared_cols,.fn=function(x){paste0(x,"_starfusion")})

fc_metadata = patient_fusioncatcher_predictions %>% 
  dplyr::rename_with(.cols=prediction_ignore_shared_cols,.fn=function(x){paste0(x,"_fusioncatcher")})

predictions_df = predictions_df %>% left_join(sf_metadata,by=merge_cols) %>% left_join(fc_metadata,by=merge_cols)

write.table(predictions_df,merged_fusion_anno_table_path,col.names = T,row.names = F,sep="\t")


## Integrate starfusion and fusioncatcher fusion level files ----


print("Combine fusion level files")


if(length(Sys.glob(sf_fusion_level_results_path))==0 | length(Sys.glob(fc_fusion_level_results_path))==0 ) {
  #quit()
  next()
}

patient_starfusion_fusion_level = read.table(sf_fusion_level_results_path,sep="\t",header=T)
patient_starfusion_fusion_level$patient_id=patient$patient_id
patient_starfusion_fusion_level = patient_starfusion_fusion_level %>% 
  harmonize_gene_identifiers() %>%
  make_identifiers_cohort_analysis()

patient_fusioncatcher_fusion_level = read.table(fc_fusion_level_results_path,sep="\t",header=T)
patient_fusioncatcher_fusion_level$patient_id=patient$patient_id
patient_fusioncatcher_fusion_level = patient_fusioncatcher_fusion_level %>% 
  harmonize_gene_identifiers() %>%
  make_identifiers_cohort_analysis()

cols_1 = names(patient_starfusion_fusion_level)
cols_2 = names(patient_fusioncatcher_fusion_level)
shared_cols=cols_1[cols_1 %in% cols_2]


merge_cols=shared_cols[!shared_cols %in% fusion_level_ignore_shared_cols]
merged_df = merge(patient_starfusion_fusion_level[,merge_cols],
                  patient_fusioncatcher_fusion_level[,merge_cols],
                  by=merge_cols)

if(nrow(merged_df)>0)  merged_df$rna_tools = "fusioncatcher, starfusion"
patient_starfusion_fusion_level$rna_tools = "starfusion"
patient_fusioncatcher_fusion_level$rna_tools = "fusioncatcher"

#make union   
fusion_level_df = rbind(  merged_df,
                          patient_starfusion_fusion_level[!patient_starfusion_fusion_level$patient_fusion_sv %in% merged_df$patient_fusion_sv,c(merge_cols,"rna_tools")],
                          patient_fusioncatcher_fusion_level[!patient_fusioncatcher_fusion_level$patient_fusion_sv %in% merged_df$patient_fusion_sv,c(merge_cols,"rna_tools")] )


#add ignore shared cols with _starfusion and _fusioncatcher
#rest can be safely added
#remove rna tools col
sf_metadata = patient_starfusion_fusion_level %>% select(-rna_tools) %>%
  dplyr::rename_with(.cols=fusion_level_ignore_shared_cols,.fn=function(x){paste0(x,"_starfusion")})

fc_metadata = patient_fusioncatcher_fusion_level %>% select(-rna_tools) %>%
  dplyr::rename_with(.cols=fusion_level_ignore_shared_cols,.fn=function(x){paste0(x,"_fusioncatcher")})

fusion_level_df = fusion_level_df %>% left_join(sf_metadata,by=merge_cols) %>% left_join(fc_metadata,by=merge_cols)

write.table(fusion_level_df,merged_fusion_level_results_path,col.names = T,row.names = F,sep="\t")
}

