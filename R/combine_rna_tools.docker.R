## Fusion sq
## Combine RNA tools
## Last update: 2022-02-03
## 

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

#TMP local use only
patients_missing_sf =  c()
patients_missing_fc =  c()
patients_missing_sf_pred =  c()
patients_missing_fc_pred =  c()

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


if(length(Sys.glob(sf_fusion_level_results_path))==0 ) {
  print(paste0("No star fusion output found: ",sf_fusion_level_results_path))
  patients_missing_sf = c(patients_missing_sf,patient$patient_id)
  
}
if(length(Sys.glob(fc_fusion_level_results_path))==0) {
  print(paste0("No fusion catcher output found: ",fc_fusion_level_results_path))
  patients_missing_fc = c(patients_missing_fc,patient$patient_id)
}

if(length(Sys.glob(sf_fusion_anno_table_path))==0 ) {
  print(paste0("No star fusion output found: ",sf_fusion_anno_table_path))
  patients_missing_sf_pred = c(patients_missing_sf_pred,patient$patient_id)
  
}
if(length(Sys.glob(fc_fusion_anno_table_path))==0) {
  print(paste0("No fusion catcher output found: ",fc_fusion_anno_table_path))
  patients_missing_fc_pred = c(patients_missing_fc_pred,patient$patient_id)
}

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

if(FALSE) {
  ##TMP debug remove later
  merge_cols=c("fusion_name","patient_fusion","patient_identifier")
  gene_anno_cols =   c("gup_gene_id", "gdw_gene_id", "gup_gene_type", "gdw_gene_type")#, "gup_ensembl_id", "gdw_ensembl_id")
  rna_cols = c("gup_sf_breakpoint","gdw_sf_breakpoint")
  
  gene_anno_cols = gene_anno_cols[gene_anno_cols %in% shared_cols]
  rna_cols = rna_cols[rna_cols %in% shared_cols]
  
  merged_df = merge(patient_starfusion_predictions[,c(merge_cols,rna_cols,gene_anno_cols)],
                    patient_fusioncatcher_predictions[,c(merge_cols,rna_cols,gene_anno_cols)],
                    by=c(merge_cols,rna_cols,gene_anno_cols))
  
  nrow(merged_df)
  
  missing_shared_cols = shared_cols[!shared_cols %in% names(merged_df)]
  missing_shared_cols
}

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

if(FALSE) {
  ##TMP debug remove later
  
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




### collect cohort merged ----

map_template_vars=c('${input_dir}'=input_dir,'${output_dir}'=output_dir,
                    '${cohort_identifier}'=cohort_identifier,'${cohort_wdir}'=cohort_wdir)

reports_dir =  stri_replace_all_fixed(reports_dir_template,names(map_template_vars), map_template_vars,vectorize=F)
base_dir =  stri_replace_all_fixed(base_dir_template,names(map_template_vars), map_template_vars,vectorize=F)

map_template_vars=c(map_template_vars,'${reports_dir}'=reports_dir,'${base_dir}'=base_dir)

##  merge fusion annotation ends up in reports dir too
## TODO implement in collect cohort file and add 'else add single file' 

#tmp dots
map_template_vars_merged = c('${analysis_type}'=".merged",
                             map_template_vars )

cohort_fusion_level_results_path = stri_replace_all_fixed(cohort_fusion_level_results_path_template,names(map_template_vars_merged), map_template_vars_merged,vectorize=F)
#cohort_results_path = stri_replace_all_fixed(cohort_results_path_template,names(map_template_vars_merged), map_template_vars_merged,vectorize=F) 
fusion_overview_path = stri_replace_all_fixed(fusion_overview_path_template,names(map_template_vars_merged), map_template_vars_merged,vectorize=F) 
#cohort_supporting_svs_path = stri_replace_all_fixed(cohort_supporting_svs_path_template,names(map_template_vars_merged), map_template_vars_merged,vectorize=F)
#cohort_pairwise_overlap_merged_path = stri_replace_all_fixed(cohort_pairwise_overlap_merged_path_template,names(map_template_vars_merged), map_template_vars_merged,vectorize=F) 
#cohort_report_path = stri_replace_all_fixed(cohort_report_path_template,names(map_template_vars_merged), map_template_vars_merged,vectorize=F)


merged_fusion_level_results = data.frame(stringsAsFactors=FALSE) 
merged_predictions= data.frame(stringsAsFactors=FALSE) 
for(id in cohort$patient_identifier) {
  patient = filter(cohort,patient_identifier==id)
  
  #note order matters, this overrides the cohort analysis type filename 
  #tmp; can be removed if dot is placed right
  #character string for patient id because could be factor
  map_template_vars_patient = c('${analysis_type}'="merged.",
                                map_template_vars,
                                '${patient_identifier}'=as.character(patient$patient_identifier) )
  
  #merged is in reports dir and normal is in base dir
  fusion_anno_table_path = stri_replace_all_fixed(fusion_anno_table_merged_path_template,names(map_template_vars_patient), map_template_vars_patient,vectorize=F)
  #fusion_anno_table_path = stri_replace_all_fixed(fusion_anno_table_path_template,names(map_template_vars_patient), map_template_vars_patient,vectorize=F)
  if(length(Sys.glob(fusion_anno_table_path))==0) { next() }
  predictions = read.table(fusion_anno_table_path,header=T,sep="\t") 
  predictions$patient_id = patient$patient_id
  merged_predictions = rbind(merged_predictions,as.data.frame(predictions))
  
  fusion_level_results_path = stri_replace_all_fixed(fusion_level_results_path_template,names(map_template_vars_patient), map_template_vars_patient,vectorize=F)
  if(length(Sys.glob(fusion_level_results_path))==0) { next() }
  fusion_level_results = read.table(fusion_level_results_path,header=T,sep="\t") 
  fusion_level_results$patient_id = patient$patient_id
  merged_fusion_level_results = rbind(merged_fusion_level_results,as.data.frame(fusion_level_results))
  
}

## TODO reannotate location precse is combine wgs
merged_fusion_level_results = merged_fusion_level_results %>% 
  mutate(location_precise_fusioncatcher = grepl(paste(proximate_bp, collapse = "|"),gup_location_fusioncatcher) & 
           grepl(paste(proximate_bp, collapse = "|"),gdw_location_fusioncatcher),
         location_precise_starfusion = grepl(paste(proximate_bp, collapse = "|"),gup_location_starfusion) & 
                  grepl(paste(proximate_bp, collapse = "|"),gdw_location_starfusion),
         precise_confident_fusioncatcher = location_precise_fusioncatcher & grepl(",",tools),
         precise_confident_starfusion = location_precise_starfusion & grepl(",",tools))

#check merged_fusion_level_results %>% filter(rna_tools=="starfusion" & precise_confident_starfusion==T & !precise_confident)

merged_fusion_level_results = merged_fusion_level_results %>% 
  mutate(precise_confident = (precise_confident_fusioncatcher==T|precise_confident_starfusion==T),
         precise_confident = ifelse(is.na(precise_confident),F,precise_confident))


#to check
merged_fusion_level_results %>% filter(is.na(precise_confident)|is.na(precise_confident_starfusion)|is.na(precise_confident_fusioncatcher)) %>% 
  select(precise_confident_fusioncatcher,precise_confident_starfusion,precise_confident) %>% unique()

write.table(merged_fusion_level_results,cohort_fusion_level_results_path,quote = FALSE,sep = "\t",row.names=FALSE)
write.table(merged_predictions,fusion_overview_path,quote = FALSE,sep = "\t",row.names=FALSE)

