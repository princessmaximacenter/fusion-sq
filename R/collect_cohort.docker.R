## Fusion-sq 
## Collect cohort 
## Last update: 2021-04-11
# Collect all patients
# Prepare for cohort-level calculations and annotation

## Update 2022-01-08 path local overrides and docker version


if(FALSE){
  #set paths for hpc
  source("/hpc/pmc_gen/ivanbelzen/github_fusion_sq/fusion-sq/R/default.conf")
  ## HPC config overrides
  source("/hpc/pmc_gen/ivanbelzen/github_fusion_sq/fusion-sq/R/default.docker.conf")
  #HPC doesnt use argparser but patient specific config instead 
  source("/hpc/pmc_gen/ivanbelzen/github_fusion_sq/fusion-sq/run/fusion_sq/fusion_sq.conf")
}


if(FALSE){
  #set paths for local use
  source("~/fusion_sq/R/default.conf")
  
  #override local for templating
  source("~/fusion_sq/R/default.docker.local.conf")
  #this is the local version of hpc config
  #source("/hpc/pmc_gen/ivanbelzen/github_fusion_sq/fusion-sq/run/fusion_sq/fusion_sq.conf")
  
  #containing the cohort level config things   
  source("~/fusion_sq/run/fusion_sq/fusion_sq.conf")
  
  #todo: make consistent with config files, this is the first time using cohort in the pipeline and patient_table_path is beter than what is now patient_table
  #output dir is where everything will be stored in this github version and also on hpc. 
  #templating should take care of differences between HPC and local => see TODO in default.conf
  
  analysis_type="fusioncatcher"
}

suppressPackageStartupMessages({
  library(tidyverse, quietly=TRUE)
  library(stringr, quietly=TRUE)
  library(stringi)
})

source(paste0(script_dir,"functions.general.R"))

if(!exists("patient_table_path") | length(Sys.glob(patient_table_path))!=1) {
  print(paste0("Cohort file not available: ",patient_table_path))
  #quit()
}


if(!exists("analysis_type")) {
  analysis_type="starfusion"
}


##  make templates where you either fill ${analysis_type}=> analysis_type,"." if fusion catcher and temporarily if star fusion just "" 
analysis_type_filename=paste0(".",analysis_type)

#order of arguments matters
map_template_vars=c('${input_dir}'=input_dir,'${output_dir}'=output_dir,
                    '${cohort_identifier}'=cohort_identifier,'${cohort_wdir}'=cohort_wdir,
                    '${analysis_type}'=analysis_type_filename)


#fill templates
reports_dir =  stri_replace_all_fixed(reports_dir_template,names(map_template_vars), map_template_vars,vectorize=F)
base_dir =  stri_replace_all_fixed(base_dir_template,names(map_template_vars), map_template_vars,vectorize=F)

map_template_vars=c(map_template_vars,'${reports_dir}'=reports_dir,'${base_dir}'=base_dir)

cohort_fusion_level_results_path = stri_replace_all_fixed(cohort_fusion_level_results_path_template,names(map_template_vars), map_template_vars,vectorize=F)
cohort_results_path = stri_replace_all_fixed(cohort_results_path_template,names(map_template_vars), map_template_vars,vectorize=F) 
fusion_overview_path = stri_replace_all_fixed(fusion_overview_path_template,names(map_template_vars), map_template_vars,vectorize=F) 
cohort_supporting_svs_path = stri_replace_all_fixed(cohort_supporting_svs_path_template,names(map_template_vars), map_template_vars,vectorize=F)
cohort_pairwise_overlap_merged_path = stri_replace_all_fixed(cohort_pairwise_overlap_merged_path_template,names(map_template_vars), map_template_vars,vectorize=F) 
cohort_report_path = stri_replace_all_fixed(cohort_report_path_template,names(map_template_vars), map_template_vars,vectorize=F)

# Exit if exists, otherwise make
if(length(Sys.glob(cohort_fusion_level_results_path))==1 ){ 
  print(paste0("File already exists: ",cohort_fusion_level_results_path))
  #quit() 
} else if(length(Sys.glob(cohort_report_path))==1){ 
  print(paste0("Clear cohort report file first, prevent inconsistencies: ",cohort_report_path))
  #quit() 
}

#for patient level output


## Read in cohort
cohort = read.table(patient_table_path,sep = "\t", header=T,stringsAsFactors = T)

## Per patient, collect output if files exist
cohort_results = data.frame(stringsAsFactors=FALSE) 
cohort_fusion_level_results = data.frame(stringsAsFactors=FALSE) 
fusion_overview = data.frame(stringsAsFactors = FALSE)
supporting_svs_df = data.frame(stringsAsFactors = FALSE)
pairwise_overlap_merged_df = data.frame(stringsAsFactors = FALSE)

## TODO adjust for running 2 rna tools as loop.
# if fusion level merge df and fusion annotation merge df are available. 
# otherwise fusion level output and predictions of either.



for(analysis_type in c("starfusion","fusioncatcher")){
  #Temporary code:
  #not needed anymore if dots are added to proper location in template
  analysis_type_filename_patient=paste0(analysis_type,".")
  ## end temp
  


for(id in cohort$patient_identifier) {
  patient = filter(cohort,patient_identifier==id)
 
  #note order matters, this overrides the cohort analysis type filename 
  #tmp; can be removed if dot is placed right
  #character string for patient id because could be factor
  map_template_vars_patient = c('${analysis_type}'=analysis_type_filename_patient,
                                map_template_vars,
                                '${patient_identifier}'=as.character(patient$patient_identifier) )
  
  fusion_anno_table_path = stri_replace_all_fixed(fusion_anno_table_path_template,names(map_template_vars_patient), map_template_vars_patient,vectorize=F)
  matching_intervals_path = stri_replace_all_fixed(matching_intervals_path_template,names(map_template_vars_patient), map_template_vars_patient,vectorize=F)
  fusion_level_results_path = stri_replace_all_fixed(fusion_level_results_path_template,names(map_template_vars_patient), map_template_vars_patient,vectorize=F)
  matching_bp_path = stri_replace_all_fixed(matching_bp_path_template,names(map_template_vars_patient), map_template_vars_patient,vectorize=F)
  pairwise_overlap_merged_path = stri_replace_all_fixed(pairwise_overlap_merged_path_template ,names(map_template_vars_patient), map_template_vars_patient,vectorize=F)
  supporting_svs_path = stri_replace_all_fixed(supporting_svs_path_template,names(map_template_vars_patient), map_template_vars_patient,vectorize=F)

  
  # to be able to annotate with overlap variable
    
  if(length(Sys.glob(fusion_anno_table_path))!=1){ next()}
  if(length(Sys.glob(matching_intervals_path))!=1){ next()}
    
  matching_intervals_table = read.table(matching_intervals_path,header=T, sep="\t")
  
  fusion_anno_table=read.table(fusion_anno_table_path,header=T,sep="\t")
  fusion_anno_table = fusion_anno_table %>% left_join(matching_intervals_table[,c("identifier","overlap_gup_gdw_adjacent_intron","overlap_gup_gdw_genebody")], by=c("identifier"))
  fusion_anno_table$patient_id = patient$patient_id
  fusion_anno_table$rna_tools=analysis_type
  fusion_overview = rbind_no_colmatch(fusion_anno_table,fusion_overview)
   
  if(length(Sys.glob(matching_bp_path))!=1){ next()}
  
  matching_bps = read.table(matching_bp_path,header=T,sep="\t") 
  if(nrow(matching_bps)>0) {
    matching_bps$patient_id = patient$patient_id 
    matching_bps$rna_tools=analysis_type
    cohort_results = rbind_no_colmatch(cohort_results,as.data.frame(matching_bps))
  
    fusion_level_results = read.table(fusion_level_results_path,header=T,sep="\t") 
    fusion_level_results$patient_id = patient$patient_id
    fusion_level_results$rna_tools=analysis_type
    cohort_fusion_level_results = rbind_no_colmatch(cohort_fusion_level_results,as.data.frame(fusion_level_results))
    
    
   if(length(Sys.glob(supporting_svs_path))!=1){next()}
    
   supporting_svs = read.table(supporting_svs_path,header=T,sep="\t",stringsAsFactors = F) 
   supporting_svs$patient_id = patient$patient_id
    supporting_svs$rna_tools=analysis_type
   supporting_svs_df = rbind_no_colmatch(supporting_svs_df,supporting_svs)
   
   
   
   if(length(Sys.glob(pairwise_overlap_merged_path))!=1){next()}
   
   pairwise_overlap_merged = read.table(pairwise_overlap_merged_path,header=T,sep="\t",stringsAsFactors = F) 
   pairwise_overlap_merged$patient_id = patient$patient_id
   pairwise_overlap_merged$rna_tools=analysis_type
   pairwise_overlap_merged_df = rbind_no_colmatch(pairwise_overlap_merged_df,pairwise_overlap_merged)
   
  }
}
}

## Unique identifier and sanity checks 
cohort_fusion_level_results = cohort_fusion_level_results %>%
  mutate(fusion_name = str_replace(str_replace(fusion_name,fixed("IGH@"),"IGH"),fixed("IGH-@-ext"),"IGH"),
         patient_fusion =  paste(patient_id,fusion_name,sep="_"),
         patient_fusion_sv =  paste(patient_id,gup_sv_merged,gdw_sv_merged,fusion_name,sep="_"),
         patient_sv_id =  paste(patient_id,gup_sv_merged,gdw_sv_merged,sep="_"))

if(nrow(cohort_fusion_level_results %>% filter(gup_sv_merged!=gdw_sv_merged & svtype!= "CTX" & specific_sv == T & grepl(",",tools)))>0) {
  print("WARNING: unexpected gup/gdw sv merged difference (non CTX and multi tools)")
  cohort_fusion_level_results %>% filter(gup_sv_merged!=gdw_sv_merged & svtype!= "CTX" & specific_sv == T & grepl(",",tools)) %>% 
    select(patient_fusion,gup_sv_merged,gdw_sv_merged,svtype,tools)
}

cohort_results = cohort_results %>%
  mutate(fusion_name = str_replace(str_replace(fusion_name,fixed("IGH@"),"IGH"),fixed("IGH-@-ext"),"IGH"),
       patient_fusion =  paste(patient_id,fusion_name,sep="_"),
       patient_fusion_sv =  paste(patient_id,gup_sv_merged,gdw_sv_merged,fusion_name,sep="_"),
       patient_sv_id =  paste(patient_id,gup_sv_merged,gdw_sv_merged,sep="_"))


cohort_results = cohort_results %>% dplyr::rename(tool= gup_tool) %>% select(-gdw_tool)

cohort_results = cohort_results %>% group_by(patient_fusion_sv) %>% mutate(svtype=toString(unique(sort(c(gup_svtype,gdw_svtype)))))

if(nrow(cohort_results %>% filter(gup_svtype != gdw_svtype & gup_location !="composite"))>0){
  print("WARNING: gup/gdw svtypes dont match")
  print(cohort_results %>% filter(gup_svtype != gdw_svtype& gup_location !="composite")) %>% ungroup() %>% select(patient_fusion,gup_location,gup_svtype,gdw_svtype,svtype)
}



fusion_overview = fusion_overview %>% 
  mutate(fusion_name = str_replace(str_replace(fusion_name,fixed("IGH@"),"IGH"),fixed("IGH-@-ext"),"IGH"),
         patient_fusion =  paste(patient_id,fusion_name,sep="_"))


write.table(cohort_fusion_level_results,cohort_fusion_level_results_path,quote = FALSE,sep = "\t",row.names=FALSE)

write.table(cohort_results,cohort_results_path,quote = FALSE,sep = "\t",row.names=FALSE)

write.table(fusion_overview,fusion_overview_path,quote = FALSE,sep = "\t",row.names=FALSE)

supporting_svs_df$bp_name = paste0(supporting_svs_df$patient_id,"_",supporting_svs_df$bp_name)
write.table(supporting_svs_df,cohort_supporting_svs_path,quote = F,sep = "\t",row.names=F,col.names = T)

write.table(pairwise_overlap_merged_df,cohort_pairwise_overlap_merged_path,quote = F,sep = "\t",row.names=F,col.names = T)



