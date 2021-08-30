## Fusion-sq 
## Collect cohort 
## Last update: 2021-04-11
# Collect all patients
# Prepare for cohort-level calculations and annotation

suppressPackageStartupMessages({
  library(tidyverse, quietly=TRUE)
  library(stringr, quietly=TRUE)
})

source("R/default.conf") #in script dir
source(paste0(script_dir,"functions.general.R"))

## Inputs defined in config file

## Output:
cohort_fusion_level_results_path = paste0(reports_dir,cohort_fusion_level_results_outfile)
cohort_results_path = paste0(reports_dir,cohort_matching_results_outfile)
fusion_overview_path = paste0(reports_dir,fusion_overview_outfile)
cohort_supporting_svs_path = paste0(reports_dir,cohort_supporting_svs_outfile)
cohort_pairwise_overlap_merged_path = paste0(reports_dir,cohort_pairwise_overlap_merged_outfile)

#to check if exists:
cohort_report_path = paste0(reports_dir,cohort_report_outfile)

# Exit if exists, otherwise make
if(length(Sys.glob(cohort_fusion_level_results_path))==1 ){ 
  print(paste0("File already exists: ",cohort_fusion_level_results_path))
  #quit() 
} else if(length(Sys.glob(cohort_report_path))==1){ 
  print(paste0("Clear cohort report file first, prevent inconsistencies: ",cohort_report_path))
  #quit() 
}

## Read in cohort
patient_metadata = read.table(patient_table,sep = "\t", header=T,stringsAsFactors = T)

## Per patient, collect output if files exist
cohort_results = data.frame(stringsAsFactors=FALSE) 
cohort_fusion_level_results = data.frame(stringsAsFactors=FALSE) 
fusion_overview = data.frame(stringsAsFactors = FALSE)
supporting_svs_df = data.frame(stringsAsFactors = FALSE)
pairwise_overlap_merged_df = data.frame(stringsAsFactors = FALSE)

for(id in patient_metadata$patient_identifier) {
  patient = filter(patient_metadata,patient_identifier==id)
  
  # to be able to annotate with overlap variable
  matching_intervals_path = paste0(base_dir,matching_intervals_outfile,patient$patient_identifier,".tsv")
  matching_intervals_table = read.table(matching_intervals_path,header=T, sep="\t")
  
  fusion_anno_table_path = paste0(base_dir,fusion_annotation_outfile,patient$patient_identifier,".tsv")
  fusion_anno_table=read.table(fusion_anno_table_path,header=T,sep="\t")
  fusion_anno_table = fusion_anno_table %>% left_join(matching_intervals_table[,c("identifier","overlap_gup_gdw_adjacent_intron","overlap_gup_gdw_genebody")], by=c("identifier"))
  fusion_anno_table$patient_id = patient$patient_id
  fusion_overview = rbind(fusion_anno_table,fusion_overview)
   
  matching_bp_path = paste0(reports_dir,matching_results_outfile,patient$patient_identifier,".tsv")
  if(length(Sys.glob(matching_bp_path))!=1){ next()}
  
  matching_bps = read.table(matching_bp_path,header=T,sep="\t") 
  if(nrow(matching_bps)>0) {
    matching_bps$patient_id = patient$patient_id 
    cohort_results = rbind(cohort_results,as.data.frame(matching_bps))
  
    fusion_level_results_path = paste0(reports_dir,fusion_level_results_outfile,patient$patient_identifier,".tsv")
    fusion_level_results = read.table(fusion_level_results_path,header=T,sep="\t") 
    fusion_level_results$patient_id = patient$patient_id
    cohort_fusion_level_results = rbind(cohort_fusion_level_results,as.data.frame(fusion_level_results))
    
    
   supporting_svs_path = paste0(reports_dir,supporting_svs_outfile,patient$patient_identifier,".tsv")
        
   if(length(Sys.glob(supporting_svs_path))!=1){next()}
    
   supporting_svs = read.table(supporting_svs_path,header=T,sep="\t",stringsAsFactors = F) 
   supporting_svs$patient_id = patient$patient_id
        
   supporting_svs_df = rbind(supporting_svs_df,supporting_svs)
   
   
   
   pairwise_overlap_merged_path = paste0(reports_dir,pairwise_overlap_merged_outfile,patient$patient_identifier,".tsv")
   
   if(length(Sys.glob(pairwise_overlap_merged_path))!=1){next()}
   
   pairwise_overlap_merged = read.table(pairwise_overlap_merged_path,header=T,sep="\t",stringsAsFactors = F) 
   pairwise_overlap_merged$patient_id = patient$patient_id
   
   pairwise_overlap_merged_df = rbind(pairwise_overlap_merged_df,pairwise_overlap_merged)
   
  }
}

## Unique identifier and sanity checks 
cohort_fusion_level_results$patient_fusion = paste0(cohort_fusion_level_results$patient_id,"_",cohort_fusion_level_results$fusion_name)
cohort_fusion_level_results = cohort_fusion_level_results %>% mutate(patient_fusion_sv =  paste(patient_id,gup_sv_merged,gdw_sv_merged,fusion_name,sep="_"))


if(nrow(cohort_fusion_level_results %>% filter(gup_sv_merged!=gdw_sv_merged & svtype!= "CTX" & specific_sv == T & grepl(",",tools)))>0) {
  print("WARNING: unexpected gup/gdw sv merged difference (non CTX and multi tools)")
  cohort_fusion_level_results %>% filter(gup_sv_merged!=gdw_sv_merged & svtype!= "CTX" & specific_sv == T & grepl(",",tools)) %>% 
    select(patient_fusion,gup_sv_merged,gdw_sv_merged,svtype,tools)
}

cohort_results = cohort_results %>% mutate(patient_fusion_sv = paste(patient_id,gup_sv_merged,gdw_sv_merged,fusion_name,sep="_"))
cohort_results$patient_fusion = paste0(cohort_results$patient_id,"_",cohort_results$fusion_name)

if(nrow(cohort_results %>% filter(gup_tool != gdw_tool))>0) {
  print("WARNING: gup/gdw tools dont match")
  print(cohort_results %>% filter(gup_tool != gdw_tool))
}
cohort_results = cohort_results %>% dplyr::rename(tool= gup_tool) %>% select(-gdw_tool)

cohort_results = cohort_results %>% group_by(patient_fusion_sv) %>% mutate(svtype=toString(unique(sort(c(gup_svtype,gdw_svtype)))))

if(nrow(cohort_results %>% filter(gup_svtype != gdw_svtype & gup_location !="composite"))>0){
  print("WARNING: gup/gdw svtypes dont match")
  print(cohort_results %>% filter(gup_svtype != gdw_svtype& gup_location !="composite")) %>% ungroup() %>% select(patient_fusion,gup_location,gup_svtype,gdw_svtype,svtype)
}


write.table(cohort_fusion_level_results,cohort_fusion_level_results_path,quote = FALSE,sep = "\t",row.names=FALSE)

write.table(cohort_results,cohort_results_path,quote = FALSE,sep = "\t",row.names=FALSE)

fusion_overview$patient_fusion = paste0(fusion_overview$patient_id,"_",fusion_overview$fusion_name)

write.table(fusion_overview,fusion_overview_path,quote = FALSE,sep = "\t",row.names=FALSE)

supporting_svs_df$bp_name = paste0(supporting_svs_df$patient_id,"_",supporting_svs_df$bp_name)
write.table(supporting_svs_df,cohort_supporting_svs_path,quote = F,sep = "\t",row.names=F,col.names = T)

write.table(pairwise_overlap_merged_df,cohort_pairwise_overlap_merged_path,quote = F,sep = "\t",row.names=F,col.names = T)


