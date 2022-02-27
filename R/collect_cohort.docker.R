## Fusion-sq 
## Collect cohort 
## Last update: 2021-04-11
# Collect all patients
# Prepare for cohort-level calculations and annotation

## Update 2022-01-08 path local overrides and docker version
## update 2022-02-26 removed tmp code

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


## Read in cohort
cohort = read.table(patient_table_path,sep = "\t", header=T,stringsAsFactors = T)


#order of arguments matters
map_template_vars=c('${input_dir}'=input_dir,'${output_dir}'=output_dir,
                    '${cohort_identifier}'=cohort_identifier,'${cohort_wdir}'=cohort_wdir)


#fill templates
reports_dir =  stri_replace_all_fixed(reports_dir_template,names(map_template_vars), map_template_vars,vectorize=F)
base_dir =  stri_replace_all_fixed(base_dir_template,names(map_template_vars), map_template_vars,vectorize=F)

map_template_vars=c(map_template_vars,'${reports_dir}'=reports_dir,'${base_dir}'=base_dir)

map_template_vars_merged = c('${analysis_type}'=".merged",
                             map_template_vars )

cohort_fusion_level_results_path = stri_replace_all_fixed(cohort_fusion_level_results_path_template,names(map_template_vars_merged), map_template_vars_merged,vectorize=F)
cohort_results_path = stri_replace_all_fixed(cohort_results_path_template,names(map_template_vars_merged), map_template_vars_merged,vectorize=F) 
fusion_overview_path = stri_replace_all_fixed(fusion_overview_path_template,names(map_template_vars_merged), map_template_vars_merged,vectorize=F) 
cohort_report_path = stri_replace_all_fixed(cohort_report_path_template,names(map_template_vars_merged), map_template_vars_merged,vectorize=F)

cohort_supporting_svs_path = stri_replace_all_fixed(cohort_supporting_svs_path_template,names(map_template_vars_merged), map_template_vars_merged,vectorize=F)
cohort_pairwise_overlap_merged_path = stri_replace_all_fixed(cohort_pairwise_overlap_merged_path_template,names(map_template_vars_merged), map_template_vars_merged,vectorize=F) 

# Exit if exists, otherwise make
if(length(Sys.glob(cohort_fusion_level_results_path))==1 ){ 
  print(paste0("File already exists: ",cohort_fusion_level_results_path))
  #quit() 
} else if(length(Sys.glob(cohort_report_path))==1){ 
  print(paste0("Clear cohort report file first, prevent inconsistencies: ",cohort_report_path))
  #quit() 
}

## Collect merged fusion data ----

# Per patient, collect output if files exist
cohort_fusion_level_results = data.frame(stringsAsFactors=FALSE) 
cohort_predictions = data.frame(stringsAsFactors = FALSE)

for(id in cohort$patient_identifier) {
  patient = filter(cohort,patient_identifier==id)
  
  # if fusion level merge df and fusion annotation merge df are available. 
  # else fusion level output and predictions of either.
  
  #tmp; can be removed if dot is placed right
  #character string for patient id because could be factor
  map_template_vars_patient = c(map_template_vars,
                                '${patient_identifier}'=as.character(patient$patient_identifier) )
  
  
  map_template_vars_analysis = c(map_template_vars_patient,'${analysis_type}'="merged.")
  merged_fusion_level_results_path = stri_replace_all_fixed(fusion_level_results_path_template,names(map_template_vars_analysis), map_template_vars_analysis,vectorize=F)
  merged_fusion_anno_table_path = stri_replace_all_fixed(fusion_anno_table_merged_path_template,names(map_template_vars_analysis), map_template_vars_analysis,vectorize=F)
  #Note that merged fusion_anno_table is in reports dir and normal is in base dir
  
  map_template_vars_analysis = c(map_template_vars_patient,'${analysis_type}'="starfusion.")
  sf_fusion_level_results_path = stri_replace_all_fixed(fusion_level_results_path_template,names(map_template_vars_analysis), map_template_vars_analysis,vectorize=F)
  sf_fusion_anno_table_path = stri_replace_all_fixed(fusion_anno_table_path_template,names(map_template_vars_analysis), map_template_vars_analysis,vectorize=F)
  
  map_template_vars_analysis = c(map_template_vars_patient,'${analysis_type}'="fusioncatcher.")
  fc_fusion_level_results_path = stri_replace_all_fixed(fusion_level_results_path_template,names(map_template_vars_analysis), map_template_vars_analysis,vectorize=F)
  fc_fusion_anno_table_path = stri_replace_all_fixed(fusion_anno_table_path_template,names(map_template_vars_analysis), map_template_vars_analysis,vectorize=F)
  
  
  if(length(Sys.glob(merged_fusion_anno_table_path))==1) { 
    predictions = read.table(merged_fusion_anno_table_path,header=T,sep="\t") 
    predictions$patient_id = patient$patient_id
    cohort_predictions = rbind(cohort_predictions,as.data.frame(predictions))
    
  } else {
    
    if(length(Sys.glob(sf_fusion_anno_table_path))==1) {
      patient_starfusion_predictions = read.table(sf_fusion_anno_table_path,sep="\t",header=T)
      patient_starfusion_predictions$patient_id=patient$patient_id
      patient_starfusion_predictions = patient_starfusion_predictions %>% 
        harmonize_gene_identifiers() %>%
        make_identifiers_cohort_analysis()
      
      patient_starfusion_predictions$rna_tools = "starfusion"
      patient_starfusion_predictions = patient_starfusion_predictions %>% 
        dplyr::rename_with(.cols=prediction_ignore_shared_cols,.fn=function(x){paste0(x,"_starfusion")})
      
      cohort_predictions = rbind_no_colmatch(cohort_predictions,as.data.frame(patient_starfusion_predictions))
    }
    if(length(Sys.glob(fc_fusion_anno_table_path))==1 ) {
      patient_fusioncatcher_predictions = read.table(fc_fusion_anno_table_path,sep="\t",header=T)
      patient_fusioncatcher_predictions$patient_id=patient$patient_id
      patient_fusioncatcher_predictions = patient_fusioncatcher_predictions %>% 
        harmonize_gene_identifiers() %>%
        make_identifiers_cohort_analysis()
      patient_fusioncatcher_predictions$rna_tools = "fusioncatcher"
      
      patient_fusioncatcher_predictions = patient_fusioncatcher_predictions %>% 
        dplyr::rename_with(.cols=prediction_ignore_shared_cols,.fn=function(x){paste0(x,"_fusioncatcher")})
      
      cohort_predictions = rbind_no_colmatch(cohort_predictions,as.data.frame(patient_fusioncatcher_predictions))
      
    }
  }
  
  if(length(Sys.glob(merged_fusion_level_results_path))==1) { 
    fusion_level_results = read.table(merged_fusion_level_results_path,header=T,sep="\t") 
    fusion_level_results$patient_id = patient$patient_id
    cohort_fusion_level_results = rbind(cohort_fusion_level_results,as.data.frame(fusion_level_results))
    
  } else {
    
    if(length(Sys.glob(sf_fusion_level_results_path))==1) {
      patient_starfusion_fusion_level = read.table(sf_fusion_level_results_path,sep="\t",header=T)
      patient_starfusion_fusion_level$patient_id=patient$patient_id
      patient_starfusion_fusion_level = patient_starfusion_fusion_level %>% 
        harmonize_gene_identifiers() %>%
        make_identifiers_cohort_analysis()
      
      patient_starfusion_fusion_level$rna_tools = "starfusion"
      patient_starfusion_fusion_level = patient_starfusion_fusion_level %>%
        dplyr::rename_with(.cols=fusion_level_ignore_shared_cols,.fn=function(x){paste0(x,"_starfusion")})
      
      cohort_fusion_level_results = rbind_no_colmatch(cohort_fusion_level_results,as.data.frame(patient_starfusion_fusion_level))
      
    }
    
    if(length(Sys.glob(fc_fusion_level_results_path))==1) {
      patient_fusioncatcher_fusion_level = read.table(fc_fusion_level_results_path,sep="\t",header=T)
      patient_fusioncatcher_fusion_level$patient_id=patient$patient_id
      patient_fusioncatcher_fusion_level = patient_fusioncatcher_fusion_level %>% 
        harmonize_gene_identifiers() %>%
        make_identifiers_cohort_analysis()
      
      patient_fusioncatcher_fusion_level$rna_tools = "fusioncatcher"
      patient_fusioncatcher_fusion_level = patient_fusioncatcher_fusion_level %>% 
        dplyr::rename_with(.cols=fusion_level_ignore_shared_cols,.fn=function(x){paste0(x,"_fusioncatcher")})
      
      cohort_fusion_level_results = rbind_no_colmatch(cohort_fusion_level_results,as.data.frame(patient_fusioncatcher_fusion_level))
      
    }
    
  }
}

## TODO reannotate location precise is combine wgs
cohort_fusion_level_results = cohort_fusion_level_results %>% 
  mutate(location_precise_fusioncatcher = grepl(paste(proximate_bp, collapse = "|"),gup_location_fusioncatcher) & 
           grepl(paste(proximate_bp, collapse = "|"),gdw_location_fusioncatcher),
         location_precise_starfusion = grepl(paste(proximate_bp, collapse = "|"),gup_location_starfusion) & 
           grepl(paste(proximate_bp, collapse = "|"),gdw_location_starfusion),
         precise_confident_fusioncatcher = location_precise_fusioncatcher & grepl(",",tools),
         precise_confident_starfusion = location_precise_starfusion & grepl(",",tools))

#check cohort_fusion_level_results %>% filter(rna_tools=="starfusion" & precise_confident_starfusion==T & !precise_confident)

cohort_fusion_level_results = cohort_fusion_level_results %>% 
  mutate(precise_confident = (precise_confident_fusioncatcher==T|precise_confident_starfusion==T),
         precise_confident = ifelse(is.na(precise_confident),F,precise_confident))


#to check
cohort_fusion_level_results %>% filter(is.na(precise_confident)|is.na(precise_confident_starfusion)|is.na(precise_confident_fusioncatcher)) %>% 
  select(precise_confident_fusioncatcher,precise_confident_starfusion,precise_confident) %>% unique()

cohort %>% filter(!patient_id %in% cohort_fusion_level_results$patient_id)

write.table(cohort_fusion_level_results,cohort_fusion_level_results_path,quote = FALSE,sep = "\t",row.names=FALSE)
write.table(cohort_predictions,fusion_overview_path,quote = FALSE,sep = "\t",row.names=FALSE)


## Collects SVs ----

supporting_svs_df = data.frame(stringsAsFactors = FALSE)
pairwise_overlap_merged_df = data.frame(stringsAsFactors = FALSE)
cohort_results = data.frame(stringsAsFactors=FALSE) 

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
    
    pairwise_overlap_merged_path = stri_replace_all_fixed(pairwise_overlap_merged_path_template ,names(map_template_vars_patient), map_template_vars_patient,vectorize=F)
    supporting_svs_path = stri_replace_all_fixed(supporting_svs_path_template,names(map_template_vars_patient), map_template_vars_patient,vectorize=F)
    matching_bp_path = stri_replace_all_fixed(matching_bp_path_template,names(map_template_vars_patient), map_template_vars_patient,vectorize=F)
    
    
    # to be able to annotate with overlap variable
      
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
      
      if(length(Sys.glob(matching_bp_path))!=1){ next()}
      matching_bps = read.table(matching_bp_path,header=T,sep="\t") 
      if(nrow(matching_bps)>0) {
        matching_bps$patient_id = patient$patient_id 
        matching_bps$rna_tools=analysis_type
        cohort_results = rbind_no_colmatch(cohort_results,as.data.frame(matching_bps))
      }
  }
}

#are all collected?
#check 
#cohort %>% filter(!patient_id %in% supporting_svs_df$patient_id)
#cohort %>% filter(!patient_id %in% cohort_results$patient_id)

supporting_svs_df$patient_sv_name = paste0(supporting_svs_df$patient_id,"_",supporting_svs_df$sv_name)

write.table(supporting_svs_df,cohort_supporting_svs_path,quote = F,sep = "\t",row.names=F,col.names = T)

write.table(pairwise_overlap_merged_df,cohort_pairwise_overlap_merged_path,quote = F,sep = "\t",row.names=F,col.names = T)


cohort_results = cohort_results  %>% 
  harmonize_gene_identifiers() %>%
  make_identifiers_cohort_analysis()

## TODO check if necessary
cohort_results = cohort_results %>% dplyr::rename(tool= gup_tool) %>% select(-gdw_tool)

write.table(cohort_results,cohort_results_path,quote = FALSE,sep = "\t",row.names=FALSE)

