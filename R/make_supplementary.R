## Fusion-sq 
## Make supplementary files 
## Uses input from Annotate cohort report
## Last update: 2021-08

## Output
# unique fusions and unique gene pairs (distinct fusions)
# patient-oriented overview
# recurrence analyses

# Changelog
## 2022-02-26
##  path local overrides and docker version
## goal is to get uq fusions / uq gene pairs df that also selects best rna-supported prediction 
## to remake tables and know what sv to present


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
}

if(!exists("patient_table_path") | length(Sys.glob(patient_table_path))!=1) {
  print(paste0("Cohort file not available: ",patient_table_path))
  #quit()
}

suppressPackageStartupMessages({
  library(tidyverse, quietly=TRUE)
  library(stringr, quietly=TRUE)
  library(dplyr)
  library(stringi)
})

source(paste0(script_dir,"functions.general.R")) 

map_template_vars=c('${resources_dir}'=resources_dir,'${input_dir}'=input_dir,'${output_dir}'=output_dir,
                    '${cohort_identifier}'=cohort_identifier,'${cohort_wdir}'=cohort_wdir)

reports_dir =  stri_replace_all_fixed(reports_dir_template,names(map_template_vars), map_template_vars,vectorize=F)
base_dir =  stri_replace_all_fixed(base_dir_template,names(map_template_vars), map_template_vars,vectorize=F)

map_template_vars=c(map_template_vars,'${reports_dir}'=reports_dir,'${base_dir}'=base_dir)

#tmp dots
map_template_vars_merged = c('${analysis_type}'=".merged", map_template_vars )

map_template_vars = map_template_vars_merged

#input
cohort_report_path = stri_replace_all_fixed(cohort_report_path_template,names(map_template_vars), map_template_vars,vectorize=F)
fusion_overview_anno_path = stri_replace_all_fixed(fusion_overview_anno_path_template,names(map_template_vars), map_template_vars,vectorize=F) 

fga_anno = read.table("~/fusion_sq/fga_anno.tmp.tsv",sep="\t",header=T)

cohort = read.table(patient_table_path,sep = "\t", header=T,stringsAsFactors = F)
cohort$supergroup=cohort$domain
cohort$supergroup_label=cohort$domain_label



## Output 

##TODO templates
round_digits=2

#cohort analyses
patient_oriented_table_path = paste0(reports_dir,patient_oriented_table_outfile)
recurrence_table_path = paste0(reports_dir,recurrence_table_outfile)

#Fusion tables
uq_fusions_path = paste0(reports_dir,uq_fusions_outfile)
uq_gene_pairs_path =  paste0(reports_dir,uq_gene_pairs_outfile)

#manuscript additional files:
manuscript_table1_path = paste0(reports_dir,"table1_uq_fusions.tsv")
manuscript_table2_path =   paste0(reports_dir,"table2_distinct_fusions_uq_gene_pairs.tsv")

#other for manual checks
selected_duplicates_hc_path = paste0(reports_dir,"uq_patient_fusions.highconf.selected_duplicates.tsv")
multi_tools_not_same_sv_path=paste0(reports_dir,"multi_tools_not_same_sv.tsv")

### START 

test_data=F

## Input 
cohort_report =  read.table(cohort_report_path,header=T,sep="\t")

## na consensus =>already done in anno cohort report and fusion anno, here just candidate repoting
if(FALSE) {

## add helper vars for more sophisticated analyses
cohort_report = cohort_report %>% 
  mutate(helper_reciprocal_fusion_name = paste0(gdw_gene_id,"--",gup_gene_id),
         helper_reciprocal_patient_fusion = paste0(patient_id,"_",helper_reciprocal_fusion_name),
         helper_reciprocal_patient_sv_id = paste0(patient_id,"_",gdw_sv_merged,"_",gup_sv_merged),
         helper_reciprocal_patient_fusion_sv = paste0(patient_id,"_",gdw_sv_merged,"_",gup_sv_merged,"_",helper_reciprocal_fusion_name))

#remove if another patient-fusion-sv is validated by both tools
helper_both = cohort_report  %>% filter(rna_tools=="fusioncatcher, starfusion")
helper_sf_only = cohort_report %>% filter(rna_tools=="starfusion" & !patient_fusion %in% helper_both$patient_fusion)
helper_fc_only = cohort_report %>% filter(rna_tools=="fusioncatcher" & !patient_fusion %in% helper_both$patient_fusion)
helper_fc_any = cohort_report %>% filter(grepl("fusioncatcher",rna_tools))
helper_sf_any = cohort_report %>% filter(grepl("starfusion",rna_tools))


# reciprocal fusion with identical SV
#always select "patient_fusion_sv %in% helper ...$patient_fusion_sv" to get unique  
reciprocal_fusion_same_sv = cohort_report %>% 
  filter( (patient_fusion_sv %in% helper_sf_only$patient_fusion_sv & 
             patient_fusion %in% helper_fc_only$helper_reciprocal_patient_fusion & 
             patient_sv_id %in% helper_fc_only$helper_reciprocal_patient_sv_id) | 
            (patient_fusion_sv %in% helper_fc_only$patient_fusion_sv & 
               patient_fusion %in% helper_sf_only$helper_reciprocal_patient_fusion & 
               patient_sv_id %in% helper_sf_only$helper_reciprocal_patient_sv_id))

reciprocal_fusion_same_sv %>% select(patient_fusion,rna_tools,helper_reciprocal_patient_fusion)

#For the ones below: by filtering out the patient_fusion_sv from above, prevent overlap between these dataframes

#both but different sv => check manually
same_fusion_different_sv = cohort_report %>% 
  filter( !patient_fusion_sv %in% reciprocal_fusion_same_sv$patient_fusion_sv) %>%
  filter( (patient_fusion_sv %in% helper_sf_only$patient_fusion_sv & 
             patient_fusion %in% helper_fc_only$patient_fusion ) | 
            (patient_fusion_sv %in% helper_fc_only$patient_fusion_sv & 
               patient_fusion %in% helper_sf_only$patient_fusion)) 

# TODO write df for user
merging_issue = same_fusion_different_sv %>%
  filter(
    (rna_tools=="fusioncatcher" & patient_fusion_sv %in% filter(same_fusion_different_sv,rna_tools=="starfusion")$patient_fusion_sv) |
      (rna_tools=="starfusion" & patient_fusion_sv %in% filter(same_fusion_different_sv,rna_tools=="fusioncatcher")$patient_fusion_sv) ) %>%
  select(patient_fusion,rna_tools,contains('breakpoint'),patient_fusion_sv,sv_names) %>% arrange(patient_fusion_sv)
# TODO write df for user

merging_issue


cohort_report = cohort_report %>% 
  mutate(rna_consensus = (patient_fusion_sv %in% c(helper_both$patient_fusion_sv,
                                                             reciprocal_fusion_same_sv$patient_fusion_sv,
                                                             reciprocal_fusion_same_sv$helper_reciprocal_patient_fusion_sv)))

cohort_report = cohort_report %>% mutate(high_confidence_rna_dna = rna_consensus & precise_confident)
cohort_report = cohort_report %>% mutate(not_high_confidence_rna_dna = !patient_fusion %in% filter(cohort_report,high_confidence_rna_dna)$patient_fusion)

}
#### 
## uq fusions ----

if(uq_fusions(cohort_report)!=(uq_fusions(filter(cohort_report,high_confidence_rna_dna))+uq_fusions(filter(cohort_report,not_high_confidence_rna_dna)))) {
  print("Warning: high and low conf sets could not be split")
} else {
  cohort_report_hc = annotate_variant_class_fractions(filter(cohort_report,high_confidence_rna_dna))
  cohort_report_lc = annotate_variant_class_fractions(filter(cohort_report,not_high_confidence_rna_dna))
  cohort_report = rbind(cohort_report_hc,cohort_report_lc) %>% unique()
}

## Make uq patient fusion ----
make_uq_patient_fusion_df = function(cohort_report) {
  
  uq_fusions_df = cohort_report
  ### Duplicate rows with same patient_fusion 
  #Impact: will not affects counts but can affect tables and SV properties
  
  #most relevant for the high conf WGS & RNA 
  ## for WGS:  multiple tools, precise location  (precise_confident from starfusion OR fusioncatcher)
  ### previous version with single RNA tool has code to select by tumor-AF
  ## for RNA: same SV supporting same gene pair of two tools, joining predictions by the patient_fusion_sv match
  ### in some cases multiple combinations => if both tools have a precise and less precise combi then can get 4 rows of those different combinations TT,FF,TF,FT.
  
  duplicate_fusions = uq_fusions_df %>% filter(patient_fusion %in% uq_fusions_df[duplicated(uq_fusions_df$patient_fusion),]$patient_fusion)
  
  #Selection:
  ## 1) if exists precise_confident by both RNA tools 
  ## 2) highest tumor AF 
  ## Chose for selecting single one instead of merging AF/length because there was a reason they were not merged before
  
  both_rna_tools_precise_confident = duplicate_fusions %>% filter(precise_confident_fusioncatcher&precise_confident_starfusion) %>% 
    group_by(patient_fusion) %>% summarize(max_tumor_af = max(tumor_af_mean,na.rm=T))
  
  single_rna_tools_precise_confident = duplicate_fusions %>% filter(!patient_fusion %in% both_rna_tools_precise_confident$patient_fusion) 
  
  fc_precise_confident = single_rna_tools_precise_confident %>% filter(precise_confident_fusioncatcher) %>% 
    group_by(patient_fusion) %>% summarize(max_tumor_af = max(tumor_af_mean,na.rm=T))
  sf_precise_confident = single_rna_tools_precise_confident %>% filter(precise_confident_starfusion) %>% 
    group_by(patient_fusion) %>% summarize(max_tumor_af = max(tumor_af_mean,na.rm=T))
  
  neither_precise_confident =  duplicate_fusions %>% filter(!precise_confident_fusioncatcher&!precise_confident_starfusion) %>% 
    group_by(patient_fusion) %>% summarize(max_tumor_af = max(tumor_af_mean,na.rm=T))
   
  both_rna_tools_precise_confident_keep = uq_fusions_df %>% merge(both_rna_tools_precise_confident,by.x=c("patient_fusion","tumor_af_mean"),by.y=c("patient_fusion","max_tumor_af"))
  fc_precise_confident_keep = uq_fusions_df %>% merge(fc_precise_confident,by.x=c("patient_fusion","tumor_af_mean"),by.y=c("patient_fusion","max_tumor_af"))
  sf_precise_confident_keep = uq_fusions_df %>% merge(sf_precise_confident,by.x=c("patient_fusion","tumor_af_mean"),by.y=c("patient_fusion","max_tumor_af"))
  neither_precise_confident_keep = uq_fusions_df %>% merge(neither_precise_confident,by.x=c("patient_fusion","tumor_af_mean"),by.y=c("patient_fusion","max_tumor_af"))
  
  
  duplicate_fusions_keep = rbind(both_rna_tools_precise_confident_keep,fc_precise_confident_keep,sf_precise_confident_keep,neither_precise_confident_keep) %>% unique()
  
  
  still_duplicated = duplicate_fusions_keep %>% filter(patient_fusion %in% duplicate_fusions_keep[duplicated(duplicate_fusions_keep$patient_fusion),]$patient_fusion)
  
  if(nrow(still_duplicated)!=0) {
    print("WARNING duplicates not removed")
    #print(still_duplicated %>% select(patient_fusion_sv,contains('precise'),tumor_af_mean,contains('location')))
  #TODO: 
    #location: intron-intron takes priority over intron-consensus I thnk. 
    #composite cannot be resolved
  }
  
  uq_fusions_df = uq_fusions_df %>% filter(!patient_fusion %in% duplicate_fusions$patient_fusion)
  uq_fusions_df=rbind(duplicate_fusions_keep,uq_fusions_df)
  
  
  #check if still the same
  if(uq_fusions_df %>% uq_fusions() != cohort_report %>% uq_fusions()) {
    print("ERROR fusions not identical numbers as start, returning original df")
    return(cohort_report)
  }
  
  #reannotate, ambiguous should be removed now
  uq_fusions_df = annotate_variant_class_fractions(uq_fusions_df)
  
  if(nrow(filter(uq_fusions_df,ambiguous))!=0) {  
    print("WARNING ambiguous not removed")
   # print(filter(uq_fusions_df,ambiguous))
  }
  
  
  #after selecting higest tumor AF only one patient-fusion-sv remains
  
  ## Apply labels
  uq_fusions_df = annotate_labels_variant_type(uq_fusions_df)
  uq_fusions_df = annotate_labels_cancer_common(uq_fusions_df)
  
  ##Adjust sv type if complex
  uq_fusions_df = uq_fusions_df %>% mutate(svtype_label = ifelse(grepl(", ",svtype),"complex",svtype))
  uq_fusions_df$svtype_label = factor(uq_fusions_df$svtype_label)
  
  return(uq_fusions_df)
}

if(FALSE) {
uq_fusions_df$helper_patient_fusion_sv_rna = paste0(uq_fusions_df$patient_fusion_sv,"_",
                                                    uq_fusions_df$fusion_predictions_fusioncatcher,"_",
                                                    uq_fusions_df$fusion_predictions_starfusion)

## TODO problem: 
uq_fusions_df[duplicated(uq_fusions_df$helper_patient_fusion_sv_rna),]

uq_fusions_hc %>% filter(grepl("IGH",patient_fusion)) 
}


uq_fusions_hc = make_uq_patient_fusion_df(filter(cohort_report,high_confidence_rna_dna))
uq_fusions_lc = make_uq_patient_fusion_df(filter(cohort_report,not_high_confidence_rna_dna))



if(!test_data) {
  #manual steps
  ## add complex case of 
  add_fusion_manually = cohort_report[cohort_report$fusion_name =="ASPSCR1--TFE3",]
  ## Apply labels
  add_fusion_manually = annotate_labels_variant_type(add_fusion_manually)
  add_fusion_manually = annotate_labels_cancer_common(add_fusion_manually)
  
  ##Adjust sv type if complex
  add_fusion_manually = add_fusion_manually %>% mutate(svtype_label = ifelse(grepl(", ",svtype),"complex",svtype))
  add_fusion_manually$svtype_label = factor(add_fusion_manually$svtype_label)
  add_fusion_manually$high_confidence=F
}

if(nrow(uq_fusions_hc)>0){
  uq_fusions_hc$high_confidence_rna_dna = TRUE
  uq_fusions_hc$high_confidence = TRUE
  
}
if(nrow(uq_fusions_lc)>0){
  uq_fusions_lc$high_confidence_rna_dna = FALSE
  uq_fusions_lc$high_confidence = FALSE
  
}

uq_fusions_df = rbind(uq_fusions_hc,uq_fusions_lc)

if(!test_data){
uq_fusions_df = rbind(uq_fusions_df,add_fusion_manually)
}

uq_fusions_df = annotate_labels_genes(uq_fusions_df)
uq_fusions_df = annotate_labels_cancer_common(uq_fusions_df)

uq_fusions_df = unique(uq_fusions_df)

uq_fusions_df = uq_fusions_df %>% mutate(anno_clinically_relevant = clinically_validated | clinically_validated_reciprocal)

selected_duplicates_hc = cohort_report_hc %>% filter(patient_fusion %in% cohort_report_hc[duplicated(cohort_report_hc$patient_fusion),c("patient_fusion")]) %>%
  mutate(selected = patient_fusion_sv %in% uq_fusions_hc$patient_fusion_sv) %>%
  select(patient_fusion,tumor_af_mean,normal_af_mean,svtype,tools,selected) 

write.table(uq_fusions_df,uq_fusions_path,quote=F,row.names = F,col.names = T,sep="\t")
write.table(selected_duplicates_hc,selected_duplicates_hc_path,quote=F,row.names = F,col.names = T,sep="\t")


## Make unique gene pairs df ----

## distinct fusions
## for HC from the uq fusions hc (selected fusions), and from all for LC fusions
## Flag to distinguish high_confidence=T set  

make_uq_gene_pairs = function(uq_fusions_df) {
  
  uq_gene_pairs_df = uq_fusions_df  %>% 
    group_by(fusion_name,somatic_variant,germline_variant,low_af) %>% 
    summarize(svtype = toString(unique(sort(unlist(strsplit(svtype,", "))))),
              tumor_af = mean(tumor_af_mean,na.rm = T),
              normal_af = mean(normal_af_mean,na.rm = T),
              svlen = mean(svlen,na.rm=T),
              patient_ids=toString(unique(sort(patient_id))),
              patient_cnt=length(unique(patient_id)),
              .groups="keep")
  
  ## duplicate fusions are labelled as ambiguous
  duplicate_fusions = uq_gene_pairs_df[duplicated(uq_gene_pairs_df$fusion_name),]
  duplicate_fusions = uq_fusions_df %>% filter(fusion_name %in% duplicate_fusions$fusion_name) %>% 
    group_by(fusion_name) %>% 
    summarize(svtype = toString(unique(sort(unlist(strsplit(svtype,", "))))),
              tumor_af = mean(tumor_af_mean,na.rm = T),
              normal_af = mean(normal_af_mean,na.rm = T),
              svlen = mean(svlen,na.rm=T),
              patient_ids=toString(unique(sort(patient_id))),
              patient_cnt=length(unique(patient_id)),
              .groups="keep")
  
  duplicate_fusions = duplicate_fusions %>% mutate(somatic_variant=F,germline_variant=F,low_af=F)
  uq_gene_pairs_df = uq_gene_pairs_df %>% filter(!fusion_name %in% duplicate_fusions$fusion_name)
  uq_gene_pairs_df=rbind(duplicate_fusions,uq_gene_pairs_df)
  
  if(nrow(uq_gene_pairs_df[duplicated(uq_gene_pairs_df$fusion_name),])!=0) {
    print("WARNING duplicates not removed")
    print(uq_gene_pairs_df[duplicated(uq_gene_pairs_df$fusion_name),])
  }
  
  
  ## Apply labels
  uq_gene_pairs_df = annotate_labels_variant_type(uq_gene_pairs_df)
  
  ## add annotation after grouping
  annotation_cols = c("anno_healthy_chimera","anno_healthy_chimera_all","anno_cancer_chimera","anno_has_onco_or_tsg","anno_clinically_relevant","
                      anno_has_oncogene","anno_has_tsg","anno_has_kinase","anno_cancer_gene_db")
  annotation_cols = annotation_cols[annotation_cols %in% names(uq_fusions_df)]
  fusion_name_labels = uq_fusions_df %>% select(fusion_name, all_of(annotation_cols)) %>% unique()
  uq_gene_pairs_df = uq_gene_pairs_df %>% left_join(fusion_name_labels,by="fusion_name")
  
  #add sv level annotation
  fusions_svs_map = uq_fusions_df %>% select(fusion_name,patient_fusion_sv,anno_sv_population)
  uq_gene_pairs_df = uq_gene_pairs_df %>% mutate(anno_sv_population = fusion_name %in% 
                                                   filter(fusions_svs_map,anno_sv_population)$fusion_name)
  
  #will not add clinically validated annotation since not present
  uq_gene_pairs_df = annotate_labels_cancer_common(uq_gene_pairs_df)
  
  ##Adjust sv type if complex
  uq_gene_pairs_df = uq_gene_pairs_df %>% mutate(svtype_label = ifelse(grepl(", ",svtype),"complex",svtype))
  uq_gene_pairs_df$svtype_label = factor(uq_gene_pairs_df$svtype_label)
  
  return(uq_gene_pairs_df)
}

uq_gene_pairs_hc = make_uq_gene_pairs(uq_fusions_df %>% filter(high_confidence_rna_dna))
uq_gene_pairs_lc = make_uq_gene_pairs(filter(cohort_report,not_high_confidence_rna_dna))

#these sets can overlap, but should not be the same patients
#also see later df multi-tool support
fusions_hc_and_lc = intersect(uq_gene_pairs_hc$fusion_name,uq_gene_pairs_lc$fusion_name)


uq_gene_pairs_hc$high_confidence = TRUE
uq_gene_pairs_hc$high_confidence_rna_dna = TRUE
uq_gene_pairs_lc$high_confidence = FALSE
uq_gene_pairs_lc$high_confidence_rna_dna = FALSE

uq_gene_pairs_df = rbind(uq_gene_pairs_hc,uq_gene_pairs_lc)

if(uq_uq_fusions(cohort_report)!=length(unique(uq_gene_pairs_df$fusion_name))) {
  print("WARNING: different fusion counts between cohort report and uq gene pairs df")
  print(paste0(uq_uq_fusions(cohort_report)," cohort report, and ",length(unique(uq_gene_pairs_df$fusion_name))," uq gene pairs df"))
}

#mutate clinrel/fusion positve patient not part of function
# NOTE: This flags all fusions that occur in clinrel, not necessarily the right one!
#uq_gene_pairs_df = uq_gene_pairs_df %>% mutate(patient_clinrel_fusion = fusion_name %in%
#                                                 filter(cohort_report,patient_clinrel_fusion)$fusion_name)

uq_gene_pairs_df = uq_gene_pairs_df %>% mutate(patient_clinrel_fusion = grepl(paste0(filter(cohort,fusion_status)$patient_id,collapse = "|"),patient_ids))

uq_gene_pairs_df = uq_gene_pairs_df %>% mutate(clinically_validated = 
                                                 fusion_name %in% filter(cohort_report,clinically_validated|clinically_validated_reciprocal)$fusion_name,
                                               anno_clinically_relevant = clinically_validated)
## redo to get the clin.rel too
uq_gene_pairs_df = annotate_labels_cancer_common(uq_gene_pairs_df)

write.table(uq_gene_pairs_df,uq_gene_pairs_path,quote=F,row.names = F,col.names = T,sep="\t")

if(FALSE) {

## few fusions are identified confidently in some patients and not in others. 
## reasons: single tool, multi tools but bp did not overlap often in SD/ALU for patients low confidence
#or inexact location like for the LINC01237--LINC01881
#AC073529.1--MID1 multiple tools for PMCID894AAJ but indersects with L1
#AL162231.2--AL589645.1 as low-AF and segdup  => tumo-specific / low af 
## most are germline

print("Few fusions are identified confidently in some patients and not in others")
cohort_report %>% filter(fusion_name %in% fusions_hc_and_lc) %>% 
  select(patient_id,fusion_name,tools,repeat_family,segdup,tumor_af_mean,normal_af_mean,gup_location,gdw_location) %>%
  unique() %>% arrange(fusion_name) #%>% View()


#Multi tool not precise confident: due to non-exact location (fullgene,genebody), composite we cannot resolve automatically, or issues calling accurate breakpoints in repeat/SD regions, or due to CNA instability. Also some very low tumor AF which can throw the mapping off? 

multi_tools_not_same_sv = cohort_report %>% filter(grepl(",",tools_any_wgs)&!grepl(",",tools)&not_precise_confident) %>%
  group_by(patient_id,fusion_name,svtype,svlen,tools,anno_sv_population,anno_healthy_chimera,gup_location,gdw_location) %>%
  summarize(repeat_family=toString(repeat_family),seg_dup=toString(segdup),
            gup_copy_ratio_l2fc_mean=mean(gup_copy_ratio_l2fc_mean,na.rm=T),
            gdw_copy_ratio_l2fc_mean=mean(gdw_copy_ratio_l2fc_mean,na.rm=T),
            gup_start_copy_ratio_l2fc=mean(gup_start_copy_ratio_l2fc,na.rm=T),
            gdw_end_copy_ratio_l2fc=mean(gdw_end_copy_ratio_l2fc,na.rm=T),
            tumor_af_mean=mean(tumor_af_mean),normal_af_mean=mean(normal_af_mean)) #%>% View()

write.table(multi_tools_not_same_sv,multi_tools_not_same_sv_path,col.names=T,quote=F,row.names = F,sep="\t")

}

### Patient count table ----
## Make patient count table of fusion burden and selected fusions per patient 
fusion_overview =  read.table(fusion_overview_anno_path,header=T,sep="\t")

#prev code
if(FALSE) {
patient_fusion_cnt_table = cnt_fusions(fusion_overview,"predicted_cnt") %>% 
  left_join(cnt_fusions(fusion_overview %>% filter(FFPM>0.1),"predicted_ffpm_cnt"),by="patient_id") %>%
  left_join(cnt_fusions(cohort_report,"validated_cnt"),by = "patient_id") %>%
  left_join(cnt_fusions(cohort_report %>% filter(precise_confident),"precise_confident_cnt"),by = "patient_id") %>%
  left_join(cnt_fusions(cohort_report %>% filter(precise_confident&somatic_variant_only),"somatic_hc_cnt"),by = "patient_id") %>%
  left_join(cnt_fusions(cohort_report %>% filter(precise_confident&low_af_only),"low_af_hc_cnt"),by = "patient_id") %>%
  left_join(cnt_fusions(cohort_report %>% filter(precise_confident&germline_variant_only),"germline_hc_cnt"),by = "patient_id") %>%
  left_join(cnt_fusions(cohort_report %>% filter(precise_confident&ambiguous),"ambiguous_hc_cnt"),by = "patient_id") %>%
  left_join(cnt_fusions(cohort_report %>% filter(precise_confident&somatic_variant_only&(anno_has_oncogene|anno_has_tsg)),"somatic_hc_oncotsg_cnt"),by = "patient_id") %>%
  
  left_join(cnt_fusions(cohort_report %>% filter(not_precise_confident),"not_precise_confident_cnt"),by = "patient_id") %>%
  left_join(cnt_fusions(cohort_report %>% filter(not_precise_confident&somatic_variant_only),"somatic_not_hc_cnt"),by = "patient_id") %>%
  left_join(cnt_fusions(cohort_report %>% filter(not_precise_confident&low_af_only),"low_af_not_hc_cnt"),by = "patient_id") %>%
  left_join(cnt_fusions(cohort_report %>% filter(not_precise_confident&germline_variant_only),"germline_not_hc_cnt"),by = "patient_id") %>%
  left_join(cnt_fusions(cohort_report %>% filter(not_precise_confident&ambiguous),"ambiguous_not_hc_cnt"),by = "patient_id") %>%
  
  left_join(cnt_fusions(cohort_report %>% filter(gup_location=="intron"&gdw_location=="intron"&specific_sv),"validated_intron_cnt"),by = "patient_id") %>%
  left_join(cnt_fusions(cohort_report %>% filter(somatic_variant),"validated_somatic_cnt"),by = "patient_id") %>%
  left_join(cnt_fusions(cohort_report %>% filter(somatic_variant&specific_sv),"validated_somatic_specific_cnt"),by = "patient_id") %>%
  left_join(cnt_fusions(cohort_report %>% filter(somatic_variant & anno_in_any_cancer),"validated_somatic_cancer_cnt"),by = "patient_id") %>%
  left_join(cnt_fusions(cohort_report %>% filter(germline_variant),"validated_germline_cnt"),by = "patient_id") %>%
  left_join(cnt_fusions(cohort_report %>% filter(germline_variant&specific_sv),"validated_germline_specific_cnt"),by = "patient_id") %>%
  left_join(cnt_fusions(cohort_report %>% filter(germline_variant & anno_in_any_cancer),"validated_germline_cancer_cnt"),by = "patient_id") %>%
  left_join(cnt_fusions(cohort_report %>% filter(gup_location=="intron"&gdw_location=="intron"&specific_sv&anno_in_any_cancer), "validated_intron_cancer_cnt"),by = "patient_id") %>%
  left_join(cnt_fusions(cohort_report %>% filter(low_af),"validated_low_af_cnt"),by = "patient_id") %>%
  
  left_join(cnt_fusions(cohort_report %>% filter(clinically_validated|clinically_validated_reciprocal),"clinical_cnt"),by = "patient_id") %>%
  left_join(patient_metadata %>% select(patient_metadata_cols),by = "patient_id") %>%
  arrange(-predicted_cnt, -validated_cnt)

patient_fusion_cnt_table[is.na(patient_fusion_cnt_table)]=0

patient_fusion_cnt_table$patient_id = factor(patient_fusion_cnt_table$patient_id)

patient_fusion_cnt_table$primary_group = factor(patient_fusion_cnt_table$primary_group)
patient_fusion_cnt_table$primary_group_shorthand_label = factor(patient_fusion_cnt_table$primary_group_shorthand_label)

patient_fusion_cnt_table$supergroup = factor(patient_fusion_cnt_table$supergroup)
patient_fusion_cnt_table$supergroup_label = factor(patient_fusion_cnt_table$supergroup_label)

## to distinguish in plots
patient_fusion_cnt_table = patient_fusion_cnt_table %>% mutate(patient_clinrel_fusion = clinical_cnt>0)
## to display fusions not germline = low af + somatic
patient_fusion_cnt_table$somatic_low_af_hc_cnt = patient_fusion_cnt_table$somatic_hc_cnt+patient_fusion_cnt_table$low_af_hc_cnt
## to display also low confidence low af + somatic on top of the high conf low af + somatic
patient_fusion_cnt_table$any_somatic_low_af_cnt = patient_fusion_cnt_table$somatic_low_af_hc_cnt + patient_fusion_cnt_table$somatic_not_hc_cnt+patient_fusion_cnt_table$low_af_not_hc_cnt
}


get_patient_fusion_cnt_table_short = function(cohort_report,fusion_overview,patient_metadata) {
  patient_metadata_cols = c("patient_id","patient_label","primary_group","primary_group_shorthand_label","supergroup","supergroup_label","fga","max_cna","min_cna","fusion_positive_patient")
  patient_metadata_cols_missing=patient_metadata_cols[!patient_metadata_cols %in% names(patient_metadata)]
  patient_metadata_cols=patient_metadata_cols[patient_metadata_cols %in% names(patient_metadata)]
  
  patient_fusion_cnt_table = cnt_fusions(fusion_overview,"predicted_cnt") %>% 
    left_join(cnt_fusions(fusion_overview %>% filter(FFPM>0.1),"predicted_ffpm_cnt"),by="patient_id") %>%
    left_join(cnt_fusions(cohort_report,"validated_cnt"),by = "patient_id") %>%
    left_join(cnt_fusions(cohort_report %>% filter(precise_confident),"precise_confident_cnt"),by = "patient_id") %>%
    left_join(cnt_fusions(cohort_report %>% filter(precise_confident&somatic_variant_only),"somatic_hc_cnt"),by = "patient_id") %>%
    left_join(cnt_fusions(cohort_report %>% filter(precise_confident&low_af_only),"low_af_hc_cnt"),by = "patient_id") %>%
    left_join(cnt_fusions(cohort_report %>% filter(precise_confident&germline_variant_only),"germline_hc_cnt"),by = "patient_id") %>%
    left_join(cnt_fusions(cohort_report %>% filter(precise_confident&ambiguous),"ambiguous_hc_cnt"),by = "patient_id") %>%
    left_join(cnt_fusions(cohort_report %>% filter(precise_confident&somatic_variant_only&(anno_has_oncogene|anno_has_tsg)),"somatic_hc_oncotsg_cnt"),by = "patient_id") %>%
    
    left_join(cnt_fusions(cohort_report %>% filter(not_precise_confident),"not_precise_confident_cnt"),by = "patient_id") %>%
    left_join(cnt_fusions(cohort_report %>% filter(not_precise_confident&somatic_variant_only),"somatic_not_hc_cnt"),by = "patient_id") %>%
    left_join(cnt_fusions(cohort_report %>% filter(not_precise_confident&low_af_only),"low_af_not_hc_cnt"),by = "patient_id") %>%
    left_join(cnt_fusions(cohort_report %>% filter(not_precise_confident&germline_variant_only),"germline_not_hc_cnt"),by = "patient_id") %>%
    left_join(cnt_fusions(cohort_report %>% filter(not_precise_confident&ambiguous),"ambiguous_not_hc_cnt"),by = "patient_id") %>%
    
    left_join(cnt_fusions(cohort_report %>% filter(clinically_validated|clinically_validated_reciprocal),"clinical_cnt"),by = "patient_id") %>%
    left_join(patient_metadata %>% select(patient_metadata_cols),by = "patient_id") %>%
    arrange(-predicted_cnt, -validated_cnt)
  
  patient_fusion_cnt_table[is.na(patient_fusion_cnt_table)]=0
  
  patient_fusion_cnt_table$patient_id = factor(patient_fusion_cnt_table$patient_id)
  
  patient_fusion_cnt_table$primary_group = factor(patient_fusion_cnt_table$primary_group)
  patient_fusion_cnt_table$primary_group_shorthand_label = factor(patient_fusion_cnt_table$primary_group_shorthand_label)
  
  patient_fusion_cnt_table$supergroup = factor(patient_fusion_cnt_table$supergroup)
  patient_fusion_cnt_table$supergroup_label = factor(patient_fusion_cnt_table$supergroup_label)
  
  ## to distinguish in plots
  patient_fusion_cnt_table = patient_fusion_cnt_table %>% mutate(patient_clinrel_fusion = clinical_cnt>0)
  ## to display fusions not germline = low af + somatic
  patient_fusion_cnt_table$somatic_low_af_hc_cnt = patient_fusion_cnt_table$somatic_hc_cnt+patient_fusion_cnt_table$low_af_hc_cnt
  ## to display also low confidence low af + somatic on top of the high conf low af + somatic
  patient_fusion_cnt_table$any_somatic_low_af_cnt = patient_fusion_cnt_table$somatic_low_af_hc_cnt + patient_fusion_cnt_table$somatic_not_hc_cnt+patient_fusion_cnt_table$low_af_not_hc_cnt
  
  return(patient_fusion_cnt_table)
}


patient_fusion_cnt_table = get_patient_fusion_cnt_table_short(filter(cohort_report,rna_consensus),filter(fusion_overview,rna_consensus), cohort)
patient_fusion_cnt_table = patient_fusion_cnt_table %>% left_join(fga_anno)

## Add names of  fusions 

patient_oriented_table = patient_fusion_cnt_table

## TODO later
if(FALSE) {

patient_oriented_table = patient_oriented_table %>% left_join(cohort_report %>% filter( clinically_validated | clinically_validated_reciprocal)  %>%
                                                                group_by(patient_id,primary_group_shorthand_label) %>% 
                                                                summarize(clinrel_fusions = toString(unique(sort(fusion_name)))))
patient_oriented_table = patient_oriented_table %>% left_join(cohort_report %>% filter( somatic_variant & precise_confident & anno_cancer_chimera) %>% 
                                                                group_by(patient_id,primary_group_shorthand_label) %>% 
                                                                summarize(somatic_cancer_chimera = toString(unique(sort(fusion_name))),
                                                                          anno_chimerseq_cancer_types=toString(anno_chimerseq_cancer_types)))
patient_oriented_table = patient_oriented_table %>% left_join(cohort_report %>%
                                                                filter( somatic_variant & precise_confident & (anno_has_oncogene | anno_has_tsg) ) %>% 
                                                                group_by(patient_id,primary_group_shorthand_label) %>% 
                                                                summarize(somatic_onco_tsg_fusions = toString(unique(sort(fusion_name)))))
patient_oriented_table = patient_oriented_table %>% left_join(cohort_report %>% filter( precise_confident & anno_cancer_gene_db & !(anno_has_oncogene | anno_has_tsg) ) %>% 
                                                                group_by(patient_id,primary_group_shorthand_label) %>% 
                                                                summarize(other_cancer_fusions = toString(unique(sort(fusion_name)))))


patient_oriented_table = patient_oriented_table %>% left_join(cohort_report %>%
                                                                filter( somatic_variant & !precise_confident & (anno_cancer_gene_db | anno_cancer_chimera) ) %>% 
                                                                group_by(patient_id,primary_group_shorthand_label) %>% 
                                                                summarize(somatic_not_hc_cancer_fusions = toString(unique(sort(fusion_name)))))


## max CNA of high conf fusion 

if("gup_copy_ratio_l2fc_mean" %in% names(cohort_report)){
max_cna_fusion_hc = uq_fusions_hc %>% group_by(patient_id) %>% summarize(
  max_cna_fusion=max(c(gup_copy_ratio_l2fc_mean,gdw_copy_ratio_l2fc_mean,gup_start_copy_ratio_l2fc,gup_end_copy_ratio_l2fc,gdw_start_copy_ratio_l2fc,gdw_end_copy_ratio_l2fc),na.rm = T))
max_cna_fusion_lc = cohort_report_lc %>% filter(!patient_id %in% max_cna_fusion_hc$patient_id) %>% group_by(patient_id) %>% summarize(  
  max_cna_fusion=max(c(gup_copy_ratio_l2fc_mean,gdw_copy_ratio_l2fc_mean,gup_start_copy_ratio_l2fc,gup_end_copy_ratio_l2fc,gdw_start_copy_ratio_l2fc,gdw_end_copy_ratio_l2fc),na.rm = T))

max_cna_fusion = rbind(max_cna_fusion_hc,max_cna_fusion_lc)
max_cna_fusion[is.infinite(max_cna_fusion$max_cna_fusion),c("max_cna_fusion")]=0

#check uq_patients(max_cna_fusion)==uq_patients(cohort_report)

patient_oriented_table = patient_oriented_table %>% left_join(max_cna_fusion,by = "patient_id")
patient_oriented_table[is.na(patient_oriented_table$max_cna_fusion),c("max_cna_fusion")]=0

}
}

patient_oriented_table %>% keep(is.numeric) %>% colSums()
patient_oriented_table %>% filter(patient_clinrel_fusion) %>% keep(is.numeric) %>% colSums()

write.table(patient_oriented_table,patient_oriented_table_path,quote=F,row.names = F,sep="\t")


## Recurrence tables ----

recurrence_table = fusion_cnt_per_attr(filter(fusion_overview,rna_consensus),"patient_id","chimeric_transcript")  %>%
  left_join(fusion_cnt_per_attr(filter(cohort_report,rna_consensus),"patient_id","any_wgs"),by="fusion_name") %>%
  filter(chimeric_transcript>1|any_wgs>1) %>%
  left_join(fusion_cnt_per_attr(filter(cohort_report,rna_consensus&precise_confident),"patient_id","high_confidence"),by="fusion_name") %>%
  left_join(fusion_cnt_per_attr(filter(cohort_report,rna_consensus&somatic_variant),"patient_id","any_tumor_specific"),by="fusion_name") %>%
  left_join(fusion_cnt_per_attr(filter(cohort_report,rna_consensus&somatic_variant&precise_confident),"patient_id","high_confidence_tumor_specific"),by="fusion_name") %>%
  left_join(fusion_cnt_per_attr(filter(cohort_report,rna_consensus&germline_variant),"patient_id","any_germline"),by="fusion_name") %>%
  left_join(fusion_cnt_per_attr(filter(cohort_report,rna_consensus&germline_variant&precise_confident),"patient_id","high_confidence_germline"),by="fusion_name") %>%

  arrange(-high_confidence)

recurrence_table[is.na(recurrence_table)]=0

write.table(recurrence_table,recurrence_table_path,quote=F,row.names = F,sep="\t")


### MANUSCRIPT exports ### ----

columns_export = c("fusion_name","patient_label","variant_type",
                   "anno_has_onco_or_tsg","anno_clinically_relevant",
                   "svtype","svlen","ffpm_max",
                   "tools",
                   "tumor_af_mean","normal_af_mean","gup_cytoband","gdw_cytoband",
                   "gup_label","gdw_label","annotation",
                   "tumor_type_label","primary_group_shorthand_label","supergroup_label",
                   "gup_fpkm_zscore_supergroup","gdw_fpkm_zscore_supergroup","gup_fpkm_supergroup_pval","gdw_fpkm_supergroup_pval",
                   "gup_copy_ratio_l2fc_mean","gdw_copy_ratio_l2fc_mean",
                   "gup_start_copy_ratio_l2fc","gdw_end_copy_ratio_l2fc",
                   "gup_location","gdw_location",
                   "high_confidence",
                   "gup_sv_merged_coordinate","gdw_sv_merged_coordinate",
                   "gup_sf_breakpoint","gdw_sf_breakpoint","gup_sf_transcript","gdw_sf_transcript",
                   "gup_sv_intron","gdw_sv_intron","gup_sv_intron_tx","gdw_sv_intron_tx",
                   "anno_has_kinase","anno_cancer_chimera","anno_sv_population","anno_healthy_chimera",
                   "anno_chimerseq_cancer_types",
                   "repeat_family","segdup",
                   "patient_high_fusion_burden","patient_clinrel_fusion",
                   "gup_ensembl_id", "gdw_ensembl_id",
                   "pairwise_distance_max","pairwise_overlap_mean",
                   "sv_gup_filter","sv_gdw_filter",
                   "fusion_predictions",
                   "tx_gup_transcript_id","tx_gup_involved_fragment","tx_gup_involved_fraction_mean",
                   "tx_gdw_transcript_id","tx_gdw_involved_fragment","tx_gdw_involved_fraction_mean"
                   )

uq_fusions_df = read.table(uq_fusions_path,header=T,sep="\t")
uq_gene_pairs_df = read.table(uq_gene_pairs_path,header=T,sep="\t")
patient_fusion_cnt_table = read.table(patient_oriented_table_path,header=T,sep="\t")
patients_high_fusion_burden = patient_fusion_cnt_table[patient_fusion_cnt_table$somatic_low_af_hc_cnt>=5,c("patient_id")]
cohort_report = read.table(cohort_report_path,header=T,sep="\t")

if(!test_data) {
  #manual steps
  complex_fusions = c("ASPSCR1--TFE3","ETV6--IGL-@-ext","LINC01344--TERT")
  add_fusion_manually = cohort_report[cohort_report$fusion_name %in% complex_fusions,]
  add_fusion_manually = annotate_labels_variant_type(add_fusion_manually)
  add_fusion_manually = annotate_labels_cancer_common(add_fusion_manually)
  add_fusion_manually = add_fusion_manually %>% mutate(svtype_label = ifelse(grepl(", ",svtype),"complex",svtype))
  add_fusion_manually$svtype_label = factor(add_fusion_manually$svtype_label)
  add_fusion_manually$high_confidence=F
  add_fusion_manually = annotate_labels_genes(add_fusion_manually)
  
  uq_fusions_df = rbind(uq_fusions_df,add_fusion_manually[,names(uq_fusions_df)])
}

uq_fusions_df = uq_fusions_df %>% unique()

### 2022-02 filter rna consensus first
#exclude normal and have high conf or 2 tools + composite
table1_uq_fusions = uq_fusions_df %>% filter(rna_consensus==T) %>%
  filter(((low_af|somatic_variant)&high_confidence)|
           (somatic_variant & ( (grepl(",",tools_any_wgs_fusioncatcher)&gup_location_fusioncatcher=="composite") | 
                                  (grepl(",",tools_any_wgs_starfusion)&gup_location_starfusion=="composite"))))

table1_uq_fusions = annotate_labels_cancer_common(table1_uq_fusions)

table1_uq_fusions = table1_uq_fusions %>% mutate(
  patient_high_fusion_burden = (patient_id %in% patients_high_fusion_burden))

table1_uq_fusions = cohort[,c("patient_id","patient_label")] %>% merge(table1_uq_fusions, by = "patient_id") 

columns_export=columns_export[columns_export %in% names(table1_uq_fusions)]

table1_uq_fusions_export = table1_uq_fusions %>% 
  select(all_of(columns_export)) %>%
  mutate(across(where(is.numeric), ~ round(., digits = round_digits))) %>%
  arrange(-anno_clinically_relevant,svtype,-anno_has_onco_or_tsg,fusion_name,variant_type)

table1_uq_fusions_export = unique(table1_uq_fusions_export)

write.table(table1_uq_fusions_export ,manuscript_table1_path,col.names=T,quote=T,row.names = F,sep="\t",na="")


### uq gene pairs
## distinct fusions

#to add label instead of id 
patient_mapping = table1_uq_fusions %>% group_by(fusion_name) %>%
  summarise(patients = toString(sort(unique(patient_label))),
            cancer_supergroups=toString(sort(unique(supergroup_label))))

table2_uq_gene_pairs = uq_gene_pairs_df %>% filter(fusion_name %in% table1_uq_fusions$fusion_name) %>% 
  left_join(patient_mapping)
  
columns_export_uq_pairs = names(table2_uq_gene_pairs[names(table2_uq_gene_pairs) %in% columns_export])
columns_export_uq_pairs=c("patients","cancer_supergroups",columns_export_uq_pairs)

table2_uq_gene_pairs = table2_uq_gene_pairs %>%
  mutate(across(where(is.numeric), ~ round(., digits = round_digits))) %>%
  arrange(-anno_clinically_relevant,-patient_cnt,-anno_has_onco_or_tsg,fusion_name,variant_type)

table2_uq_gene_pairs = unique(table2_uq_gene_pairs)

write.table(table2_uq_gene_pairs[,columns_export_uq_pairs] ,manuscript_table2_path,col.names=T,quote=T,row.names = F,sep="\t",na="")

