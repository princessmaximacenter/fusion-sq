## Fusion-sq 
## Make supplementary files 
## Uses input from Annotate cohort report
## Last update: 2021-08

## Output
# unique fusions and unique gene pairs
# patient-oriented overview
# recurrence analyses

suppressPackageStartupMessages({
  library(tidyverse, quietly=TRUE)
  library(stringr, quietly=TRUE)
  library(dplyr)
})

source("R/default.conf") #in script dir
source(paste0(script_dir,"functions.general.R"))

## Input
patient_metadata = read.table(patient_table,sep = "\t", header=T,stringsAsFactors = T)
patient_labels_lookup = patient_metadata[,c("patient_id","patient_label")]
patient_metadata_cols = c("patient_id","patient_label","primary_group","primary_group_shorthand_label","supergroup","supergroup_label","fga","max_cna","min_cna")
patient_metadata_cols=patient_metadata_cols[patient_metadata_cols %in% names(patient_metadata)]

cohort_report =  read.table(paste0(reports_dir,cohort_report_outfile),header=T,sep="\t")
fusion_overview =  read.table(paste0(reports_dir,fusion_overview_anno_outfile),header=T,sep="\t")

## Output
round_digits=2

#cohort analyses
patient_oriented_table_path = paste0(reports_dir,patient_oriented_table_outfile)
recurrence_table_path = paste0(reports_dir,recurrence_table_outfile)

#Fusion tables
uq_fusions_path = paste0(reports_dir,uq_fusions_outfile)
uq_gene_pairs_path =  paste0(reports_dir,uq_gene_pairs_outfile)

#manuscript additional files:
manuscript_table1_path = paste0(reports_dir,"table1_uq_fusions.tsv")
manuscript_table2_path =   paste0(reports_dir,"table2_uq_gene_pairs.tsv")

#other for manual checks
selected_duplicates_hc_path = paste0(reports_dir,"uq_patient_fusions.highconf.selected_duplicates.tsv")
multi_tools_not_same_sv_path=paste0(reports_dir,"multi_tools_not_same_sv.tsv")

### START 

## uq fusions

if(uq_fusions(cohort_report)!=(uq_fusions(filter(cohort_report,precise_confident))+uq_fusions(filter(cohort_report,not_precise_confident)))) {
  print("Warning: high and low conf sets could not be split")
} else {
  cohort_report_hc = annotate_variant_class_fractions(filter(cohort_report,precise_confident))
  cohort_report_lc = annotate_variant_class_fractions(filter(cohort_report,not_precise_confident))
  cohort_report = rbind(cohort_report_hc,cohort_report_lc)
}

## unique gene pairs dataframe is like recurrence but with more annotation and less emphasis on numbers?
#only needed for high confidence fusions = same sv

uq_fusions_hc = make_uq_patient_fusion_df(cohort_report_hc)
uq_fusions_lc = make_uq_patient_fusion_df(cohort_report_lc)

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
  uq_fusions_hc$high_confidence = TRUE
}
if(nrow(uq_fusions_lc)>0){
  uq_fusions_lc$high_confidence = FALSE
}

uq_fusions_df = rbind(uq_fusions_hc,uq_fusions_lc)

if(!test_data){
uq_fusions_df = rbind(uq_fusions_df,add_fusion_manually)
}

uq_fusions_df = annotate_labels_genes(uq_fusions_df)
uq_fusions_df = annotate_labels_cancer_common(uq_fusions_df)

uq_fusions_df = unique(uq_fusions_df)
selected_duplicates_hc = cohort_report_hc %>% filter(patient_fusion %in% cohort_report_hc[duplicated(cohort_report_hc$patient_fusion),c("patient_fusion")]) %>%
  mutate(selected = patient_fusion_sv %in% uq_fusions_hc$patient_fusion_sv) %>%
  select(patient_fusion,tumor_af_mean,normal_af_mean,svtype,tools,selected) 

write.table(uq_fusions_df,uq_fusions_path,quote=F,row.names = F,col.names = T,sep="\t")
write.table(selected_duplicates_hc,selected_duplicates_hc_path,quote=F,row.names = F,col.names = T,sep="\t")


## Make unique gene pairs df
## for HC from the uq fusions hc (selected fusions), and from all for LC fusions
## Flag to distinguish high_confidence=T set  
uq_gene_pairs_hc = make_uq_gene_pairs(uq_fusions_hc)
uq_gene_pairs_lc = make_uq_gene_pairs(cohort_report_lc)

#these sets can overlap, but should not be the same patients
#also see later df multi-tool support
fusions_hc_and_lc = intersect(uq_gene_pairs_hc$fusion_name,uq_gene_pairs_lc$fusion_name)


uq_gene_pairs_hc$high_confidence = TRUE
uq_gene_pairs_lc$high_confidence = FALSE

uq_gene_pairs_df = rbind(uq_gene_pairs_hc,uq_gene_pairs_lc)

if(uq_uq_fusions(cohort_report)!=length(unique(uq_gene_pairs_df$fusion_name))) {
  print("WARNING: different fusion counts between cohort report and uq gene pairs df")
  print(paste0(uq_uq_fusions(cohort_report)," cohort report, and ",length(unique(uq_gene_pairs_df$fusion_name))," uq gene pairs df"))
}

#mutate clinrel/fusion positve patient not part of function
uq_gene_pairs_df = uq_gene_pairs_df %>% mutate(patient_clinrel_fusion = fusion_name %in%
                                                 filter(cohort_report,patient_clinrel_fusion)$fusion_name)
uq_gene_pairs_df = uq_gene_pairs_df %>% mutate(clinically_validated = 
                                                 fusion_name %in% filter(cohort_report,clinically_validated|clinically_validated_reciprocal)$fusion_name)

## redo to get the clin.rel too
uq_gene_pairs_df = annotate_labels_cancer_common(uq_gene_pairs_df)

write.table(uq_gene_pairs_df,uq_gene_pairs_path,quote=F,row.names = F,col.names = T,sep="\t")



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


## MAKE patient count table of fusion burden and selected fusions per patient 

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


## Add names of  fusions 
patient_oriented_table = patient_fusion_cnt_table

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

write.table(patient_oriented_table,patient_oriented_table_path,quote=F,row.names = F,sep="\t")


## Recurrence tables

recurrence_table = fusion_cnt_per_attr(fusion_overview,"patient_id","chimeric_transcript")  %>%
  filter(chimeric_transcript>1) %>%
  left_join(fusion_cnt_per_attr(cohort_report,"patient_id","any_wgs"),by="fusion_name") %>%
  left_join(fusion_cnt_per_attr(filter(cohort_report,precise_confident),"patient_id","high_confidence"),by="fusion_name") %>%
  left_join(fusion_cnt_per_attr(filter(cohort_report,somatic_variant),"patient_id","any_tumor_specific"),by="fusion_name") %>%
  left_join(fusion_cnt_per_attr(filter(cohort_report,somatic_variant&precise_confident),"patient_id","high_confidence_tumor_specific"),by="fusion_name") %>%
  left_join(fusion_cnt_per_attr(filter(cohort_report,germline_variant),"patient_id","any_germline"),by="fusion_name") %>%
  left_join(fusion_cnt_per_attr(filter(cohort_report,germline_variant&precise_confident),"patient_id","high_confidence_germline"),by="fusion_name") %>%
  
#  left_join(fusion_cnt_per_attr(filter(cohort_report,somatic_variant&grepl(",",tools)),"patient_id","validated_somatic_2tools_cnt"),by="fusion_name") %>%
#  left_join(fusion_cnt_per_attr(filter(cohort_report,germline_variant&grepl(",",tools)),"patient_id","validated_germline_2tools_cnt"),by="fusion_name") %>%
  arrange(-high_confidence)

recurrence_table[is.na(recurrence_table)]=0

write.table(recurrence_table,recurrence_table_path,quote=F,row.names = F,sep="\t")


### MANUSCRIPT exports ###

columns_export = c("fusion_name","patient_label","variant_type",
                   "anno_has_onco_or_tsg","anno_clinically_relevant",
                   "svtype","svlen","ffpm_max",
                   "tools",
                   "tumor_af_mean","normal_af_mean","gup_cytoband","gdw_cytoband",
                   "gup_label","gdw_label","annotation",
                   "tumor_type_label","primary_group_shorthand_label","supergroup_label",
                   "gup_fpkm_zscore_supergroup","gdw_fpkm_zscore_supergroup","gup_fpkm_supergroup_pval","gdw_fpkm_supergroup_pval",
                   "gup_copy_ratio_l2fc_mean","gdw_copy_ratio_l2fc_mean",
                   #"gup_start_copy_ratio_l2fc","gdw_end_copy_ratio_l2fc",
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

uq_fusions_df = read.table(paste0(reports_dir,uq_fusions_outfile),header=T,sep="\t")
uq_gene_pairs_df = read.table(paste0(reports_dir,uq_gene_pairs_outfile),header=T,sep="\t")
patient_fusion_cnt_table = read.table(paste0(reports_dir,patient_oriented_table_outfile),header=T,sep="\t")
patients_high_fusion_burden = patient_fusion_cnt_table[patient_fusion_cnt_table$somatic_low_af_hc_cnt>=7,c("patient_id")]
cohort_report = read.table(paste0(reports_dir,cohort_report_outfile),header=T,sep="\t")

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

columns_export=columns_export[columns_export %in% names(uq_fusions_df)]

#exclude normal and have high conf or 2 tools + composite
table1_uq_fusions = uq_fusions_df %>% filter(((low_af|somatic_variant)&high_confidence)|
                                                        (somatic_variant&(grepl(",",tools_any_wgs)&gup_location=="composite")))

table1_uq_fusions = annotate_labels_cancer_common(table1_uq_fusions)

table1_uq_fusions = table1_uq_fusions %>% mutate(
  patient_high_fusion_burden = (patient_id %in% patients_high_fusion_burden))

table1_uq_fusions = patient_labels_lookup %>% merge(table1_uq_fusions, by = "patient_id") 

table1_uq_fusions_export = table1_uq_fusions %>% 
  select(all_of(columns_export)) %>%
  mutate(across(where(is.numeric), ~ round(., digits = round_digits))) %>%
  arrange(-anno_clinically_relevant,svtype,-anno_has_onco_or_tsg,fusion_name,variant_type)

table1_uq_fusions_export = unique(table1_uq_fusions_export)

write.table(table1_uq_fusions_export ,manuscript_table1_path,col.names=T,quote=T,row.names = F,sep="\t",na="")


### uq gene pairs

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

