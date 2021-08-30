## Fusion sq
## Gene expression analysis
## Last update: 2021-08-13
## Calculates z-scores of log2 FPKMs for full cohort (pan cancer) and per group (primary cancer type)

## Only fusion-associated and cancer-related genes are used at the moment, speeds up calculation
### PAR_Y genes are removed to prevent duplicates after version removal 

## The FPKM z score analysis needs to be rerun if: patient diagnosis changes; cohort changes

suppressPackageStartupMessages({
  library(tidyverse, quietly=TRUE)
  library(stringr, quietly=TRUE)
})
wdir="~/PycharmProjects/structuralvariation/fusion_sq/"
script_dir = paste0(wdir,"R/")
source(paste0(script_dir,"default.conf"))
source(paste0(script_dir,"functions.general.R")) 

patient_metadata = read.table(patient_table,sep = "\t", header=T,stringsAsFactors = T)


print("Running gene expression zscore analysis")
print(paste0("Patient count: ",nrow(patient_metadata)))

if(nrow(patient_metadata)<5) {
  print("Less than 5 patients, probably insufficient data. Please double check your input")
  quit()
}

gene_expression_analysis_path = paste0(reports_dir,gene_expression_zcore_outfile)

## subset to gene fusions
cohort_fusion_level_results_path = paste0(reports_dir,cohort_fusion_level_results_outfile)

##Prevent overwrite
if(length(Sys.glob(gene_expression_analysis_path))==1){
  #gene_expression_analysis = read.table(gene_expression_analysis_path,sep="\t",header = T)
  print(paste0("File already exists: ",gene_expression_analysis_path))
  #quit()
} 


## Subset genes 
if(length(Sys.glob(cohort_fusion_level_results_path))!=1){
  print(paste0("Needs cohort fusion-level results: ",cohort_fusion_level_results_path))
  #quit()
} 
cohort_fusion_report = read.table(cohort_fusion_level_results_path,sep = "\t", header=T)

cancer_genes = get_cancer_genes()
cancer_genes = cancer_genes[grepl("ENSG",cancer_genes$gene_id),c("gene_id")]
cancer_genes = remove_version_from_id(cancer_genes)

focus_gene_list = unique(c(cohort_fusion_report$gup_ensembl_id,cohort_fusion_report$gdw_ensembl_id,cancer_genes))
focus_gene_list = focus_gene_list[grepl("ENSG",focus_gene_list)]


## Make overview dataframe

  summary_gene_expression = data.frame(stringsAsFactors = F)
  for(id in patient_metadata$patient_identifier) {
    patient = filter(patient_metadata,patient_identifier==id)
    
    gene_expr_filepath = paste0(expression_data_dir,patient$rna_id,expression_gene_file_ext)
    gene_expression_df = read.table(gene_expr_filepath,sep="\t",header = F)
    colnames(gene_expression_df) = c("ensembl_id","counts","cpm","fpkm","chr","start","end","strand","length","ensembl_id_copy","gene_name","transcript_id")
    
    ##remove  PAR_Y genes before removing version 
    gene_expression_df = gene_expression_df[grepl("ENSG",gene_expression_df$ensembl_id)&!grepl("PAR_Y",gene_expression_df$ensembl_id),]
    gene_expression_df$ensembl_id = remove_version_from_id(gene_expression_df$ensembl_id)
    
    ##subset to focus genes
    if(length(focus_gene_list)>0) {
      gene_expression_df = gene_expression_df[gene_expression_df$ensembl_id %in% focus_gene_list,]
    } 
    
    gene_expression_df$patient_id = patient$patient_id
    
    summary_gene_expression = rbind(summary_gene_expression,gene_expression_df[,c("ensembl_id","counts","cpm","fpkm","patient_id")])
  }
  
  
  summary_gene_expression$ensembl_id = factor(summary_gene_expression$ensembl_id)
  summary_gene_expression$patient_id = factor(summary_gene_expression$patient_id)
  

## START z score analysis

##subset to patients, just to be sure, add metadata
summary_gene_expression = summary_gene_expression %>% filter(patient_id %in% patient_metadata$patient_id)
summary_gene_expression = summary_gene_expression %>% left_join(patient_metadata[,c("patient_id","supergroup_label","primary_group_shorthand_label")]
                                                                ,by="patient_id")

summary_gene_expression = summary_gene_expression %>% group_by(ensembl_id) %>% mutate(fpkm_log = log(fpkm+0.001))

summary_gene_expression = summary_gene_expression %>% group_by(ensembl_id) %>% mutate(fpkm_mean = mean(fpkm_log,na.rm = T), 
                                                                                      fpkm_sd = sd(fpkm_log,na.rm = T), 
                                                                                      fpkm_zscore = (fpkm_log - fpkm_mean )/ fpkm_sd) %>% ungroup() %>% unique()

gene_level_normal_distribution = summary_gene_expression %>% group_by(ensembl_id) %>% summarize(shapiro_pval = shapiro.test(fpkm_log)$p.value,
                                                                                                fpkm_log_normal_dist=(shapiro_pval>=0.05))

summary_gene_expression = summary_gene_expression %>% mutate(fpkm_log_normal_dist = 
                                                               ensembl_id %in% filter(gene_level_normal_distribution,fpkm_log_normal_dist)$ensembl_id)

#the grouped dataframe also contains the data from summary gene expression
# sum counts group is needed to check if not all patients have 0 counts which would error the shapiro test
summary_gene_expression_grouped = summary_gene_expression

summary_gene_expression_grouped = summary_gene_expression_grouped %>% group_by(ensembl_id,primary_group_shorthand_label) %>% 
  mutate(sum_counts_group = sum(counts), fpkm_mean_group = mean(fpkm_log,na.rm = T), fpkm_sd_group = sd(fpkm_log,na.rm = T), 
         fpkm_zscore_group = (fpkm_log - fpkm_mean_group )/ fpkm_sd_group) %>% ungroup()

summary_gene_expression_grouped = summary_gene_expression_grouped %>% group_by(ensembl_id,supergroup_label) %>% 
  mutate(sum_counts_supergroup = sum(counts), fpkm_mean_supergroup = mean(fpkm_log,na.rm = T), fpkm_sd_supergroup = sd(fpkm_log,na.rm = T), 
         fpkm_zscore_supergroup = (fpkm_log - fpkm_mean_supergroup )/ fpkm_sd_supergroup) %>% ungroup()

prigroup_normal_distribution = summary_gene_expression_grouped %>% ungroup() %>% group_by(ensembl_id,primary_group_shorthand_label) %>%
  summarize(primary_group_size=n(), shapiro_pval = ifelse((primary_group_size>2 & sum_counts_group>0),shapiro.test(fpkm_log)$p.value,0),
            fpkm_log_normal_dist_group=(shapiro_pval>=0.05)) %>% unique()

supergroup_normal_distribution = summary_gene_expression_grouped %>% ungroup() %>% group_by(ensembl_id,supergroup_label) %>%
  summarize(supergroup_size=n(), shapiro_pval = ifelse((supergroup_size>2 & sum_counts_supergroup>0),shapiro.test(fpkm_log)$p.value,0),
            fpkm_log_normal_dist_supergroup=(shapiro_pval>=0.05)) %>% unique()

summary_gene_expression_grouped = summary_gene_expression_grouped %>% 
  left_join(prigroup_normal_distribution[,c("ensembl_id","primary_group_shorthand_label","fpkm_log_normal_dist_group","primary_group_size")],
            by=c("ensembl_id","primary_group_shorthand_label")) %>% 
  left_join(supergroup_normal_distribution[,c("ensembl_id","supergroup_label","fpkm_log_normal_dist_supergroup","supergroup_size")],
            by=c("ensembl_id","supergroup_label"))


write.table(summary_gene_expression_grouped,gene_expression_analysis_path,quote = FALSE,sep = "\t",row.names=FALSE)

