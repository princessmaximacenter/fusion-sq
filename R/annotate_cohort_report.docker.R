## Fusion-sq 
## Annotate cohort report
## Last update: 2022-02-27 
## Implement output from SV pipeline 

# Input: cohort fusions  (from collect_cohort.R)
## 
# Annotate gene fusions with cancer related genes, healthy chimera, chromosome bands/cytobands, 
# Annotate gene fusions with properties of underlying SVs (see utils)
## SV overlaps with population SV databases, repeats, segmental duplications, introns

#Optional - integration with different modalities: 
## Annotate with CNAs and SNVs (see utils)
## Gene expression & z scores (from gene_expression_zscores.R), wilcox test for fusion carrying patients

# Changelog
## 2022-01-12
##  path local overrides and docker version
## 2022-02-17 
## Adjust for merged fusions, focus on the gene level annotation. 
## TODO later implement output from SV pipeline 



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
  library(GenomicRanges)
  library(dplyr)
  library(stringi)
})

source(paste0(script_dir,"functions.general.R")) 

map_template_vars=c('${resources_dir}'=resources_dir,'${input_dir}'=input_dir,'${output_dir}'=output_dir,
                    '${cohort_identifier}'=cohort_identifier,'${cohort_wdir}'=cohort_wdir)

reports_dir =  stri_replace_all_fixed(reports_dir_template,names(map_template_vars), map_template_vars,vectorize=F)
base_dir =  stri_replace_all_fixed(base_dir_template,names(map_template_vars), map_template_vars,vectorize=F)

map_template_vars=c(map_template_vars,'${reports_dir}'=reports_dir,'${base_dir}'=base_dir)

##  merge fusion annotation ends up in reports dir too

#tmp dots
map_template_vars_merged = c('${analysis_type}'=".merged",
                             map_template_vars )

map_template_vars = map_template_vars_merged

## Inputs defined in config file
## Collected for cohort:

cohort_fusion_level_results_path = stri_replace_all_fixed(cohort_fusion_level_results_path_template,names(map_template_vars_merged), map_template_vars_merged,vectorize=F)
#cohort_results_path = stri_replace_all_fixed(cohort_results_path_template,names(map_template_vars_merged), map_template_vars_merged,vectorize=F) 
fusion_overview_path = stri_replace_all_fixed(fusion_overview_path_template,names(map_template_vars_merged), map_template_vars_merged,vectorize=F) 


## Output
cohort_report_path = stri_replace_all_fixed(cohort_report_path_template,names(map_template_vars), map_template_vars,vectorize=F)
fusion_overview_anno_path = stri_replace_all_fixed(fusion_overview_anno_path_template,names(map_template_vars), map_template_vars,vectorize=F) 

## supporting SVs with annotation
cohort_supporting_svs_anno_path = stri_replace_all_fixed(cohort_supporting_svs_anno_path_template,names(map_template_vars), map_template_vars,vectorize=F)


#Gene expression & z scores 
#nb typo
#gene_expression_analysis_path = paste0(reports_dir,gene_expression_zcore_outfile)
gene_expression_zscore_outfile="gene_expression_zscore.tsv"
gene_expression_analysis_path_template = paste0("${reports_dir}",gene_expression_zscore_outfile) 
gene_expression_analysis_path = stri_replace_all_fixed(gene_expression_analysis_path_template,names(map_template_vars), map_template_vars,vectorize=F)


rt_sv_bp_overlap_path_template = "${reports_dir}cohort_supporting_svs.rt_sv_bp_overlap${analysis_type}.tsv"
rt_sv_bp_overlap_path = stri_replace_all_fixed(rt_sv_bp_overlap_path_template,names(map_template_vars), map_template_vars,vectorize=F)


# SNV files => see config

#External resources => see config for files
## TODO replace by the get_cancer_genes() function
if(FALSE) {
grobner_recurrent_path = paste0(resources_dir,grobner_recurrent_file)
oncokb_path = paste0(resources_dir,oncokb_file)
}

cosmic_path = paste0(resources_dir,cosmic_file) #seperate for fusions
grobner_druggable_path = paste0(resources_dir,grobner_druggable_file)

cancer_related_genes = get_cancer_genes(resources_dir)

## TODO: use the templates
chimerseq_path = paste0(resources_dir,chimerseq_file)
mitelman_path = paste0(resources_dir,mitelman_file)

kinases_path = paste0(resources_dir,kinases_file)
transcriptionfactors_path = paste0(resources_dir,transcriptionfactors_file)



#clinically_relevant_fusion_genes_path = paste0(resources_dir,"InclusionList_fusion_genes_Homo_sapiens_GRCh38_BioMart_v1.8.txt")
clinically_relevant_fusion_genes_path = paste0("~/Documents/resources/InclusionList_fusion_genes_Homo_sapiens_GRCh38_BioMart_v1.8.txt")


### Settings ----
## deprecated:
svs_extra_cols=c("FILTER","overlap_distance","overlap_set1_set2","overlap_set1_merged",
                 "start_copy_ratio_log2_mean","end_copy_ratio_log2_mean","copy_ratio_log2_mean","copy_ratio_log2_sd",
                 "repeat_family_start","repeat_family_end","segdup_start","segdup_end",
                 "intron_overlap","intron_tx",
                 "anno_sv_population","anno_sv_db_nstd166","anno_sv_db_nstd186","anno_sv_db_dgv")
## TODO: 
#after refactor:
sv_cna_cols = c("patient_sv_name","cr_l2fc_50_max","maf_50","cr_l2fc_50","start_cr_l2fc_50", "end_cr_l2fc_50", "start_call","end_call")
sv_pop_db_cols = c("flag_sv_db_nstd166","sv_db_nstd166_overlaps","flag_sv_db_nstd186","sv_db_nstd186_overlaps","flag_sv_db_dgv","sv_db_dgv_overlaps","flag_sv_population")
sv_repeat_cols = c("start_repeat_class", "end_repeat_class", "flag_repeat","flag_repeat_both_bp","flag_segdup","flag_segdup_both_bp")

#also used for tx selection any time merging on that level is required
patient_sv_summary_cols=c("fusion_name", "gup_sv_merged", "gdw_sv_merged", "specific_sv", "patient_id")
transcript_cols = c("transcript_type","involved_fragment","involved_fraction","involved_exon","max_exon")

### Resources ----
chromosome_bands_path = stri_replace_all_fixed(chromosome_bands_path_template,names(map_template_vars), map_template_vars,vectorize=F)
chromosome_bands = read.table(chromosome_bands_path,header=T,sep="\t",comment.char = "")
chromosome_bands = chromosome_bands %>% dplyr::rename("seqnames"="X.chrom","start"="chromStart","end"="chromEnd","cytoband"="name") 
chromosome_bands$cytoband = paste0(chromosome_bands$seqnames,chromosome_bands$cytoband)
chromosome_bands = GRanges(chromosome_bands)



## Start annotating fusion overview and then merge with fusion level results to make the report


## START: 
## Read in cohort
cohort = read.table(patient_table_path,sep = "\t", header=T,stringsAsFactors = F)

## TODO I think that is what was meant...
cohort$supergroup=cohort$domain
cohort$supergroup_label=cohort$domain_label

# Exit if exists, otherwise make
if(length(Sys.glob(cohort_report_path))==1){ 
  print(paste0("File already exists: ",cohort_report_path))
  #quit() 
} 
if(length(Sys.glob(cohort_fusion_level_results_path))!=1 ){ 
  print(paste0("Requires input file from collect_cohort: ",cohort_fusion_level_results_path))
  #quit() 
}


cohort_fusion_level_results = read.table(cohort_fusion_level_results_path,sep = "\t", header=T)
cohort_fusion_level_results = cohort_fusion_level_results %>% filter(patient_id %in% cohort$patient_id)

fusion_overview = read.table(fusion_overview_path,sep = "\t", header=T)
fusion_overview = fusion_overview  %>% filter(patient_id %in% cohort$patient_id)


### Fusion overview to cytobands ----
fusion_overview = fusion_overview %>%
  left_join(cohort %>% select(patient_id, tumor_type_label, primary_group, primary_group_label, 
                                        supergroup, supergroup_label,primary_group_shorthand_label), 
            by="patient_id") 



fusions_gup = GRanges(fusion_overview$gup_sf_breakpoint)
mcols(fusions_gup)$patient_fusion = fusion_overview$patient_fusion
mcols(fusions_gup)$side = "gup"

fusions_gdw = GRanges(fusion_overview$gdw_sf_breakpoint)
mcols(fusions_gdw)$patient_fusion = fusion_overview$patient_fusion
mcols(fusions_gdw)$side = "gdw"

fusions_sf_bp = c(fusions_gup,fusions_gdw)

fusions_overlap = findOverlaps(fusions_sf_bp,chromosome_bands)
if(length(fusions_overlap[duplicated(fusions_overlap@from)])==0) {
  fusions_cytoband_mapping = as.data.frame(fusions_sf_bp[fusions_overlap@from]$patient_fusion)
  colnames(fusions_cytoband_mapping)=c("patient_fusion")
  fusions_cytoband_mapping$side = fusions_sf_bp[fusions_overlap@from]$side
  fusions_cytoband_mapping$cytoband = chromosome_bands[fusions_overlap@to]$cytoband
}

fusions_cytoband_mapping =fusions_cytoband_mapping %>% pivot_wider(names_from = side,
                                                                   values_from=cytoband,values_fn=function(x){toString(unique(sort(x)))}) %>%
  dplyr::rename(gup_cytoband=gup,gdw_cytoband=gdw)


fusion_overview = fusion_overview %>% left_join(fusions_cytoband_mapping,by="patient_fusion")


## SV level properties ----
## supporting SVs with  annotation
if(length(Sys.glob(cohort_supporting_svs_anno_path))!=1 ){ 
  print(paste0("No annotated SVs found: ",cohort_supporting_svs_anno_path))
  #quit() 
} else {
  supporting_sv_anno_df = read.table(cohort_supporting_svs_anno_path,sep = "\t",header=T)

  #Do I need cohort results intermediate or can I just merge the SV df and left join to the fusions?
  
  supporting_sv_anno_df$repeat_class = paste0("start:",supporting_sv_anno_df$start_repeat_class,";end:",supporting_sv_anno_df$end_repeat_class)
  supporting_sv_anno_df$repeat_class=str_replace(supporting_sv_anno_df$repeat_class,"start:NA;","")
  supporting_sv_anno_df$repeat_class=str_replace(supporting_sv_anno_df$repeat_class,"end:NA","")
  
  ## First merge over 'merged svs'
  ## those which are not merged have their bp name in that field
  ## next, join on gup_ and gdw sv merged and rename the columns with gup/gdw
  # for the _both_bp flags... => used all 
  #colnames should correspoond to fusion sq pipeline code 
  merged_svs =   supporting_sv_anno_df %>% 
    group_by(sv_merged,patient_id) %>% 
    summarize(
              start_copy_ratio_l2fc = mean(start_cr_l2fc_50,na.rm=T), 
              end_copy_ratio_l2fc = mean(end_cr_l2fc_50,na.rm=T),
              copy_ratio_l2fc_mean = mean(cr_l2fc_50,na.rm=T), 
              maf_mean = mean(maf_50,na.rm=T),
             
              distance_max=max(overlap_distance,na.rm = T), 
              distance_mean=mean(overlap_distance,na.rm=T),
              overlap_min=min(overlap_set1_set2,na.rm = T),
              overlap_mean=mean(overlap_set1_set2,,na.rm = T),
              merged_overlap_min=min(overlap_set1_merged,na.rm = T),
              merged_overlap_mean=mean(overlap_set1_merged,na.rm = T),
              
              repeat_class=toString(sort(unique(repeat_class))),

              flag_repeat=any(flag_repeat), flag_repeat_both_bp=all(flag_repeat_both_bp),
              flag_segdup=any(flag_segdup), flag_segdup_both_bp=all(flag_segdup_both_bp),
              
              flag_sv_population=any(flag_sv_population), flag_sv_db_nstd166=any(flag_sv_db_nstd186),
              flag_sv_db_nstd186=any(flag_sv_db_nstd186), flag_sv_db_dgv=any(flag_sv_db_dgv),
              
              flag_sv_population_all=all(flag_sv_population),
              
              #sv_intron = toString(sort(unique(intron_overlap))),sv_intron_tx = toString(sort(unique(intron_tx))),

              sv_filter=toString(sort(unique(FILTER)))
              )
  
  merged_svs$sv_intron=""
  merged_svs$sv_intron_tx=""
  
    
  svs_anno_cols =  names(merged_svs)
  
  #subset to existing cols
  #svs_anno_cols_missing=svs_anno_cols[!svs_anno_cols %in% names(merged_svs)]
  #svs_anno_cols=svs_anno_cols[svs_anno_cols %in% names(merged_svs)]
  
  svs_anno_cols_wo_merge_cols = svs_anno_cols[!svs_anno_cols %in% c("sv_merged","patient_id","patient_sv_merged")]
  
  fusion_level_svs_cols = c("patient_fusion_sv","patient_sv_id","gup_sv_merged","gdw_sv_merged","patient_id")
  
  fusion_level_svs = cohort_fusion_level_results[,fusion_level_svs_cols] %>% 
    left_join(merged_svs[,svs_anno_cols],by=c("gup_sv_merged"="sv_merged","patient_id")) %>%
    rename_at(svs_anno_cols_wo_merge_cols,function(x){paste0("gup_",x)}) %>% 
    left_join(merged_svs[,svs_anno_cols],by=c("gdw_sv_merged"="sv_merged","patient_id")) %>%
    rename_at(svs_anno_cols_wo_merge_cols,function(x){paste0("gdw_",x)})
  
  ### Summary of supporting SVs per merged

  fusion_level_extra_sv_properties = fusion_level_svs %>% 
    group_by(across(all_of(fusion_level_svs_cols))) %>% 
    summarize(
              pairwise_distance_max=max(c(gup_distance_max,gdw_distance_max),na.rm = T), 
              pairwise_distance_mean=mean(c(gup_distance_mean,gdw_distance_mean),na.rm = T), 
              pairwise_overlap_min=min(c(gup_overlap_min,gdw_overlap_min),na.rm = T),
              pairwise_overlap_mean=mean(c(gup_overlap_mean,gdw_overlap_mean),na.rm = T),
              merged_overlap_min=min(c(gup_merged_overlap_min,gdw_merged_overlap_min),na.rm = T),
              merged_overlap_mean=mean(c(gup_merged_overlap_mean,gdw_merged_overlap_mean),na.rm = T),
              
              repeat_class=toString(sort(unique(paste0("gup:",gup_repeat_class,";gdw:",gdw_repeat_class,";")))),
              
              # segdup=toString(sort(unique(paste0("gup:",gup_segdup,";gdw:",gdw_segdup,";")))),

              flag_repeat=any(c(gup_flag_repeat,gdw_flag_repeat)), 
              flag_segdup=any(c(gup_flag_segdup,gdw_flag_segdup)), 
              flag_repeat_both_bp=all(c(gup_flag_repeat_both_bp,gdw_flag_repeat_both_bp)), 
              flag_segdup_both_bp=all(c(gup_flag_segdup_both_bp,gdw_flag_segdup_both_bp)), 
              
              flag_sv_population=any(c(gup_flag_sv_population,gdw_flag_sv_population)), 
              flag_sv_db_nstd186=any(c(gup_flag_sv_db_nstd186,gdw_flag_sv_db_nstd186)), 
              flag_sv_db_nstd166=any(c(gup_flag_sv_db_nstd166,gdw_flag_sv_db_nstd166)), 
              flag_sv_db_dgv=any(c(gup_flag_sv_db_dgv,gdw_flag_sv_db_dgv)), 
              
              flag_sv_population_all=all(c(gup_flag_sv_population_all,gdw_flag_sv_population_all))
    )

   
   #prev
  if(FALSE) {

## summarized by SV merged
fusion_level_extra_sv_properties = cohort_results %>% 
  group_by(across(all_of(patient_sv_summary_cols))) %>% 
  summarize(gup_start_copy_ratio_l2fc = mean(gup_start_copy_ratio_log2_mean,na.rm=T), gup_end_copy_ratio_l2fc = mean(gup_end_copy_ratio_log2_mean,na.rm=T),
            gdw_start_copy_ratio_l2fc = mean(gdw_start_copy_ratio_log2_mean,na.rm=T), gdw_end_copy_ratio_l2fc = mean(gdw_end_copy_ratio_log2_mean,na.rm=T),
            gup_copy_ratio_l2fc_mean = mean(gup_copy_ratio_log2_mean,na.rm=T), gup_copy_ratio_l2fc_sd = mean(gup_copy_ratio_log2_sd,na.rm=T),
            gdw_copy_ratio_l2fc_mean = mean(gdw_copy_ratio_log2_mean,na.rm=T), gdw_copy_ratio_l2fc_sd = mean(gdw_copy_ratio_log2_sd,na.rm=T),
            pairwise_distance_max=max(c(gup_overlap_distance,gdw_overlap_distance),na.rm = T), 
            pairwise_distance_mean=mean(c(gup_overlap_distance,gdw_overlap_distance),na.rm=T),
            pairwise_overlap_min=min(c(gup_overlap_set1_set2,gdw_overlap_set1_set2),na.rm = T),
            pairwise_overlap_mean=mean(c(gup_overlap_set1_set2,gdw_overlap_set1_set2),na.rm = T),
            merged_overlap_min=min(c(gup_overlap_set1_merged,gdw_overlap_set1_merged),na.rm = T),
            merged_overlap_mean=mean(c(gup_overlap_set1_merged,gdw_overlap_set1_merged),na.rm = T),
            repeat_family=toString(sort(unique(paste0("gup:",gup_repeat_family,";gdw:",gdw_repeat_family,";")))),
            segdup=toString(sort(unique(paste0("gup:",gup_segdup,";gdw:",gdw_segdup,";")))),
            gup_sv_intron = toString(sort(unique(gup_intron_overlap))),gup_sv_intron_tx = toString(sort(unique(gup_intron_tx))),
            gdw_sv_intron = toString(sort(unique(gdw_intron_overlap))),gdw_sv_intron_tx = toString(sort(unique(gdw_intron_tx))),
            anno_sv_population=any(gup_anno_sv_population,gdw_anno_sv_population), anno_sv_db_nstd166=any(gup_anno_sv_db_nstd166,gdw_anno_sv_db_nstd186),
            anno_sv_db_nstd186=any(gup_anno_sv_db_nstd186,gdw_anno_sv_db_nstd186) ,anno_sv_db_dgv=any(gup_anno_sv_db_dgv,gdw_anno_sv_db_dgv),
            sv_gup_filter=toString(sort(unique(gup_FILTER))),sv_gdw_filter=toString(sort(unique(gdw_FILTER))))
           
  }
  
 fusion_level_extra_sv_properties$repeat_class = str_replace(fusion_level_extra_sv_properties$repeat_class,"gup:;","")
 fusion_level_extra_sv_properties$repeat_class = str_replace(fusion_level_extra_sv_properties$repeat_class,"gdw:;","")
# fusion_level_extra_sv_properties$segdup = str_replace(fusion_level_extra_sv_properties$segdup,"gup:;","")
# fusion_level_extra_sv_properties$segdup = str_replace(fusion_level_extra_sv_properties$segdup,"gdw:;","")
  
  #from sv pipeline colnames to fusion sq pipeline
 fusion_level_extra_sv_properties = fusion_level_extra_sv_properties %>% mutate(
    anno_sv_population = flag_sv_population_all,
    repeat_family = repeat_class
    )
    
  

cohort_fusion_level_results = cohort_fusion_level_results %>% 
  left_join(fusion_level_svs,by=fusion_level_svs_cols) %>%
  left_join(fusion_level_extra_sv_properties,by=fusion_level_svs_cols)



}


### Replication timing domains ----

if(length(Sys.glob(rt_sv_bp_overlap_path))!=1 ){ 
  print(paste0("No annotation for replication timing domains found: ",rt_sv_bp_overlap_path))
  #quit() 
} else {
  rt_sv_bp_overlaps = read.table(rt_sv_bp_overlap_path,sep="\t",header=T)
  rt_sv_bp_wide = rt_sv_bp_overlaps %>% 
    dplyr::select(patient_id,sv_merged,patient_sv_name,sv_breakpoint_orientation,rt_anno) %>% 
    unique() %>% pivot_wider(names_from = sv_breakpoint_orientation,values_from=rt_anno) %>%
    dplyr::rename(rt_start=start,rt_end=end)
  
  ## Same as above: First merge over 'merged svs'
  ## next, join on gup_ and gdw sv merged and rename the columns with gup/gdw

  merged_rt_sv_bp = rt_sv_bp_wide %>% 
    group_by(sv_merged,patient_id) %>% 
    summarize(
      rt_sv_cnt=n(),
      rt_start_switch = sum(rt_start=="S"),
      rt_start_early = sum(rt_start=="CE"),
      rt_start_late = sum(rt_start=="CL"),
      rt_end_switch = sum(rt_end=="S"),
      rt_end_early = sum(rt_end=="CE"),
      rt_end_late = sum(rt_end=="CL"),
      rt_switch = sum(c(rt_start_switch,rt_end_switch)),
      rt_early = sum(c(rt_start_early,rt_end_early)),
      rt_late = sum(c(rt_start_late,rt_end_late)),
      rt_anno_start = toString(unique(sort(rt_start))),
      rt_anno_end = toString(unique(sort(rt_end)))
    ) 
  
  rt_anno_cols= names(merged_rt_sv_bp)[grepl("rt_",names(merged_rt_sv_bp))]  #c("rt_anno_start","rt_anno_end")
  fusion_level_svs_cols = c("patient_fusion_sv","patient_sv_id","gup_sv_merged","gdw_sv_merged","patient_id")
  
  ## merged_rt_sv_bp %>% filter(sv_cnt==rt_early) means it is early but can also be sv_cnt==rt_switch at the same time if start/end bp are different
  ## would you call then that sv early or late? Or count as 1 early bp and 1 late bp 
  
  gup_gdw_rt_anno = cohort_fusion_level_results[,c(fusion_level_svs_cols,"svtype")] %>% unique() %>% 
    left_join(merged_rt_sv_bp,by=c("gup_sv_merged"="sv_merged","patient_id")) %>%
    rename_at(rt_anno_cols,function(x){paste0("gup_",x)}) %>% 
    left_join(merged_rt_sv_bp,by=c("gdw_sv_merged"="sv_merged","patient_id")) %>%
    rename_at(rt_anno_cols,function(x){paste0("gdw_",x)})
  
  ### Summary of supporting SVs per merged
  # for intra-chr the gup/gdw will be the same and you want gup_start and end_gdw. for CTX they are different but can still do that?
 
   fusion_level_rt_anno = gup_gdw_rt_anno %>% 
    group_by(across(all_of(fusion_level_svs_cols))) %>% 
    summarize(
      #rt_anno_start=toString(sort(unique(paste0("gup:",gup_rt_anno_start,";gdw:",gdw_rt_anno_start,";")))),
      #rt_anno_end=toString(sort(unique(paste0("gup:",gup_rt_anno_end,";gdw:",gdw_rt_anno_end,";")))),
      
      gup_rt_anno=toString(sort(unique(paste0("start:",gup_rt_anno_start,";end:",gup_rt_anno_end,";")))),
      gdw_rt_anno=toString(sort(unique(paste0("start:",gdw_rt_anno_start,";end:",gdw_rt_anno_end,";")))),
      rt_switch = sum(c(gup_rt_start_switch,gup_rt_end_switch,gdw_rt_start_switch,gdw_rt_end_switch)),
      rt_early = sum(c(gup_rt_start_early,gup_rt_end_early,gdw_rt_start_early,gdw_rt_end_early)),
      rt_late = sum(c(gup_rt_start_late,gup_rt_end_late,gdw_rt_start_late,gdw_rt_end_late))
    )
  
  
   gup_gdw_rt_anno = gup_gdw_rt_anno %>%
    left_join(fusion_level_rt_anno,by=fusion_level_svs_cols)
  
  #no correction for sv type needed, does its job to normalize each fusion to 2 breakpoints
   ## Note: complex/composite variants can cause fractions like 1.5 and 0.5 
  gup_gdw_rt_anno = gup_gdw_rt_anno %>% 
    mutate( 
            rt_switch_bp = rt_switch/( gup_rt_sv_cnt+gdw_rt_sv_cnt),
            rt_early_bp = rt_early/(gup_rt_sv_cnt+gdw_rt_sv_cnt),
            rt_late_bp = rt_late/(gup_rt_sv_cnt+gdw_rt_sv_cnt)
            )
 
  #interesting cols: 
  fusion_rt_anno_cols = c("gup_rt_anno","gdw_rt_anno","gup_rt_anno_start","gdw_rt_anno_end","rt_switch_bp","rt_early_bp","rt_late_bp")

  cohort_fusion_level_results = cohort_fusion_level_results %>% 
    left_join(gup_gdw_rt_anno[,c(fusion_level_svs_cols,fusion_rt_anno_cols)],by=fusion_level_svs_cols)
  }




## Transcript properties ----
#TODO: unskip
if(FALSE) {
fusion_tx_extra_cols=data.frame()
for(id in cohort$patient_identifier) {
  patient = filter(cohort,patient_identifier==id)
  
  map_template_vars_patient = c('${analysis_type}'=analysis_type_filename_patient,
                                map_template_vars,
                                '${patient_identifier}'=as.character(patient$patient_identifier) )
  
  fusion_tx_selection_path = stri_replace_all_fixed(fusion_tx_selection_path_template,names(map_template_vars_patient), map_template_vars_patient,vectorize=F)
  
  
  if( length(Sys.glob(fusion_tx_selection_path))==1 & file.size(fusion_tx_selection_path)>5 ) { 
    fusion_tx_selection_summary = read.table(fusion_tx_selection_path,header=T,sep = "\t")
    fusion_tx_selection_summary$transcript_id = remove_version_from_id(fusion_tx_selection_summary$transcript_id)
   
    
    gup_fusion_transcripts = fusion_tx_selection_summary[fusion_tx_selection_summary$upstream==T,c("fusion_id",transcript_cols)] %>%
      rename_at(all_of(transcript_cols), function(x){paste0("tx_gup_",x)})
    gdw_fusion_transcripts = fusion_tx_selection_summary[fusion_tx_selection_summary$upstream==T,c("fusion_id",transcript_cols)] %>%
      rename_at(all_of(transcript_cols), function(x){paste0("tx_gdw_",x)})
    fusion_transcripts = gup_fusion_transcripts %>% left_join(gdw_fusion_transcripts,by="fusion_id")
    
    fusion_transcripts$patient_id=patient$patient_id
    fusion_tx_properties = cohort_results[,c("fusion_id",patient_sv_summary_cols)] %>% unique() %>% left_join(fusion_transcripts,by=c("patient_id","fusion_id")) %>%
      group_by(across(all_of(patient_sv_summary_cols))) %>% summarize(across(paste0("tx_gup_",transcript_cols),toString),
                                                                    across(paste0("tx_gdw_",transcript_cols),toString),.groups="keep")
  
    fusion_tx_extra_cols  =rbind(fusion_tx_extra_cols,fusion_tx_properties)
  } else {
    print(paste0("No transcript information found for patient: ",patient$patient_identifier," ",fusion_tx_selection_path))
    
  }
} 

cohort_fusion_level_results = cohort_fusion_level_results %>% 
  left_join(fusion_tx_extra_cols,by=patient_sv_summary_cols)

}

### External annotation of chimera and cancer genes ----
## Chimera annotation ----


## ChimerSeq
chimerseq = read.table(chimerseq_path,sep = "\t", header=T)
chimerseq = chimerseq %>% separate(col="Fusion_pair",into=c("H_gene","T_gene"),sep = "--",remove="FALSE",extra="merge")
chimerseq_summary = chimerseq %>% group_by(Fusion_pair) %>% summarize(anno_chimerseq_cancer_types = toString(unique(Cancertype))) %>%
  dplyr::rename(fusion_name = Fusion_pair)

fusion_overview = fusion_overview %>% 
  mutate(anno_gup_chimerseq = gup_gene_id %in% chimerseq$H_gene) %>%
  mutate(anno_gdw_chimerseq = gdw_gene_id %in% chimerseq$T_gene) %>%
  mutate(anno_chimerseq = fusion_name %in% chimerseq$Fusion_pair)

fusion_overview = fusion_overview %>% left_join(chimerseq_summary) 


## Mitelman database
mitelman = read.table(mitelman_path,sep = "\t", header=T)
mitelman_fusions = mitelman %>% filter(grepl("/",Gene,fixed=T))
mitelman_fusions = mitelman_fusions %>% separate(col="Gene",into=c("H_gene","T_gene"),sep = "/",remove="FALSE",extra="merge")
mitelman_fusions$fusion_name = paste0(mitelman_fusions$H_gene,"--",mitelman_fusions$T_gene)
mitelman_fusions_summary = mitelman_fusions %>% group_by(fusion_name) %>% summarise(anno_mitelman_ref = toString(sort(unique(RefNo))))

fusion_overview = fusion_overview %>% 
  mutate(anno_gup_mitelman = gup_gene_id %in% mitelman_fusions$H_gene) %>%
  mutate(anno_gdw_mitelman = gdw_gene_id %in% mitelman_fusions$T_gene) %>%
  mutate(anno_mitelman = fusion_name %in% mitelman_fusions$fusion_name)

fusion_overview = fusion_overview %>% left_join(mitelman_fusions_summary) 

## COSMIC
cosmic = read.table(cosmic_path,header=T,sep="\t")
cosmic_fusions = cosmic %>% filter(grepl("fusion",Role.in.Cancer)) %>% select(Gene.Symbol,Synonyms,Tumour.Types.Somatic.) 

fusion_overview = fusion_overview %>% 
  mutate(anno_gup_cosmic_fusion = (gup_gene_id %in% cosmic_fusions$Gene.Symbol | gup_gene_id %in% cosmic_fusions$Synonyms | gup_ensembl_id %in% cosmic_fusions$Synonyms)) %>% 
  mutate(anno_gdw_cosmic_fusion = (gdw_gene_id %in% cosmic_fusions$Gene.Symbol | gdw_gene_id %in% cosmic_fusions$Synonyms | gdw_ensembl_id %in% cosmic_fusions$Synonyms)) %>% 
  mutate(anno_has_cosmic_fusion = (anno_gup_cosmic_fusion | anno_gdw_cosmic_fusion)) 

## Cancer related gene annotation ----
#for onco/TSG/fusions together with oncokb and grobner

fusion_overview = fusion_overview %>% 
  mutate(anno_gup_cosmic = (gup_gene_id %in% filter(cancer_related_genes,source=="cosmic")$gene_id
 | gup_ensembl_id %in% filter(cancer_related_genes,source=="cosmic")$gene_id)) %>% 
  mutate(anno_gdw_cosmic = (gdw_gene_id %in% filter(cancer_related_genes,source=="cosmic")$gene_id
 | gdw_ensembl_id %in% filter(cancer_related_genes,source=="cosmic")$gene_id)) %>% 
  mutate(anno_cosmic = (anno_gup_cosmic | anno_gdw_cosmic))


## Grobner pediatric cancer
fusion_overview = fusion_overview %>% 
  mutate(anno_gup_grobner_recurrent = (gup_gene_id %in% filter(cancer_related_genes,source=="grobner")$gene_id)) %>% 
  mutate(anno_gdw_grobner_recurrent = (gdw_gene_id %in% filter(cancer_related_genes,source=="grobner")$gene_id)) %>% 
  mutate(anno_grobner_recurrent = (anno_gup_grobner_recurrent | anno_gdw_grobner_recurrent)) 

grobner_druggable = read.table(grobner_druggable_path,header=T,sep="\t")
#gene snv_indel_analysis cnv_analysis
##  for now filter on cnv_analysis == "yes"
grobner_druggable = filter(grobner_druggable,cnv_analysis=="yes")$gene

fusion_overview = fusion_overview %>% 
  mutate(anno_gup_grobner_druggable = (gup_gene_id %in% grobner_druggable)) %>% 
  mutate(anno_gdw_grobner_druggable = (gdw_gene_id %in% grobner_druggable)) %>% 
  mutate(anno_grobner_druggable = (anno_gup_grobner_druggable | anno_gdw_grobner_druggable)) 

## OncoKB
fusion_overview = fusion_overview %>% 
  mutate(anno_gup_oncokb = (gup_gene_id %in% filter(cancer_related_genes,source=="oncokb")$gene_id)) %>% 
  mutate(anno_gdw_oncokb = (gdw_gene_id %in% filter(cancer_related_genes,source=="oncokb")$gene_id)) %>% 
  mutate(anno_oncokb = (anno_gup_oncokb | anno_gdw_oncokb))



## onco tsg contain both cosmic and grobner
oncogenes=filter(cancer_related_genes,oncogene==T)$gene_id
tsg=filter(cancer_related_genes,tsg==T)$gene_id

fusion_overview = fusion_overview %>%
  mutate(anno_gup_oncogene = (gup_gene_id %in% oncogenes | gup_ensembl_id %in% oncogenes)) %>%
  mutate(anno_gdw_oncogene = (gdw_gene_id %in% oncogenes | gdw_ensembl_id %in% oncogenes)) %>% 
  mutate(anno_has_oncogene = (anno_gup_oncogene | anno_gdw_oncogene)) 

fusion_overview = fusion_overview %>% 
  mutate(anno_gup_tsg = (gup_gene_id %in% tsg | gup_ensembl_id %in% tsg)) %>% 
  mutate(anno_gdw_tsg = (gdw_gene_id %in% tsg | gdw_ensembl_id %in% tsg)) %>% 
  mutate(anno_has_tsg = (anno_gup_tsg | anno_gdw_tsg)) 


## end of cancer gene annotation


## Kinases
kinases_df = read.table(kinases_path,header=T,sep="\t")
kinases = kinases_df$gene_name

fusion_overview = fusion_overview %>% 
  mutate(anno_gup_kinase = (gup_gene_id %in% kinases | gup_ensembl_id %in% kinases)) %>%
  mutate(anno_gdw_kinase = (gdw_gene_id %in% kinases | gdw_ensembl_id %in% kinases)) %>%
  mutate(anno_has_kinase = (anno_gup_kinase | anno_gdw_kinase))

## TFs Lambert et al 2018
transcriptionfactors_df = read.table(transcriptionfactors_path,header=T,sep="\t")
transcriptionfactors = transcriptionfactors_df[transcriptionfactors_df$Is.TF.=="Yes",]
transcriptionfactors = c(transcriptionfactors$Ensembl.ID,transcriptionfactors$HGNC.symbol)

fusion_overview = fusion_overview %>% 
  mutate(anno_gup_tf = (gup_gene_id %in% transcriptionfactors | gup_ensembl_id %in% transcriptionfactors)) %>%
  mutate(anno_gdw_tf = (gdw_gene_id %in% transcriptionfactors | gdw_ensembl_id %in% transcriptionfactors)) %>%
  mutate(anno_has_tf = (anno_gup_tf | anno_gdw_tf))

## Summarize cancer chimera annotation in single variables ----
# anno_cancer_chimera: one variable for if exact gup-gdw in chimerseq or mitelman, not cosmic fusion or star fusion database because those are not exact pairs
# anno_partner_cancer_chimera: either gup or gdw in mitelman chimerseq or cosmic fusion
# anno_in_any_cancer: anno partner cancer chimera or in cosmic(onco/tsg) or in star fusion cancer database
# anno_cancer_gene_db: any of the partner genes in cosmic onco kb or grobner

fusion_overview = fusion_overview %>% mutate(anno_cancer_chimera = (anno_mitelman | anno_chimerseq))
fusion_overview = fusion_overview %>% mutate(anno_partner_cancer_chimera = 
                                               (anno_gup_mitelman | anno_gdw_mitelman | 
                                                  anno_gup_chimerseq | anno_gdw_chimerseq |
                                                  anno_has_cosmic_fusion))


fusion_overview = fusion_overview %>% mutate(anno_cancer_gene_db = 
                                               anno_cosmic | anno_oncokb | anno_grobner_recurrent)

fusion_overview = fusion_overview %>% mutate(anno_in_any_cancer = 
                                               (anno_partner_cancer_chimera | anno_cancer_gene_db | anno_grobner_druggable))

## TODO LATER: properties from STAR Fusion
#  fusion_overview = fusion_overview %>%
#    mutate(anno_gup_pfam_kinase = grepl("kinase",PFAM_LEFT,fixed=T)) %>%
#    mutate(anno_gdw_pfam_kinase = grepl("kinase",PFAM_RIGHT,fixed=T))
  
#  fusion_overview = fusion_overview %>% mutate(anno_neighbors_overlap_starfusion = grepl("NEIGHBORS_OVERLAP",annots_starfusion)) 
  
  #fusion_overview = fusion_overview %>% mutate(anno_in_any_cancer = 
  #                                               (anno_in_any_cancer | anno_starfusion_cancer_db ))
    

## Annotate clinically_relevant fusions if annotation is present

if("fusion" %in% names(cohort)) {
  #fusion == true 
  fusion_overview = fusion_overview %>% mutate(patient_clinrel_fusion = patient_id %in% filter(cohort,fusion_status)$patient_id)
  #check fusion_overview %>% filter(patient_clinrel_fusion) %>% uq_patients() == filter(cohort,fusion_status) %>% uq_patients()
  
  fusion_overview$clinically_validated=F
  fusion_overview$clinically_validated_reciprocal=F
  
  #harmonize gene identifiers / names
  cohort=cohort %>% mutate(
    left_gene = str_replace(str_replace(left_gene,fixed("IGH@"),"IGH"),fixed("IGH-@-ext"),"IGH"),
    right_gene = str_replace(str_replace(right_gene,fixed("IGH@"),"IGH"),fixed("IGH-@-ext"),"IGH"))
    
  for(pid in filter(cohort,fusion_status)$patient_identifier) {
    patient = filter(cohort,patient_identifier==pid)
    
    #check first if patient is in fusion overview!
    if(nrow(fusion_overview[fusion_overview$patient_id==patient$patient_id,])==0) { next() }
    
    fusion_overview = fusion_overview %>% dplyr::mutate(clinically_validated = ifelse(
      (patient_id==patient$patient_id & gup_gene_id == patient$left_gene & gdw_gene_id == patient$right_gene),TRUE,clinically_validated))
    
    fusion_overview = fusion_overview %>% dplyr::mutate(clinically_validated_reciprocal = ifelse(
      (patient_id==patient$patient_id & gdw_gene_id == patient$left_gene & gup_gene_id == patient$right_gene),TRUE,clinically_validated_reciprocal))
                                                          
    
  }
}


## Inclusion list / 

if(length(Sys.glob(clinically_relevant_fusion_genes_path))==1){
  clinrel_fusion_genes = read.table(clinically_relevant_fusion_genes_path,header=T,sep="\t")
  clinrel_fusion_genes = c(clinrel_fusion_genes$Gene.stable.ID,clinrel_fusion_genes$Gene.name) %>% unique()
  
  clinrel_fusion_genes = str_replace(str_replace(clinrel_fusion_genes,fixed("IGH@"),"IGH"),fixed("IGH-@-ext"),"IGH")
  
  
  fusion_overview = fusion_overview %>% 
    mutate(anno_gup_clinrel_fusion_gene = (gup_gene_id %in% clinrel_fusion_genes | gup_ensembl_id %in% clinrel_fusion_genes)) %>%
    mutate(anno_gdw_clinrel_fusion_gene = (gdw_gene_id %in% clinrel_fusion_genes | gdw_ensembl_id %in% clinrel_fusion_genes)) %>%
    mutate(anno_has_clinrel_fusion_gene = (anno_gup_tf | anno_gdw_tf))
}



### Optional: expression, SNVs and CNA if available

## Gene expression analysis ----

if(length(Sys.glob(gene_expression_analysis_path))==1){
  gene_expression_analysis = read.table(gene_expression_analysis_path,sep="\t",header = T)
  
  gene_expression_cols = c("patient_id","ensembl_id","fpkm","fpkm_zscore","fpkm_log_normal_dist",
                           "fpkm_zscore_group","fpkm_log_normal_dist_group",
                           "fpkm_zscore_supergroup","fpkm_log_normal_dist_supergroup")
  
  ## Annotate with expression and z scores 
  fusion_overview = fusion_overview %>% 
    left_join(gene_expression_analysis[,gene_expression_cols],by=c("gup_ensembl_id"="ensembl_id","patient_id")) %>% 
    dplyr::rename(gup_fpkm = fpkm, gup_fpkm_zscore=fpkm_zscore,gup_fpkm_log_normal_dist=fpkm_log_normal_dist,
           gup_fpkm_zscore_group=fpkm_zscore_group, gup_fpkm_log_normal_dist_group=fpkm_log_normal_dist_group,
           gup_fpkm_zscore_supergroup=fpkm_zscore_supergroup, gup_fpkm_log_normal_dist_supergroup=fpkm_log_normal_dist_supergroup) %>% 
    left_join(gene_expression_analysis[,gene_expression_cols],by=c("gdw_ensembl_id"="ensembl_id","patient_id")) %>% 
    dplyr::rename(gdw_fpkm=fpkm,gdw_fpkm_zscore=fpkm_zscore,gdw_fpkm_log_normal_dist=fpkm_log_normal_dist,
           gdw_fpkm_zscore_group=fpkm_zscore_group, gdw_fpkm_log_normal_dist_group=fpkm_log_normal_dist_group,
           gdw_fpkm_zscore_supergroup=fpkm_zscore_supergroup, gdw_fpkm_log_normal_dist_supergroup=fpkm_log_normal_dist_supergroup) 
  
  
} 

## Add SNVs for every patient if exists ----

##TODO just load cohort snv file?
if(FALSE) {
fusion_overview$gup_snv_high_impact=NA
fusion_overview$gdw_snv_high_impact=NA
fusion_overview$gup_snv_moderate_impact=NA
fusion_overview$gdw_snv_moderate_impact=NA

cohort_snvs=data.frame()
for(id in cohort$patient_identifier) {
  patient = filter(cohort,patient_identifier==id)
  snv_gene_summary_path = paste0(snv_summary_dir,snv_gene_summary_outfile,patient$patient_identifier,".tsv")
  
  if( length(Sys.glob(snv_gene_summary_path)) == 1) {
    snvs = read.table(snv_gene_summary_path,header=T, sep="\t")
    snvs$patient_id = patient$patient_id
    
    high_impact_snvs = filter(snvs,grepl("HIGH",vep_impact))
    moderate_impact_snvs = filter(snvs,grepl("MODERATE",vep_impact))
    
    # select 1 patient only
    fusion_overview[fusion_overview$patient_id==patient$patient_id & fusion_overview$gup_ensembl_id %in% high_impact_snvs$Gene,c("gup_snv_high_impact")]=TRUE
    fusion_overview[fusion_overview$patient_id==patient$patient_id & fusion_overview$gdw_ensembl_id %in% high_impact_snvs$Gene,c("gdw_snv_high_impact")]=TRUE
    fusion_overview[fusion_overview$patient_id==patient$patient_id & fusion_overview$gup_ensembl_id %in% moderate_impact_snvs$Gene,c("gup_snv_moderate_impact")]=TRUE
    fusion_overview[fusion_overview$patient_id==patient$patient_id & fusion_overview$gdw_ensembl_id %in% moderate_impact_snvs$Gene,c("gdw_snv_moderate_impact")]=TRUE
  } else {
    fusion_overview[fusion_overview$patient_id==patient$patient_id ,c("gup_snv_high_impact")]=FALSE
    fusion_overview[fusion_overview$patient_id==patient$patient_id ,c("gdw_snv_high_impact")]=FALSE
    fusion_overview[fusion_overview$patient_id==patient$patient_id ,c("gup_snv_moderate_impact")]=FALSE
    fusion_overview[fusion_overview$patient_id==patient$patient_id ,c("gdw_snv_moderate_impact")]=FALSE
    
  }
  
}

}
### End of additional modalities

## healthy chimera from SF or FC ----
#any by default because fusion can be called by single tool. 

healthy_chimera_starfusion_lst = fusion_overview[(fusion_overview$anno_healthy_chimera_starfusion==T),c("fusion_name")]
healthy_chimera_fusioncatcher_lst = fusion_overview[(fusion_overview$anno_healthy_chimera_fusioncatcher==T),c("fusion_name")]

fusion_overview[is.na(fusion_overview$anno_healthy_chimera_starfusion),c("anno_healthy_chimera_starfusion")]=F
fusion_overview[is.na(fusion_overview$anno_healthy_chimera_fusioncatcher),c("anno_healthy_chimera_fusioncatcher")]=F

fusion_overview = fusion_overview %>%
  mutate( anno_healthy_chimera_starfusion = fusion_name %in% healthy_chimera_starfusion_lst,
          anno_healthy_chimera_fusioncatcher= fusion_name %in% healthy_chimera_fusioncatcher_lst,
          anno_healthy_chimera_all = (anno_healthy_chimera_fusioncatcher==T&anno_healthy_chimera_starfusion==T),
          anno_healthy_chimera = (anno_healthy_chimera_fusioncatcher==T|anno_healthy_chimera_starfusion==T))
        

#check fusion_overview %>% filter(is.na(anno_healthy_chimera))


## Merge cohort fusion level results with the annotation from fusion overview ----

# Make sure to not include the columns uq to either tool to prevent explosion of rows
prediction_ignore_shared_cols=c("annots","identifier","fusion_identifier","anno_healthy_chimera","fusion_bp_distance")

ignore_anno_cols = c(paste0(prediction_ignore_shared_cols,"_starfusion"),paste0(prediction_ignore_shared_cols,"_fusioncatcher"),
                     "anno_starfusion_cancer_db","anno_fusioncatcher_cancer")

fusion_overview_anno_only = fusion_overview  %>% select(-ignore_anno_cols) %>%
  select(patient_fusion, tumor_type_label, primary_group, primary_group_label, supergroup, supergroup_label, primary_group_shorthand_label,
                                                        contains("anno"),contains("fpkm"),contains("fpkm_log"),contains("snv"),contains("cytoband"),
                                                        contains("clinrel"),contains("clinically_validated")) %>% unique() 

#check 
fusion_overview_anno_only$patient_fusion %>% unique() %>% length() == fusion_overview_anno_only %>% nrow()

cohort_report = cohort_fusion_level_results %>% 
  left_join(fusion_overview_anno_only, by="patient_fusion")

cohort_report=cohort_report %>% mutate(tumor_normal_diff_mean = (tumor_af_mean-normal_af_mean))

## Summarise over rna level predictions

  #if annotated as pfam kinase any time
# cohort_report = cohort_report %>% 
#    mutate(anno_gup_pfam_kinase = patient_fusion %in% filter(fusion_overview,gup_pfam_kinase)$patient_fusion) %>% 
#    mutate(anno_gdw_pfam_kinase = patient_fusion %in% filter(fusion_overview,gdw_pfam_kinase)$patient_fusion)

cohort_report = cohort_report %>%
  left_join( cohort %>% select(patient_id,tumor_id,normal_id,rna_id), by="patient_id")


#summary labels 
cohort_report=cohort_report %>% annotate_labels_cancer_common() %>% annotate_labels_variant_type() %>% annotate_labels_genes()

##refactor labels
cohort_report=cohort_report %>% mutate(
  high_confidence=precise_confident,
  anno_has_onco_or_tsg = (anno_has_oncogene | anno_has_tsg),
  anno_clinically_relevant = (clinically_validated|clinically_validated_reciprocal)
  )

if("gup_fpkm_zscore_group" %in% names(cohort_report)){
    
  ## If expression exists..
  #this is done here because needs to be on cohort level and not general gene expression but specific to fusions
  expression_wilcoxon=data.frame()
  for(candidate_fusion in unique(cohort_report$patient_fusion)) {
    entry=c()
    candidate = cohort_report %>% filter(patient_fusion==candidate_fusion) %>% 
      select(fusion_name,patient_id,patient_fusion,gup_ensembl_id,gdw_ensembl_id,gup_gene_id,gdw_gene_id,supergroup_label,contains("fpkm")) %>% unique()
    
    entry$patient_fusion = candidate$patient_fusion
    
    if(!is.na(candidate$gup_ensembl_id)&nrow(gene_expression_analysis[gene_expression_analysis$ensembl_id==candidate$gup_ensembl_id,])>0){
      gene_exp_candiate = gene_expression_analysis[gene_expression_analysis$ensembl_id==candidate$gup_ensembl_id,]
      w=wilcox.test(gene_exp_candiate[gene_exp_candiate$patient_id!=candidate$patient_id,c("fpkm_log")],gene_exp_candiate[gene_exp_candiate$patient_id==candidate$patient_id,c("fpkm_log")],alternative = "two.sided")
      entry$gup_fpkm_pval=w$p.value
      
      
      gene_exp_candiate = gene_exp_candiate[gene_exp_candiate$supergroup_label==candidate$supergroup_label,]
      w=wilcox.test(gene_exp_candiate[gene_exp_candiate$patient_id!=candidate$patient_id,c("fpkm_log")],gene_exp_candiate[gene_exp_candiate$patient_id==candidate$patient_id,c("fpkm_log")],alternative = "two.sided")
      entry$gup_fpkm_supergroup_pval=w$p.value  
    } else {
      entry$gup_fpkm_pval=NA
      entry$gup_fpkm_supergroup_pval=NA
    }
    
    if(!is.na(candidate$gdw_ensembl_id)&nrow(gene_expression_analysis[gene_expression_analysis$ensembl_id==candidate$gdw_ensembl_id,])>0){
      gene_exp_candiate = gene_expression_analysis[gene_expression_analysis$ensembl_id==candidate$gdw_ensembl_id,]
      w=wilcox.test(gene_exp_candiate[gene_exp_candiate$patient_id!=candidate$patient_id,c("fpkm_log")],gene_exp_candiate[gene_exp_candiate$patient_id==candidate$patient_id,c("fpkm_log")],alternative = "two.sided")
      entry$gdw_fpkm_pval=w$p.value
      gene_exp_candiate = gene_exp_candiate[gene_exp_candiate$supergroup_label==candidate$supergroup_label,]
      w=wilcox.test(gene_exp_candiate[gene_exp_candiate$patient_id!=candidate$patient_id,c("fpkm_log")],gene_exp_candiate[gene_exp_candiate$patient_id==candidate$patient_id,c("fpkm_log")],alternative = "two.sided")
      entry$gdw_fpkm_supergroup_pval=w$p.value
    } else {
      entry$gdw_fpkm_pval=NA
      entry$gdw_fpkm_supergroup_pval=NA
    }
    
    expression_wilcoxon = rbind(expression_wilcoxon,entry)
  }
  
  expression_wilcoxon_outfile="expression_wilcoxon.tsv"
  expression_wilcoxon_path_template = paste0("${reports_dir}",expression_wilcoxon_outfile) 
  expression_wilcoxon_path = stri_replace_all_fixed(expression_wilcoxon_path_template,names(map_template_vars), map_template_vars_merged,vectorize=F)
  
  write.table(expression_wilcoxon,expression_wilcoxon_path,,quote = FALSE,sep = "\t",row.names=FALSE)
  
  cohort_report = cohort_report %>% left_join(expression_wilcoxon,by="patient_fusion")
  
}


write.table(cohort_report,cohort_report_path,quote = FALSE,sep = "\t",row.names=FALSE)
write.table(fusion_overview,fusion_overview_anno_path,quote = FALSE,sep = "\t",row.names=FALSE)


## TODO: this code is also in make supplementary for reporting issues, but I want to annotate cohort report also with the RNA confidence variables

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


cohort_report = cohort_report %>% 
  mutate(rna_consensus = (patient_fusion_sv %in% c(helper_both$patient_fusion_sv,
                                                   reciprocal_fusion_same_sv$patient_fusion_sv,
                                                   reciprocal_fusion_same_sv$helper_reciprocal_patient_fusion_sv)))

cohort_report = cohort_report %>% mutate(not_precise_confident = !patient_fusion %in% filter(cohort_report,precise_confident)$patient_fusion)


cohort_report = cohort_report %>% mutate(high_confidence_rna_dna = rna_consensus & precise_confident)
cohort_report = cohort_report %>% mutate(not_high_confidence_rna_dna = !patient_fusion %in% filter(cohort_report,high_confidence_rna_dna)$patient_fusion)


#### 
## uq fusions
#annotate somatic only and germline only based on if precise confident, or of no precise conf exists

if(uq_fusions(cohort_report)!=(uq_fusions(filter(cohort_report,high_confidence_rna_dna))+uq_fusions(filter(cohort_report,not_high_confidence_rna_dna)))) {
  print("Warning: high and low conf sets could not be split")
} else {
  cohort_report_hc = annotate_variant_class_fractions(filter(cohort_report,high_confidence_rna_dna))
  cohort_report_lc = annotate_variant_class_fractions(filter(cohort_report,not_high_confidence_rna_dna))
  cohort_report = rbind(cohort_report_hc,cohort_report_lc) %>% unique()
}

#similar approach for fusion predictions

predictions_helper_fc_any = fusion_overview %>% filter(grepl("fusioncatcher",rna_tools))
predictions_helper_sf_any = fusion_overview %>% filter(grepl("starfusion",rna_tools))

fusion_overview = fusion_overview %>% mutate(rna_consensus = (rna_tools=="fusioncatcher, starfusion" | 
                                              (patient_fusion %in% predictions_helper_sf_any$patient_fusion &  patient_fusion %in% predictions_helper_fc_any$patient_fusion)))


write.table(cohort_report,cohort_report_path,quote = FALSE,sep = "\t",row.names=FALSE)
write.table(fusion_overview,fusion_overview_anno_path,quote = FALSE,sep = "\t",row.names=FALSE)
