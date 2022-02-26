
  source("~/fusion_sq/R/default.conf")
  source("~/fusion_sq/R/default.docker.local.conf")
  source("~/fusion_sq/run/fusion_sq/fusion_sq.conf")
  
  
  suppressPackageStartupMessages({
    library(GenomicRanges)
    library(tidyverse, quietly=TRUE)
    library(stringr, quietly=TRUE)
    library(stringi)
    
    library(VennDetail)
    library(ggplot2)
    library(RColorBrewer)
    library(DiagrammeR)
    
  })
  
  #source("R/default.conf") #in script dir
  source(paste0(script_dir,"functions.general.R")) 
  
  source(paste0("~/PycharmProjects/structuralvariation/sv_functional_analysis/functions.general.R"))
  
  cohort = read.table(patient_table_path,sep = "\t", header=T,stringsAsFactors = T)
  cohort$fusion_positive_patient=cohort$fusion_status
  
  
  all_fusion_catcher=read.table("~/data/fusioncatcher/fusioncatcher.20220201.tsv",sep="\t",header=F)
  
  all_fusion_catcher  = all_fusion_catcher %>% separate(col="V1",sep=":",into=c("filename","gup_gene_name")) %>% 
    separate(col="filename",sep="\\.",into=c("rna_id"),extra="drop") %>% 
    dplyr::rename(gdw_gene_name=V2,gup_sf_breakpoint=V3,gdw_sf_breakpoint=V4)
  
  all_fusion_catcher = all_fusion_catcher %>% filter(!grepl("Gene_1",gup_gene_name))
  all_fusion_catcher = all_fusion_catcher %>% left_join(cohort[,c("rna_id","patient_id","fusion_positive_patient")])
  
  all_fusion_catcher$fusion_name=paste0(all_fusion_catcher$gup_gene_name,"--",all_fusion_catcher$gdw_gene_name)
  all_fusion_catcher = all_fusion_catcher %>% mutate(fusion_name=str_replace(fusion_name,fixed("IGH@"),"IGH"))
  
  all_fusion_catcher$patient_fusion=paste0(all_fusion_catcher$patient_id,"_",all_fusion_catcher$fusion_name)
  
  all_fusion_catcher$gup_sf_breakpoint=paste0("chr",all_fusion_catcher$gup_sf_breakpoint)
  all_fusion_catcher$gdw_sf_breakpoint=paste0("chr",all_fusion_catcher$gdw_sf_breakpoint)
  
  candidates = all_fusion_catcher %>% filter(fusion_positive_patient)
  
  cohort_report_old = read.table("~/PycharmProjects/wdl_pipeline/fusion_pilot/fusion_sq_grape/reports_20210411/cohort_report.tsv",sep = "\t", header=T)
  
  cohort_report_old = cohort_report_old %>% separate_rows(gup_sf_breakpoint, sep=", ")%>% separate_rows(gdw_sf_breakpoint, sep=", ") 
  cohort_report_old = cohort_report_old %>% mutate(patient_fusion=str_replace(patient_fusion,fixed("IGH-@-ext"),"IGH"))
  
  
  clinrel_fusions=cohort_report_old %>% filter(clinically_validated|clinically_validated_reciprocal)  %>%
    filter(precise_confident) %>% select(patient_id,patient_fusion,fusion_name,gup_sf_breakpoint,gdw_sf_breakpoint) %>% unique()
  
  #fusion is not there but patient is
  
  clinrel_fusions %>% filter(!patient_fusion %in% candidates$patient_fusion & 
                               patient_id %in% candidates$patient_id)
  
  same_predictions = 
    merge(clinrel_fusions[,c("patient_fusion","gup_sf_breakpoint","gdw_sf_breakpoint")],
          candidates[,c("patient_fusion","gup_sf_breakpoint","gdw_sf_breakpoint")])
  
  #always same predictions
  clinrel_fusions %>% filter(patient_fusion %in% candidates$patient_fusion & 
                               !patient_fusion %in% same_predictions$patient_fusion) 
  
  #%>%   left_join(candidates[,c("patient_fusion","gup_sf_breakpoint","gdw_sf_breakpoint")] %>% unique(),by="patient_fusion")
  
  #missing patients => look into archived files
  clinrel_fusions %>% filter(!patient_fusion %in% candidates$patient_fusion & 
                               !patient_id %in% candidates$patient_id) %>% 
    left_join(cohort[,c("rna_id","patient_id")])
  # 
  # 
  # 
  # grep KIAA1549 PMABM000CDM.final-list_candidate-fusion-genes.txt | cut -f 1,2,9,10
  #   BRAF	KIAA1549	7:140781576:-	7:138861456:-
  #   KIAA1549	BRAF	7:138861139:-	7:140787584:-[yes]
  # 
  # [ivanbelzen@gerrit archive]$ grep KIAA1549 PMABM000BHY.final-list_candidate-fusion-genes.txt | cut -f 1,2,9,10
  # KIAA1549	BRAF	7:138861139:-	7:140787584:-[yes]
  #   [ivanbelzen@gerrit archive]$ grep KIAA1549 PMABM000BEM.final-list_candidate-fusion-genes.txt | cut -f 1,2,9,10
  # BRAF	KIAA1549	7:140776912:-	7:138871362:-
  #   KIAA1549	BRAF	7:138861139:-	7:140787584:- [yes]
  #   [ivanbelzen@gerrit archive]$ grep KIAA1549 PMABM000DIC.final-list_candidate-fusion-genes.txt | cut -f 1,2,9,10
  # [ivanbelzen@gerrit archive]$ grep ETV6 PMABM000DIC.final-list_candidate-fusion-genes.txt | cut -f 1,2,9,10
  # ETV6	RUNX1	12:11869969:+	21:34892963:- [ yes]
  #   ETV6	RUNX1	12:11869969:+	21:34887096:- [yes]
  #   RUNX1	ETV6	21:35048842:-	12:11839140:+ [yes]
  #   RUNX1	ETV6	21:35048842:-	12:11826402:+
  #   RUNX1	ETV6	21:34926215:-	12:11826402:+
  #   RUNX1	ETV6	21:35048842:-	12:11825644:+
  #   RUNX1	ETV6	21:35048842:-	12:11802974:+
  #   
  #   [ivanbelzen@gerrit archive]$ grep MYC PMABM000GPN.final-list_candidate-fusion-genes.txt | grep IGH | cut -f 1,2,9,10
  # IGH@	MYC	14:105743070:+	8:127735873:+ [yes]
  #   IGH@	MYC	14:105745508:+	8:127735796:+
  #   IGH@	MYC	14:105745203:+	8:127735796:+
  #   IGH@	MYC	14:105627807:+	8:127735796:+
  #   IGH@	MYC	14:105834561:-	8:127736830:+
  #   IGH@	MYC	14:105854656:-	8:127736322:+
  #   MYC	IGH@	8:127736871:+	14:105834584:-
  #   IGH@	MYC	14:106511162:-	8:127738128:+
  # 
  #     
  #     only 1 missing
  #   [ivanbelzen@gerrit archive]$ grep KIAA1549 PMABM000CDC.final-list_candidate-fusion-genes.txt | cut -f 1,2,9,10
  #   grep: PMABM000CDC.final-list_candidate-fusion-genes.txt: No such file or directory
  #   it is there for bowtie
  #   PMABM000CDC/results/candidate_fusion_genes_summary_BOWTIE.txt:KIAA1549	BRAF	ENSG00000122778	ENSG00000157764	ENSE00001088064	ENSE00003680515	7:138861139:-	7:140787584:-	41	13	30BOWTIE	GTCCTTCTACAGCCCAGCCCAGACGGCCAACAATCCCTGCAGT*GACTTGATTAGAGACCAAGGATTTCGTGGTGATGGAGGATCAA
  #   
  #   
  
  ### high conf candidates ----
  
  uq_fusions_old = read.table("~/PycharmProjects/wdl_pipeline/fusion_pilot/fusion_sq_grape/reports_20210411/uq_patient_fusions.tsv",sep = "\t", header=T)
  
  uq_fusions_old = uq_fusions_old %>% separate_rows(gup_sf_breakpoint, sep=", ")%>% separate_rows(gdw_sf_breakpoint, sep=", ") 
  uq_fusions_old = uq_fusions_old %>% mutate(patient_fusion=str_replace(patient_fusion,fixed("IGH-@-ext"),"IGH"))
  
  high_confidence_old = uq_fusions_old %>% filter(variant_type=="somatic" & precise_confident) %>% as.data.frame()
  high_confidence_old %>% uq_fusions()
  
  high_confidence_df = merge(high_confidence_old[,c("patient_fusion","gup_sf_breakpoint","gdw_sf_breakpoint")],
                     all_fusion_catcher[,c("patient_fusion","gup_sf_breakpoint","gdw_sf_breakpoint")]) %>% as.data.frame()
  high_confidence_df$patient_fusion %>% unique() %>% length()
  
  high_confidence_old %>% dplyr::filter(patient_fusion %in% high_confidence_df$patient_fusion)  %>% uq_fusions()
  
  high_confidence_missing = high_confidence_old %>% filter(!patient_fusion %in% high_confidence_df$patient_fusion & 
                                   patient_id %in% all_fusion_catcher$patient_id) 
  high_confidence_missing$patient_fusion %>% unique() %>% length()
  
  high_confidence_missing %>% select(patient_fusion,annotation,gup_sf_breakpoint,gdw_sf_breakpoint,contains("anno")) %>% View()
  
  high_confidence_old %>% filter(patient_id=="PMCID340AAO") %>% filter(patient_fusion=="PMCID340AAO_MYC--IGH" | patient_fusion=="PMCID340AAO_IGH--MYC" )
  #PMCID340AAO_MYC--IGH	cancer	chr8:127736623:+	chr14:105862211:+
  #matching interval from star fusion 1	FALSE	FALSE	ENSG00000136997.19		chr8:127735434-127742951:+		chr8:127736624-127738291:+		chr8:127735838-127738291:+	chr14:105856220-105863197:+	chr8:127736123-127737122:+	chr14:105861711-105862710:+	chr8:127736623-127742951:+		PMCID340AAO_PMABM000HDT_1
  #intron_consensus           sj  
  #sj for igh = chr14:105856220-105863197:+ also in fusion catcher 
  # and together with adj intron probably chr8:127736624-127738291:+
  ## not that exact combo so it will be curious to see 
  all_fusion_catcher %>% filter(patient_fusion=="PMCID340AAO_MYC--IGH" | patient_fusion=="PMCID340AAO_IGH--MYC" )
  all_fusion_catcher %>% filter(patient_id=="PMCID340AAO") %>% nrow()
  	
  
  all_fusion_catcher %>% filter(patient_fusion=="PMCID003AAM_TK1--BRCA1")
  all_fusion_catcher %>% filter(patient_id=="PMCID003AAM") %>% 
    filter(grepl("BRCA",fusion_name) | grepl("TK1",fusion_name) | 
             gup_sf_breakpoint=="chr17:78175529:-" | gdw_sf_breakpoint=="chr17:43063951:-" | 
             grepl("chr17:430639",gdw_sf_breakpoint) | grepl("chr17:781755",gup_sf_breakpoint))
  
  
  
  ### compare with fusion overview of star fusion 
  
  fusions_overview_old = read.table("~/PycharmProjects/wdl_pipeline/fusion_pilot/fusion_sq_grape/reports_20210411/fusion_overview.anno.tsv",sep = "\t", header=T)
  
  fusions_overview_old = fusions_overview_old %>% separate_rows(gup_sf_breakpoint, sep=", ")%>% separate_rows(gdw_sf_breakpoint, sep=", ") 
  fusions_overview_old = fusions_overview_old %>% mutate(patient_fusion=str_replace(patient_fusion,fixed("IGH-@-ext"),"IGH"))
  
  fusions_overview_old = fusions_overview_old %>% mutate(  inter_chr= grepl("INTER",annots))
  
  predictions_old = fusions_overview_old %>% 
    filter(JunctionReadCount>2&SpanningFragCount>2&FFPM>0.1) %>%
    filter((anno_cancer_chimera|anno_cancer_gene_db)&!anno_healthy_chimera&!anno_neighbors_overlap)
  
  prediction_candidates = merge(predictions_old[,c("patient_fusion","gup_sf_breakpoint","gdw_sf_breakpoint")],
                     all_fusion_catcher[,c("patient_fusion","gup_sf_breakpoint","gdw_sf_breakpoint")])
  
  
  all_tssv = merge(cohort_report_old[,c("patient_fusion","gup_sf_breakpoint","gdw_sf_breakpoint")],
                             all_fusion_catcher[,c("patient_fusion","gup_sf_breakpoint","gdw_sf_breakpoint")]) %>% as.data.frame()
  all_tssv$patient_fusion %>% unique() %>% length()
  
  all_tssv = all_tssv %>% left_join(cohort_report_old[,c("patient_fusion","gup_gene_id","gdw_gene_id","patient_id")])
  
 
  fusions_overview_old %>% 
    filter(patient_fusion %in% prediction_candidates$patient_fusion) %>% 
    filter(!patient_fusion %in% all_tssv$patient_fusion) %>% 
    anti_join(all_tssv,by=c("patient_id","gup_gene_id"="gdw_gene_id","gdw_gene_id"="gup_gene_id")) %>% 
    select(patient_fusion,inter_chr,gup_sf_breakpoint,gdw_sf_breakpoint,contains("anno")) %>% unique() %>% View()
  
  
  
  
  
  
  ### Cases ----
  
  fusioncatcher_output = read.table("~/fusion_sq/run/fusion_sq/output/fusion_level_results.fusioncatcher.PMCID203AAL_PMABM000HAT.tsv",
                                    sep="\t",header=T)
  
  fusioncatcher_output$patient_id=patient$patient_id
  fusioncatcher_output = fusioncatcher_output %>% 
    mutate(patient_fusion =  paste(patient_id,fusion_name,sep="_"),
           patient_fusion_sv =  paste(patient_id,gup_sv_merged,gdw_sv_merged,fusion_name,sep="_"),
           patient_sv_id =  paste(patient_id,gup_sv_merged,gdw_sv_merged,sep="_"))
  
  
  intersect(
    starfusion_output$patient_fusion_sv,
    fusioncatcher_output$patient_fusion_sv)
  
  
  starfusion_output %>% filter(!patient_fusion_sv %in% fusioncatcher_output$patient_fusion_sv) %>% 
    filter(fusion_name %in% fusioncatcher_output$fusion_name)
  
  
  fusioncatcher_output %>% filter(!patient_fusion_sv %in% starfusion_output$patient_fusion_sv)
  
  
  
  old_starfusion_output = read.table("~/PycharmProjects/wdl_pipeline/fusion_pilot/fusion_sq_grape/reports_20210411/fusion_level_results.PMCID203AAL_PMABM000HAT.tsv",
                                     sep="\t",header=T)
  
  old_starfusion_output = old_starfusion_output %>% mutate(precise_confident = location_precise&grepl(",",tools))
  
  compare_cols=c("fusion_name","gup_sf_breakpoint","gdw_sf_breakpoint","precise_confident",
                 "tools","sv_names","tumor_af_mean","svlen","svlen_spread")
  compare_cols=c(compare_cols,names(starfusion_output)[grepl("merged",names(starfusion_output))])
  
  combined_starfusion = starfusion_output[,compare_cols] %>% 
    merge(old_starfusion_output[,compare_cols], by=c("fusion_name","gup_sf_breakpoint","gdw_sf_breakpoint"))
  
  combined_starfusion %>% filter(sv_names.x==sv_names.y & precise_confident.x!=precise_confident.y)
  
  combined_starfusion %>% filter(sv_names.x!=sv_names.y & precise_confident.x==T) %>%
    filter(!fusion_name %in% filter(combined_starfusion,sv_names.x==sv_names.y)$fusion_name) %>%
    select(fusion_name,sv_names.x,sv_names.y,tumor_af_mean.x,tumor_af_mean.y) %>% View()
  
  
  #became false
  combined_starfusion %>% filter(precise_confident.y==T & precise_confident.x==F) %>%
    filter(!fusion_name %in% filter(combined_starfusion,precise_confident.x==T)$fusion_name) %>%
    select(fusion_name,sv_names.x,sv_names.y,tumor_af_mean.x,tumor_af_mean.y) 
  
  #became true
  combined_starfusion %>% filter(precise_confident.y==F & precise_confident.x==T) %>%
    filter(!fusion_name %in% filter(combined_starfusion,precise_confident.y==T)$fusion_name) %>%
    select(fusion_name,sv_names.x,sv_names.y,tumor_af_mean.x,tumor_af_mean.y) 
  
  #still looks fine, is because of sorting differences and in some cases what were 2 entries became 1
  #tumor af can change because of grouping differences => have rerun all 
  #no qualitative differences (T -> F or F -> T )
  
  #changes in overlap frac and bp distance?
  
  combined_starfusion %>% filter(sv_names.x!=sv_names.y & precise_confident.x==T) %>%
    filter(!fusion_name %in% filter(combined_starfusion,sv_names.x==sv_names.y)$fusion_name) %>%
    select(fusion_name,sv_names.x,sv_names.y,tumor_af_mean.x,tumor_af_mean.y,svlen_spread.x,svlen_spread.y) %>% View()
  
  
  
  

## Single patient comparison of new vs old integration of SVs => no difference so this is OK
  
  sf_fusion_level_results_path = "~/fusion_sq/run/fusion_sq/output/fusion_level_results.starfusion.PMCID203AAL_PMABM000HAT.tsv"
  patient_starfusion_fusion_level = read.table(sf_fusion_level_results_path,sep="\t",header=T)
  
  patient_starfusion_fusion_level$patient_id="PMCID203AAL"
  patient_starfusion_fusion_level = patient_starfusion_fusion_level %>% 
    mutate(patient_fusion =  paste(patient_id,fusion_name,sep="_"),
           patient_fusion_sv =  paste(patient_id,gup_sv_merged,gdw_sv_merged,fusion_name,sep="_"),
           patient_sv_id =  paste(patient_id,gup_sv_merged,gdw_sv_merged,sep="_"))
  
  
  
  patient_starfusion_fusion_level_old = read.table("~/PycharmProjects/wdl_pipeline/fusion_pilot/fusion_sq_grape/reports_20210411/fusion_level_results.PMCID203AAL_PMABM000HAT.tsv",
                                               sep="\t",header=T)
  patient_starfusion_fusion_level_old$patient_id="PMCID203AAL"
  
  patient_starfusion_fusion_level_old = patient_starfusion_fusion_level_old %>% select(-precise_confident) %>% dplyr::rename(precise_confident = precise_confident_sv)
  
  patient_starfusion_fusion_level_old = patient_starfusion_fusion_level_old %>% 
    mutate(patient_fusion =  paste(patient_id,fusion_name,sep="_"),
           patient_fusion_sv =  paste(patient_id,gup_sv_merged,gdw_sv_merged,fusion_name,sep="_"),
           patient_sv_id =  paste(patient_id,gup_sv_merged,gdw_sv_merged,sep="_"))

    
  patient_starfusion_fusion_level_old %>% filter(precise_confident & 
                                         !patient_fusion %in% filter(patient_starfusion_fusion_level,precise_confident)$patient_fusion  )

  patient_starfusion_fusion_level %>% filter(precise_confident & 
                                                   !patient_fusion %in% filter(patient_starfusion_fusion_level_old,precise_confident)$patient_fusion  )
  
  
  ## Script to integrate to seperate file 
  
  #goal is integrate
  #add column rna_tools match on patient fusion sv
  
  fc_fusion_level_results_path = "~/fusion_sq/run/fusion_sq/output/fusion_level_results.fusioncatcher.PMCID203AAL_PMABM000HAT.tsv"
  sf_fusion_level_results_path = "~/fusion_sq/run/fusion_sq/output/fusion_level_results.starfusion.PMCID203AAL_PMABM000HAT.tsv"
  fusion_level_results_path = "~/fusion_sq/run/fusion_sq/output/fusion_level_results.merged.PMCID203AAL_PMABM000HAT.tsv"
  
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
  
  #fusion_level_df %>% View()
  
  write.table(fusion_level_df,merged_fusion_level_results_path,col.names = T,row.names = F,sep="\t")
  
  
### Checks ----
  nrow(fusion_level_df)
  
  nrow(patient_starfusion_fusion_level)
  nrow(patient_fusioncatcher_fusion_level)
  
  
  #that is still paired on uq sv level means something! sometimes also subset of each other 
  fusion_level_df %>% filter(gup_sf_breakpoint_fusioncatcher!=gup_sf_breakpoint_starfusion | gdw_sf_breakpoint_fusioncatcher!=gdw_sf_breakpoint_starfusion)
  
  #map patient fusion to svs then look if differences arise there 
  
  fusion_sv_map = fusion_level_df %>% select(patient_fusion,patient_sv_id,rna_tools,gup_sv_merged,gdw_sv_merged,gup_sv_merged_coordinate,gdw_sv_merged_coordinate,sv_names) %>% unique()
  
  fusion_count_svs  = fusion_sv_map %>% group_by(patient_fusion,rna_tools) %>% summarize(svs_combi_cnt = length(unique(patient_sv_id)))
  
  fusion_count_svs %>% group_by(patient_fusion ) %>% summarize(rna_support_cnt = length(unique(rna_tools))) %>% filter(rna_support_cnt>1)
  
  
  #not sv but fusion 
  patient_fusioncatcher_fusion_level %>% filter(!patient_fusion_sv %in% merged_df$patient_fusion_sv & 
                                                  patient_fusion %in% filter(patient_starfusion_fusion_level)$patient_fusion  )
  
  patient_starfusion_fusion_level %>% filter(!patient_fusion_sv %in% merged_df$patient_fusion_sv & 
                                               patient_fusion %in% filter(patient_fusioncatcher_fusion_level)$patient_fusion  )
  
  
  
  completed_patients = read.table("~/fusion_sq/fusion_level.completed.lst")
  
  completed_patients = completed_patients %>% separate(col="V1",sep="\\.",into=c("filename","rna_tool","patient_identfier","ext"))
  completed_patients = completed_patients %>% separate(col="patient_identfier", sep="_.",into=c("patient_id","rna_id"),remove=F)
  
  cohort_report_old = read.table("~/PycharmProjects/wdl_pipeline/fusion_pilot/fusion_sq_grape/reports_20210411/cohort_report.tsv",sep = "\t", header=T)
  
  cohort_report_old = cohort_report_old %>% separate_rows(gup_sf_breakpoint, sep=", ")%>% separate_rows(gdw_sf_breakpoint, sep=", ") 
  cohort_report_old = cohort_report_old %>% mutate(patient_fusion=str_replace(patient_fusion,fixed("IGH-@-ext"),"IGH"))
  
  cohort_report_old %>% filter(!patient_id %in% filter(completed_patients,rna_tool=="starfusion")$patient_id) %>% select(patient_id) %>% unique()

  cohort = read.table(patient_table_path,sep = "\t", header=T,stringsAsFactors = T)

  cohort %>% filter(!patient_id %in% filter(completed_patients,rna_tool=="fusioncatcher")$patient_id) %>% select(patient_id) %>% unique()
  
  