## Fusion-sq 
## Annotate cohort report
## Last update: 2022-01-12
##  path local overrides and docker version
## 
# Input: cohort fusions  (from collect_cohort.R)
## 
# Annotate gene fusions with cancer related genes, healthy chimera, chromosome bands/cytobands, 
# Annotate gene fusions with properties of underlying SVs (see utils)
## SV overlaps with population SV databases, repeats, segmental duplications, introns

#Optional - integration with different modalities: 
## Annotate with CNAs and SNVs (see utils)
## Gene expression & z scores (from gene_expression_zscores.R), wilcox test for fusion carrying patients



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

if(!exists("patient_table_path") | length(Sys.glob(patient_table_path))!=1) {
  print(paste0("Cohort file not available: ",patient_table_path))
  #quit()
}


if(!exists("analysis_type")) {
  analysis_type="starfusion"
}

suppressPackageStartupMessages({
  library(tidyverse, quietly=TRUE)
  library(stringr, quietly=TRUE)
  library(GenomicRanges)
  library(dplyr)
  library(stringi)
})

#source("R/default.conf") #in script dir
source(paste0(script_dir,"functions.general.R")) 


##  make templates where you either fill ${analysis_type}=> analysis_type,"." if fusion catcher and temporarily if star fusion just "" 
analysis_type_filename=paste0(".",analysis_type)

if(analysis_type=="starfusion"){
  analysis_type_filename=""
} 

#order of arguments matters
map_template_vars=c('${input_dir}'=input_dir,'${output_dir}'=output_dir,
                    '${cohort_identifier}'=cohort_identifier,'${cohort_wdir}'=cohort_wdir,
                    '${analysis_type}'=analysis_type_filename)


#fill templates
reports_dir =  stri_replace_all_fixed(reports_dir_template,names(map_template_vars), map_template_vars,vectorize=F)
base_dir =  stri_replace_all_fixed(base_dir_template,names(map_template_vars), map_template_vars,vectorize=F)

map_template_vars=c(map_template_vars,'${reports_dir}'=reports_dir,'${base_dir}'=base_dir)

## Inputs defined in config file
## Collected for cohort:
cohort_fusion_level_results_path = stri_replace_all_fixed(cohort_fusion_level_results_path_template,names(map_template_vars), map_template_vars,vectorize=F)
cohort_results_path = stri_replace_all_fixed(cohort_results_path_template,names(map_template_vars), map_template_vars,vectorize=F) 
fusion_overview_path = stri_replace_all_fixed(fusion_overview_path_template,names(map_template_vars), map_template_vars,vectorize=F) 

## supporting SVs with gene annotation
cohort_supporting_svs_anno_path = stri_replace_all_fixed(cohort_supporting_svs_anno_path_template,names(map_template_vars), map_template_vars,vectorize=F)


#Gene expression & z scores 
gene_expression_analysis_path = paste0(reports_dir,gene_expression_zcore_outfile)

# CNA and SNV files => see config

#External resources => see config for files
## TODO replace by the get_cancer_genes() function
chimerseq_path = paste0(resources_dir,chimerseq_file)
mitelman_path = paste0(resources_dir,mitelman_file)
cosmic_path = paste0(resources_dir,cosmic_file)
kinases_path = paste0(resources_dir,kinases_file)
grobner_recurrent_path = paste0(resources_dir,grobner_recurrent_file)
grobner_druggable_path = paste0(resources_dir,grobner_druggable_file)
oncokb_path = paste0(resources_dir,oncokb_file)
chromosome_bands_path = paste0(resources_dir,chromosome_bands_file)
transcriptionfactors_path = paste0(resources_dir,transcriptionfactors_file)

#clinically_relevant_fusion_genes_path = paste0(resources_dir,"InclusionList_fusion_genes_Homo_sapiens_GRCh38_BioMart_v1.8.txt")
clinically_relevant_fusion_genes_path = paste0("~/Documents/resources/InclusionList_fusion_genes_Homo_sapiens_GRCh38_BioMart_v1.8.txt")


## Output
cohort_report_path = stri_replace_all_fixed(cohort_report_path_template,names(map_template_vars), map_template_vars,vectorize=F)
fusion_overview_anno_path = stri_replace_all_fixed(fusion_overview_anno_path_template,names(map_template_vars), map_template_vars,vectorize=F) 

#Temporary code:
#not needed anymore if dots are added to proper location in template
analysis_type_filename_patient=paste0(analysis_type,".")
if(analysis_type=="starfusion"){
  analysis_type_filename_patient=""
} 

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
cohort_results = read.table(cohort_results_path,sep = "\t", header=T)
fusion_overview = read.table(fusion_overview_path,sep = "\t", header=T)


svs_extra_cols=c("FILTER","overlap_distance","overlap_set1_set2","overlap_set1_merged",
                 "start_copy_ratio_log2_mean","end_copy_ratio_log2_mean","copy_ratio_log2_mean","copy_ratio_log2_sd",
                 "repeat_family_start","repeat_family_end","segdup_start","segdup_end",
                 "intron_overlap","intron_tx",
                 "anno_sv_population","anno_sv_db_nstd166","anno_sv_db_nstd186","anno_sv_db_dgv")

#also used for tx selection any time merging on that level is required
patient_sv_summary_cols=c("fusion_name", "gup_sv_merged", "gdw_sv_merged", "specific_sv", "patient_id")
transcript_cols = c("transcript_type","involved_fragment","involved_fraction","involved_exon","max_exon")


## Start annotating fusion overview and then merge with fusion level results to make the report

### Fusion overview to cytobands ----
fusion_overview = fusion_overview %>%
  left_join(cohort %>% select(patient_id, tumor_type_label, primary_group, primary_group_label, 
                                        supergroup, supergroup_label,primary_group_shorthand_label), 
            by="patient_id") 

cohort_fusion_level_results = cohort_fusion_level_results %>% filter(patient_id %in% cohort$patient_id)
cohort_results = cohort_results %>% filter(patient_id %in% cohort$patient_id)
fusion_overview = fusion_overview  %>% filter(patient_id %in% cohort$patient_id)

## cytobands
chromosome_bands = read.table(chromosome_bands_path,header=T,sep="\t",comment.char = "")
chromosome_bands = chromosome_bands %>% dplyr::rename("seqnames"="X.chrom","start"="chromStart","end"="chromEnd","cytoband"="name") 
chromosome_bands$cytoband = paste0(chromosome_bands$seqnames,chromosome_bands$cytoband)
chromosome_bands = GRanges(chromosome_bands)

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
## supporting SVs with gene annotation
if(length(Sys.glob(cohort_supporting_svs_anno_path))!=1 ){ 
  print(paste0("No annotated SVs found: ",cohort_supporting_svs_anno_path))
  #quit() 
} else {
  supporting_sv_anno_df = read.table(cohort_supporting_svs_anno_path,sep = "\t",header=T)



#subset to existing cols
svs_extra_cols_missing=svs_extra_cols[!svs_extra_cols %in% names(supporting_sv_anno_df)]
svs_extra_cols=svs_extra_cols[svs_extra_cols %in% names(supporting_sv_anno_df)]

cohort_results = cohort_results %>% left_join(supporting_sv_anno_df[,c("sv_name","patient_id",svs_extra_cols)],by=c("gup_sv_name"="sv_name","patient_id")) %>%
  rename_at(svs_extra_cols,function(x){paste0("gup_",x)}) %>% 
  left_join(supporting_sv_anno_df[,c("sv_name","patient_id",svs_extra_cols)],by=c("gdw_sv_name"="sv_name","patient_id")) %>%
  rename_at(svs_extra_cols,function(x){paste0("gdw_",x)})

if(length(svs_extra_cols_missing)>0) {
  cohort_results[,paste0("gup_",svs_extra_cols_missing)]=NA
  cohort_results[,paste0("gdw_",svs_extra_cols_missing)]=NA
}

### Summary of supporting SVs per merged  ----

cohort_results$gup_repeat_family = paste0("start:",cohort_results$gup_repeat_family_start,";end:",cohort_results$gup_repeat_family_end)
cohort_results$gup_repeat_family=str_replace(cohort_results$gup_repeat_family,"start:NA;","")
cohort_results$gup_repeat_family=str_replace(cohort_results$gup_repeat_family,"end:NA","")
cohort_results$gdw_repeat_family = paste0("start:",cohort_results$gdw_repeat_family_start,";end:",cohort_results$gdw_repeat_family_end)
cohort_results$gdw_repeat_family=str_replace(cohort_results$gdw_repeat_family,"start:NA;","")
cohort_results$gdw_repeat_family=str_replace(cohort_results$gdw_repeat_family,"end:NA","")

cohort_results$gup_segdup = paste0("start:",cohort_results$gup_segdup_start,";end:",cohort_results$gup_segdup_end)
cohort_results$gup_segdup=str_replace(cohort_results$gup_segdup,"start:NA;","")
cohort_results$gup_segdup=str_replace(cohort_results$gup_segdup,"end:NA","")
cohort_results$gdw_segdup = paste0("start:",cohort_results$gdw_segdup_start,";end:",cohort_results$gdw_segdup_end)
cohort_results$gdw_segdup=str_replace(cohort_results$gdw_segdup,"start:NA;","")
cohort_results$gdw_segdup=str_replace(cohort_results$gdw_segdup,"end:NA","")

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
           

fusion_level_extra_sv_properties$repeat_family = str_replace(fusion_level_extra_sv_properties$repeat_family,"gup:;","")
fusion_level_extra_sv_properties$repeat_family = str_replace(fusion_level_extra_sv_properties$repeat_family,"gdw:;","")
fusion_level_extra_sv_properties$segdup = str_replace(fusion_level_extra_sv_properties$segdup,"gup:;","")
fusion_level_extra_sv_properties$segdup = str_replace(fusion_level_extra_sv_properties$segdup,"gdw:;","")


cohort_fusion_level_results = cohort_fusion_level_results %>% 
  left_join(fusion_level_extra_sv_properties,by=patient_sv_summary_cols)

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
## Chimera annotation


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

## Cancer related gene annotation
#for onco/TSG/fusions together with oncokb and grobner
cancer_related_genes = get_cancer_genes()

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



## properties from STAR Fusion
if(analysis_type=="starfusion") {
  fusion_overview = fusion_overview %>%
    mutate(gup_pfam_kinase = grepl("kinase",PFAM_LEFT,fixed=T)) %>%
    mutate(gdw_pfam_kinase = grepl("kinase",PFAM_RIGHT,fixed=T))
  
  fusion_overview = fusion_overview %>% mutate(anno_neighbors_overlap = grepl("NEIGHBORS_OVERLAP",annots)) 
  
  fusion_overview = fusion_overview %>% mutate(anno_in_any_cancer = 
                                                 (anno_in_any_cancer | anno_starfusion_cancer_db ))
    
}

## Annotate clinically_relevant fusions if annotation is present

if("fusion" %in% names(cohort)) {
  #fusion == true 
  fusion_overview = fusion_overview %>% mutate(patient_clinrel_fusion = patient_id %in% filter(cohort,fusion_status)$patient_id)
  #check fusion_overview %>% filter(patient_clinrel_fusion) %>% uq_patients() == filter(cohort,fusion_status) %>% uq_patients()
  
  fusion_overview$clinically_validated=F
  fusion_overview$clinically_validated_reciprocal=F
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

## Merge cohort fusion level results with the annotation from fusion overview ----
fusion_overview_anno_only = fusion_overview  %>% select(patient_fusion, tumor_type_label, primary_group, primary_group_label, supergroup, supergroup_label, primary_group_shorthand_label,
                                                        contains("anno"),contains("fpkm"),contains("fpkm_log"),contains("snv"),contains("cytoband"),
                                                        contains("clinrel"),contains("clinically_validated")) %>% unique() 

cohort_report = cohort_fusion_level_results %>% 
  left_join(fusion_overview_anno_only, by="patient_fusion")

cohort_report=cohort_report %>% mutate(tumor_normal_diff_mean = (tumor_af_mean-normal_af_mean))

## Summarise over rna level predictions

if(analysis_type=="starfusion") { 
  #if annotated as pfam kinase any time
  cohort_report = cohort_report %>% 
    mutate(anno_gup_pfam_kinase = patient_fusion %in% filter(fusion_overview,gup_pfam_kinase)$patient_fusion) %>% 
    mutate(anno_gdw_pfam_kinase = patient_fusion %in% filter(fusion_overview,gdw_pfam_kinase)$patient_fusion)
}

cohort_report = cohort_report %>%
  left_join( cohort %>% select(patient_id,tumor_id,normal_id,rna_id), by="patient_id")


#summmary labels 
cohort_report=cohort_report %>% annotate_labels_cancer_common() %>% annotate_labels_variant_type() %>% annotate_labels_genes()

##refactor labels
cohort_report=cohort_report %>% mutate(
  high_confidence=precise_confident,
  anno_has_onco_or_tsg = (anno_has_oncogene | anno_has_tsg)#,
  #anno_clinically_relevant = (clinically_validated|clinically_validated_reciprocal)
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
    
    if(nrow(gene_expression_analysis[gene_expression_analysis$ensembl_id==candidate$gup_ensembl_id,])>0){
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
    
    if(nrow(gene_expression_analysis[gene_expression_analysis$ensembl_id==candidate$gdw_ensembl_id,])>0){
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
  
  cohort_report = cohort_report %>% left_join(expression_wilcoxon,by="patient_fusion")
  
}



#annotate somatic only and germline only based on if precise confident, or of no precise conf exists
cohort_report = cohort_report %>% mutate(not_precise_confident = !patient_fusion %in% filter(cohort_report,precise_confident)$patient_fusion)

if(uq_fusions(cohort_report)!=(uq_fusions(filter(cohort_report,precise_confident))+uq_fusions(filter(cohort_report,not_precise_confident)))) {
  print("Warning: high and low conf sets could not be split")
} else {
  cohort_report_hc = annotate_variant_class_fractions(filter(cohort_report,precise_confident))
  cohort_report_lc = annotate_variant_class_fractions(filter(cohort_report,not_precise_confident))
  cohort_report = rbind(cohort_report_hc,cohort_report_lc)
}


write.table(cohort_report,cohort_report_path,quote = FALSE,sep = "\t",row.names=FALSE)
write.table(fusion_overview,fusion_overview_anno_path,quote = FALSE,sep = "\t",row.names=FALSE)


