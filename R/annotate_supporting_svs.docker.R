### Supporting SV annotation 
## Bring together annotation from utils

## Last update: 2022-02-27 
## Refactor to use SVs anotated by SV pipeline


## cohort suppporting svs + pairwise overlap merged
# retrieve SVs annotated with pairwise overlaps from SV pipeline: SV CNA, population db, and repeats segdups 


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
  
  library(GenomicRanges, quietly=TRUE)
  library(VariantAnnotation, quietly=TRUE)
  library(StructuralVariantAnnotation, quietly=TRUE)
  library(rtracklayer, quietly=TRUE)
  #library(ensembldb)
  
})

source(paste0(script_dir,"functions.general.R"))
source(paste0(script_dir,"functions.svs.R"))


## TODO: not needed?
## Read in cohort 
if(FALSE) {
if(!exists("patient_table_path") | length(Sys.glob(patient_table_path))!=1) {
  print(paste0("Cohort file not available: ",patient_table_path))
  #quit()
}


cohort = read.table(patient_table_path,sep = "\t", header=T,stringsAsFactors = T)
}

#order of arguments matters
map_template_vars=c('${input_dir}'=input_dir,'${output_dir}'=output_dir,
                    '${cohort_identifier}'=cohort_identifier,'${cohort_wdir}'=cohort_wdir)


#fill templates
reports_dir =  stri_replace_all_fixed(reports_dir_template,names(map_template_vars), map_template_vars,vectorize=F)
base_dir =  stri_replace_all_fixed(base_dir_template,names(map_template_vars), map_template_vars,vectorize=F)

map_template_vars=c(map_template_vars,'${reports_dir}'=reports_dir,'${base_dir}'=base_dir)

map_template_vars_merged = c('${analysis_type}'=".merged",
                             map_template_vars )

## Input
cohort_supporting_svs_path = stri_replace_all_fixed(cohort_supporting_svs_path_template,names(map_template_vars_merged), map_template_vars_merged,vectorize=F)
cohort_pairwise_overlap_merged_path = stri_replace_all_fixed(cohort_pairwise_overlap_merged_path_template,names(map_template_vars_merged), map_template_vars_merged,vectorize=F) 

## SV annotated by SV pipeline
cohort_supporting_svs_path = "~/fusion_sq/run/fusion_sq/output/cohort_supporting_svs.merged.anno.tsv"

sv_cna_path = "~/fusion_sq/run/fusion_sq/output/sv_cna.fusion_sq_svs.tsv"

gtf_path = paste0(resources_dir,"ref_annot_",reference,".gtf.gz")
txdb_path = paste0(resources_dir,"ref_annot_",reference,".sqlite")
gencode_annotation_path= paste0(resources_dir,gencode_annotation_file)



#output:
cohort_supporting_svs_anno_path = stri_replace_all_fixed(cohort_supporting_svs_anno_path_template,names(map_template_vars_merged), map_template_vars_merged,vectorize=F)

#population_sv_database_matches_path = paste0(reports_dir,population_sv_database_matches_outfile)
#gene_tx_table_path = paste0(reports_dir,gene_tx_table_outfile)

# Read SVs ----
#refactor 2022-02 replaced  bp_name = patient_id + sv_name by patient_sv_name to make consistent with SV pipeline
# loading the annotated svs from SV pipeline with repeats/segdup/pop db 

supporting_svs_df = read.table(cohort_supporting_svs_path,header=T,sep="\t") 

supporting_svs_anno_df = supporting_svs_df

# Add pairwise overlap merged ----
#to assess quality of sv merged
pairwise_overlap_merged_df = read.table(cohort_pairwise_overlap_merged_path,sep = "\t",header = T)

#make consistent with patient_sv_name
pairwise_overlap_merged_df$set1 = paste0(pairwise_overlap_merged_df$patient_id,"_",pairwise_overlap_merged_df$set1)
pairwise_overlap_merged_df$set2 = paste0(pairwise_overlap_merged_df$patient_id,"_",pairwise_overlap_merged_df$set2)

pairwise_overlap_merged_df = pairwise_overlap_merged_df %>% 
  merge(supporting_svs_anno_df[c("patient_sv_name","coordinate","fusion_name")],by.x="set1",by.y="patient_sv_name") %>%
  dplyr::rename(set1_coordinate=coordinate) %>% 
  merge(supporting_svs_anno_df[c("patient_sv_name","coordinate")],by.x="set2",by.y="patient_sv_name") %>%
  dplyr::rename(set2_coordinate=coordinate)

pairwise_overlap_merged_df$overlap_distance = GenomicRanges::distance(GRanges(pairwise_overlap_merged_df$set1_coordinate),GRanges(pairwise_overlap_merged_df$set2_coordinate))

supporting_svs_anno_df = supporting_svs_anno_df %>% 
  left_join(pairwise_overlap_merged_df[c("set1","overlap_distance","overlap_set1_set2")],by=c("patient_sv_name"="set1")) 

#check supporting_svs_anno_df %>% filter(is.na(overlap_set1_merged) & !is.na(sv_merged) & grepl("merged",sv_merged)) 


## Copy number l2fc ----
cohort_sv_cna = read.table(sv_cna_path,sep="\t",header=T)

if(length(Sys.glob(sv_cna_path))==1){
  
  sv_cna_cols = c("patient_sv_name","cr_l2fc_50_max","maf_50","cr_l2fc_50","start_cr_l2fc_50", "end_cr_l2fc_50", "start_call","end_call")
  sv_cna_cols = sv_cna_cols[sv_cna_cols %in% names(cohort_sv_cna)]

  supporting_svs_anno_df = supporting_svs_anno_df %>% 
    left_join(cohort_sv_cna[,sv_cna_cols],by="patient_sv_name")
  
} else {
  print("No copy number data found")
}


write.table(supporting_svs_anno_df,cohort_supporting_svs_anno_path,quote = F,sep = "\t",row.names=F,col.names = T)


quit()

## Prepare -> convert to GRanges
supporting_svs_anno = GRanges(supporting_svs_anno_df)
names(supporting_svs_anno) = supporting_svs_anno$patient_sv_name


## Annotation resources ----
library(ensembldb)
gtf <- rtracklayer::import(gtf_path)
gene_properties = gtf[gtf$type=="gene"]

## Annotate genes
supporting_svs_anno = annotate_svs_genes(supporting_svs_anno,gene_properties)

# Create TxDb from gtf for the different versions if not exists 
if(length(Sys.glob(txdb_path))<1) {
  txdb = makeTxDbFromGFF(gtf_path,format = "gtf")
  saveDb(txdb, txdb_path)
} else {
  txdb=loadDb(txdb_path)
}
#exons = exonsBy(txdb, "tx", use.names=TRUE) 
introns = intronsByTranscript(txdb, use.names=TRUE) 

gencode_annotation = read.table(gencode_annotation_path,sep="\t",header = F)
tx_mane_select = gencode_annotation[grepl("MANE",gencode_annotation$V2),c("V1")]


### intron/exon annotation ----

if(FALSE) {
  ## Use mane select trancript to annotate introns
  ensembl_id_lst = unique(unlist(strsplit(supporting_svs_anno$ensembl_id,", ")))
  ensembl_id_lst=ensembl_id_lst[grepl("ENSG",ensembl_id_lst)]
  
  transcripts = get_tx_by_genes(gtf,ensembl_id_lst)
  
  ## Aim is single transcript per gene,
  if(length(Sys.glob(gene_tx_table_path))==1){
    gene_tx_table = read.table(gene_tx_table_path,sep = "\t",header=T)
    gene_tx_table = GRanges(gene_tx_table)
  } else {
    gene_tx_table = endoapply(split(transcripts,mcols(transcripts)$gene_id), function(x) filter_tx_canonical(x))
    gene_tx_table = unlist(gene_tx_table)
    
    write.table(gene_tx_table,gene_tx_table_path,quote = F,sep = "\t",row.names=F,col.names = T)
  }
  
  gene_tx_table = gene_tx_table[gene_tx_table$gene_id %in% ensembl_id_lst]
  
  
  introns=introns[names(introns) %in% gene_tx_table$transcript_id]
  introns_per_tx = get_features_per_tx(introns,gene_tx_table$transcript_id)
  metadata_introns = as.data.frame(mcols(introns_per_tx))
  metadata_introns = metadata_introns %>% left_join(as.data.frame(mcols(gene_tx_table)[,c("transcript_id","gene_name")]),by=c("tx_id"="transcript_id"))
  mcols(introns_per_tx)=metadata_introns
  
  supporting_svs_anno = annotate_svs_features(supporting_svs_anno,introns_per_tx) 
  mcols(supporting_svs_anno) = as.data.frame(mcols(supporting_svs_anno)) %>% dplyr::rename(intron_overlap=ft_overlap, intron_tx = tx_id)
}
### end of intron/exon annotation

write.table(supporting_svs_anno_df,cohort_supporting_svs_anno_path,quote = F,sep = "\t",row.names=F,col.names = T)





## March 2022 temporary code and part of SV pipeline: intersect SVs with replication domains breakpoints and ovelaps
#2022-03-01 


###Replication timing ----

get_replication_timing_df  = function(replication_timing_path) {
  replication_timing_df = read.table(replication_timing_path,header=T,sep="\t",comment.char = "")
  
  replication_timing_df = replication_timing_df %>% dplyr::rename(seqnames=CHR,start=POS,rt_anno=switch) %>% dplyr::mutate(end=start+49999) 
  replication_timing_df$rt_id= paste0("rt_",1:nrow(replication_timing_df))
  replication_timing_df$to_coordinate = paste0(replication_timing_df$seqnames,":",replication_timing_df$start,"-",replication_timing_df$end)
  
  return(replication_timing_df)
}
get_reciprocal_overlap_pairs_start_end = function(svs,properties,reciprocal_overlap=0,svtype_matching=F){
  ## starting bp 
  svs_start = svs
  end(svs_start)=start(svs_start)
  
  ## ending bp 
  svs_end = svs
  start(svs_end)=end(svs_end)
  
  start_overlaps = get_reciprocal_overlap_pairs(svs_start,properties,reciprocal_overlap = reciprocal_overlap,svtype_matching = svtype_matching)
  start_overlaps$sv_breakpoint_orientation="start"
  
  end_overlaps = get_reciprocal_overlap_pairs(svs_end,properties,reciprocal_overlap = reciprocal_overlap,svtype_matching = svtype_matching)
  end_overlaps$sv_breakpoint_orientation="end"
  
  overlaps= rbind(start_overlaps,end_overlaps)
  
  return(overlaps)
}


replication_timing_path_template = "~/Documents/resources/ReplicationDomain.20220114.hg38_labels_NatGenetics2018.Human_NormalCellTypes_48datasets_hg38_50KB-qunatilescaled_label_0.15_threshold.txt"
replication_timing_path = stri_replace_all_fixed(replication_timing_path_template,names(map_template_vars), map_template_vars,vectorize=F)
replication_timing_df = get_replication_timing_df(replication_timing_path)
replication_timing = GRanges(replication_timing_df)
names(replication_timing)=replication_timing$rt_id

replication_timing_df_cols=c("rt_id","rt_anno","to_coordinate","L","M","E")


supporting_svs_anno = GRanges(supporting_svs_anno_df)
names(supporting_svs_anno) = supporting_svs_anno$patient_sv_name


#Intersect 

rt_sv_overlaps = get_reciprocal_overlap_pairs(supporting_svs_anno,replication_timing,reciprocal_overlap = 0,svtype_matching = F)
rt_sv_overlaps = rt_sv_overlaps %>% dplyr::rename( patient_sv_name =set1,rt_id=set2)

rt_sv_overlaps = rt_sv_overlaps %>% left_join(supporting_svs_anno_df,by="patient_sv_name") %>% 
  left_join(replication_timing_df[,replication_timing_df_cols],by="rt_id")

sv_anno_cols = names(supporting_svs_anno_df)
rt_sv_overlaps_summary = rt_sv_overlaps %>% 
  #group_by(across(all_of(sv_anno_cols))) %>%
  group_by(patient_sv_name) %>%
  summarize( 
    rt_switch = sum(rt_anno=="S"),
    rt_early = sum(rt_anno=="CE"),
    rt_late = sum(rt_anno=="CL"),
    rt_unkn = sum(rt_anno=="N/A"),
    rt_min_interval=min(to_coordinate),
    rt_max_interval=max(to_coordinate),
  )

rt_sv_overlaps_summary_path_template = "${reports_dir}cohort_supporting_svs.rt_overlaps_summary${analysis_type}.tsv"

rt_sv_overlaps_summary_path = stri_replace_all_fixed(rt_sv_overlaps_summary_path_template,names(map_template_vars_merged), map_template_vars_merged,vectorize=F)

write.table(rt_sv_overlaps_summary,rt_sv_overlaps_summary_path,quote = F,sep = "\t",row.names=F,col.names = T)

## sv bp 

rt_sv_bp_overlaps = get_reciprocal_overlap_pairs_start_end(supporting_svs_anno,replication_timing,reciprocal_overlap = 0,svtype_matching = F)
rt_sv_bp_overlaps = rt_sv_bp_overlaps %>% dplyr::rename( patient_sv_name =set1,rt_id=set2)

rt_sv_bp_overlaps = rt_sv_bp_overlaps %>% left_join(supporting_svs_anno_df[,c("patient_sv_name","sv_merged","patient_id","patient_sv_merged")],by="patient_sv_name") %>% 
  left_join(replication_timing_df[,replication_timing_df_cols],by="rt_id")

rt_sv_bp_overlap_path_template = "${reports_dir}cohort_supporting_svs.rt_sv_bp_overlap${analysis_type}.tsv"
rt_sv_bp_overlap_path = stri_replace_all_fixed(rt_sv_bp_overlap_path_template,names(map_template_vars_merged), map_template_vars_merged,vectorize=F)

write.table(rt_sv_bp_overlaps,rt_sv_bp_overlap_path,quote = F,sep = "\t",row.names=F,col.names = T)

