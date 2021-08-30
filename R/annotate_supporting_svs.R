### Supporting SV annotation 
## Bring together annotation from utils
## Refactor 2021-08

suppressPackageStartupMessages({

library(tidyverse, quietly=TRUE)
library(stringr, quietly=TRUE)
library(stringdist, quietly=TRUE)

library(GenomicRanges, quietly=TRUE)
library(VariantAnnotation, quietly=TRUE)
library(StructuralVariantAnnotation, quietly=TRUE)
library(rtracklayer, quietly=TRUE)
library(ensembldb)

  })
#if("dplyr" %in% (.packages())){
#detach("package:dplyr", unload=TRUE) 
#  detach("package:plyr", unload=TRUE) 
#} 
#library(plyr)
#library(dplyr)

source("R/default.conf") #in script dir
source(paste0(script_dir,"functions.general.R")) 
source(paste0(script_dir,"functions.svs.R"))

## Input

cohort_supporting_sv_path = paste0(reports_dir,cohort_supporting_svs_outfile)
cohort_pairwise_overlap_merged_path = paste0(reports_dir,cohort_pairwise_overlap_merged_outfile)

## Output
cohort_supporting_sv_anno_path = paste0(reports_dir,cohort_supporting_sv_anno_outfile)
population_sv_database_matches_path = paste0(reports_dir,population_sv_database_matches_outfile)

gene_tx_table_path = paste0(reports_dir,gene_tx_table_outfile)

## Annotation resources

gtf_path = paste0(resources_dir,"ref_annot_",reference,".gtf.gz")
txdb_path = paste0(resources_dir,"ref_annot_",reference,".sqlite")
gencode_annotation_path= paste0(resources_dir,gencode_annotation_file)

# Create TxDb from gtf for the different versions if not exists 
if(length(Sys.glob(txdb_path))<1) {
  txdb = makeTxDbFromGFF(gtf_path,format = "gtf")
  saveDb(txdb, txdb_path)
} else {
  txdb=loadDb(txdb_path)
}
#exons = exonsBy(txdb, "tx", use.names=TRUE) 
introns = intronsByTranscript(txdb, use.names=TRUE) 

gtf <- rtracklayer::import(gtf_path)
gene_properties = gtf[gtf$type=="gene"]

gencode_annotation = read.table(gencode_annotation_path,sep="\t",header = F)
tx_mane_select = gencode_annotation[grepl("MANE",gencode_annotation$V2),c("V1")]

## ENDOF annotation resources

patient_metadata = read.table(patient_table,sep = "\t", header=T,na.strings = c("NA","N/A","","NULL"),stringsAsFactors = T)


supporting_sv_df = read.table(cohort_supporting_sv_path,header=T,sep="\t") 
#bp_name = patient_id + sv_name
supporting_sv_anno_df = supporting_sv_df

#Add pairwise overlap merged
pairwise_overlap_merged_df = read.table(cohort_pairwise_overlap_merged_path,sep = "\t",header = T)

#make consistent with bp_name
pairwise_overlap_merged_df$set1 = paste0(pairwise_overlap_merged_df$patient_id,"_",pairwise_overlap_merged_df$set1)
pairwise_overlap_merged_df$set2 = paste0(pairwise_overlap_merged_df$patient_id,"_",pairwise_overlap_merged_df$set2)

pairwise_overlap_merged_df = pairwise_overlap_merged_df %>% 
  merge(supporting_sv_anno_df[c("bp_name","coordinate","fusion_name")],by.x="set1",by.y="bp_name") %>%
  dplyr::rename(set1_coordinate=coordinate) %>% 
  merge(supporting_sv_anno_df[c("bp_name","coordinate")],by.x="set2",by.y="bp_name") %>%
  dplyr::rename(set2_coordinate=coordinate)

pairwise_overlap_merged_df$overlap_distance = GenomicRanges::distance(GRanges(pairwise_overlap_merged_df$set1_coordinate),GRanges(pairwise_overlap_merged_df$set2_coordinate))

supporting_sv_anno_df = supporting_sv_anno_df %>% 
  left_join(pairwise_overlap_merged_df[c("set1","overlap_distance","overlap_set1_set2")],by=c("bp_name"="set1")) 


## Prepare -> convert to GRanges
supporting_sv_anno = GRanges(supporting_sv_anno_df)
names(supporting_sv_anno) = supporting_sv_anno$bp_name

## Annotate genes
supporting_sv_anno = annotate_svs_genes(supporting_sv_anno,gene_properties)

## Use mane select trancript to annotate introns
ensembl_id_lst = unique(unlist(strsplit(supporting_sv_anno$ensembl_id,", ")))
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
  
supporting_sv_anno = annotate_svs_features(supporting_sv_anno,introns_per_tx) 
mcols(supporting_sv_anno) = as.data.frame(mcols(supporting_sv_anno)) %>% dplyr::rename(intron_overlap=ft_overlap, intron_tx = tx_id)


##Optional: repeats and segmental duplications
repeatmasker_path=paste0(resources_dir,repeatmasker_file)

if(length(Sys.glob(repeatmasker_path))>0){

  repeatmasker = read.table(repeatmasker_path,header=T,sep="\t",comment.char = "")
  
  repeatmasker = repeatmasker %>% dplyr::rename("seqnames"="genoName","start"="genoStart","end"="genoEnd")
  repeatmasker = GRanges(repeatmasker)
  names(repeatmasker) = paste0("rep_",1:length(repeatmasker))
  
  repeats = repeatmasker[repeatmasker$repClass %in% c("LINE","SINE","LTR")]
  repeats=repeats[abs(repeats$repLeft)<50]
  
  ## starting bp 
  svs_start = supporting_sv_anno
  end(svs_start)=start(svs_start)
  svs_start=annotate_svs_properties(svs_start,repeats,c("repFamily"))
  
  ## ending bp 
  svs_end = supporting_sv_anno
  start(svs_end)=end(svs_end)
  svs_end=annotate_svs_properties(svs_end,repeats,c("repFamily"))
  
  supporting_sv_anno$repeat_family_start=mcols(svs_start)$repFamily
  supporting_sv_anno$repeat_family_end=mcols(svs_end)$repFamily


} else {
  print("No repeatmasker file found. Skipping analysis")  
  print(repeatmasker_path)
}

segmental_duplications_path=paste0(resources_dir,segmental_duplications_file)
if(length(Sys.glob(segmental_duplications_path))>0){
  segmental_duplications = read.table(segmental_duplications_path,header=T,sep="\t",comment.char = "")
  segmental_duplications = segmental_duplications %>% dplyr::rename("seqnames"="chrom","start"="chromStart","end"="chromEnd","sd_name"="name")
  segmental_duplications$sd_name=paste0("sd_",segmental_duplications$sd_name)
  segdup = GRanges(segmental_duplications)
  names(segdup) = paste0("sd_",1:length(segdup))
  
  ## starting bp 
  svs_start = supporting_sv_anno
  end(svs_start)=start(svs_start)
  svs_start =annotate_svs_properties(svs_start,segdup,c("sd_name"))
  
  ## ending bp 
  svs_end = supporting_sv_anno
  start(svs_end)=end(svs_end)
  svs_end =annotate_svs_properties(svs_end,segdup,c("sd_name"))
  
  supporting_sv_anno$segdup_start=mcols(svs_start)$sd_name
  supporting_sv_anno$segdup_end=mcols(svs_end)$sd_name

}else {
  print("No segmental duplications file found. Skipping analysis")  
  print(segmental_duplications_path)
}

names(supporting_sv_anno)=rep(1:length(supporting_sv_anno))
supporting_sv_anno_df = as.data.frame(supporting_sv_anno)
supporting_sv_anno_df = unique(supporting_sv_anno_df)

## Optional: copy number
### Add CNA l2fc if exists
svs_cna_mapping_df=data.frame()
for(id in patient_metadata$patient_identifier) {
  patient = dplyr::filter(patient_metadata,patient_identifier==id)
  svs_cna_mapping_path = paste0(cna_reports_dir,svs_cna_outfile,patient$patient_identifier,".tsv")
  if( length(Sys.glob(svs_cna_mapping_path)) == 1) {
    
    svs_cna_mapping = read.table(svs_cna_mapping_path,header=T,sep="\t")
    svs_cna_mapping$patient_id = patient$patient_id
    svs_cna_mapping_df = rbind(svs_cna_mapping_df,svs_cna_mapping)
  }
}
if(nrow(svs_cna_mapping_df)>0){
  svs_cna_mapping_df$bp_name = paste0(svs_cna_mapping_df$patient_id,"_",svs_cna_mapping_df$bp_name)
  supporting_sv_anno_df = supporting_sv_anno_df %>% left_join(svs_cna_mapping_df,by=c("bp_name","sv_name","sv_merged","patient_id"))
} else {
  print("No copy number data found")
}

## Optional database overlaps

databases_lst = c("nstd166","nstd186","dgv")
sv_database_matches = data.frame()

for( sv_database_identifier in databases_lst) {
  database_overlaps_path = paste0(reports_dir,database_overlaps_outfile,sv_database_identifier,".tsv")
  if(length(Sys.glob(database_overlaps_path))==0) { 
    print("Database overlaps not found. Skipping analysis")
    print(database_overlaps_path)
    next() 
  }
  database_overlaps = read.table(database_overlaps_path,header = T,sep = "\t")
  database_overlaps = database_overlaps[,c("set1","set2","overlap_set1_set2","overlap_set2_set1","set1_svtype","set2_svtype")]
  
  overlaps_colname=paste0("anno_sv_db_",sv_database_identifier)
  supporting_sv_anno_df[,overlaps_colname]=F
  
  matching_pop_svs = supporting_sv_anno_df[,c("bp_name","svtype","patient_id")] %>%
    merge(database_overlaps, by.x = "bp_name", by.y= "set2") %>%
    dplyr::rename(db_svtype = set1_svtype)
  
  matching_pop_svs = matching_pop_svs %>% 
    group_by(patient_id,bp_name,svtype,db_svtype) %>% 
    summarize(recip_overlap = mean(c(overlap_set1_set2,overlap_set2_set1)), db_variant=toString(unique(set1)),.groups="keep")
  
  matching_pop_svs$sv_database=sv_database_identifier
  
  sv_database_matches = rbind(sv_database_matches,matching_pop_svs)
  
  #match on type in hindsight but report all
  supporting_sv_anno_df[supporting_sv_anno_df$bp_name %in% dplyr::filter(matching_pop_svs,svtype==db_svtype)$bp_name,overlaps_colname]=T
}

supporting_sv_anno_df = supporting_sv_anno_df %>% mutate(anno_sv_population = (anno_sv_db_nstd166 | anno_sv_db_nstd186 | anno_sv_db_dgv))

write.table(sv_database_matches,population_sv_database_matches_path,quote = FALSE,sep = "\t",row.names=FALSE)

write.table(supporting_sv_anno_df,cohort_supporting_sv_anno_path,quote = F,sep = "\t",row.names=F,col.names = T)


