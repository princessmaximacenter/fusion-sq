## Fusion sq
## Pairwise overlap with Population SV database
## Last update: 2021-03-15
## 
## Small refactoring 2021-08-18
## External databases
## Use supporting svs output of collect cohort
suppressPackageStartupMessages({
  
library(tidyverse, quietly=TRUE)
library(stringr, quietly=TRUE)
library(stringdist, quietly=TRUE)
library(GenomicRanges, quietly=TRUE)
library(VariantAnnotation, quietly=TRUE)
library(StructuralVariantAnnotation, quietly=TRUE)
library(rtracklayer, quietly=TRUE)
library(dplyr)
})

wdir="~/fusion_sq/"
script_dir = paste0(wdir,"R/")
source(paste0(script_dir,"default.conf"))
source(paste0(script_dir,"functions.general.R")) 
source(paste0(script_dir,"functions.svs.R"))

utils_dir=paste0(wdir,"utils/")
source(paste0(utils_dir,"functions.population_svs.R"))

### TEST DATA
source("/Users/ianthevanbelzen/fusion_sq/test_data/test_data_override.conf")


## Input

patient_metadata = read.table(patient_table,sep = "\t", header=T,na.strings = c("NA","N/A","","NULL"),stringsAsFactors = T)

## Output
cohort_supporting_svs_path = paste0(reports_dir,cohort_supporting_svs_outfile)

supporting_sv_ranges_cols = c("bp_name","from_coordinate","svLen","insLen","tumor_af","normal_af",
                              "tool", "patient_id", "tumor_normal_ratio", "somatic_variant",
                              "germline_variant", "low_af", "fusion_predictions","fusion_name")


if(length(Sys.glob(cohort_supporting_svs_path))==1) {
  supporting_svs_df = read.table(cohort_supporting_svs_path,header=T,sep="\t") 
  
} else {
  print(paste0("File missing: ",cohort_supporting_svs_path))
  quit() 
}


# Prepare supporting svs

supporting_svs_df$from_coordinate = paste0(supporting_svs_df$seqnames,":",supporting_svs_df$start,"-",supporting_svs_df$end)
supporting_bp = GRanges(supporting_svs_df)
names(supporting_bp)=supporting_bp$bp_name
seqlevelsStyle(supporting_bp)="Ensembl"

## For the CTX and other single breakpoints: resize if <30 
supporting_bp[width(supporting_bp)<30] = flank(supporting_bp[width(supporting_bp)<30],width = 15,both=T)

#####
## Per datdabase find and save overlaps 

## NCBI multiple databases with standarized formats
databases_lst = c("nstd166","nstd186")

#nstd186 (NCBI Curated Common Structural Variants)
#nstd166 (gnomAD Structural Variants)

sv_database_df_cols = c("bp_name","to_coordinate","sv_db_svlen","sv_db_af")

for( sv_database_identifier in databases_lst) {
  sv_database_path = paste0(resources_dir,sv_database_identifier,".GRCh38.variant_call.vcf.gz")
  sv_database_gr = load_sv_database_vcf(sv_database_path)
  
  ##prepare sv database
  mcols(sv_database_gr)$svtype = sv_database_gr$SVTYPE
  sv_database_ranges = make_range_sv_database(sv_database_gr) 
  
  ## cant do SV type matching during overlaps because not exactly the same
  database_overlaps = get_reciprocal_overlap_pairs(sv_database_ranges,supporting_bp,reciprocal_overlap = 0.5,svtype_matching = FALSE)

  database_overlaps$sv_database=sv_database_identifier
  database_overlaps_path = paste0(reports_dir,database_overlaps_outfile,sv_database_identifier,".tsv")
  write.table(database_overlaps,database_overlaps_path,quote = F,sep = "\t",row.names=F,col.names = T)
  
  ## some annotation to make interpretation easier
  
  sv_database_df = as.data.frame(sv_database_ranges[unique(database_overlaps$set1)])
  sv_database_df$to_coordinate = paste0("chr",sv_database_df$seqnames,":",sv_database_df$start,"-",sv_database_df$end)
  sv_database_df = sv_database_df %>% dplyr::rename(sv_db_svlen = SVLEN, sv_db_af = AF )
  
  database_overlaps = database_overlaps %>%  
    left_join(sv_database_df[,sv_database_df_cols],by=c("set1"="bp_name")) %>%
    left_join(supporting_svs_df[,supporting_sv_ranges_cols],by=c("set2"="bp_name")) %>% unique() 
  
  database_overlaps_path = paste0(reports_dir,database_overlaps_outfile,sv_database_identifier,".anno.tsv")
  write.table(database_overlaps,database_overlaps_path,quote = F,sep = "\t",row.names=F,col.names = T)

}


## DGV database

## Another type of datbase
sv_database_identifier ="dgv"
dgv_db_path = paste0(resources_dir,"dgv_GRCh38_hg38_variants_2020-02-25.tsv")
dgv_db = read.table(dgv_db_path,header=T,sep="\t")
dgv_db = dgv_db %>% dplyr::rename(bp_name = variantaccession, seqnames=chr,svtype=varianttype)

dgv_db[grepl("duplication",dgv_db$variantsubtype),c("svtype")]="DUP"
dgv_db[dgv_db$variantsubtype %in% c("gain"),c("svtype")]="DUP"

dgv_db[grepl("deletion",dgv_db$variantsubtype),c("svtype")]="DEL"
dgv_db[dgv_db$variantsubtype %in% c("loss"),c("svtype")]="DEL"

dgv_db[grepl("insertion",dgv_db$variantsubtype),c("svtype")]="INS"
dgv_db[grepl("inversion",dgv_db$variantsubtype),c("svtype")]="INV"

dgv_db[dgv_db$variantsubtype %in% c("gain+loss","complex"),c("svtype")]="other"

dgv_gr = GRanges(dgv_db[dgv_db$seqnames != "",])
names(dgv_gr) = dgv_gr$bp_name

database_overlaps = get_reciprocal_overlap_pairs(dgv_gr,supporting_bp,reciprocal_overlap = 0.5,svtype_matching = FALSE)
database_overlaps$sv_database=sv_database_identifier
database_overlaps_path = paste0(reports_dir,database_overlaps_outfile,sv_database_identifier,".tsv")
write.table(database_overlaps,database_overlaps_path,quote = F,sep = "\t",row.names=F,col.names = T)


dgv_db = dgv_db[dgv_db$bp_name %in% unique(database_overlaps$set1),]
dgv_db$to_coordinate = paste0("chr",dgv_db$seqnames,":",dgv_db$start,"-",dgv_db$end)
dgv_db = dgv_db %>% dplyr::rename(sv_db_reference = reference )

dgv_db_cols = c("bp_name","to_coordinate","sv_db_reference")

database_overlaps = database_overlaps %>%  
  left_join(dgv_db[,dgv_db_cols],by=c("set1"="bp_name")) %>%
  left_join(supporting_svs_df[,supporting_sv_ranges_cols],by=c("set2"="bp_name")) %>% unique() 

database_overlaps_path = paste0(reports_dir,database_overlaps_outfile,sv_database_identifier,".anno.tsv")

write.table(database_overlaps,database_overlaps_path,quote = F,sep = "\t",row.names=F,col.names = T)


