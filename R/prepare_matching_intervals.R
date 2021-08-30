## Fusion-sq
# Per patient, prepare matching with WGS by deriving the genomic matching intervals
# prerequisite to running run_match_wgs.R
# Transcript table is used also for transcript selection based on SV breakpoints and important also for exon imbalance calculations
## Tx table is made here because the reference database is loaded. 

### Input
## References: .gtf.gz, .sqlite and splice junction database from star fusion 
# Reference versions depend on the star fusion version (gencodev31 by default)
## star fusion: predicted fusions 

### Output
## >> file per patient
## fusion_anno_table (row per fusion, tsv)
## matching_intervals (row per fusion, tsv)
## total_matching_intervals (bed)
## transcript_table (row per transcript, tsv)

# patient_identifier = patient_id + _ + rna_id
# fusion_identifier = patient_identifier + lead number


suppressPackageStartupMessages({
  library(GenomicRanges, quietly=TRUE)
  library(AnnotationDbi, quietly=TRUE)
  library(rtracklayer, quietly=TRUE)
  library(tidyverse, quietly=TRUE)
  library(stringr, quietly=TRUE)
  library(argparser)
})

source("R/default.conf") #in script dir
source(paste0(script_dir,"functions.general.R")) 
source(paste0(script_dir,"functions.prepare_matching.R")) 

patient_metadata = read.table(patient_table,sep = "\t", header=T) 

## One patient mode 

argp = arg_parser("Prepare matching intervals")
argp = add_argument(argp, "--patient_identifier", help="Patient identifier corresponding to patient table")
argv = parse_args(argp)

if(!is.null(argv$patient_identifier) & !is.na(argv$patient_identifier)) {
  patient = dplyr::filter(patient_metadata, patient_identifier==argv$patient_identifier)
} else {
  print("Need patient identifier")
  quit()
}

print("Prepare matching intervals")
print(paste0("Running: patient: ",patient$patient_identifier))


## Reference databases
gtf_path = paste0(resources_dir,"ref_annot_",reference,".gtf.gz")
txdb_path = paste0(resources_dir,"ref_annot_",reference,".sqlite")
sjdb_path = paste0(resources_dir,"sjdb_",reference,".bed")

if(length(Sys.glob(gtf_path))<1) {
  quit("No GTF file found")
}
  
# Create TxDb from gtf for the different versions if not exists 
if(length(Sys.glob(txdb_path))<1) {
  txdb = makeTxDbFromGFF(gtf_path,format = "gtf")
  saveDb(txdb, txdb_path)
} else {
  txdb=loadDb(txdb_path)
}

gtf <- rtracklayer::import(gtf_path)
genes = gtf[gtf$type=="gene"]
genes = data.frame(genes) #annotation 
#genes$gene_id = remove_version_from_id(genes$gene_id) #optional; turned off to prevent version mismatch
gr_genes = GRanges(genes)

exons = exonsBy(txdb, "tx", use.names=TRUE) 
introns = intronsByTranscript(txdb, use.names=TRUE) 
transcripts = transcriptsBy(txdb, by = "gene")
  
sjdb_granges = import(sjdb_path,format="bed") 

transcript_df = gtf[gtf$type=="transcript"]
transcript_df = data.frame(transcript_df)

## END reference loading


##START



  sf_files = Sys.glob(paste0(starfusion_dir,patient$rna_id,"*.tsv"))
  sf_files = sf_files[grep("annotated",sf_files,invert=T)] #should only be the unnannotated now.
  
  if(length(sf_files)<1) { 
    print(paste0("EXIT: ",patient$patient_identifier,": missing star fusion file"))
    quit()
  }
  
  if(length(sf_files)>1) {
    print(paste0("WARNING: multiple star fusion files"))
  }
  starfusion_file = sf_files[1]
  metadata= get_metadata_starfusion(starfusion_file)
  if(metadata[metadata$key=="starfusion_version",c("value")]!="STAR-Fusion/1.8.0"){
    print(paste0("WARNING: ",patient$patient_identifier,": wrong star fusion version (",metadata[metadata$key=="starfusion_version",c("value")],")"))
  }

  # Make annotation table from StarFusion file
  fusion_anno_table_filepath = paste0(base_dir,fusion_annotation_outfile,patient$patient_identifier,".tsv")
  if(length(Sys.glob(fusion_anno_table_filepath))<1) {
    fusion_anno_table = make_fusion_anno_table(starfusion_file)
    ## Adjust fusion anno columns - from STAR fusion names to names used in pipeline
    fusion_anno_table =  rename_fusion_anno_columns(fusion_anno_table)
  
    write.table(fusion_anno_table,fusion_anno_table_filepath,row.names = FALSE, sep="\t",quote=FALSE)
  
  } else {
    fusion_anno_table = read.table(fusion_anno_table_filepath,header = T, sep="\t")
  }
  #endof make fusion_anno_table

  # From fusion anno table: make matching intervals (consensus) and transcript table with intervals and involved fragments
  matching_intervals_filepath = paste0(base_dir,matching_intervals_outfile,patient$patient_identifier,".tsv")
  transcript_table_filepath = paste0(base_dir,transcript_table_outfile,patient$patient_identifier,".tsv")
  total_intervals_filepath = paste0(base_dir,total_matching_intervals_outfile,patient$patient_identifier,".bed")
  
  if(length(Sys.glob(matching_intervals_filepath))>0) {
    print(paste0("EXIT: ",patient$patient_identifier,": already has matching intervals file (",matching_intervals_filepath,")"))
    #quit()
  }
  
  ## Make transcript table and consensus matching intervals
  ## Use these to collect total intervals
  total_intervals = GRanges()
  matching_intervals_table = data.frame(stringsAsFactors=FALSE)
  transcript_table = data.frame(stringsAsFactors=FALSE)
    
  for(i in 1:nrow(fusion_anno_table)){
    fusion = fusion_anno_table[i,]
    
    partner_gene_gup = fusion[,c("gup_ensembl_version","gup_sf_breakpoint","gup_gene_id")]
    colnames(partner_gene_gup)  = c("ensembl_id","sf_breakpoint","gene_name")
    partner_gene_gup$coordinate = as.numeric(strsplit(partner_gene_gup$sf_breakpoint,":")[[1]][2])
    partner_gene_gup$chromosome = strsplit(partner_gene_gup$sf_breakpoint,":")[[1]][1]
    partner_gene_gup$forward_strand = grepl("+",partner_gene_gup$sf_breakpoint,fixed=TRUE)
      
    partner_gene_gdw = fusion[,c("gdw_ensembl_version","gdw_sf_breakpoint","gdw_gene_id")]
    colnames(partner_gene_gdw)  = c("ensembl_id","sf_breakpoint","gene_name")
    partner_gene_gdw$coordinate = as.numeric(strsplit(partner_gene_gdw$sf_breakpoint,":")[[1]][2])
    partner_gene_gdw$chromosome = strsplit(partner_gene_gdw$sf_breakpoint,":")[[1]][1]
    partner_gene_gdw$forward_strand = grepl("+",partner_gene_gdw$sf_breakpoint,fixed=TRUE)
    
    gup_transcript_tbl = make_transcript_table(partner_gene_gup,TRUE)
    gdw_transcript_tbl = make_transcript_table(partner_gene_gdw,FALSE)
    
    fusion_transcript_table = rbind(gup_transcript_tbl,gdw_transcript_tbl,stringsAsFactors = FALSE)
    
    #Add to transcipt table if not empy
    if(nrow(fusion_transcript_table)>0) { 
      fusion_transcript_table$identifier = fusion$identifier
      fusion_transcript_table$fusion_identifier = fusion$fusion_identifier
      transcript_table = rbind(transcript_table,fusion_transcript_table,stringsAsFactors = FALSE)
    }
    
    ## Make matching intervals, also works with empty tx table
    
    row = c()
    row$identifier = fusion$identifier
    row$fusion_identifier = fusion$fusion_identifier
    
    gup_row = make_matching_intervals(partner_gene_gup,TRUE,gup_transcript_tbl)
    gdw_row = make_matching_intervals(partner_gene_gdw,FALSE,gdw_transcript_tbl)
    
    #tag overlapping genes so you can throw a warning also later in the pipeline
    row$overlap_gup_gdw_genebody = ifelse(length(GenomicRanges::intersect(gup_row$gene_coordinates, gdw_row$gene_coordinates))>0,T,F)
    row$overlap_gup_gdw_adjacent_intron = ifelse(length(GenomicRanges::intersect(gup_row$adjacent_intron, gdw_row$adjacent_intron))>0,T,F)
    
    # warn for overlapping adjacent introns
    if(row$overlap_gup_gdw_adjacent_intron) {
      print(paste0("Warning: ",fusion$identifier," - ",fusion$fusion_name, " has overlapping upstream downstream adjacent introns "))
    } else if(row$overlap_gup_gdw_genebody) {
      print(paste0("Warning: ",fusion$identifier," - ",fusion$fusion_name, " has overlapping upstream downstream genes "))
    }
    
    
    #expected that everything falls within gene body and/or flanking interval if there is no gene body
    total_intervals = c(total_intervals, 
                        gup_row$gene_coordinates, gdw_row$gene_coordinates,
                        gup_row$flanking, gdw_row$flanking)
    
    # returns vector with genomic ranges => to string lapply(row, function(x) { toString(x) } )
    gup_row = lapply(gup_row, function(x) { toString(x) }) %>% as.data.frame() %>% rename_with( ~ paste0("gup_",.))
    gdw_row = lapply(gdw_row, function(x) { toString(x) }) %>% as.data.frame() %>% rename_with( ~ paste0("gdw_",.))
    
    row = c(row,gup_row,gdw_row)
    ## convert to strings for saving  
    matching_intervals_table = rbind(matching_intervals_table,row, stringsAsFactors = FALSE)
  
  }

  matching_intervals_table$identifier = as.integer(matching_intervals_table$identifier)
    
  #parse total intervals: reduce overlaps, make larger, remove mitochondrial
  total_intervals = GenomicRanges::reduce(total_intervals)
  total_intervals = GenomicRanges::resize(total_intervals,width(total_intervals)+1000,fix="center")
  total_intervals = total_intervals[!seqnames(total_intervals)=="chrM"]

  write.table(matching_intervals_table,matching_intervals_filepath,row.names = FALSE, sep="\t",quote=FALSE)
  write.table(transcript_table,transcript_table_filepath,row.names = FALSE, sep="\t",quote=FALSE)
  write.table(total_intervals,total_intervals_filepath,row.names = FALSE, sep="\t",quote=FALSE)
  
  ##endof making matching intervals, transcript tables
  