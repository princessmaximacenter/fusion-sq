library(chimeraviz)
library(GenomicRanges)
library(VariantAnnotation, quietly=TRUE)
library(StructuralVariantAnnotation, quietly=TRUE)
library(rtracklayer, quietly=TRUE)
library(tidyverse, quietly=TRUE)
library(stringr, quietly=TRUE)
library(testthat, quietly=TRUE)
library(stringdist, quietly=TRUE)
library(BSgenome.Hsapiens.UCSC.hg38, quietly=TRUE)

wdir="~/PycharmProjects/wdl_pipeline/fusion_pilot/"
output_dir = paste0(wdir,"fusion_sq_grape/")
source(paste0(output_dir,"default.conf"))

base_dir = paste0(wdir,"fusion_sq_grape/base_20201109/")
source(paste0(wdir,"fusion_sq_grape/fusion_functions.grape.R"))

resources_dir="~/Documents/resources/"
starfusion_dir = paste0(data_dir,"starfusion_18/")
reference = "GRCh38_gencode_v31_CTAT_lib_Oct012019"

options(ucscChromosomeNames=FALSE)

#reference = "GRCh38_gencode_v31_CTAT_lib_Oct012019"

gtf_path = paste0(resources_dir,"Homo_sapiens.GRCh38.78.gtf.gz")
edbSqliteFile = paste0(resources_dir,"Homo_sapiens.GRCh38.97.sqlite")
#edbSqlite <- ensDbFromGtf(gtf = gtf_path,outfile="/Users/ianthevanbelzen/Documents/resources/Homo_sapiens.GRCh38.97.sqlite")
edb <- ensembldb::EnsDb(edbSqliteFile)

starfusion_dir = paste0(data_dir,"starfusion_20210116/")

patient_table = paste0(wdir,"cohort.v3.metadata.20201231.tsv")
patient_metadata = read.table(patient_table,sep = "\t", header=T,stringsAsFactors = T)

reports_dir="~/PycharmProjects/wdl_pipeline/fusion_pilot/fusion_sq_grape/reports_20210108/"
cohort_report =  read.table(paste0(reports_dir,cohort_report_outfile,".tsv"),header=T,sep="\t")
cohort_report$gup_sf_tx = remove_version_from_id(cohort_report$gup_sf_tx)
cohort_report$gdw_sf_tx = remove_version_from_id(cohort_report$gdw_sf_tx)


#for coverage plots
options(ucscChromosomeNames=FALSE) 


## functions
get_transcript_obj = function(transcripts,transcript_id) {
  #metadata columns are added to the transcript GRanges object
  #transcript ID as $transcript and $symbol
  #$feature utr5/prot.coding/utr3
  #$transcript_category exonBoundary or not
  #$coding
  #make into GRanges list and set fusion@gene_upstream@transcripts <- grangeslist_upstream
  
  tra=transcripts[[transcript_id]]
  mcols(tra)$transcript = transcript_id
  mcols(tra)$symbol = transcript_id
  tra= split_on_utr_and_add_feature(tra)
  tra_lst = GRangesList(tra)
  mcols(tra_lst)$transcript_category = "exonBoundary" #
  #mcols(tra_lst)$transcript_category = decide_transcript_category(tra,target)
  mcols(tra_lst)$coding = TRUE
  #mcols(tra_lst)$coding = !all(is.na(mcols(tra)$tx_cds_seq_start))
  mcols(tra_lst)$transcript = transcript_id
  names(tra_lst)=transcript_id
  return(tra_lst)
}


select_tx_fusion = function(edb,selected_fusion) {
  transcripts = exonsBy(edb,filter = list(
    AnnotationFilter::GeneIdFilter(
      c(
        selected_fusion@gene_upstream@ensembl_id,
        selected_fusion@gene_downstream@ensembl_id))),
    columns = c(
      "gene_id",
      "gene_name",
      "tx_id",
      "tx_cds_seq_start",
      "tx_cds_seq_end",
      "exon_id"))
  return(transcripts)
}


##Notes
#delete, useless  #library(PFAM.db)

## For protein domains 
#BiocManager::install("AnnotationHub")
#BiocManager::install("biomaRt")

library(AnnotationHub)
hub = AnnotationHub()
ahDb = query(hub, c("EnsDb", "apiens", "97"))
ahEdb = ahDb[[1]]
#listTables(ahEdb)

library(biomaRt)
ensembl = useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",version=97)

#listAttributes(ensembl) %>% dplyr::filter(grepl("interpro",name))
#listAttributes(ensembl) %>% dplyr::filter(grepl("transcript",name))
#listAttributes(ensembl)

# getBM() function is the main query function in biomaRt. It has four main arguments:
#   
# attributes: is a vector of attributes that one wants to retrieve (= the output of the query).
# filters: is a vector of filters that one wil use as input to the query.
# values: a vector of values for the filters. In case multple filters are in use, the values argument requires a list of values where each position in the list corresponds to the position of the filters in the filters argument (see examples below).
# mart: is an object of class Mart, which is created by the useMart() function.

### endof notes

#### START PLOTS ####

## PMCID308AAA_RB1--DGKB

selected_patient = "PMCID308AAA"
patient = patient_metadata %>% filter(patient_id==selected_patient)

sf_files = Sys.glob(paste0(starfusion_dir,"sf18/",patient$rna_id,"*.tsv"))
sf_files = sf_files[grep("annotated",sf_files,invert=T)] #should only be the unnannotated now.
sf_files[1] #= "/Users/ianthevanbelzen/data/starfusion_20210116/PMABM000BGM_PMCRZ962VUW_RNA-Seq.star-fusion_predicted.tsv"

fusions = import_starfusion(sf_files[1],'hg38') 

selected_fusion_wgs = cohort_report %>% filter(patient_fusion=="PMCID308AAA_RB1--DGKB")
selected_fusion = get_fusion_by_id(fusions,1)

transcripts = select_tx_fusion(edb,selected_fusion)


transcript_table =  read.table(paste0(base_dir,transcript_table_outfile,patient$patient_identifier,".tsv"),header=T,sep="\t") %>% rename(fusion_id = identifier)

transcript_table %>% filter(fusion_id==selected_fusion_wgs$fusion_predictions & grepl(selected_fusion_wgs$gup_sf_tx,transcript_id))
transcript_table %>% filter(fusion_id==selected_fusion_wgs$fusion_predictions & grepl(selected_fusion_wgs$tx_gup_transcript_id,transcript_id))

#two different ones try both
gup_tx = selected_fusion_wgs$gup_sf_tx 
gup_tx = selected_fusion_wgs$tx_gup_transcript_id #seems to be default transcipt as well 

#overlaps for gdw
#selected_fusion_wgs$gdw_sf_tx
#selected_fusion_wgs$tx_gdw_transcript_id
gdw_tx = selected_fusion_wgs$gdw_sf_tx

selected_fusion@gene_upstream@transcripts = get_transcript_obj(transcripts,gup_tx)
selected_fusion@gene_downstream@transcripts = get_transcript_obj(transcripts,gdw_tx)

plot_transcripts(selected_fusion,edb=edb,reduce_transcripts = FALSE) 

plot_fusion_separate(selected_fusion,
                     edb=edb,ylim = c(0,250),
                     bamfile = "~/data/cram_partials/PMABM000BGM_PMCRZ600VOJ_RNA-Seq_RB1_DGKB.bam",non_ucsc = FALSE,
                     reduce_transcripts = T)

plot_fusion_separate(selected_fusion,
                     edb=edb,ylim = c(0,200),
                     bamfile = "~/data/cram_partials/PMABM000BGM_PMCRZ600VOJ_RNA-Seq_RB1_DGKB.fwd.bam",non_ucsc = FALSE,
                     reduce_transcripts = T)


plot_fusion_separate(selected_fusion,
                     edb=edb,ylim = c(0,200),
                     bamfile = "~/data/cram_partials/PMABM000BGM_PMCRZ600VOJ_RNA-Seq_RB1_DGKB.rev.bam",non_ucsc = FALSE,
                     reduce_transcripts = T)


plot_transcripts(selected_fusion,
                     edb=edb,ylim = c(0,200),
                     bamfile = "~/data/cram_partials/PMABM000BGM_PMCRZ600VOJ_RNA-Seq_RB1_DGKB.fwd.bam",non_ucsc = FALSE,
                     reduce_transcripts = T)


plot_transcripts(selected_fusion,
                     edb=edb,ylim = c(0,30),
                     bamfile = "~/data/cram_partials/PMABM000BGM_PMCRZ600VOJ_RNA-Seq_RB1_DGKB.rev.bam",non_ucsc = FALSE,
                     reduce_transcripts = T)

#Domains still broken
if(FALSE){
  
  protein_domains = proteins(ahEdb, filter = ~tx_id == gup_tx ,
                             columns = c("interpro_accession","protein_domain_source","protein_domain_id", "prot_dom_start",
                                         "prot_dom_end","protein_domain_source"))
  protein_domains
interpro_domains = getBM(mart=ensembl, attributes =c("ensembl_transcript_id",
                                  "interpro","interpro_short_description","interpro_description","interpro_start","interpro_end"),
      filters="ensembl_transcript_id",
      values=c(gup_tx,gdw_tx)) 
interpro_domains
colnames(interpro_domains) = c("Transcript_id","Pfam_id","Domain_name_abbreviation","Domain_name_full","Start","End")
#reorder for bed file
interpro_domains = interpro_domains[,c("Transcript_id","Pfam_id","Start","End","Domain_name_abbreviation","Domain_name_full")]
write.table(interpro_domains,paste0(output_dir,"interpro_domains_RB1.bed"),col.names = T,row.names = F,quote = T,sep="\t")

#bedfile with 
#Transcript_id	Pfam_id	Start	End	Domain_name_abbreviation	Domain_name_full

plot_fusion_transcript_with_protein_domain(selected_fusion,
                                           edb=edb,bedfile=paste0(output_dir,"interpro_domains_RB1.bed"),
                                           gene_upstream_transcript = gup_tx,
                                           gene_downstream_transcript = gdw_tx)

plot_fusion(selected_fusion,
            edb=edb,bedfile=paste0(output_dir,"interpro_domains_RB1.bed"))
}

## merge overlapping domains? 
#Idea; also annotate the breakpoints as domains?
#selected_fusion_wgs$gup_coordinate
#selected_fusion_wgs$gdw_coordinate

### ENDOF Patient 308AAA


### Patient HOX
selected_patient = "PMCID453AAA"
patient = patient_metadata %>% filter(patient_id==selected_patient)

sf_files = Sys.glob(paste0(starfusion_dir,"sf18/",patient$rna_id,"*.tsv"))
sf_files = sf_files[grep("annotated",sf_files,invert=T)] #should only be the unnannotated now.
if(length(sf_files)>1) {
  print(sf_files)
}

fusions = import_starfusion(sf_files[1],'hg38') 

selected_fusion_wgs = cohort_report %>% filter(patient_fusion=="PMCID453AAA_MED14--HOXA9")
selected_fusion_wgs$fusion_predictions

selected_fusion = get_fusion_by_id(fusions,selected_fusion_wgs$fusion_predictions)
plot_transcripts(selected_fusion,edb=edb,reduce_transcripts = FALSE,bamfile = "~/data/cram_partials/PMABM000CNL_PMCRZ036GTD_RNA-Seq_MED14_HOXA9.bam",non_ucsc = FALSE)

transcripts = select_tx_fusion(edb,selected_fusion)

#upstream tx is allright
selected_fusion_wgs$gup_sf_tx == selected_fusion_wgs$tx_gup_transcript_id

#downstream is different
selected_fusion_wgs$gdw_sf_tx
selected_fusion_wgs$tx_gdw_transcript_id

transcript_table =  read.table(paste0(base_dir,transcript_table_outfile,patient$patient_identifier,".tsv"),header=T,sep="\t") %>% rename(fusion_id = identifier)

transcript_table %>% filter(fusion_id==selected_fusion_wgs$fusion_predictions & upstream==FALSE)

transcript_table %>% filter(fusion_id==selected_fusion_wgs$fusion_predictions & grepl(selected_fusion_wgs$gdw_sf_tx,transcript_id))
transcript_table %>% filter(fusion_id==selected_fusion_wgs$fusion_predictions & grepl("ENST00000487384",transcript_id))
transcript_table %>% filter(fusion_id==selected_fusion_wgs$fusion_predictions & grepl("ENST00000489695",transcript_id))


gup_tx = selected_fusion_wgs$gup_sf_tx 

#nomal protein but exon disrupted which gives issues with plot
gdw_tx = selected_fusion_wgs$gdw_sf_tx
gdw_tx="ENST00000396345"

selected_fusion@gene_upstream@transcripts = get_transcript_obj(transcripts,gup_tx)
selected_fusion@gene_downstream@transcripts = get_transcript_obj(transcripts,gdw_tx)

plot_transcripts(selected_fusion,edb=edb,reduce_transcripts = FALSE) 
plot_transcripts(selected_fusion,edb=edb,reduce_transcripts = FALSE,
                 bamfile = "~/data/cram_partials/PMABM000CNL_PMCRZ036GTD_RNA-Seq_MED14_HOXA9.bam",non_ucsc = FALSE,
                 ylim = c(0,300))

plot_transcripts(selected_fusion,edb=edb,reduce_transcripts = TRUE,
                 bamfile = "~/data/cram_partials/PMABM000CNL_PMCRZ036GTD_RNA-Seq_MED14_HOXA9.bam",non_ucsc = FALSE,
                 ylim = c(0,300))


plot_fusion_separate(selected_fusion,
                     edb=edb,ylim = c(0,200),
                     bamfile = "~/data/cram_partials/PMABM000CNL_PMCRZ036GTD_RNA-Seq_MED14_HOXA9.bam",non_ucsc = FALSE,
                     reduce_transcripts = T)

interpro_domains = getBM(mart=ensembl, attributes =c("ensembl_transcript_id",
                                                     "interpro","interpro_short_description","interpro_description","interpro_start","interpro_end"),
                         filters="ensembl_transcript_id",
                         values=c(gup_tx,gdw_tx)) 
interpro_domains

#replace by identifie alt protein
gdw_tx="ENST00000396345"
interpro_domains[interpro_domains$ensembl_transcript_id=="ENST00000343483",c("ensembl_transcript_id")]=gdw_tx


colnames(interpro_domains) = c("Transcript_id","Pfam_id","Domain_name_abbreviation","Domain_name_full","Start","End")
#reorder for bed file
interpro_domains = interpro_domains[,c("Transcript_id","Pfam_id","Start","End","Domain_name_abbreviation","Domain_name_full")]
write.table(interpro_domains,paste0(output_dir,"interpro_domains_MED14_HOX.bed"),col.names = T,row.names = F,quote = T,sep="\t")

selected_fusion@gene_downstream@transcripts = get_transcript_obj(transcripts,gdw_tx)

## doesnt work
plot_fusion_transcript_with_protein_domain(selected_fusion,
                                           edb=edb,bedfile=paste0(output_dir,"interpro_domains_MED14_HOX.bed"),
                                           gene_upstream_transcript = gup_tx,
                                           gene_downstream_transcript = gdw_tx,plot_downstream_protein_domains_if_fusion_is_out_of_frame = TRUE) 

## note does work if you use unsupported tx gdw_tx="ENST00000396345"
 
## so it is a bit fake but seems like homeobox is included 
## TODO expression

### ENDOF patient


### Patient ERBB
selected_patient = "PMCID912AAJ"
patient = patient_metadata %>% filter(patient_id==selected_patient)

sf_files = Sys.glob(paste0(starfusion_dir,"sf18/",patient$rna_id,"*.tsv"))
sf_files = sf_files[grep("annotated",sf_files,invert=T)] #should only be the unnannotated now.
if(length(sf_files)>1) {
  print(sf_files)
}

fusions = import_starfusion(sf_files[1],'hg38') 

selected_fusion_wgs = cohort_report %>% filter(patient_fusion=="PMCID912AAJ_ERBB4--LINC01807")
selected_fusion_wgs$fusion_predictions

selected_fusion = get_fusion_by_id(fusions,selected_fusion_wgs$fusion_predictions)

transcripts = select_tx_fusion(edb,selected_fusion)

## no tx selected
selected_fusion_wgs$gup_sf_tx 
selected_fusion_wgs$tx_gup_transcript_id
selected_fusion_wgs$gdw_sf_tx
selected_fusion_wgs$tx_gdw_transcript_id


transcript_table =  read.table(paste0(base_dir,transcript_table_outfile,patient$patient_identifier,".tsv"),header=T,sep="\t") %>% rename(fusion_id = identifier)

transcript_table %>% filter(fusion_id==selected_fusion_wgs$fusion_predictions & upstream)
transcript_table %>% filter(fusion_id==selected_fusion_wgs$fusion_predictions & upstream==FALSE)

gup_tx="ENST00000402597"
gdw_tx="ENST00000665053"

selected_fusion = get_transcripts_ensembl_db(selected_fusion, edb)

plot_transcripts(selected_fusion,edb=edb)

selected_fusion@gene_upstream@transcripts = get_transcript_obj(transcripts,gup_tx)
selected_fusion@gene_downstream@transcripts <- select_transcript(downstream_partner_gene(selected_fusion),
                                                                 which_transcripts =gdw_tx)
                                                                 #which_transcripts = "ENST00000432481")

plot_transcripts(selected_fusion,edb=edb,reduce_transcripts = T,ylim = c(0,120),
            bamfile = "~/data/cram_partials/PMABM000DCY_PMCRZ315XDV_RNA-Seq_ERBB4.rev.smaller.subsample_10.bam",non_ucsc = FALSE)

plot_transcripts(selected_fusion,edb=edb,reduce_transcripts = T,ylim = c(0,60),
                 bamfile = "~/data/cram_partials/PMABM000DCY_PMCRZ315XDV_RNA-Seq_ERBB4.rev.smaller.subsample_10.bam",non_ucsc = FALSE)


interpro_domains = getBM(mart=ensembl, attributes =c("ensembl_transcript_id",
                                                     "interpro","interpro_short_description","interpro_description","interpro_start","interpro_end"),
                         filters="ensembl_transcript_id",
                         values=c(gup_tx,gdw_tx)) 
interpro_domains
colnames(interpro_domains) = c("Transcript_id","Pfam_id","Domain_name_abbreviation","Domain_name_full","Start","End")
#reorder for bed file
interpro_domains = interpro_domains[,c("Transcript_id","Pfam_id","Start","End","Domain_name_abbreviation","Domain_name_full")]
write.table(interpro_domains,paste0(output_dir,"interpro_domains_ERBB4.bed"),col.names = T,row.names = F,quote = T,sep="\t")

plot_fusion_transcript_with_protein_domain(selected_fusion,
                                           edb=edb,bedfile=paste0(output_dir,"interpro_domains_ERBB4.bed"),
                                           gene_upstream_transcript = gup_tx,
                                           gene_downstream_transcript = gdw_tx)


### end of patienet



### patient qki-map


selected_patient = "PMCID203AAL"
patient = patient_metadata %>% filter(patient_id==selected_patient)

sf_files = Sys.glob(paste0(starfusion_dir,"sf18/",patient$rna_id,"*.tsv"))
sf_files = sf_files[grep("annotated",sf_files,invert=T)] #should only be the unnannotated now.
if(length(sf_files)>1) {
  print(sf_files)
}

fusions = import_starfusion(sf_files[1],'hg38') 

selected_fusion_wgs = cohort_report %>% filter(patient_fusion=="PMCID203AAL_QKI--MAP3K4")
selected_fusion_wgs$fusion_predictions

selected_fusion = get_fusion_by_id(fusions,selected_fusion_wgs$fusion_predictions)

transcripts = select_tx_fusion(edb,selected_fusion)

selected_fusion_wgs$gup_sf_tx 
strsplit(selected_fusion_wgs$tx_gup_transcript_id,", ")[[1]]
selected_fusion_wgs$gdw_sf_tx
selected_fusion_wgs$tx_gdw_transcript_id


transcript_table =  read.table(paste0(base_dir,transcript_table_outfile,patient$patient_identifier,".tsv"),header=T,sep="\t") %>% rename(fusion_id = identifier)
transcript_table$transcript_id=remove_version_from_id(transcript_table$transcript_id)
transcript_table %>% filter(fusion_id==selected_fusion_wgs$fusion_predictions & upstream & transcript_id %in% strsplit(selected_fusion_wgs$tx_gup_transcript_id,", ")[[1]]
)
transcript_table %>% filter(fusion_id==selected_fusion_wgs$fusion_predictions & upstream==FALSE & transcript_id %in% strsplit(selected_fusion_wgs$tx_gdw_transcript_id,", ")[[1]]
)

gup_tx= "ENST00000361752"

gdw_tx=selected_fusion_wgs$gdw_sf_tx

selected_fusion@gene_upstream@transcripts = get_transcript_obj(transcripts,gup_tx)
selected_fusion@gene_downstream@transcripts = get_transcript_obj(transcripts,gdw_tx)

plot_transcripts(selected_fusion,edb=edb,reduce_transcripts = TRUE,ylim = c(0,800),
                 bamfile = "~/data/cram_partials/PMABM000HAT_PMCRZ290IRG_RNA-Seq_QKI_MAP3K4.fwd.bam",non_ucsc = FALSE)


interpro_domains = getBM(mart=ensembl, attributes =c("ensembl_transcript_id",
                                                     "interpro","interpro_short_description","interpro_description","interpro_start","interpro_end"),
                         filters="ensembl_transcript_id",
                         values=c(strsplit(selected_fusion_wgs$tx_gup_transcript_id,", ")[[1]],strsplit(selected_fusion_wgs$tx_gdw_transcript_id,", ")[[1]])) 
interpro_domains

colnames(interpro_domains) = c("Transcript_id","Pfam_id","Domain_name_abbreviation","Domain_name_full","Start","End")
#reorder for bed file
interpro_domains = interpro_domains[,c("Transcript_id","Pfam_id","Start","End","Domain_name_abbreviation","Domain_name_full")]

write.table(interpro_domains,paste0(output_dir,"interpro_domains_QKI_MAP3K4.bed"),col.names = T,row.names = F,quote = T,sep="\t")

selected_fusion = get_fusion_by_id(fusions,selected_fusion_wgs$fusion_predictions)

plot_fusion_transcript_with_protein_domain(selected_fusion,edb=edb,
                                           gene_upstream_transcript = gup_tx,
                                           gene_downstream_transcript = "ENST00000366919",
                                           bedfile = paste0(output_dir,"interpro_domains_QKI_MAP3K4.bed"),plot_downstream_protein_domains_if_fusion_is_out_of_frame = T)

transcript_table %>% filter(fusion_id==selected_fusion_wgs$fusion_predictions & upstream==FALSE & transcript_id %in% strsplit(selected_fusion_wgs$tx_gdw_transcript_id,", ")[[1]]
)


### endof patient



## circos plots for 288AAK



selected_patient = "PMCID288AAK"
patient = patient_metadata %>% filter(patient_id==selected_patient)

sf_files = Sys.glob(paste0(starfusion_dir,"sf18/",patient$rna_id,"*.tsv"))
sf_files = sf_files[grep("annotated",sf_files,invert=T)] #should only be the unnannotated now.
if(length(sf_files)>1) {
  print(sf_files)
}

fusions = import_starfusion(sf_files[1],'hg38') 

selected_fusion_wgs = cohort_report %>% filter(patient_id==selected_patient & somatic_variant_low_af & grepl(",",tool)) 
#fusion_id_lst =as.numeric(unlist(strsplit(selected_fusion_wgs$fusion_predictions,", ")))
fusion_id_lst = as.numeric(selected_fusion_wgs$fusion_predictions[!grepl(",",selected_fusion_wgs$fusion_predictions)])
selected_fusion_wgs$fusion_predictions[grepl(",",selected_fusion_wgs$fusion_predictions)]
fusion_id_lst = c(fusion_id_lst,43,18,7)
plot_circle(fusions[fusion_id_lst])

### for patient with tp53 and many 

selected_patient = "PMCID389AAA"
patient = patient_metadata %>% filter(patient_id==selected_patient)

sf_files = Sys.glob(paste0(starfusion_dir,"sf18/",patient$rna_id,"*.tsv"))
sf_files = sf_files[grep("annotated",sf_files,invert=T)] #should only be the unnannotated now.
if(length(sf_files)>1) {
  print(sf_files)
}

fusions = import_starfusion(sf_files[1],'hg38') 

selected_fusion_wgs = cohort_report %>% filter(patient_id==selected_patient & somatic_variant_low_af & grepl(",",tool))
fusion_id_lst =as.numeric(unlist(strsplit(selected_fusion_wgs$fusion_predictions,", ")))

plot_circle(fusions[fusion_id_lst])

selected_fusion_wgs = cohort_report %>% filter(patient_id==selected_patient & somatic_variant & grepl(",",tool))
fusion_id_lst =as.numeric(unlist(strsplit(selected_fusion_wgs$fusion_predictions,", ")))

plot_circle(fusions[fusion_id_lst])

##

