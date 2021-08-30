## Fusion Pilot
## Grape

# Run match WGS
# per patient, per tool
# Configuration:
# Tool: manta, delly, gridss

### Output
## >> file per patient, per tool 
## linking table (row per fusion, tsv)
## supporting_breakpoints (bed)

# patient_identifier = patient_id + _ + rna_id
# fusion_identifier = patient_identifier + lead number


source("R/default.conf") #in script dir
source(paste0(script_dir,"functions.general.R"))
source(paste0(script_dir,"functions.match_wgs.R")) 
source(paste0(script_dir,"functions.get_vcf_filepath.R")) ## adjust if needed


library(argparser, quietly=TRUE)

argp = arg_parser("Match SV breakpoints to fusions for a single patient")
argp = add_argument(argp, "--patient_identifier", help="Patient identifier corresponding to patient table")
argp = add_argument(argp, "--tool", help="SV tool (manta/delly/gridss)")

argp = add_argument(argp, "--no_composite", flag = TRUE, help="Turn off search for composite fusions - requires more resources")

argp = add_argument(argp, "--config", default="", help="Analysis configuration file -optional")

## defaults are filled based on patient ID if not provided
argp = add_argument(argp, "--vcf", default="", help="VCF output file from Delly/GRIDSS (unfiltered) or the Manta diploid")
argp = add_argument(argp, "--vcf_somatic", default="", help="VCF somatic output file, Manta only")

argp = add_argument(argp, "--matching_intervals", default="", help="Genomic matching intervals file from base")
argp = add_argument(argp, "--total_intervals", default="", help="Total genomic intervals file from base")                    
argp = add_argument(argp, "--fusion_anno_table", default="", help="Fusion annotation table from base")                    

argv = parse_args(argp)

early_exit=TRUE
run_composite=!(argv$no_composite)

if(file.exists(argv$config)) {
  source(argv$config)
}

##
suppressPackageStartupMessages({
  library(GenomicRanges, quietly=TRUE)
  library(AnnotationDbi, quietly=TRUE)
  library(VariantAnnotation, quietly=TRUE)
  library(StructuralVariantAnnotation, quietly=TRUE)
  library(rtracklayer, quietly=TRUE)
  library(tidyverse, quietly=TRUE)
  library(stringr, quietly=TRUE)
})
### 

patient_metadata = read.table(patient_table,sep = "\t", header=T)
patient_metadata$tumor_label = paste0(patient_metadata$tumor_id,"_WGS")
patient_metadata$normal_label = paste0(patient_metadata$normal_id,"_WGS")

## subset the patient_metadata
if(!is.null(argv$patient_identifier) & !is.na(argv$patient_identifier)) {
  patient = filter(patient_metadata, patient_identifier==argv$patient_identifier)
} else {
  print("Need patient identifier")
  quit()
}

print("Match WGS")
print(paste0("Running: patient: ",patient$patient_identifier," composite: ",run_composite))


## START 

## Default paths are used in case no path is provided
# Read matching intervals and total intervals 
# Read annotation table

if(argv$matching_intervals=="") {
  argv$matching_intervals = paste0(base_dir,matching_intervals_outfile,patient$patient_identifier,".tsv")
}
if(argv$total_intervals=="") {
  argv$total_intervals = paste0(base_dir,total_matching_intervals_outfile,patient$patient_identifier,".bed")
}
if(argv$fusion_anno_table=="") {
  argv$fusion_anno_table = paste0(base_dir,fusion_annotation_outfile,patient$patient_identifier,".tsv")
}


if(length(Sys.glob(argv$matching_intervals))<1 | length(Sys.glob(argv$total_intervals))<1  | 
   length(Sys.glob(argv$fusion_anno_table))<1){
  print("Missing matching intervals, total intervals, fusion table")
    quit()
}
  
# Load WGS breakpoints

if(argv$vcf_somatic=="") {
  argv$vcf_somatic = get_vcf_filepath(patient$tumor_id,argv$tool,somatic=T)
}
if(argv$vcf=="") {
  argv$vcf = get_vcf_filepath(patient$tumor_id,argv$tool,somatic=F)
}

if(argv$vcf=="" ) {#|| argv$vcf_somatic == "") {
  print(paste0("MISSING ",patient$patient_identifier," ",argv$tool," ",patient$basename," ",(argv$vcf!="")," ",(argv$vcf_somatic!="")))# One or more missing files germline: ",argv$vcf, " somatic: ", argv$vcf_somatic))
  quit()
} 


matching_intervals_table = read.table(argv$matching_intervals,header=T, sep="\t")
fusion_anno_table=read.table(argv$fusion_anno_table,header=T,sep="\t")

if(!run_composite){
  total_intervals = read.table(argv$total_intervals,header=T, sep="\t")
  total_intervals = GRanges(total_intervals)
}  else {
  total_intervals = GRanges()
}

## For each tool different function needed to process genomic ranges

if(argv$tool == "manta") {
  all_gr = read_manta_sv_vcf(argv$vcf, argv$vcf_somatic)
}

if(argv$tool == "delly") {
  all_gr = read_delly_sv_vcf(argv$vcf)
}

if(argv$tool =="gridss"){
  all_gr = read_gridss_sv_vcf(argv$vcf)
}


if(length(all_gr)<1) quit("No breakpoints found")

  # Per fusion: matching to WGS breakpoints
  ## Retrieve intervals 
  ## Overlap breakpoints per matching interval and find partned ones
  ## If match found: update all_gr and fusion anno table accordingly and exit,
  ### early exit after found match in tier to prevent accumulation of rows in linking table,
  ## Else: continue with next interval
  ### annotate table with SV type 
  ### Maybe include in future: check breakpoint count regardless of partnered of adjacent intron and flanking
  ## NOTE: matching is not strand specific because  STAR fusion strand  is indicative of gene and and SV strand of SV type 

  fusion_id_lst = matching_intervals_table$identifier
  
  linking_table = data.frame(stringsAsFactors=FALSE) 
  linking_table_composite = data.frame(stringsAsFactors=FALSE) 
  additional_bp_table = data.frame(stringsAsFactors=FALSE) 
  
  
  match_bool = F 
  for (fusion_id in fusion_id_lst) {  
    linking_table_entry = c()
    linking_table_entry_composite = c()
    additional_bp_entry = c()
    
    #matching intervals strings to GRanges
    intervals = lapply(matching_intervals_table[fusion_id,intervals_lst], 
                 function(x) { if(x!="") GRanges(str_split(x,", ")[[1]]) else GRanges() })
  
    #find breakpoints in intervals
    
    ## Tier 1 adjacent intron
    gup_match_intron_bp = subsetByOverlaps(all_gr,intervals$gup_adjacent_intron, ignore.strand=T)
    gdw_match_intron_bp = subsetByOverlaps(all_gr,intervals$gdw_adjacent_intron, ignore.strand=T)
    
    partner_bp = check_partnered_bp(gup_match_intron_bp,gdw_match_intron_bp)
    linking_table_entry = get_linking_table(fusion_id,partner_bp,"intron","intron")
    
    if(length(linking_table_entry)>0) {
      linking_table = rbind(linking_table,linking_table_entry)
      match_bool=T
    }


    ## Tier 2: Alt. splice junctions
    gup_match_sj_bp = subsetByOverlaps(all_gr, intervals$gup_sjrange, ignore.strand=T) 
    gdw_match_sj_bp = subsetByOverlaps(all_gr, intervals$gdw_sjrange, ignore.strand=T)

    if(!early_exit || !match_bool ){
    ## intron and alt. splice junction (1,2) and vice versa (2,1)
    partner_bp = check_partnered_bp(gup_match_intron_bp,gdw_match_sj_bp)
    linking_table_entry = get_linking_table(fusion_id,partner_bp,"intron","sj")
    if(length(linking_table_entry)>0) {
      linking_table = rbind(linking_table,linking_table_entry)
      match_bool=T
    } }
    
    
    if(!early_exit || !match_bool ){
    partner_bp = check_partnered_bp(gup_match_sj_bp,gdw_match_intron_bp)
    linking_table_entry = get_linking_table(fusion_id,partner_bp,"sj","intron")
    if(length(linking_table_entry)>0) {
      linking_table = rbind(linking_table,linking_table_entry)
      match_bool=T
    } }
    
    ## tier 2,2
    if(!early_exit || !match_bool ) {
    partner_bp = check_partnered_bp(gup_match_sj_bp,gdw_match_sj_bp)
    linking_table_entry = get_linking_table(fusion_id,partner_bp,"sj","sj")
    if(length(linking_table_entry)>0) {
      linking_table = rbind(linking_table,linking_table_entry)
      match_bool=T
    }  }
    
    ## Tier 3 flanking region
    gup_match_flank_bp =subsetByOverlaps(all_gr, intervals$gup_flanking,ignore.strand=T)
    gdw_match_flank_bp = subsetByOverlaps(all_gr, intervals$gdw_flanking, ignore.strand=T)
    
    ##report bp counts regardless of partnered 
    #row$gup_tier_flank_cnt = length(gup_match_flank_bp)
    #row$gdw_tier_flank_cnt = length(gdw_match_flank_bp)
    
    if(!early_exit || !match_bool ){
    ## Intron with flanking and vice versa: (1,fl) and (fl,1)
    partner_bp = check_partnered_bp(gup_match_intron_bp,gdw_match_flank_bp)
    
    linking_table_entry = get_linking_table(fusion_id,partner_bp,"intron","flank")
    if(length(linking_table_entry)>0) {
      linking_table = rbind(linking_table,linking_table_entry)
      match_bool=T
    }  }
    
    if(!early_exit || !match_bool ){
    partner_bp = check_partnered_bp(gup_match_flank_bp,gdw_match_intron_bp)
    linking_table_entry = get_linking_table(fusion_id,partner_bp,"flank","intron")
    if(length(linking_table_entry)>0) {
      linking_table = rbind(linking_table,linking_table_entry)
      match_bool=T
    } }
    
    ## SJ with flanking
    if(!early_exit || !match_bool ){
    partner_bp = check_partnered_bp(gup_match_sj_bp,gdw_match_flank_bp)
    linking_table_entry = get_linking_table(fusion_id,partner_bp,"sj","flank")
    if(length(linking_table_entry)>0) {
      linking_table = rbind(linking_table,linking_table_entry)
      match_bool=T
    } }
    
    if(!early_exit || !match_bool ){
    partner_bp = check_partnered_bp(gup_match_flank_bp,gdw_match_sj_bp)
    linking_table_entry = get_linking_table(fusion_id,partner_bp,"flank","sj")
    if(length(linking_table_entry)>0) {
      linking_table = rbind(linking_table,linking_table_entry)
      match_bool=T
    } }  
      
    ##Tier flanking:
    if(!early_exit || !match_bool ){
    partner_bp = check_partnered_bp(gup_match_flank_bp,gdw_match_flank_bp)
    linking_table_entry = get_linking_table(fusion_id,partner_bp,"flank","flank")
    if(length(linking_table_entry)>0) {
      linking_table = rbind(linking_table,linking_table_entry)
      match_bool=T
    } }   

    if(!early_exit || !match_bool ){
  
    ## Tier 4: adjacent gene body 
    gup_match_gb_bp = subsetByOverlaps(all_gr, intervals$gup_adjacent_genebody, ignore.strand=T) 
    gdw_match_gb_bp = subsetByOverlaps(all_gr, intervals$gdw_adjacent_genebody, ignore.strand=T)
    
    ## intron and gene body 
    partner_bp = check_partnered_bp(gup_match_intron_bp,gdw_match_gb_bp)
    linking_table_entry = get_linking_table(fusion_id,partner_bp,"intron","genebody")
    if(length(linking_table_entry)>0) {
      linking_table = rbind(linking_table,linking_table_entry)
      match_bool=T
    } } 
    
    if(!early_exit || !match_bool ){
    partner_bp = check_partnered_bp(gup_match_gb_bp,gdw_match_intron_bp)
    linking_table_entry = get_linking_table(fusion_id,partner_bp,"genebody","intron")
    if(length(linking_table_entry)>0) {
      linking_table = rbind(linking_table,linking_table_entry)
      match_bool=T
    } } 
    
    ## flank and gene body 
    if(!early_exit || !match_bool ){
    partner_bp = check_partnered_bp(gup_match_flank_bp,gdw_match_gb_bp)
    linking_table_entry = get_linking_table(fusion_id,partner_bp,"flank","genebody")
    if(length(linking_table_entry)>0) {
      linking_table = rbind(linking_table,linking_table_entry)
      match_bool=T
    } } 
    
  if(!early_exit || !match_bool ){
    partner_bp = check_partnered_bp(gup_match_gb_bp,gdw_match_flank_bp)
    linking_table_entry = get_linking_table(fusion_id,partner_bp,"genebody","flank")
    if(length(linking_table_entry)>0) {
      linking_table = rbind(linking_table,linking_table_entry)
      match_bool=T
    } } 
    
    ## sj and gene body 
    if(!early_exit || !match_bool ){
    partner_bp = check_partnered_bp(gup_match_sj_bp,gdw_match_gb_bp)
    linking_table_entry = get_linking_table(fusion_id,partner_bp,"sj","genebody")
    if(length(linking_table_entry)>0) {
      linking_table = rbind(linking_table,linking_table_entry)
      match_bool=T
    } } 

    if(!early_exit || !match_bool ){
    partner_bp = check_partnered_bp(gup_match_gb_bp,gdw_match_sj_bp)
    linking_table_entry = get_linking_table(fusion_id,partner_bp,"genebody","sj")
    if(length(linking_table_entry)>0) {
      linking_table = rbind(linking_table,linking_table_entry)
      match_bool=T
    } } 
    ## genebody 
    
    if(!early_exit || !match_bool ){
    partner_bp = check_partnered_bp(gup_match_gb_bp,gdw_match_gb_bp)
    linking_table_entry = get_linking_table(fusion_id,partner_bp,"genebody","genebody")
    if(length(linking_table_entry)>0) {
      linking_table = rbind(linking_table,linking_table_entry)
      match_bool=T
    } } 
    
    ## Tier 5 full gene 
    gup_match_fullgene_bp =subsetByOverlaps(all_gr, intervals$gup_gene_coordinates,ignore.strand=T)
    gdw_match_fullgene_bp = subsetByOverlaps(all_gr, intervals$gdw_gene_coordinates, ignore.strand=T)
    
    ##report bp counts regardless of partnered 
    #row$gup_tier_fullgene_cnt = length(gup_match_fullgene_bp)
    #row$gdw_tier_fullgene_cnt = length(gdw_match_fullgene_bp)
    
    if(!early_exit || !match_bool ){
    partner_bp = check_partnered_bp(gup_match_fullgene_bp,gdw_match_fullgene_bp)
    linking_table_entry = get_linking_table(fusion_id,partner_bp,"fullgene","fullgene")
    if(length(linking_table_entry)>0) {
      linking_table = rbind(linking_table,linking_table_entry)
      match_bool=T
    } } 
    
    #matching_results_table = rbind(matching_results_table, row, stringsAsFactors = FALSE)
    
    
    ## composite here
    if(run_composite){
      if(!early_exit || !match_bool ){
        #collect partners of breakpoints in close intervals (and names) 
        gup_close_partners = c(gup_match_intron_bp,gup_match_sj_bp,gup_match_flank_bp)
        gup_partner_lst = gup_close_partners$partner
        
        gdw_close_partners = c(gdw_match_intron_bp,gdw_match_sj_bp,gdw_match_flank_bp) 
        gdw_partner_lst = gdw_close_partners$partner
        
        gup_partner_lst = gup_partner_lst[!is.na(gup_partner_lst) & gup_partner_lst %in% names(all_gr) 
                                          & !gup_partner_lst %in% gdw_partner_lst]
        gdw_partner_lst = gdw_partner_lst[!is.na(gdw_partner_lst) & gdw_partner_lst %in% names(all_gr)
                                          & !gdw_partner_lst %in% gup_partner_lst]
        #get the  breakpoints from names
        gup_partners = all_gr[gup_partner_lst,]
        gdw_partners = all_gr[gdw_partner_lst,]
        
        #make flank intervals of the partners and look if they overlap
        #then check the actual distance
        allowed_distance=5000 #10kb
        gup_partners_flank = flank(all_gr[gup_partner_lst,],allowed_distance,both = TRUE)
        gdw_partners_flank = flank(all_gr[gdw_partner_lst,],allowed_distance,both = TRUE)
        
        #this returns if the ends (partners) overlap of the bps initiated in matching interval
        gup_composite = unique(subsetByOverlaps(gup_partners_flank,gdw_partners_flank))
        if(length(gup_composite)>0) {
          #match! and look the other way to get bp names
          gdw_composite = unique(subsetByOverlaps(gdw_partners_flank,gup_partners_flank))
          
          #use the distance between these "ends" 
          gup_matching_bp_partner = gup_partners[names(gup_composite),]
          gdw_matching_bp_partner = gdw_partners[names(gdw_composite),]
          
          gup_matching_bp_partner = gup_matching_bp_partner[!duplicated(gup_matching_bp_partner$sourceId),]
          gup_matching_bp_partner = gup_matching_bp_partner[!names(gup_matching_bp_partner) %in% names(gdw_matching_bp_partner)]
          
          gdw_matching_bp_partner = gdw_matching_bp_partner[!duplicated(gdw_matching_bp_partner$sourceId),]
          gdw_matching_bp_partner = gdw_matching_bp_partner[!names(gdw_matching_bp_partner) %in% names(gup_matching_bp_partner)]
          
          composite_distances = distanceToNearest(gup_matching_bp_partner,gdw_matching_bp_partner,select="all")
          composite_distances = as.data.frame(composite_distances)
          composite_distances = composite_distances[composite_distances$distance==min(composite_distances$distance),]
          
          gup_matching_bp_partner= gup_matching_bp_partner[composite_distances$queryHits]
          gdw_matching_bp_partner= gdw_matching_bp_partner[composite_distances$subjectHits]
          
          #get orginal bp references back
          # close partners made from intron sj flank see above
          gup_matching_bp = gup_close_partners[gup_matching_bp_partner$partner]
          gdw_matching_bp = gdw_close_partners[gdw_matching_bp_partner$partner]
          
          #remove duplicates 
          #gup_matching_bp = gup_matching_bp[!names(gup_matching_bp) %in% names(gdw_matching_bp)]
          #gup_matching_bp = gup_matching_bp[!duplicated(gup_matching_bp$sourceId)]
          
          #gdw_matching_bp = gdw_matching_bp[!names(gdw_matching_bp) %in% names(gup_matching_bp)]
          #gdw_matching_bp = gdw_matching_bp[!duplicated(gdw_matching_bp$sourceId)]
        
          if(length(gup_matching_bp)>0 & length(gdw_matching_bp)>0) {
            partner_bp=GRangesList("gup"=gup_matching_bp, "gdw"=gdw_matching_bp)
            linking_table_entry = get_linking_table(fusion_id,partner_bp,"composite","composite")
            
            #assumes partnered so only one SV type and now no way to tell apart the partners of each SV -> need a more elaborate table
            #note that gup bp name and gdw bp name are NOT partnered 
            #also if not same number gup/gdw then it crashes 
            partner_bp=GRangesList("gup"=gup_matching_bp_partner, "gdw"=gdw_matching_bp_partner)
            linking_table_entry_composite = get_linking_table(fusion_id,partner_bp,"composite","composite")
            linking_table_entry_composite$composite_distance = composite_distances$distance
          
            linking_table = rbind(linking_table,linking_table_entry)
            linking_table_composite = rbind(linking_table_composite,linking_table_entry_composite)
            
          }
        }
      }
      ### ENDOF COMPOSITE
    }
  }

  # Output fusion supporting gr
  ## Output linking table

  if(nrow(linking_table)>0) {
  fusion_supporting_gr = all_gr[c(linking_table$gup_bp_name,linking_table$gdw_bp_name)]
  write.table(fusion_supporting_gr,paste0(analysis_dir,supporting_breakpoints_outfile,patient$patient_identifier,".",argv$tool,".bed"),quote = FALSE,sep = "\t",row.names=FALSE)

   #annotate with fusion name to make it clearer
  linking_table = linking_table %>% left_join(fusion_anno_table[,c("identifier","fusion_name")],by=c("fusion_id"="identifier"))
  }
  
  if(nrow(linking_table_composite)>0) {
  fusion_supporting_gr_composite = all_gr[c(linking_table_composite$gup_bp_name,linking_table_composite$gdw_bp_name)]
  write.table(fusion_supporting_gr_composite,paste0(analysis_dir,supporting_breakpoints_composite_outfile,patient$patient_identifier,".",argv$tool,".bed"),quote = FALSE,sep = "\t",row.names=FALSE)
  #annotate with fusion name to make it clearer
  linking_table_composite = linking_table_composite %>% left_join(fusion_anno_table[,c("identifier","fusion_name")],by=c("fusion_id"="identifier"))
  write.table(linking_table_composite,paste0(analysis_dir,linking_table_composite_outfile,patient$patient_identifier,".",argv$tool,".tsv"),quote = FALSE,sep = "\t",row.names=FALSE)
  
  }
  
  

#linking_table[is.na(linking_table)]=0
write.table(linking_table,paste0(analysis_dir,linking_table_outfile,patient$patient_identifier,".",argv$tool,".tsv"),quote = FALSE,sep = "\t",row.names=FALSE)

