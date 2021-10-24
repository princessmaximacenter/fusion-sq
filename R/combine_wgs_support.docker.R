## Fusion sq
## Combine WGS support from Delly, Manta, GRIDSS
## Last update: 2021-08-03
## 
## Update 2021-08-03: Slight refactor and include test data
## Update 2021-04-06: Merged ranges first pintersect
## Update 2021-03-30: SV type function applied to the output of every tool. Towards improved SV  integration
## Tx selection and remove imprecise matching bp if precise exist
### remove overlapping adjacent introns if cannot be resolved, also overlapping gene body if not in nearby intervals sv 
## Supporting SVs: using ranges use 50% reciprocal overlap & SV type match to determine if SVs are the same between tools
### Merged SVs are made for grouping but not used at the moment, stick to the individul SVs called by the tools
## Fusion-level summary from matching bps data
## Remove overlapping gene body

if(FALSE){
  #set paths for hpc
  source("/hpc/pmc_gen/ivanbelzen/github_fusion_sq/fusion-sq/R/default.conf")
  ## HPC config overrides
  source("/hpc/pmc_gen/ivanbelzen/github_fusion_sq/fusion-sq/R/hpc.default.conf")
  #HPC doesnt use argparser but patient specific config instead 
  source("/hpc/pmc_gen/ivanbelzen/structuralvariation/sv_functional_analysis/run/wilms_v2_20210923/wilms_v2_20210923.conf")
  source("/hpc/pmc_gen/ivanbelzen/structuralvariation/sv_functional_analysis/run/wilms_v2_20210923/wilms_v2_20210923.PMCID418AAA.conf")
}

suppressPackageStartupMessages({
  library(tidyverse, quietly=TRUE)
  library(stringr, quietly=TRUE)
  library(stringdist, quietly=TRUE)
  #library(argparser, quietly=TRUE)
  library(GenomicRanges, quietly=TRUE)
  library(dplyr)
  library(stringi)
  
})

source(paste0(script_dir,"functions.general.R"))
source(paste0(script_dir,"functions.svs.R"))


if(patient$patient_id =="") {
  print("patient$patient_id needs to be specified")
  #quit()
}

#order of arguments matters
map_template_vars=c('${input_dir}'=input_dir,'${output_dir}'=output_dir,'${cohort_identifier}'=cohort_identifier,'${cohort_wdir}'=cohort_wdir,'${patient_basename}'=patient$basename)

base_dir = stri_replace_all_fixed(base_dir_template,names(map_template_vars), map_template_vars,vectorize=F)
analysis_dir = stri_replace_all_fixed(analysis_dir_template,names(map_template_vars), map_template_vars,vectorize=F)
#need analysis dir for linking tables
reports_dir =  stri_replace_all_fixed(reports_dir_template,names(map_template_vars), map_template_vars,vectorize=F)

print("Combine WGS support")
print(paste0("Running: patient: ",patient$patient_identifier))


## INPUT
fusion_anno_table_path = paste0(base_dir,fusion_annotation_outfile,patient$patient_identifier,".tsv")
transcript_table_path = paste0(base_dir,transcript_table_outfile,patient$patient_identifier,".tsv")
matching_intervals_path = paste0(base_dir,matching_intervals_outfile,patient$patient_identifier,".tsv")

## OUTPUT
fusion_level_results_path = paste0(reports_dir,fusion_level_results_outfile,patient$patient_identifier,".tsv")
matching_bp_path = paste0(reports_dir,matching_results_outfile,patient$patient_identifier,".tsv")
fusion_tx_selection_path = paste0(reports_dir,fusion_tx_selection_outfile,patient$patient_identifier,".tsv")
pairwise_overlap_merged_path = paste0(reports_dir,pairwise_overlap_merged_outfile,patient$patient_identifier,".tsv")
supporting_svs_path = paste0(reports_dir,supporting_svs_outfile,patient$patient_identifier,".tsv")

## Prevent overwriting
#if(length(Sys.glob(fusion_level_results_path)==1)){ next()}


## Matching bps dataframe ##
# Load fusion anno table 

## read table to append to matching bp and allow for summarizing properties but not saved here
fusion_anno_table=read.table(fusion_anno_table_path,header=T,sep="\t",stringsAsFactors = F) 

## Load linking tables and supporting bp from tools

## from the individual tools match wgs output linking tables
## composite would need also those supporting bp because if unpartered => complex sv and cannot be resolved
matching_bps = data.frame(stringsAsFactors=FALSE)
matching_bps_composite= data.frame(stringsAsFactors=FALSE)
supporting_bps = data.frame(stringsAsFactors=FALSE)

for (tool in c("manta","delly","gridss")) {
  linking_table_path = paste0(analysis_dir,linking_table_outfile,patient$patient_identifier,".",tool,".tsv")
  if( length(Sys.glob(linking_table_path))<1 ) {
    next()
  } else if (file.size(linking_table_path)<5 ) { 
    next() 
  }  
  linking_table = read.table(linking_table_path,header=T,sep="\t",stringsAsFactors = F)
  if(nrow(linking_table)==0) { next() }
  
  linking_table$tool = tool    
  matching_bps = rbind(matching_bps, as.data.frame(linking_table))
  
  supporting_breakpoints_path = paste0(analysis_dir,supporting_breakpoints_outfile,patient$patient_identifier,".",tool,".bed")
  if( length(Sys.glob(supporting_breakpoints_path))==0  ) { next() }
  supporting_bp_tool = read.table(supporting_breakpoints_path,header=T,sep="\t",stringsAsFactors = F) 
  if(tool!="manta") {  supporting_bp_tool[, c("somatic")]=NA }
  supporting_bps = rbind(supporting_bps,supporting_bp_tool)
  
  ## Composite below:
  supporting_breakpoints_composite_path = paste0(analysis_dir,supporting_breakpoints_composite_outfile,patient$patient_identifier,".",tool,".bed")
  if( length(Sys.glob(supporting_breakpoints_composite_path))==0  ) { next() }
  supporting_bp_composite_tool = read.table(supporting_breakpoints_composite_path,header=T,sep="\t",stringsAsFactors = F) 
  if(tool!="manta") {  supporting_bp_tool[, c("somatic")]=NA }
  
  supporting_bps = rbind(supporting_bps,supporting_bp_composite_tool)
  
  linking_table_composite_path = paste0(analysis_dir,linking_table_composite_outfile,patient$patient_identifier,".",tool,".tsv")
  if( length(Sys.glob(linking_table_composite_path))==0  ) { next() }
  linking_table_composite = read.table(linking_table_composite_path,header=T,sep="\t",stringsAsFactors = F)
  
  matching_bps_composite = rbind(matching_bps_composite, as.data.frame(linking_table_composite))

  
}

##  patient has no fusions 
if(nrow(matching_bps)==0) { 
  print("Reporting summary: patient has no fusions")
  #next()
  quit()
}

## Annotate Matching bps 

# annotate with overlap variable
matching_intervals_table = read.table(matching_intervals_path,header=T, sep="\t",stringsAsFactors = F)
matching_bps = matching_bps %>% left_join(matching_intervals_table[,c("identifier","overlap_gup_gdw_adjacent_intron","overlap_gup_gdw_genebody")], by=c("fusion_id"="identifier"))

## Annotate with fusion properties from SF (Fusion overview)
matching_bps = matching_bps %>% left_join( fusion_anno_table[, c("identifier",
                                                                 "gup_ensembl_id","gdw_ensembl_id",
                                                                 "gup_gene_type","gdw_gene_type")],
                                           by=c("fusion_id"="identifier")) 


### Match transcripts to SVs / Transcript selection
## selection of SV and fusion prediction are intertwined with selection of transcript

## For intron validated fusions, predictions get the exact SV match (per fusion name)
### make sure it matches the adjacent intron interval of an exon boundary matching tx for both gup/gdw (strand does not matter)
# check for overlap between adjacent introns and exclude these bp
# keep the normal matching bps if no match could be found but annotate with specific_sv = FALSE and drop to gene body instead of intron
## dont add back those annotated as overlapping adj introns

# Dont consider as criteria factors like AF or bp distance 
# transcript_type can also be retained intron or lnc or ... dont filter on protein coding

transcript_table =  read.table(transcript_table_path,header=T,sep="\t",stringsAsFactors = F) %>%dplyr::rename(fusion_id = identifier)
transcript_table = transcript_table %>% filter(exon_boundary & !is.na(adjacent_intron)) 

#temporary identifier
matching_bps = matching_bps %>% mutate(fusion_tool = paste0(fusion_name,"_",tool))

fusions_to_specify = matching_bps %>% filter( (gup_location=="intron") | 
                                                (gdw_location=="intron"))

matching_bps_selection =  data.frame(stringsAsFactors=FALSE)
fusion_tx_selection_summary =  data.frame(stringsAsFactors=FALSE)


for(test_fusion_name in unique(fusions_to_specify$fusion_name)) {
  test_fusion = matching_bps %>% filter(fusion_name==test_fusion_name)
  
  gup_fusion_tx_selection = filter_tx_sv(transcript_table,test_fusion,upstream=T)
  gdw_fusion_tx_selection = filter_tx_sv(transcript_table,test_fusion,upstream=F)
  
  gup_fusion_tx_selection = filter_tx_properties(gup_fusion_tx_selection)
  gdw_fusion_tx_selection = filter_tx_properties(gdw_fusion_tx_selection)
  
  fusion_prediction_selection = intersect(gup_fusion_tx_selection$fusion_id,gdw_fusion_tx_selection$fusion_id)
  
  if(is.null(fusion_prediction_selection) | length(fusion_prediction_selection)==0) { next() }
  gup_fusion_tx_selection = gup_fusion_tx_selection %>% filter(fusion_id %in% fusion_prediction_selection)
  gdw_fusion_tx_selection = gdw_fusion_tx_selection %>% filter(fusion_id %in% fusion_prediction_selection)
  
  ## Filter SVs to only fall in these transcripts
  test_fusion = test_fusion %>% filter(fusion_id %in% fusion_prediction_selection)
  
  gup_sv = GRanges(test_fusion$gup_coordinate)
  gup_sv$bp_name = test_fusion$gup_bp_name
  gup_sv = subsetByOverlaps(gup_sv,GRanges(gup_fusion_tx_selection$adjacent_intron),ignore.strand=T)
  
  gdw_sv = GRanges(test_fusion$gdw_coordinate)
  gdw_sv$bp_name = test_fusion$gdw_bp_name
  gdw_sv = subsetByOverlaps(gdw_sv,GRanges(gdw_fusion_tx_selection$adjacent_intron),ignore.strand=T)
 
  ## remove overlapping adjacent interval
  if(any(as.character(gup_sv@seqnames)==as.character(gdw_sv@seqnames))) {
    overlap = subsetByOverlaps(GRanges(gup_fusion_tx_selection$adjacent_intron),GRanges(gdw_fusion_tx_selection$adjacent_intron))
    gup_sv = gup_sv[!gup_sv$bp_name %in% subsetByOverlaps(gup_sv,overlap)$bp_name]
    gdw_sv = gdw_sv[!gdw_sv$bp_name %in% subsetByOverlaps(gdw_sv,overlap)$bp_name]
  }
  
  
  # Only  keep the matching bps that comply to both, to keep them partnered
  
  test_fusion = test_fusion %>% filter(gup_bp_name %in% gup_sv$bp_name & gdw_bp_name %in% gdw_sv$bp_name)
  
  if(nrow(test_fusion)==0) { next() }
  matching_bps_selection = rbind(matching_bps_selection,test_fusion)
  
  ## Output the transcript selection
  # we sometimes cant further distinguish based on SV and thats OK
  # could look at TSL but then the fusion can also disturb splicing ofcourse
  
  gup_fusion_tx_selection$ensembl_id = unique(test_fusion$gup_ensembl_id)
  gup_fusion_tx_selection$fusion_name = test_fusion_name
  gdw_fusion_tx_selection$ensembl_id = unique(test_fusion$gdw_ensembl_id)
  gdw_fusion_tx_selection$fusion_name = test_fusion_name
  
  gup_fusion_tx_selection$upstream = TRUE
  gdw_fusion_tx_selection$upstream = FALSE
  
  fusion_tx_selection = rbind(gup_fusion_tx_selection,gdw_fusion_tx_selection)
  fusion_tx_selection_summary = rbind(fusion_tx_selection_summary,fusion_tx_selection)
}



write.table(fusion_tx_selection_summary,fusion_tx_selection_path, quote = FALSE,sep = "\t",row.names=FALSE)


## Add back fusions without specific SV too 
if(nrow(matching_bps_selection)>0) { 
  matching_bps_selection$specific_sv = TRUE
}

## remove overlapping gene body if not in nearby intervals sv 
matching_bps_residual =  matching_bps %>% filter(!fusion_tool %in% matching_bps_selection$fusion_tool) %>% 
  filter(overlap_gup_gdw_adjacent_intron==FALSE) %>% 
  filter(!overlap_gup_gdw_genebody | (gup_location %in% close_intervals & gdw_location %in% close_intervals))


# if intron/intron coding/coding not specific sv then set the validation interval to intron_consensus
matching_bps_residual[matching_bps_residual$fusion_tool %in% filter(fusions_to_specify,gup_gene_type=="protein_coding")$fusion_tool, c("gup_location")]="intron_consensus" 
matching_bps_residual[matching_bps_residual$fusion_tool %in% filter(fusions_to_specify,gdw_gene_type=="protein_coding")$fusion_tool, c("gdw_location")]="intron_consensus" 


matching_bps = rbind(matching_bps_selection,  
                     matching_bps_residual %>% mutate(specific_sv = FALSE))

#Cleanup: remove helper variable
matching_bps = matching_bps %>% select(-fusion_tool)

if(nrow(matching_bps)==0) { 
  print("Reporting summary: patient has no fusions")
  #next()
  quit()
}


## endof matching bps and tx selection


## Analyse underlying SV

#SV has 2 paradigms: bp and ranges
#1) Partner ranges for intra chr DEL DUP INV, flanking 50 bp region for CTX bps
#2) find overlaps partner ranges
#3) Merge overlapping to new ranges and merge metadata => potentially save
#4) Group by merge ID to know which variants are the same between tools
  
  
  matching_bps_sv_range_cols = c("gup_bp_name","gdw_bp_name","fusion_id","fusion_name","gup_location","gdw_location")

  #Store properties for annotation of SVs later
  ## Add composite df which contains the partners, matching wgs only does composite if nothing else is found
  if(nrow(matching_bps_composite)>0){
    matching_bps_sv_anno = rbind(matching_bps[,matching_bps_sv_range_cols],matching_bps_composite[,matching_bps_sv_range_cols])
  } else {
    matching_bps_sv_anno = matching_bps[,matching_bps_sv_range_cols]
  }

  #Subset svs to the bps in matching bps
  ## these are the best matching SVs based on the transcripts 
  ## Note: sv type is reannotated here
  supporting_bps = supporting_bps %>% filter(bp_name %in% matching_bps_sv_anno$gup_bp_name | bp_name %in% matching_bps_sv_anno$gdw_bp_name)
  supporting_bps[supporting_bps$tool!="manta",c("somatic")] = NA
  supporting_bps = unique(supporting_bps)
  supporting_bps[,c("svtype","bp_name","partner","tool")]=endoapply(supporting_bps[,c("svtype","bp_name","partner","tool")],as.character)
  supporting_bp_gr = df_to_gr(supporting_bps)
  supporting_bp_gr$svtype = get_svtype(supporting_bp_gr)
  
  # Loop through supporting_bp_gr and make partner ranges out of it - these get a new ID
  ## Partner ranges are defined based on min start/ max end of bp and its partner. 
  ## Ignore unpartnered and interchromosomal events -> return unchanged

  ## metadata columns which are kept 
  #partner field contains both after range making
  
  supporting_sv_metadata_cols =  c("sourceId",  "svtype", "svLen", "partner","FILTER",
                                          "insLen",  "tumor_af", "normal_af",  "somatic", "tool")
  supporting_svs = GRanges()
  for(gr_id in names(supporting_bp_gr)){
    supporting_bp_target = supporting_bp_gr[gr_id]
    range = make_partner_range(supporting_bp_target,supporting_bp_gr,supporting_sv_metadata_cols)
    supporting_svs[range$bp_name] = range
  }
  
  ## For the CTX and other single breakpoints: resize if <30 
  ## rename to supporting_svs_ranges because want to get back to the original ones later
  supporting_svs_ranges = supporting_svs
  supporting_svs_ranges[width(supporting_svs_ranges)<30] = flank(supporting_svs_ranges[width(supporting_svs_ranges)<30],width = 15,both=T)
  
  #map partnered ranges to bp 
  map_sv_range_bp_name= as.data.frame(supporting_svs) %>% select(bp_name) %>% 
    separate(col=bp_name,into=c("bp_name_head","bp_name_tail"),sep = "--",remove = F,fill = "left") %>% 
    select(bp_name,bp_name_head,bp_name_tail) %>%dplyr::rename(sv_name = bp_name) %>%
    gather(key = "bp_name_orient",value="bp_name",-sv_name) %>% select(sv_name,bp_name)
    
  map_sv_range_bp_name[is.na(map_sv_range_bp_name$bp_name),c("bp_name")]=map_sv_range_bp_name[is.na(map_sv_range_bp_name$bp_name),c("sv_name")]
  map_sv_range_bp_name = unique(map_sv_range_bp_name)
  
  
  # Find Same SVs
  ## which of the SVs overlap >50% and make a merged/reduced SV 

  #returns pairwise overlaps between SVs and the sv merged they have in common
  # also contains sv merged coordinate and the overlap fractions between each with merged and with eachother
  ## NB: currently not using the merged SV itself in this pipeline, but the identifier to group similar/same SV events
  overlap_merged = find_same_sv(supporting_svs_ranges,
                                supporting_svs_ranges,
                                reciprocal_overlap = 0.5,svtype_matching = T,ignore_strand = F)
  ##If a range is by itsef and has no overlap => automatically added back because we take the supporting_svs_df as base
  ## so it is NOT in the overlap_merged set => if sv_merged  == NA then set to bp_name and overlap_merged_bp == keep NA 
  
  
  
  ## Build supporting SV dataframe
  if(length(overlap_merged)>0){
    ## check if sv only assigned once 
    if(nrow(unique(overlap_merged[,c("set1","sv_merged")]))!=length(unique(overlap_merged$set1))) {
      uq_merged = unique(overlap_merged[,c("set1","sv_merged")])
      print( overlap_merged[overlap_merged$set1 %in% uq_merged[duplicated(uq_merged$set1),c("set1")],])
      
      print("WARNING SV assigned multiple times")
      print(patient$patient_identifier)
      break
    }
    
    write.table(overlap_merged,pairwise_overlap_merged_path,quote = FALSE,sep = "\t",row.names=FALSE)
    
    supporting_svs_df = as.data.frame(supporting_svs) %>% 
      left_join(overlap_merged[,c("set1","sv_merged","sv_merged_coordinate","overlap_merged_set1","overlap_set1_merged")], by=c("bp_name"="set1")) 
    
    } else {
    supporting_svs_df = as.data.frame(supporting_svs)
    supporting_svs_df$sv_merged = NA
    supporting_svs_df$overlap_merged_set1 = NA
    supporting_svs_df$overlap_set1_merged = NA
    supporting_svs_df$sv_merged_coordinate = NA
  }
  
  supporting_svs_df = unique(supporting_svs_df)
  supporting_svs_df$coordinate = paste0(supporting_svs_df$seqnames,":",supporting_svs_df$start,"-",supporting_svs_df$end,":",supporting_svs_df$strand)
  
  
  ## set NA ranges to bp name
  supporting_svs_df[is.na(supporting_svs_df$sv_merged),c("sv_merged")]=
    supporting_svs_df[is.na(supporting_svs_df$sv_merged),c("bp_name")]
  supporting_svs_df[is.na(supporting_svs_df$sv_merged_coordinate),c("sv_merged_coordinate")]=
    supporting_svs_df[is.na(supporting_svs_df$sv_merged_coordinate),c("coordinate")]
  
  supporting_svs_df = annotate_variant_af_class(supporting_svs_df)
  
  
  ## Add annotation to SVs 
  
  ## Aim: get fusion name and predictions per SV (bp pair)
  ## First fusion name per bp /side
  matching_bp_long = matching_bps_sv_anno %>%
    gather(key = "bp_name_orient",value="bp_name",-fusion_id,-fusion_name,-gup_location,-gdw_location)  
    
  ## map sv range => bp name. Join with bpname level annotation and summarize
  ## group by sv name, so both partner halves come together again, 
  ## Also works if no merged ranges found
  ## note, not yet for location
  map_sv_range_anno = map_sv_range_bp_name %>% left_join(matching_bp_long,by="bp_name") %>% group_by(sv_name) %>%
    summarize(fusion_name = toString(sort(unique(fusion_name))), fusion_predictions = toString(sort(unique(fusion_id))))
  
  #add the 'sv name' also explicitly
  supporting_svs_df$sv_name=supporting_svs_df$bp_name
  
  supporting_svs_df = 
    supporting_svs_df %>% left_join(map_sv_range_anno, by=c("sv_name"))
  
  
  write.table(supporting_svs_df,supporting_svs_path,quote = FALSE,sep = "\t",row.names=FALSE)



  ## Fusion level summary
  
  # use supporting sv properties (sv type, AF, somatic/germline/low
  # use matching bps for gene id gene type etc gup location     
  
  #properties that belong to the fusion prediction
  matching_bps_reporting_cols = c("fusion_id","fusion_name",
                                  "gup_bp_name", "gdw_bp_name",# "gup_distance", "gdw_distance",
                                  "gup_location","gdw_location","overlap_gup_gdw_genebody","specific_sv", 
                                  "gup_ensembl_id","gdw_ensembl_id",
                                  "gup_gene_type","gdw_gene_type")
                                                
  
  supporting_sv_fusion_cols = c("svtype", "svLen", "tumor_af", "normal_af", "tool",
                                   "sv_merged", "sv_merged_coordinate", "sv_name","coordinate")
  
  
  ## SV length
  ## Note: Manta and GRIDSS have  negative svLen for deletions, Delly not. Also Delly can have svlen for CTX? => set to 0 explicitly 
  supporting_svs_df[supporting_svs_df$svtype=="CTX",c("svLen")]=NA
  supporting_svs_df$svLen = abs(supporting_svs_df$svLen)
  
  ## harmonize and merge composite if exists
  if(nrow(matching_bps_composite)>0){
    matching_bps_composite =  matching_bps_composite %>% left_join( fusion_anno_table[, c("identifier","gup_ensembl_id","gdw_ensembl_id",
                                                                "gup_gene_type","gdw_gene_type")],
                                                            by=c("fusion_id"="identifier")) 
    matching_bps_composite[,names(matching_bps)[!names(matching_bps) %in% names(matching_bps_composite)]]=NA
    matching_bps_composite$specific_sv = FALSE
    matching_bps_composite = matching_bps_composite[,names(matching_bps_composite) %in% names(matching_bps)]
    matching_bps_composite$gup_distance=NA
    matching_bps_composite$gdw_distance=NA
  }
  matching_bps2 = rbind(matching_bps,matching_bps_composite)
  
  #remove columns better filled with SV properties
  matching_bps2 = matching_bps2[,matching_bps_reporting_cols]
  
  ## Annotate with fusion properties from SF (Fusion overview)
  matching_bps2 = matching_bps2 %>% left_join( fusion_anno_table[, c("identifier","FFPM",
                                                                     "gup_gene_id","gdw_gene_id",
                                                                     "gup_sf_breakpoint","gdw_sf_breakpoint",
                                                                     "gup_sf_transcript","gdw_sf_transcript",
                                                                     "predicted_frame")],
                                               by=c("fusion_id"="identifier")) 
  
  ## Merge with supporting svs
  matching_bps2 = matching_bps2 %>% 
    left_join(map_sv_range_bp_name,by=c("gup_bp_name"="bp_name")) %>% 
    left_join(supporting_svs_df[,supporting_sv_fusion_cols],by=c("sv_name")) %>%
    dplyr::rename_at(supporting_sv_fusion_cols,function(x){paste0("gup_",x)}) %>%
    left_join(map_sv_range_bp_name,by=c("gdw_bp_name"="bp_name")) %>% 
    left_join(supporting_svs_df[,supporting_sv_fusion_cols],by=c("sv_name")) %>%
    dplyr::rename_at(supporting_sv_fusion_cols,function(x){paste0("gdw_",x)})
  
  
  
  ## distance should be based on merged coordinate and sf_breakpoint
  matching_bps2$gup_gdw_sv_distance = distance(GRanges(matching_bps2$gup_sv_merged_coordinate),
                                               GRanges(matching_bps2$gdw_sv_merged_coordinate),ignore.strand=T)
  
  matching_bps2$gup_distance = distance(GRanges(matching_bps2$gup_sv_merged_coordinate),
                                        GRanges(matching_bps2$gup_sf_breakpoint),ignore.strand=T)
  
  matching_bps2$gdw_distance = distance(GRanges(matching_bps2$gdw_sv_merged_coordinate),
                                        GRanges(matching_bps2$gdw_sf_breakpoint),ignore.strand=T)
                                        
 
  # Per fusion know supporting SV and if that is found by multiple tools (as that SV type as well), group by merged sv
  ## SV properties if not composite than these should match between gup/gdw
  ## merge tumor gup/gdw normal gup/gdw AF, tools, sv type, sv length
  ## merging takes care of grouping so do not explicitly merge first per tool anymore. 

fusion_level_svs = matching_bps2 %>% group_by(fusion_name,gup_sv_merged,gdw_sv_merged,gup_sv_merged_coordinate,gdw_sv_merged_coordinate,specific_sv) %>% 
  summarize(gup_sf_breakpoint = toString(unique(sort(gup_sf_breakpoint))), gdw_sf_breakpoint = toString(unique(sort(gdw_sf_breakpoint))),
            gup_sf_transcript = toString(unique(sort(gup_sf_transcript))),gdw_sf_transcript = toString(unique(sort(gdw_sf_transcript))),
            gup_gene_id = toString(unique(gup_gene_id)), gdw_gene_id = toString(unique(gdw_gene_id)),
            gup_gene_type = toString(unique(gup_gene_type)), gdw_gene_type = toString(unique(gdw_gene_type)),
            gup_ensembl_id=unique(gup_ensembl_id),gdw_ensembl_id=unique(gdw_ensembl_id),
            ffpm_mean=mean(FFPM), ffpm_max=max(FFPM),
            predicted_frame = toString(sort(unique(predicted_frame))),
            fusion_predictions = toString(unique(sort(fusion_id))),
            gup_location=toString(unique(sort(gup_location))), gdw_location=toString(unique(sort(gdw_location))),
            overlap_gup_gdw_genebody = any(overlap_gup_gdw_genebody),
            
            sv_names = toString(unique(sort(c(as.character(gup_sv_name),as.character(gdw_sv_name))))),
            gup_coordinate = toString(unique(sort(gup_coordinate))), gdw_coordinate = toString(unique(sort(gdw_coordinate))),
            tools=toString(unique(sort(c(gup_tool,gdw_tool)))),
            tumor_af=mean(c(gup_tumor_af,gdw_tumor_af),na.rm=T), normal_af=mean(c(gup_normal_af,gdw_normal_af),na.rm=T),
            tumor_af_spread=(max(tumor_af,na.rm=T)-min(tumor_af,na.rm=T)), normal_af_spread=(max(normal_af,na.rm=T)-min(normal_af,na.rm=T)), 
            svtype=toString(unique(sort(c(gup_svtype,gdw_svtype)))), 
            svlen=ifelse(!grepl("CTX",svtype),mean(c(gup_svLen,gdw_svLen),na.rm=T),NA),
            svlen_spread=ifelse(!is.na(svlen),(max(c(gup_svLen,gdw_svLen),na.rm=T)-min(c(gup_svLen,gdw_svLen),na.rm = T)),NA),
            
            gup_gdw_sv_distance_mean = mean(gup_gdw_sv_distance,na.rm=T),
            gup_distance_mean = mean(gup_distance,na.rm=T), gdw_distance_mean = mean(gdw_distance,na.rm=T),
            
            .groups="keep") %>% ungroup() 

  
  fusion_level_svs$predicted_frame = trimws(gsub(pattern = "., ", replacement = "",gsub(pattern = ", .", replacement = "", 
                                                                                        fusion_level_svs$predicted_frame,fixed = T),fixed = T))
  
  
  #rename tumor/normal_af_mean afterwards but needed for  annotate somatic/germline/low_af
  ## if after mean still NAs then only NAs present, replace to prevent NAs with variant classification
  fusion_level_svs[is.na(fusion_level_svs$tumor_af),c("tumor_af")]=0
  fusion_level_svs[is.na(fusion_level_svs$normal_af),c("normal_af")]=0
  
  fusion_level_svs = annotate_variant_af_class(fusion_level_svs)
  
  fusion_level_svs = fusion_level_svs %>%dplyr::rename(tumor_af_mean = tumor_af, normal_af_mean = normal_af)
  
  ## Annotate as precise and confident
  fusion_level_svs = fusion_level_svs %>% mutate(location_precise = (gup_location %in% proximate_bp & gdw_location %in% proximate_bp) )
  
  #location and keep all other columns
  cohort_tools = fusion_level_svs %>% group_by(fusion_name) %>% summarize(tools_any_wgs = toString(unique(sort(tools))))
  fusion_level_svs = fusion_level_svs %>% left_join(cohort_tools)
  fusion_level_svs = fusion_level_svs %>% mutate(precise_confident = location_precise&grepl(",",tools))
  fusion_level_svs = fusion_level_svs %>% mutate(not_precise_confident = !fusion_name %in% filter(fusion_level_svs,precise_confident)$fusion_name)

write.table(matching_bps2,matching_bp_path,quote = FALSE,sep = "\t",row.names=FALSE)
  
write.table(fusion_level_svs,fusion_level_results_path,quote = FALSE,sep = "\t",row.names=FALSE)

