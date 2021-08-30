
#Tx selection: filter on transcript properties
filter_tx_properties = function(tx_table) {
  if(nrow(tx_table[tx_table$transcript_type=="protein_coding",])>0) {
    tx_table = tx_table %>% dplyr::filter(transcript_type=="protein_coding")
  }  
  
  if(nrow(tx_table[grepl("CCDS|basic",tx_table$tag),])>0) {
    tx_table = tx_table %>% dplyr::filter(grepl("CCDS|basic",tag))
  }  
  return(tx_table)
}

filter_tx_sv = function(transcript_table,test_fusion,upstream=T){
  fusion_tx = transcript_table %>% dplyr::filter(fusion_id %in% test_fusion$fusion_id)
  if(upstream) {
    fusion_tx = fusion_tx %>% dplyr::filter(upstream)
    sv_coordinate = test_fusion$gup_coordinate
  } else {
    fusion_tx = fusion_tx %>% dplyr::filter(!upstream)
    sv_coordinate = test_fusion$gdw_coordinate
    
  }
  
  #make into GRanges object for matching 
  adjacent_intron = GRanges(fusion_tx$adjacent_intron)
  adjacent_intron$transcript_id = fusion_tx$transcript_id
  adjacent_intron_subset = subsetByOverlaps(adjacent_intron, GRanges(sv_coordinate),ignore.strand=T)
  
  #if matching succeeded -> tx selection
  if(length(adjacent_intron_subset)>0){
    fusion_tx_selection = fusion_tx[fusion_tx$transcript_id %in% adjacent_intron_subset$transcript_id,]
    return(fusion_tx_selection)
  } else {
    return(data.frame())
  }
}


## Functions to load vcf

load_manta_diploidSV = function(vcf_path,sample) {
  
  germline_vcf = readVcf(vcf_path, "hg38")
  
  
  ## join and then use both for calculation of AF 
  germline_vcf_geno_PR = manta_metadata(germline_vcf, sample, "PR", somatic=FALSE)
  germline_vcf_geno_SR = manta_metadata(germline_vcf, sample, "SR", somatic=FALSE)
  
  germline_vcf_geno = germline_vcf_geno_SR %>% left_join(germline_vcf_geno_PR,by="sourceId")
  
  ## lacking tumor AF so adding dummy variable for merging
  germline_vcf_geno$tumor_af = NA
  germline_vcf_geno$somatic = NA
  germline_vcf_geno$normal_af = (as.integer(germline_vcf_geno$normal_VAR_PR) + as.integer(germline_vcf_geno$normal_VAR_SR)) / 
    ( as.integer(germline_vcf_geno$normal_REF_SR) + as.integer(germline_vcf_geno$normal_REF_PR) +
        as.integer(germline_vcf_geno$normal_VAR_SR) + as.integer(germline_vcf_geno$normal_VAR_PR) )
  
  #breakpoints need to append and annotate with allele frequencies, tool, somatic or not 
  germline_gr = breakpointRanges(germline_vcf)
  
  ## Note: manta sometimes 0s for normal AF in diploid model, also somtimes normal > tumor AF in somatic 
  
  ##  filtering the sourceIds for "manta" to remove "svrecord" since they ofen overlap with normal variants, but without metadata
  
  # in case the source Id and bp name do not match then the somatic VCF geno "bp name" matches the sourceId
  germline_gr = germline_gr[grepl("Manta",germline_gr$sourceId),]
  germline_vcf_geno$sourceId= gsub(".",":",germline_vcf_geno$sourceId,fixed=T)
  germline_gr = annotate_metadata(germline_gr,germline_vcf_geno)
  return(germline_gr)
}

###


get_svtype <- function(gr) {
  #Function copied from GRIDSS example code
  # CTX if seq names are the same
  ## return as complex if unpartnered
  
  return_gr = gr
  unpartnered_gr=gr[!gr$partner %in% names(gr)]
  gr = gr[gr$partner %in% names(gr)]
  
  svtype = ifelse(seqnames(gr) != seqnames(gr[gr$partner]), "CTX", # inter-chromosomosal
                  ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS", # TODO: improve classification of complex events
                         ifelse(strand(gr) == strand(gr[gr$partner]), "INV",
                                ifelse(xor(start(gr) < start(gr[gr$partner]), strand(gr) == "-"), "DEL",
                                       "DUP"))))
  
  if(length(unpartnered_gr)>0){
    return_gr[names(unpartnered_gr)]$svtype="complex"
  }
  return_gr[names(gr)]$svtype=svtype
  return(  return_gr$svtype)
}


annotate_variant_af_class = function(summary_df) {
  summary_df = summary_df %>% 
    mutate(tumor_normal_diff = (tumor_af-normal_af)) %>%
    mutate(tumor_normal_ratio = (tumor_normal_diff/normal_af)) 
  
  summary_df = summary_df %>% mutate(somatic_variant = (tumor_normal_diff>0.05 & tumor_normal_ratio>1.5))
  summary_df = summary_df %>% mutate(germline_variant = (normal_af>0.05 & tumor_normal_ratio<1.1))
  summary_df = summary_df %>% mutate(low_af = (tumor_normal_diff<0.05 & normal_af<0.05))
  
  summary_df = summary_df %>% mutate(somatic_variant = ifelse( (is.na(tumor_af)&!is.na(normal_af)), FALSE,somatic_variant))
  
  summary_df = summary_df %>% dplyr::select(-tumor_normal_diff)
  return(summary_df)
}

df_to_gr = function(set) {
  if(is.data.frame(set)) { 
    set = unique(set)
    row.names(set) = set$bp_name
    set =GRanges(set) 
  }
  return(set)
}

make_partner_range = function(supporting_bp_target, supporting_bp_gr, metadata_cols= c("sourceId",  "svtype", "svLen", "insLen","partner")) {
  ## function for single gr, first test if partner exists
  if(supporting_bp_target$partner %in% names(supporting_bp_gr)) {
    partner_gr = supporting_bp_gr[supporting_bp_target$partner]
  
    if(as.character(seqnames(supporting_bp_target))==as.character(seqnames(partner_gr))) {
    # for gridss the AF will not be the same and use the mean
    metadata = mcols(supporting_bp_target)[metadata_cols] 
    if(supporting_bp_target$tool=="gridss"){
      if("tumor_af" %in% metadata_cols) { metadata$tumor_af = mean(c(metadata$tumor_af,mcols(partner_gr)$tumor_af),na.rm=T) }
      if("normal_af" %in% metadata_cols) { metadata$normal_af = mean(c(metadata$normal_af,mcols(partner_gr)$normal_af),na.rm=T) }
    }
    both_names = sort(c(supporting_bp_target$bp_name,supporting_bp_target$partner))
    metadata$partner = toString(both_names)
    ## name the range to both grs in the order of the new range?  => partner also contains original bp
    
    min_start = min(start(supporting_bp_target),start(partner_gr))
    max_end = max(end(supporting_bp_target),end(partner_gr))
    
    range = GRanges(seqnames = seqnames(supporting_bp_target),
                    IRanges(start = min_start, end = max_end))
    mcols(range) =  metadata
    range$bp_name = paste0(both_names,collapse ="--")
    names(range) = range$bp_name
    return(range)
  } 
  }
  
    #else if partner not found or interchromosomal
    range = supporting_bp_target
    mcols(range)=mcols(supporting_bp_target)[metadata_cols]
    range$bp_name = supporting_bp_target$bp_name
    names(range)=range$bp_name
    return(range)
  
}


make_range_svs = function(svs_gr,sv_metadata_cols= c("sourceId",  "svtype", "svLen", "partner","insLen")){
  #sv_metadata_cols =  c("sourceId",  "svtype", "svLen", "partner", "insLen",  "tumor_af", "normal_af",  "somatic", "tool")
  
  svs = GRanges()
  for(gr_id in names(svs_gr)){
    supporting_bp_target = svs_gr[gr_id]
    range = make_partner_range(supporting_bp_target,svs_gr,sv_metadata_cols)
    svs[range$bp_name] = range
  }
  mcols(svs)$coordinate = paste0(seqnames(svs),":",start(svs),"-",end(svs))
  return(svs)
}

overlap_fraction_gr = function(set1,set2,ignore_strand=TRUE){
  overlaps_intersect = pintersect(set1, set2,ignore.strand=ignore_strand)
  overlap_fraction = width(overlaps_intersect) / width(set1)
  return(overlap_fraction)
}

get_reciprocal_overlap_pairs = function(set1,set2,reciprocal_overlap=0.5,svtype_matching=TRUE,ignore_strand=TRUE){  
  if(is.null(names(set1))) {
    names(set1)=paste0("set1_",1:length(set1))
  }
  if(is.null(names(set2))) {
    names(set2)=paste0("set2_",1:length(set2))
  }
  # make pairs, so note that this makes the lists larger
  overlap_pairs = findOverlaps(set1,set2,type="any",minoverlap = 1,ignore.strand=ignore_strand)
  
  if(svtype_matching==TRUE){
    svtype_match = set1[overlap_pairs@from]$svtype==set2[overlap_pairs@to]$svtype
    overlap_pairs = overlap_pairs[svtype_match==TRUE]
  }
  
  #remove if not at least overlaps => reciprocal overlap
  overlap_fraction = overlap_fraction_gr(set1[overlap_pairs@from], set2[overlap_pairs@to],ignore_strand) #width of intersect/width of first
  overlap_pairs = overlap_pairs[overlap_fraction > reciprocal_overlap]
  overlap_fraction = overlap_fraction_gr(set2[overlap_pairs@to], set1[overlap_pairs@from],ignore_strand)
  overlap_pairs = overlap_pairs[overlap_fraction > reciprocal_overlap]
  
  #length of intersect divided by length of first in function 
  overlap_fraction_metadata = overlap_fraction_gr(set1[overlap_pairs@from],set2[overlap_pairs@to],ignore_strand) 
  overlap_fraction_metadata_2 = overlap_fraction_gr(set2[overlap_pairs@to],set1[overlap_pairs@from],ignore_strand)
  overlap_pairs_df = as.data.frame(overlap_pairs)
  colnames(overlap_pairs_df) = c("from","to")
    
  overlap_pairs_df$set1 = names(set1[overlap_pairs@from])
  if("svtype" %in% names(mcols(set1))) {
    overlap_pairs_df$set1_svtype = set1[overlap_pairs@from]$svtype
  }
  
  overlap_pairs_df$set2 = names(set2[overlap_pairs@to])
  if("svtype" %in% names(mcols(set2))) {
    overlap_pairs_df$set2_svtype = set2[overlap_pairs@to]$svtype
  }
  overlap_pairs_df$overlap_set1_set2 = overlap_fraction_metadata
  overlap_pairs_df$overlap_set2_set1 = overlap_fraction_metadata_2
  
  #remove overlaps with self
  overlap_pairs_df = overlap_pairs_df[overlap_pairs_df$set1!=overlap_pairs_df$set2,]
  
  return(overlap_pairs_df)
}

find_same_sv = function(set1,set2,reciprocal_overlap=0.5,svtype_matching=TRUE,ignore_strand=FALSE){  
  if(is.data.frame(set1)) { set1 =df_to_gr(set1) }
  if(is.data.frame(set2)) { set2 =df_to_gr(set2) }
  
  overlap_pairs_df = get_reciprocal_overlap_pairs(set1,set2,reciprocal_overlap=reciprocal_overlap,
                                               svtype_matching=svtype_matching,ignore_strand=ignore_strand)
  overlap_pairs_bk = overlap_pairs_df
  #if same sets then make sure not the same overlap twice
  same_sets = FALSE
  if(all(set1$bp_name==set2$bp_name)){
    same_sets = TRUE
  }
  
  ## Merge SVs: reduce to the overlapping GRanges
  ## next map back to SV bps to group them, but make this first because sets are going to have their name reset
  # construct the dataframe from set1/2 intersection and then map those intersections to the collapsed ones to prevent multimapping rows
  
  sv_merged = GRanges()
  for(query in unique(overlap_pairs_df$from)) {
    hits = overlap_pairs_df[overlap_pairs_df$from==query,c("to")]
    if(length(hits)==0){next()}
    #sv_merged_entry = GenomicRanges::reduce(c(set1[query],set2[hits]),ignore.strand=ignore_strand)
    sv_merged_entry = GenomicRanges::intersect(set1[query],set2[hits],ignore.strand=ignore_strand)
    mcols(sv_merged_entry)$svtype= unique(set1[query]$svtype)
    sv_merged = c(sv_merged,sv_merged_entry)
    
    if(same_sets) {
      overlap_pairs_df = overlap_pairs_df[ !(overlap_pairs_df$from %in% hits & overlap_pairs_df$to==query),]
    }
  }
  
  
  if(length(sv_merged)==0) { return(GRanges())} 
  
  ## Only reduce the merged ranges that have high overlap with each other
  names(sv_merged)=paste0("merged_tmp_",1:length(sv_merged))
  
  overlap_between_merged = get_reciprocal_overlap_pairs(sv_merged,sv_merged,reciprocal_overlap=reciprocal_overlap,
                               svtype_matching=svtype_matching,ignore_strand=ignore_strand)
  overlapping_sv_merged_names = unique(c(overlap_between_merged$set1,overlap_between_merged$set2))
  compact_sv_merged = GRanges()
  for(query in unique(overlap_between_merged$from)) {
    hits = overlap_between_merged[overlap_between_merged$from==query,c("to")]
    if(length(hits)==0){next()}
    sv_merged_entry = GenomicRanges::reduce(c(sv_merged[query],sv_merged[hits]),ignore.strand=ignore_strand)
    mcols(sv_merged_entry)$svtype= unique(sv_merged[query]$svtype)
    compact_sv_merged = c(compact_sv_merged,sv_merged_entry)
    
    overlap_between_merged = overlap_between_merged[ !(overlap_between_merged$from %in% hits & overlap_between_merged$to==query),]
  }
 
  #add the unmapped back too
  sv_merged = unique(c(compact_sv_merged,sv_merged[!names(sv_merged) %in% overlapping_sv_merged_names]))
  sv_merged = trim(sv_merged)
  names(sv_merged)=paste0("merged_",1:length(sv_merged))
  
  # Map SVs to merged
  ## require reciprocal overlap with merged range
  ## require SV type match
  
  map_sv_merged_set1_df = get_reciprocal_overlap_pairs(sv_merged,set1,reciprocal_overlap=reciprocal_overlap,
                               svtype_matching=svtype_matching,ignore_strand=ignore_strand)
  
  
  #overlap merged to set 1 is larger than set 1 to merged, because merged ranges are larger than the original ones
  map_sv_merged_set1_df = map_sv_merged_set1_df %>% select(-to,-from) %>%
    dplyr::rename(sv_merged = set1,
           sv_merged_svtype = set1_svtype,
           set1 = set2,
           set1_svtype = set2_svtype,
           overlap_merged_set1 = overlap_set1_set2,
           overlap_set1_merged = overlap_set2_set1)
  
  
  ## assign bp to sv range that is best overlapping
  map_sv_merged_set1_overlaps =  map_sv_merged_set1_df%>% group_by(set1,sv_merged) %>% summarize(overlap_mean =mean(c(overlap_merged_set1,overlap_set1_merged)))
  map_sv_merged_set1_best = map_sv_merged_set1_overlaps %>% group_by(set1) %>% summarize(max_overlap = max(as.numeric(overlap_mean)))
  map_sv_merged_set1_uq = map_sv_merged_set1_overlaps %>% merge(map_sv_merged_set1_best,by.x=c("set1","overlap_mean"), by.y=c("set1","max_overlap") )
  
  map_sv_merged_set1_df = map_sv_merged_set1_df %>% merge(map_sv_merged_set1_uq[,c("set1","sv_merged")],by=c("set1","sv_merged"))
 
  # duplicates => remove because identical
  map_sv_merged_set1_df=map_sv_merged_set1_df[!duplicated(map_sv_merged_set1_df$set1),]
  #check no duplicates map_sv_merged_set1_df[map_sv_merged_set1_df$set1 %in% map_sv_merged_set1_df[duplicated(map_sv_merged_set1_df$set1),c("set1")],]
  
  if(!same_sets){
    map_sv_merged_set2_df = get_reciprocal_overlap_pairs(sv_merged,set2,reciprocal_overlap=reciprocal_overlap,
                                                       svtype_matching=svtype_matching,ignore_strand=ignore_strand)
  
    #overlap merged to set 1 is larger than set 1 to merged, because merged ranges are larger than the original ones
    map_sv_merged_set2_df = map_sv_merged_set2_df %>% select(-to,-from) %>%
      dplyr::rename(sv_merged = set1,
             sv_merged_svtype = set1_svtype,
             overlap_merged_set2 = overlap_set1_set2,
             overlap_set2_merged = overlap_set2_set1)
    
    map_sv_merged_set2_overlaps =  map_sv_merged_set2_df%>% group_by(set2,sv_merged) %>% summarize(overlap_mean =mean(c(overlap_merged_set2,overlap_set2_merged)))
    map_sv_merged_set2_best = map_sv_merged_set2_overlaps %>% group_by(set2) %>% summarize(max_overlap = max(as.numeric(overlap_mean)))
    map_sv_merged_set2_uq = map_sv_merged_set2_overlaps %>% merge(map_sv_merged_set2_best,by.x=c("set2","overlap_mean"), by.y=c("set2","max_overlap") )
    map_sv_merged_set2_df = map_sv_merged_set2_df %>% merge(map_sv_merged_set2_uq[,c("set2","sv_merged")],by=c("set2","sv_merged"))
    
    # duplicates => remove because identical
    map_sv_merged_set2_df=map_sv_merged_set2_df[!duplicated(map_sv_merged_set2_df$set2),]
    
  } else {
    map_sv_merged_set2_df = map_sv_merged_set1_df %>% rename_with(function(x){str_replace(x,"set1","set2")})
  }
  
  ## Make metadata frame
  ##map bp to overlaps and use that to build the metadata dataframe
  ## then you can merge the rows based on the reduced genomic ranges for certain properties like AF 
  # go back to the pairs that overlapped >50% amd sv type match
  
  bp_mapping = overlap_pairs_bk %>% select(-to,-from)
  
    #this second intersection is for the breakpoint mapping between the tools because does not have to be one to one 
  bp_mapping = bp_mapping %>% 
    left_join(map_sv_merged_set1_df[,c("sv_merged","sv_merged_svtype","set1","overlap_merged_set1","overlap_set1_merged")],by="set1") %>% 
    dplyr::rename(sv_merged_set1 = sv_merged)
  bp_mapping = unique(bp_mapping)
  
  bp_mapping = bp_mapping %>% 
    left_join(map_sv_merged_set2_df[,c("sv_merged","set2","overlap_merged_set2","overlap_set2_merged")],by="set2") %>%
    dplyr::rename(sv_merged_set2 = sv_merged)
  
  bp_mapping = unique(bp_mapping)
  
  ## should be identical otherwise not merged
  bp_mapping = bp_mapping %>% mutate(sv_merged = ifelse((sv_merged_set1==sv_merged_set2),sv_merged_set1, NA))
  bp_mapping = bp_mapping %>% dplyr::filter(!is.na(sv_merged)) %>% select(-sv_merged_set1,-sv_merged_set2)
  
  
  ## prevent same bp to multi ranges / select best matching group
  multimapping_merged = bp_mapping %>% group_by(sv_merged) %>%
    summarize(svs = toString(unique(sort(set1))),overlap_mean =mean(c(overlap_merged_set1,overlap_set1_merged)))
  sv_groups = multimapping_merged %>% group_by(svs) %>% summarize(max_overlap = max(overlap_mean))
  multimapping_merged_uq = multimapping_merged %>% merge(sv_groups, by.x=c("svs","overlap_mean"), by.y=c("svs","max_overlap") )
  
  #if still the same left then remove because identical
  multimapping_merged_uq = multimapping_merged_uq[!duplicated(multimapping_merged_uq$svs),]
  multimapping_merged_uq$newnames = paste0("merged_",1:nrow(multimapping_merged_uq))
  
  bp_mapping = bp_mapping %>% merge(multimapping_merged_uq[,c("newnames","sv_merged")],by=c("sv_merged")) 
  
  # add coordinates of merged
  sv_merged_coordinates = as.data.frame(sv_merged)
  sv_merged_coordinates$bp_name=row.names(sv_merged_coordinates)
  sv_merged_coordinates$sv_merged_coordinate = paste0(sv_merged_coordinates$seqnames,":",sv_merged_coordinates$start,"-",sv_merged_coordinates$end)
  if(!ignore_strand) {
    sv_merged_coordinates$sv_merged_coordinate = paste0(sv_merged_coordinates$sv_merged_coordinate,":",sv_merged_coordinates$strand)
  }
  bp_mapping= bp_mapping %>% left_join(sv_merged_coordinates[,c("bp_name","sv_merged_coordinate")],by=c("sv_merged"="bp_name"))
  
  #use newnames of selected ranges to fill out the gaps
  bp_mapping = bp_mapping %>% select(-sv_merged) %>% dplyr::rename(sv_merged = newnames)
  
  return(bp_mapping)
}


annotate_svs_genes = function(gr,gene_properties) {
  sv_genes = findOverlaps(gr,gene_properties,ignore.strand=T)
  gr$gene_name=NA
  gr$ensembl_id=NA
  
  multiple_genes = sv_genes[duplicated(queryHits(sv_genes))]
  single_gene = sv_genes[!queryHits(sv_genes) %in% queryHits(multiple_genes)]
  
  gr[queryHits(single_gene)]$gene_name = gene_properties[subjectHits(single_gene)]$gene_name
  gr[queryHits(single_gene)]$ensembl_id = gene_properties[subjectHits(single_gene)]$gene_id
  
  for(query in unique(queryHits(multiple_genes))) {
    hits = sv_genes[queryHits(sv_genes)==query]
    gr[query]$gene_name =  toString(unique(gene_properties[subjectHits(hits)]$gene_name))
    gr[query]$ensembl_id = toString(unique(gene_properties[subjectHits(hits)]$gene_id))
  }
  
  return(gr)
}


annotate_svs_properties = function(svs,properties,cols,use_mean=F){
  properties_sv_hits = findOverlaps(properties,svs,ignore.strand=T)
  mcols(svs)[,cols]=NA
  for(query in unique(subjectHits(properties_sv_hits))) {
    hits = properties_sv_hits[subjectHits(properties_sv_hits)==query]
    hits = properties[queryHits(hits)]
    if(use_mean==F){
      mcols(svs[query])[,cols] = toString(sort(unique(mcols(hits)[,cols])))
    } else{
      mcols(svs[query])[,cols] = mean(mcols(hits)[,cols],na.rm=T)
    }
  }
  return(svs)
}


get_tx_by_genes = function(gtf,ensembl_id_lst) {
  transcripts = gtf[gtf$type=="transcript"& gtf$gene_id %in% ensembl_id_lst]
  metadata_tx = as.data.frame(mcols(transcripts))
  metadata_tx$mane_select=F
  metadata_tx[metadata_tx$transcript_id %in% tx_mane_select,c("mane_select")]=TRUE
  
  protein_max_exon = as.data.frame(gtf[!is.na(gtf$protein_id)&gtf$gene_id %in% ensembl_id_lst]) %>% group_by(protein_id) %>%
    summarize(max_exon = max(as.numeric(exon_number),na.rm=T))
  protein_max_cds = as.data.frame(gtf[!is.na(gtf$protein_id)&gtf$gene_id %in% ensembl_id_lst & gtf$type=="CDS"]) %>% group_by(protein_id) %>%
    summarize(max_cds_width = sum(width,na.rm=T))
  
  metadata_tx = metadata_tx %>% left_join(protein_max_exon,by="protein_id")
  metadata_tx = metadata_tx %>% left_join(protein_max_cds,by="protein_id")
  
  mcols(transcripts)=metadata_tx
  return(transcripts)
}

filter_tx_canonical = function(tx_table) {
  
  if(length(tx_table[tx_table$mane_select])>0) {
    tx_table = tx_table[tx_table$mane_select]
  }
  
  if(length(tx_table[grepl("basic|CCDS|apris",tx_table$tag)])>0) {
    tx_table = tx_table[grepl("basic|CCDS",tx_table$tag)]
  }
  
  if(length(tx_table[tx_table$transcript_type=="protein_coding",])>0) {
    tx_table = tx_table[tx_table$transcript_type=="protein_coding"]
  }
  
  if(length(tx_table)==1) return(tx_table)
  
  max_cds_width = as.numeric(max(tx_table$max_cds_width,na.rm=T))
  if(length(tx_table[!is.na(tx_table$max_cds_width) & tx_table$max_cds_width > max_cds_width*0.9])>0) {
    tx_table = tx_table[!is.na(tx_table$max_cds_width) & tx_table$max_cds_width > max_cds_width*0.9]
  }
  
  if(length(tx_table)==1) return(tx_table)
  
  max_exons = as.numeric(max(tx_table$max_exon,na.rm=T))*0.8
  if(length(tx_table[!is.na(tx_table$max_exons) & tx_table$max_exons > max_exons])>0) {
    tx_table = tx_table[!is.na(tx_table$max_exons) & tx_table$max_exons > max_exons]
  }
  
  if(length(tx_table)==1) return(tx_table)
  
  min_tsl = min(tx_table$transcript_support_level,na.rm=T)
  if(length(tx_table[!is.na(tx_table$transcript_support_level) & tx_table$transcript_support_level == min_tsl])>0) {
    tx_table = tx_table[!is.na(tx_table$transcript_support_level) & tx_table$transcript_support_level == min_tsl]
  }
  
  
  if(length(tx_table)==1) return(tx_table)
  
  tx_table$tx_width = width(tx_table@ranges)
  max_width = max(tx_table$tx_width,na.rm=T)
  if(length(tx_table[tx_table$tx_width == max_width])>0) {
    tx_table = tx_table[tx_table$tx_width == max_width]
  }
  mcols(tx_table) = mcols(tx_table)[names(mcols(tx_table)) != "tx_width"]
  
  
  if(length(tx_table[!is.na(tx_table$max_cds_width) & tx_table$max_cds_width == max_cds_width])>0) {
    tx_table = tx_table[!is.na(tx_table$max_cds_width) & tx_table$max_cds_width == max_cds_width]
  }
  
  return(tx_table)
}

get_features_per_tx = function(features, transcript_lst){
  
  features_per_tx=GRanges()
  #use names because sometimes empty
  for(tx_id in transcript_lst) {
    tx_features = features[[tx_id]]
    if(length(tx_features)==0){next()}
    tx_features$tx_id =tx_id
    
    tx_features$ft_rank = rank(tx_features)
    tx_features$ft_max = max(tx_features$ft_rank,na.rm=T)
    if(any(strand(tx_features)=="-")) {
      tx_features$ft_rank = (tx_features$ft_max-tx_features$ft_rank)+1
    }
  
    features_per_tx = c(features_per_tx,tx_features)
  }
  return(features_per_tx)
}

annotate_svs_features = function(gr,features_per_tx) {
  
  sv_features= findOverlaps(gr,features_per_tx,ignore.strand=T)
  gr$ft_overlap = NA
  gr$tx_id = NA
  
  multiple_features=sv_features[duplicated(queryHits(sv_features))]
  single_feature = sv_features[!queryHits(sv_features) %in% queryHits(multiple_features)]
  if(length(single_feature)>0){
    gr[queryHits(single_feature)]$ft_overlap = paste0( features_per_tx[subjectHits(single_feature)]$gene_name,":",
                                                          features_per_tx[subjectHits(single_feature)]$ft_rank, "/", features_per_tx[subjectHits(single_feature)]$ft_max)
    gr[queryHits(single_feature)]$tx_id = features_per_tx[subjectHits(single_feature)]$tx_id
  }
  
  if(length(multiple_features)>0){
    range= gr[unique(queryHits(multiple_features))]
    end(range)=start(range)
    
    sv_features_start = findOverlaps(range,features_per_tx,ignore.strand=T)
    for(query in unique(queryHits(sv_features_start))) {
      hits = sv_features_start[queryHits(sv_features_start)==query]
      
      match = features_per_tx[subjectHits(hits)]
      for(match_gene_name in unique(match$gene_name)){
        match_gene = match[match$gene_name==match_gene_name]
        if(length(unique(match_gene$tx_id))>1) { #multiple tx?? shouldt happen
          print("WARNING: multiple tx available")
          print(match_gene)
          match_gene=match_gene[match_gene$tx_id==unique(match_gene$tx_id)[1],]
        }
        if(length(match_gene)>1) { #multiple times same gene
          match_gene_ranks = paste0(min(match_gene$ft_rank),"-",max(match_gene$ft_rank))
        } else {
          match_gene_ranks = unique(match_gene$ft_rank)
        }
        
        match_feature = paste0(", start:",  match_gene_name,":",  match_gene_ranks, "/", unique(match_gene$ft_max))
        
        range[query]$ft_overlap=  paste0(range[query]$ft_overlap, match_feature)
        range[query]$tx_id = paste0(range[query]$tx_id,"start:",unique(match_gene$tx_id))
      }
    }
    gr[unique(queryHits(multiple_features))]$ft_overlap=range$ft_overlap
    gr[unique(queryHits(multiple_features))]$tx_id=range$tx_id
      
    
    
    range= gr[unique(queryHits(multiple_features))]
    start(range)=end(range)
    sv_features_end = findOverlaps(range,features_per_tx,ignore.strand=T)
    for(query in unique(queryHits(sv_features_end))) {
      hits = sv_features_end[queryHits(sv_features_end)==query]
      
      match = features_per_tx[subjectHits(hits)]
      for(match_gene_name in unique(match$gene_name)){
        match_gene = match[match$gene_name==match_gene_name]
        if(length(unique(match_gene$tx_id))>1) { #multiple tx?? shouldt happen
          print("WARNING: multiple tx available")
          print(match_gene)
          match_gene=match_gene[match_gene$tx_id==unique(match_gene$tx_id)[1],]
        }
        if(length(match_gene)>1) { #multiple times same gene
          match_gene_ranks = paste0(min(match_gene$ft_rank),"-",max(match_gene$ft_rank))
        } else {
          match_gene_ranks = unique(match_gene$ft_rank)
        }
      
        match_feature = paste0(", end:",  match_gene_name,":",  match_gene_ranks, "/", unique(match_gene$ft_max))
        
        range[query]$ft_overlap=  paste0(range[query]$ft_overlap, match_feature)
        range[query]$tx_id = paste0(range[query]$tx_id,"end:",unique(match_gene$tx_id))
      }
    }  
    
    gr[unique(queryHits(multiple_features))]$ft_overlap=range$ft_overlap
    gr[unique(queryHits(multiple_features))]$tx_id=range$tx_id
    
  }
  
  
  gr[!is.na(gr$ft_overlap)]$ft_overlap = str_replace(gr[!is.na(gr$ft_overlap)]$ft_overlap ,"NA, ", "")
  gr[!is.na(gr$tx_id)]$tx_id = str_replace(gr[!is.na(gr$tx_id)]$tx_id ,"NA, ", "")
  
  return(gr)
}
