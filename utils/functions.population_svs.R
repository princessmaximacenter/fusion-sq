
## Load SV database (VCF)
# and annotate with its metadata
load_sv_database_vcf = function(sv_database_path) {
  sv_database_vcf = readVcf(sv_database_path , "hg38")
  sv_database_gr = rowRanges(sv_database_vcf)
  sv_database_gr$bp_name = names(sv_database_gr)
  sv_database_gr = sv_database_gr[order(sv_database_gr$bp_name)]
  
  
  sv_database_meta = as.data.frame(info(sv_database_vcf))
  sv_database_meta$bp_name = rownames(sv_database_meta)
  sv_database_meta = sv_database_meta[order(sv_database_meta$bp_name),]
  
  meta_tmp = as.data.frame(mcols(sv_database_gr))
  meta_tmp$bp_name = rownames(meta_tmp)
  meta_tmp = meta_tmp %>% left_join(sv_database_meta,by=c("bp_name")) 
  meta_tmp = meta_tmp[order(meta_tmp$bp_name),]
  
  mcols(sv_database_gr) = meta_tmp
  
  return(sv_database_gr)
}

make_range_sv_database = function(sv_database_gr){
  end(sv_database_gr[!is.na(sv_database_gr$END)])=sv_database_gr[!is.na(sv_database_gr$END)]$END
  return(sv_database_gr)
  
  #need to make a range from the sv_database_gr using the end coordinate
  keep_without_end = sv_database_gr[is.na(sv_database_gr$END)]
  
  sv_database_ranges = GRanges()
  for(bp_name in names(sv_database_gr[!is.na(sv_database_gr$END)])) {
    sv = sv_database_gr[bp_name]
    sv_metadata = mcols(sv)
    
    range = GRanges(seqnames = seqnames(sv),
                    IRanges(start = start(sv), end = sv_metadata$END))
    mcols(range)=sv_metadata
    names(range)=bp_name
    sv_database_ranges = c(sv_database_ranges,range)
  }
  
  sv_database_ranges = c(sv_database_ranges,keep_without_end)
  return(sv_database_ranges)
}

