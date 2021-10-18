
manta_metadata = function(vcf, patient, support_type, somatic=TRUE) {
  #vcf object as input
  vcf_geno_df = as.data.frame(geno(vcf)[[support_type]])
  vcf_geno_df$sourceId = rownames(vcf_geno_df)
  
  vcf_geno_df[,as.character(patient$normal_id)] = 
    gsub("[c()]", "", vcf_geno_df[,as.character(patient$normal_id)])
  
  vcf_geno_df[,as.character(patient$normal_id)] = 
    gsub(":", ", ", vcf_geno_df[,as.character(patient$normal_id)])
  
  vcf_geno_df[vcf_geno_df == "integer0"]=NA
  
  vcf_geno_df = separate(vcf_geno_df, col = as.character(patient$normal_id), 
                         into = c(paste0("normal_REF_",support_type),paste0("normal_VAR_",support_type)), sep = ", ")
  
  #remove the ":0" if couldnt be split 
  #vcf_geno_df[,paste0("normal_REF_",support_type)] = sub(":0", "", vcf_geno_df[,paste0("normal_REF_",support_type)])
  
  if(somatic){
    vcf_geno_df[,as.character(patient$tumor_id)] = 
      gsub("[c()]", "", vcf_geno_df[,as.character(patient$tumor_id)])
    
    vcf_geno_df[,as.character(patient$tumor_id)] = 
      gsub(":", ", ", vcf_geno_df[,as.character(patient$tumor_id)])
    
    vcf_geno_df[vcf_geno_df == "integer0"]=NA
    
    vcf_geno_df = separate(vcf_geno_df, col = as.character(patient$tumor_id), 
                           into = c(paste0("tumor_REF_",support_type),paste0("tumor_VAR_",support_type)), sep = ", ")
    
    #remove the ":0" if couldnt be split 
    # vcf_geno_df[,paste0("tumor_REF_",support_type)] = sub(":0", "", vcf_geno_df[,paste0("tumor_REF_",support_type)])
  }
  
  
  #vcf_geno_df[is.na(vcf_geno_df)]=0
  
  return(vcf_geno_df)
}  

annotate_metadata = function(all_gr,metadata_df) {
  all_gr=all_gr[order(names(all_gr))]
  metadata = mcols(all_gr,use.names = T)  
  metadata$bp_name = rownames(metadata)
  
  if("sourceId" %in% names(metadata_df) && nrow(dplyr::filter(metadata_df,sourceId %in% metadata$sourceId))>0) {
    metadata = metadata %>% merge(metadata_df %>% dplyr::select(sourceId,tumor_af,normal_af,somatic),by="sourceId",all.x=T) %>% unique()
  } 
  else if("bp_name" %in% names(metadata_df) && nrow(dplyr::filter(metadata_df,bp_name %in% metadata$bp_name))>0) {
    metadata = metadata %>% merge(metadata_df %>% dplyr::select(bp_name,tumor_af,normal_af,somatic),by="bp_name",all.x=T) %>% unique()
  }
  metadata=metadata[order(metadata$bp_name),]
  rownames(metadata) = metadata$bp_name  
  
  mcols(all_gr) <- metadata
  return(all_gr)
}

## for each tool different function needed to load wgs
read_manta_sv_vcf = function(vcf_germline_path, vcf_somatic_path) {
  
  # SOMATIC and GERMLINE are separate files for somatic and germline breakpoints
  if(!run_composite){
    somatic_vcf = readVcf(vcf_somatic_path, "hg38",param=total_intervals) 
  } else {
    somatic_vcf = readVcf(vcf_somatic_path, "hg38")
  }
  #SR and PR separate dataframes at first
  somatic_vcf_geno_PR = manta_metadata(somatic_vcf, patient, "PR", somatic=TRUE)
  somatic_vcf_geno_SR = manta_metadata(somatic_vcf, patient, "SR", somatic=TRUE)
  
  ## join and then use both for calculation of AF 
  somatic_vcf_geno = somatic_vcf_geno_SR %>% left_join(somatic_vcf_geno_PR,by="sourceId")
  somatic_vcf_geno$somatic=TRUE
  
  somatic_vcf_geno$tumor_af = (as.integer(somatic_vcf_geno$tumor_VAR_PR) + as.integer(somatic_vcf_geno$tumor_VAR_SR)) / 
    ( as.integer(somatic_vcf_geno$tumor_REF_SR) + as.integer(somatic_vcf_geno$tumor_REF_PR) +
        as.integer(somatic_vcf_geno$tumor_VAR_SR) + as.integer(somatic_vcf_geno$tumor_VAR_PR) )
  
  somatic_vcf_geno$normal_af = (as.integer(somatic_vcf_geno$normal_VAR_PR) + as.integer(somatic_vcf_geno$normal_VAR_SR)) / 
    ( as.integer(somatic_vcf_geno$normal_REF_SR) + as.integer(somatic_vcf_geno$normal_REF_PR) +
        as.integer(somatic_vcf_geno$normal_VAR_SR) + as.integer(somatic_vcf_geno$normal_VAR_PR) )
  
  #likely normal AF  0 but not all 
  #somatic_vcf_geno %>% filter(normal_af > tumor_af)
  
  ## germline has only normal sample 
  if(!run_composite){
    germline_vcf = readVcf(vcf_germline_path, "hg38",param=total_intervals) 
  } else {
    germline_vcf = readVcf(vcf_germline_path, "hg38")
  }
  ## join and then use both for calculation of AF 
  
  germline_vcf_geno_PR = manta_metadata(germline_vcf, patient, "PR", somatic=FALSE)
  germline_vcf_geno_SR = manta_metadata(germline_vcf, patient, "SR", somatic=FALSE)
  
  germline_vcf_geno = germline_vcf_geno_SR %>% left_join(germline_vcf_geno_PR,by="sourceId")
  germline_vcf_geno$somatic=FALSE
  
  ## lacking tumor AF so adding dummy variable for merging
  germline_vcf_geno$tumor_af = NA
  
  germline_vcf_geno$normal_af = (as.integer(germline_vcf_geno$normal_VAR_PR) + as.integer(germline_vcf_geno$normal_VAR_SR)) / 
    ( as.integer(germline_vcf_geno$normal_REF_SR) + as.integer(germline_vcf_geno$normal_REF_PR) +
        as.integer(germline_vcf_geno$normal_VAR_SR) + as.integer(germline_vcf_geno$normal_VAR_PR) )
  
  
  #breakpoints need to append and annotate with allele frequencies, tool, somatic or not 
  somatic_gr = breakpointRanges(somatic_vcf)
  germline_gr = breakpointRanges(germline_vcf)
  
  ## Note: manta sometimes 0s for normal AF in diploid model, also somtimes normal > tumor AF in somatic 
  
  ##  filtering the sourceIds for "manta" to remove "svrecord" since they ofen overlap with normal variants, but without metadata
  somatic_gr = somatic_gr[grepl("Manta",somatic_gr$sourceId),]
  somatic_vcf_geno$sourceId= gsub(".",":",somatic_vcf_geno$sourceId,fixed=T)
  
  # in case the source Id and bp name do not match then the somatic VCF geno "bp name" matches the sourceId
  ## somatic_vcf_geno[!somatic_vcf_geno$sourceId %in% somatic_gr$sourceId & !somatic_vcf_geno$sourceId %in% names(somatic_gr),]
  ### sometimes another does match with the same info like MantaDUP:TANDEM:35919:0:1:0:0:0 :1/2/3 etc
  ##somatic_vcf_geno[somatic_vcf_geno$sourceId %in% (somatic_gr[names(somatic_gr) != somatic_gr$sourceId ])$sourceId,]
  
  somatic_gr = annotate_metadata(somatic_gr,somatic_vcf_geno)
  
  germline_gr = germline_gr[grepl("Manta",germline_gr$sourceId),]
  germline_vcf_geno$sourceId= gsub(".",":",germline_vcf_geno$sourceId,fixed=T)
  germline_gr = annotate_metadata(germline_gr,germline_vcf_geno)
  
  all_gr = c(somatic_gr,germline_gr)
  #add another column for analysis type
  #annotate with tool
  mcols(all_gr)[["tool"]]="manta"
  return(all_gr) 
}
read_delly_sv_vcf = function(vcf_path) {
  ##NOTE: somatic disabled because  DV/DR/RV/RR and allele frequencies are these same
  ## Function is available for tagging based on name but not currently used by Fusion-sq 
  
  germline_vcf = readVcf(vcf_path, "hg38") 
  if(!run_composite){ 
    germline_vcf = subsetByOverlaps(germline_vcf,total_intervals)
  }
  #prevent error
  germline_vcf = germline_vcf[!isNA(info(germline_vcf)$SVTYPE)]
  info(germline_vcf)[info(germline_vcf)$SVTYPE=="BND",c("SVTYPE")]="TRA"
  
  ## Make dataframe based on this
  all_gr = breakpointRanges(germline_vcf) 
  
  ## sourceId for row names, which splits bp_nm into _bp1 and _bp2
  gr_df = as.data.frame(all_gr) %>% mutate(bp_name = rownames(.)) %>% select(sourceId)
  
  # #sourceid for somatic  
  germline_vcf_geno = as.data.frame(info(germline_vcf)) %>% mutate(sourceId = rownames(.)) %>% select(sourceId, PRECISE)
  
  #filter based on if occurs in germline_gr => some unpaired were removed during loading
  germline_vcf_geno = germline_vcf_geno %>% filter(sourceId %in% gr_df$sourceId)
  
  ## Calculate AF: how to do this depends on precise/imprecise, which is in info tag 
  #AF = RV/(RR+RV) for precise variants
  #AF = DV/(DR+DV) for imprecise variants
  germline_vcf_geno_precise = 
    filter(germline_vcf_geno,PRECISE) %>% 
    left_join (
      as.data.frame(geno(germline_vcf)$RV) %>% mutate(sourceId = rownames(.)) %>% 
        rename("tumor_RV" := !!patient$tumor_id) %>% rename("normal_RV" := !!patient$normal_id)
      , by="sourceId") %>%
    left_join(  
      as.data.frame(geno(germline_vcf)$RR) %>% mutate(sourceId = rownames(.)) %>%
        rename("tumor_RR" := !!patient$tumor_id) %>% rename("normal_RR" := !!patient$normal_id)
      , by="sourceId"
    ) %>%
    mutate(tumor_af = (tumor_RV/(tumor_RR+tumor_RV))) %>%
    mutate(normal_af = (normal_RV/(normal_RR+normal_RV)))
  
  #if both are 0 then drop the NAs, they can be rescued by DVs or otherwise we just dont have data for that variant
  germline_vcf_geno_precise = drop_na(germline_vcf_geno_precise)
  
  # calculate based on DV for all to also rescue the NAs with RV
  germline_vcf_geno_all = 
    germline_vcf_geno %>% 
    left_join (
      as.data.frame(geno(germline_vcf)$DV) %>% mutate(sourceId = rownames(.)) %>% 
        rename("tumor_DV" := !!patient$tumor_id) %>% rename("normal_DV" := !!patient$normal_id)
      , by="sourceId") %>%
    left_join(  
      as.data.frame(geno(germline_vcf)$DR) %>% mutate(sourceId = rownames(.)) %>%
        rename("tumor_DR" := !!patient$tumor_id) %>% rename("normal_DR" := !!patient$normal_id)
      , by="sourceId"
    ) %>%
    mutate(tumor_af = (tumor_DV/(tumor_DR+tumor_DV))) %>%
    mutate(normal_af = (normal_DV/(normal_DR+normal_DV)))
  
  germline_vcf_geno_all = drop_na(germline_vcf_geno_all)
  
  #if AF in precise table, then exclude it here 
  germline_vcf_af = rbind(germline_vcf_geno_precise[,c("sourceId","tumor_af","normal_af")], 
                          germline_vcf_geno_all[!germline_vcf_geno_all$sourceId %in% germline_vcf_geno_precise$sourceId,c("sourceId","tumor_af","normal_af")])
  
  
  ## Currently not used:
  gr_df$somatic = NA
  
  #merge allele frequencies, 
  gr_df = gr_df %>% left_join(germline_vcf_af, by="sourceId")
  
  #add the annotation to breakpoints
  all_gr = annotate_metadata(all_gr,gr_df)
  
  #annotate with tool
  mcols(all_gr)[["tool"]]="delly"
  
  return(all_gr)
}
read_gridss_sv_vcf = function(vcf_path) {
  ##NOTE: somatic disabled because allele frequencies are these same
  ## Function is available for tagging based on name but not currently used by Fusion-sq 
  
  if(!run_composite) {
    germline_vcf = readVcf(vcf_path, "hg38", param=total_intervals) 
  } else {
    germline_vcf = readVcf(vcf_path, "hg38")
  }
  
  all_gr = breakpointRanges(germline_vcf)
  elementMetadata(all_gr)["svtype"] = simpleEventType(all_gr)
  
  #this removes unpartnered, need to subset afterwards
  germline_vcf = germline_vcf[names(all_gr)]
  
  germline_vcf_geno = geno(germline_vcf)
  germline_gr_df = as.data.frame(all_gr) %>% mutate(bp_name = rownames(.))
  
  ## Allele frequency calculation
  #how to calculate AF depends on variant size 
  # svlen <1000 dont use the REFPAIR otherwise do 
  #uses the label of tumor/normal not the id.
  
  germline_anno = germline_gr_df %>% select(bp_name,svLen) %>% left_join(
    as.data.frame(germline_vcf_geno$VF) %>% mutate(bp_name = rownames(.)) %>% 
      rename("tumor_VF" := !!patient$tumor_label) %>% rename("normal_VF" := !!patient$normal_label)
    , by="bp_name") %>% left_join(
      as.data.frame(germline_vcf_geno$REF) %>% mutate(bp_name = rownames(.)) %>% 
        rename("tumor_REF" := !!patient$tumor_label) %>% rename("normal_REF" := !!patient$normal_label)
      , by="bp_name") %>% left_join(
        as.data.frame(germline_vcf_geno$REFPAIR) %>% mutate(bp_name = rownames(.)) %>% 
          rename("tumor_REFPAIR" := !!patient$tumor_label) %>% rename("normal_REFPAIR" := !!patient$normal_label)
        , by="bp_name") %>%
    mutate(tumor_af = ifelse(!is.na(svLen)&svLen<1000, (tumor_VF/(tumor_VF+tumor_REF)), (tumor_VF/(tumor_VF+tumor_REF+tumor_REFPAIR)))) %>%
    mutate(normal_af = ifelse(!is.na(svLen)&svLen<1000, (normal_VF/(normal_VF+normal_REF)), (normal_VF/(normal_VF+normal_REF+normal_REFPAIR))))
  
  germline_anno$somatic=NA
  
  #add the annotation to breakpoints
  all_gr = annotate_metadata(all_gr,germline_anno)
  
  #annotate with tool
  mcols(all_gr)[["tool"]]="gridss"
  return(all_gr)
}

check_partnered_bp = function(gup_matching_bp,gdw_matching_bp){
  if(length(gup_matching_bp) > 0 & length(gdw_matching_bp) > 0){  
    ##is partnered?
    gup_in_gdw = gup_matching_bp[names(gup_matching_bp) %in% gdw_matching_bp$partner,]
    gdw_in_gup = gdw_matching_bp[names(gdw_matching_bp) %in% gup_matching_bp$partner,]
    
    ##remove overlaps
    gup_in_gdw = gup_in_gdw[!names(gup_in_gdw) %in% names(gdw_in_gup)]
    gdw_in_gup = gdw_in_gup[!names(gdw_in_gup) %in% names(gup_in_gdw)]

    gup_in_gdw = gup_in_gdw[names(gup_in_gdw) %in% gdw_in_gup$partner]
    gdw_in_gup = gdw_in_gup[names(gdw_in_gup) %in% gup_in_gdw$partner]
   
    
    if(length(gup_in_gdw)>0 & length(gdw_in_gup)>0) {
      return(GRangesList("gup"=gup_in_gdw, "gdw"=gdw_in_gup))
    }  
  }
  return(GRangesList())
}

get_linking_table = function(fusion_id,partner_bp,gup_location="",gdw_location="") {
  if(length(partner_bp)==0) return( c() )
  ## for each entry in partner bp => row in linking table 
  
  #todo: how to know that bps are paired up correctly?
  #sourceId does not always work (not manta and not gridss likely)
  #check if bp name is set for all tools 
  
  ## remove row names and change col names with appending "gup/gdw" 
  ## AF separate because of GRIDSS which has different AF for breakpoints => average at patient report generation
  
  linking_table_entry = c()
  linking_table_entry$fusion_id = fusion_id
  linking_table_entry = cbind( linking_table_entry,  as.data.frame(mcols(partner_bp$gup)[,c(columns_bplevel,"somatic","svtype")])  %>% 
                                 rename_at(vars(columns_bplevel), 
                                           function(x){paste0("gup_",x)}) )
  
  rownames(linking_table_entry) = c()
  linking_table_entry = cbind( linking_table_entry,  as.data.frame(mcols(partner_bp$gdw)[,c(columns_bplevel)])  %>% 
                                 rename_at(vars(columns_bplevel), 
                                           function(x){paste0("gdw_",x)}) )
  
  rownames(linking_table_entry) = c()
  
  linking_table_entry = cbind( linking_table_entry, gup_coordinate=strsplit(toString(GRanges(seqnames=partner_bp$gup@seqnames,ranges = partner_bp$gup@ranges, strand=partner_bp$gup@strand)), ", ")[[1]])
  linking_table_entry = cbind( linking_table_entry, gdw_coordinate=strsplit(toString(GRanges(seqnames=partner_bp$gdw@seqnames,ranges = partner_bp$gdw@ranges, strand=partner_bp$gdw@strand)), ", ")[[1]])
  
  linking_table_entry = cbind(linking_table_entry, gup_distance = distance(all_gr[linking_table_entry$gup_bp_name],
                                                                           GRanges(dplyr::filter(fusion_anno_table,identifier==fusion_id)$gup_sf_breakpoint),ignore.strand=T))
  linking_table_entry = cbind(linking_table_entry, gdw_distance = distance(all_gr[linking_table_entry$gdw_bp_name],
                                                                           GRanges(dplyr::filter(fusion_anno_table,identifier==fusion_id)$gdw_sf_breakpoint),ignore.strand=T))
  
  linking_table_entry$gup_location = gup_location
  linking_table_entry$gdw_location = gdw_location
  
  return(linking_table_entry)
}


#Function copied from GRIDSS example.
#my edit: keep BND for same strand instead of inv
simpleEventType <- function(gr) {
  return(ifelse(seqnames(gr) != seqnames(gr[gr$partner]), "ITX", # inter-chromosomosal
                ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS", # TODO: improve classification of complex events
                       ifelse(strand(gr) == strand(partner(gr)), "BND",
                              ifelse(xor(start(gr) < start(partner(gr)), strand(gr) == "-"), "DEL",
                                     "DUP")))))
}

