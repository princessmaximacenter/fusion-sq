remove_version_from_id = function(identifier) {
  return(sub("\\..*","",identifier))
}

## Adjust fusion anno columns - from STAR fusion names to names used in pipeline
#because I dont save it, I do this at make patient report and collect cohort
rename_fusion_anno_columns = function(fusion_anno_table) {
  fusion_anno_table$gup_ensembl_id = remove_version_from_id(fusion_anno_table$left_ensembl_id)
  fusion_anno_table$gdw_ensembl_id = remove_version_from_id(fusion_anno_table$right_ensembl_id)
  ## add these columns to also keep versioned 
  
  fusion_anno_table = fusion_anno_table %>% 
    dplyr::rename(predicted_frame=PROT_FUSION_TYPE) %>%
    dplyr::rename(gup_sf_breakpoint = LeftBreakpoint, gdw_sf_breakpoint=RightBreakpoint) %>% 
    dplyr::rename(gup_sf_transcript = CDS_LEFT_ID, gdw_sf_transcript = CDS_RIGHT_ID) %>%
    dplyr::rename(gup_gene_id = left_gene_id) %>% dplyr::rename(gdw_gene_id = right_gene_id) %>% 
    dplyr::rename(gup_ensembl_version = left_ensembl_id) %>% dplyr::rename(gdw_ensembl_version = right_ensembl_id) %>% 
    dplyr::rename(gup_gene_type = left_gene_type) %>% dplyr::rename(gdw_gene_type = right_gene_type)
  

return(fusion_anno_table)
}



rbind_no_colmatch = function (df1,df2)  {
  if(nrow(df1)==0) {
    return(df2)
  }
  if(nrow(df2)==0) {
    return(df1)
  }
  cols_1 = names(df1)
  cols_2 = names(df2)
  col_diff_1 = cols_2[!cols_2 %in% cols_1]
  col_diff_2 = cols_1[!cols_1 %in% cols_2]
  df1[,col_diff_1]=NA
  df2[,col_diff_2]=NA
  
  join_df=rbind(df1,df2)
  return(join_df)
}

## harmonisation of identifiers

harmonize_gene_identifiers = function(df) {
  #fusioncatcher has custom ids for igh locus, but not used in my pipeline -> needs matching version 
  #remove ensembl ids for those that have no matching _version to resolve that 
  #https://github.com/ndaniel/fusioncatcher/blob/master/bin/add_custom_gene.py
  df = df %>% mutate(
    gup_gene_id = str_replace(str_replace(gup_gene_id,fixed("IGH@"),"IGH"),fixed("IGH-@-ext"),"IGH"),
    gdw_gene_id = str_replace(str_replace(gdw_gene_id,fixed("IGH@"),"IGH"),fixed("IGH-@-ext"),"IGH"),
    fusion_name = str_replace(str_replace(fusion_name,fixed("IGH@"),"IGH"),fixed("IGH-@-ext"),"IGH"))
  
  if("gup_ensembl_version" %in% names(df) & "gdw_ensembl_version" %in% names(df) ) {
    df = df %>% mutate(
      gup_ensembl_version = ifelse(grepl("ENS",gup_ensembl_version),gup_ensembl_version,NA),
      gdw_ensembl_version = ifelse(grepl("ENS",gdw_ensembl_version),gdw_ensembl_version,NA),
      gup_ensembl_id = ifelse(!is.na(gup_ensembl_version),gup_ensembl_id,NA),
      gdw_ensembl_id = ifelse(!is.na(gdw_ensembl_version),gdw_ensembl_id,NA))
    
  } else {
    #for the IGH 
    df = df %>% mutate(
      gup_ensembl_id = ifelse(gup_gene_id=="IGH",NA,gup_ensembl_id),
      gdw_ensembl_id = ifelse(gdw_gene_id=="IGH",NA,gdw_ensembl_id))
    
  }
  return(df)        
}

make_identifiers_cohort_analysis = function(df) {
  df = df %>% mutate(patient_fusion =  paste(patient_id,fusion_name,sep="_"))
  
  if("gup_sv_merged" %in% names(df) & "gdw_sv_merged" %in% names(df) ) {
    df = df %>% mutate(
      patient_fusion_sv =  paste(patient_id,gup_sv_merged,gdw_sv_merged,fusion_name,sep="_"),
      patient_sv_id =  paste(patient_id,gup_sv_merged,gdw_sv_merged,sep="_"))
  }
  
  return(df)
}
## Counting functions

uq_fusions = function(fusion_summary){
  return(nrow(unique(fusion_summary[,c("fusion_name","patient_id")])))
}

uq_uq_fusions = function(fusion_summary){
  return(length(unique(fusion_summary[,c("fusion_name")])))
}
uq_patients = function(fusion_summary){
  return(nrow(unique(fusion_summary %>% select(patient_id))))
}

### Make count overviews / recurrence / promiscuiety

fusion_cnt_per_attr = function(fusion_summary, attr = "patient_id", attr_name=NULL) {
  cnt_table = unique(fusion_summary[,c(attr,"fusion_name")]) %>% group_by(fusion_name) %>% count() 
  
  if(!is.null(attr_name)) {
    colnames(cnt_table) = c("fusion_name",attr_name)
  }
  return(cnt_table)
}

cnt_fusions = function(fusion_summary, attr_name=NULL) {
  cnt_table = fusion_summary %>% group_by(patient_id) %>% select(patient_id,fusion_name) %>% unique() %>% summarise(fusion_cnt = n(),.groups="keep") %>% ungroup()
  if(!is.null(attr_name)) {
    colnames(cnt_table) = c("patient_id",attr_name)
  }
  return(cnt_table)
}

fusion_partner_cnt_per_patient = function(fusion_summary,attr_name=NULL){
  cnt_table = unique(fusion_summary %>% select(gup_gene_id,gdw_gene_id, patient_id) %>%
                       gather(key="gene_orient",value="gene_name",gup_gene_id,gdw_gene_id) %>%
                       select(patient_id,gene_name)) %>% group_by(gene_name) %>% count() 
  
  if(!is.null(attr_name)) {
    colnames(cnt_table) = c("gene_name",attr_name)
  }
  return(cnt_table)
}

fusion_promiscuous = function(fusion_summary,attr_name=NULL){
  gup_cnt = unique(fusion_summary %>% select(gup_gene_id,gdw_gene_id)) %>% 
    group_by(gup_gene_id) %>% count() %>% rename(gup_cnt =n) %>% rename(gene_name = gup_gene_id)
  gdw_cnt = unique(fusion_summary %>% select(gup_gene_id,gdw_gene_id)) %>% 
    group_by(gdw_gene_id) %>% count() %>% rename(gdw_cnt=n) %>% rename(gene_name = gdw_gene_id)
  
  cnt_table = gup_cnt %>% merge(gdw_cnt,all=TRUE)
  cnt_table[is.na(cnt_table$gup_cnt),c("gup_cnt")]=0
  cnt_table[is.na(cnt_table$gdw_cnt),c("gdw_cnt")]=0
  cnt_table$total_cnt = cnt_table$gup_cnt + cnt_table$gdw_cnt
  
  if(!is.null(attr_name)) {
    cnt_table = cnt_table %>% select("gene_name","total_cnt")
    colnames(cnt_table) = c("gene_name",attr_name)
  }
  return(cnt_table)
}

## External resources
get_cosmic_genes = function() {
  cosmic = read.table(cosmic_path,header=T,sep="\t")
  
  cosmic_genes = strsplit(paste(cosmic$Gene.Symbol,cosmic$Synonyms,sep=",",collapse = ","),",")[[1]]
  cosmic_genes = remove_version_from_id(cosmic_genes)
  
  oncogenes = strsplit(paste( dplyr::filter(cosmic,grepl("oncogene",Role.in.Cancer))$Gene.Symbol,
                              dplyr::filter(cosmic,grepl("oncogene",Role.in.Cancer))$Synonyms,sep=",",collapse = ","),",")[[1]]
  oncogenes = remove_version_from_id(oncogenes)
  
  
  tsg = strsplit(paste( dplyr::filter(cosmic,grepl("TSG",Role.in.Cancer))$Gene.Symbol,
                        dplyr::filter(cosmic,grepl("TSG",Role.in.Cancer))$Synonyms,sep=",",collapse = ","),",")[[1]]
  tsg = remove_version_from_id(tsg)
  
  cosmic_genes = as.data.frame(x = cosmic_genes)
  colnames(cosmic_genes)=c("gene_id")
  cosmic_genes = cosmic_genes %>%  dplyr::mutate(oncogene = gene_id %in% oncogenes, tsg = gene_id %in% tsg)
  return(cosmic_genes)
}                          
get_cancer_genes = function() {
  cosmic_genes=get_cosmic_genes()
  cosmic_genes$source="cosmic"
  
  cancer_genes = cosmic_genes #gene id, oncogene, tsg
  
  oncokb = read.table(oncokb_path,header=T,sep="\t")
  oncokb = oncokb[,c("gene_name","oncogene","tsg")] %>% 
    mutate(oncogene=ifelse(oncogene=="Yes",T,F), tsg=ifelse(tsg=="Yes",T,F),source="oncokb") %>% dplyr::rename(gene_id=gene_name)
  ## TODO: ensembl tx to gene name, annotate databases that it occurs in?
  # names(oncokb)
  
  cancer_genes=rbind(cancer_genes,oncokb)
  
  grobner_recurrent = read.table(grobner_recurrent_path,header=T,sep="\t") 
  grobner_recurrent=grobner_recurrent %>% dplyr::mutate(gene_id = gene, source="grobner")
  
  grobner_onco = grobner_recurrent[grobner_recurrent$alteration_type=="Amplification",]
  grobner_onco$oncogene=T
  grobner_onco$tsg=NA 
  
  
  grobner_tsg = grobner_recurrent[grobner_recurrent$alteration_type=="Deletion" | 
                                    grobner_recurrent$alteration_type=="Gene-disrupting structural variant" ,]
  grobner_tsg$oncogene=NA
  grobner_tsg$tsg=T
  
  cancer_genes=rbind(cancer_genes, grobner_onco[,names(cancer_genes)], grobner_tsg[,names(cancer_genes)])
  
  
  return(cancer_genes)
}

annotate_variant_class_fractions = function(cohort_report) {
  cohort_report = cohort_report %>% mutate(somatic_variant_only=NA,germline_variant_only=NA,low_af_only=NA,ambiguous=NA)
  
  cohort_report = cohort_report %>% mutate(somatic_variant_only = (somatic_variant & !patient_fusion %in%
                                                                     filter(cohort_report,germline_variant|low_af)$patient_fusion),
                                           germline_variant_only = (germline_variant & !patient_fusion %in% filter(cohort_report,somatic_variant|low_af)$patient_fusion),
                                           low_af_only = (low_af & !patient_fusion %in% filter(cohort_report,somatic_variant|germline_variant)$patient_fusion))
  
  #ambiguous needs to be done afterwards
  cohort_report = cohort_report %>% mutate(ambiguous = (!patient_fusion %in% filter(cohort_report,somatic_variant_only|germline_variant_only|low_af_only)$patient_fusion))
  
  if(uq_fusions(cohort_report) == sum(uq_fusions(filter(cohort_report,somatic_variant_only)), uq_fusions(filter(cohort_report,germline_variant_only)), uq_fusions(filter(cohort_report,low_af_only)), uq_fusions(filter(cohort_report,ambiguous)))
  ){
  #  print("PASSED sanity check")
  } else {
    print("WARNING sanity check")
    
  }
  return(cohort_report)
}

## from cohort report (or part of) select best matching SV 
# used for high confidence fusions SV typing

annotate_labels_cancer_common = function(df) { 
  required_cols =  c("anno_sv_population","anno_healthy_chimera","anno_cancer_chimera","anno_cancer_gene_db")
  if(length(required_cols[required_cols %in% names(df)]) != length(required_cols)) {
    print(paste0("WARNING: Missing columns necessary for annotation label: ", required_cols[!required_cols %in% names(df)]))
    return(df)
    
  }
  
  df$annotation = ""
  df$annotation = as.character(df$annotation)
  df = df %>% mutate(annotation = ifelse( ((anno_sv_population|anno_healthy_chimera) &
                                             (anno_cancer_chimera|anno_cancer_gene_db)),"both",
                                          ifelse( (anno_sv_population|anno_healthy_chimera),"common",
                                                  ifelse( (anno_cancer_chimera|anno_cancer_gene_db) ,"cancer",""))))
  
  if("clinically_validated" %in% names(df)) {
    df = df %>% mutate(annotation = ifelse(clinically_validated,"clinical",annotation))
  }
  
  df$annotation = factor(df$annotation)
  
  return(df)
}


annotate_labels_variant_type = function(df) {
  df = df %>% mutate(variant_type = ifelse(somatic_variant,"tumor_specific",ifelse(germline_variant,"germline",
                                                                            ifelse(low_af,"low_af","ambiguous"))))
  df$variant_type = factor(df$variant_type)
  return(df)
}

annotate_labels_genes = function(cohort_report) {
  cohort_report$gup_label = ""
  cohort_report$gdw_label = ""
  cohort_report = cohort_report %>% 
    mutate(gup_label = ifelse(anno_gup_oncogene,paste0("oncogene,",gup_label),gup_label)) %>%
    mutate(gup_label = ifelse(anno_gup_tsg,paste0("tsg,",gup_label),gup_label)) %>%
    mutate(gup_label = ifelse(anno_gup_kinase,paste0("kinase,",gup_label),gup_label)) %>%
    mutate(gdw_label = ifelse(anno_gdw_oncogene,paste0("oncogene,",gdw_label),gdw_label)) %>%
    mutate(gdw_label = ifelse(anno_gdw_tsg,paste0("tsg,",gdw_label),gdw_label)) %>%
    mutate(gdw_label = ifelse(anno_gdw_kinase,paste0("kinase,",gdw_label),gdw_label))
  return(cohort_report)
}


make_uq_patient_fusion_df  = function(cohort_report) {
  if(nrow(cohort_report)==0) {return(data.frame())}
  uq_fusions_df = cohort_report
  
  #if multiple then select highest tumor-af
  #could be multiple sv types but all by multiple tools
  
  ## Chose for selecting instead of merging AF/length because there was a reason they were not merged before
  # its only for the typing and counting

  duplicate_fusions = uq_fusions_df[duplicated(uq_fusions_df$patient_fusion),]
  duplicate_fusions = uq_fusions_df %>% filter(patient_fusion %in% duplicate_fusions$patient_fusion) %>% 
    group_by(patient_fusion) %>% summarize(max_tumor_af = max(tumor_af_mean,na.rm=T))
  duplicate_fusions_keep = uq_fusions_df %>% merge(duplicate_fusions,by.x=c("patient_fusion","tumor_af_mean"),by.y=c("patient_fusion","max_tumor_af"))
  
  uq_fusions_df = uq_fusions_df %>% filter(!patient_fusion %in% duplicate_fusions$patient_fusion)
  uq_fusions_df=rbind(duplicate_fusions_keep,uq_fusions_df)
  
  ## check if unique
  if(nrow(uq_fusions_df[duplicated(uq_fusions_df$patient_fusion),])!=0) {
    print("WARNING duplicates not removed")
    #print(uq_fusions_df[duplicated(uq_fusions_df$fusion_name),])
  }
  
  #reannotate, ambiguous should be removed now
  uq_fusions_df = annotate_variant_class_fractions(uq_fusions_df)
  
  if(nrow(filter(uq_fusions_df,ambiguous))!=0) {  
    print("WARNING ambiguous not removed")
    #print(uq_fusions_df[duplicated(uq_fusions_df$fusion_name),])
  }
  
  
  #after selecting higest tumor AF only one patient-fusion-sv remains
  
  ## Apply labels
  uq_fusions_df = annotate_labels_variant_type(uq_fusions_df)
  uq_fusions_df = annotate_labels_cancer_common(uq_fusions_df)
  
  ##Adjust sv type if complex
  uq_fusions_df = uq_fusions_df %>% mutate(svtype_label = ifelse(grepl(", ",svtype),"complex",svtype))
  uq_fusions_df$svtype_label = factor(uq_fusions_df$svtype_label)
  
  return(uq_fusions_df)
}


## from the unique fusions go to unique gene pairs
make_uq_gene_pairs = function(uq_fusions_df) {
  
  uq_gene_pairs_df = uq_fusions_df  %>% 
    group_by(fusion_name,somatic_variant,germline_variant,low_af) %>% 
    summarize(svtype = toString(unique(sort(unlist(strsplit(svtype,", "))))),
              tumor_af = mean(tumor_af_mean,na.rm = T),
              normal_af = mean(normal_af_mean,na.rm = T),
              svlen = mean(svlen,na.rm=T),
              patient_ids=toString(unique(sort(patient_id))),
              patient_cnt=length(unique(patient_id)),
              .groups="keep")
  
  ## duplicate fusions are labelled as ambiguous
  duplicate_fusions = uq_gene_pairs_df[duplicated(uq_gene_pairs_df$fusion_name),]
  duplicate_fusions = uq_fusions_df %>% filter(fusion_name %in% duplicate_fusions$fusion_name) %>% 
    group_by(fusion_name) %>% 
    summarize(svtype = toString(unique(sort(unlist(strsplit(svtype,", "))))),
              tumor_af = mean(tumor_af_mean,na.rm = T),
              normal_af = mean(normal_af_mean,na.rm = T),
              svlen = mean(svlen,na.rm=T),
              patient_ids=toString(unique(sort(patient_id))),
              patient_cnt=length(unique(patient_id)),
              .groups="keep")
  
  duplicate_fusions = duplicate_fusions %>% mutate(somatic_variant=F,germline_variant=F,low_af=F)
  uq_gene_pairs_df = uq_gene_pairs_df %>% filter(!fusion_name %in% duplicate_fusions$fusion_name)
  uq_gene_pairs_df=rbind(duplicate_fusions,uq_gene_pairs_df)
  
  if(nrow(uq_gene_pairs_df[duplicated(uq_gene_pairs_df$fusion_name),])!=0) {
    print("WARNING duplicates not removed")
    print(uq_gene_pairs_df[duplicated(uq_gene_pairs_df$fusion_name),])
  }
  
  
  ## Apply labels
  uq_gene_pairs_df = annotate_labels_variant_type(uq_gene_pairs_df)
                                                
  ## add annotation after grouping
  fusion_name_labels = uq_fusions_df %>% select(fusion_name, anno_healthy_chimera,anno_cancer_chimera,anno_has_onco_or_tsg,anno_clinically_relevant,
                                                anno_has_oncogene,anno_has_tsg,anno_has_kinase,anno_cancer_gene_db) %>% unique()
  uq_gene_pairs_df = uq_gene_pairs_df %>% left_join(fusion_name_labels,by="fusion_name")
  
  #add sv level annotation
  fusions_svs_map = uq_fusions_df %>% select(fusion_name,patient_fusion_sv,anno_sv_population)
  uq_gene_pairs_df = uq_gene_pairs_df %>% mutate(anno_sv_population = fusion_name %in% 
                                                   filter(fusions_svs_map,anno_sv_population)$fusion_name)

  #will not add clinically validated annotation since not present
  uq_gene_pairs_df = annotate_labels_cancer_common(uq_gene_pairs_df)
  
  ##Adjust sv type if complex
  uq_gene_pairs_df = uq_gene_pairs_df %>% mutate(svtype_label = ifelse(grepl(", ",svtype),"complex",svtype))
  uq_gene_pairs_df$svtype_label = factor(uq_gene_pairs_df$svtype_label)
  
  return(uq_gene_pairs_df)
}

make_gene_exp_long = function(cohort_report) {
  if(!"gup_label" %in% names(cohort_report)) {
    cohort_report = annotate_labels_genes(cohort_report)
  }
  genes_up = cohort_report  %>%
    select(fusion_name,patient_id,patient_fusion,patient_fusion_sv,patient_clinrel_fusion,svtype,gup_gene_id,svtype,
           gup_fpkm,gup_fpkm_zscore,gup_fpkm_pval,gup_fpkm_log_normal_dist,gup_fpkm_zscore_group,
           gup_fpkm_zscore_supergroup,gup_fpkm_log_normal_dist_supergroup,gup_fpkm_supergroup_pval,
           gup_copy_ratio_l2fc_mean,gup_start_copy_ratio_l2fc,gup_label,gup_cytoband) %>% unique() %>%
    dplyr::rename(fpkm=gup_fpkm,
                  fpkm_zscore = gup_fpkm_zscore,
                  fpkm_zscore_dist=gup_fpkm_log_normal_dist,
                  fpkm_pval=gup_fpkm_pval,
                  fpkm_zscore_group = gup_fpkm_zscore_group,
                  fpkm_zscore_supergroup = gup_fpkm_zscore_supergroup,
                  fpkm_zscore_supergroup_dist=gup_fpkm_log_normal_dist_supergroup,
                  fpkm_supergroup_pval=gup_fpkm_supergroup_pval,
                  cna_l2fc=gup_copy_ratio_l2fc_mean,
                  cna_l2fc_bp=gup_start_copy_ratio_l2fc,
                  label=gup_label,
                  cytoband=gup_cytoband,
                  gene_name=gup_gene_id)
  
  genes_up$side="upstream"
  
  genes_down = cohort_report  %>%
    select(fusion_name,patient_id,patient_fusion,patient_fusion_sv,patient_clinrel_fusion,gdw_gene_id,svtype,
           gdw_fpkm,gdw_fpkm_zscore,gdw_fpkm_pval,gdw_fpkm_log_normal_dist,gdw_fpkm_zscore_group,
           gdw_fpkm_zscore_supergroup,gdw_fpkm_log_normal_dist_supergroup,gdw_fpkm_supergroup_pval,
           gdw_copy_ratio_l2fc_mean,gdw_end_copy_ratio_l2fc,gdw_label,gdw_cytoband) %>% unique() %>%
    dplyr::rename(fpkm=gdw_fpkm,
                  fpkm_zscore = gdw_fpkm_zscore,
                  fpkm_zscore_dist=gdw_fpkm_log_normal_dist,
                  fpkm_pval=gdw_fpkm_pval,
                  fpkm_zscore_group = gdw_fpkm_zscore_group,
                  fpkm_zscore_supergroup = gdw_fpkm_zscore_supergroup,
                  fpkm_zscore_supergroup_dist=gdw_fpkm_log_normal_dist_supergroup,
                  fpkm_supergroup_pval=gdw_fpkm_supergroup_pval,
                  cna_l2fc=gdw_copy_ratio_l2fc_mean,
                  cna_l2fc_bp=gdw_end_copy_ratio_l2fc,
                  label=gdw_label,
                  cytoband=gdw_cytoband,
                  gene_name=gdw_gene_id)
  
  genes_down$side="downstream"
  
  genes_long = rbind(genes_up,genes_down)
  genes_long$side = factor(genes_long$side,levels=c("upstream","downstream"))
  genes_long[genes_long$label=="",c("label")]=NA
  genes_long$label = substr(genes_long$label,1,nchar(genes_long$label)-1)
  genes_long$label = factor(genes_long$label,levels = c("oncogene","tsg","tsg,oncogene","kinase","kinase,tsg","kinase,oncogene","kinase,tsg,oncogene"))
  
  return(genes_long)  
}


