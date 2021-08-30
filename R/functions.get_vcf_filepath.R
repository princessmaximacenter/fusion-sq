
## Adjust to suit you own file naming convention
get_vcf_filepath = function(identifier,tool,somatic=FALSE) {
  #delly needs _WGS in the back to prevent grepping somatic too 
  if(!somatic) {
    #germline
    filepath = switch(tool,
                      manta={paste0(manta_dir,"*",identifier,"*_WGS.diploidSV*.vcf.gz")},
                      delly={paste0(delly_dir,"*",identifier,"*_WGS.delly.*.vcf.gz")},
                      gridss={paste0(gridss_dir,"*",identifier,"*_WGS.gridss.*.vcf.gz")}
    )
  } else {
    if(tool == "manta") {
      filepath = paste0(manta_dir,identifier,"*.somaticSV*.vcf.gz")
    } 
    else {
      ##for Delly and GRIDSS: no somatic in use
      filepath =""
      return()
    }
  }
  
  if(filepath==""){return()}
  
  vcf_file = Sys.glob(filepath)
  if(length(vcf_file)<1){
    print(paste0("File missing: ",filepath))
    return("")
  } else if(length(vcf_file)>1) {
    print(print(paste0("Multiple files: ",vcf_file)))
    vcf_file=vcf_file[1]
  }
  return(vcf_file)  
}
