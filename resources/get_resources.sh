#Get resources

# Gene annotation
wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/__genome_libs_StarFv1.8/GRCh38_gencode_v31_CTAT_lib_Oct012019.source.tar.gz
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.metadata.Annotation_remark.gz

## retrieve 
#gtf ref_annot_GRCh38_gencode_v31_CTAT_lib_Oct012019.gtf.gz
#sqlite ref_annot_GRCh38_gencode_v31_CTAT_lib_Oct012019.sqlite
#sjdb sjdb_GRCh38_gencode_v31_CTAT_lib_Oct012019.bed

#population database variants
wget https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd166.GRCh38.variant_call.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd166.GRCh38.variant_call.vcf.gz.tbi
wget https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd186.GRCh38.variant_call.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd186.GRCh38.variant_call.vcf.gz.tbi
wget https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd102.GRCh38.variant_call.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd102.GRCh38.variant_call.vcf.gz.tbi
wget http://dgv.tcag.ca/dgv/docs/GRCh38_hg38_variants_2020-02-25.txt -O dgv_GRCh38_hg38_variants_2020-02-25.tsv

