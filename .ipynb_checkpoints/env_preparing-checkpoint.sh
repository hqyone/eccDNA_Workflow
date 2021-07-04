parent_dir="/home/hqyone/mnt/2t/eccDNA/genome"
genome_name="hg38_p13"

REF_FA_URL="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh38.p13.genome.fa.gz"
REF_GTF_URL="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.chr_patch_hapl_scaff.annotation.gtf.gz"

[ ! -d "$parent_dir/$genome_name"] && mkdir "$parent_dir/$genome_name"
if [ ! -f "$parent_dir/$genome_name/GRCh38.p13.genome.fa" ]; then
    wget -O $REF_FA_URL -P "$parent_dir/$genome_name"
    unzip -d "$parent_dir/$genome_name/GRCh38.p13.genome.fa.gz"
    samtools faidx "$parent_dir/$genome_name/GRCh38.p13.genome.fa" -@ 3
    bwa index "$parent_dir/$genome_name/GRCh38.p13.genome.fa"
fi

if [ ! -f "$parent_dir/$genome_name/gtf/gencode.v37.chr_patch_hapl_scaff.annotation.gtf" ]; then 
    [ ! -d "$parent_dir/$genome_name/gtf"] && mkdir "$parent_dir/$genome_name/gtf"
    wget -O $REF_GTF_URL -P "$parent_dir/$genome_name/gtf"
    unzip -d "$parent_dir/$genome_name/gtf/gencode.v37.chr_patch_hapl_scaff.annotation.gtf.gz"
    sort -k1,1 -k4,4n "$parent_dir/$genome_name/gtf/gencode.v37.chr_patch_hapl_scaff.annotation.gtf"> "$parent_dir/$genome_name/gtf/gencode.v37.all.sort.gtf"
    bgzip "$parent_dir/$genome_name/gtf/gencode.v37.all.sort.gtf"
    tabix -p gff "$parent_dir/$genome_name/gtf/gencode.v37.all.sort.gtf.gz"
