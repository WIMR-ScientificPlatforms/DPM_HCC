#this was done on per project spec. dont expect any of this to work

cd PATH/TO/ega_"$snam"
module load R
R
library(dplyr)
dir.create("freebayes/unique")
tmp<-read.table("source/sample_table.txt",sep="\t",header=T) 

sdeet<-tmp %>% mutate(id=donor_id,sample=sample_id, vcf=paste0("freebayes/",sample,"_accept1pc.vcf.gz"),cond=paste0(sample_type,".",sample_number),vcfExists=vcf %in% dir("freebayes",full.names=T))
 table(sdeet$vcfExists)
  nrow(sdeet)      
write.csv(sdeet,file="source/Intersection_dataframe.csv")

# filter for failed samples
sdeet<-subset(sdeet,id %in% c(3,60,66))

lapply(na.omit(unique(sdeet$id)), function(X) {

	message(X)
	md<-filter(sdeet,id==X)
	vcfs<-md$vcf
	stopifnot(length(vcfs)>1)
	pp<-expand.grid(vcfs,vcfs) %>% subset(!Var1==Var2)

	system(paste0("cat source/header.txt > ",paste0("source/int_and_anno_donor",X,".sh")))
	apply(pp,1, function(Y) {
	paste0(paste("bcftools isec -n~10 -c all",Y[1],Y[2],'| bcftools view -T -',Y[1],' -Oz > freebayes/unique/'),X,"_",md[md$vcf==Y[1],"cond"],"_v_",md[md$vcf==Y[2],"cond"],"_1pc.vcf.gz")
	}) %>% data.frame %>% write.table(file=paste0("source/int_and_anno_donor",X,".sh"),row.names=F,col.names=F,quote=F,append=T) 
	system(paste0('echo snam=',X, ' >> ', paste0("source/int_and_anno_donor",X,".sh")) )
	system(paste0("cat source/footer.txt >> ",paste0("source/int_and_anno_donor",X,".sh")))
	system(paste0("echo ",paste0("qsub -o logs/vep_",X,".log ",paste0("source/int_and_anno_donor",X,".sh"))," >> source/runvep_commands.sh"))
	})


# where source/header.txt is
#PBS flags

module load anaconda
conda activate /project/RDS-FMH-DPM_HCC-RW/WES/software/vep
module load tabix

wdir='PATH/TO/ega_"$snam"'
cd "$wdir"


# and  source/footer.txt is 
for i in `ls freebayes/unique/"$snam"_*gz` ; do echo "$i"

tabix -p vcf "$i"

onam=`echo "$i" | sed 's/freebayes\/unique\///g' | sed 's/.vcf.gz//g'`

bcftools filter --regions chr1,chr2,chr3,chr4,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM -i 'MIN(FMT/DP)>15 & QUAL>20' "$i" | vep --gencode_basic --sift b --polyphen b --af --af_1kg --af_gnomad --max_af --format vcf --flag_pick --no_intergenic --offline --dir_cache /project/RDS-FMH-DPM_HCC-RW/WES/annotation --buffer_size 500 --fork 3 -o vep/"$onam"_soma_dp15_qual20_vep_pic_tab.txt

bcftools filter --regions chr1,chr2,chr3,chr4,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM -i 'MIN(FMT/DP)>15 & QUAL>20' "$i" | vep --gencode_basic --sift b --polyphen b --af --af_1kg --af_gnomad --max_af --format vcf --flag_pick --no_intergenic --offline --dir_cache /project/RDS-FMH-DPM_HCC-RW/WES/annotation --buffer_size 500 --fork 3 --vcf -o vep/"$onam"_soma_dp15_qual20_vep_pic.vcf

done