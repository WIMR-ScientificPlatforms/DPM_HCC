#PBS flags

module load bwa
module load samtools
module load picard
module load freebayes
module load trimgalore
module load anaconda
module load bbmap

conda activate /project/RDS-FMH-DPM_HCC-RW/WES/software/vep

cd PATH/TO/ega_"$snam"

if [ ! -d fastq/trimgalore ]; then mkdir fastq/trimgalore ; fi

trim_galore --fastqc --paired --gzip -o fastq/trimgalore  fastq/"$snam"*.fastq.gz 

bwa mem -t 8 -T 0 PATH/TO/annotation/hg38.fa fastq/trimgalore/"$snam"*_val_1.fq.gz fastq/trimgalore/"$snam"*_val_2.fq.gz |samtools view -@ 7 -Shb |  samtools sort -@ 7 -m 7G -n - | samtools fixmate -m -@ 7 - bams/"$snam"_nam_fixmate.bam

samtools sort -@ 7 -m 7G bams/"$snam"_nam_fixmate.bam | samtools markdup -@ 7 - bams/"$snam"_trim_sort_dup.bam #fine

freebayes -f PATH/TO/annotation/hg38.fa -F 0.01 -C 2 bams/"$snam"_trim_sort_dup.bam > freebayes/"$snam"_accept1pc.vcf

bgzip freebayes/"$fqdir"/"$snam"_accept1pc.vcf
tabix -p vcf freebayes/"$fqdir"/"$snam"_accept1pc.vcf.gz