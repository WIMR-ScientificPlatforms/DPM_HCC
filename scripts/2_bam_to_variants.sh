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

if [ ! -d fastq/combined ]; then mkdir fastq/combined ; fi
if [ ! -d downloads/resorted ]; then mkdir downloads/resorted ; fi

samtools sort -@ 7 -m 7G -n -o downloads/resorted/"$snam"_resort.bam downloads/"$snam".bam 

samtools fastq -@ 8 -1 fastq/"$snam"_R1.fq -2 fastq/"$snam"_R2.fq downloads/resorted/"$snam"_resort.bam 


repair.sh in1=./fastq/"$snam"_R1.fq in2=./fastq/"$snam"_R2.fq out1=./fastq/fixed_"$snam"_R1.fq out2=./fastq/fixed_"$snam"_R2.fq outsingle=./fastq/fixed_"$snam"_singleton.fq

trim_galore --fastqc --paired --gzip -o fastq/trimgalore  fastq/fixed_"$snam"_R*.fq.gz 

#combine fastqs
cat fastq/"$snam"_*_1.clean.fq.gz > fastq/combined/"$snam"_R1.fq.gz &
cat fastq/"$snam"_*_2.clean.fq.gz > fastq/combined/"$snam"_R2.fq.gz 
wait

#multiyhread freebayes
bwa mem -t 8 -T 0 PATH/TO/annotation/hg38.fa fastq/combined/"$snam"_R1.fq.gz fastq/combined/"$snam"_R2.fq.gz |samtools view -@ 7 -Shb |  samtools sort -@ 7 -m 7G -n - | samtools fixmate -m -@ 7 - bams/"$snam"_nam_fixmate.bam

samtools sort -@ 7 -m 7G bams/"$snam"_nam_fixmate.bam | samtools markdup -@ 7 - bams/"$snam"_trim_sort_dup.bam #fine

freebayes -f PATH/TO/annotation/hg38.fa -F 0.01 -C 2 bams/"$snam"_trim_sort_dup.bam > freebayes/"$snam"_accept1pc.vcf

bgzip freebayes/"$fqdir"/"$snam"_accept1pc.vcf
tabix -p vcf freebayes/"$fqdir"/"$snam"_accept1pc.vcf.gz