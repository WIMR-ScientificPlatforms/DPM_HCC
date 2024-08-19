# create env
conda create --prefix PATH/TO/DIR/vep python=3.6 
conda activate PATH/TO/DIR/vep
conda install -c bioconda ensembl-vep
pip3 install update pyega3 -U

# get annotaion
cd PATH/TO/annotation
curl -O -C - http://ftp.ensembl.org/pub/release-104/variation/vep/homo_sapiens_vep_104_GRCh38.tar.gz
tar xzf homo_sapiens_vep_104_GRCh38.tar.gz

# dbgap annotations
cd PATH/TO/project/
module load sratoolkit
vdb-config --import downloads/key.ngc 
vdb-config -i

# get files request and cart from dbgap and cp path to aspc correctly, thne
cd downloads 
vdb-decrypt DOWNLOAD

# get metadata and cart files from sra, save to downloads
prefetch downloads/cart_DAR110735_202205200000.krt #going

# sra_download_script
cd PATH/TO/project/sra/
module load sratoolkit
fasterq-dump "$snam".sra -O cd PATH/TO/project/fastq/ -e 2
gzip cd PATH/TO/project//fastq/"$snam"*fastq

# run
for i in `cut -d "," -f 1 source/SraRunTable.txt` ; do echo "$i"_unpack_slow >> jerbs.txt ; qsub -v "snam="$i"" source/sra_download_script.sh | tee -a jobs.txt ; done


# EGA
module load anaconda
conda activate /project/RDS-FMH-DPM_HCC-RW/WES/software/vep

if [ ! -d /scratch/DPM_HCC/ega_"$snam" ]; then mkdir PATH/TO/ega_"$snam" ; fi
cd PATH/TO/ega_"$snam"
if [ ! -f project_files_"$snam".txt ]; then pyega3 -cf PATH/TO/credentials.json files "$snam" > project_files_"$snam".txt ; fi
pyega3 -c 56 -cf PATH/TO/credentials.json fetch "$snam"

cd PATH/TO/ega_"$snam"
mkdir {bams,downloads,fastq,freebayes,logs,source,vep}
mkdir fastq/trimgalore
mv EGAF*/*fastq.gz fastq

# run 
for i in `cat source/samples.txt` ; do echo "$i"_alignVar >> jerbs.txt 
    qsub -o logs/varcall_"$i".log -v "snam="$i"" source/varcall.sh | tee -a jerbs.txt 
done