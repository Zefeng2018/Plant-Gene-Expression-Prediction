# ChIP-seq analysis workflow

A general framework to analysis the ChIP-Seq data.It mainly includes the download of sra data, fastq file extracting, reads quality filtering, reads mapping, uniqe mapped reads extract, PCR duplicated reads removal, peak calling (optional).

## 0. Preparation 
    mkdir work_directory
    cd work_directory
    mkdir 1fastq 2fastqc 3filtered_fq 4sam 5uniq_sam 6sorted_bam 7dedup_bam 8macs2

## 1. Download data from sra list FILE
    awk 'type=substr($0,1,3){pre_id=substr($0,1,6)}split($1,m,"."){print "/sra/sra-instant/reads/ByRun/sra/"type"/"pre_id"/"m[1]"/"$1".sra"}' SraAccList.txt > download.txt 
    for m in $(cat download.txt); do if [ ! -e $(basename $m) ]; then echo $m; ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -k 1 -QT -l 200m anonftp@ftp-trace.ncbi.nlm.nih.gov:$m .; fi; done

NOTE: the download.txt file should contain the sra id list, with each sra id is a sperate line, e.g.,
SRR1509474
SRR1509472
SRR1509477
SRR2637923
SRR2637924
SRR2637925
SRR1848405

## 2. Extract fastq files from sra data
    for sra in $(ls *.sra); do fastq-dump $m --split-files -O 1fastq;done   #used for paired-end sequenced data
    for sra in $(ls *.sra); do fastq-dump $m -O 1fastq;done                 #used for single-end sequenced data

## 3. Fastqc inspecting
    cd 1fastq
    for m in $(ls *.fastq); do echo $m; fastqc $m  -o ../2fastqc; done

## 4. Reads filter by Trimmomatic tool.
### paired-end reads
    cd ../1fastq 
    for m in $(ls *.fastq ); do java -jar ~/mysoft/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 7 $m_1.fastq $m_2.fastq  ../3filtered_fq/$m_1.paired.fastq ../3filtered_fq/$m_1.unpaired.fastq ../3filtered_fq/$m_2.paired.fastq ../3filtered_fq/$m_2.unpaired.fastq ILLUMINACLIP:~/mysoft/Trimmomatic-0.36/adapters/my_adapter.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:100
    
### (or) single-end reads
    for m in $(ls *.fastq ); do java -jar ~/mysoft/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 7 $m  ../3filtered_fq/$m_filtered.fastq  ILLUMINACLIP:/home/wuzefeng/mysoft/Trimmomatic-0.36/adapters/TruSeq9-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:100

NOTE: one may select an alternative tool in this step.

## 5. Bowtie2 mapping
    cd ../3filtered_fq/
### single-end reads
    for m in $(ls *.filtered.fastq); do bowtie2 -x ~/genome.index -1 $m_1.filtered.fastq -2  $m_2.filtered.fastq  -S ../4sam/$m.sam --no-mixed  --no-discordant --no-unal -p 7 ;done
### (or) paired-end reads
    for m in $(ls *.filtered.fastq); do bowtie2 -x ~/gebnome.index -U $m  -S ../4sam/$m.sam  -p 7  --no-unal; done

## 6. Get unique mapped reads
    cd ../4sam/
    for m in $(ls *.sam); do samtools view -Sh $m | grep -e "^@" -e "XM:i:[012][^0-9]" | grep -v "XS:i:" > ../5uniq_sam/$m.uniq.sam;done

## 7. Sorting alignemnt (.sam) file and convert to (.bam) files
    cd ../5uniq_sam
    for  m in $(ls *.sam); do  picard-tools SortSam  INPUT=$m OUTPUT=../6sorted_bam/$m.sorted.bam SO=coordinate; done

## 8. Remove PCR duplicated reads 
    cd ../6sorted_bam
    for m in $(ls *.bam); do samtools rmdup -s $m ../7dedup_bam/$m.dedup.bam ; done

## 9. Peak calling using MACS2 tool (optional)
    cd ../7dedup_bam
    macs2 callpeak -t H3K4me3.bam -c input.bam -f BAM -g 1.2e8 -n H3K4me3 -B --SPMR -q  0.01 --outdir ../8macs2/narrow/  # narrow peak
    macs2 callpeak -t H3K4me3.bam -c input.bam -f BAM -g 1.2e8 -n H3K4me3.bam -B --SPMR  --broad --broad-cutoff 0.1  --outdir ../8macs2/broad  # broad peaks


### Vistor statistics

<script type="text/javascript" src="//ra.revolvermaps.com/0/0/7.js?i=0ypfp1eocyh&amp;m=0&amp;c=ff0000&amp;cr1=ffffff&amp;sx=0" async="async"></script>


