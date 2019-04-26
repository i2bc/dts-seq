# Dts-seq

**Computanional workflow for Dts-seq analysis**

You will find here the computational protocole for the analysis of the Dts-seq data related to the publication [here](https:??).

**Contact**

- Ji Wang (<??>)
- Forence Lorieux (<??>)
- Claire Toffano-Nioche (<claire.toffano-nioche@u-psud.fr>)
- Daniel Gautheret (<??>)
- Jean Lehmann (<??>)

## RNAseq: from fastq to reads coverage

### Repository architechture 

* Protocole: creating a repository for the analysis: 
* Code:
```bash
mkdir Dts-seq ; 
cd Dts-seq ; 
mkdir 1_rawData 2_mapping 3_counting 4_?? 5_?? 6_tRNA_modification
```
* Result: the architechture of `Dts-seq` repository

### Data

#### Genome & annotations 

* Protocole: downloaded from the ncbi, acces: GCF_000005845.2_ASM584v2
* Code : 
```bash
cd 1_rawData
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz
gunzip GCF_000005845.2_ASM584v2_genomic.*.gz
```
* Result files: 
  * GCF_000005845.2_ASM584v2_genomic.fna
  * GCF_000005845.2_ASM584v2_genomic.gff

#### RNAseq Data 

- Protocole: R2 files were downloaded from ENA (acces: ??) into the local `Dts-seq/1_rawData` repository.
- Code to list all samples files:
```bash
rep in A B C ; do for samples in 3-D 4-NT 5-nD ; do ls ${rep}${samples}_R2_fastq.gz ; done ; done
```
- Result files: 15 `*_R2_fastq.gz`

|read number |     A5-nD, B5-nD, C5-nD      |       A4-NT, B4-NT, C4-NT    |       A3-D, B3-D, C3-D       |
|------------|:----------------------------:|:----------------------------:|:----------------------------:|
|QC-passed   | 11796813, 13114509, 11455875 | 11680073, 13517533, 10775965 | 15388902, 12562074, 10691188 |

### Cleanning step

- Protocol: each raw fastq file was cleaned (cutadapt software) with a specific polyA adapter following the wet protocole and reads shorter than 10 bp after polyA trimming were discarded. Quality control were done (FastQC software). 
- Code:
```bash
nohup bash -c 'for i in *_R2.fastq.gz ; do SampleName=`basename -s .fastq.gz ${i}` ; docker run -v /root/mydisk/data/JW:/data:rw -w /data docker-registry.genouest.org/bioconda/cutadapt cutadapt --adapter=AAAAAAAAAAAAAAA --minimum-length=10 --output=${SampleName}_noPolyA.fastq ${i} ; mv nohup.out cutadapt_${SampleName}.log ; done' &
nohup bash -c 'for i in *_noPolyA.fastq ; do SampleName=`basename -s _noPolyA.fastq ${i}` ; docker run -v /root/mydisk/data/JW:/data:rw -w /data docker-registry.genouest.org/ifb/fastqc fastqc -o FastQC ${i} 2> fastqc_${SampleName}.log ; done' &
cd ..
```
- Result files:
  - 15 output fastq files (`*_noPolyA.fastq`) 
  - 15 log files (`cutadapt_*.log`)
  - 15 repositories containing the quality control results (`*.html`) 

|read number |     A5-nD, B5-nD, C5-nD      |       A4-NT, B4-NT, C4-NT    |       A3-D, B3-D, C3-D       | 
|------------|:----------------------------:|:----------------------------:|:----------------------------:|
|trimmed     |  9736417,  8639091,  9522279 | 10257977, 10586494, 10181006 | 12595828, 10640726,  9146981 | 

### Mapping step

- Protocol: Reads mapping was done with bowtie2 with the local mapping option to maximise the alignment length.
- Code:
```bash
# create genomic index file for bowtie2
cd 2_mapping
docker run -v /root/mydisk/data/JW:/data:rw -w /data docker-registry.genouest.org/bioconda/bowtie2 bowtie2-build ../1_rawData/GCF_000005845.2_ASM584v2_genomic.fna index_bt2x/NC_00913
# bowtie2 run
rm nohup.out ; nohup bash -c 'for i in ../1_rawData/*_noPolyA.fastq ; do sample=`basename $i _noPolyA.fastq` ; echo "------------" $sample "-------------" ; docker run -v Dts-seq/2_mapping/:/data:rw -w /data docker-registry.genouest.org/bioconda/bowtie2 bowtie2 -x index_bt2x/NC_00913 --phred33 --local $i > ${sample}.sam ; done '
```
- Result files (into `2_mapping` repository):
  - 6 index files for bowtie2 (`*.btz2` into `index_bt2x` repository)
  - 15 *.sam

|read number |     A5-nD, B5-nD, C5-nD      |       A4-NT, B4-NT, C4-NT    |       A3-D, B3-D, C3-D       |
|------------|:----------------------------:|:----------------------------:|:----------------------------:|
|mapped      |  7209763,  5765983,  7528570 |  7829777,  7574852,  8282069 |  9771220,  8149220,  7276564 | 

### Alignments selection step

- Protocol: selection of mapped reads on the genomic sequence (with samtools view) and presenting a complete 3' end, ie. either CCA or TGG depending on the DNA strand (awk). The resulting alignment files were sorted by increasing locations and the associated index files for binary management were created.
- Code:
```bash
rep in A B C ; do for s in "3-D" "4-N"T "5-nD" ; do
   samtools view -h 2_mapping/${rep}${s}.sam NC_000913.3 | awk 'BEGIN{FS="\t";OFS="\t"}{if( ($0~/@/)||((($2==16)&&($10~/CCA$/))||(($2==0)&&($10~/^TGG/))) ){print $0}}' | samtools view -hb -o 2_mapping/${rep}${s}_CCATGG_unsort.bam - ; 
   samtools sort 2_mapping/${rep}${s}_CCATGG_unsort.bam > 2_mapping/${rep}${s}_CCATGG.bam 
   samtools index 2_mapping/${rep}${s}_CCATGG.bam ; 
done ; done ;
```
- Result files (into `2_mapping` repository): 
  - 15 *_CCATGG.bam 
  - 15 *_CCATGG.bam.bai 

|read number |     A5-nD, B5-nD, C5-nD      |       A4-NT, B4-NT, C4-NT    |       A3-D, B3-D, C3-D       | 
|------------|:----------------------------:|:----------------------------:|:----------------------------:|
|remainning  |  6494515,  5150196,  6941871 |  7229758,  6803833,  7551225 |  9045057,  7457325,  6771946 | 

### Reads coverage computation step

- Protocole: creation of coverage files (both format wig and 2 columns) with strand separation. As alignments came from R2 paire, exchange of reverse and forward strands (join).
- Code:
```bash
for i in 2_mapping/*_CCATGG.bam ; do
   sample=`basename $i .bam` ;
   # coverage for reverse strand
   samtools view -h -b -f 16 ${sample}.bam NC_000913.3 | samtools depth -m 10000000 -a - > ${sample}_depth_rev.txt ; 
   # coverage for forward strand
   samtools view -h -b -F 0x14 ${sample}.bam NC_000913.3 | samtools depth -m 10000000 -a - > ${sample}_depth_for.txt ; 
   # strands association 
   join -t $'\t' -12 -22 -o 1.3,2.3 ${sample}_depth_rev.txt ${sample}_depth_for.txt > ${sample}_depth_fr.txt ; 
   # coverage: log computation => *_covlog.wig
   awk 'BEGIN{FS="\t";print "variableStep chrom=NC_000913.3"}{if($1==0){logF=0}else{logF=log($1)/log(2)}; if($2==0){logR=0}else{logR=log($2)/log(2)};printf "%d %2.2f\n%d -%2.2f\n",NR,logF,NR,logR}' ${sample}_depth_fr.txt > ${sample}_covlog.wig ; 
done ;
```
- Result files:
  - 45 *_depth_*.txt (raw count)
  - 15 *_covlog.wig for (igv visualisation of log2 coverage) ?? not use ??

### count in 3prime of tRNA ??

```bash
awk 'BEGIN{FS="\t"}{if($3=="tRNA"){if($7~"+"){posEnd=$5}else{posEnd=$4};print "NC_000913.3:"posEnd"-"posEnd}}' NC_000913.gff > NC_000913_tRNA_3prime.list
for rep in "A" "B" "C" ; do for sample in "3-D" "4-NT" "5-nD" ; do 
   rm ${sample}_m${map}_tRNA_3prime_count.txt ; 
   for t in `more NC_000913_tRNA_3prime.list` ; do 
      samtools view 2_mapping/${rep}${sample}_CCATGG.bam ${t} | wc -l >> ${rep}${sample}_tRNA_3prime_count.txt ;
   done ; 
done ; done
```

## Modified bases: from modomics DB to gff

### Data

- tRNA with modified bases
get from modomics, acces: ?? (copy/paste in text format)

### location in genomic the sequence

blast
blast analysis
automatisation

## Terminason signal: from read coverage to ts-jump

## Tools version

- fastqc: from docker-registry.genouest.org/ifb/fastqc, version ??
- cutadapt: from docker-registry.genouest.org/bioconda/cutadapt, version 1.11
- bowtie2: from docker-registry.genouest.org/bioconda/bowtie2, bowtie2-align-s version 2.2.8 
- samtools: version 1.4

