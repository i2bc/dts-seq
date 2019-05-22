# Dts-seq

**Computational workflow for Dts-seq analysis**

You will find here the computational protocol for the analysis of Dts-seq data.

**Contact**

- Claire Toffano-Nioche (<claire.toffano-nioche@u-psud.fr>)
- Daniel Gautheret (<daniel.gautheret@u-psud.fr>)
- Jean Lehmann (<jean.lehmann@u-psud.fr>)

## RNAseq: from fastq to read coverage

### Repository architecture 

* Protocol: creating a repository for the analysis: 
* Code:
```bash
mkdir Dts-seq ; 
cd Dts-seq ; 
mkdir 1_rawData 2_processedData 3_mapping 4_selection 5_coverage 6_tRNA_modification
mkdir 2_processedData/FastQC 3_mapping/index_bt2x
```
* Result: the architecture of `Dts-seq` repository

### Data

#### Genome & annotations 

* Protocol: Genome downloaded from ncbi, accession: GCF_000005845.2_ASM584v2
* Code : 
```bash
cd 1_raw_data
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz
gunzip GCF_000005845.2_ASM584v2_genomic.*.gz
cd ..
```
* Result files: 
  * GCF_000005845.2_ASM584v2_genomic.fna
  * GCF_000005845.2_ASM584v2_genomic.gff

#### RNAseq Data 

- Protocol: R2 files were downloaded from ENA (access: ??) into the local `Dts-seq/1_rawData` repository.
- Code to list all samples files:
```bash
for rep in A B C ; do for samples in 3-D 4-NT 5-nD ; do ls 1_rawData/${rep}${samples}_R2.fastq.gz ; done ; done
```
- Result files: 9 `*_R2_fastq.gz`

|read number |     A5-nD, B5-nD, C5-nD      |       A4-NT, B4-NT, C4-NT    |       A3-D, B3-D, C3-D       |
|------------|:----------------------------:|:----------------------------:|:----------------------------:|
|QC-passed   | 11796813, 13114509, 11455875 | 11680073, 13517533, 10775965 | 15388902, 12562074, 10691188 |

### Cleaning 

- Protocol: each raw fastq file was cleaned (cutadapt) with a specific polyA adapter following the wet protocol and reads shorter than 10 bp after polyA trimming were discarded. Quality control was performed (FastQC software). 
- Code:
```bash
DTSSEQDIR="/path/to/Dts-seq/repository/from/the/root" ;
rm 2_processedData/cutadapt.log ; rm 2_processedData/fastqc.log ;
for i in `ls 1_rawData/*_R2.fastq.gz` ; do 
   SampleName=`basename -s _R2.fastq.gz ${i}` ; 
   docker run -v ${DTSSEQDIR}:/data:rw -w /data docker-registry.genouest.org/bioconda/cutadapt cutadapt --adapter=AAAAAAAAAAAAAAA --minimum-length=10 --output=2_processedData/${SampleName}_noPolyA.fastq ${i} >> 2_processedData/cutadapt.log 2>&1 ;
   docker run -v ${DTSSEQDIR}:/data:rw -w /data docker-registry.genouest.org/ifb/fastqc fastqc -o 2_processedData/FastQC 2_processedData/${SampleName}_noPolyA.fastq >> 2_processedData/fastqc.log 2>&1 ;
done'
cd ..
```
- Result files (into `2_processedData` repository):
  - 9 output fastq files (`*_noPolyA.fastq`) 
  - 2 log files (`cutadapt.log`, `fastqc.log`)
  - 9 repositories containing the quality control results (`*.html`) 

|read number |     A5-nD, B5-nD, C5-nD      |       A4-NT, B4-NT, C4-NT    |       A3-D, B3-D, C3-D       | 
|------------|:----------------------------:|:----------------------------:|:----------------------------:|
|trimmed     |  9736417,  8639091,  9522279 | 10257977, 10586494, 10181006 | 12595828, 10640726,  9146981 | 

### Mapping step

- Protocol: Reads mapping was done with bowtie2 with the local mapping option to maximise the alignment length.
- Code:
```bash
DTSSEQDIR="/path/to/Dts-seq/repository/from/the/root"
# create genome index file for bowtie2
rm 3_mapping/bowtie2-build.log ; rm 3_mapping/bowtie2-align.log ;
docker run -v ${DTSSEQDIR}:/data:rw -w /data docker-registry.genouest.org/bioconda/bowtie2 bowtie2-build 1_rawData/GCF_000005845.2_ASM584v2_genomic.fna 3_mapping/index_bt2x/NC_00913 >> 3_mapping/bowtie2-build.log 2>&1
# bowtie2 run
for i in 2_processedData/*_noPolyA.fastq ; do 
   sample=`basename $i _noPolyA.fastq` ; 
   docker run -v ${DTSSEQDIR}:/data:rw -w /data docker-registry.genouest.org/bioconda/bowtie2 bowtie2 -x 3_mapping/index_bt2x/NC_00913 --phred33 --local $i > 3_mapping/${sample}.sam 2>> 3_mapping/bowtie2-align.log ; 
done '
```
- Result files (into `2_mapping` repository):
  - 6 index files for bowtie2 (`*.btz2` into `index_bt2x` repository)
  - 9 *.sam
  - log files: `bowtie2-align.log`, `bowtie2-build.log`

|read number |     A5-nD, B5-nD, C5-nD      |       A4-NT, B4-NT, C4-NT    |       A3-D, B3-D, C3-D       |
|------------|:----------------------------:|:----------------------------:|:----------------------------:|
<--|mapped chr+spike     |  7209763,  5765983,  7528570 |  7829777,  7574852,  8282069 |  9771220,  8149220,  7276564 | -->
|mapped |   7091514, 5644988, 7434767  |   7760868, 7436857, 8097497  |   9705213, 8024150, 7164553  |



### Alignment selection

- Protocol: selection of mapped reads on the genomic sequence and presenting a complete 3' end, ie. either CCA or TGG depending on the DNA strand (awk). The resulting alignment files were sorted by increasing locations and the associated index files for binary management were created.
- Code:
```bash
rm 4_selection/selection.log
for rep in A B C ; do for s in "3-D" "4-N"T "5-nD" ; do
   awk 'BEGIN{FS="\t";OFS="\t"}{if( ($0~/@/)||((($2==16)&&($10~/CCA$/))||(($2==0)&&($10~/^TGG/))) ){print $0}}' 3_mapping/${rep}${s}.sam | samtools view -hu - | samtools sort - > 4_selection/${rep}${s}_CCATGG.bam 2>> 4_selection/selection.log ; 
   samtools index 4_selection/${rep}${s}_CCATGG.bam 2>> 4_selection/selection.log ;
done ; done ;
```
- Result files (into `3_mapping` repository): 
  - 9 *_CCATGG.bam 
  - 9 *_CCATGG.bam.bai 

|read number |     A5-nD, B5-nD, C5-nD      |       A4-NT, B4-NT, C4-NT    |       A3-D, B3-D, C3-D       | 
|------------|:----------------------------:|:----------------------------:|:----------------------------:|
|remaining   |  6494515,  5150196,  6941871 |  7229758,  6803833,  7551225 |  9045057,  7457325,  6771946 | 

### Read coverage computation

- Protocol: creation of coverage files (both format wig and 2 columns) with strand separation. As alignments came from R2 reads, exchange of reverse and forward strands (join).
- Code:
```bash
for i in 4_selection/*_CCATGG.bam ; do 
   sample=`basename $i .bam` ; 
   # coverage of reverse strand
   samtools view -h -b -f 16 ${i} NC_000913.3 | samtools depth -d 10000000 -a - > 5_coverage/${sample}_depth_rev.txt ; 
   # coverage for forward strand
   samtools view -h -b -F 0x14 ${i} NC_000913.3 | samtools depth -d 10000000 -a - > 5_coverage/${sample}_depth_for.txt ; 
   # strands association 
   join -t $'\t' -12 -22 -o 1.3,2.3 5_coverage/${sample}_depth_rev.txt 5_coverage/${sample}_depth_for.txt > 5_coverage/${sample}_depth_fr.txt ;
done 
```
- Result files (into `5_coverage` repository):
  - 9 *_depth_for.txt (coverage on forward strand)
  - 9 *_depth_rev.txt (coverage on reverse strand)
  - 9 *_depth_fr.txt (coverage on both strands)

### read count in tRNA 3' regions

```bash
awk 'BEGIN{FS="\t"}{if($3=="tRNA"){if($7~"+"){posEnd=$5}else{posEnd=$4};print "NC_000913.3:"posEnd"-"posEnd}}' NC_000913.gff > NC_000913_tRNA_3prime.list
for rep in "A" "B" "C" ; do for sample in "3-D" "4-NT" "5-nD" ; do 
   rm ${sample}_m${map}_tRNA_3prime_count.txt ; 
   for t in `more NC_000913_tRNA_3prime.list` ; do 
      samtools view 2_mapping/${rep}${sample}_CCATGG.bam ${t} | wc -l >> ${rep}${sample}_tRNA_3prime_count.txt ;
   done ; 
done ; done
```

## tRNA modified bases: conversion from Modomics DB to gff

### Data

- Protocol: Get sequences with modified bases from [modomics DB](http://modomics.genesilico.pl/sequences/list/tRNA/) for the *Escherichia coli* specie, acces: clic on "Display as ASCII" buton and copy/paste in text format file. Manually apply 2 modifications: i) deduplicate 4 tRNA names for Ini_CAU, Thr_GGU, Tyr_QUA, Val_GAC), and ii) duplicate the "_" character of selC following the footnote of the Modomics page. Create 2 fasta files from `bmModomics_nov17.txt`: i) without any bases but modified ones (`bmModomics_nov17_noBM.fasta`) and ii) without modified bases (`bmModomics_nov17_seqU.fasta`).
- Code:
```bash
sed 'n;s/[AGCU_]/ /g' bmModomics_nov17_info.fasta > bmModomics_nov17_noBM.fasta
sed 's/-//g;s/> tRNA/>tRNA/g;s/ | Escherichia coli | prokaryotic cytosol//g;s/ | /_/g;' bmModomics_nov17.txt > bmModomics_nov17_seqU.fasta
```
- Result files: 
  - 43 tRNA sequences with ?? knowed modified bases, `bmModomics_nov17.txt`
  - without any bases but modified ones (`bmModomics_nov17_noBM.fasta`)
  - without modified bases (`bmModomics_nov17_seqU.fasta`)

### Genomic coordinates of modified bases

- Protocol: alignment of modomics sequences to genomic sequence (blastn)
- Code:
```bash
blast
blast analysis
```
- Result file: `tRNA_feature.txt` 

## Termination signal: from read coverage to ts-jump

## Software version used

- fastqc: from docker-registry.genouest.org/ifb/fastqc, version 0.11.5
- cutadapt: from docker-registry.genouest.org/bioconda/cutadapt, version 1.11
- bowtie2: from docker-registry.genouest.org/bioconda/bowtie2, bowtie2-align-s version 2.2.8 
- samtools: version 1.4

