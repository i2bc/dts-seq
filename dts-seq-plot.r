# Script for Termination Signal (ts) calculation and plots
# Daniel Gautheret, 2018

rm(list=ls())
library (Biostrings)

# ---------- INPUT FILES ---------------

datadir="./"
covdir="5_coverage/"

# tRNA coordinates and names
trnacoordF<-paste0(datadir,"6_tRNA_modification/NC_000913_tRNA.tsv")
# coli fasta sequence imported from: 
#   ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
#   (unzipped)
genomeF<-paste0(datadir, "1_rawData/GCF_000005845.2_ASM584v2_genomic.fna")
# Modomics modified bases (Fasta-like but shows only mod bases)
modomicsF<-paste0(datadir, "6_tRNA_modification/bmModomics_with_tmRNA.fasta")
# list of tRNAs for printing and names of "lead" tRNAs
trnaF=paste0(datadir,"7_termination_signal/trnaprint.txt")

# Genome mapping coverage files (BigWig-like but w/ 2 strands)
#   A 2-column file with one line per genome position
#   column 1: coverage for + strand
#   column 2: coverage for minus strand
covF=c(
  "A4-NT_CCATGG_depth_fr.txt",
  "A3-D_CCATGG_depth_fr.txt",
  "A5-nD_CCATGG_depth_fr.txt",
  "B4-NT_CCATGG_depth_fr.txt", 
  "B3-D_CCATGG_depth_fr.txt",
  "B5-nD_CCATGG_depth_fr.txt",
  "C4-NT_CCATGG_depth_fr.txt",
  "C3-D_CCATGG_depth_fr.txt",
  "C5-nD_CCATGG_depth_fr.txt"
)

# Color plalette and other graphical parameters
# A (NT D nD) = red     B (NT D nD) =  green     C (NT D nD) = blue
mypalette=c("red","red","red","#AAFF00","#AAFF00","#AAFF00","blue","blue","blue")
maxx=90


# ---------- FUNCTIONS ---------------

# returns tRNA sequence from coordinates and genome (5' to 3') 
trnaseq5 <- function (genome, t) {
  if (t$strand=="+"){strand="plus"}
  if (t$strand=="-"){strand="minus"}
  v=genome[t$start:t$end]
  if (strand=="minus") {
    v=reverseComplement(v)
  }  
  return(as.character(v))
}

# compute jumps (ts) for 1 trna sequence
trnajump <- function (cov, from, to, strand) {
  v=cov[from:to,strand]
  if (strand=="plus") {
    v=rev(v)
  }  
  n=length(v)-1
  vj=c(0) # sets first jump at 0 
  for (i in 1:n){
    if (v[i]>99) {
      if (v[i]<v[i+1]) {
        jump=0
      }else{
        jump=(v[i]-v[i+1])/v[i]*100
      }
      vj=c(vj,jump)
    } else {
      vj=c(vj,NA)
    }
  }
  if (length(vj)==0){vj=c(0)} # never return empty vector
  return(vj)
}

# plot jumps for a list of tRNAs
# jmp1: array of jump values for drawing colored marks
# jmP2: vector of jump values for drawing bars
# tr: tRNA list
# genome: genome sequence
# allmods: strings of modified bases for all tRNAs
# pal: color palette for tick marks

plot_jumps <- function (jmp1, jmp2, tr, genome, allmods, pal, maxx) {  
  fname=paste0("jmp-", substr(tr[1,"name"],1,3),".pdf")
  pdf(fname, paper="a4",  width=8, height=11) 
  # fixed frame size
  par(mfrow = c(6,1))
  #  par(mfrow = c(nrow(tr),1)) 
  for (i in (1:nrow(tr))) {
    t=tr[i,]
    tname=as.character(t$name)
    title=tr[i,"title"]
    maxx=90
    maxy=100
    nbrep=length(allcov)
    if (t$strand=="+"){strand="plus"}
    if (t$strand=="-"){strand="minus"}
    # compute jumps @ each position & replicates
    jrep=c()
    for (j in (1:nbrep)) {
      cov=allcov[[j]]
      jmp=rev(trnajump(cov,t$start,t$end,strand))
      jmpl=length(jmp)
      # padding to maxx
      fill=rep(NA,max(0,maxx-jmpl))
      jmp=c(fill, jmp)
      print(length(jmp))
      jrep=cbind(jrep, jmp)
    } 
    # mean jump
    meanj=apply(jrep,1,mean,na.rm=T)
    seq=unlist(strsplit(trnaseq5(genome,t),split=''))
    mods=unlist(allmods[tname])
    # padd the sequence too
    fill=rep(" ",max(0,maxx-(length(seq))))
    seq=c(fill,seq)
    mods=c(fill,mods)
    par(mar=c(2, 2, 2, 0))
    bp=barplot (meanj, border=NA, main=title, col="grey",ylim=c(0, maxy), xlim=c(0,110),axes=F)
    # trick to allow all sequence labels to be displayed  
    impair=seq(1,length(seq),by=2)
    pair=seq(0,length(seq),by=2)
    axis(1, at=bp[impair], label=seq[impair], cex.axis=0.7, line = -1, tick=F)
    axis(1, at=bp[pair], label=seq[pair], cex.axis=0.7, line = -1, tick=F)
    axis(1, at=bp[impair], label=mods[impair], cex.axis=0.7, line = -0.4, tick=F, font=2)
    axis(1, at=bp[pair], label=mods[pair], cex.axis=0.7, line = -0.4, tick=F, font=2)
    
    yt=seq(0,maxy,by=10)
    axis(2, cex.axis=0.7, at=yt, label=yt,line=-1)
    for (j in (1:nbrep)) {
      jmp=jrep[,j]
      #points(x=bp, y=jmp, col=rainbow(nbrep)[j], pch="-", font=2)
      points(x=bp, y=jmp, col=mypalette[j], pch="-", font=2)
    } 
  }
  dev.off()
}

# Compute trna jumps (ts) for tRNA list with replicates, 5' to 3'
# tr: a trna list
# allcov: array with all genome-wide coverages
# returns array with jumps for each replicate
# jumps are right justified as here (T=tRNA, J=jumps):
#    <---------maxx-------->
# 5' --TTTTTTTTTTTTTTTTTTTTT 3'
#    ---------JJJJJJJJJJJJJJ 

trnajump_list <- function (allcov, tr, maxx) {  
  for (i in (1:nrow(tr))) {
    t=tr[i,]
    print (t$name)
    nbrep=length(allcov)
    if (t$strand=="+"){strand="plus"}
    if (t$strand=="-"){strand="minus"}
    # compute jumps @ each position & replicates
    jrep=c()
    for (j in (1:nbrep)) {
      cov=allcov[[j]]
      jmp=rev(trnajump(cov,t$start,t$end,strand))
      jmpl=length(jmp)
      # padding to maxx
      fill=rep(NA,max(0,maxx-jmpl))
      jmp=c(fill, jmp)
      jrep=cbind(jrep, jmp)
    } 
    return(jrep)
  }
}

# mean jump

trnajump_mean <- function (allcov, tr, maxx) {  
  jrep<-trnajump_list (allcov, tr, maxx) 
  meanj=apply(jrep,1,mean,na.rm=T)
  return(meanj)
}


# ------------------------
#          MAIN
# ------------------------

#---------- read tRNA names and coordinates

alltrnas=read.table(trnacoordF, header=F, sep="\t")
trnanames=as.character(alltrnas$name)
colnames(alltrnas)=c("start","end","strand","name")
alltrnas=alltrnas[order(alltrnas$name),]

#---------- read genome

allChr<-readDNAStringSet(genomeF)
genome<-allChr[[1]]
genome=RNAString(genome)

#-------- read modomics strings
# beware: result is an array of lists

trnalines = readLines(modomicsF)
allmods<-c()
tnames<-c()
for (l in trnalines) {
  print (nchar(l))
  if (grepl (">", l)) {
    tnames=c(tnames,substr(l,2,5))
  }else{
    #      v<-unlist(strsplit(l,""))
    v<-strsplit(l,"")
    allmods <- cbind(allmods,v)
  }
}
allmods=array(allmods)
rownames(allmods)=tnames

#-------- read coverage files

allcov <-list()
for (f in covF) {
  f2=paste0(covdir,f)
  cov=read.table(f2, header=F, sep="\t")
  colnames(cov)=c("plus","minus")
  allcov=append(allcov,list(cov))
  print (nrow(cov)) 
}
names(allcov)=sub("(.*)\\..*", "\\1", covF)   # keeps file prefix as lib name
# then for accessing one cov vector: allcov[[i]]  (1 coverage track across whole genome)
# example: allcov[[1]][1:10,"minus"]

# read  list of tRNAs to be printed
trnalines = readLines(trnaF)

# ------- Plot ts based on the sum of all coverages


printlist=c()
for (l in trnalines) {
  v<-unlist(strsplit(l," "))
  if (v[1]!="--") {
    tr=alltrnas[grep(v[1], alltrnas$name),]
    tr$title=l
    printlist=rbind(printlist,tr)
  } else {
    print (printlist)
    printlist=c()
  }
}

# cimput sum of all coverages
covp=c()
covm=c()
nbrep=length(allcov)
for (i in (1:nbrep)) {
  covp=cbind(covp,allcov[[i]][,"plus"])
  covm=cbind(covm,allcov[[i]][,"minus"])
}
covp1=rowSums(covp)
covm1=rowSums(covm)
sumcov <- data.frame(covp1,covm1)
colnames(sumcov)=c("plus","minus")

#----------
trnajump_list(list(sumcov), alltrnas[1,], maxx)


