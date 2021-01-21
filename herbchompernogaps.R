### HerbChomper, through-version Beta 0.1 - E. Gardner - January 21, 2021

### This version of HerbChomper trims an entire sequence end-to-end within a sliding window.

### Usage: Rscript herbchomper.R -a [alignment in] -o [alignment out] -t [sequence to trim] -r [reference sequence or "auto"] -w [size of sliding window] -i [identity cutoff] 

### The reference can be user specified or chosen automatically by specifying "-a auto". In the latter case, the sequence with the highest similarity to the target will be chosen, ignoring gaps and undetermined characters.

#get arguments
args = commandArgs(trailingOnly=TRUE)
#alignment in
seqfile<-args[grep("-a",args)+1]
#output alignment
outfile<-args[grep("-o",args)+1]
#sequence to trim
target<-args[grep("-t",args)+1]
#reference sequence
ref<-args[grep("-r",args)+1]
#sliding window size
slidingwindow<-as.numeric(args[grep("-w",args)+1])
#identity cutoff
identity<-as.numeric(args[grep("-i",args)+1])
#size of gap required to restart trimming

###### To use in the R environment, uncomment this section, set these seven variables and then run from here down
#working directory
setwd("~/scratch")
#alignment
seqfile<-("alignment.fasta")
#output alignment
outfile<-"alignment.chomped"
#sequence to trim
target<-"1"
#reference sequence
ref<-c("2")
#sliding window size
slidingwindow<-10
#identity cutoff
identity<-0.8

#check dependency and load alignments
if ("seqinr" %in% rownames(installed.packages()) == FALSE) {install.packages("seqinr")}
library(seqinr)
read.fasta(paste(seqfile))->gene
w<-slidingwindow-1

print(paste("Input file:",seqfile,collapse=" "))

#count the number of non-gap characters in the target sequence:
length1<-sum(gene[[target]][1:length(gene[[target]])]!="-")

#find where to start - forward
for (i in 1:length(gene[[target]])) {
	if (gene[[target]][i]=="-") {
		next
	}
	else {
		posF<-i
		break
	}
}

#find where to start - reverse
for (i in length(gene[[target]]):1) {
	if (gene[[target]][i]=="-") {
		next
	}
	else {
		posR<-i
		break
	}
}

#select the right reference if -a was set to "auto"

if (ref == "auto") {

#cut columns that are undetermined in the target sequence
posF->posKeep
for (i in (posF+1):posR) {
	if (gene[[target]][i] != "-" & gene[[target]][i] != "n") {
		posKeep <- c(posKeep,i)
	}
}
gene->geneCut
for (i in 1:length(geneCut)) {
geneCut[[i]]<-geneCut[[i]][posKeep]	
}
write.fasta(geneCut,names=names(gene),file=paste(seqfile,".cut.tmp",collapse="",sep=""))
read.alignment(paste(seqfile,".cut.tmp",collapse="",sep=""),"fasta")->alignment
system((paste("rm ", seqfile,".cut.tmp",collapse="",sep="")))

#calculate distances and pick the closest sequence as the reference. If multiple sequences are tied, one will be chosen arbitrarily.
as.matrix(dist.alignment(alignment, matrix = "identity",gap=1))->distances
sort(distances[target,colnames(distances)!=target])->targetDist
names(targetDist[targetDist==min(targetDist)])[1]->ref

#older slower verion of reference picking
#	similarity<-matrix(0,nrow=length(gene)-1,ncol=length(gene[[target]]))
#	rownames(similarity)<-names(gene[names(gene)!=target])
#	for (i in posF:posR) {
#		for (j in 1:nrow(similarity)) {
#			if (gene[[target]][i]==gene[[rownames(similarity)[j]]][i] & gene[[target]][i] != "-" & #gene[[target]][i] != "n" ) {
#				similarity[j,i]<-1
#			}
#			else {
#				next
#			}
#		}
#	}
#	rep(0,nrow(similarity))->simScores
#	rownames(similarity)->names(simScores)
#	for (i in 1:nrow(similarity)){
#		simScores[i]<-sum(similarity[i,])
#	}
#	names(sort(simScores))[length(simScores)]->ref

print(paste("Reference automatically set to ",ref,collapse=""))
}

print("CHOMPING... yum yum yum")

#scan forward
cutF<-0
for (i in posF:min(posR, length(gene[[target]])-slidingwindow)) {
	if (sum(gene[[target]][i:(i+w)]==gene[[ref]][i:(i+w)]) < (slidingwindow*identity)) {
		c(cutF,i)->cutF
		}
	else {
		next
		}
}

#scan reverse
cutR<--1
for (i in posR:max(posF, slidingwindow)) {
	if (sum(gene[[target]][i:(i-w)]==gene[[ref]][i:(i-w)]) < (slidingwindow*identity)) {
		c(cutR,i)->cutR
		}
	else {
		next
		}
}

#now cut the positions flagged in both the forward and reverse scans
gene[[1]][intersect(cutF,cutR)]<-"-"

print(paste(length(intersect(cutF,cutR))," characters out of ",length1," removed for target sequence ",target,sep="",collapse=""))
write.fasta(gene,names=names(gene),file=paste(outfile))
