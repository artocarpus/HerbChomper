### HerbChomper Beta 0.21 - for review - E. Gardner - May 26, 2020

### Usage: Rscript herbchomper.R -a [alignment in] -o [alignment out] -t [sequence to trim] -r [reference sequence or "auto"] -w [size of sliding window] -i [identity cutoff] -g [gap size necessary to restart trimming]

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
gapsize<-as.numeric(args[grep("-g",args)+1])

###### To use in the R environment, uncomment this section, set these seven variables and then run from here down
##working directory
#setwd("~/Downloads")
##alignment
#seqfile<-("test.fasta")
##sequence to trim
#target<-"A_jarrettiae_SAN120933"
##reference sequence
#ref<-c("A_elasticus_EG87")
##sliding window size
#slidingwindow<-10
##identity cutoff
#identity<-0.8
##size of gap required to restart trimming
#gapsize<-20

#check dependency and load alignments
if ("seqinr" %in% rownames(installed.packages()) == FALSE) {install.packages("seqinr")}
library(seqinr)
read.fasta(paste(seqfile))->gene
w<-slidingwindow-1
g<-gapsize-1

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
	similarity<-matrix(0,nrow=length(gene)-1,ncol=length(gene[[target]]))
	rownames(similarity)<-names(gene[names(gene)!=target])
	for (i in posF:posR) {
		for (j in 1:nrow(similarity)) {
			if (gene[[target]][i]==gene[[rownames(similarity)[j]]][i] & gene[[target]][i] != "-" & gene[[target]][i] != "N" ) {
				similarity[j,i]<-1
			}
			else {
				next
			}
		}
	}
	rep(0,nrow(similarity))->simScores
	rownames(similarity)->names(simScores)
	for (i in 1:nrow(similarity)){
		simScores[i]<-sum(similarity[i,])
	}
	names(sort(simScores))[length(simScores)]->ref
print(paste("Reference automatically set to ",ref,collapse=""))
}

#trim forward end
for (i in posF:posR) {
	if (sum(gene[[target]][i:(i+w)]==gene[[ref]][i:(i+w)]) < (slidingwindow*identity)) {
		gene[[target]][i]<-"-"
		}
	else {
		break
		}
}

#trim reverse
for (i in posR:posF) {
	if (sum(gene[[target]][i:(i-w)]==gene[[ref]][i:(i-w)]) < (slidingwindow*identity)) {
		gene[[target]][i]<-"-"
		}
	else {
		break
		}
}

#fix up starting and ending positions so that we don't get errors by running up against the ends
posF<-max(posF,gapsize)
posR<-min(posR,length(gene[[target]])-gapsize-1)

print("CHOMPING... yum yum yum")

#trim forward ends of gaps
for (i in posF:posR) {
	if (sum(gene[[target]][i:(i+g)]=="-")==gapsize) {
		i+g->gapF	## position to start trimming after gap
		gapF<-min(gapF,posR)
			for (j in gapF:posR) {
			if (sum(gene[[target]][j:(j+w)]==gene[[ref]][j:(j+w)]) < (slidingwindow*identity)) {
			gene[[target]][j]<-"-"
			}
			else {
				break
			}
		}
		}
	else {
		next
		}
}


#trim reverse ends of gaps
for (i in posR:posF) {
	if (sum(gene[[target]][i:(i-g)]=="-")==gapsize) {
		i-g->gapR	## position to start trimming after gap
		gapR<-max(gapR,posF)
			for (j in gapR:posF) {
			if (sum(gene[[target]][j:(j-w)]==gene[[ref]][j:(j-w)]) < (slidingwindow*identity)) {
			gene[[target]][j]<-"-"
			}
			else {
				break
			}
		}
		}
	else {
		next
		}
}

length2<-sum(gene[[target]][1:length(gene[[target]])]!="-")
print(paste(length1-length2," non-gap characters out of ",length1," removed for target sequence ",target,sep="",collapse=""))
write.fasta(gene,names=names(gene),file=paste(outfile))
