#Get command arguments - species prefix
library(protr)
library(seqinr)

args<-commandArgs(TRUE)

geneFile <- args[1]

proteinFile <- args[2]

fileName <- args[3]

peptide_list_essential=readFASTA(proteinFile)
#load gene sequence
coding_list_essential=readFASTA(geneFile)
#check sequence
peptide_list_essential <- peptide_list_essential[(sapply(peptide_list_essential, protcheck))]
#write the new file check
write.fasta(sequences=peptide_list_essential,names=names(peptide_list_essential),file.out=paste(fileName,"sequenceAA.fasta",sep=""))

library("rDNAse")
#check sequences
coding_list_essential = coding_list_essential[(sapply(coding_list_essential, dnacheck))]
#write the new file check
write.fasta(sequences=coding_list_essential,names=names(coding_list_essential),file.out=paste(fileName,"sequenceNT.fasta",sep=""))
