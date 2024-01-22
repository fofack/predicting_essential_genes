#Get command arguments - species prefix
library(protr)
library("rDNAse")

args<-commandArgs(TRUE)

geneFile <- args[1]

proteinFile <- args[2]

peptide_list_essential=readFASTA(proteinFile)
coding_list_essential=readFASTA(geneFile)



peptide_list_essential <- peptide_list_essential[(sapply(peptide_list_essential, protcheck))]
print(peptide_list_essential)
#peptide_list_nessential <- peptide_list_nessential[(sapply(peptide_list_nessential, protcheck))]
#peptide_list_cessential <- peptide_list_cessential[(sapply(peptide_list_cessential, protcheck))]
