#Get command arguments - species prefix
library(protr)

args<-commandArgs(TRUE)

geneFile <- args[1]

proteinFile <- args[2]

featureFile <- args[3]

peptide_list_essential=readFASTA(proteinFile)
coding_list_essential=readFASTA(geneFile)



peptide_list_essential <- peptide_list_essential[(sapply(peptide_list_essential, protcheck))]
#print(peptide_list_essential)
print("Begin Extract feature to AA sequence")
#peptide_list_nessential <- peptide_list_nessential[(sapply(peptide_list_nessential, protcheck))]
#peptide_list_cessential <- peptide_list_cessential[(sapply(peptide_list_cessential, protcheck))]


a1 = t(sapply(peptide_list_essential, extractAAC))
#a2 = t(sapply(peptide_list_nessential, extractAAC))

#c1 = t(sapply(peptide_list_cessential, extractAAC))

colnames(a1) <- paste("AAC", colnames(a1), sep = "_")
#colnames(a2) <- paste("AAC", colnames(a2), sep = "_")

#colnames(c1) <- paste("AAC", colnames(c1), sep = "_")

a3 = t(sapply(peptide_list_essential, extractDC))
#a4 = t(sapply(peptide_list_nessential, extractDC))

#c2 = t(sapply(peptide_list_cessential, extractDC))

colnames(a3) <- paste("DC", colnames(a3), sep = "_")
#colnames(a4) <- paste("DC", colnames(a4), sep = "_")

#colnames(c2) <- paste("DC", colnames(c2), sep = "_")

a5 = t(sapply(peptide_list_essential, extractTC))
#a6 = t(sapply(peptide_list_nessential, extractTC))

#c3 = t(sapply(peptide_list_cessential, extractTC))

colnames(a5) <- paste("TC", colnames(a5), sep = "_")
#colnames(a6) <- paste("TC", colnames(a6), sep = "_")

#colnames(c3) <- paste("TC", colnames(c3), sep = "_")

a7 = t(sapply(peptide_list_essential, extractCTriad))
#a8 = t(sapply(peptide_list_nessential, extractCTriad))

#c4 = t(sapply(peptide_list_cessential, extractCTriad))

colnames(a7) <- paste("CTriad", colnames(a7), sep = "_")
#colnames(a8) <- paste("CTriad", colnames(a8), sep = "_")

#colnames(c4) <- paste("CTriad", colnames(c4), sep = "_")

a9 = t(sapply(peptide_list_essential, extractCTDC))
#a10 = t(sapply(peptide_list_nessential, extractCTDC))

#c5 = t(sapply(peptide_list_cessential, extractCTDC))

colnames(a9) <- paste("CTDC", colnames(a9), sep = "_")
#colnames(a10) <- paste("CTDC", colnames(a10), sep = "_")

#colnames(c5) <- paste("CTDC", colnames(c5), sep = "_")

a11 = t(sapply(peptide_list_essential, extractCTDT))
#a12 = t(sapply(peptide_list_nessential, extractCTDT))

#c6 = t(sapply(peptide_list_cessential, extractCTDT))

colnames(a11) <- paste("CTDT", colnames(a11), sep = "_")
#colnames(a12) <- paste("CTDT", colnames(a12), sep = "_")

#colnames(c6) <- paste("CTDT", colnames(c6), sep = "_")

a13 = t(sapply(peptide_list_essential, extractCTDD))
#a14 = t(sapply(peptide_list_nessential, extractCTDD))

#c7 = t(sapply(peptide_list_cessential, extractCTDD))

colnames(a13) <- paste("CTDD", colnames(a13), sep = "_")
#colnames(a14) <- paste("CTDD", colnames(a14), sep = "_")

#colnames(c7) <- paste("CTDD", colnames(c7), sep = "_")

a15 = t(sapply(peptide_list_essential, extractMoreauBroto, nlag=30))
#a16 = t(sapply(peptide_list_nessential, extractMoreauBroto, nlag=10))

#c8 = t(sapply(peptide_list_cessential, extractMoreauBroto, nlag=10))

colnames(a15) <- paste("Moreau", colnames(a15), sep = "_")
#colnames(a16) <- paste("Moreau", colnames(a16), sep = "_")

#colnames(c8) <- paste("Moreau", colnames(c8), sep = "_")

a17 = t(sapply(peptide_list_essential, extractMoran, nlag=30))
#a18 = t(sapply(peptide_list_nessential, extractMoran, nlag=10))

#c9 = t(sapply(peptide_list_cessential, extractMoran, nlag=10))

colnames(a17) <- paste("Moran", colnames(a17), sep = "_")
#colnames(a18) <- paste("Moran", colnames(a18), sep = "_")

#colnames(c9) <- paste("Moran", colnames(c9), sep = "_")

a19 = t(sapply(peptide_list_essential, extractGeary, nlag=30))
#a20 = t(sapply(peptide_list_nessential, extractGeary, nlag=10))

#c10 = t(sapply(peptide_list_cessential, extractGeary, nlag=10))

colnames(a19) <- paste("Geary", colnames(a19), sep = "_")
#colnames(a20) <- paste("Geary", colnames(a20), sep = "_")

#colnames(c10) <- paste("Geary", colnames(c10), sep = "_")

a21 = t(sapply(peptide_list_essential, extractSOCN, nlag=30))
#a22 = t(sapply(peptide_list_nessential, extractSOCN, nlag=10))

#c11 = t(sapply(peptide_list_cessential, extractSOCN, nlag=10))

colnames(a21) <- paste("SOCN", colnames(a21), sep = "_")
#colnames(a22) <- paste("SOCN", colnames(a22), sep = "_")

#colnames(c11) <- paste("SOCN", colnames(c11), sep = "_")

a23 = t(sapply(peptide_list_essential, extractQSO,nlag=30))
#a24 = t(sapply(peptide_list_nessential, extractQSO,nlag=10))

#c12 = t(sapply(peptide_list_cessential, extractQSO,nlag=10))

colnames(a23) <- paste("QSO", colnames(a23), sep = "_")
#colnames(a24) <- paste("QSO", colnames(a24), sep = "_")

#colnames(c12) <- paste("QSO", colnames(c12), sep = "_")

a25 = t(sapply(peptide_list_essential, extractPAAC,lambda=30))
#a26 = t(sapply(peptide_list_nessential, extractPAAC,lambda=10))

#c13 = t(sapply(peptide_list_cessential, extractPAAC,lambda=10))

colnames(a25) <- paste("PAAC", colnames(a25), sep = "_")
#colnames(a26) <- paste("PAAC", colnames(a26), sep = "_")

#colnames(c13) <- paste("PAAC", colnames(c13), sep = "_")

a27 = t(sapply(peptide_list_essential, extractAPAAC,lambda=30))
#a28 = t(sapply(peptide_list_nessential, extractAPAAC,lambda=10))

#c14 = t(sapply(peptide_list_cessential, extractAPAAC,lambda=10))

colnames(a27) <- paste("APAAC", colnames(a27), sep = "_")
#colnames(a28) <- paste("APAAC", colnames(a28), sep = "_")

#colnames(c14) <- paste("APAAC", colnames(c14), sep = "_")


print("Fin extraction Protr")
##rDNAse sequence features
library("rDNAse")

# check if the sequence is 
coding_list_essential = coding_list_essential[(sapply(coding_list_essential, dnacheck))]
#coding_list_essential <- geneSequence
#coding_list_nessential$FBgn0032224 <- NULL
#coding_list_cessential$FBgn0030501 <- NULL

a29 = t(sapply(coding_list_essential, kmer, 2))
#a30 = t(sapply(coding_list_nessential, kmer, 2))

#c15 = t(sapply(coding_list_cessential, kmer, 2))

colnames(a29) <- paste("kmer_2", colnames(a29), sep = "_")
#colnames(a30) <- paste("kmer_2", colnames(a30), sep = "_")

#colnames(c15) <- paste("kmer_2", colnames(c15), sep = "_")

a31 = t(sapply(coding_list_essential, kmer, 3))
#a32 = t(sapply(coding_list_nessential, kmer, 3))

#c16 = t(sapply(coding_list_cessential, kmer, 3))

colnames(a31) <- paste("kmer_3", colnames(a31), sep = "_")
#colnames(a32) <- paste("kmer_3", colnames(a32), sep = "_")

#colnames(c16) <- paste("kmer_3", colnames(c16), sep = "_")

a33 = t(sapply(coding_list_essential, extrDCC))
#a34 = t(sapply(coding_list_nessential, extrDCC))

#c17 = t(sapply(coding_list_cessential, extrDCC))

colnames(a33) <- paste("DCC", colnames(a33), sep = "_")
#colnames(a34) <- paste("DCC", colnames(a34), sep = "_")

#colnames(c17) <- paste("DCC", colnames(c17), sep = "_")

a35 = t(sapply(coding_list_essential, extrDACC))
#a36 = t(sapply(coding_list_nessential, extrDACC))

#c18 = t(sapply(coding_list_cessential, extrDACC))

colnames(a35) <- paste("DACC", colnames(a35), sep = "_")
#colnames(a36) <- paste("DACC", colnames(a36), sep = "_")

#colnames(c18) <- paste("DACC", colnames(c18), sep = "_")

a37 = t(sapply(coding_list_essential, extrTCC))
#a38 = t(sapply(coding_list_nessential, extrTCC))

#c19 = t(sapply(coding_list_cessential, extrTCC))

colnames(a37) <- paste("TCC", colnames(a37), sep = "_")
#colnames(a38) <- paste("TCC", colnames(a38), sep = "_")

#colnames(c19) <- paste("TCC", colnames(c19), sep = "_")

a39 = t(sapply(coding_list_essential, extrTACC))
#a40 = t(sapply(coding_list_nessential, extrTACC))

#c20 = t(sapply(coding_list_cessential, extrTACC))

colnames(a39) <- paste("TACC", colnames(a39), sep = "_")
#colnames(a40) <- paste("TACC", colnames(a40), sep = "_")

#colnames(c20) <- paste("TACC", colnames(c20), sep = "_")

a41 = t(sapply(coding_list_essential, extrPseDNC))
#a42 = t(sapply(coding_list_nessential, extrPseDNC))

#c21 = t(sapply(coding_list_cessential, extrPseDNC))

colnames(a41) <- paste("PseDNC", colnames(a41), sep = "_")
#colnames(a42) <- paste("PseDNC", colnames(a42), sep = "_")

#colnames(c21) <- paste("PseDNC", colnames(c21), sep = "_")

a43 = t(sapply(coding_list_essential, extrPseKNC, 3))
#a44 = t(sapply(coding_list_nessential, extrPseKNC, 3))

#c22 = t(sapply(coding_list_cessential, extrPseKNC, 3))

colnames(a43) <- paste("PseKNC_3", colnames(a43), sep = "_")
#colnames(a44) <- paste("PseKNC_3", colnames(a44), sep = "_")

#colnames(c22) <- paste("PseKNC_3", colnames(c22), sep = "_")

a45 = t(sapply(coding_list_essential, extrPseKNC, 5))
#a46 = t(sapply(coding_list_nessential, extrPseKNC, 5))

#c23 = t(sapply(coding_list_cessential, extrPseKNC, 5))

colnames(a45) <- paste("PseKNC_5", colnames(a45), sep = "_")
#colnames(a46) <- paste("PseKNC_5", colnames(a46), sep = "_")

#colnames(c23) <- paste("PseKNC_5", colnames(c23), sep = "_")

print("Fin extraction Nucleotide")

###############################Loading Feature vectors##################################

x1=cbind(a1,a3)
#x2=cbind(a2,a4)

#x3=cbind(c1,c2)

x1=cbind(x1,a5)
#x2=cbind(x2,a6)

#x3=cbind(x3,c3)

x1=cbind(x1,a7)
#x2=cbind(x2,a8)

#x3=cbind(x3,c4)

x1=cbind(x1,a9)
#x2=cbind(x2,a10)

#x3=cbind(x3,c5)

x1=cbind(x1,a11)
#x2=cbind(x2,a12)

#x3=cbind(x3,c6)

x1=cbind(x1,a13)
#x2=cbind(x2,a14)

#x3=cbind(x3,c7)

x1=cbind(x1,a15)
#x2=cbind(x2,a16)

#x3=cbind(x3,c8)

x1=cbind(x1,a17)
#x2=cbind(x2,a18)

#x3=cbind(x3,c9)

x1=cbind(x1,a19)
#x2=cbind(x2,a20)

#x3=cbind(x3,c10)

x1=cbind(x1,a21)
#x2=cbind(x2,a22)

#x3=cbind(x3,c11)

x1=cbind(x1,a23)
#x2=cbind(x2,a24)

#x3=cbind(x3,c12)

x1=cbind(x1,a25)
#x2=cbind(x2,a26)

#x3=cbind(x3,c13)

x1=cbind(x1,a27)
#x2=cbind(x2,a28)

#x3=cbind(x3,c14)

##

x1=cbind(x1,a29)
#x2=cbind(x2,a30)

#x3=cbind(x3,c15)

x1=cbind(x1,a31)
#x2=cbind(x2,a32)

#x3=cbind(x3,c16)

x1=cbind(x1,a33)
#x2=cbind(x2,a34)

#x3=cbind(x3,c17)

x1=cbind(x1,a35)
#x2=cbind(x2,a36)

#x3=cbind(x3,c18)

x1=cbind(x1,a37)
#x2=cbind(x2,a38)

#x3=cbind(x3,c19)

x1=cbind(x1,a39)
#x2=cbind(x2,a40)

#x3=cbind(x3,c20)

x1=cbind(x1,a41)
#x2=cbind(x2,a42)

#x3=cbind(x3,c21)

x1=cbind(x1,a43)
#x2=cbind(x2,a44)

#x3=cbind(x3,c22)

x1=cbind(x1,a45)
#x2=cbind(x2,a46)

#x3=cbind(x3,c23)

print("Fin extraction all and fusion")
print("generation des feature avec condonw")
codonwFile=paste(featureFile,"codonwFeature.out",sep="")
Essential_Codon=read.table(codonwFile,sep="\t",header=T)
#print(Essential_Codon)
x1=cbind(x1,as.matrix(as.data.frame(lapply(Essential_Codon[,3:dim(Essential_Codon)[2]],as.numeric))))
#print(x1)
print(dim(x1))

featureFileTXT=paste(featureFile,"Feature_collection_Sequence.txt",sep="")
featureFileCSV=paste(featureFile,"Feature_collection_Sequence.csv",sep="")
write.table(x1,featureFileTXT,quote=FALSE,sep="\t")
write.table(x1,featureFileCSV,quote=FALSE,sep=",")
#write.table(x2,"Feature_collection_SequenceNEssential.txt",quote=FALSE,sep="\t")
#write.table(x3,"Feature_collection_SequenceCEssential.txt",quote=FALSE,sep="\t")

#print(dim(x2))
#print(dim(x3))
