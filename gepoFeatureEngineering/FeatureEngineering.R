library("protr")


#read fasta file
essential = readFASTA("D_melanogaster_Essential_PROTEINS.txt")
essential = essential[(sapply(essential, protcheck))]
#extract AAC Amino Acid Composotion
a1 = t(sapply(essential, extractAAC))
colnames(a1) <- paste("AAC", colnames(a1), sep = "_")
#print(colnames(a1))
#extract DC Diptide composition 
a2 =t(sapply(essential, extractDC))
colnames(a2) <- paste("DC", colnames(a2), sep = "_")
#Tripetide 
a3 =t(sapply(essential, extractTC))
colnames(a3) <- paste("TC", colnames(a3), sep = "_")
#print(a3)
#CTriad
a4 =t(sapply(essential, extractCTriad))
colnames(a4) <- paste("CTriad", colnames(a4), sep = "_")
#print(dim(a4))

a5 = t(sapply(essential, extractCTDC))
colnames(a5) <- paste("CTDC", colnames(a5), sep = "_")
#print(a5)
#CTDT
a6 =t(sapply(essential, extractCTDT))
colnames(a6) <- paste("CTDT", colnames(a6), sep = "_")
#print(a6)
a7 =t(sapply(essential, extractCTDD))
colnames(a7) <- paste("CTDD", colnames(a7), sep = "_")
#print(a7)
a8 =t(sapply(essential, extractMoreauBroto))
colnames(a8) <- paste("Moreau", colnames(a8), sep = "_")
#print(a8)
a9 = t(sapply(essential, extractMoran)) # t(extractMoran(protSeq,nlag=30))
colnames(a9) <- paste("Moran", colnames(a9), sep = "_")
#print(a9)
a10 =t(sapply(essential, extractGeary))
colnames(a10) <- paste("Geary", colnames(a10), sep = "_")
#print(a10)
a11 = t(sapply(essential, extractSOCN))
colnames(a11) <- paste("SOCN", colnames(a11), sep = "_")
#print(a11)
a12 = t(sapply(essential, extractQSO)) #t( extractQSO(protSeq,nlag=30))
colnames(a12) <- paste("QSO", colnames(a12), sep = "_")
#print(a12)
a13 = t(sapply(essential, extractPAAC)) #t(extractPAAC(protSeq,lambda=30))
colnames(a13) <- paste("PAAC", colnames(a13), sep = "_")
#print(a13)
a14 = t(sapply(essential, extractAPAAC)) # t(extractAPAAC(protSeq,lambda=30))
colnames(a14) <- paste("APAAC", colnames(a14), sep = "_")
print(a14)
