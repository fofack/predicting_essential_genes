from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
#load dataset
#aa
fastaSequences = SeqIO.parse(open("protein_sequence.fasta"), 'fasta')
seqaaDict_E=dict()
for fasta in fastaSequences:
  name, sequence = fasta.id, str(fasta.seq)
  seqaaDict_E[name]=sequence
#print(seqaaDict_E)

#aa
fastaSequences = SeqIO.parse(open("DNA_sequence.fasta"), 'fasta')
seqntDict_E=dict()
for fasta in fastaSequences:
  name, sequence = fasta.id, str(fasta.seq)
  seqntDict_E[name]=sequence
#print(seqntDict_E)
index=0
fileSeqAA=open("EssentialsequenceAA-clean.fasta", 'a')
fileSeqNT=open("EssentialsequenceNT-clean.fasta", 'a')
for cle in seqntDict_E.keys():
	if cle in seqaaDict_E:
		print(cle)
		sequenceProt=seqaaDict_E[cle].upper()
		sequenceGene=seqntDict_E[cle].upper()
		gene_Locus=cle
		recordAA=SeqRecord(Seq(sequenceProt), id=str(gene_Locus),description="")
		recordNT=SeqRecord(Seq(sequenceGene), id=str(gene_Locus),description="")
		SeqIO.write(recordAA, fileSeqAA, "fasta")
		SeqIO.write(recordNT, fileSeqNT, "fasta")
		index=index+1
