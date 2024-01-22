from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Utils import Util


class Clean:

    # fusion of the two functions of cleaning gene and protein sequences
    def main_clean_sequences(path_dir_input, path_dir_out):
        path_data = Util.find_fasta_filenames(path_dir_input)
        for i in range(len(path_data)):
            if (path_data[i] == "dna_sequence.fasta"):
                print(path_dir_out+path_data[i])
                Util.clean_sequence_dna(
                    path_dir_input+path_data[i], path_dir_out+"dna_sequence.fasta")
            else:
                Util.clean_sequence_protein(
                    path_dir_input+path_data[i], path_dir_out+"protein_sequence.fasta")

    # function to retain only sequences that are in both genes and proteins
    def clean_sequences(gene_sequence, protein_sequence, path_out):
        # load dataset
        # aa
        fastaSequences = SeqIO.parse(open(protein_sequence), 'fasta')
        seqaaDict_E = dict()
        for fasta in fastaSequences:
            name, sequence = fasta.id, str(fasta.seq)
            seqaaDict_E[name] = sequence
        # print(seqaaDict_E)
        # aa
        fastaSequences = SeqIO.parse(open(gene_sequence), 'fasta')
        seqntDict_E = dict()
        for fasta in fastaSequences:
            name, sequence = fasta.id, str(fasta.seq)
            seqntDict_E[name] = sequence
            # print(seqntDict_E)
        index = 0
        fileSeqAA = open(path_out+"EssentialsequenceAA-clean.fasta", 'a')
        fileSeqNT = open(path_out+"EssentialsequenceNT-clean.fasta", 'a')
        for cle in seqntDict_E.keys():
            if cle in seqaaDict_E:
                # print(cle)
                sequenceProt = seqaaDict_E[cle].upper()
                sequenceGene = seqntDict_E[cle].upper()
                gene_Locus = cle
                recordAA = SeqRecord(Seq(sequenceProt), id=str(
                    gene_Locus), description="")
                recordNT = SeqRecord(Seq(sequenceGene), id=str(
                    gene_Locus), description="")
                SeqIO.write(recordAA, fileSeqAA, "fasta")
                SeqIO.write(recordNT, fileSeqNT, "fasta")
                index = index+1
