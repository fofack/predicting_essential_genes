
#!bin/bash
#echo -e "Enter organism abreviation:\n"
#read name
mkdir organism_sequences
mkdir folder_tmp
mkdir folder_tmp/id
wget -nc https://rest.kegg.jp/list/$1 --output-document=folder_tmp/$1.txt

#extraction of gens numbers
pwd
cut -f 1 folder_tmp/$1.txt >> folder_tmp/id/$1.txt
mkdir organism_sequences/$1

while read line 
do gene=$(wget -qO- https://rest.kegg.jp/get/$line/NTSEQ)
   protein=$(wget -qO- https://rest.kegg.jp/get/$line/AASEQ)
   echo "$gene" >> organism_sequences/$1/dna_sequence.fasta
   echo "$protein" >> organism_sequences/$1/protein_sequence.fasta	
done < folder_tmp/id/$1.txt
