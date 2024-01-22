
from Bio import SeqIO
import os
from os import listdir
import re
import numpy as np
import csv
import pandas as pd
from sklearn import preprocessing
import tkinter as tk
from tkinter import filedialog


class Util:

    # function to retrieve the elements of a directory
    def find_dir(self, path_to_dir):
        dirnames = listdir(path_to_dir)
        return [dirname for dirname in dirnames if dirname.endswith]

    # function to retrieve fasta files from a directory
    def find_fasta_filenames(path_to_dir, suffix=".fasta"):
        filenames = listdir(path_to_dir)
        return [filename for filename in filenames if filename.endswith(suffix)]

    # cleaning of gene sequences to remove sequences <=30 characters and descriptions in front of identifiers
    def clean_sequence_dna(path_input, path_out):
        with open(path_input) as handle:
            with open(path_out, 'a') as clean_file:
                for record in SeqIO.parse(handle, "fasta"):
                    if (len(record.seq) > 30):
                        clean_file.write('>' + "%s" %
                                         record.id[4:].upper()+'\n')
                        clean_file.write("%s" % record.seq+'\n')

    # cleanup of protein sequences to remove sequences <=30 characters and descriptions in front of identifiers
    def clean_sequence_protein(path_input, path_out):
        with open(path_input) as handle:
            with open(path_out, 'a') as clean_file:
                for record in SeqIO.parse(handle, "fasta"):
                    if (len(record.seq) > 30):
                        clean_file.write('>' + "%s" %
                                         record.id[4:].upper()+'\n')
                        clean_file.write("%s" % record.seq+'\n')

    # cleaning of sequences with certain codons that do not correspond to any amino acid
    def delete_sequence_id(path_input, path_out, tab_id):
        with open(path_input) as handle:
            with open(path_out, 'a') as clean_file:
                for record in SeqIO.parse(handle, "fasta"):
                    operator = record.id in tab_id
                    if (operator == False):
                        clean_file.write('>' + "%s" % record.id+'\n')
                        clean_file.write("%s" % record.seq+'\n')

    # function to transform an array of array into a simple array
    def transform_in_simple_array(tab):
        probleme_codonw = []
        p = re.compile(
            "[A-Z]+[0-9]*[\_]?[A-Z]*[0-9]+[\_]?[0-9]*[\-]?[A-Z]*[0-9]*")
        for i in range(len(tab)):
            if (len(p.findall(tab[i]))):
                probleme_codonw.append(p.findall(tab[i]))
        id_probleme = np.array(probleme_codonw).flatten()
        return id_probleme

    def split_data(self, path_dir, path_out):
        file = self.find_dir(path_dir)
        # file.remove('.ipynb_checkpoints')
        print("split operation beging")
        for i in range(len(file)):
            with open(path_dir+file[i], "r") as handle:
                data = csv.reader(handle)
                rows = [r for r in data]
                rows[0].pop(0)
                with open(path_out+file[i], "w") as f:
                    writer = csv.writer(f)
                    for j in range(len(rows)):
                        writer.writerow(rows[j])

    def labeling(epath_file, data_csv, group_name, organism):
        df = pd.read_csv(data_csv)
        df['essentiality'] = -1
        id_kegg = df.index
        epath = pd.read_excel(epath_file, sheet_name='Sheet1')
        epath['Gene_Locus'] = epath['Gene_Locus'].str.upper()
        epath.set_index("Gene_Locus", inplace=True)
        gene_locus = epath.index
        for i in range(len(id_kegg)):
            if (id_kegg[i] in np.array(gene_locus)):
                df.at[id_kegg[i], 'essentiality'] = epath.at[id_kegg[i],
                                                             'Gene_essentaility']
        df.drop(df[df.essentiality == -1].index, inplace=True)
        mappings = {'NE': 0, 'E': 1}
        classes = df.pop('essentiality')
        essential_labels = classes.map(mappings)
        df['essentiality_status'] = essential_labels
        df = df.reset_index(drop=True)
        features = df.columns[:-1]
        target = ["essentiality_status"]
        print("Traitement du fichier:" + "  " + data_csv)
        df_pr = preprocessing.normalize(df[features])
        df_d = pd.DataFrame(df_pr, columns=features)
        target_f = pd.DataFrame(df["essentiality_status"], columns=target)
        dataset_norm = pd.concat([df_d, target_f], axis=1)
        dataset_norm.to_csv(
            'datasets/'+organism+'/'+'EssentialFeature_collection_'+group_name+'.csv', index=False)
    
    def choice_folder(text):
        root = tk.Tk()
        root.withdraw()  # Cacher la fenêtre principale

        # Demander à l'utilisateur de choisir un dossier
        choice_folder = filedialog.askdirectory(title=text)

        return choice_folder
    
    def list_folder_items(folder):
        if os.path.exists(folder):
            # Utiliser os.listdir pour obtenir la liste des éléments dans le dossier
            elements = os.listdir(folder)
            return elements
        else:
            print("Le dossier spécifié n'existe pas.")
            return []

