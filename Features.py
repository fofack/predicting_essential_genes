import iFeatureOmegaCLI as iF
import os
from Utils import Util
from Cleaning import Clean


class Feature:

    def get_feature_dna(selft, path_input_sequence, path_out):
        dna_descriptors = ["NAC", "ANF", "ENAC", "MMI", "LPDF", "EIIP", "DPCP", "TPCP",
                           "Kmer type 1", "RCKmer type 1", "PseEIIP", "PseKNC",
                           "PCPseDNC", "PCPseTNC", "SCPseDNC", "SCPseTNC", "PSTNPss",
                           "PSTNPds", "CKSNAP type 1", "DAC", "DACC", "DCC",
                           "TAC", "TACC", "TCC",]
        dna = iF.iDNA(path_input_sequence)
        for i in range(len(dna_descriptors)):
            dna.get_descriptor(dna_descriptors[i])
            if (dna_descriptors[i] == "Kmer type 1"):
                dna.to_csv(path_out+"EssentialFeature_collection_" +
                           "Kmer1"+".csv", index=True, header=True)
            elif (dna_descriptors[i] == "RCKmer type 1"):
                dna.to_csv(path_out+"EssentialFeature_collection_" +
                           "RCKmer1"+".csv", index=True, header=True)
            elif (dna_descriptors[i] == "CKSNAP type 1"):
                dna.to_csv(path_out+"EssentialFeature_collection_" +
                           "CKSNAP"+".csv", index=True, header=True)
            else:
                dna.to_csv(path_out+"EssentialFeature_collection_" +
                           dna_descriptors[i]+".csv", index=True, header=True)

    def get_feature_peptide(selft, path_file_sequence, path_out):
        protein_descriptors = ["CTDC", "CTDD", "CTDT", "EGAAC", "GAAC", "AC", "CC",
                               "ACC", "CKSAAP type 1", "KSCTriad", "DistancePair",
                               "CKSAAGP type 1", "GDPC type 1", "GTPC type 1", "ASDC",
                               "DPC type 1", "DPC type 2", "ZScale", "AAIndex",
                               "TPC type 1", "TPC type 2"]
        protein = iF.iProtein(path_file_sequence)
        for i in range(len(protein_descriptors)):
            protein.get_descriptor(protein_descriptors[i])
            if (protein_descriptors[i] == "CKSAAP type 1"):
                protein.to_csv(path_out+"EssentialFeature_collection_" +
                               "CKSAAP"+".csv", index=True, header=True)
            elif (protein_descriptors[i] == "DPC type 1"):
                protein.to_csv(path_out+"EssentialFeature_collection_" +
                               "DPC1"+".csv", index=True, header=True)
            elif (protein_descriptors[i] == "DPC type 2"):
                protein.to_csv(path_out+"EssentialFeature_collection_" +
                               "DPC2"+".csv", index=True, header=True)
            elif (protein_descriptors[i] == "TPC type 1"):
                protein.to_csv(path_out+"EssentialFeature_collection_" +
                               "TPC1"+".csv", index=True, header=True)
            elif (protein_descriptors[i] == "TPC type 2"):
                protein.to_csv(path_out+"EssentialFeature_collection_" +
                               "TPC2"+".csv", index=True, header=True)
            elif (protein_descriptors[i] == "GDPC type 1"):
                protein.to_csv(path_out+"EssentialFeature_collection_" +
                               "GDPC1"+".csv", index=True, header=True)
            elif (protein_descriptors[i] == "GTPC type 1"):
                protein.to_csv(path_out+"EssentialFeature_collection_" +
                               "GTPC1"+".csv", index=True, header=True)
            else:
                protein.to_csv(path_out+"EssentialFeature_collection_" +
                               protein_descriptors[i]+".csv", index=True, header=True)

    def get_features_ifeature(self, input_dna, out_dna, input_peptide, out_peptide):
        self.get_feature_dna(input_dna, out_dna)
        self.get_feature_peptide(input_peptide, out_peptide)

    def main(self, path_input, path_out, org):
        print("################" + " " +
              "Feature de l'organisme" + " " + org + "################")
        if not os.path.exists(path_out):
            os.mkdir(path_out)
        if os.path.exists(path_out+'/.ipynb_checkpoints'):
            os.system('rm -r {}'.format(path_out+'/.ipynb_checkpoints'))
        if not os.path.exists(path_out+'/'+org):
            os.mkdir(path_out+'/'+org)
        Clean.main_clean_sequences(
            path_input+"/", path_out+"/"+org+"/")
        Clean.clean_sequences(path_out+"/"+org+"/"+"dna_sequence.fasta",
                              path_out+"/"+org+"/"+"protein_sequence.fasta", path_out+"/"+org+"/")
        if not os.path.exists('results'):
            os.mkdir('results')
        if os.path.exists('results/.ipynb_checkpoints'):
            os.system('rm -r {}'.format('results/.ipynb_checkpoints'))
        if not os.path.exists('results/'+org):
            os.mkdir('results/'+org)
        pb = os.system(
            'codonw %s/%s/EssentialsequenceNT-clean.fasta results/%s/EssentialcodonwFeature.out results/%s/EssentialcodonwFeature.blk -nomenu -nowarn -all_indices -human' % (str(path_out), str(org), str(org), str(org)))
        print(pb)

        id_probleme = 0
        prefix = "Essential"
        if (pb != 0 and len(pb) > 0):
            id_probleme = Util.transform_in_simple_array(pb)
        if (id_probleme != 0 and len(id_probleme) > 0):
            if not os.path.exists('clean_sequence_copy'):
                os.mkdir('clean_sequence_copy')
            if not os.path.exists('clean_sequence_copy/'+org):
                os.mkdir('clean_sequence_copy/'+org)
            os.system(
                'mv sequences_clean/'+org+'/dna_sequence.fasta clean_sequence_copy/'+org)
            os.system(
                'mv sequences_clean/'+org+'/protein_sequence.fasta clean_sequence_copy/'+org)
            os.system(
                'rm sequences_clean/'+org+'/EssentialsequenceNT-clean.fasta')
            os.system(
                'rm sequences_clean/'+org+'/EssentialsequenceAA-clean.fasta')
            os.system('rm results/'+org+'/EssentialcodonwFeature.out')
            os.system('rm results/'+org+'/EssentialcodonwFeature.blk')
            Util.delete_sequence_id("clean_sequence_copy/"+org+"/dna_sequence.fasta",
                                    path_out+"/"+org+"/dna_sequence.fasta", id_probleme)
            Util.delete_sequence_id("clean_sequence_copy/"+org+"/protein_sequence.fasta",
                                    path_out+"/"+org+"/protein_sequence.fasta", id_probleme)
            Clean.clean_sequences(path_out+"/"+org+"/dna_Sequence.fasta",
                                  path_out+"/"+org+"/protein_sequence.fasta", path_out+"/"+org+"/")
            pb1 = os.system(
                'codonw %s/%s/EssentialsequenceNT-clean.fasta results/%s/EssentialcodonwFeature.out results/%s/EssentialcodonwFeature.blk -nomenu -nowarn -all_indices -human' % (str(path_out), str(org), str(org), str(org)))
            print(pb1)
            print(
                "************* Debut extraction des features de l'organisme" + " " + org + "**********")
            # extract sequence feature
            os.system(
                'Rscript gepoFeatureEngineering/scriptFeatureBySeq.R %s/%s/EssentialsequenceNT-clean.fasta %s/%s/EssentialsequenceAA-clean.fasta results/EssentialcodonwFeature.out results/%s/%s' % (str(path_out), str(org), str(path_out), str(org), str(org), str(prefix)))
        else:
            print(
                "************* Debut extraction des features de l'organisme" + " " + org + "**********")
            # extract sequence feature
            os.system(
                'Rscript gepoFeatureEngineering/scriptFeatureBySeq.R %s/%s/EssentialsequenceNT-clean.fasta %s/%s/EssentialsequenceAA-clean.fasta results/EssentialcodonwFeature.out results/%s/%s' % (str(path_out), str(org), str(path_out), str(org), str(org), str(prefix)))

        # iFeature
        if not os.path.exists('results_iFeat_tmp'):
            os.mkdir('results_iFeat_tmp')
        if not os.path.exists('results_iFeat_tmp/'+org):
            os.mkdir('results_iFeat_tmp/'+org)
        self.get_features_ifeature('sequences_clean/'+org+'/EssentialsequenceNT-clean.fasta', 'results_iFeat_tmp/',
                                   'sequences_clean/'+org+'/EssentialsequenceAA-clean.fasta', 'results_iFeat_tmp/')
        Util.split_data(Util(), 'results_iFeat_tmp/' +
                        org+'/', 'results/'+org+'/')
        os.system('rm -r results_iFeat_tmp/')

        os.system('rm results/'+org+'/EssentialcodonwFeature.out')
        os.system('rm results/'+org+'/EssentialcodonwFeature.blk')
        os.system('rm results/'+org+'/EssentialFeature_collection_Sequence.txt')
        if os.path.exists('results/'+org+'/EssentialFeature_collection_Nucleotide.csv'):
            os.system('rm results/'+org +
                      '/EssentialFeature_collection_Nucleotide.csv')
        if os.path.exists('results/'+org+'/EssentialFeature_collection_Sequence.csv'):
            os.system('rm results/'+org +
                      '/EssentialFeature_collection_Sequence.csv')

        print(
            "************* Fin extraction des features de l'organisme" + "**********")
