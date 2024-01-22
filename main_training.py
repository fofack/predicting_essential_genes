from Utils import Util
from Training import Training
import pandas as pd
import os


def main():
    subGroup = ["AAC", "AAIndex", "AC", "ACC", "ANF", "APAAC", "ASDC",
                "CC", "CKSAAGP", "CKSAAP", "CKSNAP", "CTDC", "CTDD",
                "CTDT", "Ctriad", "DAC", "DACC", "DCC", "DDE", "DPC1",
                "DPC2", "DPCP", "Dipetite", "DistancePair", "EGAAC",
                "EIIP", "ENAC", "GAAC", "GDPC1", "GTPC1", "Geary", "KSCTriad",
                "Kmer1", "LPDF", "MMI", "Moran", "Moreau", "NAC", "Nucleotide",
                "PAAC", "PCPseDNC", "PCPseTNC", "PSTNPds", "PSTNPss", "Protein",
                "PseDNC", "PseEIIP", "PseKNC", "PseKNC_3", "PseKNC_5", "QSO",
                "RCKmer1", "SCPseDNC", "SCPseTNC", "SOCN", "Sequence", "TAC",
                "TACC", "TCC", "TPC1", "TPC2", "TPCP", "Tripetite", "ZScale",
                "kmer_2", "kmer_3"]

    if not os.path.exists('data'):
        os.mkdir('data')
    if not os.path.exists('Models'):
        os.mkdir('Models')
    for i in range(len(subGroup)):
        Training.dataset_construction('datasets/', subGroupName=subGroup[i])

    # datasets = Util.find_dir(Util(), 'data/')

    # Training
    # print('debut entraînement')
    # for j in range(len(datasets)):
    #     dataset = pd.read_csv('data/'+datasets[i])
    #     subGroupName = datasets[i].split(".")[0]
    #     Training.Training(dataset=dataset, group_name=subGroupName)


# if __name__ == "__main__":
#     main()
