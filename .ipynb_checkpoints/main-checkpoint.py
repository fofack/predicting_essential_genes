from Features import Feature
import sys
from Utils import Util
import os


def main():
    organisms = Util.find_dir(Util(), 'organism_sequences')
    for o in range(len(organisms)):
        Feature.main(Feature(), 'organism_sequences/' +
                     organisms[o], 'sequences_clean', organisms[0])
    # data_epath = Util.find_dir(Util(), 'epathFile')
    # orgData = Util.find_dir(Util(), 'results')
    # if not os.path.exists('datasets'):
    #     os.mkdir('datasets')
    # for i in range(len(data_epath)):
    #     for j in range(len(orgData)):
    #         if (str(orgData[j]) == data_epath[i].split("_")[0]):
    #             if not os.path.exists('datasets/'+orgData[j]):
    #                 os.mkdir('datasets/'+orgData[j])
    #             data_csv = Util.find_dir(Util(), 'results/'+orgData[j])
    #             for s in range(len(data_csv)):
    #                 group_name = data_csv[s].split("_")[2].split(".")[0]
    #                 Util.labeling(
    #                     'epathFile/'+data_epath[i], 'results/'+orgData[j]+'/'+data_csv[s], group_name, orgData[j])


if __name__ == "__main__":
    main()
