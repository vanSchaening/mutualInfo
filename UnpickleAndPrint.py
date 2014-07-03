from pprint import pprint
import cPickle as pickle

f = "/home/axolotl/Data/conserved_site_targets_cs_0.0_hsa.pkl"
#TFtargetInteractions/s0.8/TCGA_BRCA_DNase_TF_target_mapping.pkl"
pkl_file = open(f, 'rb')

data = pickle.load(pkl_file)
print data.keys()

pkl_file.close()
