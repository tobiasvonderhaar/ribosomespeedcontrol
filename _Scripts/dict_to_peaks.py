from scipy.signal import find_peaks_cwt
import numpy as np
import pickle

#load the dictionary holding the footprint traces
with open(r"UTR_SSU_prints.p", "rb") as input_file:
    footprint_dict = pickle.load(input_file)

#prepare a new dictionary for holding peaks
peak_dict = {}

#go through each item of footprints_dict and determine peaks
for gene, this_footprint in footprint_dict.items():
    data = this_footprint.flatten()
    this_test = find_peaks_cwt(data, widths=np.arange(1,50)) - 200
    peak_dict[gene] = this_test

#save the peaks dictionary using pickle
with open(r"UTR_SSU_peaks.p", "wb") as output_file:
    pickle.dump(peak_dict, output_file)
