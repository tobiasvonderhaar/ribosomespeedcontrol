import pandas as pd
import numpy as np
import pickle

#load data for the SSU footprints and sacCer3 gene coordinates
SSU = pd.read_csv('Archer SSU footprints.csv')
coordinates = pd.read_csv('sacCer3 genes.tsv', sep='\t')
#prepare variables for the processing loop
chromosomes = coordinates['chrom'].unique()[0:15]
results = {}
#process each chromosome in turn
for chromosome in chromosomes:
    this_SSU_set = SSU.loc[SSU['Chromosome']==chromosome]
    this_coordinate_set = coordinates.loc[coordinates['chrom']==chromosome]
    this_coordinate_set = this_coordinate_set.reset_index(drop=True)
    #process each gene on this chromosome
    for geneno in range(this_coordinate_set.shape[0]):
        #prepare an array for holding the footprint data for this gene
        this_gene = this_coordinate_set['name'][geneno]
        this_array = np.zeros((1,250))
        #data have to be processed differently depending on which strand the gene is on
        if this_coordinate_set['strand'][geneno] == '+':
            #calculate the start and stop coordinates of the +50/-200 nt window
            gene_start = this_coordinate_set.iloc[geneno]['txStart']-200
            gene_stop = gene_start + 250
            #now select only those SSUs within the start and stop coordinates
            larger_set = this_SSU_set.loc[this_SSU_set['aln_start'] >= gene_start]
            between_set = larger_set.loc[larger_set['aln_stop'] < gene_stop]
            #go through the individual footprints for this gene starting with the first
            for fp in range(between_set.shape[0]):
                fp_start = between_set.iloc[fp]['aln_start'] - gene_start
                fp_stop = between_set.iloc[fp]['aln_stop'] - gene_start
                this_array[0,fp_start:fp_stop] = between_set.iloc[fp]['counts']
        elif this_coordinate_set['strand'][geneno] == '-':
            #calculate the start and stop coordinates of the +50/-200 nt window
            gene_start = this_coordinate_set.iloc[geneno]['txEnd']+200
            gene_stop = gene_start - 250
            #now select only those SSUs within the start and stop coordinates
            larger_set = this_SSU_set.loc[this_SSU_set['aln_stop'] <= gene_start]
            between_set = larger_set.loc[larger_set['aln_start'] > gene_stop]
            #go through the individual footprints for this gene starting with the last
            for fp in range(between_set.shape[0]-1,0,-1):
                #re-calculate the footprint coordinates relative to the 250 nt window
                fp_start = gene_stop - between_set.iloc[fp]['aln_stop']
                fp_stop = gene_stop - between_set.iloc[fp]['aln_start']
                this_array[0,fp_start:fp_stop] = between_set.iloc[fp]['counts']

        results[this_gene] = this_array

with open(r"UTR_SSU_prints.p", "wb") as output_file:
    pickle.dump(results, output_file)
