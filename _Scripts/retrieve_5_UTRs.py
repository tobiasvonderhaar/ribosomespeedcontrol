import pandas as pd
import requests
import time

#Add uORF information from Ingolia et al. to the basic gene dataset
uORF_data = pd.read_csv('_Data/Lit_Ingolia_2009_uORFs.csv')
dataset = pd.read_csv('_Data/An_firstTenNew_SSU.csv',index_col='ORF')
uORFs=[]

for gene in dataset.index:
    if gene in uORF_data.CDS.values:
        uORFs.append(True)
    else:
        uORFs.append(False)

dataset['uORF'] = uORFs

#isolate those entries having BOTH waiting SSUs and uORFS
dataset_2SSU_only = dataset.loc[dataset['second_SSU_peak'] == True]
dataset_both_only = dataset_2SSU_only.loc[dataset_2SSU_only['uORF'] == True]

#Look up the 5'-UTR sequences from SGD and determine the position of the first AUG
uORF_pos = []
for lineno in range(dataset_both_only.shape[0]):
    #construct a search term to retrieve the gene sequence including 5'-UTR from SGD
    search_term = 'https://www.yeastgenome.org/run_seqtools?format=fasta&type=genomic&genes='
    search_term += dataset_both_only.index[lineno]
    search_term += '&strains=S288C&up='
    search_term += str(dataset_both_only.loc[dataset_both_only.index[lineno]]['5_UTR'])
    #retrieve the query from SGD
    response = requests.get(search_term)
    data = response.text
    #remove the header line and newline charcters
    seq = ''
    seq = seq.join(data.split('\n')[1:])
    #find location of first start codon
    first_start_location = seq.find('ATG') - dataset_both_only.loc[dataset_both_only.index[lineno]]['5_UTR']
    uORF_pos.append(first_start_location)
    time.sleep(5)
uORF_positions = pd.DataFrame({'Gene':dataset_both_only.index,'uORF_pos':uORF_pos})

uORF_positions.to_csv('_Data/An_uORF_positions.csv')
