def determine_polyx(seq,test_string):
    counter = 0
    if re.search(test_string,seq) != None:
        extended_string = test_string
        counter += 1
        while re.search(extended_string,seq) != None:
            extended_string += test_string
            counter += 1
        counter -= 1

    return counter

##############Extract Infromation from sequence data###########################

from Bio import SeqIO
import pandas as pd
import numpy as np
from scipy.stats.mstats import gmean
import pickle
import re

######process basic protein sequence features
print('Processing Protein Sequence Features')
Gene,Length,poly_P,poly_KR = [],[],[],[]
for record in SeqIO.parse('_Data/genome_orf_trans.fasta','fasta'):
    Gene.append(record.name)
    this_sequence = str(record.seq).replace('*','')
    #anal = ProteinAnalysis(this_sequence)
    #pI.append(anal.isoelectric_point())
    Length.append(len(this_sequence))
    poly_P.append(determine_polyx(this_sequence,'P'))
    poly_KR.append(determine_polyx(this_sequence, '[KR]'))

#assemble everything into the base dataframe
feature_frame = pd.DataFrame({'Gene': Gene,'Length':Length,'poly_P':poly_P,'poly_KR':poly_KR})

##############Add information in introns############################
print('Processing Intron Presence')
genomic_records = list(SeqIO.parse('_Data/genome_orf_genomic.fasta','fasta'))
coding_records = list(SeqIO.parse('_Data/genome_orf_coding.fasta','fasta'))

IntGene,has_intron = [],[]

for no in range(len(coding_records)):
    gen = str(genomic_records[no].seq)
    cod = str(coding_records[no].seq)
    IntGene.append(coding_records[no].id)
    if len(gen) > len(cod):
        has_intron.append(1)
    else:
        has_intron.append(0)

Intron=[]

for Gene in feature_frame['Gene']:
    if Gene in IntGene:
        Intron.append(has_intron[IntGene.index(Gene)])
    else:
        Intron.append(0)
feature_frame['Intron'] = Intron

################read in UTR and uORF information
print('Processing UTR length and uORF presence')
#import the literature data tables
Miura = pd.read_csv('_Data/Lit_Miura_2006_5_UTRs.csv')
Naga = pd.read_csv('_Data/Lit_Nagalakshmi_2008_UTRs.csv').fillna(0)
Waern5 = pd.read_csv('_Data/Lit_Waern_Snyder_2013_5_UTR.csv')
Waern3 = pd.read_csv('_Data/Lit_Waern_Snyder_2013_3_UTR.csv')
Ingolia = pd.read_csv('_Data/Lit_Ingolia_2009_uORFs.csv')

#look up the 5'-UTR annotations for each gene and record the highest value
UTR5,UTR3,uORF = [],[],[]

for Gene in feature_frame['Gene']:
    M5,N5,W5=0,0,0
    N3,W3=0,0
    if Gene in Miura['ORF'].values:
        M5 = int(max(Miura['5_UTR'].loc[Miura['ORF']==Gene].values))
    if Gene in Naga['ORF'].values:
        N5 = int(max(Naga['5_UTR'].loc[Naga['ORF']==Gene].values))
        N3 = int(max(Naga['3_UTR'].loc[Naga['ORF']==Gene].values))
    if Gene in Waern5['ORF'].values:
        W5 = int(max(Waern5['5_UTR'].loc[Waern5['ORF']==Gene].values))
    if Gene in Waern3['ORF'].values:
        W3 = int(max(Waern3['3_UTR'].loc[Waern3['ORF']==Gene].values))
    UTR5.append(max([M5,N5,W5]))
    UTR3.append(max([N3,W3]))

    #look up the uORF annotations and record whether or not uORFs are present
    if Gene in Ingolia['CDS'].values:
        uORF.append(1)
    else:
        uORF.append(0)

feature_frame['5_UTR'] = UTR5
feature_frame['3_UTR'] = UTR3
feature_frame['uORF'] = uORF

############Add structure scores#########################
print('Processing structure scores')
#extract pars scores from supplemental data file and convert to gene scores
G,pars=[],[]
with open('_Data/Lit_Kertesz_2010_parse_scores.tab') as f:
    for line in f:
        this_line = line.replace('\n','').split('\t')
        G.append(this_line[0])
        Len = float(this_line[1])
        Scores = np.array([float(i) for i in this_line[2].split(';')])
        pars.append(np.trapz(Scores)/Len)

pars_frame = pd.DataFrame({'Gene':G,'Structure':pars})

Structure = []
#Add pars scores to feature_frame as 'Structure'
for Gene in feature_frame['Gene']:
    if Gene in pars_frame['Gene'].values:
        Structure.append(pars_frame.loc[pars_frame['Gene'] == Gene]['Structure'].values[0])
    else:
        Structure.append(0)
feature_frame['Structure'] = Structure

####################Calculate tAI values#######################
print('Processing tAI Values')

#read in the omegavalues for each codon
omegas = pd.read_csv('_Data/Lit_Sabi_2016_tAI.csv', index_col='codon')

tAIGene,tAIcalc = [],[]

for no in range(len(coding_records)):
    tAIGene.append(str(coding_records[no].id))
    this_nt_seq = str(coding_records[no].seq)
    #split sequence into codons
    this_codon_seq = []
    for index in range(0, len(this_nt_seq), 3):
        this_codon_seq.append(this_nt_seq[index:index+3])
    tAIcalc.append(gmean(np.array([omegas['omega'][codon] for codon in this_codon_seq[:-1]])))

tAI=[]

for Gene in feature_frame['Gene']:
    if Gene in tAIGene:
        tAI.append(tAIcalc[tAIGene.index(Gene)])
    else:
        tAI.append(0)
feature_frame['tAI'] = tAI

#################calculate decoding speeds#######################
print('Processing Translation Rates')

#load the table with the decoding time for each codon from Chu et al 2014
codonTimes = pd.read_csv('_Data/Lit_Chu_2014_codon_decoding_times.csv', index_col='codon')
#############coding sequence list should stillbe in memory from earlier data
speedGene,AUGSpeed,totalSpeed,slowestSpeed = [],[],[],[]
#process each sequence
for no in range(len(coding_records)):
    speedGene.append(str(coding_records[no].id))
    this_nt_seq = str(coding_records[no].seq)
    #split sequence into codons
    this_codon_seq = []
    for index in range(0, len(this_nt_seq), 3):
        this_codon_seq.append(this_nt_seq[index:index+3])
    this_speed_seq = np.array([codonTimes['decoding.time'][codon] for codon in this_codon_seq[:-1]])
    #record the average decoding time of the first ten codons
    AUGSpeed.append(np.mean(this_speed_seq[:10]))
    #record the average decoding time of the entire sequence
    totalSpeed.append(np.mean(this_speed_seq))
    #record the average deocoding time of the slowest 20 codon window (if sequence > 60 nt/20 codons)
    if(len(this_speed_seq) > 60):
        slowestSpeed.append(max([np.mean(this_speed_seq[n:n+20]) for n in range(len(this_speed_seq)-20)]))
    else:
        slowestSpeed.append(0)

speed5,speedAv,speedmin = [],[],[]

for Gene in feature_frame['Gene']:
    if Gene in speedGene:
        speed5.append(AUGSpeed[speedGene.index(Gene)])
        speedAv.append(totalSpeed[speedGene.index(Gene)])
        speedmin.append(slowestSpeed[speedGene.index(Gene)])
    else:
        speed5.append(0)
        speedAv.append(0)
        speedmin.append(0)

feature_frame['5Speed'] = speed5
feature_frame['AvSpeed'] = speedAv
feature_frame['minSpeed'] = speedmin


#Import Translation initiation rates from von der Haar 2008
initRatesFrame = pd.read_csv('_Data/Lit_vonderHaar_2008.csv')
initRatesFrame.columns = ['Gene','GeneName','MW','L','ProtAbundance','Decay','RNAAbundance']
initGene = list(initRatesFrame['Gene'])
initRates = [row['ProtAbundance'] * row['Decay'] / row['RNAAbundance'] for index,row in initRatesFrame.iterrows()]
init = []

for Gene in feature_frame['Gene']:
    if Gene in initGene:
        init.append(initRates[initGene.index(Gene)])
    else:
        init.append(0)
feature_frame['InitH'] = init

#Import translation initiation rates from Ciandrini 2013

InitC = pd.read_csv('_Data/Lit_Ciandrini_2013_Init_Rates.csv')[['Gene name (systematic)','alpha_phys']]
InitC.columns = ['Gene','InitC']
feature_frame = feature_frame.merge(InitC, how = 'left', on = 'Gene')


#################Retrieve SSU Peak information
##Calculate SSU peak ratios
print('Processing SSU Peak Ratios')

with open(r"_Data/An_UTR_SSU_prints.p", "rb") as input_file:
    footprint_dict = pickle.load(input_file)

peakGene, peakRatio = [],[]
for key,value in footprint_dict.items():
    peakGene.append(key)
    #discard genes without any peaks
    if max(value[0]) < 2:
        peakRatio.append(0)
    else:
        #if there is no AUG peak, replace this with a small value to avoid divide-by-zero errors
        if max(value[0][180:230]) <= 0:
            value[0][200] = 2
        this_ratio = max(value[0][130:180]) / max(value[0][180:230])
        peakRatio.append(this_ratio)

#do an inner merge with the feature dataframe

ratios = pd.DataFrame({'Gene':peakGene,'PeakRatio':peakRatio})

feature_frame = feature_frame.merge(ratios,how='inner',on='Gene')

##retrieve SSU classifications
peakFrame = pd.read_csv('_Data/An_firstTenNew_SSU.csv')
peakCalls = peakFrame[['ORF','second_SSU_peak']]
peakCalls.columns=['Gene','PeakCall']
feature_frame = feature_frame.merge(peakCalls,how='inner',on='Gene')

######reorder the features by theme

feature_order=['Gene','PeakRatio','PeakCall','5_UTR','uORF','Structure','InitH','InitC','Length','Intron','5Speed','AvSpeed','minSpeed','tAI','3_UTR','poly_P','poly_KR']

feature_frame = feature_frame[feature_order]
feature_frame = feature_frame.set_index('Gene')

print('Dataframe Assembly Complete')

feature_frame.to_csv('_Data/An_mRNA_Feature_Compilation.csv')
