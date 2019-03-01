

def process_blot(blot_data, control_strain = 'wt', control_Ab = 'Pgk1'):

    import pandas as pd
    import numpy as np

    #calculate band ratio control signal
    Sample, Replicate, Strain, Ab_ratio = [],[],[],[]
    #process data for each sample
    for sample in blot_data.Sample.unique():
        this_sample = blot_data.loc[blot_data.Sample == sample]
        #split this_sample into two tables depending on the antibody used
        this_test_sample = this_sample.loc[this_sample.Antibody != control_Ab]
        this_control_sample = this_sample.loc[this_sample.Antibody == control_Ab]
        #process data for each strain
        for strain in this_test_sample.Strain.unique():
            this_strain_test = this_test_sample.loc[this_test_sample.Strain == strain]
            this_strain_control = this_control_sample.loc[this_control_sample.Strain == strain]
            #process data for each replicate
            for replicate in this_strain_test.Replicate.unique():
                Sample.append(sample)
                Replicate.append(replicate)
                Strain.append(strain)
                ratio = this_strain_test.loc[this_strain_test.Replicate == replicate].Signal.values[0] / this_strain_control.loc[this_strain_control.Replicate == replicate].Signal.values[0]
                Ab_ratio.append(ratio)

    ratios = pd.DataFrame({'Sample': Sample,'Replicate':Replicate,'Strain':Strain, 'Ab_ratio':Ab_ratio})

    #calculate rexpression ratio opt to wt strain
    Sample, Replicate, fold_change = [],[],[]

    #process data for each sample
    for sample in ratios.Sample.unique():
        this_sample = ratios.loc[ratios.Sample == sample]
        #split this_sample into two tables depending on strain
        this_test_sample = this_sample.loc[this_sample.Strain != control_strain]
        this_control_sample = this_sample.loc[this_sample.Strain == control_strain]
        #process data for each replicate
        for replicate in this_test_sample.Replicate.unique():
            Sample.append(sample)
            Replicate.append(replicate)
            fo_ch = this_test_sample.loc[this_test_sample.Replicate == replicate].Ab_ratio.values[0] / np.mean(this_control_sample.Ab_ratio)
            fold_change.append(fo_ch)
    return pd.DataFrame({'Sample':Sample,'Replicate':Replicate,'protein_fold_change':fold_change})
