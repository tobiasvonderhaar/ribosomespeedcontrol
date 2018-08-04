

def process_qPCR(qPCR,control_primer='LEU3',control_strain='wt'):

    import pandas as pd
    import numpy as np

    ###############average technical replicates
    Clone,Primer,Replicate,Ct = [],[],[],[]
    #process each clone
    for clone in qPCR.clone.unique():
        this_clone = qPCR.loc[qPCR.clone == clone]
        #process each primer for this clone
        for primer in this_clone.primers.unique():
            this_primer = this_clone.loc[this_clone.primers == primer]
            #process each biological replicate for this clone-primer combination
            #create a single line with averaged Ct values
            for replicate in this_primer.replicate.unique():
                Clone.append(clone)
                Primer.append(primer)
                Replicate.append(replicate)
                Ct.append(np.mean(this_primer.loc[this_primer.replicate == replicate].Ct))

    tech_averaged = pd.DataFrame({'Clone':Clone,'Primer':Primer,'Replicate':Replicate,'Ct':Ct})

    ###############for each non-LEU3 Ct value, calculate the delta Ct with the matching LEU3 Ct value
    Clone,Replicate, Primer, Ct_LEU3 = [],[],[],[]
    #go through each clone
    for clone in tech_averaged.Clone.unique():
        this_clone = tech_averaged.loc[tech_averaged.Clone == clone]
        #go through each primer for this clone (required because multiple primers for wt)
        for primer in tech_averaged.Primer.unique():
            if primer != 'LEU3':
                this_specific_set = this_clone.loc[this_clone.Primer == primer]
                this_control_set = this_clone.loc[this_clone.Primer == control_primer]
            #go through each replicate
            for replicate in this_specific_set.Replicate.unique():
                if primer != 'LEU3':
                    Clone.append(clone)
                    Primer.append(primer)
                    Replicate.append(replicate)
                    dCt = this_control_set.loc[this_control_set.Replicate == replicate].Ct.values[0] - this_specific_set.loc[this_specific_set.Replicate == replicate].Ct.values[0]
                    Ct_LEU3.append(dCt)

    delta_Ct = pd.DataFrame({'Clone':Clone,'Primer':Primer,'Replicate':Replicate,'deltaCt_LEU3':Ct_LEU3})

    ###############calculate the ddCt value for optimised genes and the average of the corresponding wt control
    Clone, Replicate, deltadeltaCt_wt = [],[],[]
    #split delta_Ct into two tables for clones and wt control
    delta_clones = delta_Ct.loc[delta_Ct.Clone != control_strain]
    delta_wt = delta_Ct.loc[delta_Ct.Clone == control_strain]
    #go through each clone in the clones table
    for clone in delta_clones.Clone.unique():
        this_clone = delta_clones.loc[delta_clones.Clone == clone]
        this_wt = delta_wt.loc[delta_wt.Primer == this_clone.Primer.iloc[0]]
        #go through each replicate for this clone
        for replicate in this_clone.Replicate.unique():
            Clone.append(clone)
            Replicate.append(replicate)
            ddCt = np.mean(this_wt.deltaCt_LEU3) - this_clone.loc[this_clone.Replicate == replicate].deltaCt_LEU3.values[0]
            deltadeltaCt_wt.append(ddCt)

    #convert ddCt to fold change opt vs wt
    fold_change_wt = 2**-np.array(deltadeltaCt_wt)
    return pd.DataFrame({'Clone':Clone,'Replicate':Replicate,'ddCt_wt':deltadeltaCt_wt, 'fold_change_wt': fold_change_wt})
