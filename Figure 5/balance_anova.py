def balance_anova(luc_data, control='C', test_levels=[0.05,0.01,0.001]):

    import pandas as pd
    from statsmodels.stats.multicomp import pairwise_tukeyhsd
    from statsmodels.stats.multicomp import MultiComparison

    #if the variables for storing ANOVA results are already in locals(), delete them
    if 'outputI' in locals():
        del outputI
    if 'outputII' in locals():
        del outputII

    #conduct statistical tests for each of the significnce levels desginaed in the variable 'levels'
    for level in test_levels:

        #calculte test statistics for RLucI
        RLucI_data = luc_data.loc[luc_data['Factor 2'] == 'RLucI']
        mcI = MultiComparison(RLucI_data['norm2'], RLucI_data['Factor 1'])
        resultI = mcI.tukeyhsd(alpha = level)
        #convert results to data frame
        dfI = pd.DataFrame(data=resultI._results_table.data[1:], columns=resultI._results_table.data[0])
        #go through each row of dfI, and test if control exists in comparison group2. If so swap group2 and group1
        for row in range(dfI.shape[0]):
            if dfI.iloc[row].group2 == control:
                dfI.iloc[row, dfI.columns.get_loc('group2')] = dfI.iloc[row].group1
                dfI.iloc[row, dfI.columns.get_loc('group1')]  = control
        #retain only rows comparing values to the control
        dfI = dfI.loc[dfI['group1']==control]
        #retain only columns identifying the factor compared to control, and the result of the test
        dfI = dfI.loc[:,['group2','reject']]
        dfI = dfI.rename(index=str, columns = {'reject':str(level)})
        #either record the result as outputI, or merge the results with an existing outputI
        if 'outputI' in locals():
            outputI = pd.merge(outputI,dfI,on='group2')
        else:
            outputI = dfI

        #test statistics for RLucII
        RLucII_data = luc_data.loc[luc_data['Factor 2'] == 'RLucII']
        mcII = MultiComparison(RLucII_data['norm2'], RLucII_data['Factor 1'])
        resultII = mcII.tukeyhsd(alpha = level)
        dfII = pd.DataFrame(data=resultII._results_table.data[1:], columns=resultII._results_table.data[0])
        for row in range(dfII.shape[0]):
            if dfII.iloc[row].group2 == control:
                dfII.iloc[row, dfII.columns.get_loc('group2')] = dfII.iloc[row].group1
                dfII.iloc[row, dfII.columns.get_loc('group1')]  = control
        dfII = dfII.loc[dfII['group1']==control]
        dfII = dfII.loc[:,['group2','reject']]
        dfII = dfII.rename(index=str, columns = {'reject':str(level)})
        if 'outputII' in locals():
            outputII = pd.merge(outputII,dfII,on='group2')
        else:
            outputII = dfII

    outputI = outputI.rename(index = str, columns={'group2':'FlucI'})
    outputII = outputII.rename(index = str, columns={'group2':'FlucII'})

    print(outputI)
    print(outputII)
