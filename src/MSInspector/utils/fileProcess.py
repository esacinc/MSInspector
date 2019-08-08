import os
import pandas as pd
import random

def judgeColumnNonblank(df_test, columnName):
    # The columnName must be in the column list of df_test
    sampledListTest = df_test[columnName].values.tolist()
    if df_test.shape[0] >= 50:
        sampledListTest = random.sample(sampledListTest, 50)
    columnNonblankStatus = False
    if any([item != '' for item in sampledListTest]):
        columnNonblankStatus = True
    return columnNonblankStatus

def judgeSampleGroupInExp2(df_test):
    sampledListTest = df_test['SampleGroup'].values.tolist()
    SampleGroupInExp2Status = False
    if all([item.lower()[:1] in ['l', 'm', 'h'] for item in sampledListTest if item != '']):
        SampleGroupInExp2Status = True
    return SampleGroupInExp2Status

def judgeInteger(df_test, columnName):
    sampledListTest = df_test[columnName].values.tolist()
    integerStatus = False
    if all([item.isdigit() for item in sampledListTest if item != '']):
        integerStatus = True
    return integerStatus

def headerModify(inf_tmp, outf_tmp, experiment_type_tmp, skyzip_file_dir_baseName, fileNameInferredExpStatusDic, fileNameEmptyStatusDic):
    # This function is to remove space in the column name at the first row of inf_tmp.
    outf = open(outf_tmp, 'w')
    with open(inf_tmp, 'r') as inf:
        line = inf.readline()
        outf.write(line.replace(' ', ''))
        line = inf.readline()
        while line != '':
            # replace #N/A with '', so that when R read this file, it will transform '' into NA automatically.
            outf.write(line.replace('#N/A', ''))
            line = inf.readline()
    outf.close()
    os.remove(inf_tmp)
    # Inferred experiment type from outf_tmp.
    # Read outf_tmp into a dataframe and keep all the values in the format of string.
    df_tmp = pd.read_csv(outf_tmp, sep='\t', header=0, converters={i: str for i in range(0, 100)})
    if df_tmp.shape[0] == 0:
        fileNameEmptyStatusDic.update({skyzip_file_dir_baseName:True})
        fileNameInferredExpStatusDic.update({skyzip_file_dir_baseName:None})
        skylineTemp_type = None
    else:
        fileNameEmptyStatusDic.update({skyzip_file_dir_baseName:False})
        # Start to infer the experiment type. 
        exp3Stutus = judgeColumnNonblank(df_tmp, 'Exp3SampleGroup')
        exp4Stutus = judgeColumnNonblank(df_tmp, 'Exp4SampleGroup')
        exp2NewTemplateSatuts = judgeColumnNonblank(df_tmp, 'Exp2SampleGroup')
        dayStatus = judgeColumnNonblank(df_tmp, 'Day')
        analyteConcentrationStatus = judgeColumnNonblank(df_tmp, 'AnalyteConcentration')
        internalStandardConcentrationStatus = judgeColumnNonblank(df_tmp, 'InternalStandardConcentration')
        if exp3Stutus:
            inferredExpType = 'exp3'
            skylineTemp_type = None
        elif exp4Stutus:
            inferredExpType = 'exp4'
            skylineTemp_type = None
        elif (not exp2NewTemplateSatuts) and (not exp3Stutus) and (not exp4Stutus) and dayStatus:
            inferredExpType = 'exp5'
            skylineTemp_type = None
        elif exp2NewTemplateSatuts:
            inferredExpType = 'exp2'
            skylineTemp_type = 'new'
        elif (not exp2NewTemplateSatuts) and (not exp3Stutus) and (not exp4Stutus) and analyteConcentrationStatus and internalStandardConcentrationStatus:
            inferredExpType = 'exp1'
            skylineTemp_type = 'new'
        else:
            # The skyline file must be from exp1 and exp2 of old templates which have the same columns names.
            # Infer the experiment type base the content of column SampleGroup and Concentration.
            sampleGroupStatus = judgeSampleGroupInExp2(df_tmp)
            integerStatus = judgeInteger(df_tmp, 'Concentration')
            if sampleGroupStatus and integerStatus:
                inferredExpType = 'exp2'
                skylineTemp_type = 'old'
            else:
                inferredExpType = 'exp1'
                skylineTemp_type = 'old'
        if inferredExpType != experiment_type_tmp:
            fileNameInferredExpStatusDic.update({skyzip_file_dir_baseName:True})
        else:
            fileNameInferredExpStatusDic.update({skyzip_file_dir_baseName:False})
    # Judge is the skyline template old or new if experiment_type_tmp is exp1 or exp2
    #if experiment_type_tmp in ["exp1", "exp2"]:
    #    if df_tmp.shape[0] < 10:
    #        sampledList1 = df_tmp['AnalyteConcentration'].values.tolist()
    #        sampledList2 = df_tmp['InternalStandardConcentration'].values.tolist()
    #    else:
    #        sampledList1 = random.sample(df_tmp['AnalyteConcentration'].values.tolist(), 10)
    #        sampledList2 = random.sample(df_tmp['InternalStandardConcentration'].values.tolist(), 10)
    #    if any([item != '' for item in sampledList1]) and any([item != '' for item in sampledList2]):
    #        skylineTemp_type = 'new'
    #    else:
    #        skylineTemp_type = 'old'
    #else:
    #    skylineTemp_type = None
    return skylineTemp_type