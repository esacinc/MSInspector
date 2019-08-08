import os
import sys
import time
import pandas as pd
import numpy as np
import xml.etree.ElementTree as ET
import re
import subprocess

def select_rows_old(df,search_strings):
	# This function is to return the rows which contain all of the strings in search_strings
	unq,IDs = np.unique(df,return_inverse=True)
	unqIDs = np.searchsorted(unq,search_strings)
	df_out1 = df[((IDs.reshape(df.shape) == unqIDs[:,None,None]).any(-1)).all(0)]
	df_out2 = df[~((IDs.reshape(df.shape) == unqIDs[:,None,None]).any(-1)).all(0)]
	return df_out1, df_out2

def is_number(s):
	try:
		f = float(s)
		if f!=f or f == float('inf') or f == float('-inf'):
			return False
		return True
	except ValueError:
		return False

def infer_dataType(s):
	if s.isdigit():
		if int(s) == 0:
			data_type = 'zero integer'
		else:
			data_type = 'integer'
	elif is_number(s):
		if float(s) >= 0:
			data_type = 'number'
		else:
			data_type = 'negative number'
	else:
		data_type = 'text'
	return data_type

def string_finder(row, words, waived_list):
	# row is a type of Series.
	row1 = row.tolist()
	row2 = [v for u, v in enumerate(row1) if u not in waived_list]
	if any(word == field for field in row2 for word in words):
		return True
	return False

def select_rows(df, match, waived_col_id_list):
	ids = df.apply(string_finder, words=match, waived_list=waived_col_id_list, axis=1)
	df_out1 = df[ids]
	df_out2 = df[~ids]
	return df_out1, df_out2

def check_missing(listToCheck):
	# If status is True, that means listToCheck has '' or missing values exist.
	status = False
	for item in listToCheck:
		if item == '':
			status = True
			break
	return status

def locate_missing(valueList, colNameList, match, waived_list):
	outList = [colNameList[index] for index, value in enumerate(valueList) if (value in match) and (index not in waived_list)]
	return "Error", "Attribute", "Essential attribute(s) has(have) missing values, including "+'; '.join(outList)+'.'

def identify_uniProtKB_entryID(proteinName):
# This function is to extract uniProtKB_entryID from the protein name.
	if '|' in proteinName:
		tmp = proteinName.split('|')
		uniProtKB_entryID_tmp = tmp[1]
		# Judge whether uniProtKB_entryID is a legal uniProtKB entry ID based on its pattern using regular expression. Please refer to https://www.uniprot.org/help/accession_numbers
		pattern = re.compile('[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}')
		match = pattern.match(uniProtKB_entryID_tmp)
		if match:
			uniProtKB_entryID = uniProtKB_entryID_tmp
		else:
			uniProtKB_entryID = proteinName
	else:
		uniProtKB_entryID = proteinName
	return uniProtKB_entryID

def qcAnalysisGlobal(experiment_type, skyTsvDirList, fileNameList, required_col_dic, waived_col_dic, fileNameSkylineTmpTypeDic, col_dataType_dic):
	# Check all the *.tsv files in skyTsvDirList to make sure all of them exist,
	# Otherwise, it means some *.sky.zip files are not successfully transformed into *.tsv file via SkylineCMd and an error will be yielded.
	#failedList = [fileNameList[i] for i, item in enumerate(skyTsvDirList) if not os.path.isfile(item)]
	#if len(failedList) >0:
	#	print >> sys.stderr, "The file(s) which is(are): %s can't be opened by Skyline. Please check the version compatibility between Skyline document and Skyline program."%(", ".join(failedList))
	#	sys.exit(1)
	# Only keep the required columns based on required_col_dic
	if experiment_type in ['exp1', 'exp2']:
		fileNameSelceted = fileNameSkylineTmpTypeDic.keys()[0]
		skylineTemp_type =  fileNameSkylineTmpTypeDic[fileNameSelceted]
		required_col_list = required_col_dic[experiment_type][skylineTemp_type]
		waived_col_list = waived_col_dic[experiment_type][skylineTemp_type]
		col_dataType_forCheck_dic = col_dataType_dic[experiment_type][skylineTemp_type]
	else:
		required_col_list = required_col_dic[experiment_type]
		waived_col_list = waived_col_dic[experiment_type]
		col_dataType_forCheck_dic = col_dataType_dic[experiment_type]
	# touch the first skyFile
	dfTemplate = pd.read_csv(skyTsvDirList[0], sep='\t', header=0, converters={i: str for i in range(0, 100)})
	dfTemplate = dfTemplate[required_col_list]
	# Because the exported .tsv files are based on the skyline report template, they must have the same column information.
	# If one sky document file lacks one column data, the exported .tsv file will assign missing values to this this column.
	errorDf = pd.DataFrame(columns=['SkyDocumentName', 'IssueType', 'IssueSubtype', 'IssueReason']+list(dfTemplate.columns.values))
	normalDf = pd.DataFrame(columns=['SkyDocumentName']+list(dfTemplate.columns.values))
	peptide_excluded_in_Rscript_Df = pd.DataFrame(columns=['peptide', 'precursorCharge', 'isotopeLabelType', 'transition', 'uniProtKBID', 'proteinName', 'SkyDocumentName'])
	# Step 1: detect missing values "" without regard to the waived_col of the specific the experiment type
	search_strings = ['']
	#waived_col_list = waived_col_dic[experiment_type]
	for i, skyFileDir in enumerate(skyTsvDirList):
		# Read the *.tsv file into dataframe and keep all the values in the format of string. 
		df = pd.read_csv(skyFileDir, sep='\t', header=0, converters={j: str for j in range(0, 100)})
		df = df[required_col_list]
		# Add 'SkyDocumentName' into df as the first column
		col_name1 = df.columns.tolist()
		col_name1.insert(0, 'SkyDocumentName')
		df = df.reindex(columns=col_name1)
		df['SkyDocumentName'] = fileNameList[i]
		waived_col_id_list = [s for s, item in enumerate(df.columns.tolist()) if item in waived_col_list]
		#print waived_col_id_list
		dfTmp1, dfTmp2 = select_rows(df, search_strings, waived_col_id_list)
		#print dfTmp1.shape
		removedPeptideList = []
		if dfTmp1.shape[0] > 0:
			col_name = dfTmp1.columns.tolist()
			col_name.insert(1,'IssueType')
			col_name.insert(2,'IssueSubtype')
			col_name.insert(3,'IssueReason')
			dfTmp1 = dfTmp1.reindex(columns=col_name)
			waived_col_id_list_update = [s for s, item in enumerate(dfTmp1.columns.tolist()) if item in waived_col_list]
			#print dfTmp1
			#print waived_col_id_list_update
			resultTmp= dfTmp1.apply(lambda row: locate_missing([row[colName] for colName in dfTmp1.columns.values], list(dfTmp1.columns.values), search_strings, waived_col_id_list_update), axis=1)
			# resultTmp is a Series
			indexTmp1 = 0
			# suppress SettingWithCopyWarning
			pd.options.mode.chained_assignment = None
			for indexTmp, itemTmp in resultTmp.iteritems():
				dfTmp1['IssueType'][dfTmp1.index[indexTmp1]] = itemTmp[0]
				dfTmp1['IssueSubtype'][dfTmp1.index[indexTmp1]] = itemTmp[1]
				dfTmp1['IssueReason'][dfTmp1.index[indexTmp1]] = itemTmp[2]
				indexTmp1 = indexTmp1 + 1
			errorDf = pd.concat([errorDf, dfTmp1],  ignore_index=True)
			# Deduplicate the 'PeptideModifiedSequence' with the PrecursorCharge
			peptideList = list(set(dfTmp1['PeptideModifiedSequence']+'$$$$'+dfTmp1['PrecursorCharge']))
			#peptideList = list(set(dfTmp1['PeptideModifiedSequence']))
			removedPeptideList = removedPeptideList + peptideList
		# The peptides in dfTmp1[''] will be removed from dfTmp2, because the peptide information is 
		# incomplete and these peptides should not be exported into *.tsv file for downstream R code.
		# Pay attention: If part of data for one peptide sequence have missing value and the rest part don't have missing value, all the rows of this peptide with specific precursor charge will be removed.
		#
		#for removedPeptide in removedPeptideList:
		#	removedPeptide1 = removedPeptide.split('$$$$')[0]
		#	removedPeptide2 = removedPeptide.split('$$$$')[1]
		#	dfTmp2 = dfTmp2[~((dfTmp2['PeptideModifiedSequence'] == removedPeptide1) & (dfTmp2['PrecursorCharge'] == removedPeptide2))]
		#
		# Look into dfTmp2 to detect potential missing values and incorrect data type for some columns.
		fileNameSelceted = fileNameList[i]
		skylineTemp_type =  fileNameSkylineTmpTypeDic[fileNameSelceted]
		# Traverse PeptideModifiedSequence in dfTmp2
		PeptideModifiedSequenceList = set(dfTmp2['PeptideModifiedSequence'].tolist())
		for peptideSeq in PeptideModifiedSequenceList:
			dfTmp2_1 = dfTmp2[dfTmp2['PeptideModifiedSequence'] == peptideSeq]
			# Traverse PrecursorCharge in dfTmp2_1
			PrecursorChargeList = set(dfTmp2_1['PrecursorCharge'].tolist())
			for precursorCharge in PrecursorChargeList:
				dfTmp2_2 = dfTmp2_1[dfTmp2_1['PrecursorCharge'] == precursorCharge]
				#print dfTmp2_2
				# Step 1: If experiment_type is 'exp1', check whether there are missing values in ISSpike or PeptideConcentrationIS, Concentration or PeptideConcentration and MultiplicationFactor
				if experiment_type == 'exp1' and skylineTemp_type == 'old':
					iSSpikeStatus = check_missing(dfTmp2_2['ISSpike'].tolist())
					peptideConcentrationISStatus = check_missing(dfTmp2_2['PeptideConcentrationIS'].tolist())
					concentrationStatus = check_missing(dfTmp2_2['Concentration'].tolist())
					peptideConcentrationStatus = check_missing(dfTmp2_2['PeptideConcentration'].tolist())
					multiplicationFactorStatus = check_missing(dfTmp2_2['MultiplicationFactor'].tolist())
					if (iSSpikeStatus and peptideConcentrationISStatus) or (concentrationStatus and (peptideConcentrationStatus or multiplicationFactorStatus)):
						errorDfTmp = pd.DataFrame(columns=['SkyDocumentName', 'IssueType', 'IssueSubtype', 'IssueReason']+list(dfTemplate.columns.values))
						skyDocumentName = dfTmp2_2['SkyDocumentName'].tolist()[0]
						issueType = "Error"
						issueSubtype = "Attribute"
						proteinName = dfTmp2_2['ProteinName'].tolist()[0]
						columnsWithIssue = []
						if iSSpikeStatus and peptideConcentrationISStatus:
							columnsWithIssue.append('ISSpike or PeptideConcentrationIS')
						if concentrationStatus and (peptideConcentrationStatus or multiplicationFactorStatus):
							columnsWithIssue.append('Concentration or PeptideConcentration and MultiplicationFactor')
						issueReason = "Essential attribute(s) has(have) missing values, including "+'; '.join(columnsWithIssue)+'.'
						errorDfTmp.loc[len(errorDfTmp)] = [skyDocumentName, issueType, issueSubtype, issueReason, proteinName, peptideSeq, '', precursorCharge] + ['']*(errorDfTmp.shape[1]-8)
						errorDf = pd.concat([errorDf, errorDfTmp],  ignore_index=True)
						# Deduplicate the 'PeptideModifiedSequence' with the PrecursorCharge
						peptideList = list(set(dfTmp2_2['PeptideModifiedSequence']+'$$$$'+dfTmp2_2['PrecursorCharge']))
						#peptideList = list(set(dfTmp2_2['PeptideModifiedSequence']))
						for item in peptideList:
							if item not in removedPeptideList:
								removedPeptideList.append(item)
				# Step 2: Check the data type of some required columns. This need to be done later.
				col_dataType_forCheck_list = col_dataType_forCheck_dic.keys()
				col_with_wrong_dataType_list = []
				for col_name in col_dataType_forCheck_list:
					col1 = set([infer_dataType(item) for item in dfTmp2_2[col_name].tolist() if item != ''])
					col2 = col_dataType_forCheck_dic[col_name]
					if all([item in col2  for item in col1]):
						pass
					else:
						# This column has unqualified data type.
						unqualified_dataType = [item for item in col1 if item not in col2]
						unqualified_dataType_infor = 'Attribute "' + col_name + '" is annotated using ' + ', '.join(unqualified_dataType)
						col_with_wrong_dataType_list.append(unqualified_dataType_infor)
				if len(col_with_wrong_dataType_list) > 0:
					errorDfTmp = pd.DataFrame(columns=['SkyDocumentName', 'IssueType', 'IssueSubtype', 'IssueReason']+list(dfTemplate.columns.values))
					skyDocumentName = dfTmp2_2['SkyDocumentName'].tolist()[0]
					issueType = "Error"
					issueSubtype = "Attribute"
					proteinName = dfTmp2_2['ProteinName'].tolist()[0]
					issueReason = "Essential attribute(s) has(have) unqualified data types: " + '; '.join(col_with_wrong_dataType_list)+'.'
					errorDfTmp.loc[len(errorDfTmp)] = [skyDocumentName, issueType, issueSubtype, issueReason, proteinName, peptideSeq, '', precursorCharge] + ['']*(errorDfTmp.shape[1]-8)
					errorDf = pd.concat([errorDf, errorDfTmp],  ignore_index=True)
					# Deduplicate the 'PeptideModifiedSequence' with the PrecursorCharge
					peptideList = list(set(dfTmp2_2['PeptideModifiedSequence']+'$$$$'+dfTmp2_2['PrecursorCharge']))
					for item in peptideList:
						if item not in removedPeptideList:
							removedPeptideList.append(item)
		for removedPeptide in removedPeptideList:
			removedPeptide1 = removedPeptide.split('$$$$')[0]
			removedPeptide2 = removedPeptide.split('$$$$')[1]
			dfTmp2 = dfTmp2[~((dfTmp2['PeptideModifiedSequence'] == removedPeptide1) & (dfTmp2['PrecursorCharge'] == removedPeptide2))]
		# Add the dfTmp2 into normalDf
		normalDf = pd.concat([normalDf, dfTmp2], ignore_index=True)
		# Extract the information of the removed peptides.
		for removedPeptide in removedPeptideList:
			removedPeptide1 = removedPeptide.split('$$$$')[0]
			removedPeptide2 = removedPeptide.split('$$$$')[1]
			dfTmp3 = df[(df['PeptideModifiedSequence']==removedPeptide1) & (df['PrecursorCharge']==removedPeptide2)]
			dfTmp3['fragment_ion_complete'] = dfTmp3['FragmentIon']+" ("+dfTmp3['ProductCharge']+"+)"
			peptide_list = dfTmp3['PeptideModifiedSequence'].unique()
			for input_peptide_sequence in peptide_list:
				dfTmp4 = dfTmp3[(dfTmp3['PeptideModifiedSequence']==input_peptide_sequence)]
				protein_list = dfTmp3['ProteinName'].unique()
				protein_uniProtID_list = [identify_uniProtKB_entryID(protein_tmp) for protein_tmp in protein_list]
				for indexLabel, protein_tmp in enumerate(protein_list):
					protein_uniProtID = protein_uniProtID_list[indexLabel]
					dfTmp5 = dfTmp4[(dfTmp4['ProteinName']==protein_tmp)]
					for precursorchargeTmp in dfTmp5['PrecursorCharge'].unique():
						dfTmp6 = dfTmp5[(dfTmp5['PrecursorCharge']==precursorchargeTmp)]
						isotopeLabelType_list = dfTmp6['IsotopeLabelType'].unique()
						isotopeLabelType_list.sort()
						isotopeLabelTypeTmp = '|'.join(isotopeLabelType_list)
						transitionTmp = []
						for isotopeLabelTypeSubtmp in isotopeLabelType_list:
							fragmentIon_list = dfTmp6[(dfTmp6['IsotopeLabelType']==isotopeLabelTypeSubtmp)]['fragment_ion_complete'].unique()
							fragmentIon_list.sort()
							transitionTmp.append(':'.join([isotopeLabelTypeSubtmp, '|'.join(fragmentIon_list)]))
						transitionTmp = ';'.join(transitionTmp)
						peptide_excluded_in_Rscript_Df_tmp = pd.DataFrame({'peptide':[input_peptide_sequence], 'precursorCharge':[precursorchargeTmp], 'isotopeLabelType':[isotopeLabelTypeTmp], 'transition':[transitionTmp], 'uniProtKBID':[protein_uniProtID], 'proteinName':[protein_tmp], 'SkyDocumentName':[fileNameList[i]]}, columns=['peptide', 'precursorCharge', 'isotopeLabelType', 'transition', 'uniProtKBID', 'proteinName', 'SkyDocumentName'])
						peptide_excluded_in_Rscript_Df = pd.concat([peptide_excluded_in_Rscript_Df, peptide_excluded_in_Rscript_Df_tmp],  ignore_index=True)
	return errorDf, normalDf, peptide_excluded_in_Rscript_Df

def detectIS(skyFileDir, fileName, experiment_type, error_report_path, errorDfColNumber):
	# Parse the *.sky file whose format is XML.
	for event, elem in ET.iterparse(skyFileDir):
			if event == 'end':
				if elem.tag == 'peptide_modifications':
					internal_standard= elem.get('internal_standard', default=None)
					if internal_standard == 'none':
					# This means that the user doesn't set the Internal Standard Type when preparing data for upload.
					# This situation works for experiment 1, 2, 3, 4 and 5.
					 	#print >> sys.stderr, "Internal_standard value in the peptide_modifications underneath peptide_settings of the *.sky file of %s is unset. Please check it."%(skyFileDir)
						if experiment_type=='exp1' or experiment_type=='exp2':
							errorInfor = '\t'.join([os.path.basename(fileName), 'Error', 'Internal standard','Internal standard type in the peptide_modifications underneath peptide_settings is set to be none.']+['']*(errorDfColNumber-4))+'\n'
						else:
							errorInfor = '\t'.join([os.path.basename(fileName), 'Error', 'Internal standard','Internal standard type in the peptide_modifications underneath peptide_settings is set to be none. Please set it to be heavy.']+['']*(errorDfColNumber-4))+'\n'
						with open(error_report_path, 'a') as outfTmp:
							outfTmp.write(errorInfor)
						internal_standard_type = 'none'
					elif internal_standard is None:
						# This means that internal_standard is set to be heavy by default.
						internal_standard_type = 'heavy'
					else:
						internal_standard_type = internal_standard
					# In exp3,exp4 and exp5, the default internal standard should always be heavy. A check will be performed in experiment3_qc.R, experiment4_qc.R and experiment5_qc.R
					#if experiment_type=='exp3' or experiment_type=='exp4' or experiment_type == 'exp5':
						# If internal_standard is not heavy, that means the user forgets to set it to be a correct value when preparing data for upload. 
					#	if internal_standard != 'heavy':
					#		errorInfor = '\t'.join([os.path.basename(fileName), 'Error', 'Internal standard','Internal standard type in the peptide_modifications underneath peptide_settings is incorrect. Please set it to be heavy']+['']*(errorDfColNumber-4))+'\n'
					#	with open(error_report_path, 'a') as outfTmp:
					#		outfTmp.write(errorInfor)
			elem.clear()
	return internal_standard_type

def qcAnalysisRcode(experiment_type, error_report_path, dataset_path, fileList_path, mypeptideType_file_path, RscriptBinary, rScript, plot_output, plot_output_dir):
	if experiment_type == 'exp2' or experiment_type == 'exp5' or experiment_type == 'exp3' or experiment_type == 'exp4':
		#os.system('"%s" %s %s %s %s %s >> %s'%(RscriptBinary, rScript, dataset_path, fileList_path, plot_output, plot_output_dir, error_report_path))
		subprocess.call([RscriptBinary, rScript, dataset_path, fileList_path, str(plot_output), plot_output_dir, ">>", error_report_path])
	if experiment_type == 'exp1':
		#os.system('"%s" %s %s %s %s %s %s >> %s'%(RscriptBinary, rScript, dataset_path, fileList_path, plot_output, plot_output_dir, mypeptideType_file_path, error_report_path))
		subprocess.call([RscriptBinary, rScript, dataset_path, fileList_path, str(plot_output), plot_output_dir, mypeptideType_file_path, ">>", error_report_path])