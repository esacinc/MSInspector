import pandas as pd
import os
import base64
import re
import sys

def deduplicateList(list_input):
	list_ouput = []
	for item in list_input:
		if item not in list_ouput:
			list_ouput.append(item)
	return list_ouput

def countPeptide(inputList):
	# The elements in inputList is (peptideSeq, precursorCharge) or (peptideSeq)
	peptideList = []
	if len(inputList) == 0:
		countOutput = 0
	else:
		for item in inputList:
			if isinstance(item, tuple):
				peptide = item[0]
			else:
				peptide = item
			if peptide not in peptideList:
				peptideList.append(peptide)
		countOutput = len(peptideList)
	return countOutput

def png2base64String(pngsrc):
	with open(pngsrc, "rb") as imageFile:
		pngString = base64.b64encode(imageFile.read())
	return pngString

def string2escape(inputstr):
	outstr = ''
	for i in inputstr:
		if i in ['[', ']', '+', '.']:
			i = '\\'+ i
		outstr = outstr + i
	return outstr

def locateFigure(peptideSeq, precursorCharge, fileList, experiment_type, dirName):
	if experiment_type == 'exp1':
		# There will be three figures plotting in the html.
		pattern1 = re.compile(r"^%s_linear_(.*)_ResponseCurveQuery.response_curve.png$"%(string2escape(peptideSeq))) 
		pattern2 = re.compile(r"^%s_log_(.*)_ResponseCurveQuery.response_curve.png$"%(string2escape(peptideSeq)))
		pattern3 = re.compile(r"^%s_residual_(.*)_ResponseCurveQuery.response_curve.png$"%(string2escape(peptideSeq)))
		f1 = ''
		f2 = ''
		f3 = ''
		for fileTmp in fileList:
			match = pattern1.match(fileTmp)
			if match:
				f1 = match.group()
				break
		for fileTmp in fileList:
			match = pattern2.match(fileTmp)
			if match:
				f2 = match.group()
				break
		for fileTmp in fileList:
			match = pattern3.match(fileTmp)
			if match:
				f3 = match.group()
				break
		f1 = os.path.join(dirName, f1)
		f2 = os.path.join(dirName, f2)
		f3 = os.path.join(dirName, f3)
		
		f1_string = 'data:image/png;base64,'+png2base64String(f1)
		f2_string = 'data:image/png;base64,'+png2base64String(f2)
		f3_string = 'data:image/png;base64,'+png2base64String(f3)
		return ([f1, f2, f3], [f1_string, f2_string, f3_string])
	if  experiment_type == 'exp2' or experiment_type == 'exp5' or experiment_type == 'exp3' or experiment_type == 'exp4':
		pattern1 = re.compile(r"^%s_%s_(.*).png$"%(string2escape(peptideSeq), string2escape(precursorCharge))) 
		f1 = ''
		for fileTmp in fileList:
			match = pattern1.match(fileTmp)
			if match:
				f1 = match.group()
				break
		f1 = os.path.join(dirName, f1)
		f1_string = 'data:image/png;base64,'+png2base64String(f1)
		return ([f1], [f1_string])

def locateTable(peptideSeq, precursorCharge, fileList, experiment_type, dirName):
	if experiment_type == 'exp1' or experiment_type == 'exp3':
		# There will be two tables shown in the html.
		if experiment_type == 'exp1':
			pattern1 = re.compile(r"^%s_(.*)_ResponseCurveAnalysis.LODCTable.tsv$"%(string2escape(peptideSeq)))
			pattern2 = re.compile(r"^%s_(.*)_ResponseCurveAnalysis.fitTable.tsv$"%(string2escape(peptideSeq)))
		else:
			pattern1 = re.compile(r"^%s_%s_(.*)_ave_values_for_spike_levels.tsv$"%(string2escape(peptideSeq), string2escape(precursorCharge)))
			pattern2 = re.compile(r"^%s_%s_(.*)_summary_table.tsv$"%(string2escape(peptideSeq), string2escape(precursorCharge)))
		f1 = ''
		f2 = ''
		for fileTmp in fileList:
			match = pattern1.match(fileTmp)
			if match:
				f1 = match.group()
				break
		for fileTmp in fileList:
			match = pattern2.match(fileTmp)
			if match:
				f2 = match.group()
				break
		f1 = os.path.join(dirName, f1)
		f2 = os.path.join(dirName, f2)
		f1_df = pd.read_csv(f1, sep='\t', header=0, converters={i: str for i in range(0, 100)})
		f2_df = pd.read_csv(f2, sep='\t', header=0, converters={i: str for i in range(0, 100)})
		return([f1, f2], [f1_df, f2_df])
	if experiment_type == 'exp2' or experiment_type == 'exp4' or experiment_type == 'exp5':
		pattern1 = re.compile(r"^%s_%s_(.*)_CV_results.tsv$"%(string2escape(peptideSeq), string2escape(precursorCharge)))
		f1 = ''
		for fileTmp in fileList:
			match = pattern1.match(fileTmp)
			if match:
				f1 = match.group()
				break
		f1 = os.path.join(dirName, f1)
		f1_df = pd.read_csv(f1, sep='\t', header=0, converters={i: str for i in range(0, 100)})
		return ([f1], [f1_df])

def peptide_infor_parse(is_infor_file, peptide_infor_file, peptide_excluded_in_Rscript_infor_file, assayFileList, assayInforDic, peptideTrackDic, peptideSeqChargeIsotopeDic, experiment_type):
	is_infor_df = pd.read_csv(is_infor_file, sep='\t', header=0, converters={i: str for i in range(0, 100)})
	#print is_infor_df
	peptide_infor_df = pd.read_csv(peptide_infor_file, sep='\t', header=0, converters={i: str for i in range(0, 100)})
	peptide_excluded_in_Rscript_infor_df = pd.read_csv(peptide_excluded_in_Rscript_infor_file, sep='\t', header=0, converters={i: str for i in range(0, 100)})
    
    # For the peptides in peptide_excluded_in_Rscript_infor_file.
	for item in peptide_excluded_in_Rscript_infor_df.SkyDocumentName.unique():
		assayFileList.append(item)
		assayInforDic.update({item:{}})
		peptideTrackDic.update({item:[]})
		peptideSeqChargeIsotopeDic.update({item:{}})
		# Add status, protein number, peptide number, precursor number, transition number, internal_standard, experimentType into assayInforDic[item]
		assayInforDic[item]['status'] = False
		peptide_infor_df_tmp = peptide_excluded_in_Rscript_infor_df[peptide_excluded_in_Rscript_infor_df['SkyDocumentName'] == item]
		assayInforDic[item]['proteinNum'] = len(peptide_infor_df_tmp['proteinName'].unique())
		assayInforDic[item]['peptideNum'] = len(peptide_infor_df_tmp['peptide'].unique())
		precursorNum = 0
		transitionNum = 0
		for index, row in peptide_infor_df_tmp.iterrows():
			precursorNum = precursorNum + len(row['isotopeLabelType'].split('|'))
			peptide = row['peptide']
			precursorCharge = row['precursorCharge']
			isotopeLableSet = row['isotopeLabelType'].split('|')
			try:
				peptideSeqChargeIsotopeDic[item][peptide][precursorCharge] = isotopeLableSet
			except KeyError:
				peptideSeqChargeIsotopeDic[item].update({peptide:{precursorCharge:isotopeLableSet}})
			if experiment_type == 'exp1':
				if (peptide) not in peptideTrackDic[item]:
					peptideTrackDic[item].append((peptide))
			elif experiment_type == 'exp2':
				if (peptide, precursorCharge) not in peptideTrackDic[item]:
					peptideTrackDic[item].append((peptide, precursorCharge))
			elif experiment_type == 'exp3':
				if (peptide, precursorCharge) not in peptideTrackDic[item]:
					peptideTrackDic[item].append((peptide, precursorCharge))
			elif experiment_type == 'exp4':
				if (peptide, precursorCharge) not in peptideTrackDic[item]:
					peptideTrackDic[item].append((peptide, precursorCharge))
			else:
				if (peptide, precursorCharge) not in peptideTrackDic[item]:
					peptideTrackDic[item].append((peptide, precursorCharge))
			for item1 in row['transition'].split(';'):
				transitionNum = transitionNum + len(item1.split(':')[1].split('|'))
		assayInforDic[item]['precursorNum'] = precursorNum
		assayInforDic[item]['transitionNum'] = transitionNum
		assayInforDic[item]['internalStandard'] = is_infor_df[is_infor_df['SkyDocumentName']==item]['internal_standard'].values[0]
		if experiment_type == 'exp1':
			assayInforDic[item]['experimentType'] = 'Response Curve'
		elif experiment_type == 'exp2':
			assayInforDic[item]['experimentType'] = 'Repeatability'
		elif experiment_type == 'exp3':
			assayInforDic[item]['experimentType'] = 'Selectivity'
		elif experiment_type == 'exp4':
			assayInforDic[item]['experimentType'] = 'Stability'
		else:
			assayInforDic[item]['experimentType'] = 'Endogenous'
	# Update the rest peptides in peptide_infor_file. These peptides (normal_data.tsv) are going to be analyzed in the R script. 
	for item in peptide_infor_df.SkyDocumentName.unique():
		peptide_infor_df_tmp = peptide_infor_df[peptide_infor_df['SkyDocumentName'] == item]
		if item not in assayFileList:
			assayFileList.append(item)
			assayInforDic.update({item:{}})
			peptideTrackDic.update({item:[]})
			peptideSeqChargeIsotopeDic.update({item:{}})
		assayInforDic[item]['status'] = True
		# Add protein number, peptide number, precursor number, transition number, internal_standard, experimentType into assayInforDic[item]
		try:
			assayInforDic[item]['proteinNum'] = assayInforDic[item]['proteinNum'] + len(peptide_infor_df_tmp['proteinName'].unique())
		except KeyError:
			assayInforDic[item]['proteinNum'] = len(peptide_infor_df_tmp['proteinName'].unique())
		try:
			assayInforDic[item]['peptideNum'] = assayInforDic[item]['peptideNum'] + len(peptide_infor_df_tmp['peptide'].unique())
		except KeyError:
			assayInforDic[item]['peptideNum'] = len(peptide_infor_df_tmp['peptide'].unique())
		precursorNum = 0
		transitionNum = 0
		for index, row in peptide_infor_df_tmp.iterrows():
			precursorNum = precursorNum + len(row['isotopeLabelType'].split('|'))
			peptide = row['peptide']
			precursorCharge = row['precursorCharge']
			isotopeLableSet = row['isotopeLabelType'].split('|')
			try:
				peptideSeqChargeIsotopeDic[item][peptide][precursorCharge] = isotopeLableSet
			except KeyError:
				peptideSeqChargeIsotopeDic[item].update({peptide:{precursorCharge:isotopeLableSet}})
			if experiment_type == 'exp1':
				if (peptide) not in peptideTrackDic[item]:
					peptideTrackDic[item].append((peptide))
			elif experiment_type == 'exp2':
				if (peptide, precursorCharge) not in peptideTrackDic[item]:
					peptideTrackDic[item].append((peptide, precursorCharge))
			elif experiment_type == 'exp3':
				if (peptide, precursorCharge) not in peptideTrackDic[item]:
					peptideTrackDic[item].append((peptide, precursorCharge))
			elif experiment_type == 'exp4':
				if (peptide, precursorCharge) not in peptideTrackDic[item]:
					peptideTrackDic[item].append((peptide, precursorCharge))
			else:
				if (peptide, precursorCharge) not in peptideTrackDic[item]:
					peptideTrackDic[item].append((peptide, precursorCharge))
			for item1 in row['transition'].split(';'):
				transitionNum = transitionNum + len(item1.split(':')[1].split('|'))
		try:
			assayInforDic[item]['precursorNum'] = assayInforDic[item]['precursorNum'] + precursorNum
		except:
			assayInforDic[item]['precursorNum'] = precursorNum
		try:
			assayInforDic[item]['transitionNum'] = assayInforDic[item]['transitionNum'] + transitionNum
		except:
			assayInforDic[item]['transitionNum'] = transitionNum
		if 'internalStandard' not in assayInforDic[item].keys():
			assayInforDic[item]['internalStandard'] = is_infor_df[is_infor_df['SkyDocumentName']==item]['internal_standard'].values[0]
		if 'experimentType' not in assayInforDic[item].keys():
			if experiment_type == 'exp1':
				assayInforDic[item]['experimentType'] = 'Response Curve'
			elif experiment_type == 'exp2':
				assayInforDic[item]['experimentType'] = 'Repeatability'
			elif experiment_type == 'exp3':
				assayInforDic[item]['experimentType'] = 'Selectivity'
			elif experiment_type == 'exp4':
				assayInforDic[item]['experimentType'] = 'Stability'
			else:
				assayInforDic[item]['experimentType'] = 'Endogenous'

def qc_report_infor_parse(qc_report_file, is_inferred_file, assayInforDic, peptideTrackDic, peptideOutputDic, assayFileList, experiment_type):
	dir_tmp = os.path.join(os.path.dirname(qc_report_file), 'figures_tables.tmp')
	file_list = os.listdir(dir_tmp)
	#print file_list
	qc_report_infor_df = pd.read_csv(qc_report_file, sep='\t', header=0, converters={i: str for i in range(0, 100)})
	is_inferred_infor_df = pd.read_csv(is_inferred_file, sep='\t', header=0, converters={i: str for i in range(0, 100)})
	for item in qc_report_infor_df.SkyDocumentName.unique():
		peptideOutputDic.update({item:{}})
		qc_report_infor_df_tmp = qc_report_infor_df[qc_report_infor_df['SkyDocumentName'] == item]
		is_inferred_tmp = is_inferred_infor_df[is_inferred_infor_df['SkyDocumentName'] == item]['internal_standard'].iloc[0]
		if not assayInforDic[item]['status']:
			# It means that all of the peptides in this sky document have missing attribute issues which will cause errors.
			assayInforDic[item]['isQuality'] = "Internal standard type can't be inferred. All the peptides have missing values or incorrect data types in some essential attributes."
			assayInforDic[item]['peptideSeqErrors'] = peptideTrackDic[item]
			assayInforDic[item]['peptideSeqWarnings'] = []
			assayInforDic[item]['peptideSeqWithoutIssues'] = []
			for index, row in qc_report_infor_df_tmp.iterrows():
				peptide = row['PeptideModifiedSequence']
				precursorCharge = row['PrecursorCharge']
				issueType = row['IssueType']
				issueSubtype = row['IssueSubtype']
				issueReason = row['IssueReason']
				if experiment_type == 'exp1':
					peptideTerm = (peptide)
				elif experiment_type == 'exp2':
					peptideTerm = (peptide, precursorCharge)
				elif experiment_type == 'exp3':
					peptideTerm = (peptide, precursorCharge)
				elif experiment_type == 'exp4':
					peptideTerm = (peptide, precursorCharge)
				else:
					peptideTerm = (peptide, precursorCharge)
				if peptideTerm not in peptideOutputDic[item].keys():
					peptideOutputDic[item].update({peptideTerm:{}})
				if 'issueReason' not in peptideOutputDic[item][peptideTerm].keys():
					peptideOutputDic[item][peptideTerm].update({'issueReason':[issueReason]})
				else:
					if issueReason not in peptideOutputDic[item][peptideTerm]['issueReason']:
						peptideOutputDic[item][peptideTerm]['issueReason'].append(issueReason)
		else:
			if "Internal standard" in qc_report_infor_df_tmp['IssueSubtype'].values and qc_report_infor_df_tmp['IssueType'].values[0] == 'Error':
				# It means that Internal Standard has issues for this skyDocumentName
				# It could be none or incorrect value
				assayInforDic[item]['isQuality'] = qc_report_infor_df_tmp[qc_report_infor_df_tmp['IssueSubtype']=="Internal standard"]['IssueReason'].values[0]+' Errors happen for all the peptides.'
			else:
				assayInforDic[item]['isQuality'] = 'Correct'
			# the rows whose column 'IssueSubtype' are 'Internal standard' need to be deleted from qc_report_infor_df_tmp
			qc_report_infor_df_tmp = qc_report_infor_df_tmp[qc_report_infor_df_tmp['IssueSubtype']!="Internal standard"]
			if assayInforDic[item]['isQuality'] == 'Correct':
				for index, row in qc_report_infor_df_tmp.iterrows():
					peptide = row['PeptideModifiedSequence']
					precursorCharge = row['PrecursorCharge']
					issueType = row['IssueType']
					issueSubtype = row['IssueSubtype']
					issueReason = row['IssueReason']
					if experiment_type == 'exp1':
						peptideTerm = (peptide)
					elif experiment_type == 'exp2':
						peptideTerm = (peptide, precursorCharge)
					elif experiment_type == 'exp3':
						peptideTerm = (peptide, precursorCharge)
					elif experiment_type == 'exp4':
						peptideTerm = (peptide, precursorCharge)
					else:
						peptideTerm = (peptide, precursorCharge)
					
					if peptideTerm not in peptideOutputDic[item].keys():
						peptideOutputDic[item].update({peptideTerm:{}})
						if issueType == 'Error':
							try:
								assayInforDic[item]['peptideSeqErrors'].append(peptideTerm)
							except KeyError:
								assayInforDic[item].update({'peptideSeqErrors': [peptideTerm]})
						else:
							# Pay attention: 
							# if issueSubtype != "Internal Standard spike peptide concentration" and experiment_type == "exp1"
							# This kind of warning will be treated as peptide without issue. But the issueReason: "All of the concentrations of the internal standard peptide are zero." will be shown in the output html.
							# So for the peptideTerm in assayInforDic[item]['peptideSeqWarnings'], the issueReason will evaluated later.
							try:
								assayInforDic[item]['peptideSeqWarnings'].append(peptideTerm)
							except KeyError:
								assayInforDic[item].update({'peptideSeqWarnings': [peptideTerm]})
							# Update the 'imageSrc', 'imagePNGbase64', 'table' in peptideOutputDic[item][peptideTerm]
							# For different experiments, the naming of the figures and tables are different. So process them individually.
							# Using regular expression to search file_list
							figure_src_list, base64String_list = locateFigure(peptide, precursorCharge, file_list, experiment_type, dir_tmp)
							table_src_list, df_list = locateTable(peptide, precursorCharge, file_list, experiment_type, dir_tmp)
							if 'imageSrc' not in peptideOutputDic[item][peptideTerm].keys():
								peptideOutputDic[item][peptideTerm].update({'imageSrc':figure_src_list})
								peptideOutputDic[item][peptideTerm].update({'imagePNGbase64':base64String_list})
								peptideOutputDic[item][peptideTerm].update({'table':df_list})
					try:
						peptideOutputDic[item][peptideTerm]['issueReason'].append(issueReason)
					except KeyError:
						peptideOutputDic[item][peptideTerm].update({'issueReason':[issueReason]})
			else:
				assayInforDic[item]['peptideSeqErrors'] = peptideTrackDic[item]
				# peptideTerm in assayInforDic[item]['peptideSeqErrors'] will not be updated into peptideOutputDic.
				# so that these peptideTerm will not be displayed in keptAssayFileErrorTable
				# however, if there are some other peptides with error information in qc_report_infor_df_tmp,
				# these peptides will be updated into peptideOutputDic
				# Pay attention, if isQuality in not correct, there won't be peptides with warnings in qc_report_infor_df_tmp
				# Keep parsing qc_report_infor_df_tmp to extract error information.
				for index, row in qc_report_infor_df_tmp.iterrows():
					peptide = row['PeptideModifiedSequence']
					precursorCharge = row['PrecursorCharge']
					issueType = row['IssueType']
					issueSubtype = row['IssueSubtype']
					issueReason = row['IssueReason']
					if experiment_type == 'exp1':
						peptideTerm = (peptide)
					elif experiment_type == 'exp2':
						peptideTerm = (peptide, precursorCharge)
					elif experiment_type == 'exp3':
						peptideTerm = (peptide, precursorCharge)
					elif experiment_type == 'exp4':
						peptideTerm = (peptide, precursorCharge)
					else:
						peptideTerm = (peptide, precursorCharge)

					if issueType == 'Error':
						if peptideTerm not in peptideOutputDic[item].keys():
							peptideOutputDic[item].update({peptideTerm:{}})
						try:
							peptideOutputDic[item][peptideTerm]['issueReason'].append(issueReason)
						except KeyError:
							peptideOutputDic[item][peptideTerm].update({'issueReason':[issueReason]})
				assayInforDic[item]['peptideSeqWarnings'] = []
				assayInforDic[item]['peptideSeqWithoutIssues'] = []
			
			# Find the peptides without issues and add them into assayInforDic and peptideOutputDic.
			petideTermWithoutIssues = []
			#print assayInforDic
			#print peptideTrackDic[item]
			#print peptideOutputDic[item]
			#print '&&&&&&&&&&&'
			for subitem in peptideTrackDic[item]:
				if 'peptideSeqErrors' in assayInforDic[item].keys():
					peptideSeqErrorsListTmp = set(assayInforDic[item]['peptideSeqErrors'])
				else:
					peptideSeqErrorsListTmp = set([])
				if 'peptideSeqWarnings' in assayInforDic[item].keys():
					peptideSeqWarningsListTmp = set(assayInforDic[item]['peptideSeqWarnings'])
				else:
					peptideSeqWarningsListTmp = set([])
				if subitem not in peptideOutputDic[item].keys():
					if subitem not in peptideSeqErrorsListTmp and subitem not in peptideSeqWarningsListTmp:
						petideTermWithoutIssues.append(subitem)
						try:
							assayInforDic[item]['peptideSeqWithoutIssues'].append(subitem)
						except KeyError:
							assayInforDic[item].update({'peptideSeqWithoutIssues': [subitem]})
						
			# if the issueReason of peptideTerm in assayInforDic[item]['peptideSeqWarnings'] is only "All of the concentrations of the internal standard peptide are zero. The internal standard is assumed to be the endogenous peptide.",
			# move that peptideTerm into assayInforDic[item]['peptideSeqWithoutIssues']
			if 'peptideSeqErrors' not in assayInforDic[item].keys():
				assayInforDic[item]['peptideSeqErrors'] = []
			if 'peptideSeqWarnings' not in assayInforDic[item].keys():
				assayInforDic[item]['peptideSeqWarnings'] = []
			if 'peptideSeqWithoutIssues' not in assayInforDic[item].keys():
				assayInforDic[item]['peptideSeqWithoutIssues'] = []
			
			if len(assayInforDic[item]['peptideSeqWarnings']) > 0:
				keptPeptideSeqWarningsList = []
				for peptideTermTmp in assayInforDic[item]['peptideSeqWarnings']:
					if len(peptideOutputDic[item][peptideTermTmp]['issueReason']) == 1 and peptideOutputDic[item][peptideTermTmp]['issueReason'][0] == "All of the concentrations of the internal standard peptide are zero. The internal standard is assumed to be the endogenous peptide.":
						try:
							assayInforDic[item]['peptideSeqWithoutIssues'].append(peptideTermTmp)
						except KeyError:
							assayInforDic[item].update({'peptideSeqWithoutIssues': [peptideTermTmp]})
					else:
						keptPeptideSeqWarningsList.append(peptideTermTmp)
				assayInforDic[item]['peptideSeqWarnings'] = keptPeptideSeqWarningsList

			#print assayInforDic
			#print petideTermWithoutIssues
			# Update 
			for peptideTermTmp in petideTermWithoutIssues:
				if isinstance(peptideTermTmp,tuple):
					peptide = peptideTermTmp[0]
					precursorCharge = peptideTermTmp[1]
				else:
					peptide = peptideTermTmp
					precursorCharge = ''
				figure_src_list, base64String_list = locateFigure(peptide, precursorCharge, file_list, experiment_type, dir_tmp)
				table_src_list, df_list = locateTable(peptide, precursorCharge, file_list, experiment_type, dir_tmp)
				if peptideTermTmp not in peptideOutputDic[item].keys():
					peptideOutputDic[item].update({peptideTermTmp:{}})
				if 'imageSrc' not in peptideOutputDic[item][peptideTermTmp].keys():
					peptideOutputDic[item][peptideTermTmp].update({'imageSrc':figure_src_list})
					peptideOutputDic[item][peptideTermTmp].update({'imagePNGbase64':base64String_list})
					peptideOutputDic[item][peptideTermTmp].update({'table':df_list})
		# Count the peptides.
		assayInforDic[item]['peptideWithErrors'] = countPeptide(assayInforDic[item]['peptideSeqErrors'])
		assayInforDic[item]['peptideWithWarnings'] = countPeptide(assayInforDic[item]['peptideSeqWarnings'])
		assayInforDic[item]['peptideWithoutIssues'] = countPeptide(assayInforDic[item]['peptideSeqWithoutIssues'])
		# Deduplicate the issueReason in peptideOutputDic
		for peptideTermTmp2 in peptideOutputDic[item].keys():
			if 'issueReason' in peptideOutputDic[item][peptideTermTmp2].keys():
				issueReasonList = deduplicateList(peptideOutputDic[item][peptideTermTmp2]['issueReason'])
				peptideOutputDic[item][peptideTermTmp2]['issueReason'] = issueReasonList

	assayFileWithoutIssueList = []
	for item in assayFileList:
		if item not in qc_report_infor_df.SkyDocumentName.unique():
			assayFileWithoutIssueList.append(item)
	for item in assayFileWithoutIssueList:
		peptideOutputDic.update({item:{}})
		assayInforDic[item]['isQuality'] = 'Correct'
		assayInforDic[item]['peptideSeqErrors'] = []
		assayInforDic[item]['peptideSeqWarnings'] = []
		assayInforDic[item]['peptideSeqWithoutIssues'] = peptideTrackDic[item]
		#print peptideTrackDic[item]
		#print peptideTrackDic[item][0]
		for peptideTermTmp in assayInforDic[item]['peptideSeqWithoutIssues']:
			if isinstance(peptideTermTmp,tuple):
				peptide = peptideTermTmp[0]
				precursorCharge = peptideTermTmp[1]
			else:
				peptide = peptideTermTmp
				precursorCharge = ''
			figure_src_list, base64String_list = locateFigure(peptide, precursorCharge, file_list, experiment_type, dir_tmp)
			table_src_list, df_list = locateTable(peptide, precursorCharge, file_list, experiment_type, dir_tmp)
			if peptideTermTmp not in peptideOutputDic[item].keys():
				peptideOutputDic[item].update({peptideTermTmp:{}})
			if 'imageSrc' not in peptideOutputDic[item][peptideTermTmp].keys():
				peptideOutputDic[item][peptideTermTmp].update({'imageSrc':figure_src_list})
				peptideOutputDic[item][peptideTermTmp].update({'imagePNGbase64':base64String_list})
				peptideOutputDic[item][peptideTermTmp].update({'table':df_list})
		# Count the peptides.
		if 'peptideSeqErrors' not in assayInforDic[item].keys():
			assayInforDic[item]['peptideSeqErrors'] = []
		if 'peptideSeqWarnings' not in assayInforDic[item].keys():
			assayInforDic[item]['peptideSeqWarnings'] = []
		if 'peptideSeqWithoutIssues' not in assayInforDic[item].keys():
			assayInforDic[item]['peptideSeqWithoutIssues'] = []
		assayInforDic[item]['peptideWithErrors'] = countPeptide(assayInforDic[item]['peptideSeqErrors'])
		assayInforDic[item]['peptideWithWarnings'] = countPeptide(assayInforDic[item]['peptideSeqWarnings'])
		assayInforDic[item]['peptideWithoutIssues'] = countPeptide(assayInforDic[item]['peptideSeqWithoutIssues'])