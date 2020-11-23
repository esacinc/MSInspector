# CHANGE LOG:
# 

__author__ = 'Yin Lu'
__copyright__ = 'Copyright 2018, ESAC, Inc'
__credits__ = ['Yin Lu']
__license__ = 'GPL'
__version__ = '1.1'
__maintainer__ = 'Yin Lu'
__email__ = 'yin.lu@esacinc.com'

# System imports
import os
import sys
import argparse
import time
import shutil
import warnings
from distutils.spawn import find_executable
from zipfile import ZipFile
import subprocess
try:
	import pandas as pd
except ImportError:
    print 'The package named pandas is mising. Please refer to the tutorial to install it.'
    sys.exit(1)
try:
	import jinja2
except ImportError:
    print 'The package named Jinja2 is mising. Please refer to the tutorial to install it.'
    sys.exit(1)

from MSInspector.utils.parseArgument import *
from MSInspector.utils.qcAnalysis import *
from MSInspector.utils.fileProcess import *
from MSInspector.utils.fileParse import *

def main():
	warnings.filterwarnings("ignore")
	# Parse the input arguments
	args = parseArgument()
	# Preprocess the inout arguments
	experiment_type = args.experiment_type
	experiment_type = experiment_type.lower()
	skyrTempsDir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'skyrTemps')
	rScriptsDir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'rScripts')
	htmlTempsDir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'htmlTemps')

	if experiment_type not in ['exp1', 'exp2', 'exp3', 'exp4', 'exp5']:
		print >> sys.stderr, "Invalid experiment type. Please check it."
		sys.exit(1)
	skyrTemp = os.path.join(skyrTempsDir, "MSInspector_report.skyr")
	reportName = "MSInspector_report"

	# The required columns for exp1, exp2, exp3, exp4 and exp5 are defined in a dictionary.
	common_col_list = ['ProteinName', 'PeptideModifiedSequence', 'IsotopeLabelType', 'PrecursorCharge', 'ProductCharge', 'FragmentIon', 'Area', 'ReplicateName']
	required_col_dic = {'exp1':{'old':common_col_list+['Replicate', 'Background', 'SampleGroup', 'ISSpike', 'PeptideConcentrationIS', 'Concentration', 'PeptideConcentration', 'MultiplicationFactor', 'donotuse'],
							'new':common_col_list+['ReplicateNumber', 'Background', 'SampleGroup', 'InternalStandardConcentration', 'AnalyteConcentration', 'ConcentrationMultiplier', 'donotuse']},
					'exp2':{'old':common_col_list+['Replicate', 'Concentration', 'SampleGroup'],
						    'new':common_col_list+['ReplicateNumber', 'Day', 'Exp2SampleGroup']},
					'exp3':common_col_list+['ReplicateNumber', 'AnalyteConcentration', 'Exp3SampleGroup'],
					'exp4':common_col_list+['ReplicateNumber', 'Exp4SampleGroup', 'Temperature', 'Time', 'TimeUnits', 'FreezeThawCycles'],
					'exp5':common_col_list+['ReplicateNumber', 'Day', 'SampleGroup']
					}
	# A dictionary about the data types for several required columns is defined. The '' is ignored automatically when judging the data type.
	# The data types of the required columns will be checked subsequently.
	# If the content of the column is in a controlled vocabulary, the values will be checked as well.
	col_dataType_dic = {'exp1':{'old':{'Replicate':('integer'), 'ISSpike':('integer', 'zero integer', 'number'), 'PeptideConcentrationIS':('integer', 'zero integer', 'number'), 'Concentration':('integer', 'zero integer', 'number'), 'PeptideConcentration':('integer', 'zero integer', 'number'), 'MultiplicationFactor':('integer', 'zero integer', 'number')},
							'new':{'ReplicateNumber':('integer'), 'InternalStandardConcentration':('integer', 'zero integer', 'number'), 'AnalyteConcentration':('integer', 'zero integer', 'number'), 'ConcentrationMultiplier':('integer', 'zero integer', 'number')}
							},
					'exp2':{'old':{'Replicate':('integer'), 'Concentration':('integer')}, 
					       'new':{'ReplicateNumber':('integer'), 'Day':('integer')}
					      	},
					'exp3':{'ReplicateNumber':('integer'), 'AnalyteConcentration':('integer', 'zero integer', 'number')},
					'exp4':{'ReplicateNumber':('integer'), 'Temperature':('number', 'negative number', 'integer', 'zero integer'), 'Time':('zero integer', 'number', 'integer'), 'FreezeThawCycles':('integer', 'zero integer')},
					'exp5':{'ReplicateNumber':('integer'), 'Day':('integer')}
					}
	# In the exported *.tsv from skyline document, there will be some columns where missing values are allowed.
	# A dictionary is defined for exp1, exp2, exp3, exp4 and exp5
	waived_col_dic = {'exp1':{'old':['Area','Background','Concentration', 'PeptideConcentration', 'ISSpike', 'PeptideConcentrationIS','MultiplicationFactor', 'donotuse'],
				              'new':['Area','Background','ConcentrationMultiplier', 'donotuse']},
					'exp2':{'old':['Area'],
						    'new':['Area']},
					'exp3':['Area'],
					'exp4':['Area', 'Time', 'TimeUnits'],
					'exp5':['Area'],
					}
	
	mypeptideType_file = args.mypeptideType_file
	if experiment_type == 'exp1':
		# check whether mypeptideType_file is a file.
		# mypeptideType_file is compulsory when experiment_type == 'exp1'.
		if os.path.isfile(mypeptideType_file):
			pass
		elif mypeptideType_file == 'Null':
			print >> sys.stderr, "Since -e is exp1, the directory of the file whose first column is the *.sky.zip and second column is peptide type must be indicated. Please check it"
			sys.exit(1)
		else:
			print >> sys.stderr, "the file with peptide types which is located at %s can't be found. Please check it."%(mypeptideType_file)
			sys.exit(1)
	# if mypeptideType_file == 'Null' and  experiment_type != 'exp1', it means that this file won't be considered for exp2, exp3, exp4 and exp5
	
	# Judge the mypeptideType in mypeptideType_file.
	mypeptideType_Dic = {}
	if experiment_type == 'exp1':
		with open(mypeptideType_file, 'r') as inf:
			line = inf.readline()
			if line.strip() != 'SkyDocumentName\tpeptide_standard_purity':
				print >> sys.stderr, 'The first row of mypeptideType_file should be "SkyDocumentName\tpeptide_standard_purity". Please check it.' %(mypeptideType_file)
				sys.exit(1)
			line = inf.readline()
			while line != "":
				contentTmp = line.strip().split('\t')
				fileName = contentTmp[0]
				mypeptideType = contentTmp[1]
				if mypeptideType.lower() not in ["purified", "crude"]:
					print >> sys.stderr, 'mypeptideType of %s is not purified or crude. Please check it.' %(mypeptideType_file)
					sys.exit(1)
				mypeptideType_Dic.update({fileName:mypeptideType.lower()})
				line = inf.readline()

	dir = os.path.abspath(args.output_dir)
	# set the outout directory
	tStamp = time.strftime('%Y-%m-%d_%H-%M-%S',time.localtime(time.time()))
	outputdir = os.path.join(dir, experiment_type+'.'+tStamp)
	if os.path.isdir(outputdir):
		shutil.rmtree(outputdir, ignore_errors=True)
	skyFileTmpdir = os.path.join(outputdir, 'skyFiles.tmp')
	# outputdir stores all the output files and skyFileTmpdir stores the unzipped *.sky
	os.mkdir(outputdir)
	os.mkdir(skyFileTmpdir)
	
	#plot_output = args.plot_output
	plot_output = True
	if plot_output:
		plot_output_dir = os.path.join(outputdir, 'figures_tables.tmp')
		if os.path.isdir(plot_output_dir):
			shutil.rmtree(plot_output_dir, ignore_errors=True)
		os.mkdir(plot_output_dir)
	else:
		plot_output_dir = None

	skylineCmdLog = os.path.join(outputdir, 'skylineCmdRun.log')
	skylineCmdLogOutf = open(skylineCmdLog, "w")
	outf1 = os.path.join(outputdir, 'QC_report.tsv')
	outf2 = os.path.join(outputdir, 'normal_data.tsv')
	outf3 = os.path.join(outputdir, 'file_namelist_IS.tsv')
	outf5 = os.path.join(outputdir, 'peptide_excluded_in_Rscript_infor.tsv')

	SkylineCmdBinary = args.SkylineCmdBinary
	# Check the binaries directory to make sure it can work
	#SkylineCmdBinary = os.path.join(SkylineCmdBinary, 'SkylineCmd.exe')
	executable1 = find_executable(SkylineCmdBinary)
	if not executable1:
		print >> sys.stderr, "SkylineCmd.exe can't be found in %s. Please check it."%(SkylineCmdBinary)
		sys.exit(1)
	
	RscriptBinary = args.RscriptBinary
	executable2 = find_executable(RscriptBinary)
	if not executable2:
		print >> sys.stderr, "Rscript.exe can't be found in %s. Please check it."%(RscriptBinary)
		sys.exit(1)
	
	# Judge the input type: multiple *.sky.zip files or a directory.
	input = args.input
	skyzip_file_dir_list = []
	if len(input) == 1 and os.path.isdir(os.path.abspath(input[0])):
		# the input is a directory and all the files with suffix of *.sky.zip will be traversed.
		dirTmp = os.path.abspath(input[0])
		for subfile in os.listdir(dirTmp):
			if subfile[-8:] == ".sky.zip":
				skyzip_file_dir_list.append(os.path.join(dirTmp, subfile))
		if len(skyzip_file_dir_list) == 0:
			print >> sys.stderr, "There is no *.sky.zip file in the directory of %s. Please check it."%(input[0])
			sys.exit(1)
	else:
		if all(os.path.basename(inputTmp)[-8:]=='.sky.zip' for inputTmp in input):
			for inputTmp in input:
				skyzip_file_dir_list.append(os.path.abspath(inputTmp))
			if len(skyzip_file_dir_list) == 0:
				print >> sys.stderr, "The input files are not valid *.sky.zip files. Please check it."
				sys.exit(1)
		else:
			invalidList = [os.path.basename(inputTmp) for inputTmp in input if os.path.basename(inputTmp)[-8:] != '.sky.zip']
			print >> sys.stderr, "The input files %s are not valid *.sky.zip files. Please check it." %('; '.join(invalidList))
			sys.exit(1)
	# Judge whether there are duplicated skyzip file in the input list
	# Pay attention to the format of name of *.skyzip file, usually it's suffilx is *_2017-02-03_18-00-08.sky.zip. After being unzipped, _2017-02-03_18-00-08 will be removed.
	skyzip_file_dir_basename_list_raw = [os.path.basename(item) for item in skyzip_file_dir_list]
	skyzip_file_dir_basename_list = []
	for item in skyzip_file_dir_list:
		subitem = os.path.basename(item)
		#if subitem[:-8].split('_') >= 2:
		if len(subitem[:-8].split('_')) >= 2:
			listTmp = [subitem3 for subitem2 in subitem[:-8].split('_')[-2:] for subitem3 in subitem2.split('-')]
			if all(subitem4.isdigit() for subitem4 in listTmp):
				skyzip_file_dir_basename_list.append('_'.join(subitem[:-8].split('_')[:-2])+'_*_*.sky.zip')
			else:
				skyzip_file_dir_basename_list.append(subitem)
		else:
			skyzip_file_dir_basename_list.append(subitem)
	skyzip_file_dir_duplicated_list = [item for item in set(skyzip_file_dir_basename_list) if skyzip_file_dir_basename_list.count(item) > 1]
	if len(skyzip_file_dir_duplicated_list) > 0:
		print >> sys.stderr, "There are duplicated *.sky.zip files in the input which are: %s. Please check it." %('; '.join(skyzip_file_dir_duplicated_list))
		sys.exit(1)
	
	# If experiment_type == 'exp1', judge whether all the *.sky.zip files have peptide_standard_purity types by comparing with the mypeptideType_file.
	if experiment_type == 'exp1':
		if len(skyzip_file_dir_basename_list_raw) != len(mypeptideType_Dic):
			print >> sys.stderr, "The number of *.sky.zip files in %s is not equal to those in the input of the parameter -i. Please check it." %(mypeptideType_file)
			sys.exit(1)
		for fileName in skyzip_file_dir_basename_list_raw:
			if fileName not in mypeptideType_Dic.keys():
				print >> sys.stderr, "The  *.sky.zip file from the input of the parameter -i which is %s can't be found in %s. Please check it." %(fileName, mypeptideType_file)
				sys.exit(1)
	# Check is done!!!

	# Unzip the sky.zip files in skyzip_file_dir_list one by one
	print "Start to access and parse *.sky.zip files..."
	time1 = time.time()
	skyTsvDirList = []
	skyFileDirList = []
	#pandasList = []
	fileNameList = []
	fileNameSkylineTmpTypeDic = {}
	fileNameEmptyStatusDic = {}
	fileNameInferredExpStatusDic = {}
	for skyzip_file_dir in skyzip_file_dir_list:
		zf = ZipFile(skyzip_file_dir, 'r')
		zf.extractall(skyFileTmpdir)
		zf.close()
		#  test the suffix of os.path.basename(skyzip_file_dir)[:-8]
		subitem = os.path.basename(skyzip_file_dir)
		if subitem[:-8].split('_') >= 2:
			listTmp = [subitem3 for subitem2 in subitem[:-8].split('_')[-2:] for subitem3 in subitem2.split('-')]
			if all(subitem4.isdigit() for subitem4 in listTmp):
				input_file_sky_tmp = '_'.join(subitem[:-8].split('_')[:-2])+'.sky'
			else:
				input_file_sky_tmp = subitem[:-8]+'.sky'
		else:
			input_file_sky_tmp = subitem[:-8]+'.sky'
		input_file_sky = os.path.join(skyFileTmpdir, input_file_sky_tmp)
		# Parse and transform the input_file_sky into report file based on the specific skyline document report template.
		output_file = os.path.join(skyFileTmpdir, input_file_sky_tmp[:-4]+'_tmp.tsv')
		# Run SkylineCmd.exe toward the sky document file and export a *.tsv report file
		#print input_file_sky
		#print output_file
		#print '"%s" --timestamp --in="%s" --report-file="%s" --report-format=TSV --report-name=%s --report-add=%s --report-conflict-resolution=overwrite --report-invariant >> %s'%(SkylineCmdBinary, input_file_sky, output_file, reportName, skyrTemp, skylineCmdLog)
		ps = subprocess.Popen('"%s" --timestamp --in="%s" --report-file="%s" --report-format=TSV --report-name=%s --report-add=%s --report-conflict-resolution=overwrite --report-invariant'%(SkylineCmdBinary, input_file_sky, output_file, reportName, skyrTemp),shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		#print '"%s" --timestamp --in="%s" --report-file="%s" --report-format=TSV --report-name=%s --report-add=%s --report-conflict-resolution=overwrite --report-invariant'%(SkylineCmdBinary, input_file_sky, output_file, reportName, skyrTemp)
		while True:
			buff = ps.stdout.readline()
			if buff == '' and ps.poll() != None:
				break
			else:
				skylineCmdLogOutf.write(buff.strip()+'\n')
		skylineCmdLogOutf.write('****************\n')
		#os.system('"%s" --timestamp --in="%s" --report-file="%s" --report-format=TSV --report-name=%s --report-add=%s --report-conflict-resolution=overwrite --report-invariant >> %s'%(SkylineCmdBinary, input_file_sky, output_file, reportName, skyrTemp, skylineCmdLog))
		#os.system('"%s" --timestamp --in=%s --report-file=%s --report-format=TSV --report-name=%s --report-add=%s --report-conflict-resolution=overwrite >> %s'%(SkylineCmdBinary, input_file_sky, output_file, reportName, skyrTemp, skylineCmdLog))
		output_file_new = os.path.join(skyFileTmpdir, input_file_sky_tmp[:-4]+'.tsv')
		# judge whether output_file is successfully generated or not
		if os.path.isfile(output_file):
			for rowNumber, line in enumerate(open(output_file, 'rU')):
			   pass
			rowNumber += 1
			if rowNumber < 2:
				print >> sys.stderr, "The *.sky file in %s is blank after been parsed by %s. Please check it."%(skyzip_file_dir, SkylineCmdBinary)
				sys.exit(1)
		else:
			print >> sys.stderr, "The file %s can't be accessed by %s. Please check the version compatibility between Skyline document and Skyline program."%(skyzip_file_dir, SkylineCmdBinary)
			sys.exit(1)
		
		# Firstly, modify the output_file. If the row number of output_file is 0, True will be stored in fileNameEmptyStatusDic.
		# Then, infer the experiment Type based on output_file_new,
		# if the inferred experiment type is inconsistent with the experiment type from input arguments, True will be stored in fileNameInferredExpStatusDic.
		# Meanwhile, judge which skyline template is used when experiment_type is exp1 or exp2. For example, new or old
		skylineTempType = headerModify(output_file, output_file_new, experiment_type, subitem, fileNameInferredExpStatusDic, fileNameEmptyStatusDic)
		skyTsvDirList.append(output_file_new)
		skyFileDirList.append(input_file_sky)
		fileNameList.append(os.path.basename(skyzip_file_dir))
		fileNameSkylineTmpTypeDic.update({os.path.basename(skyzip_file_dir):skylineTempType})
	# If some Skyline files are blank, the program will stop and error message will be printed.
	if any([fileNameEmptyStatusDic[fileName] for fileName in fileNameList]):
		fileNameBlankList = [fileName for fileName in fileNameList if fileNameEmptyStatusDic[fileName]]
		print >> sys.stderr, "The input *.sky.zip file(s), including %s, has(have) no data. Please check it(them) in Skyline.\n"%(", ".join(fileNameBlankList))
		sys.exit(1)
	# If the inferred experiment type is different from the input one for some Skyline files, the program will stop and error message will be printed.
	if any([fileNameInferredExpStatusDic[fileName] for fileName in fileNameList]):
		fileNameWrongExpList = [fileName for fileName in fileNameList if fileNameInferredExpStatusDic[fileName]]
		print >> sys.stderr, "The experiment type(s) of the input *.sky.zip file(s), including %s, is(are) wrong according to the annotated attributes. Please check it(them) in Skyline.\n"%(", ".join(fileNameWrongExpList))
		sys.exit(1)
	# Because for exp1 and exp2, there may be two skyline templates. All the input *.sky.zip files must share the same template, otherwise, an error will be thrown and the program will be stopped.
	if experiment_type in ['exp1', 'exp2']:
		if len(set([fileNameSkylineTmpTypeDic[fileName] for fileName in fileNameList])) > 1:
			print >> sys.stderr, "The input *.sky.zip files, including %s, are not annotated using the same Skyline template. Please check them in Skyline.\n"%(", ".join(fileNameList))
			sys.exit(1)

	time2 = time.time()
	skylineCmdLogOutf.close()
	if experiment_type == 'exp1':
		fileNameSelceted = fileNameSkylineTmpTypeDic.keys()[0]
		skylineTemp_type =  fileNameSkylineTmpTypeDic[fileNameSelceted]
		if skylineTemp_type == 'old':
			rScript = os.path.join(rScriptsDir, "experiment1_oldTemplate_qc.R")
		else:
			rScript = os.path.join(rScriptsDir, "experiment1_newTemplate_qc.R")
	elif experiment_type == 'exp2':
		fileNameSelceted = fileNameSkylineTmpTypeDic.keys()[0]
		skylineTemp_type =  fileNameSkylineTmpTypeDic[fileNameSelceted]
		if skylineTemp_type == 'old':
			rScript = os.path.join(rScriptsDir, "experiment2_oldTemplate_qc.R")
		else:
			rScript = os.path.join(rScriptsDir, "experiment2_newTemplate_qc.R")
	elif experiment_type == 'exp3':
		rScript = os.path.join(rScriptsDir, "experiment3_qc.R")
	elif experiment_type == 'exp4':
		rScript = os.path.join(rScriptsDir, "experiment4_qc.R")
		reportName="experiment4"
	else :
		rScript = os.path.join(rScriptsDir, "experiment5_qc.R")
		reportName="experiment5"
	
	print "It takes %.2fsec."%(time2-time1)
	# Step 1: QC each *.sky in skyTsvDirList
	# In this step, the peptides with missing values and duplicated peptides will be stored in errorDf, the rest peptides will be stored in normalDf
	print "QC of %s is running..."%(reportName)
	errorDf, normalDf, peptide_excluded_in_Rscript_Df = qcAnalysisGlobal(experiment_type, skyTsvDirList, fileNameList, required_col_dic, waived_col_dic, fileNameSkylineTmpTypeDic, col_dataType_dic)
	errorDfColNumber = errorDf.columns.size
	errorDf.to_csv(outf1, sep='\t', header=True, index=False)
	normalDf.to_csv(outf2, sep='\t', header=True, index=False)
	peptide_excluded_in_Rscript_Df.to_csv(outf5, sep='\t', header=True, index=False)
	# Use the Series in pandas to deduplicate fileName rapidly
	#uniqueSkyFileList = normalDf['SkyDocumentName'].value_counts().index.tolist()
	#fileNameList = [fileName for fileName in fileNameList if fileName in uniqueSkyFileList]
	skyFileDirListTmp = [os.path.basename(skyFileDir) for skyFileDir in skyFileDirList]
	if len(skyFileDirListTmp) != len(set(skyFileDirListTmp)):
		# It means there are duplicated *.sky files afer *.sky.zip files are unzipped, the program will be forced to exit.
		skyFileDirDup_list = []
		for item_tmp1 in set(skyFileDirListTmp):
			if skyFileDirListTmp.count(item_tmp1) > 1:
				skyFileDirDup_list.append(item_tmp1)
		print >> sys.stderr, "There are multiple *.sky files with the same name(s) in the input after being unzipped, which are(is) : %s. Please check it."%('; '.join(skyFileDirDup_list))
		sys.exit(1)
	
	with open(outf3, 'w') as outfTmp:
		outfTmp.write('SkyDocumentName\tinternal_standard\n')
		for index, item in enumerate(fileNameList):
			internal_standard = detectIS(skyFileDirList[index], fileNameList[index], experiment_type, outf1, errorDfColNumber)
			outfTmp.write(item+'\t'+internal_standard+'\n')
	#print fileNameList
	# Step 2: QC each fileName in fileNameList from normalDf for the specific experiment type by running the corresponding R code.
	qcAnalysisRcode(experiment_type, outf1, outf2, outf3, mypeptideType_file, RscriptBinary, rScript, plot_output, plot_output_dir)
	# Step 3: Create the QC_report.html file
	assayFileList = []
	assayFileAnchorDic = {}
	assayInforDic = {}
	peptideTrackDic = {}
	peptideTrackAnchorDic = {}
	peptideSeqChargeIsotopeDic = {}
	peptideOutputDic = {}
	peptide_infor_file = os.path.join(plot_output_dir, 'peptide_infor.tsv')
	is_inferred_file = os.path.join(plot_output_dir, 'internal_standard_inferred_infor.tsv')
	# Step 3.1: Parse peptide_infor_file and add the data into assayFileList and assayInforDic
	peptide_infor_parse(outf3, peptide_infor_file, outf5, assayFileList, assayInforDic, peptideTrackDic, peptideSeqChargeIsotopeDic, experiment_type)
	# Step 3.2: Parse QC_report.tsv and add the data into assayInforDic
	# Normally if the R script comes across errors, the function of qc_report_infor_parse will come across IOError
	"""
	try:
		qc_report_infor_parse(outf1, is_inferred_file, assayInforDic, peptideTrackDic, peptideOutputDic, assayFileList, experiment_type)
	except Exception, e:
		print e
		print '\n'
		print "Errors happen when executing the R script. You can: 1) Check the compatibility of the installed R packages in R console by referring to the tutorial; 2) If there is no problem in 1), it means that some undefined errors in the Skyline document exist and please contact Yin Lu through yin.lu@eascinc.com to report the issue."
		sys.exit(1)
	"""
	qc_report_infor_parse(outf1, is_inferred_file, assayInforDic, peptideTrackDic, peptideOutputDic, assayFileList, experiment_type)
	# For each key in peptideSeqDisplay, assign an unique anchor id
	id = 0
	for item in assayFileList:
		peptideTrackAnchorDic.update({item:{}})
		for subitem in peptideTrackDic[item]:
			id = id + 1
			peptideTrackAnchorDic[item][subitem] = "id"+str(id)
	keptAssayFileErrorTable = []
	keptAssayErrorTable = {}
	keptAssayFileWarningTable = []
	keptAssayFileWithoutIssueTable = []
	
	id1 = 0
	for item in assayFileList:
		if assayInforDic[item]['isQuality'] in ['Correct',
											 'The internal standard in the skyline file is set to be none. Errors happen for all the peptides.',
											 'The internal standard in the skyline file is set to be none. Please set it to be heavy. Errors happen for all the peptides.',
											 'The internal standard in the skyline file is set to be heavy, while the inferred internal standard is light. Errors happen for all the peptides.',
											 'The internal standard in the skyline file is set to be light, while the inferred internal standard is heavy. Errors happen for all the peptides.',
											 "Internal standard type can't be inferred. All the peptides have missing values or incorrect data types in some essential attributes."]:
			if len(assayInforDic[item]['peptideSeqErrors']) > 0:
				for subitem in assayInforDic[item]['peptideSeqErrors']:
					if subitem in peptideOutputDic[item].keys():
						if item not in keptAssayFileErrorTable:
							keptAssayFileErrorTable.append(item)
						try:
							keptAssayErrorTable[item].append(subitem)
						except:
							keptAssayErrorTable.update({item:[subitem]})
			if len(assayInforDic[item]['peptideSeqWarnings']) > 0:
				keptAssayFileWarningTable.append(item)
				id1 = id1 + 1
				if item not in assayFileAnchorDic.keys():
					assayFileAnchorDic.update({item:{}})
				assayFileAnchorDic[item]['peptideSeqWarnings'] = "skyid"+str(id1)
			if len(assayInforDic[item]['peptideSeqWithoutIssues']) > 0:
				keptAssayFileWithoutIssueTable.append(item)
				id1 = id1 + 1
				if item not in assayFileAnchorDic.keys():
					assayFileAnchorDic.update({item:{}})
				assayFileAnchorDic[item]['peptideSeqWithoutIssues'] = "skyid"+str(id1)
	
	# Step 3.3: Render the html template using the data of assayFileList.
	jinja_env = jinja2.Environment(loader = jinja2.FileSystemLoader(searchpath=htmlTempsDir))
	template_file = "report_template.html"
	template = jinja_env.get_template(template_file)
	#print keptAssayFileErrorTable
	#print "$$$$$$$$$$"
	#print keptAssayErrorTable
	#print "$$$$$$$$$$"
	#print assayInforDic
	#print "$$$$$$$$$$"
	#print peptideOutputDic
	html_rendered = template.render(assayFileList=assayFileList, assayInforDic=assayInforDic, peptideOutputDic=peptideOutputDic, peptideSeqChargeIsotopeDic=peptideSeqChargeIsotopeDic, keptAssayFileErrorTable=keptAssayFileErrorTable, keptAssayErrorTable=keptAssayErrorTable, keptAssayFileWarningTable= keptAssayFileWarningTable, keptAssayFileWithoutIssueTable=keptAssayFileWithoutIssueTable, peptideTrackAnchorDic=peptideTrackAnchorDic, experiment_type=experiment_type, assayFileAnchorDic= assayFileAnchorDic)
	outf4 = os.path.join(outputdir, 'QC_report.html') 
	with open(outf4, 'w') as handle:
		handle.write(html_rendered)
	
	# Step 3.4: Output a csv file which contains the tables for all the peptides. The tables should be like what appear on the assay portal detail pages for each experiment.
	outf6 = os.path.join(outputdir, 'statistical_tables_on_Assay_Portal.csv')
	output_statistical_tables(peptide_infor_file, plot_output_dir, outf6, experiment_type)
	time3 = time.time()
	print "It takes %.2fsec."%(time3-time2)
	

if __name__ == '__main__':
	main()
