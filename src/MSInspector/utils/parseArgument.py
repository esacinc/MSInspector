import os
import sys
import argparse

def parseArgument():
	dir_tmp = os.getcwd()
	version = '1.1'
	usageTmp = '\r{}\n\
##                                                                           ##\n\
##      MSInspector (Quality control of the five experiments in            ##\n\
##                       CPTAC assay portal)                                 ##\n\
##                                                                           ##\n\
##      last change: 01/15/2024                                              ##\n\
##                                                                           ##\n\
##                                                                           ##\n\
###############################################################################\n\
\n'.format('###############################################################################'.ljust(len('usage:')))
	#usage = usageTmp+r'Usage: MSInspector [-h] [<args>]\n\nExample 1:\nMSInspector -ps "C:\Program Files\SkylineDailyRunner\SkylineDailyRunner.exe"  -pr "C:\Program Files\R\R-3.5.2\\bin\Rscript.exe" -i "D:\Skyline_analysis\qcAssayPortal\data\UVicPC_Borchers-MousePlasma_Agilent6490_directMRM-Exp1\20160309_MouseV2B1.5_refined_2018-07-03_14-59-18.sky.zip" -e exp1 -t "D:\Skyline_analysis\qcAssayPortal\data\UVicPC_Borchers-MousePlasma_Agilent6490_directMRM-Exp1\meta_data.tsv"\nExample 2:\nMSInspector -ps "C:\Program Files\Skyline\SkylineCmd.exe" -pr "C:\Program Files\R\R-3.5.2\\bin\Rscript.exe" -i "D:\Skyline_analysis\qcAssayPortal\data\UVicPC_Borchers-MousePlasma_Agilent6490_directMRM-Exp2" -e exp2\n'
	usage = usageTmp+'Usage: MSInspector [-h] [<args>]\n\nExample 1:\nMSInspector -ps '+r'"C:\Program Files\Skyline\SkylineCmd.exe""'+' -pr '+r'"C:\Program Files\R\R-4.3.2\bin\Rscript.exe"'+' -i '+r'"D:\Skyline_analysis\qcAssayPortal\data\UVicPC_Borchers-MousePlasma_Agilent6490_directMRM-Exp1\20160309_MouseV2B1.5_refined_2018-07-03_14-59-18.sky.zip"'+' -e exp1 -t '+r'"D:\Skyline_analysis\qcAssayPortal\data\UVicPC_Borchers-MousePlasma_Agilent6490_directMRM-Exp1\meta_data.tsv"'+'\nExample 2:\nMSInspector -ps '+r'"C:\Program Files\Skyline\SkylineCmd.exe"'+' -pr '+r'"C:\Program Files\R\R-4.3.2\bin\Rscript.exe"'+' -i '+r'"D:\Skyline_analysis\qcAssayPortal\data\UVicPC_Borchers-MousePlasma_Agilent6490_directMRM-Exp2"'+' -e exp2\n'
	#parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description=usage, usage=argparse.SUPPRESS)
	parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description=usage, usage=argparse.SUPPRESS)
	parser.add_argument('-ps', required=True, dest='SkylineCmdBinary', metavar='<dir required>', help="the path to SkylineCmd.exe, SkylineRunner.exe, or SkylineDailyRunner.exe in the Windows OS")
	parser.add_argument('-pr', required=True, dest='RscriptBinary', metavar='<dir required>', help="the path to Rscript.exe in the Windows OS")
	parser.add_argument('-i', nargs='+', required=True, dest='input', metavar='input <required>', help='two options: 1. A directory where all the *.sky.zip files are located 2.*.sky.zip *.sky.zip ... (at least one input *.sky.zip)')
	parser.add_argument('-e', required=True, dest='experiment_type', metavar='string <required>', help='the experiment type. Choose one from Options: exp1 , exp2 , exp3 , exp4, exp5')
	parser.add_argument('-t', default='Null', dest='mypeptideType_file', metavar='peptide type file <required>', help='the directory of the file whose first column is the *.sky.zip and second column is peptide type ("purified" or "crude"). When -e is exp1, it must be assigned. Otherwise, please leave it blank.')
	parser.add_argument('-o', default=dir_tmp, dest='output_dir', metavar='<dir>', help='the directory of the outputs (default: current directory)')
	#parser.add_argument('-p', dest='plot_output', action = 'store_true', help='switch to plot figures and generate tables of the specific experiment type (default: off)')
	parser.add_argument('--version', action='version', version='%s'%(version))
	return parser.parse_args()