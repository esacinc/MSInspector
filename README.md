MSInspector
================================================

A Python program to evaluate the data quality of the five experiments in CPTAC assay portal, the format of input files is *.sky.zip.<br />
The experiments are:<br />
Experiment 1, response curve (or ResponseCurve)<br />
Experiment 2, evaluation of repeatability (or ValidationSamples)<br />
Experiment 3, selectivity (or Selectivity)<br />
Experiment 4, stability (Stability)<br />
Experiment 5, reproducible detection of endogenous analyte (Endogenous)<br />

The work flow of MSInspector is shown [here](https://github.com/esacinc/MSInspector/tree/master/src/MSInspector/doc/workflow.pdf).

A report file (QC_report.html) will be generated which will capture the details of errors and warnings when 
running the R codes for the individual experiment.
The details of the errors and warnings are listed [here](https://github.com/esacinc/MSInspector/tree/master/src/MSInspector/doc/issue_categories.pdf).

Documentation
-------------

* [Installation](#installation)
  * [Install MSInspector](#install-MSInspector)
* [How to use it](#how-to-use-it)
  * [Command line](#command-line)
* [Changelog](#changelog)
* [Citation](#citation)

Installation
------------

### Install MSInspector

MSInspector is implemented as a Python program running on a Windows platform that requires pre-installation of Skyline https://skyline.ms/project/home/software/Skyline/begin.view and R https://cran.r-project.org/bin/windows/base/  <br />

It has the following dependencies:
* For the Python (v3.11.\*) programming language: requires Python-related libraries, including pandas(>= v2.0.0) and Jinja2(v3.1.2). <br /><br />
* For the R (v4.3.\* is recommended) programming language: requires R-related libraries, including Cairo, evaluate, stringr, plyr, MASS, ggplot2, reshape2 and dplyr. In the R console, in order to install required libraries, please run `install.packages(c("Cairo", "evaluate", "stringr", "plyr", "MASS", "ggplot2", "reshape2", "dplyr"))`. <br /><br />
* For Skyline (v23.1 is recommended): requires pre-installed Skyline command-line interface https://skyline.ms/_webdav/home/software/Skyline/@files/docs/Skyline%20Command-Line%20Interface-3_7.pdf. Or install Skyline administrator downloaded from https://skyline.ms/wiki/home/software/Skyline/page.view?name=install-administator-64. The command-line toolSkylineCmd.exe can be found in the installation directory. <br />

MSInspector can be installed from the source code by pip:<br />
1) Download MSInspector source code from URL and unzip the zipped file folder.<br />
2) If the package of wheel is not installed, run `pip install wheel` to install it.<br />
3) Change the directory to MSInspector's directory and run `python setup.py bdist_wheel` to build a wheel file for the subsequent installation via pip.<br />
4) Run `pip install .\dist\MSInspector-1.2-py2-none-any.whl` to install MSInspector.<br />


How to use it
-------------

### Command line

    Usage: MSInspector [-h] [<args>]

    Example 1:
    MSInspector -ps "C:\Program Files\Skyline\SkylineCmd.exe"" -pr "C:\Program Files\R\R-4.3.2\bin\Rscript.exe" -i "D:\Skyline_analysis\qcAssayPortal\data\UVicPC_Borchers-MousePlasma_Agilent6490_directMRM-Exp1\20160309_MouseV2B1.5_refined_2018-07-03_14-59-18.sky.zip" -e exp1 -t "D:\Skyline_analysis\qcAssayPortal\data\UVicPC_Borchers-MousePlasma_Agilent6490_directMRM-Exp1\meta_data.tsv"
    
    Example 2:
    MSInspector -ps "C:\Program Files\Skyline\SkylineCmd.exe" -pr "C:\Program Files\R\R-4.3.2\bin\Rscript.exe" -i "D:\Skyline_analysis\qcAssayPortal\data\UVicPC_Borchers-MousePlasma_Agilent6490_directMRM-Exp2" -e exp2

    optional arguments:
      -h, --help            show this help message and exit
      -ps <dir required>    the path to SkylineCmd.exe, SkylineRunner.exe, or SkylineDailyRunner.exe in the Windows OS
      -pr <dir required>    the path to Rscript.exe in the Windows OS
      -i input <required> [input <required> ...]
                            two options: 1. A directory where all the *.sky.zip files are located 2.*.sky.zip *.sky.zip ... (at least one input *.sky.zip)
      -e string <required>  the experiment type. Choose one from Options: exp1 , exp2
      -t peptide type file <required> the directory of the file whose first column is the *.sky.zip and second column is peptide type ("purified" or "crude"). When -e is exp1, it must be assigned. Otherwise, please leave it blank.
      -o <dir>              the directory of the outputs (default: current directory)
      --version             show program's version number and exit

Changelog
---------
2019-05-31 Add a function to infer which template (new or old) is being applied to annotate the data for experiment 1 and experiment 2.<br />
2019-06-13 Add a function to QA the experiment type from the input parameter.<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Add a function to QA the data type of several required attributes.<br />
2019-06-26 Fix the bug caused by the space in the directory of the input parameters.<br />
2019-07-18 Add a function to capture the error caused by wrong annotation of IS Spike for experiment 1.<br />
2019-08-22 Fix the issue of the variable width of the horizontal lines in Response Curve for the experiment 1.<br />
2019-12-03 Add errors and warnings detecting function to experiment 3.<br />
2020-01-17 Add errors and warnings detecting function to experiment 4.<br />
2020-03-19 Add errors and warnings detecting function to experiment 5.<br />
2024-01-15 Update Python from Python 2.7 to Python 3.11<br />

Citation
--------
