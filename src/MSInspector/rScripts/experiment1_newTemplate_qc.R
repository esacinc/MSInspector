# This sample code is used to qc skyline document data from experiment 1 of Assay Portal.
# It's modified based on:
# 1) esac-panorama-master\AssayPortal\resources\reports\schemas\targetedms\ResponseCurveQuery\ResponseCurve.r
# 2) esac-panorama-master\AssayPortal\resources\reports\schemas\targetedms\ResponseCurveAnalysis\ResponseCurveAnalysis.r

suppressWarnings(suppressMessages(library(Cairo)))
#suppressWarnings(suppressMessages(library(Rlabkey)))
suppressWarnings(suppressMessages(library(evaluate)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(require(plyr)))
suppressWarnings(suppressMessages(require(ggplot2)))
suppressWarnings(suppressMessages(require(MASS)))
suppressWarnings(suppressMessages(require(reshape2)))
suppressWarnings(suppressMessages(require(dplyr)))

fitline <- function(xall, yall){
    TotalPoint <- length(xall)
    if (TotalPoint -2 >=1){
        usePoint <- TotalPoint-2
        myRSquare <- 1;
        while (myRSquare > 0.98 && usePoint >=1){
            x <- xall[usePoint : TotalPoint]
            y <- yall[usePoint : TotalPoint]
            w = 1/y;
            fit <- lm(y~x, weights = w)
            myRSquare <- summary(fit)$r.squared
            usePoint <- usePoint - 1
        }
        if (myRSquare <= 0.98){
            if (usePoint == (TotalPoint -3)) {
                usePoint <- usePoint + 1
            }else{
                usePoint <- usePoint + 2
            }
        }else{
            usePoint <- usePoint + 1
        }
        bestUsePoint <- usePoint
        bestEndPoint <- TotalPoint
        x <- xall[bestUsePoint : bestEndPoint]
        y <- yall[bestUsePoint : bestEndPoint]
        w = 1/y;
        fit <- lm(y~x, weights = w)
        bestRSquare <- summary(fit)$r.squared
        longest <- bestEndPoint - bestUsePoint + 1
    }
    #if (myRSquare < 0.98 && usePoint == (TotalPoint-3)){
    if ( TotalPoint-3 >= 1){
        usePoint <- TotalPoint-3
        myRSquare <- 1
        while (myRSquare > 0.98 && usePoint >=1 ){
            x <- xall[usePoint : (TotalPoint-1)]
            y <- yall[usePoint : (TotalPoint-1)]
            w = 1/y
            fit <- lm(y~x, weights = w)
            myRSquare <- summary(fit)$r.squared
            usePoint <- usePoint - 1
        }
        if (myRSquare <= 0.98){
            if (usePoint == (TotalPoint -4)) {
                usePoint <- usePoint + 1
            }else{
                usePoint <- usePoint + 2
            }
        }else{
            usePoint <- usePoint + 1
        }
        x <- xall[usePoint : (TotalPoint-1)]
        y <- yall[usePoint : (TotalPoint-1)]
        w = 1/y
        fit <- lm(y~x, weights = w)
        myRSquare <- summary(fit)$r.squared
        replace <- 0
        if (bestRSquare < 0.98 && myRSquare > 0.98){
            replace <- 1
        }else if ( ((bestRSquare > 0.98 && myRSquare > 0.98) || (bestRSquare < 0.98 && myRSquare < 0.98))&& longest < (TotalPoint-usePoint)) {
            replace <- 1
        }
        if (replace == 1) {
            bestUsePoint <- usePoint;
            bestEndPoint <- TotalPoint-1
            bestRSquare <- myRSquare;
            longest <- bestEndPoint - bestUsePoint + 1
        }
    }
    #if (myRSquare < 0.98 && usePoint == (TotalPoint-4)){
    if (TotalPoint-4 >= 1){
        usePoint <- TotalPoint-4
        myRSquare <- 1
        while (myRSquare > 0.98 && usePoint >=1){
            x <- xall[usePoint : (TotalPoint-2)]
            y <- yall[usePoint : (TotalPoint-2)]
            w = 1/y
            fit <- lm(y~x, weights = w)
            myRSquare <- summary(fit)$r.squared
            usePoint <- usePoint - 1
        }
        if (myRSquare <= 0.98){
            if (usePoint == (TotalPoint -5)) {
                usePoint <- usePoint + 1
            }else{
                usePoint <- usePoint + 2
            }
        }else{
            usePoint <- usePoint + 1
        }
        x <- xall[usePoint : (TotalPoint-2)]
        y <- yall[usePoint : (TotalPoint-2)]
        w = 1/y
        fit <- lm(y~x, weights = w)
        myRSquare <- summary(fit)$r.squared
        replace <- 0
        if (bestRSquare < 0.98 && myRSquare > 0.98){
            replace <- 1
        }else if (((bestRSquare > 0.98 && myRSquare > 0.98) || (bestRSquare < 0.98 && myRSquare < 0.98)) && longest < (TotalPoint-usePoint-1)) {
            replace <- 1
        }
        if (replace == 1) {
            bestUsePoint <- usePoint
            bestEndPoint <- TotalPoint-2
            bestRSquare <- myRSquare
            longest <- bestEndPoint - bestUsePoint + 1
        }
    }
    if ( TotalPoint <3 ) {
        list()
    }else{
        x <- xall[bestUsePoint : bestEndPoint]
        y <- yall[bestUsePoint : bestEndPoint]
        w <- 1/y
        fit <- lm(y~x, weights = w)
        mCoef <- coef(summary(fit))
        mResidual <- c(TotalPoint)
        for (i in 1: TotalPoint){
            mResidual[i] <- (yall[i]-(xall[i]*mCoef[2,1]+mCoef[1,1]))/(yall[i]+(xall[i]*mCoef[2,1]+mCoef[1,1]))
        }
        usedPoints = bestEndPoint - bestUsePoint + 1
        Rsquare <-summary(fit)$r.squared
        list(mCoef=mCoef, mResidual=mResidual, Rsquare=Rsquare, usedPoints = usedPoints)
    }
}

detect_outliers <- function(x, na.rm = TRUE, ...) {
    # This function is abandoned.
    qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
    H <- 1.5 * IQR(x, na.rm = na.rm)
    y <- x
    y[x < (qnt[1] - H)] <- NA
    y[x > (qnt[2] + H)] <- NA
    y
}

identify_uniProtKB_entryID  <- function(x) {
    # This function is to extract uniProtKB_entryID from the protein name.
    if (grepl('\\|', x)) {
        tmp <- strsplit(x, split = '\\|')
        uniProtKB_entryID_tmp <- tmp[[1]][2]
        # Judge whether uniProtKB_entryID is a legal uniProtKB entry ID based on its pattern using regular expression. Please refer to https://www.uniprot.org/help/accession_numbers
        if (str_detect(uniProtKB_entryID_tmp, "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}")) {
            uniProtKB_entryID <- uniProtKB_entryID_tmp
        } else {
            uniProtKB_entryID <- x
        }
    } else {
        uniProtKB_entryID <- x
    }
    uniProtKB_entryID
}

args <- commandArgs(trailingOnly = TRUE)
dataset_path <- args[1]
fileList_path <- args[2]
plot_output <- args[3]
plot_output_dir <- args[4]
mypeptideType_file_path <- args[5]

#dataset_path <- "normal_data.tsv"
#fileList_path <- "file_namelist_IS.tsv"
#plot_output <- "True"
#plot_output_dir <- "D:\\Skyline_analysis\\qcAssayPortal\\qcAssayPortal\\src\\qcAssayPortal\\rScripts\\test\\exp1_new\\tmp"
#mypeptideType_file_path <- "mypeptideType_file.tsv"

if (plot_output == 'True') {
    plot_output <- TRUE
} else {
    plot_output <- FALSE
}
RCAnalysisSwitch <- TRUE
qcSlopeSwitch <- TRUE

cv_threshold <- 0.5
pValue_threshold <- 0.05
rSquare_threshold <- 0.5

# Load data from local table
labkey.data.total <- read.table(file=dataset_path, header=TRUE, sep='\t')
fileDf <- read.table(file=fileList_path, header=TRUE, sep='\t')
mypeptideType_file_Df <- read.table(file=mypeptideType_file_path, header=TRUE, sep='\t')

colNumber <- ncol(labkey.data.total)
# Add the internal standard type and peptide type into labkey.data.total based on the column of SkyDocumentName
# A new column named internal_standard is added
labkey.data.total <- merge(labkey.data.total, fileDf,  by.x="SkyDocumentName", by.y="SkyDocumentName", all.x=TRUE)
# A new column named peptide_standard_purity is added
labkey.data.total <- merge(labkey.data.total, mypeptideType_file_Df,  by.x="SkyDocumentName", by.y="SkyDocumentName", all.x=TRUE)

thenames <- tolower(names(labkey.data.total))
names(labkey.data.total) <- thenames

# Rename some column names of labkey.data.total
colnames(labkey.data.total)[colnames(labkey.data.total)=="skydocumentname"] <- "SkyDocumentName"
colnames(labkey.data.total)[colnames(labkey.data.total)=="proteinname"] <- "protein"
colnames(labkey.data.total)[colnames(labkey.data.total)=="isotopelabeltype"] <- "isotopelabel"
#colnames(labkey.data.total)[colnames(labkey.data.total)=="replicatename"] <- "replicate"
colnames(labkey.data.total)[colnames(labkey.data.total)=="analyteconcentration"] <- "concentration"
colnames(labkey.data.total)[colnames(labkey.data.total)=="internalstandardconcentration"] <- "isspike"
colnames(labkey.data.total)[colnames(labkey.data.total)=="concentrationmultiplier"] <- "multiplicationfactor"
colnames(labkey.data.total)[colnames(labkey.data.total)=="replicatenumber"] <- "replicate"

# Reassign the data type for each column in labkey.data.total
labkey.data.total[,'SkyDocumentName'] <- as.character(labkey.data.total[,'SkyDocumentName'])
labkey.data.total[,'protein'] <- as.character(labkey.data.total[,'protein'])
labkey.data.total[,'peptidemodifiedsequence'] <- as.character(labkey.data.total[,'peptidemodifiedsequence'])
labkey.data.total[,'isotopelabel'] <- as.character(labkey.data.total[,'isotopelabel'])
labkey.data.total[,'precursorcharge'] <- as.character(labkey.data.total[,'precursorcharge'])
labkey.data.total[,'productcharge'] <- as.character(labkey.data.total[,'productcharge'])
labkey.data.total[,'fragmention'] <- as.character(labkey.data.total[,'fragmention'])
#labkey.data.total[,'replicate'] <- as.character(labkey.data.total[,'replicate'])
labkey.data.total[,'replicate'] <- as.numeric(as.character(labkey.data.total[,'replicate']))
labkey.data.total[,'concentration'] <-  as.numeric(as.character(labkey.data.total[,'concentration']))
labkey.data.total[,'samplegroup'] <-  as.character(labkey.data.total[,'samplegroup'])
labkey.data.total[,'area'] <-  as.numeric(as.character(labkey.data.total[,'area']))
labkey.data.total[,'background'] <-  as.numeric(as.character(labkey.data.total[,'background']))
labkey.data.total[,'donotuse'] <-  as.character(labkey.data.total[,'donotuse'])
labkey.data.total[,'isspike'] <-  as.numeric(as.character(labkey.data.total[,'isspike']))
#labkey.data.total[,'peptideconcentrationis'] <-  as.numeric(as.character(labkey.data.total[,'peptideconcentrationis']))
#labkey.data.total[,'peptideconcentration'] <-  as.numeric(as.character(labkey.data.total[,'peptideconcentration']))
labkey.data.total[,'multiplicationfactor'] <-  as.numeric(as.character(labkey.data.total[,'multiplicationfactor']))

# Add a new temporary column fragment_ion_complete
if (nrow(labkey.data.total) ==0) {
    labkey.data.total$fragment_ion_complete <- integer(0)
} else {
    labkey.data.total$fragment_ion_complete <- paste(labkey.data.total[ ,'fragmention'], " (", labkey.data.total[ ,'productcharge'], "+)", sep='' )
}

# Write peptide information into output file.
log_filename <- paste(plot_output_dir, "\\peptide_infor.tsv", sep='' )
logdf <- data.frame(peptide=as.character(), precursorCharge=as.character(), isotopeLabelType=as.character(), transition=as.character(), uniProtKBID=as.character(), proteinName=as.character(), SkyDocumentName=as.character())

# Separate the error detecting codes from the warning detecting codes.
# Traverse the SkyDocumentName in fileDf to detect all the possible errors.
# Create a list to store the  peptides with errors for each SkyDocumentName.
df_skydoc_error_peptide <- list()
df_skydoc_error_peptide2 <- list()
for (SkyDocumentName in as.character(fileDf[, "SkyDocumentName"])) {
    labkey.data1 <- labkey.data.total[labkey.data.total$SkyDocumentName==SkyDocumentName, ]
    mypeptideType <- labkey.data1$peptide_standard_purity[1]
    # Get a list of all peptides
    peptide_list <- unique(labkey.data1[ , 'peptidemodifiedsequence'])
    peptide_list_with_error <- c()
    peptide_list_with_error2 <- c()
    for (input_peptide_sequence in peptide_list) {
        labkey.data2 <- labkey.data1[labkey.data1$peptidemodifiedsequence==input_peptide_sequence, ]
        # Get a list of all protein names, although usually one peptide with a specific precursor charge has only one protein.
        protein_list <- unique(labkey.data2[ , 'protein'])
        protein_uniProtID_list <- sapply(protein_list, identify_uniProtKB_entryID, USE.NAMES=FALSE)
        for (indexLabel in 1:length(protein_list)) {
            input_protein_name <- protein_list[indexLabel]
            protein_uniProtID <- protein_uniProtID_list[indexLabel]
            # Choose the specific peptide with a specific precursor charge from a specific protein.
            labkey.data <- labkey.data2[labkey.data2$protein==input_protein_name, ]
            internal_standards <- unique(labkey.data$internal_standard)
            
            if (nrow(labkey.data) >= 1) {
                for (precursorchargeTmp in unique(labkey.data$precursorcharge)) {
                    labkey.data.tmp <- labkey.data[labkey.data$precursorcharge==precursorchargeTmp, ]
                    isotopeLabelTypeTmp <- paste(sort(unique(labkey.data.tmp$isotopelabel)), collapse = '|')
                    transitionTmp <- c()
                    for (isotopeLabelTypeSubtmp in sort(unique(labkey.data.tmp$isotopelabel))) {
                        transitionTmp <- c(transitionTmp, paste(isotopeLabelTypeSubtmp, paste(sort(unique(labkey.data.tmp[labkey.data.tmp$isotopelabel==isotopeLabelTypeSubtmp, ]$fragment_ion_complete)), collapse = '|'), sep=':'))
                    }
                    transitionTmp <- paste(transitionTmp, collapse = ';')
                    logdfTmp <- data.frame(peptide=input_peptide_sequence, precursorCharge=precursorchargeTmp, isotopeLabelType=isotopeLabelTypeTmp, transition= transitionTmp, uniProtKBID=protein_uniProtID, proteinName=input_protein_name, SkyDocumentName=SkyDocumentName)
                    logdf <- rbind(logdf, logdfTmp)
                }
            }
            
            # It's impossile that there will be more than 1 internal standard for one specific SkyDocumentName. So the following code will be commented.
            #if (length(internal_standards) > 1) {
            #    errorType <- "Error"
            #    errorSubtype <- "Internal standard"
            #    errorReason <- paste("Multiple internal standard types found: ", paste(internal_standards, collapse='; '), '.', sep='')
            #    errorInfor <- paste(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', '', '', '', '', '', '', '', '', '', '', '', '', '',sep='\t')
            #    #cat(errorInfor)
            #    errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence), rep('', colNumber-3)), collapse='\t')
            #    cat('\n')
            #    peptide_list_with_error <- c(peptide_list_with_error, input_peptide_sequence)
            #    next
            #}
            
            # Remove the rows whose donotuse value is set.
            labkey.data$donotuse[is.na(labkey.data$donotuse)] <- "FALSE"
            labkey.data$donotuse[tolower(labkey.data$donotuse) == "x"] <- "TRUE"
            labkey.data <- labkey.data[tolower(labkey.data$donotuse) != "true",]
                 
            if (!is.na(labkey.data$multiplicationfactor[1])){
                labkey.data$concentration <- labkey.data$concentration * labkey.data$multiplicationfactor
            }
            
            if(length(unique(labkey.data$concentration)) <2) {
                errorType <- "Error"
                errorSubtype <- "Concentration"
                errorReason <- paste("More than one concentration levels are needed: ", paste(unique(labkey.data$concentration), collapse='; '), ".", sep='')
                #errorInfor <- paste(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', '', '', '', '', '', '', '', '', '', '', '', '', '',sep='\t')
                errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence), rep('', colNumber-3)), collapse='\t')
                cat(errorInfor)
                cat('\n')
                peptide_list_with_error <- c(peptide_list_with_error, input_peptide_sequence)
                peptide_list_with_error2 <- c(peptide_list_with_error2, input_peptide_sequence)
                next
            }
            
            #if (is.na(labkey.data$isspike[1])){
            #    labkey.data$isspike <- labkey.data$peptideconcentrationis
            #}
            
            labkey.data$area[is.na(labkey.data$area)] <- 0 
            if (length(grep("background", names(labkey.data))) >0){
                labkey.data$background[is.na(labkey.data$background)] <- 0 
                labkey.data$rawArea <- labkey.data$area + labkey.data$background 
            }else{
                labkey.data$rawArea <- labkey.data$area 
            }
            
            # Evaluate the function of dcast to capture the warning message: Aggregation function missing: defaulting to length
            evaOut <- evaluate("dcast(labkey.data, protein + peptidemodifiedsequence + precursorcharge + productcharge + fragmention + replicate + concentration + samplegroup + isspike ~ isotopelabel, value.var='rawArea')")
            if (length(evaOut) == 3) {
                # In this condition, warning message "Aggregation function missing: defaulting to length" is printed
                errorType <- "Error"
                errorSubtype <- "Area values of heavy or light Isotope"
                errorReason <- paste("More than one area values exist for the combination of protein, peptidemodifiedsequence, precursorcharge, productcharge, fragmention, replicate, concentration, samplegroup, isspike and isotopelabel.")
                #errorInfor <- paste(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', '', '', '', '', '', '', '', '', '', '', '', '', '',sep='\t')
                errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence), rep('', colNumber-3)), collapse='\t')
                cat(errorInfor)
                cat('\n')
                peptide_list_with_error <- c(peptide_list_with_error, input_peptide_sequence)
                peptide_list_with_error2 <- c(peptide_list_with_error2, input_peptide_sequence)
                next
            }
            
            df <- dcast(labkey.data, protein + peptidemodifiedsequence + precursorcharge + productcharge + fragmention + replicate + concentration + samplegroup + isspike ~ isotopelabel, value.var="rawArea")
            colnames(df) <-  c("ProteinName", "PeptideModifiedSequence", "PrecursorCharge", "ProductCharge", "FragmentIon", "Replicate", "Concentration", "SampleGroup", "ISSpike", "heavyArea", "lightArea")

            df$lightArea[is.na(df$heavyArea)] <- NA
            df$heavyArea[is.na(df$lightArea)] <- NA
            df$heavyArea[df$lightArea==0] <- NA

            if (internal_standards[1] != "heavy"){
               df$HLRatio <- df$heavyArea/df$lightArea
            }
            if (internal_standards[1] == "heavy"){
               df$HLRatio <- df$lightArea/df$heavyArea
            }
            
            df <- df[is.finite(df$HLRatio),]
            
            # Calculate the number of precursorcharge.fragmention.productcharge in df, if it is less than 2, throw an error. There should be at least two fragment ions kept in the response curve.
            df$Fragment_Ion <- paste(df[ ,'PrecursorCharge'], df[ ,'FragmentIon'], df[,'ProductCharge'], sep=".")
            labkey.data$fragment_ion <- paste(labkey.data[ ,'precursorcharge'], labkey.data[ ,'fragmention'], labkey.data[,'productcharge'], sep=".")
            s1 <- unique(df$Fragment_Ion)
            s2 <- unique(labkey.data$fragment_ion)
            if (length(s1) == 0) {
                errorType <- "Error"
                errorSubtype <- "Fragment ion"
                errorReasonTmp <- "In response curve, no fragment ion with both heavy and light isotope exists."
                omittedFragementIon <- c()
                omittedFragementIonIsotope <- c()
                for (item in s2) {
                    omittedFragementIon <- c(omittedFragementIon, item)
                    omittedFragementIonIsotope <- c(omittedFragementIonIsotope, unique(labkey.data[labkey.data$fragment_ion == item, ]$isotopelabel)[1])
                }
                if (length(omittedFragementIon) == 0) {
                    errorReason <- errorReasonTmp
                } else if (length(omittedFragementIon) == 1) {
                    errorReason <- paste(errorReasonTmp, "The other fragment ion with one type of isotope is: ", paste(paste(omittedFragementIon, omittedFragementIonIsotope, sep='-'), collapse = ","), '.', sep="")
                } else if (length(omittedFragementIon) > 1) {
                    errorReason <- paste(errorReasonTmp, "The other fragment ions with one type of isotope are: ", paste(paste(omittedFragementIon, omittedFragementIonIsotope, sep='-'), collapse = ","), '.', sep="")
                }
                #errorInfor <- paste(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', '', '', '', '', '', '', '', '', '', '', '', '', '',sep='\t')
                errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence), rep('', colNumber-3)), collapse='\t')
                cat(errorInfor)
                cat('\n')
                peptide_list_with_error <- c(peptide_list_with_error, input_peptide_sequence)
                peptide_list_with_error2 <- c(peptide_list_with_error2, input_peptide_sequence)
                next
            }
            # delete column of Fragment_Ion in df and column of fragment_ion in labkey.data
            df[c('Fragment_Ion')] <- list(NULL)
            labkey.data[c('fragment_ion')] <- list(NULL)
            
            df2 <- dcast(labkey.data, protein + peptidemodifiedsequence + precursorcharge + productcharge + fragmention + replicate + concentration + samplegroup + isspike ~ isotopelabel, value.var="rawArea")
            colnames(df2) <-  c("ProteinName", "PeptideModifiedSequence", "PrecursorCharge", "ProductCharge", "FragmentIon", "Replicate", "Concentration", "SampleGroup", "ISSpike", "heavyArea", "lightArea")
            df2$lightArea[df2$heavyArea ==0] <- NA;
            df2$heavyArea[df2$lightArea ==0] <- NA;
            
            if (internal_standards[1]  != "heavy"){
               df2$Ratio <- df2$heavyArea/df2$lightArea
            }
            if (internal_standards[1] == "heavy"){
               df2$Ratio <- df2$lightArea/df2$heavyArea
            }
            
            df2 <- df2[is.finite(df2$Ratio),]
            
            if (internal_standards[1] != "heavy"){
                TSum= ddply(df2, .(ProteinName, PeptideModifiedSequence, Concentration, SampleGroup, Replicate, ISSpike), summarize, Ratio = sum(heavyArea, na.rm=TRUE)/sum(lightArea, na.rm=TRUE), heavyArea=sum(heavyArea, na.rm=TRUE), lightArea=sum(lightArea, na.rm=TRUE))
            }
            if (internal_standards[1] == "heavy"){
                TSum = ddply(df2, .(ProteinName, PeptideModifiedSequence, Concentration, SampleGroup, Replicate, ISSpike), summarize, Ratio = sum(lightArea, na.rm=TRUE)/sum(heavyArea, na.rm=TRUE), heavyArea=sum(heavyArea, na.rm=TRUE), lightArea=sum(lightArea, na.rm=TRUE))
            }
            
            TSum$FragmentIon <- "SUM"
            TSum$PrecursorCharge <- ""
            TSum$ProductCharge <- ""
            
            TSum <- TSum[, names(df2)]
            
            df2 <- rbind(df2, TSum)
            
            if (all(df2$ISSpike == 0)) {
                df2$MeasuredConcentration <- df2$Ratio
            } else {
                df2$MeasuredConcentration <- df2$Ratio * df2$ISSpike
            }

            # Set the PrecursorCharge of the SUM in df2.
            precursorChargeSet <- c()
            for (precursorChargeTmp in unique(df2$PrecursorCharge)) {
              if (precursorChargeTmp != '') {
                precursorChargeSet <- c(precursorChargeSet, precursorChargeTmp) 
              }
            }
            precursorChargeSet <- precursorChargeSet[1]
            df2[df2$FragmentIon == 'SUM', ]$PrecursorCharge <- precursorChargeSet
            
            curveDataIndex <- with(df2,  order(ProteinName, PeptideModifiedSequence, PrecursorCharge, FragmentIon, ProductCharge, Concentration, Replicate))
            thisPeptide <- df2[curveDataIndex,]
            thisPeptide$PrecursorCharge <- substr(thisPeptide$PrecursorCharge,1,1)
            thisPeptide$ProductCharge <- substr(thisPeptide$ProductCharge,1,1)
            uniquePeptide <- unique(thisPeptide$PeptideModifiedSequence)
            
            # Calculate LOD/LOQ. This part is omitted.
            thisPeptide$FragmentIon <- paste(thisPeptide$PrecursorCharge, thisPeptide$FragmentIon, thisPeptide$ProductCharge, sep=".") 
            # Keep the fragment ions in keptFragmentIon
            result= ddply(df2, .(ProteinName, PeptideModifiedSequence, PrecursorCharge, FragmentIon, ProductCharge, Concentration, SampleGroup), 
                summarize, MedianR=median(Ratio, na.rm= TRUE), MinR = min(Ratio, na.rm=TRUE), MaxR=max(Ratio, na.rm=TRUE), CVR=sd(Ratio, na.rm= TRUE)/mean(Ratio, na.rm= TRUE), MedianMeasuredC=median(MeasuredConcentration, na.rm= TRUE), SDMeasuredC=sd(MeasuredConcentration, na.rm= TRUE), MinMeasuredC=min(MeasuredConcentration, na.rm= TRUE), MaxMeasuredC=max(MeasuredConcentration, na.rm= TRUE))

            result$Median[is.na(result$MedianR)] <- 0
            
            curveDataIndex <- with( result,  order(ProteinName, PeptideModifiedSequence, PrecursorCharge, FragmentIon, ProductCharge, Concentration))
            thisPeptide <- result[curveDataIndex,]
            thisPeptide$PrecursorCharge <- substr(thisPeptide$PrecursorCharge,1,1)
            thisPeptide$ProductCharge <- substr(thisPeptide$ProductCharge,1,1)
            uniquePeptide <- unique(thisPeptide$PeptideModifiedSequence)
            thisPeptide$FragmentIon <- paste(thisPeptide$PrecursorCharge, thisPeptide$FragmentIon, thisPeptide$ProductCharge, sep=".")
            thisPeptide$MaxMeasuredC[thisPeptide$MaxMeasuredC > max(thisPeptide$MedianMeasuredC)*2] <- max(thisPeptide$MedianMeasuredC)

            #Calculate the robust linear fitting model of concentration of analyte vs.  median of measured concentration of analyte.
            uniqueT <- unique(thisPeptide$FragmentIon)
            thisPeptide <- thisPeptide[thisPeptide$Concentration >0,]
            thisPeptide <- thisPeptide[with(thisPeptide, order(FragmentIon, Concentration)),]
            thisPeptideR <- thisPeptide[(thisPeptide$MedianMeasuredC != 0 & is.finite(thisPeptide$MedianMeasuredC)),c("FragmentIon", "Concentration", "MedianMeasuredC")]
            fitR <- NULL
            for (j in 1:length(uniqueT)){
                thisFragmentIon <- thisPeptideR[(thisPeptideR$FragmentIon == uniqueT[j]),]
                x <- thisFragmentIon$Concentration
                y <- thisFragmentIon$MedianMeasuredC
                w <- 1/(y)^2
                tfit <- try(rlm(x~y, weights=w, method="MM", maxit=1000), silent=TRUE)
                if (!inherits ( tfit, "try-error")){
                    mCoef <- coef(summary(tfit))
                    r2 <- (cor(x,y))^2
                    fitR <- rbind(fitR, c(uniqueT[j], mCoef[2,1], mCoef[1,1], mCoef[2,2], mCoef[1,2],r2))
                }else{
                   #print("Regression fit failed")
                   invisible()
                }
            }
            if (is.null(fitR)) {
                # It means regression fit of all the fragment ions fail.
                # However, the information of the peptide is sufficient to infer internal standard.
                errorType <- "Error"
                errorSubtype <- "Internal standard peptide concentration"
                errorReason <- "The annotation of the concentration of the internal standard peptide has issues. Please check the annotation of IS Spike or PeptideConcentrationIS"
                errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence), rep('', colNumber-3)), collapse='\t')
                cat(errorInfor)
                cat('\n')
                peptide_list_with_error <- c(peptide_list_with_error, input_peptide_sequence)
                next
            }
        }
    }
    peptide_list_with_error <- unique(peptide_list_with_error)
    peptide_list_with_error2 <- unique(peptide_list_with_error2)
    # if peptide_list_with_error is c(), the SkyDocumentName will not be inserted.
    df_skydoc_error_peptide[[SkyDocumentName]] <- peptide_list_with_error
    df_skydoc_error_peptide2[[SkyDocumentName]] <- peptide_list_with_error2
}

if (plot_output) {
    write.table(logdf, file=log_filename, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
}

# Infer the internal standard type for each SkyDocumentName by randomly sampled 5 peptides.
is_inferred_filename <- paste(plot_output_dir, "\\internal_standard_inferred_infor.tsv", sep='' )
df_internal_standard_inferred <- data.frame(SkyDocumentName=as.character(), internal_standard=as.character())
for (SkyDocumentName in as.character(fileDf[, "SkyDocumentName"])) {
    labkey.data1 <- labkey.data.total[labkey.data.total$SkyDocumentName==SkyDocumentName, ]
    peptide_list <- unique(labkey.data1[ , 'peptidemodifiedsequence'])
    # Remove the peptides with errors.
    if (SkyDocumentName %in% names(df_skydoc_error_peptide2)) {
        peptide_list <- setdiff(peptide_list,df_skydoc_error_peptide2[[SkyDocumentName]])
    }
    # Randomly sampling 5 peptides. If the number of the peptides is less than 5, keep all of them.
    if (length(peptide_list) >= 5) {
        peptide_list_sampled <- sample(peptide_list, 5, replace = FALSE)
    } else {
        peptide_list_sampled <- peptide_list
    }
    # Set internal_standard_hypothetical to be light
    if (length(peptide_list_sampled) == 0) {
        internal_standard_inferred <- "can't be inferred"
    } else {
        internal_standard_hypothetical <- 'light'
        value1 <- 0
        value2 <- 0
        
        for (input_peptide_sequence in peptide_list_sampled) {
            labkey.data2 <- labkey.data1[labkey.data1$peptidemodifiedsequence==input_peptide_sequence, ]
            protein_list <- unique(labkey.data2[ , 'protein'])
            protein_uniProtID_list <- sapply(protein_list, identify_uniProtKB_entryID, USE.NAMES=FALSE)
            for (indexLabel in 1:length(protein_list)) {
                input_protein_name <- protein_list[indexLabel]
                protein_uniProtID <- protein_uniProtID_list[indexLabel]
                labkey.data <- labkey.data2[labkey.data2$protein==input_protein_name, ]
                internal_standards <- internal_standard_hypothetical
                # Remove the rows whoes donotuse value is set.
                labkey.data$donotuse[is.na(labkey.data$donotuse)] <- "FALSE"
                labkey.data$donotuse[tolower(labkey.data$donotuse) == "x"] <- "TRUE"
                labkey.data <- labkey.data[tolower(labkey.data$donotuse) != "true",]
                
                if (!is.na(labkey.data$multiplicationfactor[1])){
                    labkey.data$concentration <- labkey.data$concentration * labkey.data$multiplicationfactor
                }
                
                #if (is.na(labkey.data$isspike[1])){
                #   labkey.data$isspike <- labkey.data$peptideconcentrationis
                #}
                labkey.data$area[is.na(labkey.data$area)] <- 0 
                if (length(grep("background", names(labkey.data))) >0){
                    labkey.data$background[is.na(labkey.data$background)] <- 0 
                    labkey.data$rawArea <- labkey.data$area + labkey.data$background 
                }else{
                    labkey.data$rawArea <- labkey.data$area 
                }
                
                df <- dcast(labkey.data, protein + peptidemodifiedsequence + precursorcharge + productcharge + fragmention + replicate + concentration + samplegroup + isspike ~ isotopelabel, value.var="rawArea")
                colnames(df) <-  c("ProteinName", "PeptideModifiedSequence", "PrecursorCharge", "ProductCharge", "FragmentIon", "Replicate", "Concentration", "SampleGroup", "ISSpike", "heavyArea", "lightArea")
    
                df$lightArea[is.na(df$heavyArea)] <- NA
                df$heavyArea[is.na(df$lightArea)] <- NA
                df$heavyArea[df$lightArea==0] <- NA
                # Because internal_standard_hypothetical is set to be light, HLRatio = heavyArea/lightArea
                df$HLRatio <- df$heavyArea/df$lightArea
                
                temp <- ddply(df, .(ProteinName, PeptideModifiedSequence, PrecursorCharge, FragmentIon, ProductCharge), summarize, medianArea=median(heavyArea, na.rm=TRUE))
                orderT <- with(temp,  order(ProteinName, PeptideModifiedSequence, -medianArea))
                # Pay attention: Only top three fragment ions with the largest medianArea are kept.
                df <- merge(df, temp[orderT[1:3], ])
                result= ddply(df, .(ProteinName, PeptideModifiedSequence, PrecursorCharge, FragmentIon, ProductCharge, Concentration, SampleGroup), 
                    summarize, Median=median(HLRatio, na.rm= TRUE), Min = min(HLRatio, na.rm=TRUE), Max=max(HLRatio, na.rm=TRUE))
                result$Median[is.na(result$Median)] <- 0
                
                curveDataIndex <- with( result,  order(ProteinName, PeptideModifiedSequence, PrecursorCharge, FragmentIon, ProductCharge, Concentration))
                thisPeptide <- result[curveDataIndex,]
                thisPeptide$PrecursorCharge <- substr(thisPeptide$PrecursorCharge,1,1)
                thisPeptide$ProductCharge <- substr(thisPeptide$ProductCharge,1,1)
                uniquePeptide <- unique(thisPeptide$PeptideModifiedSequence)
                thisPeptide$FragmentIon <- paste(thisPeptide$PrecursorCharge, thisPeptide$FragmentIon, thisPeptide$ProductCharge, sep=".")
                
                fragmentIon_list <- unique(thisPeptide$FragmentIon)
                for (fragmentIon_tmp in fragmentIon_list) {
                    thisFragmentIon <- thisPeptide[thisPeptide$FragmentIon == fragmentIon_tmp, ]
                    minConcentration <- min(thisFragmentIon$Concentration)
                    maxConcentration <- max(thisFragmentIon$Concentration)
                    medHLRation_min <- median(thisFragmentIon[thisFragmentIon$Concentration == minConcentration, ]$Median, na.rm=TRUE)
                    medHLRation_max <- median(thisFragmentIon[thisFragmentIon$Concentration == maxConcentration, ]$Median, na.rm=TRUE)
                    if (is.na(medHLRation_min) | is.na(medHLRation_max)) {
                        invisible()
                    } else if (medHLRation_max >= medHLRation_min) {
                        value1 <- value1 + 1
                    } else {
                        value2 <- value2 + 1
                    }
                }
            }
        }
        #print(paste('value1:', value1, 'value2:', value2, sep=' '))
        if (value1==0 & value2==0) {
            internal_standard_inferred <- "can't be inferred"
        } else if (value1 >= value2) {
            internal_standard_inferred <- 'light'
        } else {
            internal_standard_inferred <- 'heavy'
        }
    }
    df_internal_standard_inferred_tmp <- data.frame(SkyDocumentName=SkyDocumentName, internal_standard=internal_standard_inferred)
    df_internal_standard_inferred <- rbind(df_internal_standard_inferred, df_internal_standard_inferred_tmp)
}

write.table(df_internal_standard_inferred, file=is_inferred_filename, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

# sequece <- 'AAAAAAGAGPEMVR'
# Traverse the SkyDocumentName in fileDf to detect all the possible warnings.
for (SkyDocumentName in as.character(fileDf[, "SkyDocumentName"])) {
    # Evaluate the internal_standard, if the internal standard is wrong, errors will arise.
    original_internal_standard <- as.character(fileDf[fileDf$SkyDocumentName == SkyDocumentName, ]$internal_standard)
    inferred_internal_standard <- as.character(df_internal_standard_inferred[df_internal_standard_inferred$SkyDocumentName == SkyDocumentName, ]$internal_standard)
    if (inferred_internal_standard[1] == "can't be inferred") {
        # In this situation, use the original internal standard for calculations.
        #next
        invisible()
    }
    if (original_internal_standard[1] == 'none') {
        # Just jump out of the loop. Don't print the errorInfor, because it has already be printed in the function of detectIS in qcAnalysis.py
        next
    }
    
    if ((inferred_internal_standard[1] != "can't be inferred") & (original_internal_standard[1] != inferred_internal_standard[1])) {
        errorType <- "Error"
        errorSubtype <- "Internal standard"
        errorReason <- paste('The internal standard in the skyline file is set to be ', original_internal_standard, ', while the inferred internal standard is ', inferred_internal_standard, '.', sep='')
        #errorInfor <- paste(SkyDocumentName, errorType, errorSubtype, errorReason, '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', sep='\t')
        errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason), rep('', colNumber-1)), collapse='\t')
        cat(errorInfor)
        cat('\n')
        next
    }

    labkey.data1 <- labkey.data.total[labkey.data.total$SkyDocumentName==SkyDocumentName, ]
    mypeptideType <- labkey.data1$peptide_standard_purity[1]
    # Get a list of all peptides
    peptide_list <- unique(labkey.data1[ , 'peptidemodifiedsequence'])
    # Remove the peptides with errors.
    if (SkyDocumentName %in% names(df_skydoc_error_peptide)) {
        peptide_list <- setdiff(peptide_list,df_skydoc_error_peptide[[SkyDocumentName]])
    }
    for (input_peptide_sequence in peptide_list) {
        labkey.data2 <- labkey.data1[labkey.data1$peptidemodifiedsequence==input_peptide_sequence, ]
        # Get a list of all protein names, although usually one peptide with a specific precursor charge has only one protein.
        protein_list <- unique(labkey.data2[ , 'protein'])
        protein_uniProtID_list <- sapply(protein_list, identify_uniProtKB_entryID, USE.NAMES=FALSE)
        for (indexLabel in 1:length(protein_list)) {
            input_protein_name <- protein_list[indexLabel]
            protein_uniProtID <- protein_uniProtID_list[indexLabel]
            # Choose the specific peptide with a specific precursor charge from a specific protein.
            labkey.data <- labkey.data2[labkey.data2$protein==input_protein_name, ]
            internal_standards <- unique(labkey.data$internal_standard)
                        
            # Remove the rows whoes donotuse value is set.
            labkey.data$donotuse[is.na(labkey.data$donotuse)] <- "FALSE"
            labkey.data$donotuse[tolower(labkey.data$donotuse) == "x"] <- "TRUE"
            labkey.data <- labkey.data[tolower(labkey.data$donotuse) != "true",]
            
            if (!is.na(labkey.data$multiplicationfactor[1])){
                labkey.data$concentration <- labkey.data$concentration * labkey.data$multiplicationfactor
            }
            
            #if (is.na(labkey.data$isspike[1])){
            #    labkey.data$isspike <- labkey.data$peptideconcentrationis
            #}
            labkey.data$area[is.na(labkey.data$area)] <- 0 
            if (length(grep("background", names(labkey.data))) >0){
                labkey.data$background[is.na(labkey.data$background)] <- 0 
                labkey.data$rawArea <- labkey.data$area + labkey.data$background 
            }else{
                labkey.data$rawArea <- labkey.data$area 
            }
            
            df <- dcast(labkey.data, protein + peptidemodifiedsequence + precursorcharge + productcharge + fragmention + replicate + concentration + samplegroup + isspike ~ isotopelabel, value.var="rawArea")
            colnames(df) <-  c("ProteinName", "PeptideModifiedSequence", "PrecursorCharge", "ProductCharge", "FragmentIon", "Replicate", "Concentration", "SampleGroup", "ISSpike", "heavyArea", "lightArea")

            df$lightArea[is.na(df$heavyArea)] <- NA
            df$heavyArea[is.na(df$lightArea)] <- NA
            df$heavyArea[df$lightArea==0] <- NA

            if (internal_standards[1] != "heavy"){
               df$HLRatio <- df$heavyArea/df$lightArea
            }
            if (internal_standards[1] == "heavy"){
               df$HLRatio <- df$lightArea/df$heavyArea
            }
            
            df <- df[is.finite(df$HLRatio),]
            
            # Calculate the number of precursorcharge.fragmention.productcharge in df, if it is less than 2, throw an warning. There should be at least two fragment ions kept in the response curve.
            df$Fragment_Ion <- paste(df[ ,'PrecursorCharge'], df[ ,'FragmentIon'], df[,'ProductCharge'], sep=".")
            labkey.data$fragment_ion <- paste(labkey.data[ ,'precursorcharge'], labkey.data[ ,'fragmention'], labkey.data[,'productcharge'], sep=".")
            s1 <- unique(df$Fragment_Ion)
            s2 <- unique(labkey.data$fragment_ion)
            if (length(s1) < 3) {
                errorType <- "Warning"
                errorSubtype <- "Fragment ion less than three"
                if (length(s1) == 1) {
                    errorReasonTmp <- paste("In response curve, only one fragment ion ", s1[1] , " (with both heavy and light isotopes) exists.", sep="")
                } else {
                    errorReasonTmp <- paste("In response curve, only two fragment ions ", s1[1] , ", " , s1[2] , " (with both heavy and light isotopes) exist.", sep="")
                }
                omittedFragementIon <- c()
                omittedFragementIonIsotope <- c()
                for (item in s2) {
                    if (!(item %in% s1)) {
                        omittedFragementIon <- c(omittedFragementIon, item)
                        omittedFragementIonIsotope <- c(omittedFragementIonIsotope, unique(labkey.data[labkey.data$fragment_ion == item, ]$isotopelabel)[1])
                    }
                }
                if (length(omittedFragementIon) == 0) {
                    errorReason <- errorReasonTmp
                } else if (length(omittedFragementIon) == 1) {
                    errorReason <- paste(errorReasonTmp, "The other fragment ion with one type of isotope is: ", paste(paste(omittedFragementIon, omittedFragementIonIsotope, sep='-'), collapse = ","), '.', sep="")
                } else if (length(omittedFragementIon) > 1) {
                    errorReason <- paste(errorReasonTmp, "The other fragment ions with one type of isotope are: ", paste(paste(omittedFragementIon, omittedFragementIonIsotope, sep='-'), collapse = ","), '.', sep="")
                }
                #errorInfor <- paste(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', '', '', '', '', '', '', '', '', '', '', '', '', '',sep='\t')
                errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence), rep('', colNumber-3)), collapse='\t')
                cat(errorInfor)
                cat('\n')
            }
            
            # delete column of Fragment_Ion in df and column of fragment_ion in labkey.data
            df[c('Fragment_Ion')] <- list(NULL)
            labkey.data[c('fragment_ion')] <- list(NULL)
            
            if (internal_standards[1] != "heavy"){
                TSum = ddply(df, .(ProteinName, PeptideModifiedSequence, PrecursorCharge, Concentration, SampleGroup, Replicate), summarize, HLRatio = sum(heavyArea, na.rm=TRUE)/sum(lightArea, na.rm=TRUE))
            }
            if (internal_standards[1] == "heavy"){
                TSum = ddply(df, .(ProteinName, PeptideModifiedSequence, PrecursorCharge, Concentration, SampleGroup, Replicate), summarize, HLRatio = sum(lightArea, na.rm=TRUE)/sum(heavyArea, na.rm=TRUE))
            }
            
            TSum$FragmentIon <- "SUM"
            TSum$ProductCharge <- ""
            
            resultSum = ddply(TSum, .(ProteinName, PeptideModifiedSequence, PrecursorCharge, FragmentIon, ProductCharge, Concentration, SampleGroup), 
                summarize, Median=median(HLRatio, na.rm= TRUE), Min = min(HLRatio, na.rm=TRUE), Max=max(HLRatio, na.rm=TRUE), CV=sd(HLRatio, na.rm= TRUE)/mean(HLRatio, na.rm= TRUE), LOD=mean(HLRatio, na.rm= TRUE)+ 3* sd(HLRatio, na.rm= TRUE)) 
            
            if (internal_standards[1] != "heavy"){
                temp = ddply(df, .(ProteinName, PeptideModifiedSequence, PrecursorCharge, FragmentIon, ProductCharge), summarize, medianArea=median(heavyArea, na.rm=TRUE))
            }

            if (internal_standards[1] == "heavy"){
                temp = ddply(df, .(ProteinName, PeptideModifiedSequence, PrecursorCharge, FragmentIon, ProductCharge), summarize, medianArea=median(lightArea, na.rm=TRUE))
            }
            
            orderT <- with(temp,  order(ProteinName, PeptideModifiedSequence, -medianArea))
            # Pay attention: Only top three fragment ions with the largest medianArea are kept.
            # However, in the fitTable and LODTable, all the fragment ions are listed there, so it's better to keep the top three fragment ions in the tables.
            df <- merge(df, temp[orderT[1:3], ])
            resultT= ddply(df, .(ProteinName, PeptideModifiedSequence, PrecursorCharge, FragmentIon, ProductCharge, Concentration, SampleGroup), 
                summarize, Median=median(HLRatio, na.rm= TRUE), Min = min(HLRatio, na.rm=TRUE), Max=max(HLRatio, na.rm=TRUE), CV=sd(HLRatio, na.rm= TRUE)/mean(HLRatio, na.rm= TRUE), LOD=mean(HLRatio, na.rm= TRUE)+ 3* sd(HLRatio, na.rm= TRUE))

            result <- rbind(resultT, resultSum)
            result$Median[is.na(result$Median)] <- 0
            
            curveDataIndex <- with( result,  order(ProteinName, PeptideModifiedSequence, PrecursorCharge, FragmentIon, ProductCharge, Concentration))
            thisPeptide <- result[curveDataIndex,]
            thisPeptide$PrecursorCharge <- substr(thisPeptide$PrecursorCharge,1,1)
            thisPeptide$ProductCharge <- substr(thisPeptide$ProductCharge,1,1)
            uniquePeptide <- unique(thisPeptide$PeptideModifiedSequence)
            
            if (plot_output) {
                thisPeptide$FragmentIon <- paste(thisPeptide$PrecursorCharge, thisPeptide$FragmentIon, thisPeptide$ProductCharge, sep=".")
                keptFragmentIon <- unique(thisPeptide$FragmentIon)
                samePeptideLength <- dim(thisPeptide)[1]
                mProtein <- unlist(strsplit(as.character(thisPeptide$ProteinName[1]),"[.]"))[1]
                #mTitle <- paste("Analyte: ", mProtein, ".", uniquePeptide[1], "\n", sep="")
                mTitle <- paste("Analyte: ", protein_uniProtID, ".", uniquePeptide[1], "\n", sep="")
                thisPeptide$Max[thisPeptide$Max > max(thisPeptide$Median)*2] <- max(thisPeptide$Median)
                for (myplotType in c("linear", "log", "residual")) {
                    if ( tolower(myplotType) == "linear") {
                        mxlabel <- "\nTheoretical Concentration (fmol/uL)"
                        mcolor <- "black"
                        if (tolower(mypeptideType) == "crude"){
                            mxlabel <- "\nTheoretical Concentration (fmol/uL)\nEstimated from unpurified peptide"
                            mcolor <- "red"
                        }
                        CairoPNG(filename=paste(plot_output_dir, "\\", input_peptide_sequence, "_", myplotType, '_', indexLabel, '_ResponseCurveQuery.response_curve.png', sep=''), width=800, height=600, bg="white")
                        p <- ggplot(data=thisPeptide, aes(x=Concentration, y=Median, color=FragmentIon)) + geom_errorbar(aes(ymin=Min, ymax=Max), width=.8) + geom_smooth(method=lm, se=FALSE) +geom_point(size=2) + xlab(mxlabel) + ylab("Peak Area Ratio") + theme(title=element_text(size=18, colour="black"), axis.text=element_text(size=16), axis.text.x=element_text(colour=mcolor), axis.title=element_text(size=20) , axis.title.x=element_text(colour=mcolor), legend.position = c(0.15, 0.85), legend.title = element_text(size=14), legend.text = element_text(size=14))+ labs(title=mTitle)+ scale_colour_discrete(name = "Transition")
                        print(p)
                        dev.off()
                        
                        # Extract the slopes for each FragmentIon
                        if (qcSlopeSwitch) {
                            lm_data <- as.data.frame(thisPeptide %>%
                                                          group_by(FragmentIon) %>%
                                                          do({
                                                            mod = lm(Median ~ Concentration, data = .)
                                                            data.frame(Intercept = coef(mod)[1], Slope = coef(mod)[2], rSquare = summary(mod)$r.squared , pValue = summary(mod)$coefficients[2,4])
                                                          })
                                                          )
                            # Judge whether there are outliers among lm_data$Slope according to the coefficient of variance
                            cv_slope <- sd(lm_data$Slope)/mean(lm_data$Slope)
                            if (is.nan(cv_slope)) {
                                errorType <- "Warning"
                                errorSubtype <- "Bad linear regression fitting"
                                errorReason <- paste("The slopes of the fragments ions (", paste(unique(thisPeptide$FragmentIon), collapse=', '), ") is(are) 0.", sep= '')
                                #errorInfor <- paste(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', '', '', '', '', '', '', '', '', '', '', '', '', '',sep='\t')
                                errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence), rep('', colNumber-3)), collapse='\t')
                                cat(errorInfor)
                                cat('\n')
                            } else if (cv_slope > cv_threshold) {
                                errorType <- "Warning"
                                errorSubtype <- "Bad linear regression fitting"
                                errorReason <- paste("The coefficient of variance of slopes of the fragment ions (", paste(unique(thisPeptide$FragmentIon), collapse=', '), ") is(are) larger than ", cv_threshold, ".", sep= '')
                                #errorInfor <- paste(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', '', '', '', '', '', '', '', '', '', '', '', '', '',sep='\t')
                                errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence), rep('', colNumber-3)), collapse='\t')
                                cat(errorInfor)
                                cat('\n')
                            } else {
                                invisible()
                            }
                            # Judge whether the linear fits of the fragment ions are good or not according to the slopes.
                            rSquareMin <- min(lm_data$rSquare)
                            pValueMax <- max(lm_data$pValue)
                            if ( rSquareMin < rSquare_threshold || pValueMax > pValue_threshold) {
                                errorType <- "Warning"
                                errorSubtype <- "Bad linear regression fitting"
                                errorReason <- "The quality of fit of linear model is poor due to the low R2 or large p value of F-test in the process of linear regression."
                                #errorInfor <- paste(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', '', '', '', '', '', '', '', '', '', '', '', '', '',sep='\t')
                                errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence), rep('', colNumber-3)), collapse='\t')
                                cat(errorInfor)
                                cat('\n')
                            }
                        }
                    }
                    if (tolower(myplotType) == "log"){
                        mxlabel <- "\nLog Theoretical Concentration (fmol/uL)"
                        mcolor <- "black"
                        if (tolower(mypeptideType) == "crude"){
                            mxlabel <- "\nLog Theoretical Concentration (fmol/uL)\nEstimated from unpurified peptide"
                            mcolor <- "red"
                         }       
                         CairoPNG(filename=paste(plot_output_dir, "\\", input_peptide_sequence, "_", myplotType, '_', indexLabel, '_ResponseCurveQuery.response_curve.png', sep=''), width=800, height=600, bg="white")
                         pd <- position_dodge(.05)
                         p <- ggplot(data=thisPeptide[thisPeptide$Concentration >0,], aes(x=log(Concentration,10), y=log(Median,10), color=FragmentIon)) + geom_errorbar(aes(ymin=log(Min,10), ymax=log(Max,10)), position=pd, width=.08) +geom_point(position=pd, size=2) + xlab(mxlabel) + ylab("Log Peak Area Ratio") + theme(title=element_text(size=18, colour="black"), axis.text=element_text(size=16),   axis.text.x=element_text(colour=mcolor), axis.title=element_text(size=20) , axis.title.x=element_text( colour=mcolor), legend.position = c(0.2, 0.85), legend.title = element_text(size=14), legend.text = element_text(size=14))+ labs(title=mTitle)+ scale_colour_discrete(name = "Transition")    
                         print(p)
                         dev.off()
                    }
                    if ( tolower(myplotType) == "residual"){
                        mxlabel <- "\nLog Theoretical Concentration (fmol/uL)"
                        mcolor <- "black"
                        if (tolower(mypeptideType) == "crude"){
                            mxlabel <- "\nLog Theoretical Concentration (fmol/uL)\nEstimated from unpurified peptide"
                            mcolor <- "red"
                         }       
                         uniqueT <- unique(thisPeptide$FragmentIon)
                         thisPeptide <- thisPeptide[with(thisPeptide, order(FragmentIon, Concentration)),]
                         thisPeptideR <- thisPeptide[(thisPeptide$Median != 0 & is.finite(thisPeptide$Median)),  ]      
                         for (j in 1:length(uniqueT)){
                             thisFragmentIon <- thisPeptideR[(thisPeptideR$FragmentIon == uniqueT[j]),]		
                             x <- thisFragmentIon$Concentration;
                             y <- thisFragmentIon$Median;
                             if (length(x) > 2){
                                 Lvalue <- fitline(x,y)       
                                 thisPeptideR$Residual[(thisPeptideR$FragmentIon == uniqueT[j])] <- Lvalue$mResidual;
                             }else{
                                 thisPeptideR$Residual[(thisPeptideR$FragmentIon == uniqueT[j])] <- NA
                             }
                         }
                         CairoPNG(filename=paste(plot_output_dir, "\\", input_peptide_sequence, "_", myplotType, '_', indexLabel, '_ResponseCurveQuery.response_curve.png', sep=''), width=800, height=600, bg="white")
                         p <- ggplot(data=thisPeptideR, aes(x=log(Concentration,10), y=Residual, color=FragmentIon)) +geom_point(size=2) + xlab(mxlabel) + ylab("Residual") + theme(title=element_text(size=18, colour="black"), axis.text=element_text(size=16) , axis.text.x=element_text(colour=mcolor), axis.title=element_text(size=20) , axis.title.x=element_text(colour=mcolor), legend.title = element_text(size=14), legend.text = element_text(size=14))+ labs(title=mTitle)+ scale_colour_discrete(name = "Transition")
                         print(p)
                         dev.off()
                    }
                }
                #write.table(thisPeptide, file = paste(plot_output_dir, "\\", input_peptide_sequence, "_", gsub(">", "", gsub("\\|", "-", protein_uniProtID)), "_", SkyDocumentName, '.tsv', sep=''), sep = "\t", qmethod = "double", col.names=NA)
            }
            
            # Run ResponseCurveAnalysis 
            if (RCAnalysisSwitch) {
                df2 <- dcast(labkey.data, protein + peptidemodifiedsequence + precursorcharge + productcharge + fragmention + replicate + concentration + samplegroup + isspike ~ isotopelabel, value.var="rawArea")
                colnames(df2) <-  c("ProteinName", "PeptideModifiedSequence", "PrecursorCharge", "ProductCharge", "FragmentIon", "Replicate", "Concentration", "SampleGroup", "ISSpike", "heavyArea", "lightArea")
                df2$lightArea[df2$heavyArea ==0] <- NA;
                df2$heavyArea[df2$lightArea ==0] <- NA;
                
                if (internal_standards[1]  != "heavy"){
                   df2$Ratio <- df2$heavyArea/df2$lightArea
                }
                if (internal_standards[1] == "heavy"){
                   df2$Ratio <- df2$lightArea/df2$heavyArea
                }
                
                df2 <- df2[is.finite(df2$Ratio),]
                
                if (internal_standards[1] != "heavy"){
                    TSum= ddply(df2, .(ProteinName, PeptideModifiedSequence, Concentration, SampleGroup, Replicate, ISSpike), summarize, Ratio = sum(heavyArea, na.rm=TRUE)/sum(lightArea, na.rm=TRUE), heavyArea=sum(heavyArea, na.rm=TRUE), lightArea=sum(lightArea, na.rm=TRUE))
                }
                if (internal_standards[1] == "heavy"){
                    TSum = ddply(df2, .(ProteinName, PeptideModifiedSequence, Concentration, SampleGroup, Replicate, ISSpike), summarize, Ratio = sum(lightArea, na.rm=TRUE)/sum(heavyArea, na.rm=TRUE), heavyArea=sum(heavyArea, na.rm=TRUE), lightArea=sum(lightArea, na.rm=TRUE))
                }
                
                TSum$FragmentIon <- "SUM"
                TSum$PrecursorCharge <- ""
                TSum$ProductCharge <- ""
                
                TSum <- TSum[, names(df2)]
                
                df2 <- rbind(df2, TSum)
                
                
                if (all(df2$ISSpike == 0)) {
                    errorType <- "Warning"
                    errorSubtype <- "Internal Standard spike peptide concentration"
                    errorReason <- "All of the concentrations of the internal standard peptide are zero. The internal standard is assumed to be the endogenous peptide."
                    #errorInfor <- paste(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', '', '', '', '', '', '', '', '', '', '', '', '', '',sep='\t')
                    errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence), rep('', colNumber-3)), collapse='\t')
                    cat(errorInfor)
                    cat('\n')
                    df2$MeasuredConcentration <- df2$Ratio
                } else {
                    df2$MeasuredConcentration <- df2$Ratio * df2$ISSpike
                }

                # Set the PrecursorCharge of the SUM in df2.
                precursorChargeSet <- c()
                for (precursorChargeTmp in unique(df2$PrecursorCharge)) {
                  if (precursorChargeTmp != '') {
                    precursorChargeSet <- c(precursorChargeSet, precursorChargeTmp) 
                  }
                }
                precursorChargeSet <- precursorChargeSet[1]
                df2[df2$FragmentIon == 'SUM', ]$PrecursorCharge <- precursorChargeSet
                
                curveDataIndex <- with(df2,  order(ProteinName, PeptideModifiedSequence, PrecursorCharge, FragmentIon, ProductCharge, Concentration, Replicate))
                thisPeptide <- df2[curveDataIndex,]
                thisPeptide$PrecursorCharge <- substr(thisPeptide$PrecursorCharge,1,1)
                thisPeptide$ProductCharge <- substr(thisPeptide$ProductCharge,1,1)
                uniquePeptide <- unique(thisPeptide$PeptideModifiedSequence)
                
                if (plot_output) {
                    # Calculate LOD/LOQ
                    thisPeptide$FragmentIon <- paste(thisPeptide$PrecursorCharge, thisPeptide$FragmentIon, thisPeptide$ProductCharge, sep=".") 
                    # Keep the fragment ions in keptFragmentIon
                    
                    lowConc <- sort(unique(thisPeptide$Concentration[thisPeptide$Concentration>0]))[1]
                    usedData <- thisPeptide[,c("FragmentIon", "Concentration", "Replicate", "Ratio")]
                    usedData <- usedData[(usedData$Concentration == 0 | usedData$Concentration == lowConc ), ]

                    LODData <- ddply(usedData,  .(FragmentIon, Concentration), 
                        summarize, mmean=mean(Ratio, na.rm= TRUE), msd=sd(Ratio, na.rm=TRUE), mqt=suppressWarnings(qt(0.95,sum(!is.na(Ratio))-1)))
                    
                    methods <- c("blank+low_conc", "blank_only", "rsd_limit")
                    #first method
                    LODData[LODData$Concentration == lowConc,'mmean'] <- 0
                    LOD1 <- ddply(LODData, .(FragmentIon), summarize, LOD=sum(mmean+msd*mqt))
                    LOD1$LOQ <- 3*LOD1$LOD
                    names(LOD1) <- c("FragmentIon", paste(methods[1],"LOD", sep="_"), paste(methods[1],"LOQ", sep="_"))
                    usedData <- usedData[(usedData$Concentration == 0) ,]
                    if (dim(usedData)[1] >0) {
                        LOD2 <- ddply(usedData, .(FragmentIon), 
                        summarize, LOD=mean(Ratio, na.rm= TRUE)+ 3* sd(Ratio, na.rm= TRUE), LOQ=mean(Ratio, na.rm= TRUE)+ 10 * sd(Ratio, na.rm= TRUE) )
                    
                        names(LOD2) <- c("FragmentIon", paste(methods[2],"LOD", sep="_"), paste(methods[2],"LOQ", sep="_"))
                    } else {
                       LOD2 = LOD1;
                       LOD2[,2:3] <- NA
                       names(LOD2) <- c("FragmentIon", paste(methods[2],"LOD", sep="_"), paste(methods[2],"LOQ", sep="_"))
                    }
                    usedData <- thisPeptide[,c("FragmentIon", "Concentration", "Replicate", "Ratio")]
                    usedData <- usedData[(usedData$Concentration>0 & !is.na(usedData$Ratio)),  ]
                    usedData$estimateConc <- usedData$Concentration * usedData$Ratio;
                    LODData <- ddply(usedData, .(FragmentIon, Concentration), summarize, rsd=sd(estimateConc, na.rm= TRUE)/ mean(estimateConc))
                    LOD3 <- NULL 
                    rsd.max <- 0.15
                    temp <- by(LODData, LODData[,c( "FragmentIon")],
                            function(x){
                                tfit <- try (fit <- lm (log(x$rsd) ~ log (x$Concentration)), silent=TRUE)
                                if (!inherits (tfit, "try-error")) {
                                    # record regression results and LOD/LOQ
                                    intercept <- fit$coefficients[1]
                                    slope <- fit$coefficients[2]
                                    loq <- exp ( (log (rsd.max) - intercept) / slope )
                                    lod <- loq/3
                                } else {
                                    lod <- loq <- NA
                                }
                                LOD3 <<- rbind(LOD3, c(as.character(x[1,1]),lod, loq))
                            })
                    LOD3 <- data.frame(LOD3)
                    names(LOD3) <- c("FragmentIon", paste(methods[3],"LOD", sep="_"), paste(methods[3],"LOQ", sep="_"))
                    LOD3 <- LOD3[with(LOD3, order(FragmentIon)),]
                    LOD3[,2] <- as.numeric(as.character(LOD3[,2]))
                    LOD3[,3] <- as.numeric(as.character(LOD3[,3]))
                    
                    LODTable <- merge(LOD1, LOD2, by="FragmentIon", all=T)
                    LODTable <- merge(LODTable, LOD3,  by="FragmentIon", all=T)           

                    result= ddply(df2, .(ProteinName, PeptideModifiedSequence, PrecursorCharge, FragmentIon, ProductCharge, Concentration, SampleGroup), 
                        summarize, MedianR=median(Ratio, na.rm= TRUE), MinR = min(Ratio, na.rm=TRUE), MaxR=max(Ratio, na.rm=TRUE), CVR=sd(Ratio, na.rm= TRUE)/mean(Ratio, na.rm= TRUE), MedianMeasuredC=median(MeasuredConcentration, na.rm= TRUE), SDMeasuredC=sd(MeasuredConcentration, na.rm= TRUE), MinMeasuredC=min(MeasuredConcentration, na.rm= TRUE), MaxMeasuredC=max(MeasuredConcentration, na.rm= TRUE))

                    result$Median[is.na(result$MedianR)] <- 0
                    
                    curveDataIndex <- with( result,  order(ProteinName, PeptideModifiedSequence, PrecursorCharge, FragmentIon, ProductCharge, Concentration))
                    thisPeptide <- result[curveDataIndex,]
                    thisPeptide$PrecursorCharge <- substr(thisPeptide$PrecursorCharge,1,1)
                    thisPeptide$ProductCharge <- substr(thisPeptide$ProductCharge,1,1)
                    uniquePeptide <- unique(thisPeptide$PeptideModifiedSequence)
                    
                    thisPeptide$FragmentIon <- paste(thisPeptide$PrecursorCharge, thisPeptide$FragmentIon, thisPeptide$ProductCharge, sep=".")
   
                    samePeptideLength <- dim(thisPeptide)[1]
                    mProtein <- unlist(strsplit(as.character(thisPeptide$ProteinName[1]),"[.]"))[1]
                    #mTitle <- paste("Analyte: ", mProtein, ".", uniquePeptide[1], "\n", sep="")
                    mTitle <- paste("Analyte: ", protein_uniProtID, ".", uniquePeptide[1], "\n", sep="")
                    thisPeptide$MaxMeasuredC[thisPeptide$MaxMeasuredC > max(thisPeptide$MedianMeasuredC)*2] <- max(thisPeptide$MedianMeasuredC)
                    
                    # Annotate the following code to inhibit the generation of figures of *._ResponseCurveAnalysis.response_curve.png
                    if (FALSE) {
                        for (myplotType in c("linear", "log")) {
                            if ( tolower(myplotType) == "linear") {
                                mxlabel <- "\nTheoretical Concentration (fmol/ug)"
                                mcolor <- "black"
                                if (mypeptideType == "crude"){
                                    mxlabel <- "\nTheoretical Concentration (fmol/ug)\nEstimated from unpurified peptide"
                                    mcolor <- "red"
                                }
                                CairoPNG(filename=paste(plot_output_dir, "\\", input_peptide_sequence, "_", myplotType, '_', indexLabel, '_ResponseCurveAnalysis.response_curve.png', sep=''), width=800, height=600, bg="white")
                                p<- ggplot(data=thisPeptide, aes(x=Concentration, y=MedianMeasuredC, color=FragmentIon)) + geom_errorbar(aes(ymin=MinMeasuredC, ymax=MaxMeasuredC), width=.8) + geom_smooth(method=lm, se=FALSE) +geom_point(size=2) + xlab(mxlabel) + ylab("Measured Concentration (fmol/ug") + theme(title=element_text(size=18, colour="black"), axis.text=element_text(size=16), axis.text.x=element_text(colour=mcolor), axis.title=element_text(size=20), axis.title.x=element_text(colour=mcolor), legend.position = c(0.15, 0.85), legend.title = element_text(size=14), legend.text = element_text(size=14))+ labs(title=mTitle)+ scale_colour_discrete(name = "Transition")
                                print(p)
                                dev.off()
                            }
                            if (tolower(myplotType) == "log"){
                                mxlabel <- "\nLog Theoretical Concentration (fmol/ug)"
                                mcolor <- "black"
                                if (mypeptideType == "crude"){
                                    mxlabel <- "\nLog Theoretical Concentration (fmol/ug)\nEstimated from unpurified peptide"
                                  mcolor <- "red"
                                }
                                CairoPNG(filename=paste(plot_output_dir, "\\", input_peptide_sequence, "_", myplotType, '_', indexLabel, '_ResponseCurveAnalysis.response_curve.png', sep=''), width=800, height=600, bg="white")
                                pd <- position_dodge(.05)
                                p<- ggplot(data=thisPeptide[thisPeptide$Concentration >0,], aes(x=log(Concentration,10), y=log(MedianMeasuredC,10), color=FragmentIon)) + geom_errorbar(aes(ymin=log(MinMeasuredC,10) , ymax=log(MaxMeasuredC, 10)),position=pd, width=.08)+geom_point(position=pd, size=2) + xlab(mxlabel) + ylab("Log Measured Concentration (fmol/ug") + theme(title=element_text(size=18, colour="black"), axis.text=element_text(size=16), axis.text.x=element_text( colour=mcolor), axis.title=element_text(size=20), axis.title.x=element_text(colour=mcolor), legend.position = c(0.15, 0.85), legend.title = element_text(size=14), legend.text = element_text(size=14))+ labs(title=mTitle)+ scale_colour_discrete(name = "Transition")
                                print(p)
                                dev.off()
                            }
                        }
                    }

                    #Calculate the robust linear fitting model of concentration of analyte vs.  median of measured concentration of analyte.
                    {
                        uniqueT <- unique(thisPeptide$FragmentIon)
                        thisPeptide <- thisPeptide[thisPeptide$Concentration >0,]
                        thisPeptide <- thisPeptide[with(thisPeptide, order(FragmentIon, Concentration)),]
                        thisPeptideR <- thisPeptide[(thisPeptide$MedianMeasuredC != 0 & is.finite(thisPeptide$MedianMeasuredC)),c("FragmentIon", "Concentration", "MedianMeasuredC")]
                        fitR <- NULL
                        for (j in 1:length(uniqueT)){
                            thisFragmentIon <- thisPeptideR[(thisPeptideR$FragmentIon == uniqueT[j]),]
                            x <- thisFragmentIon$Concentration
                            y <- thisFragmentIon$MedianMeasuredC
                            w <- 1/(y)^2
                            tfit <- try(rlm(x~y, weights=w, method="MM", maxit=1000), silent=TRUE)
                            if (!inherits ( tfit, "try-error")){
                                mCoef <- coef(summary(tfit))
                                r2 <- (cor(x,y))^2
                                fitR <- rbind(fitR, c(uniqueT[j], mCoef[2,1], mCoef[1,1], mCoef[2,2], mCoef[1,2],r2))
                            }else{
                               #print("Regression fit failed")
                               invisible()
                            }
                        }
                        fitR <- as.data.frame(fitR)
                        names(fitR) <-  c("FragmentIon", "Slope", "Intercept", "SlopeStdErr", "InterceptStdErr", "RSquare")
                    }
                    for (i in 2: dim(fitR)[2]){
                        fitR[,i] <- as.numeric(as.character(fitR[,i]))
                    }
                    
                    # In fitR and LODTable, all the fragment ions will be kept.
                    # In fitR and LODTable, only keep the fragment ions in keptFragmentIon
                    #fitR <- fitR[fitR$FragmentIon %in% keptFragmentIon, ]
                    #LODTable <- LODTable[LODTable$FragmentIon %in% keptFragmentIon, ]
                    
                    # Write tables
                    write.table(format(fitR, digits=3), file = paste(plot_output_dir, "\\", input_peptide_sequence, '_', indexLabel, '_ResponseCurveAnalysis.fitTable.tsv', sep=''), sep = "\t", qmethod = "double", col.names=TRUE, row.names=FALSE, quote=FALSE)
                    #write.csv(format(fitR, digits=3), file=paste(plot_output_dir, "\\", input_peptide_sequence, '_ResponseCurveAnalysis.fitTable.csv', sep=''), row.names=FALSE, quote=FALSE)
                    
                    mergeTable <- merge(LODTable, fitR, by="FragmentIon", all=T)
                    mergeTable[,2:7] <-  (mergeTable[,2:7] - mergeTable[,9])  /(mergeTable[,8])
                    mergeTable <- mergeTable[,c(1,2,4,6,3,5,7)]
                    write.table(format(mergeTable, digits=3), file = paste(plot_output_dir, "\\", input_peptide_sequence, '_', indexLabel, '_ResponseCurveAnalysis.LODCTable.tsv', sep=''), sep = "\t", qmethod = "double", col.names=TRUE, row.names=FALSE, quote=FALSE)   
                    #write.csv(format(mergeTable[,1:7], digits=3), file=paste(plot_output_dir, "\\", input_peptide_sequence, '_ResponseCurveAnalysis.LODCTable.csv', sep=''), row.names=FALSE, quote=FALSE)
                }
            }
        }
    }
}
