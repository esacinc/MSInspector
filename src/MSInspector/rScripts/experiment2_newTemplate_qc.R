# This sample code is used to qc skyline document data from experiment 1 of Assay Portal.
# It's modified based on:
# esac-panorama-master\AssayPortal\resources\reports\schemas\targetedms\QCAnalysisQuery\web_portal_QC.r

# UPDATE: return top 3 plots plus sum of ions plot (Sept 2014)

# UPDATE: new sort to get top 3 plots (March 2016)

# UPDATE: in cases where there are medium labeled peptides instead of light labeled peptides, 
#   'medium' Isotope Label is converted to 'light' Isotope Label to calculate variable/constant ratios (April 2016)
# For a reverse curve: H/L = H/M. For a forward curve: L/H = M/H.

# UPDATE: convert all sample_group values to Lo, Med, Hi based on first character (L, M, H) of sample_group (April 2016)

#library(Cairo) # need for producing PNG image using Panorama
#suppressWarnings(suppressMessages(library(Rlabkey)))
suppressWarnings(suppressMessages(library(Cairo)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(evaluate)))
suppressWarnings(suppressMessages(require(reshape2)))
suppressWarnings(suppressMessages(require(dplyr, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(require(plyr)))

options(dplyr.summarise.inform = FALSE)

# ***** plot_QC function *****
plot_QC <- function(plot_fragment_ion_results, input_peptide_sequence, current_ion, days) {
    if (current_ion == 'sum'){
        current_ion <- 'sum of ions'
    }

    plot_title <- paste(input_peptide_sequence, current_ion, sep='\n')
    min_calculated_area_ratio <- min(plot_fragment_ion_results$calculated_area_ratio)
    max_calculated_area_ratio <- max(plot_fragment_ion_results$calculated_area_ratio)

    # Expand right side of clipping rect to make room for the legend
    par(xpd=TRUE, mar=par()$mar+c(0,0,0,4))
    # bty="L",
    
    suppressWarnings(plot(plot_fragment_ion_results$day, plot_fragment_ion_results$calculated_area_ratio, log="y", yaxt="n", col=ifelse(plot_fragment_ion_results$sample_group=='Hi', 'red',
                                                                                                                        ifelse(plot_fragment_ion_results$sample_group=='Med', 'blue', 'green')),
    
        pch=ifelse(plot_fragment_ion_results$replicate==1 , 1,
                    ifelse(plot_fragment_ion_results$replicate==2 , 0,
                            ifelse(plot_fragment_ion_results$replicate==3 , 2,
                                    ifelse(plot_fragment_ion_results$replicate==4 , 5,
                                        ifelse(plot_fragment_ion_results$replicate==5 , 6,
                                    4))))),
    
        cex=2, lwd=1, main=plot_title, xlab="Time (day)", ylab="Measured (area ratio) [log-scale]", cex.lab=1.75, cex.axis=1.75))

    #x_axis_values <- c(1,2,3,4,5)
    x_axis_values <- days
    #y_axis_values <- c(0.01,0.1,1,10,100,format(max_calculated_measured_concentration,digits=3))
    y_axis_values <- c(format(min_calculated_area_ratio,digits=3),format(median(plot_fragment_ion_results$calculated_area_ratio),digits=3),format(max_calculated_area_ratio,digits=3))
    #format(y_axis_values,scientific=FALSE,digits=4)
    axis(2, y_axis_values, labels=format(y_axis_values,scientific=FALSE)) # draw y axis with required labels
    #par(xpd=TRUE)
    #legend(x=4,y=max_calculated_area_ratio,legend=c("rep1","rep2","rep3","repX","Hi","Med","Lo"),pch=c(1,0,2,4,18,18,18), col=c("black","black","black","black","red","blue","green"), bty="n")
    legend(x=max(days)+0.2,y=max_calculated_area_ratio,legend=c("rep1","rep2","rep3","rep4","rep5","Hi","Med","Lo"),pch=c(1,0,2,5,6,18,18,18), col=c("black","black","black","black","black","red","blue","green"), cex=1, bty="n")
    # Restore default clipping rect
    par(mar=c(5, 4, 4, 2) + 0.1)
} # end plot_QC function

identify_uniProtKB_entryID  <- function(x) {
    # This function is to extract uniProtKB_entryID from the protein name.
    if (grepl('\\|', x)) {
        tmp <- strsplit(x, split = '\\|')
        uniProtKB_entryID_tmp <- tmp[[1]][2]
        # judge whether uniProtKB_entryID is a legal uniProtKB entry ID based on its pattern using regular expression. Please refer to https://www.uniprot.org/help/accession_numbers
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

check_completeness <- function(x) {
    # it will return the state of completeness. If it's TRUE, it means that it's complete. Otherwise, it means that it's wrongly annotated.
    x_tmp <- sort(unique(tolower(x)))
    completeness_state <- TRUE
    for (x_item in x_tmp) {
        if (!(x_item %in% c('h', 'm', 'l', 'hi', 'med', 'lo', 'high', 'medium', 'low'))) {
            completeness_state <- FALSE
            break
        }
    }
    completeness_state
}


args <- commandArgs(trailingOnly = TRUE)
dataset_path <- args[1]
fileList_path <- args[2]
plot_output <- args[3]
plot_output_dir <- args[4]


#dataset_path <- "normal_data.tsv"
#fileList_path <- "file_namelist_IS.tsv"
#plot_output <- "True"
#plot_output_dir <- "D:\\Skyline_analysis\\qcAssayPortal\\qcAssayPortal\\src\\qcAssayPortal\\rScripts\\test\\exp2_new_error\\tmp"



if (plot_output == 'True') {
    plot_output <- TRUE
} else {
    plot_output <- FALSE
}

cv_threshold_all <- 60.0
cv_threshold_individual <- 40.0
cv_threshold <- 20.0
allInjectionSwitch <- TRUE
pValue_threshold <- 0.05
relative_difference_threshold <- 0.2

# Load data from local table
QC_set_total <- read.table(file=dataset_path, header=TRUE, sep='\t')
fileDf <- read.table(file=fileList_path, header=TRUE, sep='\t')

colNumber <- ncol(QC_set_total)
thenames <- tolower(names(QC_set_total))
names(QC_set_total) <- thenames

# Rename columns in QC_set_total dataframe (replace Panorama names with new names used by R script)
colnames(QC_set_total)[colnames(QC_set_total)=="skydocumentname"] <- "SkyDocumentName"
colnames(QC_set_total)[colnames(QC_set_total)=="proteinname"] <- "protein_name"
colnames(QC_set_total)[colnames(QC_set_total)=="peptidemodifiedsequence"] <- "peptide"
colnames(QC_set_total)[colnames(QC_set_total)=="isotopelabeltype"] <- "isotope_label_type"     # light, heavy
colnames(QC_set_total)[colnames(QC_set_total)=="precursorcharge"] <- "precursor_charge"
colnames(QC_set_total)[colnames(QC_set_total)=="productcharge"] <- "product_charge"
colnames(QC_set_total)[colnames(QC_set_total)=="fragmention"] <- "fragment_ion_only"
colnames(QC_set_total)[colnames(QC_set_total)=="area"] <- "area"
colnames(QC_set_total)[colnames(QC_set_total)=="replicatename"] <- "replicate_name"
colnames(QC_set_total)[colnames(QC_set_total)=="replicatenumber"] <- "replicate"
colnames(QC_set_total)[colnames(QC_set_total)=="day"] <- "day"      # day
colnames(QC_set_total)[colnames(QC_set_total)=="exp2samplegroup"] <- "sample_group"     # Lo, Med, Hi

if (nrow(QC_set_total) ==0) {
  QC_set_total$fragment_ion <- integer(0)
} else {
  QC_set_total$fragment_ion <- paste(QC_set_total[ ,'fragment_ion_only'], " (", QC_set_total[ ,'product_charge'], "+)", sep='' )
}
  
# Preprocess the columns in QC_set_total
# convert columns from character to numeric
QC_set_total[,'day'] <- as.numeric(as.character(QC_set_total[,'day']))
QC_set_total[,'replicate'] <- as.numeric(as.character(QC_set_total[,'replicate']))
QC_set_total[,'area'] <- as.numeric(as.character(QC_set_total[,'area']))

# remove factor version
QC_set_total[,'SkyDocumentName'] <- as.character(QC_set_total[,'SkyDocumentName'])
QC_set_total[,'protein_name'] <- as.character(QC_set_total[,'protein_name'])
QC_set_total[,'peptide'] <- as.character(QC_set_total[,'peptide'])
QC_set_total[,'precursor_charge'] <- as.character(QC_set_total[,'precursor_charge'])
QC_set_total[,'product_charge'] <- as.character(QC_set_total[,'product_charge'])
QC_set_total[,'isotope_label_type'] <- as.character(QC_set_total[,'isotope_label_type'])
QC_set_total[,'fragment_ion_only'] <- as.character(QC_set_total[,'fragment_ion_only'])
QC_set_total[,'replicate_name'] <- as.character(QC_set_total[,'replicate_name'])
QC_set_total[,'fragment_ion'] <- as.character(QC_set_total[,'fragment_ion'])
QC_set_total[,'sample_group'] <- as.character(QC_set_total[,'sample_group'])
#QC_set_total[,'sample_group'] <- as.character(QC_set_total[,'sample_group'])

# Convert all sample_group values to Lo, Med, Hi based on first character of sample_group
if (nrow(QC_set_total) > 0) {
  QC_set_total[toupper(substr(QC_set_total$sample_group, 1, 1))=="L","sample_group"] <- "Lo"
  QC_set_total[toupper(substr(QC_set_total$sample_group, 1, 1))=="M","sample_group"] <- "Med"
  QC_set_total[toupper(substr(QC_set_total$sample_group, 1, 1))=="H","sample_group"] <- "Hi"
}


# Initialize an empty data.frame QC_set_filtered with the column names same as QC_set_total's.
#QC_set_filtered = data.frame()
#for (k in colnames(QC_set_total)) {
#    QC_set_filtered[[k]]<-as.character()
#}

# Write peptide information into output file.
log_filename <- paste(plot_output_dir, "\\peptide_infor.tsv", sep='' )
logdf <- data.frame(peptide=as.character(), precursorCharge=as.character(), isotopeLabelType=as.character(), transition=as.character(), uniProtKBID=as.character(), proteinName=as.character(), SkyDocumentName=as.character())

# Separate the error detecting codes from the warning detecting codes.
# Traverse the SkyDocumentName in fileDf to detect all the possible errors.
# Create a list to store the  peptides with errors for each SkyDocumentName.
# df_skydoc_error_peptide stores the peptides with errors. But there my be one special situation that one specific precursor charge of the peptide has errors but the others don't have error.
# Therefore, df_skydoc_error_peptide_precursorCharge is used to store peptide plus precursor charge.
df_skydoc_error_peptide <- list()
df_skydoc_error_peptide_precursorCharge <- data.frame(SkyDocumentName=as.character(), protein_name=as.character(), peptide=as.character(), precursorCharge=as.character())
for (SkyDocumentName in as.character(fileDf[, "SkyDocumentName"])) {
    QC_set1 <- QC_set_total[QC_set_total$SkyDocumentName==SkyDocumentName, ]
    # Get a list of all peptides
    peptide_list <- unique(QC_set1[ , 'peptide'])
    peptide_list_with_error <- c()
    for (input_peptide_sequence in peptide_list) {
        QC_set2 <- QC_set1[QC_set1$peptide==input_peptide_sequence, ]
        # Get a list of all precursor_charges
        precursor_charge_list <- unique(QC_set2[ , 'precursor_charge'])
        for (input_precursor_charge in precursor_charge_list) {
            QC_set3 <- QC_set2[QC_set2$precursor_charge==input_precursor_charge, ]
            # Get a list of all protein names, although usually one peptide with a specific precursor charge has only one protein.
            protein_list <- unique(QC_set3[ , 'protein_name'])
            protein_uniProtID_list <- sapply(protein_list, identify_uniProtKB_entryID, USE.NAMES=FALSE)
            for (indexLabel in 1:length(protein_list)) {
                input_protein_name <- protein_list[indexLabel]
                protein_uniProtID <- protein_uniProtID_list[indexLabel]
                # Choose the specific peptide with a specific precursor charge from a specific protein.
                QC_setTmp <- QC_set3[QC_set3$protein_name==input_protein_name, ]
                if (nrow(QC_setTmp) >= 1) {
                    isotopeLabelTypeTmp <- paste(sort(unique(QC_setTmp$isotope_label_type)), collapse = '|')
                    transitionTmp <- c()
                    for (isotopeLabelTypeSubtmp in sort(unique(QC_setTmp$isotope_label_type))) {
                        transitionTmp <- c(transitionTmp, paste(isotopeLabelTypeSubtmp, paste(sort(unique(QC_setTmp[QC_setTmp$isotope_label_type==isotopeLabelTypeSubtmp, ]$fragment_ion)), collapse = '|'), sep=':'))
                    }
                    transitionTmp <- paste(transitionTmp, collapse = ';')
                    logdfTmp <- data.frame(peptide=input_peptide_sequence, precursorCharge=input_precursor_charge, isotopeLabelType=isotopeLabelTypeTmp, transition=transitionTmp, uniProtKBID=protein_uniProtID, proteinName=input_protein_name, SkyDocumentName=SkyDocumentName)
                    logdf <- rbind(logdf, logdfTmp)
                }
                # Get a list of all unique fragment ions, unique days, unique sample groups (Lo, Med, Hi), unique replicates and unique isotope_label_type associated with current peptide
                fragment_ion_list <- unique(QC_setTmp[ , 'fragment_ion'])
                days <- sort(unique(QC_setTmp[ , 'day']))
                sample_groups <- sort(unique(QC_setTmp[ , 'sample_group']))
                replicates <- sort(unique(QC_setTmp[ , 'replicate']))
                isotope_label_types <- unique(QC_setTmp[ , 'isotope_label_type'])
                # Check the controlled vocabulary for the column sample_group
                sample_groups_for_check <- sort(unique(QC_setTmp[ , 'sample_group']))
                completeness_state <- check_completeness(sample_groups_for_check)
                if (! completeness_state) {
                    errorType <- "Error"
                    errorSubtype <- "Attribute"
                    errorReason <- paste('The controlled vocabulary of the attribute of sample group must be High, Medium, Low, while it has other annotated terms: ', paste(sample_groups_for_check, collapse = ', '), '.', sep='')
                    #errorInfor <- paste(SkyDocumentName, errorType, errorSubtype, errorReason, '', '', '', '', '', '', '', '', '', '', '', '', sep='\t')
                    #errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason), rep('', colNumber-1)), collapse='\t')
                    errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge), rep('', colNumber-5)), collapse='\t')
                    cat(errorInfor)
                    cat('\n')
                    next
                }

                # *** for medium labeled peptides ***
                if(('light' %in% isotope_label_types) & ('medium' %in% isotope_label_types)) {
                    errorType <- "Error"
                    errorSubtype <- "Light and Medium isotope"
                    errorReason <- "Both light and medium isotope labels found in the peptide with a specific charge."
                    #errorInfor <- paste(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge, '', '', '', '', '', '', '', '',sep='\t')
                    errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge), rep('', colNumber-5)), collapse='\t')
                    cat(errorInfor)
                    cat('\n')
                    peptide_list_with_error <- c(peptide_list_with_error, input_peptide_sequence)
                    df_skydoc_error_peptide_precursorChargeTmp <- data.frame(SkyDocumentName=SkyDocumentName, protein_name=input_protein_name, peptide=input_peptide_sequence, precursorCharge=input_precursor_charge)
                    df_skydoc_error_peptide_precursorCharge <- rbind(df_skydoc_error_peptide_precursorCharge, df_skydoc_error_peptide_precursorChargeTmp)
                    next
                } else {
                    QC_setTmp$isotope_label_type[QC_setTmp$isotope_label_type == "medium"] <- "light"
                }

                # Make judgement whether there are multiple heavy or light area for the combination of fragment_ion, replicate, day and sample_group.
                # If it happens, traverse fragment_ion_list, days, sample_groups and replicates to evaluate the fragment_ion under the specific combination of day, sample_group, and replicate.
                # The reason to this error is that the annotation of column 
                evaOut1 <- evaluate("dcast(QC_setTmp, protein_name + peptide + precursor_charge + fragment_ion + replicate + day + sample_group ~ isotope_label_type, value.var='area')")
                evaOut2 <- evaluate("dcast(QC_setTmp, protein_name + peptide + precursor_charge + fragment_ion + replicate_name + day + sample_group ~ isotope_label_type, value.var='area')")
                
                if (length(evaOut1) == 3) {
                    # In this condition, some replicate information is wrong for some combinations of fragment_ion, day and sample_group.
                    # The wrongly annotated replicate need to be generated.
                    df1 <- suppressMessages(dcast(QC_setTmp, protein_name + peptide + precursor_charge + fragment_ion + replicate + day + sample_group ~ isotope_label_type, value.var='area'))
                    # Evaluate the fragment_ion under the specific combination of day, sample_group, and replicate.
                    # df1 can be used to extract the combinations
                    errorReasonTmp <- c()
                    for (index in 1:nrow(df1)) {
                        current_set <- QC_setTmp[QC_setTmp$fragment_ion==df1[index, ]$fragment_ion & QC_setTmp$day==df1[index, ]$day & QC_setTmp$sample_group==df1[index, ]$sample_group & QC_setTmp$replicate==df1[index, ]$replicate, ]
                        lightCount <- nrow(current_set[current_set$isotope_label_type=='light', ])
                        heavyCount <- nrow(current_set[current_set$isotope_label_type=='heavy', ])
                        
                        if (lightCount != 1 | heavyCount != 1) {
                            fragment_ion_error <- current_set$fragment_ion[1]
                            if (lightCount != 1) {
                                light_error <- paste(lightCount, ' light isotopes')
                            } else {
                                light_error <- ''
                            }
                            if (heavyCount != 1) {
                                heavy_error <- paste(heavyCount, ' heavy isotopes')
                            } else {
                                heavy_error <- ''
                            }
                            errorReason_item1 <- paste(paste('For ', fragment_ion_error, ': ', sep=''), paste(heavy_error, light_error, sep=' '), sep='')
                            if (length(evaOut2) == 2) {
                                errorReason_item2 <- paste(unique(current_set$replicate_name), collapse = ' | ')
                                errorReason_item <- paste(errorReason_item1, ' due to multiple values in attributes: replicate_name (', errorReason_item2, ')', sep='')
                            } else {
                                errorReason_item <- paste(errorReason_item1, ' due to wrongly annotated values in attributes: replicate, day or sample group', sep='')
                            }
                            if (!(errorReason_item %in% errorReasonTmp)) {
                                errorReasonTmp <- c(errorReasonTmp, errorReason_item)
                            }
                        } else {
                            invisible() 
                        }
                    }
                    if (length(errorReasonTmp) > 0) {
                        errorType <- "Error"
                        errorSubtype <- "Area values of heavy or light Isotope"
                        errorReason <- paste(paste(errorReasonTmp, collapse='. '), '.', sep='')
                        #errorInfor <- paste(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge, '', '', '', '', '', '', '', '',sep='\t')
                        errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge), rep('', colNumber-5)), collapse='\t')
                        cat(errorInfor)
                        cat('\n')
                        peptide_list_with_error <- c(peptide_list_with_error, input_peptide_sequence)
                        df_skydoc_error_peptide_precursorChargeTmp <- data.frame(SkyDocumentName=SkyDocumentName, protein_name=input_protein_name, peptide=input_peptide_sequence, precursorCharge=input_precursor_charge)
                        df_skydoc_error_peptide_precursorCharge <- rbind(df_skydoc_error_peptide_precursorCharge, df_skydoc_error_peptide_precursorChargeTmp)
                        next
                    }
                } else {
                    invisible()
                }
            }
        }
    }
    peptide_list_with_error <- unique(peptide_list_with_error)
    df_skydoc_error_peptide[[SkyDocumentName]] <- peptide_list_with_error
}

if (plot_output) {
    write.table(logdf, file=log_filename, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
}

df_skydoc_error_peptide_precursorCharge$SkyDocumentName <-  as.character(df_skydoc_error_peptide_precursorCharge$SkyDocumentName)
df_skydoc_error_peptide_precursorCharge$protein_name <-  as.character(df_skydoc_error_peptide_precursorCharge$protein_name)
df_skydoc_error_peptide_precursorCharge$peptide <-  as.character(df_skydoc_error_peptide_precursorCharge$peptide)
df_skydoc_error_peptide_precursorCharge$precursorCharge <- as.character(df_skydoc_error_peptide_precursorCharge$precursorCharge)

# Infer the internal standard type for each SkyDocumentName by randomly sampled 5 peptides.
if (FALSE) {
    df_internal_standard_inferred <- data.frame(SkyDocumentName=as.character(), internal_standard=as.character())
    for (SkyDocumentName in as.character(fileDf[, "SkyDocumentName"])) {
        df_internal_standard_inferred_tmp <- data.frame(SkyDocumentName=SkyDocumentName, internal_standard='light')
        df_internal_standard_inferred <- rbind(df_internal_standard_inferred, df_internal_standard_inferred_tmp)
    }
}

# Write peptide information into output file.
is_inferred_filename <- paste(plot_output_dir, "\\internal_standard_inferred_infor.tsv", sep='' )
df_internal_standard_inferred <- data.frame(SkyDocumentName=as.character(), internal_standard=as.character())
for (SkyDocumentName in as.character(fileDf[, "SkyDocumentName"])) {
    QC_set_1 <- QC_set_total[QC_set_total$SkyDocumentName==SkyDocumentName, ]
    peptide_list <- unique(QC_set_1[ , 'peptide'])
    
    # Remove the peptides with errors.
    if (SkyDocumentName %in% names(df_skydoc_error_peptide)) {
        peptide_list <- setdiff(peptide_list,df_skydoc_error_peptide[[SkyDocumentName]])
    }
    
    # Randomly sampling 5 peptides. If the number of the peptides is less than 5, keep all of them.
    if (length(peptide_list) >= 5) {
        peptide_list_sampled <- sample(peptide_list, 5, replace = FALSE)
    } else {
        peptide_list_sampled <- peptide_list
    }
    
    if (length(peptide_list_sampled) == 0) {
        internal_standard_inferred <- "can't be inferred"
    } else {
        # If the internal_standard_hypothetical is light, the curve_type will be reverse.
        # If the internal_standard_hypothetical is heavy, the curve_type will be forward.
        internal_standard_hypothetical <- 'light'
        curve_type_hypothetical <- 'reverse'
        value1 <- 0
        value2 <- 0
        
        for (input_peptide_sequence in peptide_list_sampled) {
            QC_set_2 <- QC_set_1[QC_set_1$peptide==input_peptide_sequence, ]
            precursor_charge_list <- unique(QC_set_2[ , 'precursor_charge'])
            for (input_precursor_charge in precursor_charge_list) {
                QC_set_3 <- QC_set_2[QC_set_2$precursor_charge==input_precursor_charge, ]
                protein_list <- unique(QC_set_3[ , 'protein_name'])
                protein_uniProtID_list <- sapply(protein_list, identify_uniProtKB_entryID, USE.NAMES=FALSE)
                for (indexLabel in 1:length(protein_list)) {
                    input_protein_name <- protein_list[indexLabel]
                    protein_uniProtID <- protein_uniProtID_list[indexLabel]
                    QC_set <- QC_set_3[QC_set_3$protein_name==input_protein_name, ]
                    
                    fragment_ion_list <- unique(QC_set[ , 'fragment_ion'])
                    days <- sort(unique(QC_set[ , 'day']))
                    sample_groups <- sort(unique(QC_set[ , 'sample_group']))
                    replicates <- sort(unique(QC_set[ , 'replicate']))
                    isotope_label_types <- unique(QC_set[ , 'isotope_label_type'])
                    QC_set$isotope_label_type[QC_set$isotope_label_type == "medium"] <- "light"
                    
                    fragment_ion_results <- data.frame()
                    for (current_ion in fragment_ion_list) {
                        for (current_day in days) {
                            for (current_sample in sample_groups) {
                                for (current_rep in replicates) {
                                    current_set <- QC_set[QC_set$fragment_ion==current_ion & QC_set$day==current_day & QC_set$sample_group==current_sample & QC_set$replicate==current_rep, ]
                                    current_set_count <- nrow(current_set)
                                    if (current_set_count == 2) {
                                        light_area <- current_set[current_set$isotope_label_type=='light', 'area' ]
                                        heavy_area <- current_set[current_set$isotope_label_type=='heavy', 'area' ]
                                        theoretical_area <- light_area
                                        measured_area <- heavy_area

                                        if(theoretical_area==0 | is.na(theoretical_area) | is.na(measured_area)){
                                            calculated_area_ratio = NA
                                        }
                                        else {
                                            calculated_area_ratio <- measured_area/theoretical_area
                                        }
                                        fragment_ion_results <-  rbind(fragment_ion_results, data.frame(Peptide = input_peptide_sequence, Protein_Name = input_protein_name, Precursor_Charge = input_precursor_charge,
                                                                                                      fragment_ion = current_ion, day = current_day, sample_group = current_sample, replicate = current_rep, light_area=light_area, heavy_area=heavy_area,
                                                                                                      theoretical_area=theoretical_area, measured_area=measured_area,
                                                                                                      calculated_area_ratio=calculated_area_ratio, ion_category='individual') )
                                    }
                                }
                            }
                        }
                    }
                    
                    # Selected the top 3 fragment ions which are ranked by med_total_CV.
                    # Step 1:
                    # ***** calculate CV (Coefficient of Variation) *****
                    LoMedHi <- c('Lo', 'Med', 'Hi')
                    ions <- fragment_ion_list
                    # make CV summary data frame
                    CV_results <- data.frame()
                    for (current_ion in ions) {
                      # define columns (make row for each fragment ion)
                        CV_results <-  rbind(CV_results, data.frame(fragment_ion=current_ion, low_intra_CV=NA, med_intra_CV=NA, high_intra_CV=NA,
                                                                  low_inter_CV=NA, med_inter_CV=NA, high_inter_CV=NA,
                                                                  low_total_CV=NA, med_total_CV=NA, high_total_CV=NA,
                                                                  low_count=NA, med_count=NA, high_count=NA))
                    }
                    
                    # Step 2:
                    # *** intra-assay CV ***
                    for (current_ion in ions) {
                        for (current_LoMedHi in LoMedHi) {
                            avg_intra_assay_CV <-0
                            individual_intra_assay_CVs <- c()
                            for (current_day in days) {
                                #for (current_rep in replicates) {
                                current_set <- fragment_ion_results[fragment_ion_results$fragment_ion==current_ion & fragment_ion_results$day==current_day & fragment_ion_results$sample_group==current_LoMedHi, ]
                                
                                # remove rows with a value of NA for calculated_area_ratio
                                current_set <- current_set[complete.cases(current_set[ , 'calculated_area_ratio' ]),]
                                
                                if (nrow(current_set) <= 1 ){
                                    percent_CV <- NA
                                }
                                else{
                                    percent_CV <- (sd(current_set$calculated_area_ratio))/(mean(current_set$calculated_area_ratio)) * 100
                                    individual_intra_assay_CVs <- c(individual_intra_assay_CVs, percent_CV)
                                }
                                #  } # end current_rep
                            } # end current_day

                            if (length(individual_intra_assay_CVs)==0){
                                avg_CV <- NA
                                count <- 0
                            }
                            else{
                                avg_CV <- mean(individual_intra_assay_CVs, na.rm = TRUE)
                            }

                            if (current_LoMedHi=='Lo'){
                                CV_results[CV_results$fragment_ion==current_ion, 'low_intra_CV'] <- round(avg_CV, digits=1)
                            }
                            else if(current_LoMedHi=='Med'){
                                CV_results[CV_results$fragment_ion==current_ion, 'med_intra_CV'] <- round(avg_CV, digits=1)
                            }
                            else if(current_LoMedHi=='Hi'){
                                CV_results[CV_results$fragment_ion==current_ion, 'high_intra_CV'] <- round(avg_CV, digits=1)
                            }
                        }
                    }
                    
                    # Step 3:
                    # *** inter-assay CV ***
                    for (current_ion in ions) {
                        for (current_LoMedHi in LoMedHi) {
                            avg_inter_assay_CV <-0
                            individual_inter_assay_CVs <- c()
                            
                            avg_CV  <- 0
                            # for (current_day in days) {
                            for (current_rep in replicates) {           
                                current_set <- fragment_ion_results[fragment_ion_results$fragment_ion==current_ion & fragment_ion_results$replicate==current_rep & fragment_ion_results$sample_group==current_LoMedHi, ]
                                # remove rows with a value of NA for calculated_area_ratio
                                current_set <- current_set[complete.cases(current_set[ , 'calculated_area_ratio' ]),]
                                if (nrow(current_set) <= 1 ){
                                  percent_CV <- NA
                                }
                                else{
                                  percent_CV <- (sd(current_set$calculated_area_ratio))/(mean(current_set$calculated_area_ratio)) * 100
                                  individual_inter_assay_CVs <- c(individual_inter_assay_CVs, percent_CV)
                                }
                            } # end current_rep
                            # } # end current_day
                            
                            if (length(individual_inter_assay_CVs)==0){
                                avg_CV <- NA
                                count <- 0
                            }
                            else{
                                avg_CV <- mean(individual_inter_assay_CVs, na.rm = TRUE)
                            }
                            if (current_LoMedHi=='Lo'){
                                CV_results[CV_results$fragment_ion==current_ion, 'low_inter_CV'] <- round(avg_CV, digits=1)
                            }
                            else if(current_LoMedHi=='Med'){
                                CV_results[CV_results$fragment_ion==current_ion, 'med_inter_CV'] <- round(avg_CV, digits=1)
                            }
                            else if(current_LoMedHi=='Hi'){
                                CV_results[CV_results$fragment_ion==current_ion, 'high_inter_CV'] <- round(avg_CV, digits=1)
                            }
                        }
                    }
                    
                    # Calculate total variability
                    CV_results[ , 'med_total_CV'] <- round(sqrt((CV_results[ , 'med_intra_CV'])*(CV_results[ , 'med_intra_CV']) + (CV_results[ , 'med_inter_CV'])*(CV_results[ , 'med_inter_CV'])), digits=1)
                    CV_results[ , 'high_total_CV'] <- round(sqrt((CV_results[ , 'high_intra_CV'])*(CV_results[ , 'high_intra_CV']) + (CV_results[ , 'high_inter_CV'])*(CV_results[ , 'high_inter_CV'])), digits=1)
                    CV_results[ , 'low_total_CV'] <- round(sqrt((CV_results[ , 'low_intra_CV'])*(CV_results[ , 'low_intra_CV']) + (CV_results[ , 'low_inter_CV'])*(CV_results[ , 'low_inter_CV'])), digits=1)
                    
                    # Determine counts
                    for (current_ion in ions) {
                        CV_results[CV_results$fragment_ion==current_ion, 'high_count'] <- nrow(fragment_ion_results[fragment_ion_results$sample_group=='Hi' & fragment_ion_results$fragment_ion==current_ion & !is.na(fragment_ion_results$calculated_area_ratio), ] )
                        CV_results[CV_results$fragment_ion==current_ion, 'med_count'] <- nrow(fragment_ion_results[fragment_ion_results$sample_group=='Med' & fragment_ion_results$fragment_ion==current_ion & !is.na(fragment_ion_results$calculated_area_ratio), ] )
                        CV_results[CV_results$fragment_ion==current_ion, 'low_count'] <- nrow(fragment_ion_results[fragment_ion_results$sample_group=='Lo' & fragment_ion_results$fragment_ion==current_ion & !is.na(fragment_ion_results$calculated_area_ratio), ] )
                    }

                    ions_to_plot <- c()
                    # Determine fragment ions to plot
                    if (length(ions) <= 3 ) {
                        ions_to_plot <- ions
                    } else {
                        results_to_plot <- CV_results[!is.na(CV_results$low_total_CV) & !is.na(CV_results$med_total_CV) & !is.na(CV_results$high_total_CV) , ]
                        # new sort to get Top 3 plots
                        results_to_plot <- results_to_plot[order(results_to_plot$med_total_CV, results_to_plot$low_total_CV, results_to_plot$high_total_CV), ]
                        three_lowest_med_total_CV <- head(results_to_plot, 3)
                        ions_to_plot <- as.character(three_lowest_med_total_CV[ , 'fragment_ion'])
                    }
                    
                    #print(input_protein_name)
                    #print(ions_to_plot)
                    # Transverse the ion in ions_to_plot. For day of Lo and Hi in each day, compare the median value of calculated_area_ratio.
                    for (ion_tmp in ions_to_plot) {
                        plot_fragment_ion_results <- fragment_ion_results[!is.na(fragment_ion_results$calculated_area_ratio), ]
                        days_list <- sort(unique(plot_fragment_ion_results$day))
                        for (day_tmp in days_list) {
                            concentration_sort <- unique(plot_fragment_ion_results$sample_group)
                            # Sort concentration_sort by the sequence of Hi, Med and Lo
                            hi_index <- which(concentration_sort=='Hi')
                            med_index <- which(concentration_sort=='Med')
                            lo_index <- which(concentration_sort=='Lo')
                            valid_index <- c()
                            for (item in c(hi_index, med_index, lo_index)) {
                                if (length(item) == 1) {
                                    valid_index <- c(valid_index, item)
                                }
                            }
                            if (length(valid_index) >= 2) {
                                hi_index_new <- concentration_sort[valid_index][1]
                                lo_index_new <- concentration_sort[valid_index][length(concentration_sort[valid_index])]
                                medHLRation_hi <- median(plot_fragment_ion_results[plot_fragment_ion_results$day == day_tmp & plot_fragment_ion_results$sample_group == hi_index_new, ]$calculated_area_ratio, na.rm=TRUE)
                                medHLRation_lo <- median(plot_fragment_ion_results[plot_fragment_ion_results$day == day_tmp & plot_fragment_ion_results$sample_group == lo_index_new, ]$calculated_area_ratio, na.rm=TRUE)
                                if (is.na(medHLRation_hi) | is.na(medHLRation_lo)) {
                                    invisible()
                                } else if (medHLRation_hi >= medHLRation_lo) {
                                    value1 <- value1 + 1
                                } else {
                                    value2 <- value2 + 1
                                }
                            }
                        }
                    
                    }
                }
            }
        }
        if (value1==0 & value2==0) {
            internal_standard_inferred <- "can't be inferred"
        } else if (value1 >= value2) {
            internal_standard_inferred <- 'light'
        } else {
            internal_standard_inferred <- 'heavy'
        }
        #print(value1)
        #print('$$$$$$$$')
        #print(value2)
    }
    df_internal_standard_inferred_tmp <- data.frame(SkyDocumentName=SkyDocumentName, internal_standard=internal_standard_inferred)
    df_internal_standard_inferred <- rbind(df_internal_standard_inferred, df_internal_standard_inferred_tmp)
}
if (plot_output) {
    write.table(df_internal_standard_inferred, file=is_inferred_filename, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
}

    
#cat(dim(QC_set_filtered))
# Traverse the SkyDocumentName in fileDf to detect all the possible warnings.
for (SkyDocumentName in as.character(fileDf[, "SkyDocumentName"])) {
    # Evaluate the internal_standard, if the internal standard is wrong, errors will arise.
    original_internal_standard <- as.character(fileDf[fileDf$SkyDocumentName == SkyDocumentName, ]$internal_standard)
    inferred_internal_standard <- as.character(df_internal_standard_inferred[df_internal_standard_inferred$SkyDocumentName == SkyDocumentName, ]$internal_standard)
    if (inferred_internal_standard[1] == "can't be inferred") {
        next
    }
    if (original_internal_standard[1] == 'none') {
        # Just jump out of the loop. Don't print the errorInfor, because it has already be printed in the function of detectIS in qcAnalysis.py
        next
    }

    if (original_internal_standard[1] != inferred_internal_standard[1]) {
        errorType <- "Error"
        errorSubtype <- "Internal standard"
        errorReason <- paste('The internal standard in the skyline file is set to be ', original_internal_standard, ', while the inferred internal standard is ', inferred_internal_standard, '.', sep='')
        #errorInfor <- paste(SkyDocumentName, errorType, errorSubtype, errorReason, '', '', '', '', '', '', '', '', '', '', '', '', sep='\t')
        errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason), rep('', colNumber-1)), collapse='\t')
        cat(errorInfor)
        cat('\n')
        next
    }

    QC_set_1 <- QC_set_total[QC_set_total$SkyDocumentName==SkyDocumentName, ]
    peptide_list <- unique(QC_set_1[ , 'peptide'])
    
    internal_standard <- inferred_internal_standard[1]
    if (internal_standard == 'light') {
        curve_type <- 'reverse'
    } else if (internal_standard == 'heavy'){
        curve_type <- 'forward'
    } else {
        next
    }
    
    # Remove the peptides plus precursorCharge with errors in df_skydoc_error_peptide_precursorCharge
    #if (SkyDocumentName %in% names(df_skydoc_error_peptide)) {
    #    peptide_list <- setdiff(peptide_list,df_skydoc_error_peptide[[SkyDocumentName]])
    #}
    for (input_peptide_sequence in peptide_list) {
        QC_set_2 <- QC_set_1[QC_set_1$peptide==input_peptide_sequence, ]
        precursor_charge_list <- unique(QC_set_2[ , 'precursor_charge'])
        for (input_precursor_charge in precursor_charge_list) {
            QC_set_3 <- QC_set_2[QC_set_2$precursor_charge==input_precursor_charge, ]
            protein_list <- unique(QC_set_3[ , 'protein_name'])
            protein_uniProtID_list <- sapply(protein_list, identify_uniProtKB_entryID, USE.NAMES=FALSE)
            for (indexLabel in 1:length(protein_list)) {
                input_protein_name <- protein_list[indexLabel]
                protein_uniProtID <- protein_uniProtID_list[indexLabel]
                # Judge whether SkyDocumentName, input_protein_name, input_peptide_sequence and input_precursor_charge exist in df_skydoc_error_peptide_precursorCharge
                if ( nrow(subset(df_skydoc_error_peptide_precursorCharge, SkyDocumentName==SkyDocumentName & protein_name==input_protein_name & peptide==input_peptide_sequence & precursorCharge==input_precursor_charge)) > 0 ) {
                    # This means that SkyDocumentName, input_protein_name, input_peptide_sequence and input_precursor_charge exist in df_skydoc_error_peptide_precursorCharge
                    next
                }
                
                QC_set <- QC_set_3[QC_set_3$protein_name==input_protein_name, ]

                fragment_ion_list <- unique(QC_set[ , 'fragment_ion'])
                days <- sort(unique(QC_set[ , 'day']))
                sample_groups <- sort(unique(QC_set[ , 'sample_group']))
                replicates <- sort(unique(QC_set[ , 'replicate']))
                isotope_label_types <- unique(QC_set[ , 'isotope_label_type'])
                QC_set$isotope_label_type[QC_set$isotope_label_type == "medium"] <- "light"
                
                # capture the warning caused by the number of fragment ions 
                if (length(fragment_ion_list) < 3) {
                    errorType <- "Warning"
                    errorSubtype <- "Fragment ion"
                    errorReason <- paste("In repeatability graph, the number of fragment ions is ",  length(fragment_ion_list), " < 3, the fragment ions is(are): ", paste(fragment_ion_list, collapse = ', '), '.', sep="")
                    #errorInfor <- paste(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge, '', '', '', '', '', '', '', '', sep='\t')
                    errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge), rep('', colNumber-5)), collapse='\t')
                    cat(errorInfor)
                    cat('\n')
                }
                
                fragment_ion_results <- data.frame()
                #now <- Sys.time()
                # # ***** prepare file to print PNG images *****
                # only plot top 3 ions plus one more for sum
                image_frame_count <- 4
                #CairoPNG(filename=paste(plot_output_dir, "\\", input_peptide_sequence, "_", input_precursor_charge, "_", gsub('\\|', '_', input_protein_name), "_", curve_type ,"_", trunc(as.numeric(now)), ".png", sep=""), width=image_frame_count*400, height=400, bg="white", units="px")
                #CairoPNG(filename=paste(plot_output_dir, "\\", input_peptide_sequence, "_", input_precursor_charge, "_", curve_type ,"_", trunc(as.numeric(now)), ".png", sep=""), width=image_frame_count*400, height=400, bg="white", units="px")
                CairoPNG(filename=paste(plot_output_dir, "\\", input_peptide_sequence, "_", input_precursor_charge, "_", curve_type, ".png", sep=""), width=image_frame_count*400, height=400, bg="white", units="px")
                par(mfrow= c(1, image_frame_count))
                
                # Step 1:
                # Traverse fragment_ion_list, days, sample_groups and replicates
                # to evaluate the fragment_ion under the specific combination of day, sample_group, and replicate.
                for (current_ion in fragment_ion_list) {
                    for (current_day in days) {
                        for (current_sample in sample_groups) {
                            for (current_rep in replicates) {
                                current_set_count <- 0
                                light_area <- 0
                                heavy_area <- 0
                                theoretical_area <- 0
                                measured_area <- 0
                                calculated_area_ratio <- 0
                                
                                current_set <- QC_set[QC_set$fragment_ion==current_ion & QC_set$day==current_day & QC_set$sample_group==current_sample & QC_set$replicate==current_rep, ]
                                current_set_count <- nrow(current_set)
                                
                                if (current_set_count == 2) {
                                    light_area <- current_set[current_set$isotope_label_type=='light', 'area' ]
                                    heavy_area <- current_set[current_set$isotope_label_type=='heavy', 'area' ]

                                    if(curve_type=='forward') {
                                        theoretical_area <- heavy_area
                                        measured_area <- light_area
                                    }
                                    else {
                                        theoretical_area <- light_area
                                        measured_area <- heavy_area
                                    }
                                    if(theoretical_area==0 | is.na(theoretical_area) | is.na(measured_area)){
                                        calculated_area_ratio = NA
                                    }
                                    else {
                                        calculated_area_ratio <- measured_area/theoretical_area
                                    }
                                    fragment_ion_results <-  rbind(fragment_ion_results, data.frame(Peptide = input_peptide_sequence, Protein_Name = input_protein_name, Precursor_Charge = input_precursor_charge,
                                                                                                  fragment_ion = current_ion, day = current_day, sample_group = current_sample, replicate = current_rep, light_area=light_area, heavy_area=heavy_area,
                                                                                                  theoretical_area=theoretical_area, measured_area=measured_area,
                                                                                                  calculated_area_ratio=calculated_area_ratio, ion_category='individual') )
                                }
                                else {
                                }
                            } # end current_rep
                        } # end current_sample
                    } # end current_day
                }  # end current_ion
                # Step 2: repeat calculations for sum of ions 
                for (current_day in days) {
                    for (current_sample in sample_groups) {
                        for (current_rep in replicates) {
                            sum_light_area <- 0
                            sum_heavy_area <- 0
                            skipped_count <- 0
                            skip_current_sample <- 'true'

                            sum_theoretical_area <- 0
                            sum_measured_area <- 0
                            for (current_ion in fragment_ion_list) {
                                current_set_count <- 0
                                calculated_area_ratio <- 0
                                light_area <- 0
                                heavy_area <- 0
                                theoretical_area <- 0
                                measured_area <- 0
                                
                                current_set <- QC_set[QC_set$fragment_ion==current_ion & QC_set$day==current_day & QC_set$sample_group==current_sample & QC_set$replicate==current_rep, ]
                                
                                current_set_count <- nrow(current_set)
                                
                                if (current_set_count == 2) {
                                    light_area <- current_set[current_set$isotope_label_type=='light', 'area' ]
                                    heavy_area <- current_set[current_set$isotope_label_type=='heavy', 'area' ]
                                    if(curve_type=='forward'){
                                      theoretical_area <- heavy_area
                                      measured_area <- light_area
                                    }
                                    else {
                                      theoretical_area <- light_area
                                      measured_area <- heavy_area
                                    }

                                    add_to_sum <- 'false'

                                    if(theoretical_area==0 | is.na(theoretical_area) | is.na(measured_area)){
                                      skipped_count <- skipped_count + 1
                                    }
                                    else {
                                      sum_theoretical_area <- sum_theoretical_area + theoretical_area
                                      sum_measured_area <- sum_measured_area + measured_area

                                      sum_light_area <- sum_light_area + light_area
                                      sum_heavy_area <- sum_heavy_area + heavy_area

                                      add_to_sum <- 'true'
                                    }
                                    skip_current_sample <- 'false'
                                } 
                                else {
                                }
                            }
                            
                            if (skip_current_sample=='false'){
                                if(sum_theoretical_area==0){
                                  calculated_area_ratio <- NA
                                }
                                else {
                                  calculated_area_ratio <- sum_measured_area/sum_theoretical_area
                                }
                                fragment_ion_results <-  rbind(fragment_ion_results, data.frame(Peptide = input_peptide_sequence, Protein_Name = input_protein_name, Precursor_Charge = input_precursor_charge,
                                                                                                fragment_ion = 'sum', day = current_day, sample_group = current_sample, replicate = current_rep, light_area=sum_light_area, heavy_area=sum_heavy_area,
                                                                                                theoretical_area=sum_theoretical_area, measured_area=sum_measured_area,
                                                                                                calculated_area_ratio=calculated_area_ratio, ion_category='all') )
                            }
                        } # end current_rep
                    } # end current_sample
                } # end current_day
                # Step 3:
                # ***** calculate CV (Coefficient of Variation) *****
                LoMedHi <- c('Lo', 'Med', 'Hi')
                ions <- c(fragment_ion_list, 'sum')
                # make CV summary data frame
                CV_results <- data.frame()
                for (current_ion in ions) {
                  # define columns (make row for each fragment ion)
                    CV_results <-  rbind(CV_results, data.frame(fragment_ion=current_ion, low_intra_CV=NA, med_intra_CV=NA, high_intra_CV=NA,
                                                              low_inter_CV=NA, med_inter_CV=NA, high_inter_CV=NA,
                                                              low_total_CV=NA, med_total_CV=NA, high_total_CV=NA,
                                                              low_count=NA, med_count=NA, high_count=NA))
                }
                # *** intra-assay CV ***
                for (current_ion in ions) {
                    for (current_LoMedHi in LoMedHi) {
                        avg_intra_assay_CV <-0
                        individual_intra_assay_CVs <- c()
                        for (current_day in days) {
                            #for (current_rep in replicates) {
                            current_set <- fragment_ion_results[fragment_ion_results$fragment_ion==current_ion & fragment_ion_results$day==current_day & fragment_ion_results$sample_group==current_LoMedHi, ]
                            
                            # remove rows with a value of NA for calculated_area_ratio
                            current_set <- current_set[complete.cases(current_set[ , 'calculated_area_ratio' ]),]
                            
                            if (nrow(current_set) <= 1 ){
                                percent_CV <- NA
                            }
                            else{
                                percent_CV <- (sd(current_set$calculated_area_ratio))/(mean(current_set$calculated_area_ratio)) * 100
                                individual_intra_assay_CVs <- c(individual_intra_assay_CVs, percent_CV)
                            }
                            #  } # end current_rep
                        } # end current_day

                        if (length(individual_intra_assay_CVs)==0){
                            avg_CV <- NA
                            count <- 0
                        }
                        else{
                            avg_CV <- mean(individual_intra_assay_CVs, na.rm = TRUE)
                        }

                        if (current_LoMedHi=='Lo'){
                            CV_results[CV_results$fragment_ion==current_ion, 'low_intra_CV'] <- round(avg_CV, digits=1)
                        }
                        else if(current_LoMedHi=='Med'){
                            CV_results[CV_results$fragment_ion==current_ion, 'med_intra_CV'] <- round(avg_CV, digits=1)
                        }
                        else if(current_LoMedHi=='Hi'){
                            CV_results[CV_results$fragment_ion==current_ion, 'high_intra_CV'] <- round(avg_CV, digits=1)
                        }
                    } # end current_LoMedHi
                } # current_ion
                # *** inter-assay CV ***
                for (current_ion in ions) {
                    for (current_LoMedHi in LoMedHi) {
                        avg_inter_assay_CV <-0
                        individual_inter_assay_CVs <- c()
                        
                        avg_CV  <- 0
                        # for (current_day in days) {
                        for (current_rep in replicates) {           
                            current_set <- fragment_ion_results[fragment_ion_results$fragment_ion==current_ion & fragment_ion_results$replicate==current_rep & fragment_ion_results$sample_group==current_LoMedHi, ]
                            # remove rows with a value of NA for calculated_area_ratio
                            current_set <- current_set[complete.cases(current_set[ , 'calculated_area_ratio' ]),]
                            if (nrow(current_set) <= 1 ){
                              percent_CV <- NA
                            }
                            else{
                              percent_CV <- (sd(current_set$calculated_area_ratio))/(mean(current_set$calculated_area_ratio)) * 100
                              individual_inter_assay_CVs <- c(individual_inter_assay_CVs, percent_CV)
                            }
                        } # end current_rep
                        # } # end current_day
                        
                        if (length(individual_inter_assay_CVs)==0){
                            avg_CV <- NA
                            count <- 0
                        }
                        else{
                            avg_CV <- mean(individual_inter_assay_CVs, na.rm = TRUE)
                        }
                        if (current_LoMedHi=='Lo'){
                            CV_results[CV_results$fragment_ion==current_ion, 'low_inter_CV'] <- round(avg_CV, digits=1)
                        }
                        else if(current_LoMedHi=='Med'){
                            CV_results[CV_results$fragment_ion==current_ion, 'med_inter_CV'] <- round(avg_CV, digits=1)
                        }
                        else if(current_LoMedHi=='Hi'){
                            CV_results[CV_results$fragment_ion==current_ion, 'high_inter_CV'] <- round(avg_CV, digits=1)
                        }
                    } # end current_LoMedHi
                } # current_ion

                # *** END: inter-assay CV ***
                # calculate total variability
                CV_results[ , 'med_total_CV'] <- round(sqrt((CV_results[ , 'med_intra_CV'])*(CV_results[ , 'med_intra_CV']) + (CV_results[ , 'med_inter_CV'])*(CV_results[ , 'med_inter_CV'])), digits=1)
                CV_results[ , 'high_total_CV'] <- round(sqrt((CV_results[ , 'high_intra_CV'])*(CV_results[ , 'high_intra_CV']) + (CV_results[ , 'high_inter_CV'])*(CV_results[ , 'high_inter_CV'])), digits=1)
                CV_results[ , 'low_total_CV'] <- round(sqrt((CV_results[ , 'low_intra_CV'])*(CV_results[ , 'low_intra_CV']) + (CV_results[ , 'low_inter_CV'])*(CV_results[ , 'low_inter_CV'])), digits=1)

                # determine counts
                for (current_ion in ions) {
                    CV_results[CV_results$fragment_ion==current_ion, 'high_count'] <- nrow(fragment_ion_results[fragment_ion_results$sample_group=='Hi' & fragment_ion_results$fragment_ion==current_ion & !is.na(fragment_ion_results$calculated_area_ratio), ] )
                    CV_results[CV_results$fragment_ion==current_ion, 'med_count'] <- nrow(fragment_ion_results[fragment_ion_results$sample_group=='Med' & fragment_ion_results$fragment_ion==current_ion & !is.na(fragment_ion_results$calculated_area_ratio), ] )
                    CV_results[CV_results$fragment_ion==current_ion, 'low_count'] <- nrow(fragment_ion_results[fragment_ion_results$sample_group=='Lo' & fragment_ion_results$fragment_ion==current_ion & !is.na(fragment_ion_results$calculated_area_ratio), ] )
                }

                ions_to_plot <- c()
                ions_in_table <- c()

                # determine fragment ions to plot
                if (length(ions) <= 4 ) {
                    results_to_plot <- CV_results[CV_results$fragment_ion!='sum' & !is.na(CV_results$low_total_CV) & !is.na(CV_results$med_total_CV) & !is.na(CV_results$high_total_CV) , ]
                    ions_list <- as.character(results_to_plot[ , 'fragment_ion'])
                    results_to_plot2 <- CV_results[CV_results$fragment_ion=='sum' & !is.na(CV_results$low_total_CV) & !is.na(CV_results$med_total_CV) & !is.na(CV_results$high_total_CV) , ]
                    ions_to_plot <- c(ions_list, as.character(results_to_plot2[ , 'fragment_ion']))
                    ions_in_table <- ions_to_plot
                } else {
                    results_to_plot <- CV_results[CV_results$fragment_ion!='sum' & !is.na(CV_results$low_total_CV) & !is.na(CV_results$med_total_CV) & !is.na(CV_results$high_total_CV) , ]
                    ions_list <- as.character(results_to_plot[ , 'fragment_ion'])
                    # new sort to get Top 3 plots
                    results_to_plot <- results_to_plot[order(results_to_plot$med_total_CV, results_to_plot$low_total_CV, results_to_plot$high_total_CV), ]
                    three_lowest_med_total_CV <- head(results_to_plot, 3)
                    three_lowest_ions <- as.character(three_lowest_med_total_CV[ , 'fragment_ion'])
                    results_to_plot2 <- CV_results[CV_results$fragment_ion=='sum' & !is.na(CV_results$low_total_CV) & !is.na(CV_results$med_total_CV) & !is.na(CV_results$high_total_CV) , ]
                    ions_to_plot  <- c(three_lowest_ions, as.character(results_to_plot2[ , 'fragment_ion']))
                    ions_in_table <- c(ions_list, as.character(results_to_plot2[ , 'fragment_ion']))
                }
                
                # Check wether ions_to_plot is blank or now. If it's blank, an error message will be thrown.
                if (length(ions_to_plot) == 0) {
                    errorType <- "Error"
                    errorSubtype <- "Missing points"
                    errorReason <- "In repeatability graph, there is no fragment ions and their CVs can't be caluclated."
                    errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge), rep('', colNumber-5)), collapse='\t')
                    cat(errorInfor)
                    cat('\n')
                    next
                }
                
                errorReason_missing_point_on_day <- c()
                errorReason_missing_concerntration_on_day <- c()
                errorReason_missing_replicate_point_on_day <- c()
                mean_for_concentration_by_ion <- data.frame(Peptide=as.character(), Protein_Name=as.character(), Precursor_Charge=as.character(), fragment_ion=as.character(), day=as.numeric(), mean_Hi=as.numeric(), mean_Med=as.numeric(), mean_Lo=as.numeric())
                for (current_plot_ion in ions_to_plot ) {
                    plot_fragment_ion_results <- fragment_ion_results[!is.na(fragment_ion_results$calculated_area_ratio), ]
                    plot_fragment_ion_results <- plot_fragment_ion_results[plot_fragment_ion_results$fragment_ion == current_plot_ion, ]
                    plot_days <- sort(unique(plot_fragment_ion_results[ , 'day']))
                    # Capture the situation where no points are observed at specific day.
                    daysTmp <- setdiff(days,plot_days)
                    if (length(daysTmp) > 0) {
                        errorReason_tmp <- paste("For fragment ion ", current_plot_ion, ", there are no points on day ",  paste(daysTmp, collapse = ', '), sep="")
                        errorReason_missing_point_on_day <- c(errorReason_missing_point_on_day, errorReason_tmp)
                    }

                    # Do not make plots for fragment ions with no data
                    if (nrow(plot_fragment_ion_results) != 0) {
                        # make QC plot for current ion
                        plot_QC(plot_fragment_ion_results, input_peptide_sequence, current_plot_ion, plot_days)
                        # check whether points are diplayed by Hi, Med, Lo from top to bottom.
                        # The strategy is to compare the mean of all the points in each concentration and the mean values should be Hi > Med > Lo for each fragment ion.
                        # calculate the mean values for the frament ion in each concentration
                        if (current_plot_ion != 'sum') {
                            tmp1 <- ddply(plot_fragment_ion_results, .(Peptide, Protein_Name, Precursor_Charge, fragment_ion, day, sample_group), 
                                      summarize, mean=mean(calculated_area_ratio, na.rm= TRUE))
                            tmp2 <- dcast(tmp1, Peptide + Protein_Name + Precursor_Charge + fragment_ion + day ~ sample_group, value.var="mean")
                            colnames(tmp2)[which(colnames(tmp2)=='Hi')] <- 'mean_Hi'
                            colnames(tmp2)[which(colnames(tmp2)=='Med')] <- 'mean_Med'
                            colnames(tmp2)[which(colnames(tmp2)=='Lo')] <- 'mean_Lo'
                            mean_for_concentration_by_ion <- suppressWarnings(dplyr::bind_rows(mean_for_concentration_by_ion, tmp2))
                        }
                    }
                    
                    # Transverse each day and count the number of replicates, if the number is less than 3, a warning will arise.
                    for (plot_day in plot_days) {
                        plot_fragment_ion_results_tmp1 <- plot_fragment_ion_results[plot_fragment_ion_results$day == plot_day, ]
                        concentration_list <- unique(plot_fragment_ion_results_tmp1$sample_group)
                        concentrationsTmp <- setdiff(sample_groups,concentration_list) 
                        if (length(concentrationsTmp) > 0) {
                            errorReason_tmp <- paste("For fragment ion ", current_plot_ion, ", not all of three concentrations ('Hi', 'Med', 'Lo') exist on day ",  plot_day, sep="")
                            errorReason_missing_concerntration_on_day <- c(errorReason_missing_concerntration_on_day, errorReason_tmp)
                        }
                        concentration_with_less_replicates <- c()
                        for (concentration_tmp in concentration_list) {
                            plot_fragment_ion_results_tmp2 <- plot_fragment_ion_results_tmp1[plot_fragment_ion_results_tmp1$sample_group == concentration_tmp, ]
                            replicate_count <- length(unique(plot_fragment_ion_results_tmp2$replicate))
                            if (replicate_count < 3) {
                                concentration_with_less_replicates <- c(concentration_with_less_replicates, concentration_tmp)
                            }
                        }
                        if (length(concentration_with_less_replicates) > 0) {
                            errorReason_tmp <- paste("For fragment ion ", current_plot_ion, ", there are less than 3 replicates in concentration(s) ", paste(concentration_with_less_replicates, collapse=', '), " on day ",  plot_day, sep="")
                            errorReason_missing_replicate_point_on_day <- c(errorReason_missing_replicate_point_on_day, errorReason_tmp)
                        }
                    }
                }
                dev.off()
                
                if (length(errorReason_missing_point_on_day) > 0) {
                    errorType <- "Warning"
                    errorSubtype <- "Missing points"
                    errorReason <- paste(paste(errorReason_missing_point_on_day, collapse='. '), '.', sep="")
                    #errorInfor <- paste(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge, '', '', '', '', '', '', '', '', sep='\t')
                    errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge), rep('', colNumber-5)), collapse='\t')
                    cat(errorInfor)
                    cat('\n')
                }
                
                if (length(errorReason_missing_concerntration_on_day) >0 ) {
                    errorType <- "Warning"
                    errorSubtype <- "Missing points"
                    errorReason <- paste(paste(errorReason_missing_concerntration_on_day, collapse='. '), '.', sep="")
                    #errorInfor <- paste(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge, '', '', '', '', '', '', '', '', sep='\t')
                    errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge), rep('', colNumber-5)), collapse='\t')
                    cat(errorInfor)
                    cat('\n')
                }
                
                if (length(errorReason_missing_replicate_point_on_day) >0 ) {
                    errorType <- "Warning"
                    errorSubtype <- "Missing points"
                    errorReason <- paste(paste(errorReason_missing_replicate_point_on_day, collapse='. '), '.', sep="")
                    #errorInfor <- paste(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge, '', '', '', '', '', '', '', '', sep='\t')
                    errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge), rep('', colNumber-5)), collapse='\t')
                    cat(errorInfor)
                    cat('\n')
                }
                
                normal_count <- 0
                abnormal_count <- 0
                abnormal_fragment_ion <- c()
                
                for (i in 1:nrow(mean_for_concentration_by_ion)) {
                    mean_Hi <- mean_for_concentration_by_ion[i, 'mean_Hi']
                    mean_Med <- mean_for_concentration_by_ion[i, 'mean_Med']
                    mean_Lo <- mean_for_concentration_by_ion[i, 'mean_Lo']
                    compare_1 <- mean_Hi >= mean_Med
                    compare_2 <- mean_Med >= mean_Lo
                    compare_1_na_status <- !is.na(compare_1)
                    compare_2_na_status <- !is.na(compare_2)
                    
                    if (compare_1_na_status | compare_2_na_status) {
                      if (compare_1_na_status & compare_2_na_status) {
                        if (compare_1 & compare_2) {
                          normal_count <- normal_count + 1
                        } else {
                          abnormal_count <- abnormal_count + 1
                          abnormal_fragment_ion <- c(abnormal_fragment_ion, mean_for_concentration_by_ion[i, 'fragment_ion'])
                        }
                      } else if (compare_1_na_status) {
                        if (compare_1) {
                          normal_count <- normal_count + 1
                        } else {
                          abnormal_count <- abnormal_count + 1
                          abnormal_fragment_ion <- c(abnormal_fragment_ion, mean_for_concentration_by_ion[i, 'fragment_ion'])
                        }
                      } else {
                        if (compare_2) {
                          normal_count <- normal_count + 1
                        } else {
                          abnormal_count <- abnormal_count + 1
                          abnormal_fragment_ion <- c(abnormal_fragment_ion, mean_for_concentration_by_ion[i, 'fragment_ion'])
                        }
                      }
                    }
                }
                abnormal_fragment_ion <- unique(abnormal_fragment_ion)
                
                if (abnormal_count > normal_count) {
                    errorType <- "Warning"
                    errorSubtype <- "Abnormal display of points"
                    errorReason <- paste('For fragment ions: ', paste(abnormal_fragment_ion, collapse=', '), ', some points are not displayed by the sequence of Hi, Med and low from top to bottom.')
                    #errorInfor <- paste(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge, '', '', '', '', '', '', '', '', sep='\t')
                    errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge), rep('', colNumber-5)), collapse='\t')
                    cat(errorInfor)
                    cat('\n')
                }
                

                # output to files
                CV_results_sorted <- CV_results[match(ions_in_table , CV_results$fragment_ion), ]
                write.table(CV_results_sorted, file=paste(plot_output_dir, "\\", input_peptide_sequence, "_", input_precursor_charge, "_", curve_type ,"_CV_results", ".tsv", sep=""), sep = "\t", qmethod = "double", col.names=TRUE, row.names=FALSE, quote=FALSE)
                #write.csv(CV_results_sorted, file=paste(plot_output_dir, "\\", input_peptide_sequence, "_", input_precursor_charge, "_", curve_type ,"_CV_results_", trunc(as.numeric(now)), ".csv", sep=""), row.names = FALSE)
                # Top three plots were selected based on med_total_CV, low_total_CV and high_total_CV.
                # med_total_CV, low_total_CV and high_total_CV of the will be used to evaluate the quality.
                CV_results_sub <- subset(CV_results_sorted, select = c(fragment_ion, low_total_CV, med_total_CV, high_total_CV))
                # Evaluate the quality based on the total_CV of the fragment_ion "sum" and the rest fragment ions in ions_to_plot
                CV_results_sub_all <- CV_results_sub[CV_results_sub$fragment_ion == "sum", ]
                CV_results_sub_individual <- CV_results_sub[CV_results_sub$fragment_ion != "sum" & CV_results_sub$fragment_ion %in% ions_to_plot, ]
                rownames(CV_results_sub_all) <- CV_results_sub_all$fragment_ion
                rownames(CV_results_sub_individual) <- CV_results_sub_individual$fragment_ion
                # For sum, low_total_CV, med_total_CV and high_total_CV will be used to compared with cv_threshold_all
                CV_results_sub_all <- CV_results_sub_all[,-c(1)]
                # For individual ions, med_total_CV and high_total_CV will be used to compared with cv_threshold_individual
                CV_results_sub_individual <- CV_results_sub_individual[, -c(1, 2)]

                if (FALSE) {
                    if (any(CV_results_sub_all > cv_threshold_all) | any(CV_results_sub_individual > cv_threshold_individual)) {
                        index_all_array <- which(CV_results_sub_all > cv_threshold_all)
                        index_individual_array <- which(matrix(CV_results_sub_individual > cv_threshold_individual, ncol=ncol(CV_results_sub_individual)), arr.ind=TRUE)
                        
                        # if index_individual_array is 1*2 matrix, index_individual_array shouldn't be ordered by the first column.
                        if (length(index_individual_array) == 2) {
                            index_individual_row <- unique(index_individual_array[,1])
                        } else {
                            index_individual_array <- index_individual_array[order(index_individual_array[,1]), ]
                            index_individual_row <- unique(index_individual_array[,1])
                        }
                        
                        concentration_cv_array_all <- c()
                        for (index_tmp in index_all_array) {
                            concentration_cv_array_all <- c(concentration_cv_array_all, paste('sum-', strsplit(colnames(CV_results_sub_all)[index_tmp], split='_')[[1]][1], ': ', paste(CV_results_sub_all[1, index_tmp], '%', sep=''), sep=''))
                        }
                        
                        concentration_cv_array_individual <- c()
                        for (index_tmp in index_individual_row) {
                            fragment_ion_selected <- rownames(CV_results_sub_individual)[index_tmp]
                            fragment_ion_selected_id <- which(rownames(CV_results_sub_individual)==fragment_ion_selected)
                            index1_array <- index_individual_array[, 2][index_individual_array[,1]== index_tmp]
                            for (index1 in index1_array) {
                                concentration_cv_array_individual <- c(concentration_cv_array_individual, paste(fragment_ion_selected, '-', strsplit(colnames(CV_results_sub_individual)[index1], split='_')[[1]][1], ': ', paste(CV_results_sub_individual[fragment_ion_selected_id, index1], '%', sep=''), sep=''))
                            }
                        }
                        if (length(concentration_cv_array_all) >= 1) {
                            errorReason_all <- paste('(', paste(concentration_cv_array_all, collapse = ', '), ' larger than the threshold of ',cv_threshold_all, '%)', sep ='')
                        } else {
                            errorReason_all <- ''
                        }
                        
                        if (length(concentration_cv_array_individual) >= 1) {
                            errorReason_individual <- paste('(', paste(concentration_cv_array_individual, collapse = ', '), ' larger than the threshold of ', cv_threshold_individual, '%)', sep='')
                        } else {
                          errorReason_individual <- ''
                        }
                        errorType <- "Warning"
                        errorSubtype <- "Bad distribution of points"
                        errorReason <- paste("The total coefficient of variation of fragment ions under concentration(s) is(are) large ", errorReason_all, ' ', errorReason_individual, '.', sep="")
                        #errorInfor <- paste(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge, '', '', '', '', '', '', '', '', sep='\t')
                        errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge), rep('', colNumber-5)), collapse='\t')
                        cat(errorInfor)
                        cat('\n')
                    }
                }
                CV_results_sub_2 <- subset(CV_results_sorted, select = c(fragment_ion, low_intra_CV, med_intra_CV, high_intra_CV, low_inter_CV, med_inter_CV, high_inter_CV, low_total_CV, med_total_CV, high_total_CV))
                CV_results_sub_2_plot <- CV_results_sub_2[CV_results_sub_2$fragment_ion  %in% ions_to_plot, ]
                rownames(CV_results_sub_2_plot) <- CV_results_sub_2_plot$fragment_ion
                CV_results_sub_2_plot <- CV_results_sub_2_plot[,-c(1)]
                if (any(CV_results_sub_2_plot > cv_threshold)) {
                    index_individual_array <- which(matrix(CV_results_sub_2_plot > cv_threshold, ncol=ncol(CV_results_sub_2_plot)), arr.ind=TRUE)
                    # if index_individual_array is 1*2 matrix, index_individual_array shouldn't be ordered by the first column.
                    if (length(index_individual_array) == 2) {
                        index_individual_row <- unique(index_individual_array[,1])
                    } else {
                        index_individual_array <- index_individual_array[order(index_individual_array[,1]), ]
                        index_individual_row <- unique(index_individual_array[,1])
                    }
                    errorReason_tmp <- c()
                    for (index_tmp in index_individual_row) {
                        fragment_ion_selected <- rownames(CV_results_sub_2_plot)[index_tmp]
                        fragment_ion_selected_id <- which(rownames(CV_results_sub_2_plot)==fragment_ion_selected)
                        index1_array <- index_individual_array[, 2][index_individual_array[,1]== index_tmp]
                        intra_assay_abonormal <- c()
                        inter_assay_abonormal <- c()
                        total_abnormal <- c()
                        for (index1 in index1_array) {
                            strsplit_tmp <- strsplit(colnames(CV_results_sub_2_plot)[index1], split='_')[[1]]
                            concentraion_term <- strsplit_tmp[1]
                            cv_type_term <- strsplit_tmp[2]
                            if (cv_type_term == "intra") {
                                intra_assay_abonormal <- c(intra_assay_abonormal, concentraion_term)
                            } else if (cv_type_term == "inter") {
                                inter_assay_abonormal <- c(inter_assay_abonormal, concentraion_term)
                            } else {
                                total_abnormal <- c(total_abnormal, concentraion_term)
                            }
                        }
                        errorReason_message <- paste('For fragment ion ', fragment_ion_selected, ',', sep='')
                        if (length(intra_assay_abonormal) > 0) {
                            errorReason_message_tmp <- paste(' intra_assay CV of concentration(s) of ', paste(intra_assay_abonormal, collapse=', '), ' is(are) larger than the threshold of ', cv_threshold, '%', sep='')
                            errorReason_message <- paste(errorReason_message, errorReason_message_tmp, ',', sep='')
                        }
                        if (length(inter_assay_abonormal) > 0) {
                            errorReason_message_tmp <- paste(' inter_assay CV of concentration(s) of ', paste(inter_assay_abonormal, collapse=', '), ' is(are) larger than the threshold of ', cv_threshold, '%',  sep='')
                            errorReason_message <- paste(errorReason_message, errorReason_message_tmp, ',', sep='')
                        }
                        if (length(total_abnormal) > 0) {
                            errorReason_message_tmp <- paste(' total CV of concentration(s) of ', paste(total_abnormal, collapse=', '), ' is(are) larger than the threshold of ', cv_threshold, '%',  sep='')
                            errorReason_message <- paste(errorReason_message, errorReason_message_tmp, '.', sep='')
                        }
                        errorReason_tmp <- c(errorReason_tmp, errorReason_message)
                    }
                    errorType <- "Warning"
                    errorSubtype <- "Bad distribution of points"
                    errorReason <- paste(errorReason_tmp, collapse =" ")
                    errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge), rep('', colNumber-5)), collapse='\t')
                    cat(errorInfor)
                    cat('\n')
                }
                if (FALSE) {
                    # This only applys to the ions to plot
                    fragment_ion_results_tmp <- fragment_ion_results[fragment_ion_results$fragment_ion %in% ions_to_plot, ]
                    fragment_ion_results_tmp2 <- as.data.frame(fragment_ion_results_tmp %>%
                                                     group_by(Peptide, Protein_Name, Precursor_Charge, fragment_ion) %>%
                                                     summarise(Ratio_mean=mean(calculated_area_ratio), Ratio_max=max(calculated_area_ratio), Ratio_min=min(calculated_area_ratio)))
                    fragment_ion_results_tmp3 <- as.data.frame(fragment_ion_results_tmp %>%
                                                                 group_by(Peptide, Protein_Name, Precursor_Charge, fragment_ion, sample_group) %>%
                                                                 summarise(Ratio_mean=mean(calculated_area_ratio), Ratio_max=max(calculated_area_ratio), Ratio_min=min(calculated_area_ratio)))
                    fragment_ion_results_tmp4 <- fragment_ion_results_tmp3[!((fragment_ion_results_tmp3$Ratio_max <= 1.3*fragment_ion_results_tmp3$Ratio_mean & fragment_ion_results_tmp3$Ratio_max >= 0.7*fragment_ion_results_tmp3$Ratio_mean) & (fragment_ion_results_tmp3$Ratio_min <= 1.3*fragment_ion_results_tmp3$Ratio_mean & fragment_ion_results_tmp3$Ratio_min >= 0.7*fragment_ion_results_tmp3$Ratio_mean)),]
                    if (nrow(fragment_ion_results_tmp4) > 0) {
                        unique_fragment_ion <- unique(fragment_ion_results_tmp4$fragment_ion)
                        errorType <- "Warning"
                        errorSubtype <- "High variance"
                        errorReasonTmp <- c()
                        for (fragment_ion_tmp in unique_fragment_ion) {
                            errorReasonTmp <- c(errorReasonTmp, paste("for fragment ion ", fragment_ion_tmp, " at the concentration(s) of ", paste(fragment_ion_results_tmp4[fragment_ion_results_tmp4$fragment_ion == fragment_ion_tmp, ]$sample_group, collapse=', '), sep = ''))
                        }
                        errorReason <- paste('When checking the samples with the concentrations of Low, Medium and High, ', paste(errorReasonTmp, collapse = ', '), ', not all of the area ratios are within 30% of the mean.')
                        errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge), rep('', colNumber-5)), collapse='\t')
                        cat(errorInfor)
                        cat('\n')
                    }
                }
                
                if (allInjectionSwitch) {
                    # This only applies to the ions to plot
                    fragment_ion_results_tmp <- fragment_ion_results[fragment_ion_results$fragment_ion %in% ions_to_plot, ]
                    fragment_ion_results_tmp_sum <- fragment_ion_results_tmp[fragment_ion_results_tmp$fragment_ion=='sum', ]
                    fragment_ion_results_tmp_idividual <- fragment_ion_results_tmp[fragment_ion_results_tmp$fragment_ion!='sum', ]
                    
                    heavyArea_transition_ratio <- c()
                    lightArea_transition_ratio <- c()
                    for (i in 1:nrow(fragment_ion_results_tmp_idividual)) {
                        row <- fragment_ion_results_tmp_idividual[i,]
                        row_sum_tmp <- fragment_ion_results_tmp_sum[fragment_ion_results_tmp_sum$day==row$day & fragment_ion_results_tmp_sum$sample_group==row$sample_group & fragment_ion_results_tmp_sum$replicate==row$replicate, ]
                        heavyArea_transition_ratio <- c(heavyArea_transition_ratio, row$heavy_area/row_sum_tmp$heavy_area)
                        lightArea_transition_ratio <- c(lightArea_transition_ratio, row$light_area/row_sum_tmp$light_area)
                    }
                    fragment_ion_results_tmp_idividual$heavyArea_transition_ratio <- heavyArea_transition_ratio
                    fragment_ion_results_tmp_idividual$lightArea_transition_ratio <- lightArea_transition_ratio
                    
                    # Step1: Check whether there are all samples above the LLOQ, no all of transition ratio are within 30% of the mean.
                    fragment_ion_results_tmp_idividual_2 <- as.data.frame(fragment_ion_results_tmp_idividual %>%
                                                            dplyr::group_by(Protein_Name, Peptide, Precursor_Charge, fragment_ion) %>%
                                                            dplyr::summarise(heavyArea_transition_ratio_mean=mean(heavyArea_transition_ratio, na.rm=TRUE), heavyArea_transition_ratio_max=max(heavyArea_transition_ratio, na.rm=TRUE), heavyArea_transition_ratio_min=min(heavyArea_transition_ratio, na.rm=TRUE),
                                                                      lightArea_transition_ratio_mean=mean(lightArea_transition_ratio, na.rm=TRUE), lightArea_transition_ratio_max=max(lightArea_transition_ratio, na.rm=TRUE), lightArea_transition_ratio_min=min(lightArea_transition_ratio, na.rm=TRUE)))
                    # Remove the rows with NAs
                    row.has.na <- apply(fragment_ion_results_tmp_idividual_2, 1, function(x){any(is.na(x))})
                    fragment_ion_results_tmp_idividual_2 <-  fragment_ion_results_tmp_idividual_2[!row.has.na,]
                    
                    if (nrow(fragment_ion_results_tmp_idividual_2) > 0) {
                        fragment_ion_results_tmp_idividual_2_heavy <- fragment_ion_results_tmp_idividual_2[!((fragment_ion_results_tmp_idividual_2$heavyArea_transition_ratio_max <= 1.3*fragment_ion_results_tmp_idividual_2$heavyArea_transition_ratio_mean & fragment_ion_results_tmp_idividual_2$heavyArea_transition_ratio_max >= 0.7*fragment_ion_results_tmp_idividual_2$heavyArea_transition_ratio_mean) & (fragment_ion_results_tmp_idividual_2$heavyArea_transition_ratio_min <= 1.3*fragment_ion_results_tmp_idividual_2$heavyArea_transition_ratio_mean & fragment_ion_results_tmp_idividual_2$heavyArea_transition_ratio_min >= 0.7*fragment_ion_results_tmp_idividual_2$heavyArea_transition_ratio_mean)),]
                        fragment_ion_results_tmp_idividual_2_light <- fragment_ion_results_tmp_idividual_2[!((fragment_ion_results_tmp_idividual_2$lightArea_transition_ratio_max <= 1.3*fragment_ion_results_tmp_idividual_2$lightArea_transition_ratio_mean & fragment_ion_results_tmp_idividual_2$lightArea_transition_ratio_max >= 0.7*fragment_ion_results_tmp_idividual_2$lightArea_transition_ratio_mean) & (fragment_ion_results_tmp_idividual_2$lightArea_transition_ratio_min <= 1.3*fragment_ion_results_tmp_idividual_2$lightArea_transition_ratio_mean & fragment_ion_results_tmp_idividual_2$lightArea_transition_ratio_min >= 0.7*fragment_ion_results_tmp_idividual_2$lightArea_transition_ratio_mean)), ]
                        
                        if (nrow(fragment_ion_results_tmp_idividual_2_heavy)+nrow(fragment_ion_results_tmp_idividual_2_light) > 0) {
                            fragment_ion_results_tmp_idividual_2_trans <- data.frame(fragment_ion=as.character(), isotope_warining=as.character(), stringsAsFactors = FALSE)
                            if (nrow(fragment_ion_results_tmp_idividual_2_heavy) > 0) {
                                thisPeptide_4_trans_tmp <- data.frame(fragment_ion=fragment_ion_results_tmp_idividual_2_heavy$fragment_ion, isotope_warining=rep("heavy", nrow(fragment_ion_results_tmp_idividual_2_heavy)))
                                fragment_ion_results_tmp_idividual_2_trans <- rbind(fragment_ion_results_tmp_idividual_2_trans, thisPeptide_4_trans_tmp)
                            }
                            if (nrow(fragment_ion_results_tmp_idividual_2_light) > 0) {
                                for (i in 1:nrow(fragment_ion_results_tmp_idividual_2_light)) {
                                if (fragment_ion_results_tmp_idividual_2_light$fragment_ion[i] %in% fragment_ion_results_tmp_idividual_2_trans$fragment_ion) {
                                    fragment_ion_results_tmp_idividual_2_trans$isotope_warining <- as.character(fragment_ion_results_tmp_idividual_2_trans$isotope_warining)
                                    fragment_ion_results_tmp_idividual_2_trans$isotope_warining[fragment_ion_results_tmp_idividual_2_trans$fragment_ion == fragment_ion_results_tmp_idividual_2_light$fragment_ion[i]] <- 'heavy or light'
                                } else {
                                    thisPeptide_4_trans_tmp <- data.frame(fragment_ion=fragment_ion_results_tmp_idividual_2_light$fragment_ion[i], isotope_warining='light')
                                    fragment_ion_results_tmp_idividual_2_trans <- rbind(fragment_ion_results_tmp_idividual_2_trans, thisPeptide_4_trans_tmp)
                                }
                                }
                            }
                            errorType <- "Warning"
                            errorSubtype <- "High variance"
                            errorReasonTmp <- c()
                            for (i in 1:nrow(fragment_ion_results_tmp_idividual_2_trans)) {
                                errorReasonTmp <- c(errorReasonTmp, paste("for fragment ion ", fragment_ion_results_tmp_idividual_2_trans[i, ]$fragment_ion, ", not all of the transition ratios from ", fragment_ion_results_tmp_idividual_2_trans[i, ]$isotope_warining, ' isotope labeled peptide are within 30% of the mean', sep = ''))
                            }
                            errorReason <- paste('When checking the samples with the concentrations of Low, Medium and High, ', paste(errorReasonTmp, collapse = ', '), '.', sep='')
                            errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge), rep('', colNumber-5)), collapse='\t')
                            cat(errorInfor)
                            cat('\n')
                        }
                        
                        # Step 2: for each ion transition, take the transformation + t-test approach for the calculation of p value to determine if there is a significant difference between the means of two groups (heavy and mean).
                        fragment_ion_results_tmp_idividual$heavyArea_transition_ratio_log <- suppressWarnings(log(fragment_ion_results_tmp_idividual$heavyArea_transition_ratio/(1-fragment_ion_results_tmp_idividual$heavyArea_transition_ratio)))
                        fragment_ion_results_tmp_idividual$lightArea_transition_ratio_log <- suppressWarnings(log(fragment_ion_results_tmp_idividual$lightArea_transition_ratio/(1-fragment_ion_results_tmp_idividual$lightArea_transition_ratio)))
                        
                        # Replace Inf and -Inf with NA
                        fragment_ion_results_tmp_idividual$heavyArea_transition_ratio_log[which(fragment_ion_results_tmp_idividual$heavyArea_transition_ratio_log == Inf | fragment_ion_results_tmp_idividual$heavyArea_transition_ratio_log == -Inf)] <- NA
                        fragment_ion_results_tmp_idividual$lightArea_transition_ratio_log[which(fragment_ion_results_tmp_idividual$lightArea_transition_ratio_log == Inf | fragment_ion_results_tmp_idividual$lightArea_transition_ratio_log == -Inf)] <- NA
                        
                        # Since t test can't be applied to the vectors with one element.
                        # Remove fragment ion with only one data point
                        fragment_ion_to_remove <- c()
                        for (fragmention_ion_to_check in unique(fragment_ion_results_tmp_idividual$fragment_ion)) {
                          if (nrow(fragment_ion_results_tmp_idividual[fragment_ion_results_tmp_idividual$fragment_ion==fragmention_ion_to_check, ]) == 1) {
                            fragment_ion_to_remove <- c(fragment_ion_to_remove, fragmention_ion_to_check)
                          }
                        }
                        fragment_ion_results_tmp_idividual <- fragment_ion_results_tmp_idividual[!(fragment_ion_results_tmp_idividual$fragment_ion %in% fragment_ion_to_remove), ]
                        if (nrow(fragment_ion_results_tmp_idividual) > 1) {
                            try ({
                                thisPeptide_5 <- as.data.frame(fragment_ion_results_tmp_idividual %>%
                                                                 dplyr::group_by(Protein_Name, Peptide, Precursor_Charge, fragment_ion) %>%
                                                                 dplyr::summarise(heavyArea_transition_ratio_log_mean=mean(heavyArea_transition_ratio_log, na.rm=TRUE), lightArea_transition_ratio_log_mean=mean(lightArea_transition_ratio_log, na.rm=TRUE), pValue.log=t.test(heavyArea_transition_ratio_log, lightArea_transition_ratio_log)$p.value, heavyArea_transition_ratio_mean=mean(heavyArea_transition_ratio, na.rm=TRUE), lightArea_transition_ratio_mean=mean(lightArea_transition_ratio, na.rm=TRUE), pValue=t.test(heavyArea_transition_ratio, lightArea_transition_ratio)$p.value))
                                thisPeptide_5$Internal_standard <- rep(internal_standard, nrow(thisPeptide_5))
                                if (internal_standard=='light') {
                                    thisPeptide_5$relative_difference <- (thisPeptide_5$heavyArea_transition_ratio_mean-thisPeptide_5$lightArea_transition_ratio_mean)/thisPeptide_5$lightArea_transition_ratio_mean
                                } else {
                                    thisPeptide_5$relative_difference <- (thisPeptide_5$lightArea_transition_ratio_mean-thisPeptide_5$heavyArea_transition_ratio_mean)/thisPeptide_5$heavyArea_transition_ratio_mean
                                }
                                
                                thisPeptide_5 <- thisPeptide_5[thisPeptide_5$pValue < pValue_threshold & abs(thisPeptide_5$relative_difference > relative_difference_threshold), ]
                                if (nrow(thisPeptide_5) > 0) {
                                    errorType <- "Warning"
                                    errorSubtype <- "High variance"
                                    errorReasonTmp <- c()
                                    for (i in 1:nrow(thisPeptide_5)) {
                                        errorReasonTmp <- c(errorReasonTmp, paste("for fragment ion ", thisPeptide_5[i, ]$fragment_ion, ", the mean of transition ratios from the heavy isotope labeled peptide (", round(thisPeptide_5[i, ]$heavyArea_transition_ratio_mean, 4), ') is significantly different from the mean of transition ratios from the light isotope labeled peptide (', round(thisPeptide_5[i, ]$lightArea_transition_ratio_mean, 4), ') and the relative difference (comparing to the internal standard type: ' ,thisPeptide_5[i, ]$Internal_standard, ') is ', paste(round(100*thisPeptide_5[i, ]$relative_difference, 4), "%", sep=""), sep = ''))
                                    }
                                    errorReason <- paste('When checking the samples with the concentrations of Low, Medium and High, ', paste(errorReasonTmp, collapse = ', '), '.', sep='')
                                    errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge), rep('', colNumber-5)), collapse='\t')
                                    cat(errorInfor)
                                    cat('\n')
                                }
                            }, silent = TRUE)
                        }
                    }
                }
            }
        }
    } 
}