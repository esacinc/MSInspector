# This sample code is used to qc skyline document data from experiment 5 of Assay Portal.
# It's modified based on:
#   esac-panorama-master\experiment-5\code\Experiment_5_for_integration_Panorama.R

# UPDATE: For experiment 5, removed all sample group info, as there should be a single sample
# In the updated skyline template of experiment 5, it should have heavy as the internal standard type. (The default internal standard type is heavy, no matter whether it's set or not)

suppressWarnings(suppressMessages(library(Cairo)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(evaluate)))
suppressWarnings(suppressMessages(require(reshape2)))

plot_QC <- function(plot_fragment_ion_results, input_peptide_sequence, current_ion, days) {


  if (current_ion == 'all'){

      current_ion <- 'sum of ions'
  }


  plot_title <- paste(input_peptide_sequence, current_ion, sep='\n')

  # In order to avoid log(0), replace 0 with NA, so that the points with y=0 will not be diplayed.
  impute_status <- FALSE
  if (any(plot_fragment_ion_results$calculated_area_ratio == 0)) {
    impute_status <- TRUE
    ids_tmp <- which(plot_fragment_ion_results$calculated_area_ratio == 0)
    #replace_value <- min(plot_fragment_ion_results$calculated_area_ratio[plot_fragment_ion_results$calculated_area_ratio != 0])/10000
    #plot_fragment_ion_results$calculated_area_ratio[ids_tmp] <- NA
    impute_value <- 
      min(plot_fragment_ion_results$calculated_area_ratio[
        plot_fragment_ion_results$calculated_area_ratio != 0])/100
    plot_fragment_ion_results$calculated_area_ratio[ids_tmp] <- impute_value
  }
  
  min_calculated_area_ratio <- min(plot_fragment_ion_results$calculated_area_ratio, na.rm = TRUE)
  max_calculated_area_ratio <- max(plot_fragment_ion_results$calculated_area_ratio, na.rm = TRUE)
  
  # Expand right side of clipping rect to make room for the legend
  par(xpd=TRUE, mar=par()$mar+c(0,0,0,4))

  # bty="L",

  ##get minimum and maximum calculated area ratios
  min_area_ratio <- min(plot_fragment_ion_results$calculated_area_ratio, na.rm = TRUE)
  max_area_ratio <- max(plot_fragment_ion_results$calculated_area_ratio, na.rm = TRUE)
  
  ##get symbols correponding to different replicates
  symbols_rep <- case_when(plot_fragment_ion_results$replicate==1 ~ 1,
                           plot_fragment_ion_results$replicate==2 ~ 0,
                           plot_fragment_ion_results$replicate==3 ~ 2,
                           plot_fragment_ion_results$replicate==4 ~ 5,
                           plot_fragment_ion_results$replicate==5 ~ 6,
                           plot_fragment_ion_results$replicate==6 ~ 3,
                           plot_fragment_ion_results$replicate==7 ~ 4,
                           plot_fragment_ion_results$replicate==8 ~ 16,
                           plot_fragment_ion_results$replicate==9 ~ 9,
                           plot_fragment_ion_results$replicate==10 ~ 17,
                           TRUE ~ 18)
  
  suppressWarnings(plot(plot_fragment_ion_results$day,
                        plot_fragment_ion_results$calculated_area_ratio,
                        log="y", yaxt="n", 
                        pch=symbols_rep,
                        cex=2, lwd=1, main=plot_title,
                        xlab="Time (day)", ylab="Measured (area ratio) [log-scale]", 
                        cex.lab=1.75, cex.axis=1.75,
                        ylim = c(min_area_ratio*0.8,max_area_ratio*1.2)))

  ##highlight points with the imputed values
  if (impute_status) {
      points(plot_fragment_ion_results$day[ids_tmp], 
             plot_fragment_ion_results$calculated_area_ratio[ids_tmp],
             pch = symbols_rep[ids_tmp],
             col="red",
             bg="red",
             cex=2,
             lwd=2)
      segments(x0 = min(plot_fragment_ion_results$day),
               x1 = max(plot_fragment_ion_results$day),
               y0 = impute_value,
               y1 = impute_value,
               col="red",
               lwd=1.5,
               lty=2)
      text(x = 1.8,
           y = impute_value*1.2,
           labels = "Imputed 0 values",
           col = "red")
  }
  x_axis_values <- days

  y_axis_values <- c(format(min_calculated_area_ratio,digits=3),format(median(plot_fragment_ion_results$calculated_area_ratio),digits=3),format(max_calculated_area_ratio,digits=3))

  suppressWarnings(axis(2, y_axis_values, labels=format(y_axis_values,scientific=FALSE))) # draw y axis with required labels

  #par(xpd=TRUE)
  #legend(x=4,y=max_calculated_area_ratio,legend=c("rep1","rep2","rep3","repX","Hi","Med","Lo"),pch=c(1,0,2,4,18,18,18), col=c("black","black","black","black","red","blue","green"), bty="n")
  
  legend_replicate <- paste("rep", unique(sort(plot_fragment_ion_results$replicate_number)), sep = "")
  pch_repicate_tmp <- unique(sort(plot_fragment_ion_results$replicate_number))
  col_repicate <- rep("black", length(pch_repicate_tmp))
  legend(x=max(days)+0.2, y=max_calculated_area_ratio,
        legend=legend_replicate, pch=symbols_rep, col=col_repicate, cex=1, bty="n")

  # Restore default clipping rect
  par(mar=c(5, 4, 4, 2) + 0.1)
} 

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

args <- commandArgs(trailingOnly = TRUE)
dataset_path <- args[1]
fileList_path <- args[2]
plot_output <- args[3]
plot_output_dir <- args[4]

#dataset_path <- "normal_data.tsv"
#fileList_path <- "file_namelist_IS.tsv"
#plot_output <- "True"
#plot_output_dir <- "D:\\Skyline_analysis\\qcAssayPortal\\qcAssayPortal\\src\\qcAssayPortal\\rScripts\\test\\improve_exp5\\tmp"

if (plot_output == 'True') {
    plot_output <- TRUE
} else {
    plot_output <- FALSE
}

cv_threshold <- 20.0

##test it out with these parameters
# dataset_path <- "CPTAC_TEST/Broad_Carr_CellLysate_TSQQuantiva_IMACMRM_TEST"
# input_protein_name <- "AURKA"
# input_peptide_sequence <- "TT[+80.0]LC[+57.0]GTLDYLPPEMIEGR"
# input_precursor_charge <- 3
# curve_type <- "forward"

# Load data from local table
QC_set_total <- read.table(file=dataset_path, header=TRUE, sep='\t')
fileDf <- read.table(file=fileList_path, header=TRUE, sep='\t')

#QC_set <- labkey.selectRows(
#  baseUrl="https://panoramaweb.org/labkey",
#  folderPath=paste("/CPTAC Assay Portal/",dataset_path,"/ValidationSamples",sep=""),
#  schemaName="targetedms",
#  queryName="QCAnalysisQuery",
#  viewName="",
#  colFilter=makeFilter(c("Protein", "EQUAL", input_protein_name),
#                       c("PeptideModifiedSequence", "EQUAL", input_peptide_sequence),
#                       c("PrecursorCharge","EQUAL",input_precursor_charge)),
#  containerFilter=NULL
#)

# Transform the column names to match those from embedded panorama query

colNumber <- ncol(QC_set_total)
thenames = tolower(names(QC_set_total))
thenames = gsub(" ","", thenames) 
names(QC_set_total) <- thenames

# sample row
#1 YARS.IPI00007074 VDAQFGGIDQR heavy 2 1 y6 3573011.0 GO_QCorig_Broad_1000ng_Interlab_092412_031 3 2 Med 1.0 2 22199 #66fba526-16af-1031-a003-


# rename columns in QC_set_total dataframe (replace Panorama names with new names used by R script)
QC_set_total <- dplyr::rename(QC_set_total,
                        SkyDocumentName = skydocumentname,
                        peptide = peptidemodifiedsequence,
                        protein_name = proteinname,
                        precursor_charge = precursorcharge,
                        product_charge = productcharge,
                        fragment_ion_only = fragmention,
                        day = day,
                        replicate_number = replicatenumber,
                        replicate_name = replicatename,
                        sample_group= samplegroup,
                        isotope_label_type = isotopelabeltype,
                        area = area)

if (nrow(QC_set_total) ==0) {
    QC_set_total$fragment_ion <- integer(0)
} else {
    QC_set_total$fragment_ion <- paste(QC_set_total[ ,'fragment_ion_only'], " (", QC_set_total[ ,'product_charge'], "+)", sep='' )
}

ion_category <- 'error'

# convert columns from character to numeric
QC_set_total[,'day'] <- as.numeric(as.character(QC_set_total[,'day']))
QC_set_total[,'replicate_number'] <- as.numeric(as.character(QC_set_total[,'replicate_number']))
QC_set_total[,'area'] <- as.numeric(as.character(QC_set_total[,'area']))
# remove factor version
QC_set_total[,'fragment_ion'] <- as.character(QC_set_total[,'fragment_ion'])
QC_set_total[,'isotope_label_type'] <- as.character(QC_set_total[,'isotope_label_type'])

# Write peptide information into output file.
log_filename <- paste(plot_output_dir, "\\peptide_infor.tsv", sep='' )
logdf <- data.frame(peptide=as.character(), precursorCharge=as.character(), isotopeLabelType=as.character(), transition=as.character(), uniProtKBID=as.character(), proteinName=as.character(), SkyDocumentName=as.character())

#########################################
# Separate the error detecting codes from the warning detecting codes.
# Traverse the SkyDocumentName in fileDf to detect all the possible errors.
# The details will be added later.
#########################################
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
            protein_list <- as.character(unique(QC_set3[ , 'protein_name']))
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
                # Get a list of all unique fragment ions, unique datas, unique replicates, unique isotope_label_type associated with current peptide.
                fragment_ion_list <- unique(QC_setTmp[ , 'fragment_ion'])
                days <- sort(unique(QC_setTmp[ , 'day']))
                replicates <- sort(unique(QC_setTmp[ , 'replicate_number']))
                isotope_label_types <- unique(QC_setTmp[ , 'isotope_label_type'])
                
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
                
                
                # Make judgement whether there are multiple heavy or light area for the combination of fragment_ion, replicate, day.
                # If it happens, traverse fragment_ion_list, days, sample_groups and replicates to evaluate the fragment_ion under the specific combination of day, sample_group, and replicate.
                # The reason to this error is that the annotation of column 
                evaOut1 <- evaluate("dcast(QC_setTmp, protein_name + peptide + precursor_charge + fragment_ion + replicate + day ~ isotope_label_type, value.var='area')")
                evaOut2 <- evaluate("dcast(QC_setTmp, protein_name + peptide + precursor_charge + fragment_ion + replicate_name + day ~ isotope_label_type, value.var='area')")
                
                if (length(evaOut1) == 3) {
                    # In this condition, some replicate information is wrong for some combinations of fragment_ion, day and sample_group.
                    # The wrongly annotated replicate need to be generated.
                    df1 <- suppressMessages(dcast(QC_setTmp, protein_name + peptide + precursor_charge + fragment_ion + replicate + day ~ isotope_label_type, value.var='area'))
                    # Evaluate the fragment_ion under the specific combination of day, sample_group, and replicate.
                    # df1 can be used to extract the combinations
                    errorReasonTmp <- c()
                    for (index in 1:nrow(df1)) {
                        current_set <- QC_setTmp[QC_setTmp$fragment_ion==df1[index, ]$fragment_ion & QC_setTmp$day==df1[index, ]$day & QC_setTmp$replicate==df1[index, ]$replicate, ]
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
                                errorReason_item <- paste(errorReason_item1, ' due to wrongly annotated values in attributes: replicate or day', sep='')
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
    df_skydoc_error_peptide[[SkyDocumentName]] <- peptide_list_with_error
}

if (plot_output) {
    write.table(logdf, file=log_filename, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
}

df_skydoc_error_peptide_precursorCharge$SkyDocumentName <-  as.character(df_skydoc_error_peptide_precursorCharge$SkyDocumentName)
df_skydoc_error_peptide_precursorCharge$protein_name <-  as.character(df_skydoc_error_peptide_precursorCharge$protein_name)
df_skydoc_error_peptide_precursorCharge$peptide <-  as.character(df_skydoc_error_peptide_precursorCharge$peptide)
df_skydoc_error_peptide_precursorCharge$precursorCharge <- as.character(df_skydoc_error_peptide_precursorCharge$precursorCharge)

#########################################
# Infer the internal standard type for each SkyDocumentName by randomly sampled 5 peptides.
# Since the internal standard type has to be predefined to be heavy in experiment, there is no need to infer the internal standard type.
# For Exp5, the internal standard type can't be inferred, so in df_internal_standard_inferred every inferred internal_stand will keep the original setting.
#########################################
df_internal_standard_inferred <- data.frame(SkyDocumentName=as.character(), internal_standard=as.character())
for (SkyDocumentName in as.character(fileDf[, "SkyDocumentName"])) {
    internal_standard_inferred <- as.character(fileDf[fileDf$SkyDocumentName == SkyDocumentName, ]$internal_standard)[1]
    df_internal_standard_inferred_tmp <- data.frame(SkyDocumentName=SkyDocumentName, internal_standard=internal_standard_inferred)
    df_internal_standard_inferred <- rbind(df_internal_standard_inferred, df_internal_standard_inferred_tmp)
}
is_inferred_filename <- paste(plot_output_dir, "\\internal_standard_inferred_infor.tsv", sep='' )
if (plot_output) {
    write.table(df_internal_standard_inferred, file=is_inferred_filename, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
}

#########################################
# Traverse the SkyDocumentName in fileDf to detect all the possible warnings and generate the images and the tables for peptide without issues.
# Detecting warnings functions will be added later. So generating images an tables for peptides without issues can be implemented firstly.
#########################################
for (SkyDocumentName in as.character(fileDf[, "SkyDocumentName"])) {
    # The defualt internal_standard should always be heavy. Evaluate the internal_standard, if the internal standard is wrong, errors will arise.
    original_internal_standard <- as.character(fileDf[fileDf$SkyDocumentName == SkyDocumentName, ]$internal_standard)
    inferred_internal_standard <- as.character(df_internal_standard_inferred[df_internal_standard_inferred$SkyDocumentName == SkyDocumentName, ]$internal_standard)
    if (original_internal_standard[1] == 'none') {
        # Just jump out of the loop. Don't print the errorInfor, because it has already be printed in the function of detectIS in qcAnalysis.py
        next
    }
    
    # This won't happen.
    if (original_internal_standard[1] != inferred_internal_standard[1]) {
        errorType <- "Error"
        errorSubtype <- "Internal standard"
        errorReason <- paste('The internal standard in the skyline file is set to be ', original_internal_standard, ' which is incorrect, please set the Internal standard type in the peptide_modifications underneath peptide_settings to be heavy.', sep='')
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
    # In fact, curve_type is always 'forward'.
    
    for (input_peptide_sequence in peptide_list) {
        QC_set_2 <- QC_set_1[QC_set_1$peptide==input_peptide_sequence, ]
        precursor_charge_list <- unique(QC_set_2[ , 'precursor_charge'])
        for (input_precursor_charge in precursor_charge_list) {
            QC_set_3 <- QC_set_2[QC_set_2$precursor_charge==input_precursor_charge, ]
            protein_list <- as.character(unique(QC_set_3[ , 'protein_name']))
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
                replicates <- sort(unique(QC_set[ , 'replicate_number']))
                isotope_label_types <- unique(QC_set[ , 'isotope_label_type'])
                QC_set$isotope_label_type[QC_set$isotope_label_type == "medium"] <- "light"
                
                if (length(fragment_ion_list) < 3) {
                    errorType <- "Warning"
                    errorSubtype <- "Fragment ion"
                    errorReason <- paste("In reproducible detection of endogenous analyte graph, the number of fragment ions is ",  length(fragment_ion_list), " < 3, the fragment ions is(are): ", paste(fragment_ion_list, collapse = ', '), '.', sep="")
                    errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge), rep('', colNumber-5)), collapse='\t')
                    cat(errorInfor)
                    cat('\n')
                }
                
                fragment_ion_results <- data.frame()
                # # ***** prepare file to print PNG images *****
                # only plot top 3 ions plus one more for sum
                image_frame_count <- 4
                CairoPNG(filename=paste(plot_output_dir, "\\", input_peptide_sequence, "_", input_precursor_charge, "_", curve_type, ".png", sep=""), width=image_frame_count*400, height=400, bg="white", units="px")
                par(mfrow= c(1, image_frame_count))
                
                # Step 1:
                # Traverse fragment_ion_list, days, sample_groups and replicates
                # to evaluate the fragment_ion under the specific combination of day and replicate_number.

                for (current_ion in fragment_ion_list) {
                    for (current_day in days) {
                        for (current_rep in replicates) {
                            current_set_count <- 0
                            light_area <- 0
                            heavy_area <- 0
                            theoretical_area <- 0
                            measured_area <- 0
                            calculated_area_ratio <- 0
                            
                            current_set <- QC_set[QC_set$fragment_ion==current_ion & QC_set$day==current_day & QC_set$replicate_number==current_rep, ]
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
                                
                                if(theoretical_area==0 | is.na(theoretical_area) | is.na(measured_area)){
                                    calculated_area_ratio = NA
                                }
                                else {
                                    calculated_area_ratio <- measured_area/theoretical_area
                                }
                                fragment_ion_results <-  rbind(fragment_ion_results, data.frame(Peptide = input_peptide_sequence, Protein_Name = input_protein_name, Precursor_Charge = input_precursor_charge,
                                                                                                  fragment_ion = current_ion, day = current_day, replicate_number = current_rep, light_area=light_area, heavy_area=heavy_area,
                                                                                                  theoretical_area=theoretical_area, measured_area=measured_area,
                                                                                                  calculated_area_ratio=calculated_area_ratio, ion_category='individual') )
                            }
                            else {
                            }
                        } # end current_rep
                    } # end current_day
                } # end current_ion
                
                # Step 2: repeat calculations for sum of ions 
                for (current_day in days) {
                    for (current_rep in replicates) {
                        sum_light_area <- 0
                        sum_heavy_area <- 0
                        skipped_count <- 0
                        skip_current_sample <- 'true'
                        
                        sum_theoretical_area <- 0
                        sum_measured_area <- 0
                        for (current_ion in fragment_ion_list)  {
                            current_set_count <- 0
                            calculated_area_ratio <- 0
                            light_area <- 0
                            heavy_area <- 0
                            theoretical_area <- 0
                            measured_area <- 0
                            
                            current_set <- QC_set[QC_set$fragment_ion==current_ion & QC_set$day==current_day & QC_set$replicate_number==current_rep, ]

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
                        } # end current_ion
                        
                        if (skip_current_sample=='false') {
                            if(sum_theoretical_area==0){
                                calculated_area_ratio <- NA
                            }
                            else {
                                calculated_area_ratio <- sum_measured_area/sum_theoretical_area
                            }
                            fragment_ion_results <-  rbind(fragment_ion_results, data.frame(Peptide = input_peptide_sequence, Protein_Name = input_protein_name, Precursor_Charge = input_precursor_charge,
                                                                                              fragment_ion = 'all', day = current_day, replicate_number = current_rep, light_area=sum_light_area, heavy_area=sum_heavy_area,
                                                                                              theoretical_area=sum_theoretical_area, measured_area=sum_measured_area,
                                                                                              calculated_area_ratio=calculated_area_ratio, ion_category='all') )
                        }
                    } # end current_rep
                } # end current_day
                # Step 3:
                # ***** calculate CV (Coefficient of Variation) *****
                ions <- c(fragment_ion_list, 'all')
                # make CV summary data frame
                CV_results <- data.frame(fragment_ion = ions,
                                         intra_CV=NA,
                                         inter_CV=NA,
                                         total_CV=NA,
                                         total_count=NA)
                # *** intra-assay CV ***
                for (current_ion in ions) {
                    avg_intra_assay_CV <-0
                    individual_intra_assay_CVs <- c()
                    for (current_day in days) {
                        current_set <- fragment_ion_results[fragment_ion_results$fragment_ion==current_ion & 
                                            fragment_ion_results$day==current_day, ]
                        # remove rows with a value of NA for calculated_area_ratio
                        current_set <- current_set[complete.cases(current_set[ , 'calculated_area_ratio' ]),]
                        if (nrow(current_set) <= 1 ){
                            percent_CV <- NA
                        }
                        else {
                            percent_CV <- (sd(current_set$calculated_area_ratio))/(mean(current_set$calculated_area_ratio)) * 100
                            individual_intra_assay_CVs <- c(individual_intra_assay_CVs, percent_CV)
                        } 
                    } # end current_day
                    
                    if (length(individual_intra_assay_CVs)==0){
                        avg_CV <- NA
                        count <- 0
                    }
                    else {
                        avg_CV <- mean(individual_intra_assay_CVs, na.rm = TRUE)
                    }
                    CV_results[CV_results$fragment_ion==current_ion, "intra_CV"] <- round(avg_CV, digits=1)
                } # end current_ion
                
                # *** END: intra-assay CV ***
                for (current_ion in ions) {
                    avg_inter_assay_CV <-0
                    individual_inter_assay_CVs <- c()
                    avg_CV  <- 0
                    for (current_rep in replicates) {
                        current_set <- fragment_ion_results[fragment_ion_results$fragment_ion==current_ion & 
                                           fragment_ion_results$replicate_number==current_rep, ]
                        # remove rows with a value of NA for calculated_area_ratio
                        current_set <- current_set[complete.cases(current_set[ , 'calculated_area_ratio' ]),]
                        if (nrow(current_set) <= 1 ){
                            percent_CV <- NA
                        }
                        else {
                            percent_CV <- (sd(current_set$calculated_area_ratio))/(mean(current_set$calculated_area_ratio)) * 100
                            individual_inter_assay_CVs <- c(individual_inter_assay_CVs, percent_CV)
                        }
                    } # end current_rep
                    if (length(individual_inter_assay_CVs)==0){
                        avg_CV <- NA
                        count <- 0
                    }
                    else {
                        avg_CV <- mean(individual_inter_assay_CVs, na.rm = TRUE)
                    }
                    CV_results[CV_results$fragment_ion==current_ion, 'inter_CV'] <- round(avg_CV, digits=1)
                } # end current_ion
                # *** END: inter-assay CV ***
                # calculate total variability
                CV_results[ , 'total_CV'] <- round(sqrt((CV_results[ , 'intra_CV'])*(CV_results[ , 'intra_CV']) + 
                                          (CV_results[ , 'inter_CV'])*(CV_results[ , 'inter_CV'])), digits=1)
                # determine counts
                for (current_ion in ions){
                    CV_results[CV_results$fragment_ion==current_ion, "total_count"] <- 
                    nrow(fragment_ion_results[fragment_ion_results$fragment_ion==current_ion & !is.na(fragment_ion_results$calculated_area_ratio), ] )
                }
                
                ions_to_plot <- c()
                ions_in_table <- c()
                
                # determine fragment ions to plot
                if (length(ions) <= 4 ) {
                    results_to_plot <- CV_results[CV_results$fragment_ion!='all' & !is.na(CV_results$intra_CV) & !is.na(CV_results$inter_CV) & !is.na(CV_results$total_CV) & !is.na(CV_results$total_count), ]
                    ions_list <- as.character(results_to_plot[ , 'fragment_ion'])
                    results_to_plot2 <- CV_results[CV_results$fragment_ion=='all' & !is.na(CV_results$intra_CV) & !is.na(CV_results$inter_CV) & !is.na(CV_results$total_CV) & !is.na(CV_results$total_count), ]
                    ions_to_plot <- c(ions_list, as.character(results_to_plot2[ , 'fragment_ion']))
                    ions_to_plot <- ions
                } else {
                    results_to_plot <- CV_results[CV_results$fragment_ion!='all' & !is.na(CV_results$intra_CV) & !is.na(CV_results$inter_CV) & !is.na(CV_results$total_CV) & !is.na(CV_results$total_count), ]
                    ions_list <- as.character(results_to_plot[ , 'fragment_ion'])
                    # new sort to get Top 3 plots
                    results_to_plot <- results_to_plot[order(results_to_plot$total_CV, results_to_plot$intra_CV, results_to_plot$inter_CV), ]
                    three_lowest_total_CV <- head(results_to_plot, 3)
                    three_lowest_ions <- as.character(three_lowest_total_CV[ , 'fragment_ion'])
                    results_to_plot2 <- CV_results[CV_results$fragment_ion=='all' & !is.na(CV_results$intra_CV) & !is.na(CV_results$inter_CV) & !is.na(CV_results$total_CV) & !is.na(CV_results$total_count), ]
                    ions_to_plot  <- c(three_lowest_ions, as.character(results_to_plot2[ , 'fragment_ion']))
                    ions_in_table <- c(ions_list, as.character(results_to_plot2[ , 'fragment_ion']))
                    #results_to_plot <- CV_results[CV_results$fragment_ion!='all' & !is.na(CV_results$total_CV), ]
                    # new sort to get Top 3 plots
                    #results_to_plot <- results_to_plot[order(results_to_plot$total_CV), ]
                    #three_lowest_total_CV <- head(results_to_plot, 3)
                    #three_lowest_ions <- as.character(three_lowest_total_CV[ , 'fragment_ion'])
                    #ions_to_plot  <- c(three_lowest_ions, 'all')
                }
                
                # Check whether ions_to_plot is blank or now. If it's blank, an error message will be thrown.
                if (length(ions_to_plot) == 0) {
                    errorType <- "Error"
                    errorSubtype <- "Missing points"
                    errorReason <- "In Reproducible detection of endogenous analyte graph, there is no fragment ions and their CVs can't be calculated."
                    errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge), rep('', colNumber-5)), collapse='\t')
                    cat(errorInfor)
                    cat('\n')
                    next
                }

                par(mfrow = c(1,length(ions_to_plot)))
                
                errorReason_missing_point_on_day <- c()
                errorReason_missing_replicate_point_on_day <- c()
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
                    }
                    
                    # Transverse each day and count the number of replicates, if the number is less than 3, a warning will arise.
                    day_with_less_replicates <- c()
                    for (plot_day in plot_days) {
                        plot_fragment_ion_results_tmp1 <- plot_fragment_ion_results[plot_fragment_ion_results$day == plot_day, ]
                      
                        #concentration_list <- unique(plot_fragment_ion_results_tmp1$sample_group)
                        #concentrationsTmp <- setdiff(sample_groups,concentration_list) 
                        #if (length(concentrationsTmp) > 0) {
                            #errorReason_tmp <- paste("For fragment ion ", current_plot_ion, ", not all of three concentrations ('Hi', 'Med', 'Lo') exist on day ",  plot_day, sep="")
                            #errorReason_missing_concerntration_on_day <- c(errorReason_missing_concerntration_on_day, errorReason_tmp)
                        #}
                        replicate_count <- length(unique(plot_fragment_ion_results_tmp1$replicate_number))
                        if (replicate_count < 3) {
                            day_with_less_replicates <- c(day_with_less_replicates, plot_day)
                        }
                    }
                    if (length(day_with_less_replicates) > 0) {
                        errorReason_tmp <- paste("For fragment ion ", current_plot_ion, ", there are less than 3 replicates on day(s): ",  paste(day_with_less_replicates, collapse=', '), sep="")
                        errorReason_missing_replicate_point_on_day <- c(errorReason_missing_replicate_point_on_day, errorReason_tmp)
                    }
                }
                dev.off()
                
                if (length(errorReason_missing_point_on_day) > 0) {
                    errorType <- "Warning"
                    errorSubtype <- "Missing points"
                    errorReason <- paste(paste(errorReason_missing_point_on_day, collapse='. '), '.', sep="")
                    errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge), rep('', colNumber-5)), collapse='\t')
                    cat(errorInfor)
                    cat('\n')
                }
                
                if (length(errorReason_missing_replicate_point_on_day) >0 ) {
                    errorType <- "Warning"
                    errorSubtype <- "Missing points"
                    errorReason <- paste(paste(errorReason_missing_replicate_point_on_day, collapse='. '), '.', sep="")
                    errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge), rep('', colNumber-5)), collapse='\t')
                    cat(errorInfor)
                    cat('\n')
                }
                write.table(CV_results, file=paste(plot_output_dir, "\\", input_peptide_sequence, "_", input_precursor_charge, "_", curve_type ,"_CV_results", ".tsv", sep=""), sep = "\t", qmethod = "double", col.names=TRUE, row.names=FALSE, quote=FALSE)
                
                if (FALSE) {
                    CV_results_sorted <- CV_results[match(ions_in_table , CV_results$fragment_ion), ]
                    CV_results_sub_2 <- subset(CV_results_sorted, select = c(fragment_ion, intra_CV, inter_CV, total_CV))
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
                            items_abnormal <- c()
                            for (index1 in index1_array) {
                                strsplit_tmp <- strsplit(colnames(CV_results_sub_2_plot)[index1], split='_')[[1]]
                                cv_type_term <- strsplit_tmp[1]
                                if (cv_type_term == "intra") {
                                    items_abnormal <- c(items_abnormal, 'intra_assay CV')
                                } else if (cv_type_term == "inter") {
                                    items_abnormal <- c(items_abnormal, 'inter_assay CV')
                                } else {
                                    items_abnormal <- c(items_abnormal, 'total CV')
                                }
                            }
                            if (length(items_abnormal) > 0) {
                                errorReason_message <- paste('For fragment ion ', fragment_ion_selected, ', ', paste(items_abnormal, collapse=', '), ' is(are) larger than the threshold of ', cv_threshold, '%.', sep='')
                                errorReason_tmp <- c(errorReason_tmp, errorReason_message)
                            }
                        }
                        errorType <- "Warning"
                        errorSubtype <- "Bad distribution of points"
                        errorReason <- paste(errorReason_tmp, collapse =" ")
                        errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge), rep('', colNumber-5)), collapse='\t')
                        cat(errorInfor)
                        cat('\n')
                    }
                }
            }
        }
    }
}