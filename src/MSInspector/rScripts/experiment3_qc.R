# This sample code is used to qc skyline document data from experiment 3 of Assay Portal.
# It's modified based on:
#   esac-panorama-master\experiment-3\code\Experiment_3_with_labkey_connection_updated_only_plot_4_ions.R

suppressWarnings(suppressMessages(library(Cairo)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(evaluate)))
suppressWarnings(suppressMessages(require(reshape2)))

##variable that specifies whether the script takes a standalone csv file vs going through Panorama
standalone <- TRUE
options(warn=-1)

# Multiple plot function - from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  suppressWarnings(suppressMessages(require(grid)))

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# ***** plot_QC function *****
plot_QC <- function(plot_fragment_ion_results, input_peptide_sequence, current_ion) {
  if (current_ion == 'all'){
    current_ion <- 'sum of ions'
  }

  plot_title <- paste(input_peptide_sequence, current_ion, sep='\n')

  ##average results for each replicate for the same cell line and spike level
  plot_grouped_reps <- plot_fragment_ion_results %>% group_by(cell_line, spike_level)
  plot_ave_reps <- plot_grouped_reps %>%
    summarize(calculated_area_ratio_ave_reps = mean(calculated_area_ratio))

  ##get minimum and maximum calculated area ratios
  min_area_ratio <- min(plot_ave_reps$calculated_area_ratio_ave_reps)
  max_area_ratio <- max(plot_ave_reps$calculated_area_ratio_ave_reps)

  plot_ave_reps$horiz_line <- min_area_ratio

  g <- ggplot(plot_ave_reps, aes(x=spike_level, y=calculated_area_ratio_ave_reps,
                                 color=cell_line, group=cell_line)) +
    geom_point(size=2.5) +
    geom_smooth(method="lm",linetype="dashed",se = FALSE) +
    ##geom_line(aes(group=cell_line),linetype="dashed") +
    ggtitle(plot_title) +
    ##coord_trans(y="log10") +
    scale_y_continuous(limits = c(min_area_ratio*0.8, max_area_ratio*1.2)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    guides(color = guide_legend(title="Cell line")) +
    xlab("Spike-in level") + ylab("Measured (area ratio) [linear scale]")
  g

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

cv_cal <- function(x) {
    sd(x)/mean(x)
}

cv_cal2 <- function(x) {
    x_min <- min(x)
    if (x_min < 0) {
       x <- x + abs(x_min)
    }
    sd(x)/mean(x)
}

cv_cal3 <- function(x) {
    x_tmp <- abs(x)
    sd(x_tmp)/mean(x_tmp)
}

intercept_scope <- function(x) {
    max(x)-min(x)
}

rmse <- function(x) {
    sqrt(sum(x**2)/length(x))
}

locate_abnormal_value <- function(x, valueType, thresold_value) {
    # x is a data.frame with the format of:
    # fragment_ion          1          2          3          4           5           6
    # y8 (1+) 0.06990185 0.10183874 0.01727739 0.07365045 0.052791834 0.026315809
    # y6 (1+) 0.08755183 0.05586970 0.07590068 0.05780173 0.012987405 0.004313037
    # valueType is rSquare or pValue
    out_list <- list()
    # Traverse each row to make a comparison except the first column.
    for (i in 1:nrow(x)) {
        for (k in 2:ncol(x)) {
            if (valueType == "rSquare") {
                if (x[i, k] < thresold_value) {
                    out_list <- c(out_list, list(c(x[i,1], colnames(x)[k])))
                }
            } else {
                if (x[i, k] > thresold_value) {
                    out_list <- c(out_list, list(c(x[i,1], colnames(x)[k])))
                }
            }
        }
    }
    out_list
}

args <- commandArgs(trailingOnly = TRUE)
dataset_path <- args[1]
fileList_path <- args[2]
plot_output <- args[3]
plot_output_dir <- args[4]

#dataset_path <- "normal_data.tsv"
#fileList_path <- "file_namelist_IS.tsv"
#plot_output <- "True"
#plot_output_dir <- "D:\\Skyline_analysis\\qcAssayPortal\\qcAssayPortal\\src\\qcAssayPortal\\rScripts\\test\\transition_ratio_exp3\\tmp"

if (plot_output == 'True') {
    plot_output <- TRUE
} else {
    plot_output <- FALSE
}

cv_threshold <- 0.5
rmse_threshold <- 1.0
pValue_threshold <- 0.05
rSquare_threshold <- 0.5
relative_difference_threshold <- 0.2
interceptSwtich <- FALSE
nonBlankConcentrationSwitch <- TRUE


# Load data from local table
QC_set_total <- read.table(file=dataset_path, header=TRUE, sep='\t')
fileDf <- read.table(file=fileList_path, header=TRUE, sep='\t')

#QC_set <- labkey.selectRows(
#    baseUrl="http://cptac-proliant-linux.esacinc.com/labkey",
#    folderPath=paste("/CPTAC Assay Portal/",dataset_path,"/Selectivity",sep=""),
#    schemaName="targetedms",
#    queryName="QCAnalysisQuery",
#    viewName="",
#    colFilter=makeFilter(c("Protein", "EQUAL", input_protein_name),
#                         c("PeptideModifiedSequence", "EQUAL", input_peptide_sequence),
#                         c("PrecursorCharge","EQUAL",input_precursor_charge)),
#    containerFilter=NULL
#)

#current_number = 1234;

# Transform the column names to match those from embedded panorama query
colNumber <- ncol(QC_set_total)
thenames <- tolower(names(QC_set_total))
thenames <- gsub(" ","", thenames)
names(QC_set_total) <- thenames

# sample row
#1 YARS.IPI00007074 VDAQFGGIDQR heavy 2 1 y6 3573011.0 GO_QCorig_Broad_1000ng_Interlab_092412_031 3 2 Med 1.0 2 22199 #66fba526-16af-1031-a003-


# rename columns in QC_set_total dataframe (replace Panorama names with new names used by R script)
colnames(QC_set_total)[colnames(QC_set_total)=="skydocumentname"] <- "SkyDocumentName"
colnames(QC_set_total)[colnames(QC_set_total)=="peptidemodifiedsequence"] <- "peptide"
colnames(QC_set_total)[colnames(QC_set_total)=="proteinname"] <- "protein_name"
colnames(QC_set_total)[colnames(QC_set_total)=="replicatename"] <- "replicate_name"
colnames(QC_set_total)[colnames(QC_set_total)=="precursorcharge"] <- "precursor_charge"
colnames(QC_set_total)[colnames(QC_set_total)=="productcharge"] <- "product_charge"
colnames(QC_set_total)[colnames(QC_set_total)=="fragmention"] <- "fragment_ion_only"
colnames(QC_set_total)[colnames(QC_set_total)=="analyteconcentration"] <- "analyte_concentration"      # three concentration levels (no spike and two analyte spikes)
colnames(QC_set_total)[colnames(QC_set_total)=="replicatenumber"] <- "replicate" # recplicate number
colnames(QC_set_total)[colnames(QC_set_total)=="samplegroup"] <- "sample_group"     # individual samples of the matrix of interest
colnames(QC_set_total)[colnames(QC_set_total)=="isotopelabeltype"] <- "isotope_label_type"     # light, heavy
colnames(QC_set_total)[colnames(QC_set_total)=="area"] <- "area"

if (nrow(QC_set_total) ==0) {
  QC_set_total$fragment_ion <- integer(0)
} else {
  QC_set_total$fragment_ion <- paste(QC_set_total[ ,'fragment_ion_only'], " (", QC_set_total[ ,'product_charge'], "+)", sep='' )
}

ion_category <- 'error'

# convert some classes for some columns
QC_set_total[,'sample_group'] <- as.character(QC_set_total[,'sample_group'])
# remove factor version
QC_set_total[,'fragment_ion'] <- as.character(QC_set_total[,'fragment_ion'])
QC_set_total[,'isotope_label_type'] <- as.character(QC_set_total[,'isotope_label_type'])
QC_set_total[,'precursor_charge'] <- as.character(QC_set_total[,'precursor_charge'])
QC_set_total[,'product_charge'] <- as.character(QC_set_total[,'product_charge'])

# convert columns from character to numeric
QC_set_total[,'replicate'] <- as.numeric(as.character(QC_set_total[,'replicate']))
QC_set_total[,'area'] <- as.numeric(as.character(QC_set_total[,'area']))
QC_set_total[,'analyte_concentration'] <- as.numeric(as.character(QC_set_total[,'analyte_concentration']))

# Write peptide information into output file.
log_filename <- paste(plot_output_dir, "\\peptide_infor.tsv", sep='' )
logdf <- data.frame(peptide=as.character(), precursorCharge=as.character(), isotopeLabelType=as.character(), transition=as.character(), uniProtKBID=as.character(), proteinName=as.character(), SkyDocumentName=as.character())

#########################################
# Separate the error detecting codes from the warning detecting codes.
# Traverse the SkyDocumentName in fileDf to detect all the possible errors.
# Create a list to store the  peptides with errors for each SkyDocumentName.
# df_skydoc_error_peptide stores the peptides with errors. But there my be one special situation that one specific precursor charge of the peptide has errors but the others don't have error.
# Therefore, df_skydoc_error_peptide_precursorCharge is used to store peptide plus precursor charge.
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
                # Get a list of all unique fragment ions, unique sample groups (cell lines), unique analyte concentrations, unique replicates and unique isotope_label_type associated with current peptide
                fragment_ion_list <- unique(QC_setTmp[ , 'fragment_ion'])
                cell_lines <- sort(unique(QC_setTmp[ , 'sample_group']))
                spike_levels <- sort(unique(QC_setTmp[ , 'analyte_concentration']))
                replicates <- sort(unique(QC_setTmp[ , 'replicate']))
                isotope_label_types <- unique(QC_setTmp[ , 'isotope_label_type'])
                # Check medium labeled peptides.
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
                # Check the analyte concentrations.
                if(length(spike_levels) <2) {
                    errorType <- "Error"
                    errorSubtype <- "Concentration"
                    errorReason <- paste("More than one analyte concentration levels are needed: ", paste(spike_levels, collapse='; '), ".", sep='')
                    errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge), rep('', colNumber-5)), collapse='\t')
                    cat(errorInfor)
                    cat('\n')
                    peptide_list_with_error <- c(peptide_list_with_error, input_peptide_sequence)
                    peptide_list_with_error2 <- c(peptide_list_with_error2, input_peptide_sequence)
                    next
                }
                # Make judgement whether there are multiple heavy or light area for the combination of fragment_ion, replicate, analyte concentration and sample_group.
                # If it happens, traverse fragment_ion_list, analyte concentrations, sample_groups and replicates to evaluate the fragment_ion under the specific combination of analyte concentration, sample_group, and replicate.
                # The reason to this error is that the annotation of column
                evaOut1 <- evaluate("dcast(QC_setTmp, protein_name + peptide + precursor_charge + fragment_ion + replicate + analyte_concentration + sample_group ~ isotope_label_type, value.var='area')")
                evaOut2 <- evaluate("dcast(QC_setTmp, protein_name + peptide + precursor_charge + fragment_ion + replicate_name + analyte_concentration + sample_group ~ isotope_label_type, value.var='area')")
                if (length(evaOut1) == 3) {
                    #cat('Error happens.............................................................')
                    # In this condition, some replicate information is wrong for some combinations of fragment_ion, analyte concentration and sample_group.
                    # The wrongly annotated replicate need to be generated.
                    df1 <- suppressMessages(dcast(QC_setTmp, protein_name + peptide + precursor_charge + fragment_ion + replicate + analyte_concentration + sample_group ~ isotope_label_type, value.var='area'))
                    # Evaluate the fragment_ion under the specific combination of analyte_concentration, sample_group, and replicate.
                    # df1 can be used to extract the combinations
                    errorReasonTmp <- c()
                    for (index in 1:nrow(df1)) {
                        current_set <- QC_setTmp[QC_setTmp$fragment_ion==df1[index, ]$fragment_ion & QC_setTmp$analyte_concentration==df1[index, ]$analyte_concentration & QC_setTmp$sample_group==df1[index, ]$sample_group & QC_setTmp$replicate==df1[index, ]$replicate, ]
                        
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
                                errorReason_item <- paste(errorReason_item1, ' due to wrongly annotated values in attributes: replicate, analyte concentration or sample group', sep='')
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

#########################################
# Infer the internal standard type for each SkyDocumentName by randomly sampled 5 peptides.
#########################################
df_internal_standard_inferred <- data.frame(SkyDocumentName=as.character(), internal_standard=as.character())
for (SkyDocumentName in as.character(fileDf[, "SkyDocumentName"])) {
    QC_set1 <- QC_set_total[QC_set_total$SkyDocumentName==SkyDocumentName, ]
    peptide_list <- unique(QC_set1[ , 'peptide'])
    # Remove the Rows in QC_set1 where the rows in df_skydoc_error_peptide_precursorCharge match.
    #df_skydoc_error_peptide_precursorCharge_toCheck <- df_skydoc_error_peptide_precursorCharge[df_skydoc_error_peptide_precursorCharge$SkyDocumentName==SkyDocumentName, ]
    #ids <- c()
    #for (row_id in 1:nrow(df_skydoc_error_peptide_precursorCharge_toCheck)) {
    #    id_tmp <- which(QC_set1$peptide==df_skydoc_error_peptide_precursorCharge_toCheck[row_id, "peptide"] & QC_set1$precursor_charge==df_skydoc_error_peptide_precursorCharge_toCheck[row_id, "precursorCharge"] & QC_set1$protein_name==df_skydoc_error_peptide_precursorCharge_toCheck[row_id, "protein_name"])
    #    ids <- c(ids, id_tmp)
    #}
    #ids <- unique(ids)
    #QC_set1 <- QC_set1[-ids, ]
    # Randomly sampling 5 peptides with specific precursor_charge and protein_name.
    #QC_set1_tmp <- QC_set1[,c('peptide', 'precursor_charge', 'protein_name')]
    #QC_set1_tmp <- unique(QC_set1_tmp)
    #random_row_id <- c()
    #if (nrow(QC_set1_tmp) >= 5) {
    #    random_row_id <- sample(1:nrow(QC_set1_tmp), 5, replace = FALSE)
    #} else {
    #    random_row_id <- 1:nrow(QC_set1_tmp)
    #}
    #Remove the peptides with errors.
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
            QC_set2 <- QC_set1[QC_set1$peptide==input_peptide_sequence, ]
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
                    # Get a list of all unique fragment ions, unique sample groups (cell lines), unique analyte concentrations, unique replicates and unique isotope_label_type associated with current peptide
                    fragment_ion_list <- unique(QC_setTmp[ , 'fragment_ion'])
                    cell_lines <- sort(unique(QC_setTmp[ , 'sample_group']))
                    spike_levels <- sort(unique(QC_setTmp[ , 'analyte_concentration']))
                    replicates <- sort(unique(QC_setTmp[ , 'replicate']))
                    isotope_label_types <- unique(QC_setTmp[ , 'isotope_label_type'])
                    QC_setTmp$isotope_label_type[QC_setTmp$isotope_label_type == "medium"] <- "light"
                    df <- dcast(QC_setTmp, protein_name + peptide + precursor_charge + fragment_ion + replicate + analyte_concentration + sample_group ~ isotope_label_type, value.var='area')
                    colnames(df) <-  c("ProteinName", "PeptideModifiedSequence", "PrecursorCharge", "FragmentIon", "Replicate", "Concentration", "SampleGroup", "heavyArea", "lightArea")
                    df$lightArea[is.na(df$heavyArea)] <- NA
                    df$heavyArea[is.na(df$lightArea)] <- NA
                    df$heavyArea[df$lightArea==0] <- NA
                    # Because internal_standard_hypothetical is set to be light, HLRatio = heavyArea/lightArea
                    df$HLRatio <- df$heavyArea/df$lightArea
                    # Get the ions that have the highest median area across the entire dataset
                    ions_total <- unique(df$FragmentIon)
                    if (length(ions_total) <= 3 ) {
                        ions_to_plot_tmp <- ions_total
                    } else {
                        median_area_per_ion_tmp <- tapply(df$HLRatio,
                                                  INDEX = df$FragmentIon,
                                                  FUN=median,
                                                  na.rm=TRUE)
                        ions_to_plot_tmp <- names(sort(median_area_per_ion_tmp, decreasing = TRUE))[1:3]
                    }
                    for (fragmentIon_tmp in ions_to_plot_tmp) {
                        thisFragmentIon_tmp <- df[df$FragmentIon == fragmentIon_tmp, ]
                        SampleGroup_list <- unique(thisFragmentIon_tmp$SampleGroup)
                        for (sampleGroup_tmp in SampleGroup_list) {
                            thisFragmentIon <- thisFragmentIon_tmp[thisFragmentIon_tmp$SampleGroup == sampleGroup_tmp, ]
                            minConcentration <- min(thisFragmentIon$Concentration)
                            maxConcentration <- max(thisFragmentIon$Concentration)
                            medHLRation_min <- median(thisFragmentIon[thisFragmentIon$Concentration == minConcentration, ]$HLRatio, na.rm=TRUE)
                            medHLRation_max <- median(thisFragmentIon[thisFragmentIon$Concentration == maxConcentration, ]$HLRatio, na.rm=TRUE)
                            if (is.na(medHLRation_min) | is.na(medHLRation_max)) {
                                invisible()
                            } else if (medHLRation_max >= medHLRation_min) {
                                value1 <- value1 + 1
                            } else {
                              value2 <- value2 + 1
                            }
                        }
                    }
                    #temp <- ddply(df, .(ProteinName, PeptideModifiedSequence, PrecursorCharge, FragmentIon, Replicate, Concentration, SampleGroup), summarize, medianArea=median(HLRatio, na.rm=TRUE))
                    #orderT <- with(temp,  order(ProteinName, PeptideModifiedSequence, -medianArea))
                    # Pay attention: Only top three fragment ions with the largest medianArea are kept.
                    #df <- merge(df, temp[orderT[1:3], ])
                    #result= ddply(df, .(ProteinName, PeptideModifiedSequence, PrecursorCharge, FragmentIon, Concentration, SampleGroup), 
                    #    summarize, Median=median(HLRatio, na.rm= TRUE), Min = min(HLRatio, na.rm=TRUE), Max=max(HLRatio, na.rm=TRUE))
                    #result$Median[is.na(result$Median)] <- 0
                    #curveDataIndex <- with( result,  order(ProteinName, PeptideModifiedSequence, PrecursorCharge,FragmentIon, Concentration, SampleGroup))
                    #thisPeptide <- result[curveDataIndex,]
                    #fragmentIon_list <- unique(thisPeptide$FragmentIon)
                    #for (fragmentIon_tmp in fragmentIon_list) {
                    #    thisFragmentIon_tmp <- thisPeptide[thisPeptide$FragmentIon == fragmentIon_tmp, ]
                    #    SampleGroup_list <- unique(thisFragmentIon_tmp$SampleGroup)
                    #    for (sampleGroup_tmp in SampleGroup_list) {
                    #        thisFragmentIon <- thisFragmentIon_tmp[thisFragmentIon_tmp$SampleGroup == sampleGroup_tmp, ]
                    #        minConcentration <- min(thisFragmentIon$Concentration)
                    #        maxConcentration <- max(thisFragmentIon$Concentration)
                    #        medHLRation_min <- median(thisFragmentIon[thisFragmentIon$Concentration == minConcentration, ]$Median, na.rm=TRUE)
                    #        medHLRation_max <- median(thisFragmentIon[thisFragmentIon$Concentration == maxConcentration, ]$Median, na.rm=TRUE)
                    #        if (is.na(medHLRation_min) | is.na(medHLRation_max)) {
                    #            invisible()
                    #        } else if (medHLRation_max >= medHLRation_min) {
                    #            value1 <- value1 + 1
                    #        } else {
                    #            value2 <- value2 + 1
                    #        }
                    #    }
                    #}
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
is_inferred_filename <- paste(plot_output_dir, "\\internal_standard_inferred_infor.tsv", sep='' )
if (plot_output) {
    write.table(df_internal_standard_inferred, file=is_inferred_filename, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
}

#########################################
# Traverse the SkyDocumentName in fileDf to detect all the possible warnings and generate the images and the tables for peptide without issues.
# Detecting warnings functions will be added later. So generating images an tables for peptides without issues can be implemented firstly.
#########################################
#rmse_list <- c()
for (SkyDocumentName in as.character(fileDf[, "SkyDocumentName"])) {
    # The defualt internal_standard should always be heavy. Evaluate the internal_standard, if the internal standard is wrong, errors will arise.
    original_internal_standard <- as.character(fileDf[fileDf$SkyDocumentName == SkyDocumentName, ]$internal_standard)
    inferred_internal_standard <- as.character(df_internal_standard_inferred[df_internal_standard_inferred$SkyDocumentName == SkyDocumentName, ]$internal_standard)
    if (original_internal_standard[1] == 'none') {
        # Just jump out of the loop. Don't print the errorInfor, because it has already be printed in the function of detectIS in qcAnalysis.py
        next
    }
    
    if (original_internal_standard[1] != inferred_internal_standard[1]) {
        errorType <- "Error"
        errorSubtype <- "Internal standard"
        #errorReason <- paste('The internal standard in the skyline file is set to be ', original_internal_standard, ' which is incorrect, please set the Internal standard type in the peptide_modifications underneath peptide_settings to be heavy.', sep='')
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
                # get a list of all unique fragment ions associated with current peptide
                fragment_ion_list <- unique(QC_set[ , 'fragment_ion'])
                #  get a list of all unique cell lines associated with current peptide
                cell_lines <- sort(unique(QC_set[ , 'sample_group']))
                # get a list of all unique spike levels (0,5,10) associated with current peptide
                spike_levels <- sort(unique(QC_set[ , 'analyte_concentration']))
                # get a list of all unique replicates associated with current peptide
                replicates <- sort(unique(QC_set[ , 'replicate']))
                # for medium labeled peptides
                isotope_label_types <- unique(QC_set[ , 'isotope_label_type'])
                QC_set$isotope_label_type[QC_set$isotope_label_type == "medium"] <- "light"
                if (length(cell_lines) <6) {
                    errorType <- "Warning"
                    errorSubtype <- "Sample group"
                    errorReason <- 'The number of sample group is less than 6.'
                    errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge), rep('', colNumber-5)), collapse='\t')
                    cat(errorInfor)
                    cat('\n')
                }
                sum_light_area <- 0
                sum_heavy_area <- 0
                theoretical_area <- 0
                measured_area <- 0
                fragment_ion_results <- data.frame()
                
                for (current_ion in fragment_ion_list) {
                    for (current_cell_line in cell_lines) {
                        for (current_spike_level in spike_levels) {
                            for (current_rep in replicates) {
                                current_set_count <- 0
                                light_area <- 0
                                heavy_area <- 0
                                theoretical_area <- 0
                                measured_area <- 0
                                calculated_area_ratio <- 0
                                
                                current_set <- QC_set[QC_set$fragment_ion==current_ion & QC_set$sample_group==current_cell_line & QC_set$analyte_concentration==current_spike_level & QC_set$replicate==current_rep, ]
                                current_set_count <- nrow(current_set)
                                
                                if (current_set_count == 2) {
                                    light_area <- current_set[current_set$isotope_label_type=='light', 'area' ]
                                    heavy_area <- current_set[current_set$isotope_label_type=='heavy', 'area' ]
                                    if(curve_type=='forward'){
                                        theoretical_area <- heavy_area
                                        measured_area <- light_area
                                    }
                                    else if (curve_type=='reverse'){
                                        theoretical_area <- light_area
                                        measured_area <- heavy_area
                                    }
                                    else{
                                        stop("invalid curve type")
                                    }
                                    
                                    if(theoretical_area==0 | is.na(theoretical_area) | is.na(measured_area)){
                                        calculated_area_ratio = NA
                                    }
                                    else {
                                        calculated_area_ratio <- measured_area/theoretical_area
                                    }
                                    
                                    fragment_ion_results <-  rbind(fragment_ion_results, data.frame(Peptide = input_peptide_sequence, Protein_Name = input_protein_name, Precursor_Charge = input_precursor_charge,
                                                                                                    fragment_ion = current_ion, cell_line = current_cell_line, spike_level = current_spike_level, replicate = current_rep, light_area=light_area, heavy_area=heavy_area,
                                                                                                    theoretical_area=theoretical_area, measured_area=measured_area,
                                                                                                    calculated_area_ratio=calculated_area_ratio, ion_category='individual'))
                                }
                                else {
                                    invisible()
                                }
                            } # end current_rep
                        } # end current_sample
                    } # end current_cell_line
                } # end current_ion
                
                # ***** repeat calculations for sum of ions *****
                for (current_cell_line in cell_lines) {
                    for (current_spike_level in spike_levels) {
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
                                
                                current_set <- QC_set[QC_set$fragment_ion==current_ion & QC_set$sample_group==current_cell_line & QC_set$analyte_concentration==current_spike_level & QC_set$replicate==current_rep, ]
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
                                    invisible()
                                }
                            } # end current_ion

                            if (skip_current_sample=='false'){
                                if(sum_theoretical_area==0){
                                    calculated_area_ratio <- NA
                                }
                                else {
                                    calculated_area_ratio <- sum_measured_area/sum_theoretical_area
                                }
                                fragment_ion_results <-  rbind(fragment_ion_results, data.frame(Peptide = input_peptide_sequence, Protein_Name = input_protein_name, Precursor_Charge = input_precursor_charge,
                                                                                                fragment_ion = 'all', cell_line = current_cell_line, spike_level = current_spike_level, replicate = current_rep, light_area=sum_light_area, heavy_area=sum_heavy_area,
                                                                                                theoretical_area=sum_theoretical_area, measured_area=sum_measured_area,
                                                                                                calculated_area_ratio=calculated_area_ratio, ion_category='all'))
                            }
                        } # end current_cell_line
                    } # end current_spike_level
                } # end current_rep

                ions <- c(fragment_ion_list, 'all')
                # make summary data tables for slopes
                summary_table <- data.frame(matrix(NA, ncol=length(cell_lines)+1, nrow=length(ions)))
                colnames(summary_table) <- c("fragment_ion", cell_lines)
                summary_table$fragment_ion <- ions
                
                # make summary data tables for intercepts
                summary_table_intercept <- data.frame(matrix(NA, ncol=length(cell_lines)+1, nrow=length(ions)))
                colnames(summary_table_intercept) <- c("fragment_ion", cell_lines)
                summary_table_intercept$fragment_ion <- ions
                
                # make summary data tables for R squares
                # make summary data tables for p values
                summary_table_rSquare <- data.frame(matrix(NA, ncol=length(cell_lines)+1, nrow=length(ions)))
                colnames(summary_table_rSquare) <- c("fragment_ion", cell_lines)
                summary_table_rSquare$fragment_ion <- ions
                summary_table_pValue <- data.frame(matrix(NA, ncol=length(cell_lines)+1, nrow=length(ions)))
                colnames(summary_table_pValue) <- c("fragment_ion", cell_lines)
                summary_table_pValue$fragment_ion <- ions
                # also save the 3 data points in a file
                values_for_spike_levels <- data.frame(fragment_ion=c(),
                                                      cell_line=c(),
                                                      spike_level=c(),
                                                      calculated_rea_ratio_ave_rep=c())
                # do this for each ion
                for(current_ion in ions) {
                    plot_fragment_ion_results <- fragment_ion_results[!is.na(fragment_ion_results$calculated_area_ratio), ]
                    plot_fragment_ion_results <- plot_fragment_ion_results[plot_fragment_ion_results$fragment_ion == current_ion, ]
                    # average results for each replicate for the same cell line and spike level
                    plot_grouped_reps <- plot_fragment_ion_results %>% group_by(cell_line, spike_level)
                    plot_ave_reps <- plot_grouped_reps %>%
                      summarize(calculated_area_ratio_ave_reps = mean(calculated_area_ratio))
                    plot_ave_reps$cell_line <- as.character(plot_ave_reps$cell_line)
                    # for each cell line, get the estimated slope
                    for(current_cell_line in cell_lines) {
                        cell_line_subset <- plot_ave_reps %>% dplyr::filter(cell_line == current_cell_line)
                        mod = lm(calculated_area_ratio_ave_reps ~ spike_level, data = cell_line_subset)
                        summary_table[summary_table$fragment_ion == current_ion,current_cell_line] <- coef(mod)["spike_level"]
                        summary_table_intercept[summary_table_intercept$fragment_ion == current_ion,current_cell_line] <- mod$coefficients[1]
                        summary_table_rSquare[summary_table_rSquare$fragment_ion == current_ion,current_cell_line] <- summary(mod)$r.squared
                        summary_table_pValue[summary_table_pValue$fragment_ion == current_ion,current_cell_line] <- summary(mod)$coefficients[2,4]
                    }
                    v_current_ion <- cbind(fragment_ion = rep(current_ion, nrow(plot_ave_reps)),
                                           plot_ave_reps)
                    values_for_spike_levels <- rbind(values_for_spike_levels,
                                                       v_current_ion)
                }
                
                # save values for spike levels
                values_for_spike_levels$calculated_area_ratio_ave_reps <- format(round(values_for_spike_levels$calculated_area_ratio_ave_reps, digits = 4))
                write.table(values_for_spike_levels, file=paste(plot_output_dir, "\\", input_peptide_sequence, "_", input_precursor_charge, "_", curve_type ,"_ave_values_for_spike_levels", ".tsv", sep=""), sep = "\t", , qmethod = "double", col.names=TRUE, row.names=FALSE, quote=FALSE)
                
                # determine the ions to plot
                ions_to_plot <- ions
                if (length(ions) <= 4 ) {
                    ions_to_plot <- ions
                } else {
                    # get the ions that have the highest median area across the entire dataset
                    # first get the median area ratio by ion
                    median_area_per_ion <- tapply(fragment_ion_results$calculated_area_ratio,
                                                  INDEX = fragment_ion_results$fragment_ion,
                                                  FUN=median,
                                                  na.rm=TRUE)
                    # remove the "all" entry (since that gets included anyway)
                    median_area_per_ion <- median_area_per_ion[names(median_area_per_ion)!="all"]
                    three_highest_area_ratios <- names(sort(median_area_per_ion, decreasing = TRUE))[1:3]
                    ions_to_plot  <- c(three_highest_area_ratios, 'all')
                }
                
                # If the number of the fragment ions is less than 3, a warning is raised.
                # But this is allowed, and it's commented.

                #if (length(ions_to_plot) <4) {
                #    errorType <- "Warning"
                #    errorSubtype <- "Fragment ion"
                #    errorReason <- paste('The number of fragment ion is less than 3.')
                #    errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge), rep('', colNumber-5)), collapse='\t')
                #    cat(errorInfor)
                #    cat('\n')
                #}
                
                # For each Fragment ion, if its spike level is less than 3, a warning message will be raised.
                # If the replicate number in any concentration is less than 2, a warning message will be raised.
                errorReason_missing_point <- c()
                for (current_plot_ion in ions_to_plot ) {
                    plot_fragment_ion_results <- fragment_ion_results[!is.na(fragment_ion_results$calculated_area_ratio), ]
                    plot_fragment_ion_results <- plot_fragment_ion_results[plot_fragment_ion_results$fragment_ion == current_plot_ion, ]
                    errorReason_tmp <- c()
                    spike_level_tmp <- unique(plot_fragment_ion_results$spike_level)
                    if (length(spike_level_tmp) < 3) {
                        errorReason_tmp1 <- "there are Less than 3 spike levels"
                        errorReason_tmp <- c(errorReason_tmp, errorReason_tmp1)
                    }
                    spike_level_with_warnings <- c()
                    spike_level_id <- which(spike_level_tmp==0)
                    spike_level_tmp <- spike_level_tmp[-spike_level_id]
                    plot_fragment_ion_results_tmp <- plot_fragment_ion_results[plot_fragment_ion_results$spike_level != 0, ]
                    for (spike_level_tmp1 in spike_level_tmp) {
                        plot_fragment_ion_results_tmp2 <- plot_fragment_ion_results_tmp[plot_fragment_ion_results_tmp$spike_level==spike_level_tmp1, ]
                        if (length(unique(plot_fragment_ion_results_tmp2$replicate)) < 2) {
                            spike_level_with_warnings <- c(spike_level_with_warnings, spike_level_tmp1)
                        } 
                    }
                    if (length(spike_level_with_warnings) > 0) {
                        errorReason_tmp2 <- paste("there are less than 2 replicates in spike level(s): ",  paste(spike_level_with_warnings, collapse = ", "), ".", sep="")
                        errorReason_tmp <- c(errorReason_tmp, errorReason_tmp2)
                    }
                    
                    if (length(errorReason_tmp) > 0) {
                        errorReason_missing_point <- c(errorReason_missing_point, paste("For fragment ion ", current_plot_ion, ", ", paste(errorReason_tmp, collapse= ", "), sep=""))
                    }
                }
                if (length(errorReason_missing_point) > 0){
                    errorType <- "Warning"
                    errorSubtype <- "Missing points"
                    errorReason <- paste(errorReason_missing_point, collapse= " ")
                    errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge), rep('', colNumber-5)), collapse='\t')
                    cat(errorInfor)
                    cat('\n')
                }
                
                image_frame_count <- 4
                CairoPNG(paste(plot_output_dir, "\\", input_peptide_sequence, "_", input_precursor_charge, "_", curve_type, ".png", sep=""), width=image_frame_count*480, height=400)
                all_plots <- list()
                for (current_plot_ion in ions_to_plot ){
                    # if there are any NAs, throw a warning
                    #if(is.na(sum(fragment_ion_results$calculated_area_ratio))) {
                    #    warning("Some are ratios could not be calculated. They are excluded from the plot and the tables.")
                    #}
                    
                    plot_fragment_ion_results <- fragment_ion_results[!is.na(fragment_ion_results$calculated_area_ratio), ]
                    plot_fragment_ion_results <- plot_fragment_ion_results[plot_fragment_ion_results$fragment_ion == current_plot_ion, ]
                    # do not make plots for fragment ions with no data
                    if (nrow(plot_fragment_ion_results) != 0) {
                        # make QC plot for current ion
                        current_plot <- plot_QC(plot_fragment_ion_results, input_peptide_sequence, current_plot_ion)
                    }
                    
                    all_plots[[current_plot_ion]] <- current_plot
                }
                multiplot(plotlist=all_plots, cols=image_frame_count)
                dev.off()
                # output to files
                write.table(format(summary_table, digits=4), file=paste(plot_output_dir, "\\", input_peptide_sequence, "_", input_precursor_charge, "_", curve_type ,"_summary_table.tsv", sep=""), sep = "\t", qmethod = "double", col.names=TRUE, row.names=FALSE, quote=FALSE)
                # Only keep the rSqares, pValues and intercepts of the ions_to_plot
                summary_table_intercept_plot <- data.frame(matrix(nrow=0, ncol=length(cell_lines)+1))
                colnames(summary_table_intercept_plot) <- c("fragment_ion", cell_lines)
                for (current_plot_ion in ions_to_plot) {
                  summary_table_intercept_plot <- rbind(summary_table_intercept_plot,
                                                        summary_table_intercept[summary_table_intercept$fragment_ion==current_plot_ion, ])
                }
                rownames(summary_table_intercept_plot) <- as.character(c(1:nrow(summary_table_intercept_plot)))
                
                summary_table_rSquare_plot <- data.frame(matrix(nrow=0, ncol=length(cell_lines)+1))
                colnames(summary_table_rSquare_plot) <- c("fragment_ion", cell_lines)
                for (current_plot_ion in ions_to_plot) {
                    summary_table_rSquare_plot <- rbind(summary_table_rSquare_plot,
                                                       summary_table_rSquare[summary_table_rSquare$fragment_ion==current_plot_ion, ])
                }
                rownames(summary_table_rSquare_plot) <- as.character(c(1:nrow(summary_table_rSquare_plot)))
                
                summary_table_pValue_plot <- data.frame(matrix(nrow=0, ncol=length(cell_lines)+1))
                colnames(summary_table_pValue_plot) <- c("fragment_ion", cell_lines)
                for (current_plot_ion in ions_to_plot) {
                    summary_table_pValue_plot <- rbind(summary_table_pValue_plot,
                                                       summary_table_pValue[summary_table_pValue$fragment_ion==current_plot_ion, ])
                }
                rownames(summary_table_pValue_plot) <- as.character(c(1:nrow(summary_table_pValue_plot)))
                # output to files
                write.table(format(summary_table_intercept_plot, digits=4), file=paste(plot_output_dir, "\\", input_peptide_sequence, "_", input_precursor_charge, "_", curve_type ,"_summary_table_intercept.tsv", sep=""), sep = "\t", qmethod = "double", col.names=TRUE, row.names=FALSE, quote=FALSE)
                write.table(format(summary_table_rSquare_plot, digits=4), file=paste(plot_output_dir, "\\", input_peptide_sequence, "_", input_precursor_charge, "_", curve_type ,"_summary_table_rSquare.tsv", sep=""), sep = "\t", qmethod = "double", col.names=TRUE, row.names=FALSE, quote=FALSE)
                write.table(format(summary_table_pValue_plot, digits=3), file=paste(plot_output_dir, "\\", input_peptide_sequence, "_", input_precursor_charge, "_", curve_type ,"_summary_table_pValue.tsv", sep=""), sep = "\t", qmethod = "double", col.names=TRUE, row.names=FALSE, quote=FALSE)
                # locate abnormal values of r square or p value.
                rSquare_abnormal_list <- locate_abnormal_value(summary_table_rSquare_plot, "rSquare", rSquare_threshold)
                #pValue_abnormal_list <- locate_abnormal_value(summary_table_pValue_plot, "rSquare", rSquare_threshold)
                if (length(rSquare_abnormal_list) > 0) {
                    errorType <- "Warning"
                    errorSubtype <- "Bad linear regression fitting"
                    rSquare_abnormal_df <- data.frame(matrix(unlist(rSquare_abnormal_list), ncol=2, byrow=T),stringsAsFactors=FALSE)                    
                    errorReason <- paste('The quality of fit of linear model is poor due to R2 < ', rSquare_threshold, ' for fragment ion ', paste(paste(rSquare_abnormal_df[,1], ' in sample group of ', rSquare_abnormal_df[,2], sep=''), collapse=", "), '.', sep='')
                    errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge), rep('', colNumber-5)), collapse='\t')
                    cat(errorInfor)
                    cat('\n')
                }
                
                # Abnormal p value need to be confirmed by Simina. Now it's commented. 
                #if (length(pValue_abnormal_list) > 0) {
                #    errorType <- "Warning"
                #    errorSubtype <- "Bad linear regression fitting"
                #    pValue_abnormal_df <- data.frame(matrix(unlist(pValue_abnormal_list), ncol=2, byrow=T),stringsAsFactors=FALSE)                    
                #    errorReason <- paste('The quality of fit of linear model is poor due to p value > ', pValue_threshold, ' for fragment ion ', paste(paste(pValue_abnormal_df[,1], ' in sample group of ', pValue_abnormal_df[,2], sep=''), collapse=", "), '.', sep='')
                #    errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge), rep('', colNumber-5)), collapse='\t')
                #    cat(errorInfor)
                #    cat('\n')
                #}

                # Calculate the CVs for the slops of fragment ions across different cell lines.
                # Calculate the CVs for the slops of cell lines across different fragment ions.
                summary_table_tmp1 <- data.frame(matrix(nrow=0, ncol=length(cell_lines)+1))
                colnames(summary_table_tmp1) <- c("fragment_ion", cell_lines)
                for (current_plot_ion in ions_to_plot) {
                    summary_table_tmp1 <- rbind(summary_table_tmp1,
                                                       summary_table[summary_table$fragment_ion==current_plot_ion, ])
                }
                rownames(summary_table_tmp1) <- as.character(c(1:nrow(summary_table_tmp1)))
                summary_table_tmp <- summary_table_tmp1[2:ncol(summary_table_tmp1)]
                cv_table <- data.frame(component=as.character(), type=as.character(), cv=as.numeric())
                row_cv <- apply(summary_table_tmp, 1, cv_cal)
                col_cv <- apply(summary_table_tmp, 2, cv_cal)
                for (indexComponent in c(1:length(ions_to_plot))) {
                    cv_table_tmp <- data.frame(component=ions_to_plot[indexComponent], type="fragment ion", cv=row_cv[indexComponent])
                    cv_table <- rbind(cv_table, cv_table_tmp)
                }
                for (indexComponent in c(1:ncol(summary_table_tmp))) {
                    cv_table_tmp <- data.frame(component=colnames(summary_table_tmp)[indexComponent], type="cell line", cv=col_cv[indexComponent])
                    cv_table <- rbind(cv_table, cv_table_tmp)
                }
                rownames(cv_table) <- c(1:nrow(cv_table))
                nan_list <- rownames(cv_table)[is.nan(cv_table$cv)]
                if (length(nan_list) > 0) {
                    errorType <- "Warning"
                    errorSubtype <- "Bad linear regression fitting"
                    errorReason <- paste(paste('The slopes of ', cv_table[nan_list, ]$type, ' of ', cv_table[nan_list, ]$component, ' is(are) 0', sep=''), collapse= '. ')
                    errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge), rep('', colNumber-5)), collapse='\t')
                    cat(errorInfor)
                    cat('\n')
                }
                large_cv_list <- rownames(cv_table)[cv_table$cv > cv_threshold]
                if (length(large_cv_list) > 0) {
                    errorType <- "Warning"
                    errorSubtype <- "Bad linear regression fitting"
                    errorReason <- paste(paste('The coefficient of variance of the slopes of ', cv_table[large_cv_list, ]$type, ' of ', cv_table[large_cv_list, ]$component, ' is ', format(cv_table[large_cv_list, ]$cv, digits=4), ' > ', cv_threshold, sep=''), collapse= '. ')
                    errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge), rep('', colNumber-5)), collapse='\t')
                    cat(errorInfor)
                    cat('\n')
                }
                write.table(format(cv_table, digits=3), file=paste(plot_output_dir, "\\", input_peptide_sequence, "_", input_precursor_charge, "_", curve_type ,"_cv_table.tsv", sep=""), sep = "\t", qmethod = "double", col.names=TRUE, row.names=FALSE, quote=FALSE)
                
                # Calculate the rmse for the intercepts of fragment ions across different cell lines.
                # Calculate the rmse for the slops of cell lines across different fragment ions.
                # If Exp 3 takes light as internal standard isotopes and spikes heavy isotopes, the expected intercept of y value should be zero.
                # rmse for the intercepts of fragment ions across different cell lines is considered.
                # If Exp 3 takes heavy as internal standard isotopes and spikes light isotopes, we can't garantee each matrix doesn't contain endogenouse proteins. The expected intercept of y values for each matrix could be different.
                # MSInspector can't report warnings about the intercepts.
                if (interceptSwtich) {
                    if (internal_standard == 'light') {
                        summary_table_tmp2 <- data.frame(matrix(nrow=0, ncol=length(cell_lines)+1))
                        colnames(summary_table_tmp2) <- c("fragment_ion", cell_lines)
                        for (current_plot_ion in ions_to_plot) {
                          summary_table_tmp2 <- rbind(summary_table_tmp2,
                                                      summary_table_intercept[summary_table_intercept$fragment_ion==current_plot_ion, ])
                        }
                        rownames(summary_table_tmp2) <- as.character(c(1:nrow(summary_table_tmp2)))
                        summary_table_tmp3 <- summary_table_tmp2[2:ncol(summary_table_tmp2)]
                        cv_table2 <- data.frame(component=as.character(), type=as.character(), rmse=as.numeric())
                        row_cv2 <- apply(summary_table_tmp3, 1, rmse)
                        col_cv2 <- apply(summary_table_tmp3, 2, rmse)
                        for (indexComponent in c(1:length(ions_to_plot))) {
                          cv_table_tmp <- data.frame(component=ions_to_plot[indexComponent], type="fragment ion", rmse=row_cv2[indexComponent])
                          cv_table2 <- rbind(cv_table2, cv_table_tmp)
                        }
                        for (indexComponent in c(1:ncol(summary_table_tmp3))) {
                          cv_table_tmp <- data.frame(component=colnames(summary_table_tmp3)[indexComponent], type="cell line", rmse=col_cv2[indexComponent])
                          cv_table2 <- rbind(cv_table2, cv_table_tmp)
                        }
                        
                        rownames(cv_table2) <- c(1:nrow(cv_table2))
                        large_cv_list2 <- rownames(cv_table2)[cv_table2$rmse > rmse_threshold]
                        #rmse_list <- c(rmse_list, cv_table2$rmse)
                        if (length(large_cv_list2) > 0) {
                          errorType <- "Warning"
                          errorSubtype <- "Bad linear regression fitting"
                          errorReason <- paste(paste('The RMSE of intercepts at y axis of ', cv_table2[large_cv_list2, ]$type, ' of ', cv_table2[large_cv_list2, ]$component, ' is ', format(cv_table2[large_cv_list2, ]$rmse, digits=4), ' > ', rmse_threshold, sep=''), collapse= '. ')
                          errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge), rep('', colNumber-5)), collapse='\t')
                          cat(errorInfor)
                          cat('\n')
                        }
                        write.table(format(cv_table2, digits=3), file=paste(plot_output_dir, "\\", input_peptide_sequence, "_", input_precursor_charge, "_", curve_type ,"_rmse_table_intercept.tsv", sep=""), sep = "\t", qmethod = "double", col.names=TRUE, row.names=FALSE, quote=FALSE)
                    }
                }
                if (FALSE) {
                    fragment_ion_results_plot <- fragment_ion_results[fragment_ion_results$fragment_ion %in% ions_to_plot, ]
                    fragment_ion_results_plot_2 <- fragment_ion_results_plot[fragment_ion_results_plot$spike_level > 0, ]
                    fragment_ion_results_plot_2 <- as.data.frame(fragment_ion_results_plot_2 %>%
                                                     group_by(Peptide, Protein_Name, Precursor_Charge, fragment_ion, spike_level) %>%
                                                     summarise(Ratio_mean=mean(calculated_area_ratio), Ratio_max=max(calculated_area_ratio), Ratio_min=min(calculated_area_ratio)))
                    fragment_ion_results_plot_3 <- fragment_ion_results_plot_2[!((fragment_ion_results_plot_2$Ratio_max <= 1.3*fragment_ion_results_plot_2$Ratio_mean & fragment_ion_results_plot_2$Ratio_max >= 0.7*fragment_ion_results_plot_2$Ratio_mean) & (fragment_ion_results_plot_2$Ratio_min <= 1.3*fragment_ion_results_plot_2$Ratio_mean & fragment_ion_results_plot_2$Ratio_min >= 0.7*fragment_ion_results_plot_2$Ratio_mean)),]
                    if (nrow(fragment_ion_results_plot_3) > 0) {
                        unique_fragment_ion <- unique(fragment_ion_results_plot_3$fragment_ion)
                        errorType <- "Warning"
                        errorSubtype <- "High variance"
                        errorReasonTmp <- c()
                        for (fragment_ion_tmp in unique_fragment_ion) {
                          errorReasonTmp <- c(errorReasonTmp, paste("for fragment ion ", fragment_ion_tmp, " at the spike level(s) of ", paste(fragment_ion_results_plot_3[fragment_ion_results_plot_3$fragment_ion == fragment_ion_tmp, ]$spike_level, collapse=', '), sep = ''))
                        }
                        errorReason <- paste('When checking the samples with the non blank spike levels, ', paste(errorReasonTmp, collapse = ', '), ', not all of the area ratios are within 30% of the mean.')
                        errorInfor <- paste(c(c(SkyDocumentName, errorType, errorSubtype, errorReason, input_protein_name, input_peptide_sequence, '', input_precursor_charge), rep('', colNumber-5)), collapse='\t')
                        cat(errorInfor)
                        cat('\n')
                    }
                }
                
                if (nonBlankConcentrationSwitch) {
                    fragment_ion_results_plot <- fragment_ion_results[fragment_ion_results$fragment_ion %in% ions_to_plot, ]
                    fragment_ion_results_plot_2 <- fragment_ion_results_plot[fragment_ion_results_plot$spike_level > 0, ]
                    
                    fragment_ion_results_tmp_sum <- fragment_ion_results_plot_2[fragment_ion_results_plot_2$fragment_ion=='all', ]
                    fragment_ion_results_tmp_idividual <- fragment_ion_results_plot_2[fragment_ion_results_plot_2$fragment_ion!='all', ]
                    
                    heavyArea_transition_ratio <- c()
                    lightArea_transition_ratio <- c()
                    for (i in 1:nrow(fragment_ion_results_tmp_idividual)) {
                      row <- fragment_ion_results_tmp_idividual[i,]
                      row_sum_tmp <- fragment_ion_results_tmp_sum[fragment_ion_results_tmp_sum$cell_line==row$cell_line & fragment_ion_results_tmp_sum$spike_level==row$spike_level & fragment_ion_results_tmp_sum$replicate==row$replicate, ]
                      heavyArea_transition_ratio <- c(heavyArea_transition_ratio, row$heavy_area/row_sum_tmp$heavy_area)
                      lightArea_transition_ratio <- c(lightArea_transition_ratio, row$light_area/row_sum_tmp$light_area)
                    }
                    fragment_ion_results_tmp_idividual$heavyArea_transition_ratio <- heavyArea_transition_ratio
                    fragment_ion_results_tmp_idividual$lightArea_transition_ratio <- lightArea_transition_ratio
                    
                    # Step1: Check whether there are all samples in non-blank levels, no all of transition ratio are within 30% of the mean.
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
                            errorReason <- paste('When checking samples with the non blank spike levels, ', paste(errorReasonTmp, collapse = ', '), '.', sep='')
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
                                  errorReason <- paste('When checking samples with the non blank spike levels, ', paste(errorReasonTmp, collapse = ', '), '.', sep='')
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