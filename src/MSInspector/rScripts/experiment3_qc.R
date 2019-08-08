# This sample code is used to qc skyline document data from experiment 3 of Assay Portal.
# It's modified based on:
#   esac-panorama-master\experiment-3\code\Experiment_3_with_labkey_connection_updated_only_plot_4_ions.R

suppressWarnings(suppressMessages(library(Cairo)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(stringr)))

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

args <- commandArgs(trailingOnly = TRUE)
dataset_path <- args[1]
fileList_path <- args[2]
plot_output <- args[3]
plot_output_dir <- args[4]

#dataset_path <- "normal_data.tsv"
#fileList_path <- "file_namelist_IS.tsv"
#plot_output <- "True"
#plot_output_dir <- "D:\\Skyline_analysis\\qcAssayPortal\\qcAssayPortal\\src\\qcAssayPortal\\rScripts\\test\\debug_exp3\\tmp"

if (plot_output == 'True') {
    plot_output <- TRUE
} else {
    plot_output <- FALSE
}

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
colnames(QC_set_total)[colnames(QC_set_total)=="exp3samplegroup"] <- "sample_group"     # individual samples of the matrix of interest
colnames(QC_set_total)[colnames(QC_set_total)=="isotopelabeltype"] <- "isotope_label_type"     # light, heavy
colnames(QC_set_total)[colnames(QC_set_total)=="area"] <- "area"

QC_set_total$fragment_ion <- paste(QC_set_total[ ,'fragment_ion_only'], " (", QC_set_total[ ,'product_charge'], "+)", sep='' )

ion_category <- 'error'

# convert some classes for some columns
QC_set_total[,'sample_group'] <- as.character(QC_set_total[,'sample_group'])
# remove factor version
QC_set_total[,'fragment_ion'] <- as.character(QC_set_total[,'fragment_ion'])
QC_set_total[,'isotope_label_type'] <- as.character(QC_set_total[,'isotope_label_type'])
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
             }
        }
    }
    df_skydoc_error_peptide[[SkyDocumentName]] <- peptide_list_with_error
}

if (plot_output) {
    write.table(logdf, file=log_filename, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
}

#########################################
# Infer the internal standard type for each SkyDocumentName by randomly sampled 5 peptides.
# Since the internal standard type has to be predefined to be heavy in experiment, there is no need to infer the internal standard type.
#########################################
df_internal_standard_inferred <- data.frame(SkyDocumentName=as.character(), internal_standard=as.character())
for (SkyDocumentName in as.character(fileDf[, "SkyDocumentName"])) {
    df_internal_standard_inferred_tmp <- data.frame(SkyDocumentName=SkyDocumentName, internal_standard="heavy")
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
                # *** This part will be moved into error detection function chunk. ***
                if(('light' %in% isotope_label_types) & ('medium' %in% isotope_label_types)) {
                    stop("both light and medium isotope label found in dataset")
                } else {
                    QC_set$isotope_label_type[QC_set$isotope_label_type == "medium"] <- "light"
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
                                # *** This part will be moved into error detection function chunk. ***
                                if(nrow(current_set[current_set$isotope_label_type=='light', ]) > 1){
                                    stop("more than one light isotope")
                                }
                                if(nrow(current_set[current_set$isotope_label_type=='heavy', ]) > 1){
                                    stop("more than one heavy isotope")
                                }
                                
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
                                # *** This part will be moved into error detection function chunk. ***
                                if(nrow(current_set[current_set$isotope_label_type=='light', ]) > 1){
                                    stop("more than one light isotope")
                                }
                                if(nrow(current_set[current_set$isotope_label_type=='heavy', ]) > 1){
                                    stop("more than one heavy isotope")
                                }
                                
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
                                                                                                fragment_ion = 'all', cell_line = current_cell_line, spike_level = current_spike_level, replicate = current_rep, light_area=light_area, heavy_area=heavy_area,
                                                                                                theoretical_area=theoretical_area, measured_area=measured_area,
                                                                                                calculated_area_ratio=calculated_area_ratio, ion_category='all'))
                            }
                        } # end current_cell_line
                    } # end current_spike_level
                } # end current_rep

                ions <- c(fragment_ion_list, 'all')
                # make summary data tables
                summary_table <- data.frame(matrix(NA, ncol=length(cell_lines)+1, nrow=length(ions)))
                colnames(summary_table) <- c("fragment_ion", cell_lines)
                summary_table$fragment_ion <- ions
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
                        summary_table[summary_table$fragment_ion == current_ion,current_cell_line] <-
                            coef(lm(calculated_area_ratio_ave_reps ~ spike_level, data=cell_line_subset))["spike_level"]
                    }
                    v_current_ion <- cbind(fragment_ion = rep(current_ion, nrow(plot_ave_reps)),
                                           plot_ave_reps)
                    values_for_spike_levels <- rbind(values_for_spike_levels,
                                                       v_current_ion)
                }
                
                # save values for spike levels
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
                
                image_frame_count <- 4
                CairoPNG(paste(plot_output_dir, "\\", input_peptide_sequence, "_", input_precursor_charge, "_", curve_type, ".png", sep=""), width=image_frame_count*480, height=400)
                all_plots <- list()
                for (current_plot_ion in ions_to_plot ){
                    # if there are any NAs, throw a warning
                    if(is.na(sum(fragment_ion_results$calculated_area_ratio))) {
                        warning("Some are ratios could not be calculated. They are excluded from the plot and the tables.")
                    }
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
            }
        }
    }
}