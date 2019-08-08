# This sample code is used to qc skyline document data from experiment 4 of Assay Portal.
# It's modified based on:
# esac-panorama-master\experiment-4\code\Experiment_4_with_labkey_connection_updated_v3_new_inter_CV_calc.R

suppressWarnings(suppressMessages(library(Cairo))) # need for producing PNG image
suppressWarnings(suppressMessages(library(Rlabkey)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(stringr)))

##variable that specifies whether the script takes a standalone csv file vs going through Panorama
options(warn=-1)

# Multiple plot function - from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  #require(grid)
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

# Sort the Time levels in the order of Control, Autosampler, Frozen, FTx1 and FTx2
sortLevel <- function(timeFactors) {
    time_points <- as.character(unique(timeFactors))
    df_tmp <- data.frame(time_point<-as.character(), sample_group<-as.character(), sample_group_type<-as.numeric(), time_value<-as.numeric())
    for (time_point in time_points) {
      sample_group <- strsplit(time_point, " ")[[1]][1]
      if (sample_group == "Control") {
        sample_group_type <- 1
      } else if (sample_group == "Frozen") {
        sample_group_type <- 3
      } else if (sample_group == "FTx1") {
        sample_group_type <- 4
      } else if (sample_group == "FTx2") {
        sample_group_type <- 5
      }
      else {
        sample_group_type <- 2
      }
      time_value <- strsplit(time_point, " ")[[1]][3]
      if (suppressWarnings(!is.na(as.numeric(time_value)))) {
        time_value <- as.numeric(time_value)
        df_tmp1 <- data.frame(time_point=time_point, sample_group=sample_group, sample_group_type=sample_group_type, time_value=time_value)
      } else {
        df_tmp1 <- data.frame(time_point=time_point, sample_group=sample_group, sample_group_type=sample_group_type, time_value=0)
      }
      df_tmp <- rbind(df_tmp, df_tmp1)
    }
    # Sort by sample_group_type, then sort by time_value
    df_tmp <- df_tmp[with(df_tmp, order(sample_group_type, time_value)), ]
    as.character(df_tmp$time_point)
}


# ***** plot_QC function *****
plot_QC <- function(plot_fragment_ion_results, input_peptide_sequence, current_ion, times) {
    if (current_ion == 'all'){
        current_ion <- 'sum of ions'
    }
    plot_title <- paste(input_peptide_sequence, current_ion, sep='\n')
    temp <- as.character(plot_fragment_ion_results$sample_group)
    plot_fragment_ion_results$Time <- paste(temp, 
                                            " (", 
                                            plot_fragment_ion_results$time, ")",
                                            sep="")
    # Sort the Time based on the order of Control, 4C, Frozen and FT.
    sortedLevels <- sortLevel(factor(plot_fragment_ion_results$Time))
    plot_fragment_ion_results$Time <- factor(plot_fragment_ion_results$Time,
                                            levels = sortedLevels)
    colnames(plot_fragment_ion_results)[colnames(plot_fragment_ion_results) == "sample_group"] <- 
        "Condition"
    ##set values that are zero to the smallest non-zero value divided by 2
    ##flag this as well
    plot_fragment_ion_results$zero_values <- FALSE
    plot_fragment_ion_results$zero_values[plot_fragment_ion_results$calculated_area_ratio == 0] <- TRUE
    plot_fragment_ion_results$calculated_area_ratio[plot_fragment_ion_results$calculated_area_ratio == 0] <-
        min(plot_fragment_ion_results$calculated_area_ratio[plot_fragment_ion_results$calculated_area_ratio > 0])/2
    ##get minimum and maximum calculated area ratios
    min_area_ratio <- min(plot_fragment_ion_results$calculated_area_ratio)
    max_area_ratio <- max(plot_fragment_ion_results$calculated_area_ratio)
    plot_fragment_ion_results$horiz_line <- min_area_ratio
    ##perform ANOVA and get p-value for hypothesis of no difference between groups
    lm_null <- lm(calculated_area_ratio ~ 1, data = plot_fragment_ion_results)
    lm_groups <- lm(calculated_area_ratio ~ Time, data = plot_fragment_ion_results)
    p_val <- signif(anova(lm_null, lm_groups)$"Pr(>F)"[2],2)
    plot_title <- paste(plot_title, "\n", "P-value from ANOVA: ", p_val, sep="")
    g <- ggplot(plot_fragment_ion_results, aes(x=Time, y=calculated_area_ratio, color=Condition,
                                                shape=Replicate)) +
        geom_point() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle(plot_title) +
        coord_trans(y="log10") +
        scale_y_continuous(limits = c(min_area_ratio*0.8,max_area_ratio*1.2)) +
        theme(plot.title = element_text(hjust = 0.5)) +
        xlab("Condition (time)") + ylab("Measured (area ratio) [log-scale]")
    ##if there were any zeros, add in a horizontal line at the smallest non-zero value divided by 2
    if(sum(plot_fragment_ion_results$zero_values) > 0)
    {
        g <- g + geom_hline(aes(yintercept=horiz_line),color="grey",linetype="dashed") 
        warning("Dashed line indicates values of 0 - They are added at the smallest non-zero value divided by 2.")
    }
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
#plot_output_dir <- "D:\\Skyline_analysis\\qcAssayPortal\\qcAssayPortal\\src\\qcAssayPortal\\rScripts\\test\\debug_exp4\\tmp"

if (plot_output == 'True') {
    plot_output <- TRUE
} else {
    plot_output <- FALSE
}

# Load data from local table
QC_set_total <- read.table(file=dataset_path, header=TRUE, sep='\t')
fileDf <- read.table(file=fileList_path, header=TRUE, sep='\t')

#QC_set <- labkey.selectRows(
#        baseUrl="http://cptac-proliant-linux.esacinc.com/labkey",
#        folderPath=paste("/CPTAC Assay Portal/",dataset_path,"/Stability",sep=""),
#        schemaName="targetedms",
#        queryName="QCAnalysisQuery",
#        viewName="",
#        colFilter=makeFilter(c("Protein", "EQUAL", input_protein_name),
#                             c("PeptideModifiedSequence", "EQUAL", input_peptide_sequence),
#                             c("PrecursorCharge","EQUAL",input_precursor_charge)),
#        containerFilter=NULL
#)


# Transform the column names to match those from embedded panorama query
colNumber <- ncol(QC_set_total)
thenames <- tolower(names(QC_set_total))
thenames <- gsub(" ","", thenames) 
names(QC_set_total) <- thenames

# rename columns in QC_set_total dataframe (replace Panorama names with new names used by R script)
colnames(QC_set_total)[colnames(QC_set_total)=="skydocumentname"] <- "SkyDocumentName"
colnames(QC_set_total)[colnames(QC_set_total)=="proteinname"] <- "protein_name"
colnames(QC_set_total)[colnames(QC_set_total)=="peptidemodifiedsequence"] <- "peptide"
colnames(QC_set_total)[colnames(QC_set_total)=="isotopelabeltype"] <- "isotope_label_type"     # light, heavy
colnames(QC_set_total)[colnames(QC_set_total)=="precursorcharge"] <- "precursor_charge"
colnames(QC_set_total)[colnames(QC_set_total)=="productcharge"] <- "product_charge"
colnames(QC_set_total)[colnames(QC_set_total)=="fragmention"] <- "fragment_ion_only"
colnames(QC_set_total)[colnames(QC_set_total)=="area"] <- "area"
colnames(QC_set_total)[colnames(QC_set_total)=="replicatename"] <- "replicate"
colnames(QC_set_total)[colnames(QC_set_total)=="replicatenumber"] <- "Replicate"
colnames(QC_set_total)[colnames(QC_set_total)=="exp4samplegroup"] <- "sample_group"     # temperature
colnames(QC_set_total)[colnames(QC_set_total)=="time"] <- "time_only"
colnames(QC_set_total)[colnames(QC_set_total)=="timeunits"] <- "time_units"
colnames(QC_set_total)[colnames(QC_set_total)=="freezethawcycles"] <- "freeze_thaw_cycles"

if (nrow(QC_set_total) ==0) {
  QC_set_total$fragment_ion <- integer(0)
} else {
  QC_set_total$fragment_ion <- paste(QC_set_total[ ,'fragment_ion_only'], " (", QC_set_total[ ,'product_charge'], "+)", sep='')
}
QC_set_total$time_units <- tolower(QC_set_total$time_units)
time <-  sapply(1:nrow(QC_set_total), function(x) ifelse(QC_set_total[x,'sample_group'] == "FT", paste(QC_set_total[x,'freeze_thaw_cycles'], 'cycles', sep=' '), paste(QC_set_total[x,'time_only'], QC_set_total[x,'time_units'], sep=' ')))
QC_set_total <- cbind(QC_set_total, time)
QC_set_total$time <- gsub("^0 cycles", "0 cycle", QC_set_total$time)
QC_set_total$time <- gsub("^1 cycles", "1 cycle", QC_set_total$time)
QC_set_total$time <- gsub("^0 hours", "0 hour", QC_set_total$time)
QC_set_total$time <- gsub("^1 hours", "1 hour", QC_set_total$time)
QC_set_total$time <- gsub("^1 days", "1 day", QC_set_total$time)
QC_set_total$time <- gsub("^1 weeks", "1 week", QC_set_total$time)

# convert columns from character to numeric
QC_set_total[,'area'] <- as.numeric(as.character(QC_set_total[,'area']))

# remove factor version
QC_set_total[,'replicate'] <- as.character(QC_set_total[,'replicate'])
QC_set_total[,'Replicate'] <- as.character(QC_set_total[,'Replicate'])
QC_set_total[,'time'] <- as.character(QC_set_total[,'time'])
QC_set_total[,'fragment_ion'] <- as.character(QC_set_total[,'fragment_ion'])
QC_set_total[,'sample_group'] <- as.character(QC_set_total[,'sample_group'])
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
                fragment_ion_list <- unique(QC_set[ , 'fragment_ion'])
                times <- sort(unique(QC_set[ , 'time']))
                sample_groups <- sort(unique(QC_set[ , 'sample_group']))
                Replicates <- sort(unique(QC_set[ , 'Replicate']))
                replicates <- sort(unique(QC_set[ , 'replicate']))
                isotope_label_types <- unique(QC_set[ , 'isotope_label_type'])
                
                # This will be moved into error detection function later.
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
                # # ***** prepare file to print PNG images *****
                # only plot top 3 ions plus one more for sum
                image_frame_count <- 4
                #par(mfrow= c(1, image_frame_count))
                #CairoPNG(filename=paste(plot_output_dir, "\\", input_peptide_sequence, "_", input_precursor_charge, "_", curve_type, ".png", sep=""), width=image_frame_count*480, height=400, bg="white", units="px")
                # Step 1:
                # Traverse fragment_ion_list, times, sample_groups and Replicates
                # to evaluate the fragment_ion under the specific combination of time, sample_group and Replicate.
                for (current_ion in fragment_ion_list) {
                    for (current_time in times) {
                        for (current_sample in sample_groups) {
                            for (current_Rep in Replicates) {
                                current_set_count <- 0
                                light_area <- 0
                                heavy_area <- 0
                                theoretical_area <- 0
                                measured_area <- 0
                                calculated_area_ratio <- 0
                                current_set <- QC_set[QC_set$fragment_ion==current_ion & QC_set$time==current_time & QC_set$sample_group==current_sample & QC_set$Replicate==current_Rep, ]
                                
                                # This will be moved into error detection part.
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
                                                                                                    fragment_ion = current_ion, time = current_time, sample_group = current_sample, Replicate = current_Rep, light_area=light_area, heavy_area=heavy_area,
                                                                                                    theoretical_area=theoretical_area, measured_area=measured_area,
                                                                                                    calculated_area_ratio=calculated_area_ratio, ion_category='individual') )
                                }
                                else {
                                }
                            } # end current_rep
                        } # end current_sample
                    } # end current_day
                } # end current_ion

                # Step 2: repeat calculations for sum of ions 
                for (current_time in times) {
                    for (current_sample in sample_groups) {
                        # Use replicate number here.
                        for (current_Rep in Replicates) {
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
                                current_set <- QC_set[QC_set$fragment_ion==current_ion & QC_set$time==current_time & QC_set$sample_group==current_sample & QC_set$Replicate==current_Rep, ]
                                # This will be moved into error detection part.
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
                                                                                                fragment_ion = 'all', time = current_time, sample_group = current_sample, Replicate = current_Rep, light_area=sum_light_area, heavy_area=sum_heavy_area,
                                                                                                theoretical_area=sum_theoretical_area, measured_area=sum_measured_area,
                                                                                                calculated_area_ratio=calculated_area_ratio, ion_category='all') )
                            }
                        
                        } # end current_rep
                        
                    } # end current_sample
                    
                } # end current_day
                
                # Step 3:
                # ***** calculate CV (Coefficient of Variation) *****
                conditions <- c("Control", "Autosampler", "Frozen", "FTx1", "FTx2")
                ##if condition happens to be something else, then replace Autosampler by that
                for(t in 4:10)
                {
                    temp <- paste(t,"C",sep="")
                    if(length(grep(temp, as.character(fragment_ion_results$sample_group))) > 0) {
                        conditions[2] <- temp
                    }
                }
                ions <- c(fragment_ion_list, 'all')
                # make CV summary data frame
                CV_results <- data.frame(fragment_ion = ions,
                         control_intra_CV = NA,
                         actual_temp_intra_CV = NA,
                         frozen_intra_CV = NA,
                         FTx1_intra_CV = NA,
                         FTx2_intra_CV = NA,
                         all_intra_CV = NA,
                         all_inter_CV = NA,
                         control_count = NA,
                         control_count_light_area_0 = NA,
                         control_count_heavy_area_0 = NA,
                         actual_temp_count = NA,
                         actual_temp_count_light_area_0 = NA,
                         actual_temp_count_heavy_area_0 = NA,
                         frozen_count = NA,
                         frozen_count_light_area_0 = NA,
                         frozen_count_heavy_area_0 = NA,
                         FTx1_count = NA,
                         FTx1_count_light_area_0 = NA,
                         FTx1_count_heavy_area_0 = NA,
                         FTx2_count = NA,
                         FTx2_count_light_area_0 = NA,
                         FTx2_count_heavy_area_0 = NA)
                # *** intra-assay CV ***
                for (current_ion in ions)  {
                    individual_intra_assay_CVs_all <- c()
                    for(current_cond in conditions) {
                        avg_intra_assay_CV <-0
                        individual_intra_assay_CVs <- c()
                        for (current_time in times) {
                            current_set <- fragment_ion_results[fragment_ion_results$fragment_ion==current_ion & 
                                                                    fragment_ion_results$time==current_time &
                                                                    fragment_ion_results$sample_group==current_cond, ]
                            # remove rows with a value of NA for calculated_area_ratio
                            current_set <- current_set[complete.cases(current_set[ , 'calculated_area_ratio' ]),]
                            if (nrow(current_set) <= 1 ){
                                percent_CV <- NA
                            }
                            else{
                                percent_CV <- (sd(current_set$calculated_area_ratio))/(mean(current_set$calculated_area_ratio)) * 100
                                individual_intra_assay_CVs <- c(individual_intra_assay_CVs, percent_CV)
                                individual_intra_assay_CVs_all <- c(individual_intra_assay_CVs_all, percent_CV)
                            }
                        } # end current_day
                        if (length(individual_intra_assay_CVs)==0){
                            avg_CV <- NA
                            count <- 0
                        }
                        else{
                            avg_CV <- mean(individual_intra_assay_CVs, na.rm = TRUE)
                        }
                        if (current_cond=="Control") {
                            CV_results[CV_results$fragment_ion==current_ion, "control_intra_CV"] <- avg_CV
                        }
                        else if(current_cond==conditions[2]){
                            CV_results[CV_results$fragment_ion==current_ion, "actual_temp_intra_CV"] <- avg_CV
                        }
                        else if(current_cond=="Frozen"){
                            CV_results[CV_results$fragment_ion==current_ion, "frozen_intra_CV"] <- avg_CV
                        }
                        else if(current_cond=="FTx1"){
                            CV_results[CV_results$fragment_ion==current_ion, "FTx1_intra_CV"] <- avg_CV
                        } else {
                            CV_results[CV_results$fragment_ion==current_ion, "FTx2_intra_CV"] <- avg_CV
                        }
                    } #end current_cond
                    ##now also get this over all conditions
                    CV_results[CV_results$fragment_ion==current_ion, "all_intra_CV"] <- mean(individual_intra_assay_CVs_all, na.rm = TRUE)
                } # current_ion
                # round for display purposes:
                CV_results[,c("control_intra_CV", "actual_temp_intra_CV", "frozen_intra_CV", "FTx1_intra_CV", "FTx2_intra_CV", "all_intra_CV")] <- 
                  sapply(CV_results[,c("control_intra_CV", "actual_temp_intra_CV", "frozen_intra_CV", "FTx1_intra_CV", "FTx2_intra_CV", "all_intra_CV")],
                         round, digits=1)
                # *** END: intra-assay CV ***
                # *** get the number of samples with areas == 0 in each condition ***
                for (current_ion in ions) {
                    for(current_cond in conditions) {
                        current_set <- fragment_ion_results[fragment_ion_results$fragment_ion==current_ion & 
                                                            fragment_ion_results$sample_group==current_cond, ]
                        if (current_cond=="Control") {
                            CV_results[CV_results$fragment_ion==current_ion,  
                                        c("control_count_light_area_0",
                                        "control_count_heavy_area_0")] <- 
                                c(sum(current_set$light_area == 0),
                                sum(current_set$heavy_area == 0))
                        }
                        if(current_cond==conditions[2]){
                            CV_results[CV_results$fragment_ion==current_ion,  
                                        c("actual_temp_count_light_area_0",
                                        "actual_temp_count_heavy_area_0")] <- 
                                c(sum(current_set$light_area == 0),
                                sum(current_set$heavy_area == 0))
                        }
                        else if(current_cond=="Frozen"){
                            CV_results[CV_results$fragment_ion==current_ion,  
                                        c("frozen_count_light_area_0",
                                        "frozen_count_heavy_area_0")] <- 
                                c(sum(current_set$light_area == 0),
                                sum(current_set$heavy_area == 0))
                        }
                        else if(current_cond=="FTx1"){
                            CV_results[CV_results$fragment_ion==current_ion,  
                                        c("FTx1_count_light_area_0",
                                        "FTx1_count_heavy_area_0")] <- 
                                c(sum(current_set$light_area == 0),
                                sum(current_set$heavy_area == 0))
                        } else {
                            CV_results[CV_results$fragment_ion==current_ion,  
                                        c("FTx2_count_light_area_0",
                                        "FTx2_count_heavy_area_0")] <- 
                                c(sum(current_set$light_area == 0),
                                sum(current_set$heavy_area == 0))
                        }
                    } # end current_cond
                } # end current_ion
                # *** END: get the number of samples with areas == 0 in each condition ***
                # *** inter-assay CV across all conditions ***
                for (current_ion in ions) {
                    current_set <- fragment_ion_results[fragment_ion_results$fragment_ion==current_ion, ]
                    ##for each replicate, get the mean and the standard deviation across all conditions
                    mean_reps <- tapply(current_set$calculated_area_ratio, INDEX=current_set$Replicate, FUN = mean)
                    sd_reps <- tapply(current_set$calculated_area_ratio, INDEX=current_set$Replicate, FUN = sd)
                    percent_CV <- sd_reps/mean_reps * 100
                    ##now calculate overall inter-assay CV!
                    ##need to remove NAs in case any of them are left
                    CV_results[CV_results$fragment_ion==current_ion, "all_inter_CV"] <- mean(percent_CV, na.rm = TRUE)
                } # current_ion
                # *** END: inter-assay CV across all conditions ***
                # round for display purposes:
                CV_results[,"all_inter_CV"] <- sapply(CV_results[,"all_inter_CV"], round, digits=1)
                # determine counts
                for (current_ion in ions){
                    CV_results[CV_results$fragment_ion==current_ion, "control_count"] <- 
                        nrow(fragment_ion_results[fragment_ion_results$sample_group=="Control"&fragment_ion_results$fragment_ion==current_ion & !is.na(fragment_ion_results$calculated_area_ratio), ] )
                    CV_results[CV_results$fragment_ion==current_ion, "actual_temp_count"] <- 
                        nrow(fragment_ion_results[fragment_ion_results$sample_group==conditions[2]&fragment_ion_results$fragment_ion==current_ion & !is.na(fragment_ion_results$calculated_area_ratio), ] )
                    CV_results[CV_results$fragment_ion==current_ion, "frozen_count"] <- 
                        nrow(fragment_ion_results[fragment_ion_results$sample_group=="Frozen"&fragment_ion_results$fragment_ion==current_ion & !is.na(fragment_ion_results$calculated_area_ratio), ] )
                    CV_results[CV_results$fragment_ion==current_ion, "FTx1_count"] <- 
                        nrow(fragment_ion_results[fragment_ion_results$sample_group=="FTx1"&fragment_ion_results$fragment_ion==current_ion & !is.na(fragment_ion_results$calculated_area_ratio), ] )
                    CV_results[CV_results$fragment_ion==current_ion, "FTx2_count"] <- 
                        nrow(fragment_ion_results[fragment_ion_results$sample_group=="FTx2"&fragment_ion_results$fragment_ion==current_ion & !is.na(fragment_ion_results$calculated_area_ratio), ] )
                }
                ions_to_plot <- c()
                # determine fragment ions to plot
                if (length(ions) <= 4 ) {
                    ions_to_plot <- ions
                } else {
                    #results_to_plot <- CV_results[CV_results$fragment_ion != 'all' & 
                                                    #!is.na(CV_results$control_total_CV) & 
                    #                                !is.na(CV_results$actual_temp_total_CV),
                                                    #!is.na(CV_results$frozen_total_CV) & 
                                                    #!is.na(CV_results$FT_total_CV) , 
                    #                                ]
                    results_to_plot <- CV_results[CV_results$fragment_ion != 'all',]
                    # new sort to get Top 3 plots
                    results_to_plot <- results_to_plot[order(results_to_plot$all_inter_CV), ]
                    three_lowest_inter_CV <- head(results_to_plot, 3)
                    three_lowest_ions <- as.character(three_lowest_inter_CV[ , 'fragment_ion'])
                    ions_to_plot  <- c(three_lowest_ions, 'all')
                }
                #graphics.off()
                CairoPNG(filename=paste(plot_output_dir, "\\", input_peptide_sequence, "_", input_precursor_charge, "_", curve_type, ".png", sep=""), width=image_frame_count*480, height=400, bg="white", units="px")
                all_plots <- list()
                for (current_plot_ion in ions_to_plot ){
                    ##if there are any NAs, throw a warning
                    # This will be moved into warning detection part later.
                    if(is.na(sum(fragment_ion_results$calculated_area_ratio))) {
                        warning("Some are ratios could not be calculated. They are excluded from the plot and the tables.")
                    }
                    plot_fragment_ion_results <- fragment_ion_results[!is.na(fragment_ion_results$calculated_area_ratio), ]
                    plot_fragment_ion_results <- plot_fragment_ion_results[plot_fragment_ion_results$fragment_ion == current_plot_ion, ]
                    plot_times <- sort(unique(plot_fragment_ion_results[ , 'time']))
                    # do not make plots for fragment ions with no data
                    if (nrow(plot_fragment_ion_results) != 0) {
                        # make QC plot for current ion
                        current_plot <- plot_QC(plot_fragment_ion_results, input_peptide_sequence, current_plot_ion, plot_times)
                    }
                    all_plots[[current_plot_ion]] <- current_plot
                }
                multiplot(plotlist=all_plots, cols=image_frame_count)
                dev.off()
                # output to files
                write.table(CV_results, file=paste(plot_output_dir, "\\", input_peptide_sequence, "_", input_precursor_charge, "_", curve_type ,"_CV_results", ".tsv", sep=""), sep = "\t", qmethod = "double", col.names=TRUE, row.names=FALSE, quote=FALSE)
            }
        }
    }
}