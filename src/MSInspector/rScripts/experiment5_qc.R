# This sample code is used to qc skyline document data from experiment 5 of Assay Portal.
# It's modified based on:
#   esac-panorama-master\experiment-5\code\Experiment_5_for_integration_Panorama.R

# UPDATE: For experiment 5, removed all sample group info, as there should be a single sample
# In the updated skyline template of experiment 5, it should have heavy as the internal standard type. (The default internal standard type is heavy, no matter whether it's set or not)

suppressWarnings(suppressMessages(library(Cairo)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(stringr)))

plot_QC <- function(plot_fragment_ion_results, input_peptide_sequence, current_ion, days) {


  if (current_ion == 'all'){

      current_ion <- 'sum of ions'
  }


  plot_title <- paste(input_peptide_sequence, current_ion, sep='\n')

  min_calculated_area_ratio <- min(plot_fragment_ion_results$calculated_area_ratio)
  max_calculated_area_ratio <- max(plot_fragment_ion_results$calculated_area_ratio)



  # Expand right side of clipping rect to make room for the legend
  par(xpd=TRUE, mar=par()$mar+c(0,0,0,4))

  # bty="L",

  ##get minimum and maximum calculated area ratios
  min_area_ratio <- min(plot_fragment_ion_results$calculated_area_ratio)
  max_area_ratio <- max(plot_fragment_ion_results$calculated_area_ratio)
  
  suppressWarnings(plot(plot_fragment_ion_results$day, plot_fragment_ion_results$calculated_area_ratio, log="y", yaxt="n", 
       pch=ifelse(plot_fragment_ion_results$replicate==1 , 1,
                  ifelse(plot_fragment_ion_results$replicate==2 , 0,
                         ifelse(plot_fragment_ion_results$replicate==3 , 2,
                                ifelse(plot_fragment_ion_results$replicate==4 , 5,
                                       ifelse(plot_fragment_ion_results$replicate==5 , 6,
                                4))))),

       cex=2, lwd=1, main=plot_title, xlab="Time (day)", ylab="Measured (area ratio) [log-scale]", 
       cex.lab=1.75, cex.axis=1.75,
       ylim = c(min_area_ratio*0.8,max_area_ratio*1.2)))


  #x_axis_values <- c(1,2,3,4,5)
  x_axis_values <- days

  #y_axis_values <- c(0.01,0.1,1,10,100,format(max_calculated_measured_concentration,digits=3))
  y_axis_values <- c(format(min_calculated_area_ratio,digits=3),format(median(plot_fragment_ion_results$calculated_area_ratio),digits=3),format(max_calculated_area_ratio,digits=3))




  #format(y_axis_values,scientific=FALSE,digits=4)

  axis(2, y_axis_values, labels=format(y_axis_values,scientific=FALSE)) # draw y axis with required labels

  #par(xpd=TRUE)
  #legend(x=4,y=max_calculated_area_ratio,legend=c("rep1","rep2","rep3","repX","Hi","Med","Lo"),pch=c(1,0,2,4,18,18,18), col=c("black","black","black","black","red","blue","green"), bty="n")

  legend(x=max(days)+0.2,y=max_calculated_area_ratio,legend=c("rep1","rep2","rep3","rep4","rep5"),pch=c(1,0,2,5,6), col=c("black","black","black","black","black"), cex=1, bty="n")

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
#plot_output_dir <- "D:\\Skyline_analysis\\qcAssayPortal\\qcAssayPortal\\src\\qcAssayPortal\\rScripts\\test\\8\\tmp"

if (plot_output == 'True') {
    plot_output <- TRUE
} else {
    plot_output <- FALSE
}

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
                        sample_group= samplegroup,
                        isotope_label_type = isotopelabeltype,
                        area = area)

QC_set_total$fragment_ion <- paste(QC_set_total[ ,'fragment_ion_only'], " (", QC_set_total[ ,'product_charge'], "+)", sep='' )

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
                days <- sort(unique(QC_set[ , 'day']))
                replicates <- sort(unique(QC_set[ , 'replicate_number']))
                isotope_label_types <- unique(QC_set[ , 'isotope_label_type'])
                
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
                    ions_to_plot <- ions
                } else {
                    results_to_plot <- CV_results[CV_results$fragment_ion!='all' & !is.na(CV_results$total_CV), ]
                    # new sort to get Top 3 plots
                    results_to_plot <- results_to_plot[order(results_to_plot$total_CV), ]
                    three_lowest_total_CV <- head(results_to_plot, 3)
                    three_lowest_ions <- as.character(three_lowest_total_CV[ , 'fragment_ion'])
                    ions_to_plot  <- c(three_lowest_ions, 'all')
                }
                
                par(mfrow = c(1,length(ions_to_plot)))
                
                for (current_plot_ion in ions_to_plot ) {
                    plot_fragment_ion_results <- fragment_ion_results[!is.na(fragment_ion_results$calculated_area_ratio), ]
                    plot_fragment_ion_results <- plot_fragment_ion_results[plot_fragment_ion_results$fragment_ion == current_plot_ion, ]
                    plot_days <- sort(unique(plot_fragment_ion_results[ , 'day']))
                    # do not make plots for fragment ions with no data
                    if (nrow(plot_fragment_ion_results) != 0) {
                        # make QC plot for current ion
                        plot_QC(plot_fragment_ion_results, input_peptide_sequence, current_plot_ion, plot_days)
                    }
                }
                dev.off()
                write.table(CV_results, file=paste(plot_output_dir, "\\", input_peptide_sequence, "_", input_precursor_charge, "_", curve_type ,"_CV_results", ".tsv", sep=""), sep = "\t", qmethod = "double", col.names=TRUE, row.names=FALSE, quote=FALSE)
            }
        }
    }
}