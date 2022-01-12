#!/usr/bin/env Rscript

#' @title ##methyl_master_formatting_custom
#' @description Formatting the custom routine data
#' This function is for specifying and running a custom data set
#' CNV methylation analysis
#' @param custom.idat.files.dir The input directory for the sesame routine
#' @param custom.output.dir The ouptut directory for the sesame routine
#' @param custom.sample.sheet.path The path to the MethylMaster sample sheet
#' @param custom.comparison the 2-element vector of Sample_Group levels to be
#' compared. First element is taken as the treatment and second as the control
#' if "reference" is set to "internal" the second element is ignored
#' @param custom.file.sep the file separator to use
#' @param custom.data.cache the sesame data cache to use
#' @param custom.data.normal the sesame normal data set to use,
#' e.g. Epic.5.Normal
#' @param custom.ref.version the sesame reference version (default is hg38)
#' @param custom.reference the sesame reference to use
#' @param custom.split.by which column, if any, to split the analyses by
#' @param custom.save.seg save the segmentation results as .RData object
#' @param custom.hm450.mean.correct apply the hm450 AutoCorrectPeak algorithm
#' @param custom.form.add.meta additional metadata column to be passed along
#' @param ... additional parameters to passs to methyl_master_sesame
#' @import data.table
#' @import dplyr
#' @import sesame
#' @import sesameData
#' @import ExperimentHub
#' @import future
#' @import profvis
#' @return A seg object for downstream analysis
#' @export
