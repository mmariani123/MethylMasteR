#!/usr/bin/env Rscript

#' @title methyl_master_load_data_objects
#' @description Possible function for loading data objects
#' required for MethylMaster
#' Load required data objects:
#' load(file=paste0(files.dir,"/","hm450.manifest.hg38.rda"))
#' annotation_df <- as.data.frame(
#'  data.table::fread(
#'    paste0(work.dir,
#'           file.sep,
#'           "HM450.hg38.manifest.tsv")
#'  )
#' )
#' @param files.dir The parameter containing the R files
#' @param file.sep The the desired file separator
#' @return loaded data object(s)
#' @export
methyl_master_load_data_objects <- function(files.dir,
                                            file.sep
                                            ){

##Need the probe450kfemanno obj for champ
load(paste0(files.dir,
            file.sep,
            "probe450kfemanno.rda"))
annotation_df <- probe450kfemanno
##rm(probe450kfemanno)
return(annotation_df)

}
