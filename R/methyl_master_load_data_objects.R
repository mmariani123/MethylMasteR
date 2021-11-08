#!/usr/bin/env Rscript

##Michael Mariani Dartmouth College 2021

##Required data objects for MethylMaster

##Load required data objects:

##load(file=paste0(files.dir,"/","hm450.manifest.hg38.rda"))
##annotation_df <- as.data.frame(
##  data.table::fread(
##    paste0(work.dir,
##           file.sep,
##           "HM450.hg38.manifest.tsv")
##  )
##)

methyl_master_load_data_objects <- function(files.dir, file.sep){

##Need the probe450kfemanno obj for champ
load(paste0(files.dir,
            file.sep,
            "probe450kfemanno.rda"))
annotation_df <- probe450kfemanno
##rm(probe450kfemanno)
return(annotation_df)

}
