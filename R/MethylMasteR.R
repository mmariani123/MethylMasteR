#!/usr/bin/env Rscript

#' Main functions:
#'
#' MultiMethylv1.0 : CNV calling from methylation data
#' Michael Mariani PhD Dartmouth College 2021
#' Possible routines:
#' "test"            ##Run a quick test
#' "process_sesame", ##Preprocess the TCGA and cord data in sesame format
#' "sesame",         ##Run Sesame  CNV calling (get segments)
#' "k450" ,          ##Run 450K    CNV calling (get segments)
#' "champ" ,         ##Run ChAMP   CNV calling (get segments)
#' "epicopy" ,       ##Run EpiCopy CNV calling (get segments)
#'
#' @param infile Path to the input file
#' @param os.type
#' @param sample.sheet.path
#' @param idat.folder.name
#' @param r.lib.path
#' @param proj
#' @param visualize
#' @param visualize.individual
#' @param weighted.mean
#' @param routine
#' @param sex
#' @param main.dir
#' @importFrom magrittr %>%
#' @import profmem
#' @import profvis
#' @return #NULL
#' @export
methyl_master <- function(os.type              = "linux",
                          sample.sheet.path    = "Sample_Sheet.csv",
                          idat.folder.name     = "blca_idat_download",
                          r.lib.path           = .libPaths()[1],
                          proj                 = "TCGA-BLCA",
                          visualize            = FALSE,
                          visualize.individual = FALSE,
                          weighted.mean        = "normal",
                          routine              = "sesame",
                          sex                  = "gender_reported",
                          main.dir             = getwd()){

#################### Load global variables ####################################

.libPaths(r.lib.path)

if(os.type=="windows"){
  file.sep="\\\\"
}else{
  file.sep="/"
}

file.spacer="_"

files.dir <- paste0(main.dir,
                   file.sep,
                   "files")

data.dir <- paste0(main.dir,
                file.sep,
                "data")

if(routine=="test"){

#################### LOAD PACKAGES ############################################

source(paste0(files.dir,
              file.sep,
              "methyl_master_load_packages.R"))

  methyl_master_load_packages(routine)

}else{

source(paste0(files.dir,
                file.sep,
                "methyl_master_load_packages.R"))

  methyl_master_load_packages(routine)

######################## LOAD FUNCTIONS #######################################

source(paste0(files.dir,
              file.sep,
              "methyl_master_helper_functions.R"))

source(paste0(files.dir,
			  file.sep,
              "methyl_master_champ_extra_functions.R"))

source(paste0(files.dir,
              file.sep,
              "methyl_master_findSegments2.R"))

source(paste0(files.dir,
              file.sep,
              "methyl_master_population_ranges.R"))

source(paste0(files.dir,
              file.sep,
              "methyl_master_CopyNumber450kCancer.R"))

######################## LOAD OBJECTS ########################################

#source(paste0(files.dir,
#              file.sep,
#              "methyl_master_load_data_objects.R"))
#
#annotation_df <- methyl_master_load_data_objects(files.dir=files.dir,
#                                                 file.sep=file.sep)

######################## SAMPLE SHEET ########################################

##Sample Sheet
sample.sheet.csv <- paste0(sample.sheet.path) %>%
                      read.csv(header = TRUE,
                               stringsAsFactors = FALSE)

}

########################### SETUP ############################################

if(routine=="test"){
  ref <- "testing"
  sub.workflow <- "tested"
}else if(routine=="download"){
  ref <- "downloading"
  sub.workflow <- "downloaded"
}else if(routine=="process_sesame"){
  ref <- "processing_sesame"
  sub.workflow <- "processed_sesame"
  cord.files.path <- paste0(data.dir,
                            file.sep,
                            sesame.cord.dir)
}else if(routine=="sesame"){
  ref <- sesame.ref
  sub.workflow <- "complete"
  if(sesame.ref=="norm"){
    ##Get sesame normal samples:
    sesameData::sesameDataCache(sesame.data.cache)
    sesame_ssets_normal <- sesameData::sesameDataGet(sesame.data.normal)
  }else if(sesame.ref=="cord"){
    cord.files.path <- paste0(data.dir,
                              file.sep,
                              sesame.cord.dir)
  }
  idat.files.dir <- paste0(data.dir,
                           file.sep,
                           idat.folder.name)
}else if(routine=="k450"){
  ref <- k450.ref
  sub.workflow <- k450.workflow
  idat.files.dir <- paste0(data.dir,
                           file.sep,
                           idat.folder.name)
}else if(routine=="champ"){
  ref <- champ.control.group
  sub.workflow <- "complete"
  idat.files.dir <- paste0(data.dir,
                           file.sep,
                           idat.folder.name)
}else if(routine=="epicopy"){
  ref <- epi.ref
  sub.workflow <- "complete"
  idat.files.dir <- paste0(data.dir,
                           file.sep,
                           idat.folder.name)
  epi.target.dir <- idat.pooled.files.dir
}else{
  stop("Error: please select a valid <routine>")
}

work.dir <- paste0(main.dir,
                   file.sep,
                   routine,
                   file.spacer,
                   ref,
                   file.spacer,
                   sub.workflow)

if(!dir.exists(work.dir)){
  dir.create(work.dir)
}else{
  print("<work.dir> exists, overwriting")
  message("<work.dir> exists, overwriting")
  unlink(work.dir,
         recursive=TRUE,
         force=TRUE)
  dir.create(work.dir)
}

if(visualize.individual==TRUE){
  individual.plots.dir=paste0(work.dir,
                              file.sep,
                              "individual_plots")
  if(!dir.exists(individual.plots.dir)){
    dir.create(individual.plots.dir)
  }else{
    print("<individual.plots.dir> exists, overwriting")
    message("<individual.plots.dir> exists, overwriting")
    unlink(individual.plots.dir,
           recursive=TRUE,
           force=TRUE)
    dir.create(individual.plots.dir)
  }
}

setwd(work.dir)

######################## Run Pipeline ##############################
####################################################################
####################################################################
####################################################################
####################################################################

##profvis.out <- profvis({

profmem.out <- profmem::profmem({

start.time <- Sys.time()

switch(routine,

######################## Run Test ##################################
####################################################################
####################################################################
####################################################################
####################################################################

test = {

  print("Hello World")
  profvis::pause(10)

},

####################### DOWNLOAD ############################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################

download = {

  source(Paste0(scripts.dir,
                file.sep,
                "methyl_master_download.R")
  )

},

########################## Proces SeSAMe #############################
######################################################################
######################################################################
######################################################################
######################################################################

process_sesame = {

  source(Paste0(scripts.dir,
                file.sep,
                "methyl_master_process_sesame.R")
         )

},

################## SeSAMe CNSegmentation #############################
######################################################################
######################################################################
######################################################################
######################################################################

sesame = {

  source(Paste0(scripts.dir,
                file.sep,
                "methyl_master_sesame.R")
  )

},

################# segmentation: cnv450k analysis ############################
#############################################################################
#############################################################################
#############################################################################
#############################################################################

k450 = {

  source(Paste0(scripts.dir,
                file.sep,
                "methyl_master_k450.R")
  )

}, ##End K450

############################# ChAMP ##########################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################

champ = {

  source(Paste0(scripts.dir,
                file.sep,
                "methyl_master_champ.R")
  )

}, ##End ChAMP

############################# Epicopy ########################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################

epicopy = {

  source(Paste0(scripts.dir,
                file.sep,
                "methyl_master_epicopy.R")
  )

} ##End Epicopy

) ##End switch routine

################### TIME/MEM PROFILE ####################
#########################################################
#########################################################
#########################################################
#########################################################

##start.time <- Sys.time()
end.time <- Sys.time()

total.time <<- end.time - start.time

}) ##End profmem

##}) ##End profvis

print(total.time)

writeLines(c("\nTotal time:\n",
             total.time,
             "\nGC output: \n",
             capture.output(gc()),
             "\nProfmem max output (bytes): \n",
             {profmem.out$bytes[!is.na(profmem.out$bytes)] %>%
                 max()}),
           paste0(work.dir,
                  file.sep,
                  "time_mem.",
                  routine,
                  ".",
                  gsub("-|\\s|:",".",Sys.time(), perl = TRUE),
                  ".txt")
)

##For profvis:
##htmlwidgets::saveWidget(profvis.out,
##                        paste0(work.dir,
##                               file.sep,
##                               routine,
##                               ".",
##                               gsub("-|\\s|:",".",
##                                    Sys.time(),
##                                    perl = TRUE),
##                               ".profvis.out.html"))

## Can open in browser from R
## browseURL(paste0(work.dir,file.sep,"profile.html"))

##gc.output <- capture.output(gc())
##print(profvis.out)

####################### Overlaps and visualize ###############################
##############################################################################
##############################################################################
##############################################################################
##############################################################################

if(visualize==TRUE){

  source(Paste0(scripts.dir,
                file.sep,
                "methyl_master_visualize.R")
  )

}

##################### Write Session Info #####################################

writeLines(capture.output(sessionInfo()),
           paste0(work.dir,
                  file.sep,
                  "sessionInfo.",
                  routine,
                  ".",
                  gsub("-|\\s|:",".",Sys.time(), perl = TRUE),
                  ".txt")
)

##################### End pipeline ###########################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################

}
