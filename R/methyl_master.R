#!/usr/bin/env Rscript

#' @title methyl_master
#' @description Main function:
#' MultiMethylv1.0 : CNV calling from methylation data
#' Michael Mariani PhD Dartmouth College 2021
#' Possible routines:
#' "test"            ##Run a quick test
#' "process_sesame", ##Preprocess the TCGA and cord data in sesame format
#' "sesame",         ##Run Sesame  CNV calling (get segments)
#' "hm450" ,         ##Run 450K    CNV calling (get segments)
#' "champ" ,         ##Run ChAMP   CNV calling (get segments)
#' "epicopy" ,       ##Run EpiCopy CNV calling (get segments)
#' "compare"         ##Run algorithm comparison functionality
#' @param input.dir
#' @param output.dir
#' @param sample.sheet.path
#' @param file.sep
#' @param r.lib.path
#' @param n.cores
#' @param os.type
#' @param proj
#' @param visualize
#' @param visualize.individual
#' @param weighted.mean
#' @param routine
#' @param reference
#' @param split.by
#' @param comparison
#' @param overlap.density
#' @param epi.run.gistic
#' @param compare.name
#' @param ...
#' @importFrom magrittr %>%
#' @import profmem
#' @import profvis
#' @return #NULL
#' @export
methyl_master <- function(input.dir            = NULL,
                          output.dir           = NULL,
                          sample.sheet.path    = NULL,
                          file.sep             = NULL,
                          r.lib.path           = .libPaths()[1],
                          n.cores              = 1,
                          os.type              = "linux",
                          proj                 = "TCGA-BLCA",
                          visualize            = FALSE,
                          visualize.individual = FALSE,
                          weighted.mean        = "normal",
                          routine              = "test",
                          reference            = NULL,
                          split.by             = NULL,
                          comparison           = NULL,
                          overlap.density      = 0.1,
                          epi.run.gistic       = FALSE,
                          compare.name         = NULL
                          ...
                          ){

#################### Load global variables ####################################

.libPaths(r.lib.path)

if(os.type=="windows"){
  file.sep="\\\\"
}else{
  file.sep="/"
}

file.spacer="_"

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

if(!dir.exists(output.dir)){
  dir.create(output.dir)
}else{
  print("<output.dir> exists, overwriting")
  message("<output.dir> exists, overwriting")
  unlink(output.dir,
         recursive=TRUE,
         force=TRUE)
  dir.create(output.dir)
}

if(visualize.individual==TRUE){
  individual.plots.dir=paste0(output.dir,
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

######################## Run Pipeline ##############################
####################################################################
####################################################################
####################################################################
####################################################################

if(routine!="compare"){

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

################## SeSAMe CNSegmentation #############################
######################################################################
######################################################################
######################################################################
######################################################################

sesame = {

  seg <- methyl_master_sesame(sesame.idat.files.dir = input.dir,
                       sesame.output.dir            = output.dir,
                       sesame.sample.sheet.path     = sample.sheet.path,
                       sesame.comparison            = comparison,
                       sesame.file.sep              = file.sep,
                       sesame.data.cache            = "EPIC",
                       sesame.data.normal           = "EPIC.5.normal",
                       sesame.ref.version           = "hg38",
                       sesame.reference             = reference,##"Sample_Group"
                       sesame.split.by              = split.by,
                       ...
                       )

},

################# segmentation: cnv450k analysis ############################
#############################################################################
#############################################################################
#############################################################################
#############################################################################

hm450 = {

  seg <- methyl_master_k450(k450.input.file.dir     = input.dir,
                            k450.output.file.dir    = output.dir,
                            k450.sample.sheet.path  = sample.sheet.path,
                            k450.comparison         = comparisons,
                            k450.split.by           = split.by,
                            k450.reference          = reference,
                            k450.sesame.data.cache  = "EPIC",
                            k450.sesame.data.normal = 'EPIC.5.normal',
                            k450.sesame.ref.version = "hg38",
                            ...)

}, ##End K450

############################# ChAMP ##########################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################

champ = {

  seg <- methyl_master_champ(champ.directory=input.file.dir,
                             champ.array.type="450K",
                             champ.batch.name=c("batch"),
                             champ.padj=0.05,
                             champ.ncores=n.cores,
                             champ.control=TRUE,
                             champ.control.group="normal",
                             ##champ.contol.group="champCtls"
                             champ.runimpute=TRUE,
                             champ.runQC=TRUE,
                             champ.runnorm=TRUE,
                             champ.runSVD=TRUE,
                             champ.runCombat=TRUE,
                             champ.runDMP=TRUE,
                             champ.runDMR=TRUE,
                             champ.runBlock=TRUE,
                             champ.runGSEA=TRUE,
                             ...
                             )

}, ##End ChAMP

############################# Epicopy ########################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################

epicopy = {

  seg <- methyl_master_epicopy(epi.target.dir=input.file.dir,
                               epi.output.dir=output.dir,
                               epi.single.file=TRUE,
                               epi.single.file.path=sample.sheet.path,
                               epi.run.gistic=epi.run.gistic,
                               epi.comparisons=comparisons,
                               epi.ncores=n.cores,
                               epi.ref="median",
                               epi.normals="Sample_Group",
                               epi.samp.names=NULL,
                               epi.qn=FALSE,
                               epi.mode.bw=0.1,
                               epi.mode.method="naive",
                               epi.normal.cnv=TRUE,
                               epi.mean.center=TRUE,
                               epi.filter.probes=FALSE,
                               epi.retained.probes=NULL,
                               epi.keepfnobj=TRUE,
                               epi.fn.output=NULL)

}, ##End Epicopy

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
           paste0(output.dir,
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

  source("R/methyl_master_visualize.R")

}

##################### Write Session Info #####################################

writeLines(capture.output(sessionInfo()),
           paste0(output.dir,
                  file.sep,
                  "sessionInfo.",
                  routine,
                  ".",
                  gsub("-|\\s|:",".",Sys.time(), perl = TRUE),
                  ".txt")
)

}else{
  
  methyl_master_compare(compare.input.dir = NULL,
                        compare.output.dir = NULL,
                        compare.output.name = NULL)
  
}

##################### End pipeline ###########################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################

}
