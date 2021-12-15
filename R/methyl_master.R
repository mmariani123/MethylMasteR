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
#' @param sesame.data.cache
#' @param sesame.data.normal
#' @param sesame.ref.version
#' @param epi.run.gistic
#' @param compare.name
#' @param olaps.split.field
#' @param estimate.recurrence
#' @param ov.less.stringent.ra.setting
#' @param ov.keep.extra.columns
#' @param ov.pvalue
#' @param ov.simplify.reduce
#' @param save.seg
#' @param create.dir
#' @param ...
#' @import CNVRanger
#' @import matter
#' @importFrom magrittr %>%
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
                          sesame.data.cache    = "EPIC",
                          sesame.data.normal   = 'EPIC.5.normal',
                          sesame.ref.version   = "hg38",
                          epi.run.gistic       = FALSE,
                          compare.name         = NULL,
                          save.seg             = FALSE,
                          olaps.split.field    = "Sample_ID",
                          estimate.recurrence  = FALSE,
                          ov.less.stringent.ra.setting = ov.less.stringent.ra.setting,
                          ov.pvalue            = ov.pvalue,
                          ov.keep.extra.columns = TRUE,
                          ov.simplify.reduce  = weightedmean,
                          create.dir           = FALSE,
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

if(create.dir==TRUE){
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

profmem.out <- matter::profmem({

start.time <- Sys.time()

switch(routine,

######################## Run Test ##################################
####################################################################
####################################################################
####################################################################
####################################################################

test = {

  print("Hello World")
  ##profvis::pause(10)

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
                       sesame.save.seg              = save.seg,
                       ...
                       )

},

################# segmentation: cnv450k analysis ############################
#############################################################################
#############################################################################
#############################################################################
#############################################################################

hm450 = {

  seg <- methyl_master_hm450(hm450.input.file.dir     = input.dir,
                             hm450.output.file.dir    = output.dir,
                             hm450.sample.sheet.path  = sample.sheet.path,
                             hm450.comparison         = comparison,
                             hm450.split.by           = split.by,
                             hm450.reference          = reference,
                             hm450.sesame.data.cache  = "EPIC",
                             hm450.sesame.data.normal = 'EPIC.5.normal',
                             hm450.sesame.ref.version = "hg38",
                             hm.450.save.seg          = save.seg,
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
                             champ.save.seg=save.seg,
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
                               epi.comparison=comparison,
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
                               epi.fn.output=NULL,
                               epi.save.seg=save.seg,
                               ...
                               )

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
             profmem.out[3]
             ##{profmem.out$bytes[!is.na(profmem.out$bytes)] %>%
             ##     max()}
             ),
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

## 0: homozygous deletion (2-copy loss)
## 1: heterozygous deletion (1-copy loss)
## 2: normal diploid state
## 3: 1-copy gain
## 4: amplification (>= 2-copy gain)

##Overall the pipeline results structure is
##a dataframe with the following fields:
##"Sample_ID"
##"chrom"
##"loc.start"
##"loc.end"
##"num.mark"
##"seg.mean"
##"state"
##"method"

################# OVERLAPS AND VISUALIZE ####################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################

##samples.list <- as.list(c(paste0("normal_",seq(1,10)),
##                            paste0("tumor_",seq(1,10))))

##seg$chrom <- as.numeric(seg$chrom)
##seg <- na.omit(seg)
##changing the chrom to number and removing NAs from X and Y doesn't seem
##to fix the population ranges problem.

#################### VISUALISE INDIVIDUAL ######################################
##load(paste0(output.dir,.Platform$file.sep,"seg.RData"))
mapply(x=seg,y=names(seg),
  function(x=x,
           y=y,
           output.dir1=output.dir,
           routine1=routine,
           olaps.split.field1=olaps.split.field,
           ov.keep.extra.columns1=ov.keep.extra.columns,
           overlap.density1=overlap.density,
           estimate.recurrence1=estimate.recurrence,
           ov.simplify.reduce1=ov.simplify.reduce,
           ov.less.stringent.ra.setting1=ov.less.stringent.ra.setting,
           ov.pvalue1=ov.pvalue){
    methyl_master_olaps_and_visualize(
       ov.seg                       = x,
       ov.name                      = y,
       ov.output.dir                = output.dir1,
       ov.routine                   = routine1,
       ov.split.field               = olaps.split.field1,
       ov.keep.extra.columns        = ov.keep.extra.columns1,
       ov.overlap.density           = overlap.density1,
       ov.estimate.recurrence       = estimate.recurrence1,
       ov.simplify.reduce           = ov.simplify.reduce1,
       ov.less.stringent.ra.setting = ov.less.stringent.ra.setting1,
       ov.pvalue                    = ov.pvalue1
       )})

if(visualize.individual==TRUE){
  ##load(paste0(output.dir,.Platform$file.sep,"seg.RData"))
  mapply(x=seg,y=names(seg),
         function(x=x,
                  y=y,
                  output.dir1=output.dir
                  ){
                 methyl_master_plot_individual(
                   pi.seg=x,
                   pi.name=y,
                   pi.output.dir=output.dir1
                 )})

}##End visualize individual

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
