#!/usr/bin/env Rscript

#' @title methyl_master
#' @description MethyMaster main function:
#' CNV calling platform from methylation data algorithms
#' and additional comparison functionality
#' Michael Mariani PhD Dartmouth College 2021-2022
#' Possible routines:
#' "test"            ##Run a quick test
#' "process_sesame", ##Preprocess the TCGA and cord data in sesame format
#' "sesame",         ##Run Sesame  CNV calling (get segments)
#' "hm450" ,         ##Run 450K    CNV calling (get segments)
#' "champ" ,         ##Run ChAMP   CNV calling (get segments)
#' "epicopy" ,       ##Run EpiCopy CNV calling (get segments)
#' "compare"         ##Run algorithm comparison functionality
#' @param input.dir Input directory
#' @param output.dir Output directory
#' @param sample.sheet.path Path to sample sheet
#' @param file.sep File separator
#' @param r.lib.path Path to r library default is .libPaths()
#' @param n.cores Number of cores to run
#' @param os.type The os type running the software
#' @param proj The project name
#' @param visualize Visualize the output
#' @param visualize.individual Visualize individual sample plots for sesame
#' @param simplify.reduce The reduction function being applied to the output
#' @param routine The specific routine to run
#' @param genome.version The genome version
#' @param reference The reference to use
#' @param reference.name The reference name, for the epicopy routine
#' @param split.by Which additional metadata column to split the analysis by
#' @param comparison Which groups from the metadata "Sample_Group" to compare
#' @param form.thresholds The thresholds for cnv state calling
#' @param overlap.density The overlap density applied to final results from
#' the modified population ranges function
#' @param sesame.data.cache The sesame data cache to use
#' @param sesame.data.normal Which sesame normal data to use e.g."EPIC.5.Normal"
#' @param hm450.workflow which hm450 workflow to use with the hm450 routine:
#' "A", "B", or "C"
#' @param champ.padj adjusted p-value threshold for the ChAMP routine
#' @param champ.control select the champ control to use
#' @param champ.run.combat run combat batch correction in the champ routine
#' @param champ.run.dmp run dmp in the champ routin
#' @param champ.run.dmr run dmr in the champ routine
#' @param champ.run.block run block in the champ routine
#' @param champ.run.gsea run gsea in the champ routine
#' @param champ.run.epimod run epimod in the champ routine
#' @param epi.run.gistic run GISTIC in the epicopy routine
#' @param compare.name the compare.name to use
#' @param olaps.split.field which field in the olaps to split the results by
#' the default (recommended) is to use "Sample_ID"
#' @param estimate.recurrence whether to use the estimate.recurrence feature
#' in populationsRanges when looking at overlaps within a particular routine
#' @param ov.less.stringent.ra.setting whether to use less stringent
#' method when comparing overlaps of cnvs from an individual routine
#' @param ov.keep.extra.columns whether to keep extra meta data columns or
#' not when collating results from an individual routine
#' @param ov.pvalue the pvalue cutoff to use within a routine
#' @param ov.simplify.reduce which reduction method to use to compare overlaps
#' within a routine
#' @param save.seg whether to save the seg output or not
#' @param create.dir create output directory when running MethylMaster,
#' if already exists it will be over-written
#' @param compare.list.files whether to list the files individually to be
#' compared with the compare routine or not, otherwise will search through
#' the input.dir path
#' @param compare.files.in vector of file paths to individual routine results
#' to use with the compare routine
#' @param compare.names the names of the files (with paths) to run with the
#' compare routine
#' @param compare.olaps.size the overlaps size to use when considering an
#' overlap hit across multiple routine results (using the compare routin)
#' @param ... additional parameters to be passed to methyl_master
#' @import CNVRanger
#' @import matter
#' @importFrom magrittr %>%
#' @return NULL
#' @export
methyl_master <- function(input.dir            = NULL,
                          output.dir           = NULL,
                          sample.sheet.path    = NULL,
                          file.sep             = "/",
                          r.lib.path           = .libPaths()[1],
                          n.cores              = 1,
                          os.type              = "linux",
                          proj                 = "TCGA-BLCA",
                          visualize            = FALSE,
                          visualize.individual = FALSE,
                          simplify.reduce      = weightedmean,
                          routine              = "test",
                          reference            = NULL,
                          reference.name       = "all",
                          split.by             = NULL,
                          comparison           = NULL,
                          form.thresholds      = NULL,
                          overlap.density      = 0.1,
                          sesame.data.cache    = "EPIC",
                          sesame.data.normal   = 'EPIC.5.normal',
                          genome.version       = "hg38",
                          hm450.workflow       = "B",
                          champ.padj           = 0.05,
                          champ.control        = FALSE,
                          champ.run.combat     = TRUE,
                          champ.run.dmp        = TRUE,
                          champ.run.dmr        = TRUE,
                          champ.run.block      = TRUE,
                          champ.run.gsea       = TRUE,
                          champ.run.epimod     = TRUE,
                          epi.run.gistic       = FALSE,
                          compare.name         = NULL,
                          save.seg             = FALSE,
                          olaps.split.field    = "Sample_ID",
                          estimate.recurrence  = FALSE,
                          ov.less.stringent.ra.setting =
                            ov.less.stringent.ra.setting,
                          ov.pvalue            = ov.pvalue,
                          ov.keep.extra.columns = TRUE,
                          ov.simplify.reduce    = weightedmean,
                          create.dir           = FALSE,
                          compare.list.files   = FALSE,
                          compare.files.in     = NULL,
                          compare.names        = NULL,
                          compare.olaps.size   = 1,
                          ...
                          ){

#################### Load global variables ####################################

##cat(
##paste(
##"###############################################################",
##"###############################################################",
##"###############################################################",
##"###############################################################",
##"###############################################################",
##" __  __      _   _           _ __  __           _       ____   ",
##"|  \/  | ___| |_| |__  _   _| |  \/  | __ _ ___| |_ ___|  _ \  ",
##"| |\/| |/ _ \ __| '_ \| | | | | |\/| |/ _` / __| __/ _ \ |_\) | ",
##"| |  | |  __/ |_| | | | |_| | | |  | |\(_| \__ \ ||  __/  _ <  ",
##"|_|  |_|\___|\__|_| |_|\__, |_|_|  |_|\__,_|___/\__\___|_| \_\ ",
##"                       |___/                                   ",
##"                                                               ",
##"###############################################################",
##"###############################################################",
##"###############################################################",
##"###############################################################",
##"###############################################################",
##sep="\n"))

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

##Individual plotting functionality has been moved to
##methyl_master_formatting_sesame for the time being
##if(visualize.individual==TRUE){
##  individual.plots.dir=paste0(output.dir,
##                              file.sep,
##                              "individual_plots")
##  if(!dir.exists(individual.plots.dir)){
##    dir.create(individual.plots.dir)
##  }else{
##    print("<individual.plots.dir> exists, overwriting")
##    message("<individual.plots.dir> exists, overwriting")
##    unlink(individual.plots.dir,
##           recursive=TRUE,
##           force=TRUE)
##    dir.create(individual.plots.dir)
##  }
##}

######################## Check any flags ###########################

if(reference=="internal" &
   routine=="epicopy" &
   !is.na(reference.name)){
  stop(paste0("ERROR: epicopy internal .rda object ",
              "no longer works, built ~ 2016. ",
              "must set <reference.name> to NA and epicopy ",
              "will use median, otherwise set <reference> ",
              "to 'comparison' etc."))
}

######################## Run Pipeline ##############################
####################################################################
####################################################################
####################################################################
####################################################################

if(routine!="compare"){

##profvis.out <- profvis({

##profmem.out <- matter::profmem({
peakram.out <- peakRAM::peakRAM({

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
                       sesame.genome.version        = genome.version,
                       sesame.reference             = reference,##"Sample_Group"
                       sesame.split.by              = split.by,
                       sesame.save.seg              = save.seg,
                       ...
                       )

}, ##End SeSAMe routine

################# segmentation: cnv450k analysis ############################
#############################################################################
#############################################################################
#############################################################################
#############################################################################

hm450 = {

  seg <- methyl_master_hm450(hm450.input.dir          = input.dir,
                             hm450.output.dir         = output.dir,
                             hm450.sample.sheet.path  = sample.sheet.path,
                             hm450.comparison         = comparison,
                             hm450.split.by           = split.by,
                             hm450.file.sep           = file.sep,
                             hm450.reference          = reference,
                             hm450.workflow           = hm450.workflow,
                             hm450.sesame.data.cache  = "EPIC",
                             hm450.sesame.data.normal = 'EPIC.5.normal',
                             hm450.genome.version     = genome.version,
                             hm.450.save.seg          = save.seg,
                             ...)

}, ##End HM450 routine

############################# ChAMP ##########################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################

champ = {

  seg <- methyl_master_champ(champ.input.dir=input.dir,
                             champ.output.dir=output.dir,
                             champ.sample.sheet = sample.sheet.path,
                             champ.array.type="450K",
                             champ.batch.name=c("batch"),
                             champ.padj=champ.padj,
                             champ.ncores=n.cores,
                             champ.control=champ.control,
                             champ.control.group="normal",
                             champ.comparison=comparison,
                             ##champ.contol.group="champCtls"
                             champ.runimpute=TRUE,
                             champ.runQC=TRUE,
                             champ.runnorm=TRUE,
                             champ.runSVD=TRUE,
                             champ.runCombat=champ.run.combat,
                             champ.runDMP=champ.run.dmp,
                             champ.runDMR=champ.run.dmr,
                             champ.runBlock=champ.run.block,
                             champ.runGSEA=champ.run.gsea,
                             champ.runEpiMod=champ.run.epimod,
                             ...
                             )

}, ##End ChAMP routine

############################# Epicopy ########################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################

epicopy = {

  seg <- methyl_master_epicopy(epi.input.dir=input.dir,
                               epi.output.dir=output.dir,
                               epi.single.file=TRUE,
                               epi.single.file.path=sample.sheet.path,
                               epi.reference.group = reference,
                               epi.run.gistic=epi.run.gistic,
                               epi.comparison=comparison,
                               epi.ncores=n.cores,
                               epi.ref="median",
                               epi.normals=reference.name,
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

}, ##End Epicopy routine

custom = {

  seg <- methyl_master_custom(custom.idat.files.dir = input.dir,
                              custom.output.dir        = output.dir,
                              custom.sample.sheet.path = sample.sheet.path,
                              custom.comparison        = comparison,
                              custom.file.sep          = file.sep,
                              custom.data.cache        = "EPIC",
                              custom.data.normal       = "EPIC.5.normal",
                              custom.genome.version    = genome.version,
                              custom.reference         = reference,
                              custom.save.seg          = save.seg,
                              ...
  )

} ##End custom routine

) ##End switch statement

################### TIME/MEM PROFILE ####################
#########################################################
#########################################################
#########################################################
#########################################################

##start.time <- Sys.time()
end.time <- Sys.time()

##total.time <<- as.character(end.time - start.time)
total.time <- as.character(end.time - start.time)

##}) ##End profmem

}) ##End peakRAM

##}) ##End profvis

##print(profmem.out[3])
print(total.time)
print(peakram.out$Peak_RAM_Used_MiB)

writeLines(c("\nTotal time:\n",
             total.time,
             "\nGC output: \n",
             capture.output(gc()),
             "\npeakRAM::peakRAM() max output (bytes): \n",
             peakram.out$Peak_RAM_Used_MiB,
             "Mebibytes"
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

####################### Formatting ###########################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################

if(routine=="sesame"){
  seg <-
    methyl_master_formatting_sesame(sesame.form.seg=seg,
              sesame.form.output.dir=output.dir,
              sesame.form.sample.sheet.path=sample.sheet.path,
              sesame.form.comparison=comparison,
              sesame.form.thresholds=form.thresholds,
              sesame.form.save.seg=save.seg,
              sesame.form.plot.individual=visualize.individual,
              sesame.form.auto.corrected=sesame.hm450.mean.correct,
              sesame.form.add.meta=sesame.form.add.meta,
              ...)
}else if(routine=="hm450"){
  seg <- methyl_master_formatting_hm450(hm450.form.seg=seg,
                                 hm450.form.output.dir=output.dir,
                                 hm450.form.sample.sheet.path=sample.sheet.path,
                                 hm450.form.workflow=hm450.workflow,
                                 hm450.form.thresholds=form.thresholds,
                                 hm450.form.comparison=comparison,
                                 hm450.form.save.seg=save.seg,
                                 hm450.form.anno.file.path=hm450.anno.file.path)
}else if(routine=="champ"){
  seg <- methyl_master_formatting_champ(champ.form.seg=seg,
                                        champ.form.output.dir=output.dir,
                                        champ.form.thresholds=form.thresholds,
                                        champ.form.save.seg=save.seg,
                                        champ.form.comparison=comparison,
                                        champ.form.padj=champ.padj)
}else if(routine=="epicopy"){
  seg <- methyl_master_formatting_epicopy(epi.form.seg=seg,
                                          epi.form.output.dir=output.dir,
                                          epi.form.thresholds=form.thresholds,
                                          epi.form.save.seg=save.seg,
                                          epi.form.comparison=comparison)
}}else if(routine=="custom"){
  methyl_master_formatting_custom(custom.form.seg=seg,
                                  custom.form.output.dir=output.dir,
                                  custom.form.sample.sheet.path=sample.sheet.path,
                                  custom.form.comparison=comparison,
                                  custom.form.thresholds=form.thresholds,
                                  custom.form.save.seg=save.seg,
                                  custom.form.plot.individual=visualize.individual,
                                  custom.form.auto.corrected=sesame.hm450.mean.correct,
                                  ...))
}else{
  stop(paste0("Error you must select a valid <routine> value"))
}

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

methyl_master_olaps_and_visualize(
       ov.seg                       = seg,
       ov.name                      = names(seg),
       ov.output.dir                = output.dir,
       ov.routine                   = routine,
       ov.split.field               = olaps.split.field,
       ov.keep.extra.columns        = ov.keep.extra.columns,
       ov.overlap.density           = overlap.density,
       ov.estimate.recurrence       = estimate.recurrence,
       ov.simplify.reduce           = ov.simplify.reduce,
       ov.less.stringent.ra.setting = ov.less.stringent.ra.setting,
       ov.pvalue                    = ov.pvalue,
       ov.plot.individual           = ov.plot.individual
       )

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

  methyl_master_compare(compare.list.files=compare.list.files,
                        compare.files.in=compare.files.in,
                        compare.input.dir=input.dir,
                        compare.output.dir=output.dir,
                        compare.names=compare.names,
                        compare.olaps.size=compare.olaps.size)

}

##################### End pipeline ###########################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################

}
