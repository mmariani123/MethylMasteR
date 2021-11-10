#!/usr/bin/env Rscript

#' @title MethylMasteR Epicopy() version
#' @description MethylMasteR Epicopy() version
#'
#' Lock updated epicopy functions in the environment:
#' load(paste0(work.dir,
#'            file.sep,
#'            "sesame_rgset_tumor.RData"))
#' load(paste0(work.dir,
#'            file.sep,
#'            "sesame_rgset_normal.RData"))
#' if(use.epi.normal==TRUE){
#'  rgnormal <- sesame::SigSetsToRGChannelSet(sesame_ssets_normal)
#' }else{
#'  load(paste0(work.dir,
#'              file.sep,
#'              "sesame_rgset_normal.RData"))
#'  rgnormal <- rgset_normal
#'  rm(rgset_normal)
#' }
#' Not run:
#' Run epicopy on included example TCGA data
#' input_loc <- system.file('extdata',
#'                         'raw_idat',
#'                         package = 'Epicopy')
#' No output directory, returns only local R object
#'
#' @param epi.target.dir
#' @param epi.output.dir
#' @param epi.ref
#' @param epi.normals
#' @param epi.samp.names
#' @param epi.qn
#' @param epi.mode.bw
#' @param epi.mode.method
#' @param epi.normal.cnv
#' @param epi.mean.center
#' @param epi.filter.probes
#' @param epi.retained.probes
#' @param epi.keepfnobj
#' @param epi.fn.output
#' @param ...
#' @import Epicopy
#' @import utils
#' @export
methyl_master_epicopy <- function(epi.target.dir=NULL,
                                  epi.output.dir=NULL,
                                  epi.single.file=TRUE,
                                  epi.single.file.path=NULL,
                                  epi.comparisons=NULL,
                                  epi.run.gistic=FALSE,
                                  epi.ncores=n.cores,
                                  epi.ref="median", #How to calculate LRR
                                  epi.less.stringent.ra=FALSE,
                                  epi.normals="Sample_Group",
                                  epi.samp.names=NULL,
                                  epi.qn=FALSE,
                                  epi.mode.bw=0.1,
                                  epi.mode.method="naive",
                                  epi.normal.cnv=TRUE,
                                  epi.mean.center=TRUE,
                                  epi.filter.probes=FALSE,
                                  epi.retained.probes=NULL,
                                  epi.keepfnobj=FALSE,
                                  epi.fn.output=NULL,
                                  ...){

#lock in modified epicopy functions to original namespace:
rlang::env_unlock(env = asNamespace('Epicopy'))
rlang::env_binding_unlock(env = asNamespace('Epicopy'))
assign('.coerce.pData', .coerce.pData.mm, envir = asNamespace('Epicopy'))
##assign('LRRtoCNA', LRRtoCNA.mm, envir = asNamespace('Epicopy'))
##assign('getLRR', getLRR.mm, envir = asNamespace('Epicopy'))
assign('epicopy', epicopy.mm, envir = asNamespace('Epicopy'))
##utils::assignInNamespace('.funnorm', .funnorm.mm,
##                         ns = asNamespace('Epicopy'))
rlang::env_binding_lock(env = asNamespace('Epicopy'))
rlang::env_lock(asNamespace('Epicopy'))

rlang::env_unlock(env = asNamespace('minfi'))
rlang::env_binding_unlock(env = asNamespace('minfi'))
utils::assignInNamespace("read.metharray", read.metharray.mm,
                         ns= asNamespace("minfi"))
utils::assignInNamespace('read.metharray.sheet', read.metharray.sheet.mm,
                         ns = asNamespace('minfi'))
rlang::env_binding_lock(env = asNamespace('minfi'))
rlang::env_lock(asNamespace('minfi'))

##source(paste0(scripts.dir,file.sep,"salas_mm_epicopy_salas.R"))

if(is.null(epi.output.dir)){
  epi.output.dir <- "."
}

epicopy_results <- Epicopy::epicopy(target_dir = epi.target.dir,
                           ##sesame_rgset_tumor,
                           ##input_loc,
                           output_dir     = epi.output.dir,
                           single.file    = epi.single.file,
                           single.file.path = epi.single.file.path,
                           comparisons    = epi.comparisons,
                           ncores         = epi.ncores,
                           Ref            = epi.ref,
                           Normals        = epi.normals,
                           ##sampNames    = "Sample_Name",
                           sampNames      = epi.samp.names,
                           QN             = epi.qn,
                           mode.bw        = epi.mode.bw,
                           mode.method    = epi.mode.method,
                           normal.cnv     = epi.normal.cnv,
                           mean.center    = epi.mean.center,
                           filterProbes   = epi.filter.probes,
                           retainedProbes = epi.retained.probes,
                           keep_fnobj     = epi.keepfnobj,
                           fn_output      = epi.fn.output,
                           ...
                           )

epicopy_results$output$ID <-
  gsub("^X","",epicopy_results$output$ID,perl = TRUE)

##intersect(rownames(tumor), rownames(normal))
##intersect(tumor$X, normal$X)
##save(epicopy_results,
##     file = paste0(work.dir,
##                   file.sep,
##                   "epicopy_example_results.RData"))

save(epicopy_results,
     file=paste0(epi.output.dir,
                 .Platform$file.sep,
                 "epicopy_results.rda")
)

write.csv(epicopy_results$output,
          file=paste0(epi.output.dir,
                      .Platform$file.sep,
                      "epicopy_results.csv"),
          row.names = FALSE,
          col.names = TRUE,
          quote = FALSE)

gistic.output.dir <-
paste0(epi.output.dir,
       .Platform$file.sep,
       "gistic_results")
message("Creating <gistic.output.dir> and running export_gistic() ...")
dir.create(gistic.output.dir)

if(epi.run.gistic==TRUE){
  Epicopy::export_gistic(epicopy_results,
                output_dir = gistic.output.dir,
                filterbycount = TRUE,
                min_probes = 50)
}

##Segmentation file written.
##Marker file written.

}
