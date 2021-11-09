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
#' @import Epicopy
#' @import utils
#' @export
methyl_master_epicopy <- function(epi.target.dir=idat.files.dir,
                                  epi.output.dir=NULL,
                                  epi.single.file=TRUE,
                                  epi.single.file.path=NULL,
                                  epi.comparisons=NULL,
                                  epi.ncores=1,
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
                                  epi.keepfnobj=TRUE,
                                  epi.fn.output=NULL){

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
  epi.output.dir <- output.dir
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
                           fn_output      = epi.fn.output)

##epicopy_results$output$Sample_ID %>% unique() %in% normal$X
##epicopy_results$output$Sample_ID %>% unique() %in% tumor$X
##Need to remove the "X" from the beginning of the sample names
##epicopy_results_fix <- epicopy_results
##epicopy_results_fix$output$Sample_ID <-
##gsub("^X","",epicopy_results$output$Sample_ID,perl = TRUE)
##epicopy_results_fix$output$Sample_ID %>% unique() %in% normal$X
##epicopy_results_fix$output$Sample_ID %>% unique() %in% tumor$X

epicopy_results$output$Sample_ID <-
  gsub("^X","",epicopy_results$output$Sample_ID,perl = TRUE)

##intersect(rownames(tumor), rownames(normal))
##intersect(tumor$X, normal$X)
##save(epicopy_results,
##     file = paste0(work.dir,
##                   file.sep,
##                   "epicopy_example_results.RData"))

## Assign sex information:
##colnames(epicopy_results$output)
epicopy_results$output$karyotype    <- ""
epicopy_results$output$sex_reported <- ""
epicopy_results$output$sex_inferred <- ""
epicopy_results$output$treatment    <- ""
for(i in 1:length(unique(epicopy_results$output$Sample_ID))){
  ##print(i)
  sample.now <- unique(epicopy_results$output$Sample_ID)[i]
  if(sample.now %in% rownames(tumor)){
    epicopy_results$output[epicopy_results$output$Sample_ID==sample.now,]$karyotype    <-
      tumor[tumor$X==sample.now, "karyotype"]
    epicopy_results$output[epicopy_results$output$Sample_ID==sample.now,]$sex_reported <-
      tumor[tumor$X==sample.now, "gender_reported"]
    epicopy_results$output[epicopy_results$output$Sample_ID==sample.now,]$sex_inferred  <-
      tumor[tumor$X==sample.now, "sex_inferred"]
    epicopy_results$output[epicopy_results$output$Sample_ID==sample.now,]$treatment <- "tumor"
  }else if(sample.now %in% rownames(normal)){
    epicopy_results$output[epicopy_results$output$Sample_ID==sample.now,]$karyotype    <-
      normal[normal$X==sample.now, "karyotype"]
    epicopy_results$output[epicopy_results$output$Sample_ID==sample.now,]$sex_reported <-
      normal[normal$X==sample.now, "gender_reported"]
    epicopy_results$output[epicopy_results$output$Sample_ID==sample.now,]$sex_inferred <-
      normal[normal$X==sample.now, "sex_inferred"]
    epicopy_results$output[epicopy_results$output$Sample_ID==sample.now,]$treatment <- "normal"
  }else{
    stop(paste0("Error: epicopy processing:",
                "samples names must belong ",
                "to either tumor or normal ",
                "treatments"))
  }
}

save(epicopy_results,
     file=paste0(work.dir,
                 file.sep,
                 "epicopy_results.rda")
)

write.csv(epicopy_results$output,
          file=paste0(work.dir,
                      file.sep,
                      "epicopy_results.csv"),
          row.names = FALSE,
          col.names = TRUE,
          quote = FALSE)

if(run.gistic==TRUE){
  export_gistic(epicopy_results,
                filterbycount = TRUE,
                min_probes = 50,
                output_dir = work.dir)
}

##Segmentation file written.
##Marker file written.

}
