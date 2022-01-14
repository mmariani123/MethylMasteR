#!/usr/bin/env Rscript

#' @title methyl_master_epicopy
#' @description MethyMaster version of Epicopy::epicopy() function
#' Epicopy by:
#' Sean, Soonweng Cho, Hyunseok Kim and Leslie Cope (2017).
#' Epicopy: Get CNVinformation from 450K array. R package version 0.99.0.
#' @param epi.input.dir The idat input dir
#' @param epi.output.dir The output dir for epicopy functionality
#' @param epi.single.file Whether a single file is used
#' @param epi.single.file.path The single file path
#' @param epic.manifest.path The EPIC manifest file path
#' @param epi.reference.group The reference group to use
#' @param epi.comparison The MehtylMaster comparison vector
#' sheet
#' @param epi.ncores The number of cores to use, multiple cores may not be
#' supported on all systems
#' @param epi.ref The Epicopy ref to use
#' @param epi.normals The Epciopy normals to specify
#' @param epi.samp.names The Epicopy sample names can be custom specified here
#' @param epi.qn The Epicopy qn parameter
#' @param epi.mode.bw The Epicopy mode.bw parameter
#' @param epi.mode.method The Epicopy mode.method parameter
#' @param epi.normal.cnv The Epicopy normal.cnv parameter
#' @param epi.mean.center Whether to perform Epicopy mean centering
#' @param epi.filter.probes The Epicopy filter.probes parameter
#' @param epi.retained.probes The Epicopy retained probes parameter
#' @param epi.keepfnobj The Epicopy keepfnobj parameter
#' @param epi.fn.output The Epicopy fn.output parameter
#' @param epi.run.gistic The Epicopy run.gistic parameter
#' (Whetehr to run GISTIC)
#' @import Epicopy
#' @return epicopy results data frame
#' @export
methyl_master_epicopy <- function(epi.input.dir,
                                  epi.output.dir,
                                  epi.single.file,
                                  epi.single.file.path,
                                  epic.manifest.path,
                                  epi.reference.group,
                                  epi.comparison,
                                  epi.ncores,
                                  epi.ref,
                                  epi.normals,
                                  ##sampNames = "Sample_Name",
                                  epi.samp.names,
                                  epi.qn,
                                  epi.mode.bw,
                                  epi.mode.method,
                                  epi.normal.cnv,
                                  epi.mean.center,
                                  epi.filter.probes,
                                  epi.retained.probes,
                                  epi.keepfnobj,
                                  epi.fn.output,
                                  epi.run.gistic=FALSE,
                                  ...
                                  ){

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
                         ns = asNamespace("minfi"))
utils::assignInNamespace('read.metharray.sheet', read.metharray.sheet.mm,
                         ns = asNamespace('minfi'))
rlang::env_binding_lock(env = asNamespace('minfi'))
rlang::env_lock(asNamespace('minfi'))

epicopy_results <- Epicopy::epicopy(target_dir = epi.input.dir,
                           output_dir          = epi.output.dir,
                           single.file         = epi.single.file,
                           single.file.path    = epi.single.file.path,
                           epic.manifest.path  = epic.manifest.path,
                           reference_group     = epi.reference.group,
                           comparison          = epi.comparison,
                           ncores              = epi.ncores,
                           Ref                 = epi.ref,
                           Normals             = epi.normals,
                           ##sampNames         = "Sample_Name",
                           sampNames           = epi.samp.names,
                           QN                  = epi.qn,
                           mode.bw             = epi.mode.bw,
                           mode.method         = epi.mode.method,
                           normal.cnv          = epi.normal.cnv,
                           mean.center         = epi.mean.center,
                           filterProbes        = epi.filter.probes,
                           retainedProbes      = epi.retained.probes,
                           keep_fnobj          = epi.keepfnobj,
                           fn_output           = epi.fn.output
                           )

  if(epi.run.gistic==TRUE){
      gistic.output.dir <-
        paste0(epi.output.dir,
              .Platform$file.sep,
              "gistic_results")
  if(!dir.exists(gistic.output.dir)){
      dir.create(gistic.output.dir)
  }else{
    print("<gistic.output.dir> exists, overwriting")
    message("<gistic.output.dir> exists, overwriting")
    unlink(gistic.output.dir,
           recursive=TRUE,
           force=TRUE)
    dir.create(gistic.output.dir)
  }
    message("Running export_gistic() ...")
    Epicopy::export_gistic(epicopy_results[[1]],
                          output_dir = gistic.output.dir,
                          filterbycount = TRUE,
                          min_probes = 50)
  }

gistic.df <- read.table(
  file=paste0(gistic.output.dir,.Platform$file.sep,"segmentation_output.txt"),
                        header=TRUE,
                        sep="\t",
                        stringsAsFactors = FALSE)

seg <- list(gistic.df)

return(seg)

}
