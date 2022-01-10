#!/usr/bin/env Rscript

#' @title methyl_master_champ
#' @description My function to call the ChAMP::champ.process() function
#' For ChAMP refs see:
#' Tian Y, Morris TJ, Webster AP, Yang Z, Beck S, Andrew F, Teschendorff AE
#' (2017). "ChAMP: updated methylation analysis pipeline for Illumina
#' BeadChips." _Bioinformatics_, btx513. doi: 10.1093/bioinformatics/btx513
#' (URL: https://doi.org/10.1093/bioinformatics/btx513).
#'
#' Morris TJ, Butcher LM, Teschendorff AE, Chakravarthy AR, Wojdacz TK, Beck S
#' (2014). "ChAMP: 450k Chip Analysis Methylation Pipeline." _Bioinformatics_,
#' *30*(3), 428-430. doi: 10.1093/bioinformatics/btt684 (URL:
#' https://doi.org/10.1093/bioinformatics/btt684).
#'
#' champ.lasso method is described in:
#'
#'   Butcher LM, Beck S (2015). "Probe Lasso: A novel method to rope in
#' differentially methylated regions with 450K DNA methylation data."
#' _Methods_, *72*, 21-28. doi: 10.1016%2Fj.ymeth.2014.10.036 (URL:
#' https://doi.org/10.1016%2Fj.ymeth.2014.10.036).
#' @param champ.input.dir The input idat files dir
#' @param champ.output.dir The output dir for ChAMP output files
#' @param champ.sample.sheet The MethylMaster sample sheet path
#' @param champ.array.type The array type being used, default is "hm450"
#' @param champ.batch.name The field in the sample sheet used as batch field
#' for batch correction with combat, default is c("Batch")
#' @param champ.padj The p.adj value to filter ChAMP result by (default is 0.05)
#' @param champ.ncores The number of cores to use with ChAMP, note that more
#' than 1, i.e. parallel processing, may not work on all systems.
#' @param champ.control The champ control to use
#' @param champ.control.group The specific ChAMP control group
#' @param champ.comparison The MethylMaster comparison 2 element vector
#' @param champ.runimpute Whether to run Champ.runimpute or not
#' @param champ.runQC Whether to run champ.runQC or not
#' @param champ.runnorm Whether to run champ.runnorm or not
#' @param champ.runSVD Whether to run champ.runSVD or not
#' @param champ.runCombat Whether to run combat in champ (champ.runCombat)
#' @param champ.runDMP Whether to run champ.runDMP
#' @param champ.runDMR Whether to run champ.DMR
#' @param champ.runBlock Whether to run champ.runBlock
#' @param champ.runGSEA Whether to run champ.runGSEA
#' @param champ.runEpiMod Whether to run champ.runEpiMod
#' @param ... Additional parameter to pass to ChAMP
#' @import ChAMP
#' @import org.Hs.eg.db
#' @return ChAMP results stored in a list
#' @export
methyl_master_champ <- function(champ.input.dir=getwd(),
                                champ.output.dir=getwd(),
                                champ.sample.sheet=NULL,
                                champ.array.type="hm450",
                                champ.batch.name=c("batch"),
                                champ.padj=0.05,
                                champ.ncores=1,
                                champ.control=TRUE,
                                champ.control.group="normal",
                                champ.comparison=c("tumor","normal"),
                                champ.split.by=NULL,
                                champ.runimpute=TRUE,
                                champ.runQC=TRUE,
                                champ.runnorm=TRUE,
                                champ.runSVD=TRUE,
                                champ.runCombat=TRUE,
                                champ.runDMP=TRUE,
                                champ.runDMR=TRUE,
                                champ.runBlock=TRUE,
                                champ.runGSEA=TRUE,
                                champ.runEpiMod=TRUE,
                                ...){

##library(org.Hs.eg.db) ##loads AnnotationDbi, contains org.Hs.egSYMBOL

##source(paste0("G:\\My Drive\\dartmouth",
##              "\\salas_lab_working\\cnv",
##              "\\cnv_testing\\files",
##              "\\methyl_master_champ_helper_functions.R"))

##ChAMP.EpiMod START requires below or throws error
##load(file=paste0("G:\\My Drive\\dartmouth",
##                 "\\salas_lab_working\\cnv",
##                 "\\cnv_testing\\probe450kfemanno.rda"))

##lock in modified epicopy functions to original namespace:
rlang::env_unlock(env = asNamespace('ChAMP'))
rlang::env_binding_unlock(env = asNamespace('ChAMP'))
assign('champ.import',  champ.import.mm,  envir = asNamespace('ChAMP'))
assign('champ.load',    champ.load.mm,    envir = asNamespace('ChAMP'))
assign('champ.process', champ.process.mm, envir = asNamespace('ChAMP'))
##assign('LRRtoCNA', LRRtoCNA.mm, envir = asNamespace('Epicopy'))
##assign('getLRR', getLRR.mm, envir = asNamespace('Epicopy'))
##assign('epicopy', epicopy.mm, envir = asNamespace('Epicopy'))
##utils::assignInNamespace('.funnorm', .funnorm.mm,
##                         ns = asNamespace('Epicopy'))
rlang::env_binding_lock(env = asNamespace('ChAMP'))
rlang::env_lock(asNamespace('ChAMP'))

##rlang::env_unlock(env = asNamespace('minfi'))
##rlang::env_binding_unlock(env = asNamespace('minfi'))
##utils::assignInNamespace("read.metharray", read.metharray.mm,
##                         ns= asNamespace("minfi"))
##utils::assignInNamespace('read.metharray.sheet', read.metharray.sheet.mm,
##                         ns = asNamespace('minfi'))
##rlang::env_binding_lock(env = asNamespace('minfi'))
##rlang::env_lock(asNamespace('minfi'))

##print(champ.sample.sheet)

champ.results <- ChAMP::champ.process(runload       = TRUE,
                                      directory     = champ.input.dir,
                                      resultsDir    = champ.output.dir,
                                      sample.sheet  = champ.sample.sheet,
                                      arraytype     = champ.array.type,
                                      ##"450K"
                                      batchname     = champ.batch.name,
                                      adjPVal       = champ.padj,
                                      cores         = champ.ncores,
                                      control       = champ.control,
                                      controlGroup  = champ.control.group,
                                      compare.group = champ.comparison,
                                      ##champ.split.by = champ.split.by,
                                      runimpute     = champ.runimpute,
                                      runQC         = champ.runQC,
                                      runnorm       = champ.runnorm,
                                      runSVD        = champ.runSVD,
                                      runCombat     = champ.runCombat,
                                      runDMP        = champ.runDMP,
                                      runDMR        = champ.runDMR,
                                      runBlock      = champ.runBlock,
                                      runGSEA       = champ.runGSEA,
                                      runEpiMod     = champ.runEpiMod,
                                      ...)

seg <- list(champ.results)

return(seg)

}
