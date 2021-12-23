#!/usr/bin/env Rscript

#' @title methyl_master_hm450
#' @description My version of hm450 analyses
#' ############ Samples and ref also need to be in same order ##########
#' "CpG probe IDs not in the same order in data and ctrl!"
#' Modified by Michael Mariani PhD Dartmouth College 2021
#' @param champ.input.dir
#' @param champ.output.dir
#' @param champ.sample.sheet
#' @param champ.array.type
#' @param champ.batch.name
#' @param champ.padj
#' @param champ.ncores
#' @param champ.control
#' @param champ.control.group
#' @param champ.comparison
#' @param champ.runimpute
#' @param champ.runQC
#' @param champ.runnorm
#' @param champ.runSVD
#' @param champ.runCombat
#' @param champ.runDMP
#' @param champ.runDMR
#' @param champ.runBlock
#' @param champ.runGSEA
#' @param champ.runEpiMod
#' @param ...
#' @import ChAMP
#' @import org.Hs.eg.db
#' @return #champ results
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
