#!/usr/bin/env Rscript

#' @title ##methyl_master_Sesame
#' @description MethylMasteR run SeSAMe function
#'
#' First part:: Very similar to PROCESS SESAME
#' command to prepare for segementation
#' this is included so that the total
#' time and mem for running sesame
#' can be recorded
#'
#' @param sesame.idat.files.dir
#' @param sesame.output.dir
#' @param sesame.sample.sheet.path
#' @param sesame.comparison
#' @param sesame.file.sep
#' @param sesame.data.cache
#' @param sesame.data.normal
#' @param sesame.ref.version
#' @param sesame.reference
#' @param sesame.split.by
#' @param sesame.save.seg
#' @param ...
#' @import data.table
#' @import dplyr
#' @import sesame
#' @import sesameData
#' @import ExperimentHub
#' @import future
#' @import profvis
#' @return A seg object for downstream analysis
#' @export

methyl_master_sesame <- function(sesame.idat.files.dir       = NULL,
                                 sesame.output.dir           = NULL,
                                 sesame.ref                  = NULL,
                                 sesame.sample.sheet.path    = NULL,
                                 sesame.comparison           = NULL,
                                 sesame.file.sep             = NULL,
                                 sesame.data.cache           = "EPIC",
                                 sesame.data.normal          = "EPIC.5.normal",
                                 sesame.ref.version          = "hg38",
                                 sesame.reference            = "internal",
                                 sesame.split.by             = NULL,
                                 sesame.save.seg             = FALSE,
                                 ...
                      ){

##Get the samples from the sample sheet:
sesame.sample.sheet.df <- read.csv(sesame.sample.sheet.path,
                                   header=TRUE,
                                   stringsAsFactors = FALSE)

if(sesame.reference=="internal"){

  if(is.null(sesame.split.by)){

    ##----------------------------------------------------------
    ##| SEnsible Step-wise Analysis of DNA MEthylation (SeSAMe)
    ##| --------------------------------------------------------
    ##| Please cache the annotation data for your array platform
    ##| (e.g. EPIC) by calling "sesameDataCache("EPIC")"
    ##| or "sesameDataCacheAll()". This needs to be done only
    ##| once per SeSAMe installation.

    ##Get sesame normal samples:
    ##Loads all refs as <sesame.data.cache> defaults to "EPIC"?
    sesameData::sesameDataCache(sesame.data.cache)
    sesame_ssets_normal <- sesameData::sesameDataGet(sesame.data.normal)
    ExperimentHub::setExperimentHubOption("CACHE", idat.files.dir)
    ExperimentHub::ExperimentHub()
    treatment_idat_prefixes <- sesame::searchIDATprefixes(idat.files.dir,
                                                        recursive=TRUE)
    sesameData::sesameDataCacheAll()

    treatment.names <- sesame.sample.sheet.df[
      sesame.sample.sheet.df$Sample_Group %in%
        sesame.comparison[1],"Sample_Name"]

    treatment_idat_prefixes <-  treatment_idat_prefixes[gsub(paste0(".*",
                                     .Platform$file.sep),
                                     "",
                                     treatment_idat_prefixes) %in%
                                     treatment.names]

    treatment.platform <- unique(sesame.sample.sheet.df[
                                 sesame.sample.sheet.df$Sample_Name %in%
                                 treatment.names, "Platform"])

    sesame_sset <- sesame::openSesame(treatment_idat_prefixes,
                            mask = TRUE,
                            sum.TypeI = TRUE,
                            platform = treatment.platform,
                            what="sigset")

    ##sesame_seg %<-% foreach(i = 1:length(names(sesame_sset))) %do% {
    sesame_seg <- foreach(i = 1:length(names(sesame_sset))) %do% {
      sesame::cnSegmentation(sesame_sset[[i]],
                            sesame_ssets_normal,
                            refversion = sesame.ref.version)
    }
    names(sesame_seg) <- names(sesame_sset)

    seg <- list(sesame_seg)

  }else{

    stop(paste0("Error: Epic.5.Normal internal reference samples are all male",
                "\ncan't split by sex!"))
    ##Get sesame normal samples:

    split.by.cat <- unique(sesame.sample.sheet.df[[sesame.split.by]])

    sesameData::sesameDataCache(sesame.data.cache)
    sesame_ssets_normal <- sesameData::sesameDataGet(sesame.data.normal)
    normal.sexes <- unlist(lapply(sesame_ssets_normal,sesame::inferSex))
    ssesame_ssets.normal.1 <-
      sesame_ssets_normal[names(which(normal.sexes==split.by.cat[1]))]
    ssesame_ssets.normal.2 <-
      sesame_ssets_normal[names(which(normal.sexes==split.by.cat[2]))]

    ExperimentHub::setExperimentHubOption("CACHE", idat.files.dir)

    ExperimentHub::ExperimentHub()

    idat_prefixes <- sesame::searchIDATprefixes(idat.files.dir,
                                                          recursive=TRUE)

    sesameData::sesameDataCacheAll()

    treatment.samples.1 <-
      sesame.sample.sheet.df[
        sesame.sample.sheet.df$Sample_Group==sesame.comparison[1] &
          sesame.sample.sheet.df[[sesame.split.by]]==split.by.cat[1],
        "Sample_Name"]

    treatment.samples.2 <-
      sesame.sample.sheet.df[
        sesame.sample.sheet.df$Sample_Group==sesame.comparison[1] &
          sesame.sample.sheet.df[[sesame.split.by]]==split.by.cat[2],
        "Sample_Name"]

    idat_prefixes.treatment.1 <-
      idat_prefixes[gsub(".*/.*_","",idat_prefixes) %in%
                      gsub(paste0(".*",sesame.file.sep,".*_"),
                           "",treatment.paths.1)]

    idat_prefixes.treatment.2 <-
      idat_prefixes[gsub(".*/.*_","",idat_prefixes) %in%
                      gsub(paste0(".*",sesame.file.sep,".*_"),
                           "",treatment.paths.2)]

    sesame_sset.1 <- sesame::openSesame(treatment_idat_prefixes.1,
                                      mask = TRUE,
                                      sum.TypeI = TRUE,
                                      platform = sesame.platform,
                                      what="sigset")

    sesame_seg.1 <- foreach(i = 1:length(names(sesame_sset.1))) %do% {
      sesame::cnSegmentation(sesame_sset.1[[i]],
                             sesame_ssets_normal.1,
                             refversion = sesame.ref.version)
    }
    names(sesame_seg.1) <- names(sesame_sset.1)

    sesame_sset.2 <- sesame::openSesame(treatment_idat_prefixes.2,
                                        mask = TRUE,
                                        sum.TypeI = TRUE,
                                        platform = sesame.platform,
                                        what="sigset")

    sesame_seg.2 <- foreach(i = 1:length(names(sesame_sset.2))) %do% {
      sesame::cnSegmentation(sesame_sset.2[[i]],
                             sesame_ssets_normal.2,
                             refversion = sesame.ref.version)
    }
    names(sesame_seg.2) <- names(sesame_sset.2)

    seg <- list(sesame_seg.1, sesame_seg.2)

  }

}else if(sesame.reference=="comparison"){

  if(is.null(sesame.split.by)){

    treatment.samples <-
      sesame.sample.sheet.df[
        sesame.sample.sheet.df$Sample_Group==sesame.comparison[1],
        "Sample_Name"]

    control.samples <-
        sesame.sample.sheet.df[
        sesame.sample.sheet.df$Sample_Group==sesame.comparison[2],
        "Sample_Name"]

      treatment.paths <- sesame.sample.sheet.df[
      sesame.sample.sheet.df$Sample_Group==sesame.comparison[1],
      "Basename"]

    control.paths   <-  sesame.sample.sheet.df[
      sesame.sample.sheet.df$Sample_Group==sesame.comparison[2],
      "Basename"]

    treatment.platform <-     unique(sesame.sample.sheet.df[
      sesame.sample.sheet.df$Sample_Group==sesame.comparison[1],
      "Platform"])

    control.platform   <-     unique(sesame.sample.sheet.df[
      sesame.sample.sheet.df$Sample_Group==sesame.comparison[2],
      "Platform"])

    ExperimentHub::setExperimentHubOption("CACHE",
                         sesame.idat.files.dir)

    ExperimentHub::ExperimentHub()

    idat_prefixes <- sesame::searchIDATprefixes(sesame.idat.files.dir,
                                              recursive=TRUE)

    idat_prefixes.treatment <-
      idat_prefixes[gsub(".*/","",idat_prefixes) %in%
                      gsub(paste0(".*",sesame.file.sep),
                           "",treatment.paths)]

    idat_prefixes.control   <-
      idat_prefixes[gsub(".*/","",idat_prefixes) %in%
                      gsub(paste0(".*",sesame.file.sep),
                           "",control.paths)]

    sesameData::sesameDataCacheAll()

    sesame_sset.treatment <- sesame::openSesame(idat_prefixes.treatment,
                            mask = TRUE,
                            sum.TypeI = TRUE,
                            platform = treatment.platform,
                            what="sigset")

    sesame_sset.control <- sesame::openSesame(idat_prefixes.control,
                            mask = TRUE,
                            sum.TypeI = TRUE,
                            platform = control.platform,
                            what="sigset")

    sesame_seg <- foreach(i = 1:length(names(sesame_sset.treatment))) %do% {
      sesame::cnSegmentation(sesame_sset.treatment[[i]],
                             sesame_sset.control,
                             refversion = sesame.ref.version)
    }
    names(sesame_seg) <- names(sesame_sset.treatment)

    seg <- list(sesame_seg)

  }else{

    split.by.cat <- unique(sesame.sample.sheet.df[[sesame.split.by]])

    treatment.samples.1 <-
      sesame.sample.sheet.df[
      sesame.sample.sheet.df$Sample_Group==sesame.comparison[1] &
      sesame.sample.sheet.df[[sesame.split.by]]==split.by.cat[1],"Sample_Name"]

    control.samples.1 <-
      sesame.sample.sheet.df[
      sesame.sample.sheet.df$Sample_Group==sesame.comparison[2] &
      sesame.sample.sheet.df[[sesame.split.by]]==split.by.cat[1],"Sample_Name"]

    treatment.paths.1 <-
      sesame.sample.sheet.df[
        sesame.sample.sheet.df$Sample_Group==sesame.comparison[1] &
        sesame.sample.sheet.df[[sesame.split.by]]==split.by.cat[1],"Basename"]

    control.paths.1   <-
      sesame.sample.sheet.df[
        sesame.sample.sheet.df$Sample_Group==sesame.comparison[2] &
        sesame.sample.sheet.df[[sesame.split.by]]==split.by.cat[1],"Basename"]

    treatment.platform.1 <-
      unique(sesame.sample.sheet.df[
        sesame.sample.sheet.df$Sample_Group==sesame.comparison[1] &
        sesame.sample.sheet.df[[sesame.split.by]]==split.by.cat[2],"Platform"])

    control.platform.1 <-
      unique(sesame.sample.sheet.df[
        sesame.sample.sheet.df$Sample_Group==sesame.comparison[2] &
        sesame.sample.sheet.df[[sesame.split.by]]==split.by.cat[2],"Platform"])

    treatment.samples.2 <-
      sesame.sample.sheet.df[
      sesame.sample.sheet.df$Sample_Group==sesame.comparison[1] &
      sesame.sample.sheet.df[[sesame.split.by]]==split.by.cat[2],"Sample_Name"]

    control.samples.2 <-
      sesame.sample.sheet.df[
      sesame.sample.sheet.df$Sample_Group==sesame.comparison[2] &
      sesame.sample.sheet.df[[sesame.split.by]]==split.by.cat[2],"Sample_Name"]

    treatment.paths.2 <-
      sesame.sample.sheet.df[
        sesame.sample.sheet.df$Sample_Group==sesame.comparison[1] &
        sesame.sample.sheet.df[[sesame.split.by]]==split.by.cat[2],"Basename"]

    control.paths.2   <-
      sesame.sample.sheet.df[
        sesame.sample.sheet.df$Sample_Group==sesame.comparison[2] &
        sesame.sample.sheet.df[[sesame.split.by]]==split.by.cat[2],"Basename"]

    treatment.platform.2 <-
      unique(sesame.sample.sheet.df[
        sesame.sample.sheet.df$Sample_Group==sesame.comparison[1] &
        sesame.sample.sheet.df[[sesame.split.by]]==split.by.cat[2],"Platform"])

    control.platform.2   <-
      unique(sesame.sample.sheet.df[
        sesame.sample.sheet.df$Sample_Group==sesame.comparison[2] &
        sesame.sample.sheet.df[[sesame.split.by]]==split.by.cat[2],"Platform"])

    ExperimentHub::setExperimentHubOption("CACHE",
                           sesame.idat.files.dir)

    ExperimentHub::ExperimentHub()

    idat_prefixes <- sesame::searchIDATprefixes(sesame.idat.files.dir,
                                                recursive=TRUE)

    idat_prefixes.treatment.1 <-
      idat_prefixes[gsub(".*/","",idat_prefixes) %in%
                    gsub(paste0(".*",sesame.file.sep),
                         "",treatment.paths.1)]

    idat_prefixes.control.1   <-
      idat_prefixes[gsub(".*/","",idat_prefixes) %in%
                    gsub(paste0(".*",sesame.file.sep),
                         "",control.paths.1)]

    idat_prefixes.treatment.2 <-
      idat_prefixes[gsub(".*/","",idat_prefixes) %in%
                    gsub(paste0(".*",sesame.file.sep),
                         "",treatment.paths.2)]

    idat_prefixes.control.2   <-
      idat_prefixes[gsub(".*/","",idat_prefixes) %in%
                    gsub(paste0(".*",sesame.file.sep),
                         "",control.paths.2)]

    sesameData::sesameDataCacheAll()

    ##future::plan("multisession")

    ##sesame_sset.treatment.1 %<-%
    sesame_sset.treatment.1 <- sesame::openSesame(idat_prefixes.treatment.1,
                                      mask = TRUE,
                                      sum.TypeI = TRUE,
                                      platform = treatment.platform.1,
                                      what="sigset")

    ##sesame_sset.control.1 %<-%
    sesame_sset.control.1 <- sesame::openSesame(idat_prefixes.control.1,
                                        mask = TRUE,
                                        sum.TypeI = TRUE,
                                        platform = control.platform.1,
                                        what="sigset")

    ##sesame_sset.treatment.2 %<-%
      sesame_sset.treatment.2 <-
      sesame::openSesame(idat_prefixes.treatment.2,
      mask = TRUE,
      sum.TypeI = TRUE,
      platform = treatment.platform.2,
      what="sigset")

    ##sesame_sset.control.2 %<-%
      sesame_sset.control.2 <-
      sesame::openSesame(idat_prefixes.control.2,
      mask = TRUE,
      sum.TypeI = TRUE,
      platform = control.platform.2,
      what="sigset")

    ##while(!all(resolved(x))){
    ##  profvis::pause(5)
    ##}

    ##sesame_seg.treatment.1 %<-%
      sesame_seg.treatment.1 <-
      foreach(i = 1:length(names(sesame_sset.treatment.1))) %do% {
      sesame::cnSegmentation(sesame_sset.treatment.1[[i]],
                             sesame_sset.control.1,
                             refversion = sesame.ref.version)
    }
    names(sesame_seg.treatment.1) <- names(sesame_sset.treatment.1)

    ##sesame_seg.treatment.2 %<-%
    sesame_seg.treatment.2 <-
      foreach(i = 1:length(names(sesame_sset.treatment.2))) %do% {
      sesame::cnSegmentation(sesame_sset.treatment.2[[i]],
                             sesame_sset.control.2,
                             refversion = sesame.ref.version)
    }
    names(sesame_seg.treatment.2) <- names(sesame_sset.treatment.2)

    seg <- list(sesame_seg.treatment.1, sesame_seg.treatment.2)

    ##load("C:\\Users\\Mike\\Desktop\\cnv_testing\\sesame_normal_sex\\seg.RData")

  } ##End split.by

} ##End reference

} ##End methyl_master_sesame

###################### Possible additional functionality ######################

##if(sesame.calc.betas==TRUE){
##  sesame_betas <- sesame::openSesame(idat_prefixes,
##                                     mask = TRUE,
##                                     sum.TypeI = TRUE,
##                                     platform = sesame.platform)
##}
##
##if(sesame.infer.sex.karyotypes==TRUE){
##  sesame.karyotype <- foreach(i = 1:length(names(sesame_sset))) %do% {
##    sesame::inferSexKaryotypes(sesame_sset[[i]])
##  }
##  names(sesame.karyotype) <- names(sesame_sset)
##  sesame.karyotype <- as.data.frame(unlist(sesame_karyotype))
##  colnames(sesame.karyotype) <- "karyotype"
##  data.table::fwrite(sesame.karyotype,
##                     file = sesame.output.dir,
##                     col.names = TRUE,
##                     row.names = FALSE,
##                     quote = FALSE)
##}
##
##if(sesame.calc.qc==TRUE){
##  sesame_qc[[i]] <- foreach(i = 1:length(names(sesame_sset))) %do% {
##    sesame::sesameQC(sesame_sset[[i]])
##  }
##  names(sesame_qc) <- names(sesame_sset)
##  sesame_qc <- as.data.frame(data.table::rbindlist(sesame_qc))
##  rownames(sesame_qc) <- names(sesame_sset)
##}
##
##save(sesame_seg_normal_normal,
##     file=paste0(work.dir,
##                 file.sep,
##                 "sesame_seg_normal_normal.RData"))
##
##save(sesame_seg_normal_tumor,
##     file=paste0(work.dir,
##                 file.sep,
##                 "sesame_seg_normal_tumor.RData"))
##
##if(sesame.calc.betas==TRUE){
##  sesame_betas <- sesame::openSesame(idat_prefixes,
##                                     mask = TRUE,
##                                     sum.TypeI = TRUE,
##                                     platform = sesame.platform)
##}
##
##if(sesame.infer.sex.karyotypes==TRUE){
##  sesame.karyotype <- foreach(i = 1:length(names(sesame_sset))) %do% {
##    sesame::inferSexKaryotypes(sesame_sset[[i]])
##  }
##  names(sesame.karyotype) <- names(sesame_sset)
##  sesame.karyotype <- as.data.frame(unlist(sesame_karyotype))
##  colnames(sesame.karyotype) <- "karyotype"
##  data.table::fwrite(sesame.karyotype,
##                     file = sesame.output.dir,
##                     col.names = TRUE,
##                     row.names = FALSE,
##                     quote = FALSE)
##}
##
##if(sesame.calc.qc==TRUE){
##  sesame_qc[[i]] <- foreach(i = 1:length(names(sesame_sset))) %do% {
##    sesame::sesameQC(sesame_sset[[i]])
##  }
##  names(sesame_qc) <- names(sesame_sset)
##  sesame_qc <- as.data.frame(data.table::rbindlist(sesame_qc))
##  rownames(sesame_qc) <- names(sesame_sset)
##}
##if(sesame.calc.betas==TRUE){
##  sesame_betas <- sesame::openSesame(treatment_idat_prefixes,
##                                     mask = TRUE,
##                                     sum.TypeI = TRUE,
##                                     platform = sesame.platform)
##}
##
##if(sesame.calc.karyotype==TRUE){
##  sesame_karyotype <- foreach(i = 1:length(names(sesame_sset))) %do% {
##    sesame::inferSexKaryotypes(sesame_sset[[i]])
##  }
##  names(sesame_karyotype) <- names(sesame_sset)
##  sesame_karyotype <- as.data.frame(unlist(sesame_karyotype))
##  colnames(sesame_karyotype) <- "karyotype"
##}
##
##if(sesame.calc.qc==TRUE){
##  sesame_qc <- foreach(i = 1:length(names(sesame_sset))) %do% {
##    sesame::sesameQC(sesame_sset[[i]])
##  }
##  names(sesame_qc) <- names(sesame_sset)
##  sesame_qc <- as.data.frame(data.table::rbindlist(sesame_qc))
##  rownames(sesame_qc) <- names(sesame_sset)
##}
##
##if(sesame.calc.rg==TRUE){
##  sesame_rgset <- sesame::SigSetsToRGChannelSet(sesame_sset)
##}
