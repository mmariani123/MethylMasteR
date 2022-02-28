#!/usr/bin/env Rscript

#' @title ##methyl_master_sesame
#' @description MethylMasteR run SeSAMe function
#' @param sesame.idat.files.dir The input directory for the sesame routine
#' @param sesame.output.dir The ouptut directory for the sesame routine
#' @param sesame.sample.sheet.path The path to the MethylMaster sample sheet
#' @param sesame.comparison the 2-element vector of Sample_Group levels to be
#' compared. First element is taken as the treatment and second as the control
#' if "reference" is set to "internal" the second element is ignored
#' @param sesame.file.sep the file separator to use
#' @param sesame.data.cache the sesame data cache to use
#' @param sesame.data.normal the sesame normal data set to use,
#' e.g. Epic.5.Normal
#' @param sesame.genome.version the sesame reference version (default is hg38)
#' @param sesame.reference the sesame reference to use
#' @param sesame.split.by which column, if any, to split the analyses by
#' @param sesame.save.seg save the segmentation results as .RData object
#' @param ... additional parameters to passs to methyl_master_sesame
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
                                 sesame.genome.version              = "hg38",
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

    ##----------------------------------------------------------
    ##| SEnsible Step-wise Analysis of DNA MEthylation (SeSAMe)
    ##| --------------------------------------------------------
    ##| Please cache the annotation data for your array platform
    ##| (e.g. EPIC) by calling "sesameDataCache("EPIC")"
    ##| or "sesameDataCacheAll()". This needs to be done only
    ##| once per SeSAMe installation.

    ##Get sesame normal samples:
    ##Loads all refs as <sesame.data.cache> defaults to "EPIC"?
    ##sesameData::sesameDataCache(sesame.data.cache)

    sesame_ssets_normal <- sesameData::sesameDataGet(sesame.data.normal)
    tmp.dir <- paste0(sesame.output.dir,
                      .Platform$file.sep,
                      "tmp")
    dir.create(path=tmp.dir)
    ExperimentHub::setExperimentHubOption("CACHE",
                                          ##sesame.idat.files.dir,
                                          tmp.dir)
    ExperimentHub::ExperimentHub()
    ##sesameData::sesameDataCache("EPIC")
    sesameData::sesameDataCacheAll()

    treatment_idat_prefixes <- sesame::searchIDATprefixes(sesame.idat.files.dir,
                                                        recursive=TRUE)
    ##sesameData::sesameDataCache("HM450")

    treatment.names <- sesame.sample.sheet.df[
      sesame.sample.sheet.df$Sample_Group %in%
        sesame.comparison[1],"Sample_Name"]

    treatment_idat_prefixes <-  treatment_idat_prefixes[gsub(".*/(?!/)",
                                     "",
                                     treatment_idat_prefixes,
                                     perl=TRUE) %in%
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
                            refversion = sesame.genome.version)
    }
    names(sesame_seg) <- names(sesame_sset)

  seg <- list(sesame_seg)

}else if(sesame.reference=="comparison"){

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

    tmp.dir <- paste0(sesame.output.dir,
                      .Platform$file.sep,
                      "tmp")
    dir.create(path=tmp.dir)
    ExperimentHub::setExperimentHubOption("CACHE",
                                          ##sesame.idat.files.dir,
                                          tmp.dir)

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
                             refversion = sesame.genome.version)
    }
    names(sesame_seg) <- names(sesame_sset.treatment)

    seg <- list(sesame_seg)

} ##End reference

##remove the temp dir for cached files:
unlink(tmp.dir, recursive = TRUE)

return(seg)

} ##End methyl_master_sesame
