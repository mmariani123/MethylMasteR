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
#' @param sesame.hm450.mean.correct
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
                                 sesame.hm450.mean.correct = FALSE,
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
    ExperimentHub::setExperimentHubOption("CACHE", sesame.idat.files.dir)
    ExperimentHub::ExperimentHub()
    treatment_idat_prefixes <- sesame::searchIDATprefixes(sesame.idat.files.dir,
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

    ##colnames(sesame_seg)
    ##"chr"
    ##"startPos"
    ##"endPos"
    ##"median"
    ##"mean"
    ##"sd"
    ##"smp"
    ##"p.val"

    if(sesame.hm450.mean.correct==TRUE){

      ##ReadData {CopyNumber450kCancer}	R Documentation
      ##Function Reads the Data (i.e. regions file and
      ##sample list file) for CopyNumber450kCancer
      ##Description
      ##The input should be two files regions file and
      ##sample list file.
      ##regions file: contains the data for all the
      ##regions/segments in all the sample
      ##sample list file: contains the number of the samples,
      ##the names of the samples and user comment.

      ##The header of the segments/regions file should
      ##be in this order and with these names:
      ##"Sample", "Chromosome", "bp.Start", "bp.End",
      ##"Num.of.Markers", "Mean".
      ##The segments file should have all the samples
      ##in one file Be carful for the dots and it is
      ##case sensitive.

      ##Sample	Sample name
      ##Chromosome	Chromosome number chr1, chr2,
      ##....., chrX, chrY
      ##bp.Start	number, start point for the segment
      ##bp.End	end point for the segment
      ##Num.of.Markers	Number of the probes or
      ##markers in the segment
      ##Mean	is the log value for the segment
      ##(the mean of the log values for all the
      ##probes in the segment, it is the same
      ##value that is used in CopyNumber450k package)
      ##The header of the sample list file should be
      ##in this order and with these names: To check
      ##if the header of the sample list file is ok
      ##"Number", "Sample", "Comment"
      ##Be carful it is case sensitive.

      ##Number	is the number of the sample 1,2,3,....
      ##Sample	the name of the samples
      ##Comment	any comment the user want to see
      ##in the reviewing step and in the QC file,
      ##(ex. karyotyping)
      ##Usage
      ##ReadData(regions_file,
      ##Sample_list,
      ##copynumber450k = FALSE)
      ##Arguments
      ##regions_file
      ##The segmentaion file (CSV file)

      ##Sample_list
      ##The CSV file that contains the names of the samples
      ##and the user comments

      ##copynumber450k
      ##True if the file is the output of copynumber450k,
      ##defualt is FALSE.

      ##Examples
      ##example
      ##the package contains example files:
      ##regions.csv and sample_list.csv
      #to load the example regions.csv
      ##and sample_list.csv files
      ##regions <- system.file("extdata",
      ##"regions.csv",
      ##package="CopyNumber450kCancer")
      ##sample_list <- system.file("extdata",
      ##"sample_list.csv",
      ##package="CopyNumber450kCancer")

      ##Create the object for the package
      ##object <- ReadData(regions,sample_list)

      ##Baseline autocorrection,
      ##this will creat different plot and QC
      ##and new regions file in the working directory
      ##object <- AutoCorrectPeak(object)

      ##For manual revision and manual baseline
      ##determination
      ##object <- ReviewPlot(object)

      ##To plot the final plots
      ##PlotCNV(object) ## to plot all the samples
      ##PlotCNV(object,
      ##select= c(1,4),
      ##comment=FALSE,
      ##cutoff=0.1,
      ##markers=20) ## to plot some samples

for(i in 1:length(names(sesame_seg))){

  signals.now <- sesame_seg[[i]]$seg.signals

  ##"ID"
  ##"chrom"
  ##"loc.start"
  ##"loc.end"
  ##"num.mark"
  ##"seg.mean"
  ##"seg.sd"
  ##"seg.median"
  ##"seg.mad"
  ##"pval"
  ##"lcl"
  ##"ucl"

  regions.df <- data.frame(Sample           = signals.now$ID,
                           Chromosome       = signals.now$chrom,
                           bp.Start         = signals.now$loc.start,
                           bp.End           = signals.now$loc.end,
                           Num.of.Markers   = signals.now$num.mark,
                           Mean             = signals.now$seg.mean,
                           stringsAsFactors = FALSE)

      sample_list.df <-
        data.frame(Number=seq(from=1,
                              by=1,
                              to=length(
                                unique(signals.now$ID))),
                   Sample=unique(signals.now$ID),
                   Comment="",
                   stringsAsFactors = FALSE)

      sample_list <- paste0(output.dir,
                            .Platform$file.sep,
                            "sample_list.csv")
      write.table(sample_list.df,
                  file = sample_list,
                  sep=",",
                  row.names = FALSE,
                  quote = FALSE)

      regions_file <- paste0(output.dir,
                             .Platform$file.sep,
                             "regions.csv")
      write.table(regions.df,
                  file = regions_file,
                  sep=",",
                  row.names = FALSE,
                  quote = FALSE)

      signals.in <-
        CopyNumber450kCancer::AutoCorrectPeak(
          ReadData(regions_file,
                   sample_list,
                   copynumber450k = FALSE))

      ## For manual revision and manual baseline determination
      ##signals.in <- CopyNumber450kCancer::ReviewPlot(signals.in)
      ##Above requires interacting with plot
      ## To plot the final plots

      CopyNumber450kCancer::PlotCNV(signals.in) ## to plot all the samples

      ##PlotCNV(signals.in,
      ##        select= c(1,4),
      ##        comment=FALSE,
      ##        cutoff=0.1,
      ##        markers=20) ## to plot some samples

    }
    }

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

    treatment.platform.1 <- unique(sesame.sample.sheet.df[
      sesame.sample.sheet.df$Sample_Name %in%
        treatment.samples.1, "Platform"])

    treatment.platform.2 <- unique(sesame.sample.sheet.df[
      sesame.sample.sheet.df$Sample_Name %in%
        treatment.samples.2, "Platform"])

    sesame_sset.1 <- sesame::openSesame(treatment_idat_prefixes.1,
                                      mask = TRUE,
                                      sum.TypeI = TRUE,
                                      platform = treatment.platform.1,
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
                                        platform = treatment.platform.2,
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
