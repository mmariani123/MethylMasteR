#!/usr/bin/env Rscript

#' @title ##methyl_master_custom
#' @description MethylMasteR run custom function
#' @param custom.idat.files.dir The input directory for the custom routine
#' @param custom.output.dir The ouptut directory for the custom routine
#' @param custom.sample.sheet.path The path to the MethylMaster sample sheet
#' @param custom.comparison the 2-element vector of Sample_Group levels to be
#' compared. First element is taken as the treatment and second as the control
#' if "reference" is set to "internal" the second element is ignored
#' @param custom.file.sep the file separator to use
#' @param custom.data.cache the custom data cache to use
#' @param custom.data.normal the custom normal data set to use,
#' e.g. Epic.5.Normal
#' @param custom.ref.version the custom reference version (default is hg38)
#' @param custom.reference the custom reference to use
#' @param custom.save.seg save the segmentation results as .RData object
#' @param ... additional parameters to passs to methyl_master_custom
#' @import data.table
#' @import dplyr
#' @import sesame
#' @import sesameData
#' @import ExperimentHub
#' @import future
#' @import profvis
#' @return A seg object for downstream analysis
#' @export
methyl_master_custom <- function(custom.idat.files.dir       = getwd(),
                                 custom.output.dir           = getwd(),
                                 custom.ref                  = NULL,
                                 custom.sample.sheet.path    = NULL,
                                 custom.comparison           = NULL,
                                 custom.file.sep             = NULL,
                                 custom.data.cache           = "EPIC",
                                 custom.data.normal          = "EPIC.5.normal",
                                 custom.ref.version          = "hg38",
                                 custom.reference            = "internal",
                                 custom.save.seg             = FALSE,
                                 ...
                      ){

##Get the samples from the sample sheet:
custom.sample.sheet.df <- read.csv(custom.sample.sheet.path,
                                   header=TRUE,
                                   stringsAsFactors = FALSE)

if(custom.reference=="internal"){

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

    sesame_ssets_normal <- sesameData::sesameDataGet(custom.data.normal)
    ExperimentHub::setExperimentHubOption("CACHE", custom.idat.files.dir)
    ExperimentHub::ExperimentHub()
    ##sesameData::sesameDataCache("EPIC")
    sesameDataCacheAll()

    treatment_idat_prefixes <- sesame::searchIDATprefixes(custom.idat.files.dir,
                                                        recursive=TRUE)
    ##sesameData::sesameDataCache("HM450")

    treatment.names <- custom.sample.sheet.df[
      custom.sample.sheet.df$Sample_Group %in%
        custom.comparison[1],"Sample_Name"]

    treatment_idat_prefixes <-  treatment_idat_prefixes[gsub(".*/(?!/)",
                                     "",
                                     treatment_idat_prefixes,
                                     perl=TRUE) %in%
                                     treatment.names]

    treatment.platform <- unique(custom.sample.sheet.df[
                                 custom.sample.sheet.df$Sample_Name %in%
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
                            refversion = custom.ref.version)
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

    ##sub.dir <- paste0(custom.output.dir,
    ##                  .Platform$file.sep,
    ##                  "auto_corrected")

    ##if(exists(sub.dir)){

    ##  unlink(sub.dir, recursive = TRUE)
    ##  dir.create(sub.dir)

    ##}else{

    ##  dir.create(sub.dir)

    ##}

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

    sesame.now <- sesame_seg
    rm(sesame_seg)

    custom_seg <- list()
    for(i in 1:length(names(sesame.now))){

      signals.now <- sesame.now[[i]]$seg.signals

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

      regions.df <- data.frame(Sample         = signals.now$ID,
                           Chromosome         = signals.now$chrom,
                           bp.Start           = signals.now$loc.start,
                           bp.End             = signals.now$loc.end,
                           Num.of.Markers     = signals.now$num.mark,
                           Mean               = signals.now$seg.mean,
                           stringsAsFactors   = FALSE)

      sample_list.df <-
        data.frame(Number=seq(from=1,
                              by=1,
                              to=length(
                                unique(signals.now$ID))),
                   Sample=unique(signals.now$ID),
                   Comment="",
                   stringsAsFactors = FALSE)

      sample_list <- paste0(custom.output.dir,
                            .Platform$file.sep,
                            "sample_list.csv")
      write.table(sample_list.df,
                  file = sample_list,
                  sep=",",
                  row.names = FALSE,
                  quote = FALSE)

      regions_file <- paste0(custom.output.dir,
                             .Platform$file.sep,
                             "regions.csv")
      write.table(regions.df,
                  file = regions_file,
                  sep=",",
                  row.names = FALSE,
                  quote = FALSE)

      custom_seg[[i]] <-
        AutoCorrectPeak.mm(
          object=CopyNumber450kCancer::ReadData(regions_file,
                                                sample_list,
                                                copynumber450k = FALSE),
          output.dir = paste0(output.dir,
                              .Platform$file.sep,
                              names(sesame.now)[i])
        )


      ## For manual revision and manual baseline determination
      ##signals.in <- CopyNumber450kCancer::ReviewPlot(signals.in)
      ##Above requires interacting with plot
      ## To plot the final plots

      ##CopyNumber450kCancer::PlotCNV(signals.in) ## to plot all the samples

      ##PlotCNV(signals.in,
      ##        select= c(1,4),
      ##        comment=FALSE,
      ##        cutoff=0.1,
      ##        markers=20) ## to plot some samples

    }

    file.remove(regions_file)
    file.remove(sample_list)

  }

names(custom_seg) <- names(sesame_sset)
seg <- list(custom_seg)

##}else if(sesame.reference=="comparison"){
##
##    treatment.samples <-
##      sesame.sample.sheet.df[
##        sesame.sample.sheet.df$Sample_Group==sesame.comparison[1],
##        "Sample_Name"]
##
##    control.samples <-
##        sesame.sample.sheet.df[
##        sesame.sample.sheet.df$Sample_Group==sesame.comparison[2],
##        "Sample_Name"]
##
##      treatment.paths <- sesame.sample.sheet.df[
##      sesame.sample.sheet.df$Sample_Group==sesame.comparison[1],
##      "Basename"]
##
##    control.paths   <-  sesame.sample.sheet.df[
##      sesame.sample.sheet.df$Sample_Group==sesame.comparison[2],
##      "Basename"]
##
##    treatment.platform <-     unique(sesame.sample.sheet.df[
##      sesame.sample.sheet.df$Sample_Group==sesame.comparison[1],
##      "Platform"])
##
##    control.platform   <-     unique(sesame.sample.sheet.df[
##      sesame.sample.sheet.df$Sample_Group==sesame.comparison[2],
##      "Platform"])
##
##    ExperimentHub::setExperimentHubOption("CACHE",
##                         sesame.idat.files.dir)
##
##    ExperimentHub::ExperimentHub()
##
##    idat_prefixes <- sesame::searchIDATprefixes(sesame.idat.files.dir,
##                                              recursive=TRUE)
##
##    idat_prefixes.treatment <-
##      idat_prefixes[gsub(".*/","",idat_prefixes) %in%
##                      gsub(paste0(".*",sesame.file.sep),
##                           "",treatment.paths)]
##
##    idat_prefixes.control   <-
##      idat_prefixes[gsub(".*/","",idat_prefixes) %in%
##                      gsub(paste0(".*",sesame.file.sep),
##                           "",control.paths)]
##
##    sesameData::sesameDataCacheAll()
##
##    sesame_sset.treatment <- sesame::openSesame(idat_prefixes.treatment,
##                            mask = TRUE,
##                            sum.TypeI = TRUE,
##                            platform = treatment.platform,
##                            what="sigset")
##
##    sesame_sset.control <- sesame::openSesame(idat_prefixes.control,
##                            mask = TRUE,
##                            sum.TypeI = TRUE,
##                            platform = control.platform,
##                            what="sigset")
##
##    sesame_seg <- foreach(i = 1:length(names(sesame_sset.treatment))) %do% {
##      sesame::cnSegmentation(sesame_sset.treatment[[i]],
##                             sesame_sset.control,
##                             refversion = sesame.ref.version)
##    }
##    names(sesame_seg) <- names(sesame_sset.treatment)
##
##    seg <- list(sesame_seg)
##
##} ##End reference

} ##End methyl_master_custom
