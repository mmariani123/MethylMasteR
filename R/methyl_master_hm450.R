#!/usr/bin/env Rscript

#' @title methyl_master_hm450
#' @description My version of hm450 analyses
#' ############ Samples and ref also need to be in same order ##########
## "CpG probe IDs not in the same order in data and ctrl!"
#' @param hm450.input.dir
#' @param hm450.output.dir
#' @param hm450.sample.sheet
#' @param hm450.reference
#' @param hm450.workflow
#' @param hm450.file.sep
#' @param hm450.comparison
#' @param hm450.split.by
#' @param hm450.sesame.data.cache
#' @param hm450.sesame.data.normal
#' @param hm450.sesame.ref.version
#' @param hm450.save.seg
#' @param hm450.run.mean.correct
#' @param ...
#' @import data.table
#' @import Biobase
#' @importFrom ExperimentHub setExperimentHubOption unbox
#' @importFrom ExperimentHub ExperimentHub unbox
#' @importFrom sesameData sesamesesameDataGet
#' @importFrom sesameData sesamesesameDataCache
#' @importFrom sesame SigSetsToRGChannel
#' @importFrom sesame openSesame
#' @return #seg.out
#' @export
methyl_master_hm450 <- function(hm450.input.dir=NULL,
                                hm450.output.dir=NULL,
                                hm450.sample.sheet.path=NULL,
                                hm450.comparison=c("tumor","normal"),
                                hm450.split.by=NULL,
                                hm450.reference="internal",
                                hm450.workflow="B",
                                hm450.file.sep="/",
                                hm450.sesame.data.cache="EPIC",
                                hm450.sesame.data.normal="EPIC.5.normal",
                                hm450.sesame.ref.version="hg38",
                                hm450.save.seg=FALSE,
                                hm450.run.mean.correct=FALSE,
                                ...
                                ){

  ##load(hm450.anno.file.path)
  ##load("G:\\My Drive\\dartmouth\\salas_lab_working\\cnv\\cnv_testing\\probe450kfemanno.rda")
  ##load("G:\\My Drive\\dartmouth\\salas_lab_working\\cnv\\cnv_testing\\hm450.manifest.hg38.rda")
  ##https://www.bioconductor.org/packages/release/BiocViews.html#___IlluminaChip

  ##candidates_data_treatment_B.1 <-
  ##readRDS(file="C:\\Users\\Mike\\Desktop\\candidates_data_treatment_B.1.RDS")
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  annotation1 <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  annotation1 <- as.data.frame(annotation1)

  reformat_hm450 <- function(df,anno){
    setDT(anno)
    anno[strand=="+",startPos:=pos,]
    anno[strand=="+",endPos:=pos+50,]
    anno[strand=="-",startPos:=pos,]
    anno[strand=="-",endPos:=pos-50,]
    anno.start <- anno[,c("Name","startPos")]
    anno.end   <- anno[,c("Name","endPos")]
    ##which(!candidates_data_treatment_B.1$startCG %in% annotation1$Name)
    ##table(c(1,2,2,3) %in% c(3,2,1))
    ##table(c(1,2,3) %in% c(3,2,11))
    ##table(annotation1$Name %in% candidates_data_treatment_B.1$startCG)
    ##table(candidates_data_treatment_B.1$startCG %in% annotation1$Name)
    ##table(annotation1$Name %in% candidates_data_treatment_B.1$endCG)
    ##table(candidates_data_treatment_B.1$endCG %in% annotation1$Name)
    anno.start$startPos <- as.character(anno.start$startPos)
    anno.end$endPos     <- as.character(anno.end$endPos)
    df <- dplyr::inner_join(df,
                            anno.start,
                            by=c("startCG"="Name"))
    df <- dplyr::inner_join(df,
                            anno.end,
                            by=c("endCG"="Name"))
    df <- df[,c(1,9,10,4,5,6,7,8)]
    return(df)
  }

  ##df.rf <- reformat_hm450(candidates_data_treatment_B.1, annotation1)

  sample.sheet.df <- read.csv(hm450.sample.sheet.path,
                              header=TRUE,
                              stringsAsFactors = FALSE)

  if(hm450.reference=="internal"){

    if(is.null(hm450.split.by)){

    sample.sheet.df <- sample.sheet.df[sample.sheet.df$Sample_Group %in%
                                                  hm450.comparison,]

    sesameData::sesameDataCache(hm450.sesame.data.cache)

    sesame_ssets_normal <- sesameData::sesameDataGet(hm450.sesame.data.normal)

    ExperimentHub::setExperimentHubOption("CACHE", hm450.input.dir)

    ExperimentHub::ExperimentHub()

    idat_prefixes <- sesame::searchIDATprefixes(hm450.input.dir,
                                                recursive=TRUE)
    sesameData::sesameDataCacheAll()

    treatment.names <- sample.sheet.df[
      sample.sheet.df$Sample_Group %in%
        hm450.comparison[1],"Sample_Name"]

    treatment.paths <-
      sample.sheet.df[
        sample.sheet.df$Sample_Group==hm450.comparison[1],
        "Basename"]

    treatment.platform <- unique(sample.sheet.df[
      sample.sheet.df$Sample_Name %in%
        treatment.names, "Platform"])

    idat_prefixes.treatment <-
      idat_prefixes[gsub(".*/","",idat_prefixes) %in%
                      gsub(paste0(".*",hm450.file.sep),
                           "",treatment.paths)]

    sesame_sset <- sesame::openSesame(idat_prefixes.treatment,
                                      mask = TRUE,
                                      sum.TypeI = TRUE,
                                      platform = treatment.platform,
                                      what="sigset")

    hm450_rgset <- sesame::SigSetsToRGChannelSet(sesame_sset)
    hm450_rgset.control <- sesame::SigSetsToRGChannelSet(sesame_ssets_normal)

    hm450_cn_methylset.treatment <-
      minfi::getCN(minfi::preprocessRaw(hm450_rgset))
    names(hm450_cn_methylset.treatment) <- names(sesame_sset)

    hm450_cn_methylset.control <-
      minfi::getCN(minfi::preprocessRaw(hm450_rgset.control))
    names(hm450_cn_methylset.control) <- names(sesame_ssets_normal)

    if(hm450.workflow=="A"){

      proc_treatment_A  <- hm450_cn_methylset.treatment
      rm(hm450_cn_methylset.treatment)

      proc_control_A <- hm450_cn_methylset.control
      rm(hm450_cn_methylset.control)
      med_control_A  <- apply(proc_control_A, 1, "median")

      candidates_data_treatment_A <-
        findSegments2(proc_treatment_A[, , drop = FALSE],
                      med_control_A,
                      proc_control_A) ##Scaled with B but not A

      candidates_data_treatment_A <-
        reformat_hm450(candidates_data_treatment_A, annotation1)

      seg <- list(candidates_data_treatment_A)

    }else if(hm450.workflow=="B"){

      ## With z-Transformation, illumina
      proc_control_B  <- hm450_cn_methylset.control
      rm(hm450_cn_methylset.control)
      proc_control_B[is.infinite(proc_control_B)] <- NA
      proc_control_B <- scale(proc_control_B) ##Z transform
      med_control_B <- apply(proc_control_B, 1, "median")

      proc_treatment_B  <- hm450_cn_methylset.treatment
      rm(hm450_cn_methylset.treatment)
      proc_treatment_B[is.infinite(proc_treatment_B)] <- NA
      proc_treatment_B <- scale(proc_treatment_B)
      med_treatment_B <- apply(proc_treatment_B, 1, "median")

      candidates_data_treatment_B <-
        findSegments2(proc_treatment_B[, , drop = FALSE],
                      med_control_B,
                      proc_control_B) ##Scaled with B but not A

      candidates_data_treatment_B <-
        reformat_hm450(candidates_data_treatment_B, annotation1)

      ##colnames(candidates_data_treatment_B)
      ##"chr"
      ##"startPos"
      ##"endPos"
      ##"median"
      ##"mean"
      ##"sd"
      ##"smp"
      ##"p.val"

      if(hm450.run.mean.correct==TRUE){

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

        candidates_data_treatment_B <-
          CopyNumber450kCancer::AutoCorrectPeak(candidates_data_treatment_B)

      }

      ## candidates_data_treatment_B.2$start <-
      ##  annotation1[candidates_data_treatment_B.2$startCG %in%
      ##                annotation$Name, annotation1$startPos]

      ##candidates_data_treatment_B.2$end <-
      ##  annotation1[candidates_data_treatment_B.2$endCG %in%
      ##                       annotation$Name, annotation1$endPos]

      seg <- list(candidates_data_treatment_B)

    }else if(hm450.workflow=="C"){

      ## Conumee-path, Illumina

      proc_control_C <- hm450_cn_methylset.control
      rm(hm450_cn_methylset.control)

      proc_treatment_C <- hm450_cn_methylset.treatment
      rm(hm450_cn_methylset.treatment)

      ##Beware of getting the following error:
      ##Beware of error:
      ##Error in cnAnalysis450k::runConumee(proc_treatment_C,
      ##                                    proc_control_C) :
      ##  CpG probe IDs not in the same order in data and ctrl!

      ##proc_treatment_sorted_C <-
      ##  proc_treatment_sorted_C[female_shared_names,]
      ##proc_control_sorted_C   <-
      ##  proc_control_sorted_C[male_shared_names,]

      ####################### RUN CONUMEE ##############################

      proc_treatment_C_ordered <-
        proc_control_C[order(rownames(proc_treatment_C)),]

      proc_control_C_ordered <-
        proc_control_C[rownames(proc_control_C) %in%
                                   rownames(proc_treatment_C_ordered),]
      proc_control_C_ordered <-
        proc_control_C[rownames(proc_treatment_C_ordered),]

      require(IlluminaHumanMethylation450kanno.ilmn12.hg19)
      require(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

      candidates_data_treatment_C <-
        cnAnalysis450k::runConumee(proc_treatment_C_ordered,
                                   proc_control_C_ordered,
                                   ##proc_treatment_sorted_C.2,
                                   ##proc_control_sorted_C.2,
                                   what="segments"
        )

      candidates_data_treatment_C <-
        reformat_hm450(candidates_data_treatment_C, annotation1)

      seg <- list(candidates_data_treatment_C.1)

    } ##End hm450 workflow

    }else{

      stop(paste0("Error: Epic.5.Normal internal ",
                  "\nreference samples are all male",
                  "\ncan't split by sex!"))

      sample.sheet.df <- sample.sheet.df[sample.sheet.df$Sample_Group %in%
                                                  hm450.comparison,]

      split.by.cat <- unique(sesame.sample.sheet.df[[sesame.split.by]])

      sesameData::sesameDataCache(hm450.sesame.data.cache)
      sesame_ssets_normal <- sesameData::sesameDataGet(hm450.sesame.data.normal)
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
        sample.sheet.df[
          sample.sheet.df$Sample_Group==hm450.comparison[1] &
            sample.sheet.df[[hm450.split.by]]==split.by.cat[1],
            "Sample_Name"]

      treatment.samples.2 <-
        sample.sheet.df[
          sample.sheet.df$Sample_Group==sesame.comparison[1] &
            sample.sheet.df[[hm450.split.by]]==split.by.cat[2],
            "Sample_Name"]

      idat_prefixes.treatment.1 <-
        idat_prefixes[gsub(".*/","",idat_prefixes) %in%
                      gsub(paste0(".*",hm450.file.sep),
                             "",treatment.paths.1)]

      idat_prefixes.treatment.2 <-
        idat_prefixes[gsub(".*/","",idat_prefixes) %in%
                      gsub(paste0(".*",hm450.file.sep),
                             "",treatment.paths.2)]

      treatment.platform.1 <- unique(sample.sheet.df[
        sample.sheet.df$Sample_Name %in%
          treatment.samples.1, "Platform"])

      treatment.platform.2 <- unique(sample.sheet.df[
        sample.sheet.df$Sample_Name %in%
          treatment.samples.2, "Platform"])

      sesame_sset.1 <- sesame::openSesame(treatment_idat_prefixes.1,
                                          mask = TRUE,
                                          sum.TypeI = TRUE,
                                          platform = treatment.platform.1,
                                          what="sigset")

      hm450_rgset.1 <- sesame::SigSetsToRGChannelSet(sesame_sset.1)

      hm450_cn_methylset_treatment.1 <-
        minfi::getCN(minfi::preprocessRaw(hm450_rgset.1))

      names(hm450_cn_methylset_treatment.1) <- names(sesame_sset.1)

      sesame_sset.2 <- sesame::openSesame(treatment_idat_prefixes.2,
                                          mask = TRUE,
                                          sum.TypeI = TRUE,
                                          platform = treatment.platform.2,
                                          what="sigset")

      hm450_rgset.2 <- sesame::SigSetsToRGChannelSet(sesame_sset.2)

      hm450_cn_methylset_treatment.2 <-
        minfi::getCN(minfi::preprocessRaw(hm450_rgset.2))

      names(hm450_cn_methylset_treatment.2) <- names(sesame_sset.2)

      seg <- list(hm450_cn_methylset_treatment.1,
                  hm450_cn_methylset_treatment.2)

    }

  }else{

    if(is.null(hm450.split.by)){

    sample.sheet.df <- sample.sheet.df[sample.sheet.df$Sample_Group %in%
                                       hm450.comparison,]

    ExperimentHub::setExperimentHubOption("CACHE", idat.files.dir)

    ExperimentHub::ExperimentHub()

    idat_prefixes <- searchIDATprefixes(idat.files.dir, recursive=TRUE)

    sesameData::sesameDataCacheAll()

    treatment.names <- sample.sheet.df[
      sample.sheet.df$Sample_Group %in%
        hm450.comparison[1],"Sample_Name"]

    treatment_idat_prefixes <-  idat_prefixes[gsub(".*",
                                "",
                                idat_prefixes) %in%
                                treatment.names]

    treatment.platform <- unique(sample.sheet.df[
      sample.sheet.df$Sample_Name %in%
        treatment.names, "Platform"])

    control.names <- sample.sheet.df[
      sample.sheet.df$Sample_Group %in%
        hm450.comparison[2],"Sample_Name"]

    control_idat_prefixes <-  idat_prefixes[gsub(".*/",
                                            "",
                                            idat_prefixes) %in%
                                            control.names]

    control.platform <- unique(sample.sheet.df[
      sample.sheet.df$Sample_Name %in%
        control.names, "Platform"])

    sesame_ssets_control <- openSesame(control_idat_prefixes,
                                   mask = TRUE,
                                   sum.TypeI = TRUE,
                                   platform = control.platform,
                                   what="sigset")

    sesame_ssets_treatment <- openSesame(treatment_idat_prefixes,
                                   mask = TRUE,
                                   sum.TypeI = TRUE,
                                   platform = treatment.platform,
                                   what="sigset")

    sesame_rgset_treatment <-
      sesame::SigSetsToRGChannelSet(sesame_ssets_treatment)

    sesame_rgset_control <-
      sesame::SigSetsToRGChannelSet(sesame_ssets_control)

    sesame_cn_methylset_treatment <-
      minfi::getCN(minfi::preprocessRaw(sesame_rgset_treatment))

    sesame_cn_methylset_control <-
      minfi::getCN(minfi::preprocessRaw(sesame_rgset_control))

    if(hm450.workflow=="A"){

      proc_treatment_A  <- sesame_cn_methylset_treatment
      rm(sesame_cn_methylset_treatment)

      proc_control_A <- sesame_cn_methylset_control
      rm(sesame_cn_methylset_control)
      med_control_A  <- apply(proc_control_A, 1, "median")

      candidates_data_treatment_A <-
        findSegments2(proc_treatment_A[, , drop = FALSE],
                      med_control_A,
                      proc_control_A) ##Scaled with B but not A

      seg <- list(candidates_data_treatment_A)

    }else if(hm450.workflow=="B"){

      ## With z-Transformation, illumina
      proc_treatment_B  <- sesame_cn_methylset_treatment
      rm(sesame_cn_methylset_treatment)
      proc_treatment_B[is.infinite(proc_treatment_B)] <- NA
      proc_treatment_B <- scale(proc_treatment_B)

      proc_control_B  <- sesame_cn_methylset_control
      rm(sesame_cn_methylset_control)
      proc_control_B[is.infinite(proc_control_B)] <- NA
      proc_control_B <- scale(proc_control_B) ##Z transform
      med_control_B <- apply(proc_control_B, 1, "median")

      candidates_data_treatment_B <-
        findSegments2(proc_treatment_B[, , drop = FALSE],
                      med_control_B,
                      proc_control_B) ##Scaled with B but not A

      seg <- list(candidates_data_treatment_B)

    }else if(hm450.workflow=="C"){

      ## Conumee-path, Illumina

      proc_treatment_C <- sesame_cn_methylset_treatment
      rm(sesame_cn_methylset_treatment)

      proc_control_C <- sesame_cn_methylset_control
      rm(sesame_cn_methylset_control)

      ##What was I doing with shared names here?
      ##proc_treatment_sorted_C.1 <-
      ##  proc_treatment_sorted_C.1[female_shared_names,]
      ##proc_control_sorted_C.1   <-
      ##  proc_control_sorted_C.1[male_shared_names,]

      ####################### RUN CONUMEE ##############################

      candidates_data_treatment_C <-
        cnAnalysis450k::runConumee(proc_treatment_C,
                                   proc_control_C
                                   ##proc_treatment_sorted_C,
                                   ##proc_control_sorted_C
                                   )

      seg <- list(candidates_data_treatment_C)

    }

    }else{

      sample.sheet.df <- sample.sheet.df[sample.sheet.df$Sample_Group %in%
                                                  hm450.comparison,]

      split.by.cat <- unique(sample.sheet.df[[hm450.split.by]])

      treatment.samples.1 <-
        sample.sheet.df[
          sample.sheet.df$Sample_Group==hm450.comparison[1] &
            sample.sheet.df[[hm450.split.by]]==split.by.cat[1],
            "Sample_Name"]

      control.samples.1 <-
        sample.sheet.df[
          sample.sheet.df$Sample_Group==hm450.comparison[2] &
            sample.sheet.df[[hm450.split.by]]==split.by.cat[1],
            "Sample_Name"]

      treatment.paths.1 <-
        sample.sheet.df[
          sample.sheet.df$Sample_Group==hm450.comparison[1] &
            sample.sheet.df[[hm450.split.by]]==split.by.cat[1],
            "Basename"]

      control.paths.1   <-
        sample.sheet.df[
          sample.sheet.df$Sample_Group==hm450.comparison[2] &
            sample.sheet.df[[hm450.split.by]]==split.by.cat[1],
            "Basename"]

      treatment.platform.1 <-
        unique(sample.sheet.df[
          sample.sheet.df$Sample_Group==hm450.comparison[1] &
            sample.sheet.df[[hm450.split.by]]==split.by.cat[2],
            "Platform"])

      control.platform.1 <-
        unique(sample.sheet.df[
          sample.sheet.df$Sample_Group==hm450.comparison[2] &
            sample.sheet.df[[hm450.split.by]]==split.by.cat[2],
            "Platform"])

      treatment.samples.2 <-
        sample.sheet.df[
          sample.sheet.df$Sample_Group==hm450.comparison[1] &
            sample.sheet.df[[hm450.split.by]]==split.by.cat[2],
            "Sample_Name"]

      control.samples.2 <-
        sample.sheet.df[
          sample.sheet.df$Sample_Group==hm450.comparison[2] &
            sample.sheet.df[[hm450.split.by]]==split.by.cat[2],
            "Sample_Name"]

      treatment.paths.2 <-
        sample.sheet.df[
          sample.sheet.df$Sample_Group==hm450.comparison[1] &
            sample.sheet.df[[hm450.split.by]]==split.by.cat[2],
            "Basename"]

      control.paths.2   <-
        sample.sheet.df[
          sample.sheet.df$Sample_Group==hm450.comparison[2] &
            sample.sheet.df[[hm450.split.by]]==split.by.cat[2],
            "Basename"]

      treatment.platform.2 <-
        unique(sample.sheet.df[
          sample.sheet.df$Sample_Group==hm450.comparison[1] &
            sample.sheet.df[[hm450.split.by]]==split.by.cat[2],
            "Platform"])

      control.platform.2   <-
        unique(sample.sheet.df[
          sample.sheet.df$Sample_Group==hm450.comparison[2] &
            sample.sheet.df[[hm450.split.by]]==split.by.cat[2],
            "Platform"])

      ExperimentHub::setExperimentHubOption("CACHE",
                                          hm450.input.dir)

      ExperimentHub::ExperimentHub()

      idat_prefixes <- sesame::searchIDATprefixes(hm450.input.dir,
                                                recursive=TRUE)

      idat_prefixes.treatment.1 <-
        idat_prefixes[gsub(".*/","",idat_prefixes) %in%
                      gsub(paste0(".*",hm450.file.sep),
                           "",treatment.paths.1)]

      idat_prefixes.treatment.2 <-
        idat_prefixes[gsub(".*/","",idat_prefixes) %in%
                      gsub(paste0(".*",hm450.file.sep),
                             "",treatment.paths.2)]

      idat_prefixes.control.1   <-
        idat_prefixes[gsub(".*/","",idat_prefixes) %in%
                      gsub(paste0(".*",hm450.file.sep),
                           "",control.paths.1)]

      idat_prefixes.control.2   <-
        idat_prefixes[gsub(".*/","",idat_prefixes) %in%
                        gsub(paste0(".*",hm450.file.sep),
                             "",control.paths.2)]

      sesameData::sesameDataCacheAll()

      sesame_sset.treatment.1 <- sesame::openSesame(idat_prefixes.treatment.1,
                                                mask = TRUE,
                                                sum.TypeI = TRUE,
                                                platform = treatment.platform.1,
                                                what="sigset")

      sesame_sset.treatment.2 <- sesame::openSesame(idat_prefixes.treatment.2,
                                                  mask = TRUE,
                                                  sum.TypeI = TRUE,
                                                platform = treatment.platform.2,
                                                  what="sigset")

      sesame_sset.control.1 <- sesame::openSesame(idat_prefixes.control.1,
                                              mask = TRUE,
                                              sum.TypeI = TRUE,
                                              platform = control.platform.1,
                                              what="sigset")

      sesame_sset.control.2 <- sesame::openSesame(idat_prefixes.control.2,
                                                mask = TRUE,
                                                sum.TypeI = TRUE,
                                                platform = control.platform.2,
                                                what="sigset")

      sesame_rgset.treatment.1 <-
        sesame::SigSetsToRGChannelSet(sesame_sset.treatment.1)

      sesame_rgset.treatment.2 <-
        sesame::SigSetsToRGChannelSet(sesame_sset.treatment.2)

      sesame_rgset.control.1 <-
        sesame::SigSetsToRGChannelSet(sesame_sset.control.1)

      sesame_rgset.control.2 <-
        sesame::SigSetsToRGChannelSet(sesame_sset.control.2)

      sesame_cn_methylset.treatment.1 <-
        minfi::getCN(minfi::preprocessRaw(sesame_rgset.treatment.1))

      sesame_cn_methylset.treatment.2 <-
        minfi::getCN(minfi::preprocessRaw(sesame_rgset.treatment.2))

      sesame_cn_methylset.control.1 <-
        minfi::getCN(minfi::preprocessRaw(sesame_rgset.control.1))

      sesame_cn_methylset.control.2 <-
        minfi::getCN(minfi::preprocessRaw(sesame_rgset.control.2))

      if(hm450.workflow=="A"){

        proc_treatment_A.1  <- sesame_cn_methylset.treatment.1
        rm(sesame_cn_methylset.treatment.1)

        proc_treatment_A.2  <- sesame_cn_methylset.treatment.2
        rm(sesame_cn_methylset.treatment.2)

        proc_control_A.1 <- sesame_cn_methylset.control.1
        rm(sesame_cn_methylset.control.1)
        med_control_A.1  <- apply(proc_control_A.1, 1, "median")

        proc_control_A.2 <- sesame_cn_methylset.control.2
        rm(sesame_cn_methylset.control.2)
        med_control_A.2 <- apply(proc_control_A.2, 1, "median")

        candidates_data_treatment_A.1 <-
          findSegments2(proc_treatment_A.1[, , drop = FALSE],
                      med_control_A.1,
                      proc_control_A.1) ##Scaled with B but not A

        candidates_data_treatment_A.2 <-
          findSegments2(proc_treatment_A.2[, , drop = FALSE],
                      med_control_A.2,
                      proc_control_A.2)

        seg <- list(candidates_data_treatment_A.1,
                    candidates_data_treatment_A.2)

      }else if(hm450.workflow=="B"){

        ## With z-Transformation, illumina
        proc_control_B.1  <- sesame_cn_methylset.control.1
        rm(sesame_cn_methylset.control.1)
        proc_control_B.1[is.infinite(proc_control_B.1)] <- NA
        proc_control_B.1 <- scale(proc_control_B.1) ##Z transform
        med_control_B.1 <- apply(proc_control_B.1, 1, "median")

        proc_control_B.2 <- sesame_cn_methylset.control.2
        rm(sesame_cn_methylset.control.2)
        proc_control_B.2[is.infinite(proc_control_B.2)] <- NA
        proc_control_B.2 <- scale(proc_control_B.2)
        med_control_B.2 <- apply(proc_control_B.2, 1, "median")

        proc_treatment_B.1  <- sesame_cn_methylset.treatment.1
        rm(sesame_cn_methylset.treatment.1)
        proc_treatment_B.1[is.infinite(proc_treatment_B.1)] <- NA
        proc_treatment_B.1 <- scale(proc_treatment_B.1)
        med_treatment_B.1 <- apply(proc_treatment_B.1, 1, "median")

        proc_treatment_B.2  <- sesame_cn_methylset.treatment.2
        rm(sesame_cn_methylset.treatment.2)
        proc_treatment_B.2[is.infinite(proc_treatment_B.2)] <- NA
        proc_treatment_B.2 <- scale(proc_treatment_B.2)
        med_treatment_B.2 <- apply(proc_treatment_B.2, 1, "median")

        candidates_data_treatment_B.1 <-
          findSegments2(proc_treatment_B.1[, , drop = FALSE],
                      med_control_B.1,
                      proc_control_B.1) ##Scaled with B but not A

        candidates_data_treatment_B.1 <-
          reformat_hm450(candidates_data_treatment_B.1, annotation1)

        candidates_data_treatment_B.2 <-
          findSegments2(proc_treatment_B.2[, , drop = FALSE],
                      med_control_B.2,
                      proc_control_B.2)

        candidates_data_treatment_B.2 <-
          reformat_hm450(candidates_data_treatment_B.2, annotation1)

       ## candidates_data_treatment_B.2$start <-
        ##  annotation1[candidates_data_treatment_B.2$startCG %in%
        ##                annotation$Name, annotation1$startPos]

        ##candidates_data_treatment_B.2$end <-
        ##  annotation1[candidates_data_treatment_B.2$endCG %in%
        ##                       annotation$Name, annotation1$endPos]

        seg <- list(candidates_data_treatment_B.1,
                    candidates_data_treatment_B.2)

      }else if(hm450.workflow=="C"){

      ## Conumee-path, Illumina

        proc_control_C.1 <- sesame_cn_methylset.control.1
        rm(sesame_cn_methylset.control.1)

        proc_control_C.2 <- sesame_cn_methylset.control.2
        rm(sesame_cn_methylset.control.2)

        proc_treatment_C.1 <- sesame_cn_methylset.treatment.1
        rm(sesame_cn_methylset.treatment.1)

        proc_treatment_C.2 <- sesame_cn_methylset.treatment.2
        rm(sesame_cn_methylset.treatment.2)

        ##What was I doing with shared names here:
        ##proc_treatment_sorted_C.1 <-
        ##  proc_treatment_sorted_C.1[female_shared_names,]
        ##proc_treatment_sorted_C.2 <-
        ##  proc_treatment_sorted_C.2[female_shared_names,]
        ##proc_control_sorted_C.1   <-
        ##  proc_control_sorted_C.1[male_shared_names,]
        ##proc_control_sorted_C.2   <-
        ##  proc_control_sorted_C.2[male_shared_names,]

      ####################### RUN CONUMEE ##############################

        candidates_data_treatment_C.1 <-
          cnAnalysis450k::runConumee(proc_treatment_C.1,
                                     proc_control_C.1
                                     ##proc_treatment_sorted_C.1,
                                     ##proc_control_sorted_C.1
                                     )

        candidates_data_treatment_C.1 <-
          cnAnalysis450k::runConumee(proc_treatment_C.2,
                                     proc_control_C.2
                                     ##proc_treatment_sorted_C.2,
                                     ##proc_control_sorted_C.2
                                     )

        seg <- list(candidates_data_treatment_C.1,
                    candidates_data_treatment_C.2)

      } ##End hm450 workflow

    } ##End hm450.split.by

  } ##End hm450.reference

  return(seg)

} ##End function
