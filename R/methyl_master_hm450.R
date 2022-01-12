#!/usr/bin/env Rscript

#' @title methyl_master_hm450
#' @description MethylMaster version of HM450 analysis
#' @param hm450.input.dir The input .idat directory
#' @param hm450.output.dir The output directory for the HM450 routine
#' @param hm450.sample.sheet The MethylMaster sample sheet path
#' @param hm450.reference The reference to use for the HM450 routine
#' @param hm450.workflow The specific HM450 workflow to use: "A", "B", or "C"
#' @param hm450.file.sep The system-specific file separator to use
#' @param hm450.comparison The two-element MethylMaster comparison vector
#' to splut up the HM450 analysis by
#' @param hm450.sesame.data.cache The sesame data control data cache: if using
#' a sesame dataset as the reference
#' @param hm450.sesame.data.normal The sesame data normal data set: if using a
#' sesame dataset as the reference
#' @param hm450.genome.version The sesame ref version (e.g. hg38) if using
#' a sesame dataset as the reference
#' @param hm450.save.seg Whether to save the HM450 routine segmentatin results
#' @param ... Additional parameters to pass to methyl_master_hm450
#' @import data.table
#' @import Biobase
#' @importFrom ExperimentHub setExperimentHubOption
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom sesameData sesameDataGet
#' @importFrom sesameData sesameDataCache
#' @importFrom sesame SigSetsToRGChannelSet
#' @importFrom sesame openSesame
#' @return HM450 routine segmentation CNV results stored in a list
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
                                hm450.genome.version="hg38",
                                hm450.save.seg=FALSE,
                                ...
                                ){

  ##load(hm450.anno.file.path)
  ##load("G:\\My Drive\\dartmouth\\salas_lab_working\\cnv\\
  ##cnv_testing\\probe450kfemanno.rda")
  ##load("G:\\My Drive\\dartmouth\\salas_lab_working\\cnv\\
  ##cnv_testing\\hm450.manifest.hg38.rda")
  ##https://www.bioconductor.org/packages/release/BiocViews.html#___IlluminaChip

  ##candidates_data_treatment_B.1 <-
  ##readRDS(file="C:\\Users\\Mike\\Desktop\\candidates_data_treatment_B.1.RDS")

  ##hm450.anno.file.path <- paste0("G:\\My Drive\\dartmouth",
  ##                                 "\\salas_lab\\cnv",
  ##                                 "\\hm450.manifest.hg38.rda")
  ##load(hm450.anno.file.path)
  ##load("G:\\My Drive\\dartmouth\\salas_lab_working\\cnv
  ##\\cnv_testing\\probe450kfemanno.rda")
  ##load("G:\\My Drive\\dartmouth\\salas_lab_working\\cnv
  ##\\cnv_testing\\hm450.manifest.hg38.rda")
  ##https://www.bioconductor.org/packages/release/BiocViews.html#___IlluminaChip
  ##annotation1 <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  ##annotation1 <- as.data.frame(annotation1)

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

      candidates_data_treatment_B <-
        findSegments2(proc_treatment_B[, , drop = FALSE],
                      med_control_B,
                      proc_control_B) ##Scaled with B but not A

      candidates_data_treatment_B <-
        reformat_hm450(candidates_data_treatment_B, annotation1)

      ## candidates_data_treatment_B.2$start <-
      ##  annotation1[candidates_data_treatment_B.2$startCG %in%
      ##                annotation$Name, annotation1$startPos]

      ##candidates_data_treatment_B.2$end <-
      ##  annotation1[candidates_data_treatment_B.2$endCG %in%
      ##                       annotation$Name, annotation1$endPos]

      seg <- list(candidates_data_treatment_B)

  }else{

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


  } ##End hm450.reference

  return(seg)

} ##End function
