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

    treatment.names <-    sesame.sample.sheet.df[
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
    treatment_idat_prefixes <- sesame::searchIDATprefixes(idat.files.dir,
                                                          recursive=TRUE)
    sesameData::sesameDataCacheAll()

    sesame_sset.1 <- sesame::openSesame(treatment_idat_prefixes,
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

    sesame_sset.2 <- sesame::openSesame(treatment_idat_prefixes,
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
        sesame.sample.sheet.df$Sample_Group==sesame.comparison[1],"Sample_Name"]

    control.samples <-
        sesame.sample.sheet.df[
        sesame.sample.sheet.df$Sample_Group==sesame.comparison[2],"Sample_Name"]

      treatment.paths <-     sesame.sample.sheet.df[
      sesame.sample.sheet.df$Sample_Group==sesame.comparison[1], "Basename"]

    control.paths   <-     sesame.sample.sheet.df[
      sesame.sample.sheet.df$Sample_Group==sesame.comparison[2], "Basename"]

    treatment.platform <-     unique(sesame.sample.sheet.df[
      sesame.sample.sheet.df$Sample_Group==sesame.comparison[1], "Platform"])

    control.platform   <-     unique(sesame.sample.sheet.df[
      sesame.sample.sheet.df$Sample_Group==sesame.comparison[2], "Platform"])

    ExperimentHub::setExperimentHubOption("CACHE",
                         sesame.idat.files.dir)

    ExperimentHub::ExperimentHub()

    idat_prefixes <- sesame::searchIDATprefixes(sesame.idat.files.dir,
                                              recursive=TRUE)

    idat_prefixes.treatment <-
      idat_prefixes[gsub(".*/.*_","",idat_prefixes) %in%
                      gsub(paste0(".*",sesame.file.sep,".*_"),
                           "",treatment.paths)]

    idat_prefixes.control   <-
      idat_prefixes[gsub(".*/.*_","",idat_prefixes) %in%
                      gsub(paste0(".*",sesame.file.sep,".*_"),
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
      idat_prefixes[gsub(".*/.*_","",idat_prefixes) %in%
                    gsub(paste0(".*",sesame.file.sep,".*_"),
                         "",treatment.paths.1)]

    idat_prefixes.control.1   <-
      idat_prefixes[gsub(".*/.*_","",idat_prefixes) %in%
                    gsub(paste0(".*",sesame.file.sep,".*_"),
                         "",control.paths.1)]

    idat_prefixes.treatment.2 <-
      idat_prefixes[gsub(".*/.*_","",idat_prefixes) %in%
                    gsub(paste0(".*",sesame.file.sep,".*_"),
                         "",treatment.paths.2)]

    idat_prefixes.control.2   <-
      idat_prefixes[gsub(".*/.*_","",idat_prefixes) %in%
                    gsub(paste0(".*",sesame.file.sep,".*_"),
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

  } ##End split.by

} ##End reference

if(sesame.save.seg==TRUE){
  save(seg, file =  paste0(sesame.output.dir, .Platform$file.sep, "seg.RData"))
}

######################## FORMATTING #########################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################

if(reference=="internal"){

  if(is.null(split.by)){

    sesame.seg <- seg[[1]]

    sesame_seg_treatment.df <- binding_frames_mm(
      sesame.seg)

    sesame_seg_treatment.df <- sesame_seg_treatement.df[
      sesame_seg_treatment.df$pval <= 0.05,]

    colnames(sesame_seg_treatment.df)[1] <- "Sample_ID"
    seg <- sesame_seg_treatment.df [,c(1,2,3,4,5,6,10)]
    seg$state <- round(2^seg$seg.mean * 2)
    seg$state[seg$state > 4] <- 4
    seg$method <- "sesame"
    row.names(seg) <- NULL
    seg <- na.omit(seg)
    seg$chrom <-
      unlist(strsplit(seg$chrom,split="chr"))[c(FALSE,TRUE)]

    seg.out <- list(seg.1)
    names(seg.out) <- "internal"

    ##For sesame we want to omit the "chr"
    ##from chrom for intial plotting
    ##then add back in downstream for ggplot

  }else{

    sesame.seg.treatment <- sesame.seg[names(sesame.seg) %in%
    sample.sheet.csv[sample.sheet.csv$Sample_Group %in% comparison[1],
    "Sample_Name"]]

    split.by.1 <- unique(sample.sheet.csv[[split.by]])[1]
    split.by.2 <- unique(sample.sheet.csv[[split.by]])[2]

    sesame_seg_treatment.1 <- sesame.seg.treatment[
      names(sesame.seg.treatment) %in%
        sample.sheet.csv[sample.sheet.csv$gender_reported %in% split.by.1,
                         "Sample_Name"]]

    sesame_seg_treatment.2 <- sesame.seg[names(sesame.seg) %in%
      sample.sheet.csv[sample.sheet.csv$gender_reported %in% split.by.2,
        "Sample_Name"]]

    sesame_seg_treatment.1.df <- binding_frames_mm(
      sesame_seg_treatement.1)

    sesame_seg_treatment.1.df <- sesame_seg_treatement.1.df[
      sesame_seg_treatment.1.df$pval <= 0.05,]

    sesame_seg_treatment.2.df <- binding_frames_mm(
      sesame_seg_treatment.2)

    sesame_seg_treatment.2.df <- sesame_seg_treatement.2.df[
      sesame_seg_treatment.2.df$pval <= 0.05,]

    colnames(sesame_seg_treatment.1.df)[1] <- "Sample_ID"
    seg.1 <- sesame_seg_treatment.1.df [,c(1,2,3,4,5,6,10)]
    colnames(sesame_seg_treatment.2.df)[1] <- "Sample_ID"
    seg.2 <- sesame_seg_treatment.2.df [,c(1,2,3,4,5,6,10)]
    seg$state <- round(2^seg$seg.mean * 2)
    seg$state[seg$state > 4] <- 4
    seg$method <- "sesame"
    row.names(seg) <- NULL
    seg <- na.omit(seg)
    seg$chrom <-
      unlist(strsplit(seg$chrom,split="chr"))[c(FALSE,TRUE)]

    seg.out <- list(seg.1,
                    seg.2)

    names(seg.out) <- c("internal.1", "internal.2")

    ##For sesame we want to omit the "chr"
    ##from chrom for intial plotting
    ##then add back in downstream for ggplot

  }

}else{

  if(is.null(split.by)){
    sesame.seg.treatment <- sesame.seg[names(sesame.seg) %in%
      sample.sheet.csv[sample.sheet.csv$Sample_Group %in% comparison[1],
        "Sample_Name"]]

    sesame.seg.control <- sesame.seg[names(sesame.seg) %in%
      sample.sheet.csv[sample.sheet.csv$Sample_Group %in% comparison[2],
        "Sample_Name"]]

    sesame_seg_treatment.df <- binding_frames_mm(
      sesame_seg_treatement)

    sesame_seg_treatment.df <- sesame_seg_treatment.df[
      sesame_seg_treatment.df$pval <= 0.05,]

    sesame_seg_control.df <- binding_frames_mm(
      sesame_seg_control)

    sesame_seg_control.df <- sesame_seg_control.df[
      sesame_seg_control.df$pval <= 0.05,]

    colnames(sesame_seg_treatment.df)[1] <- "Sample_ID"
    seg.treatment <- sesame_seg_treatment.df [,c(1,2,3,4,5,6,10)]
    seg.treatment$state <- round(2^seg$seg.mean * 2)
    seg.treatment$state[seg$state > 4] <- 4
    seg.treatment$method <- "sesame"
    row.names(seg.treatment) <- NULL
    seg.treatment <- na.omit(seg.treatment)
    seg.treatment$chrom <-
      unlist(strsplit(seg.treatment$chrom,split="chr"))[c(FALSE,TRUE)]

    colnames(sesame_seg_control.df)[1] <- "Sample_ID"
    seg.control <- sesame_seg_control.df [,c(1,2,3,4,5,6,10)]
    seg.control$state <- round(2^seg$seg.mean * 2)
    seg.control$state[seg$state > 4] <- 4
    seg.control$method <- "sesame"
    row.names(seg.control) <- NULL
    seg.control <- na.omit(seg.control)
    seg.control$chrom <-
      unlist(strsplit(seg.control$chrom,split="chr"))[c(FALSE,TRUE)]

    seg.out <- list(seg.treatment,
                    seg.control)

    names(seg.out) <- c(comparison[1], comparison[2])

    ##For sesame we want to omit the "chr"
    ##from chrom for intial plotting
    ##then add back in downstream for ggplot

  }else{

    if(is.null(split.by)){

      sesame.seg.treatment <- sesame.seg[names(sesame.seg) %in%
        sample.sheet.csv[sample.sheet.csv$Sample_Group %in% comparison[1],
          "Sample_Name"]]

      sesame.seg.control <- sesame.seg[names(sesame.seg) %in%
        sample.sheet.csv[sample.sheet.csv$Sample_Group %in% comparison[2],
          "Sample_Name"]]

      sesame_seg_treatment.df <- binding_frames_mm(
        sesame_seg_treatement)

      sesame_seg_treatment.df <- sesame_seg_treatment.df[
        sesame_seg_treatment.df$pval <= 0.05,]

      sesame_seg_control.df <- binding_frames_mm(
        sesame_seg_control)

      sesame_seg_control.df <- sesame_seg_control.df[
        sesame_seg_control.df$pval <= 0.05,]

      split.categories <- unique(sesame_seg_treatment.df[[split.by]])

      sesame_seg_treatment.df.1 <- sesame_seg_treatment.df[
        sesame_seg_treatment.df[[split.by]] %in% split.categories[1],]

      sesame_seg_treatment.df.2 <- sesame_seg_treatment.df[
        sesame_seg_treatment.df[[split.by]] %in% split.categories[2],]

      sesame_seg_control.df.1 <- sesame_seg_treatment.df[
        sesame_seg_treatment.df[[split.by]] %in% split.categories[1],]

      sesame_seg_control.df.2 <- sesame_seg_treatment.df[
        sesame_seg_treatment.df[[split.by]] %in% split.categories[2],]

      colnames(sesame_seg_treatment.df.1)[1] <- "Sample_ID"
      seg.treatment.1 <- sesame_seg_treatment.df.1[,c(1,2,3,4,5,6,10)]
      seg.treatment.1$state <- round(2^seg$seg.mean * 2)
      seg.treatment.1$state[seg$state > 4] <- 4
      seg.treatment.1$method <- "sesame"
      row.names(seg.treatment.1) <- NULL
      seg.treatment.1 <- na.omit(seg.treatment.1)
      seg.treatment.1$chrom <-
        unlist(strsplit(seg.treatment.1$chrom,split="chr"))[c(FALSE,TRUE)]

      colnames(sesame_seg_treatment.df.2)[1] <- "Sample_ID"
      seg.treatment.2 <- sesame_seg_treatment.df.2[,c(1,2,3,4,5,6,10)]
      seg.treatment.2$state <- round(2^seg$seg.mean * 2)
      seg.treatment.2$state[seg$state > 4] <- 4
      seg.treatment.2$method <- "sesame"
      row.names(seg.treatment.2) <- NULL
      seg.treatment.2 <- na.omit(seg.treatment.2)
      seg.treatment.2$chrom <-
        unlist(strsplit(seg.treatment.2$chrom,split="chr"))[c(FALSE,TRUE)]

      colnames(sesame_seg_control.df.1)[1] <- "Sample_ID"
      seg.control.1 <- sesame_seg_control.df.1[,c(1,2,3,4,5,6,10)]
      seg.control.1$state <- round(2^seg.1$seg.mean * 2)
      seg.control.1$state[seg.1$state > 4] <- 4
      seg.control.1$method <- "sesame"
      row.names(seg.control.1) <- NULL
      seg.control.1 <- na.omit(seg.control.1)
      seg.control.1$chrom <-
        unlist(strsplit(seg.control.1$chrom,split="chr"))[c(FALSE,TRUE)]

      colnames(sesame_seg_control.df.2)[1] <- "Sample_ID"
      seg.control.2 <- sesame_seg_control.df.2[,c(1,2,3,4,5,6,10)]
      seg.control.2$state <- round(2^seg.2$seg.mean * 2)
      seg.control.2$state[seg.2$state > 4] <- 4
      seg.control.2$method <- "sesame"
      row.names(seg.control.2) <- NULL
      seg.control.2 <- na.omit(seg.control.2)
      seg.control.2$chrom <-
        unlist(strsplit(seg.control.2$chrom,split="chr"))[c(FALSE,TRUE)]

      ##For sesame we want to omit the "chr"
      ##from chrom for intial plotting
      ##then add back in downstream for ggplot
      seg.out <- list(seg.treatment.1,
                      seg.treatment.2,
                      seg.control.1,
                      seg.control.2)

      names(seg.out) <- c(paste0(comparison[1],"_1"),
                          paste0(comparison[1],"_2"),
                          paste0(comparison[2],"_1"),
                          paste0(comparison[2],"_2"))

    }

  }

} ##End formating

return(seg.out)

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
