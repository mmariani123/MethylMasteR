#!/usr/bin/env Rscript

#' @title methyl_master_hm450
#' @description My version of hm450 analyses
#' ############ Samples and ref also need to be in same order ##########
## "CpG probe IDs not in the same order in data and ctrl!"
#' @param hm450.input.dir
#' @param hm450.output.dir
#' @param hm450.sample.sheet
#' @param hm450.reference
#' @param hm450.file.sep
#' @param hm450.comparison
#' @param hm450.split.by
#' @param hm450.sesame.data.cache
#' @param hm450.sesame.data.normal
#' @param hm450.sesame.ref.version
#' @param hm450.save.seg
#' @param ...
#' @import Biobase
#' @importFrom ExperimentHub setExperimentHubOption unbox
#' @importFrom ExperimentHub ExperimentHub unbox
#' @importFrom sesameData sesamesesameDataGet
#' @importFrom sesameData sesamesesameDataCache
#' @importFrom sesame SigSetsToRGChannel
#' @importFrom sesame openSesame
#' @return #seg.out
#' @export
methyl_master_hm450 <- function(hm450.idat.files.dir=NULL,
                               hm450.sample.sheet.path=NULL,
                               hm450.comparison=c("tumor","normal"),
                               hm450.split.by=NULL,
                               hm450.reference="internal",
                               hm450.sesame.data.cache="EPIC",
                               hm450.sesame.data.normal="EPIC.5.normal",
                               hm450.sesame.ref.version="hg38",
                               hm450.save.seg=FALSE,
                               ...
                               ){
  if(hm450.reference=="internal"){

    if(is.null(hm450.split.by){

    sesameData::sesameDataCache(sesame.data.cache)

    sesame_ssets_normal <- sesameData::sesameDataGet(sesame.data.normal)

    ExperimentHub::setExperimentHubOption("CACHE", idat.files.dir)

    ExperimentHub::ExperimentHub()

    treatment_idat_prefixes <- sesame::searchIDATprefixes(idat.files.dir,
                                                          recursive=TRUE)
    sesameData::sesameDataCacheAll()

    treatment.names <-    sesame.sesame.sample.sheet.df[
      sesame.sesame.sample.sheet.df$Sample_Group %in%
        sesame.comparison[1],"Sample_Name"]

    treatment_idat_prefixes <-
      treatment_idat_prefixes[gsub(paste0(".*",
                                          .Platform$file.sep),
                                          "",
                                          treatment_idat_prefixes) %in%
                                          treatment.names]

    treatment.platform <- unique(sesame.sesame.sample.sheet.df[
      sesame.sesame.sample.sheet.df$Sample_Name %in%
        treatment.names, "Platform"])

    hm450_sset <- sesame::openSesame(treatment_idat_prefixes,
                                      mask = TRUE,
                                      sum.TypeI = TRUE,
                                      platform = treatment.platform,
                                      what="sigset")

    hm450_rgset <- sesame::SigSetsToRGChannelSet(hm450_sset)

    hm450_cn_methylset_treatment <-
      minfi::getCN(minfi::preprocessRaw(hm450_rgset))

    save(hm450_cn_methylset_treatment,
         file=paste0(hm.450.output.dir,
                     .Platform$file.sep,
                     "hm450_cn_methylset_treatment.RData"))

    seg <- list(hm450_cn_methylset_treatment)

    }else{

      stop(paste0("Error: Epic.5.Normal internal ",
                  "\nreference samples are all male",
                  "\ncan't split by sex!"))

      split.by.cat <- unique(hm450.sample.sheet.df[[hm450.split.by]])

      sesameData::sesameDataCache(hm450.data.cache)

      hm450_ssets_normal <- sesameData::sesameDataGet(hm450.data.normal)

      normal.sexes <- unlist(lapply(hm450_ssets_normal,sesame::inferSex))

      hm450_ssets.normal.1 <-
        hm450_ssets_normal[names(which(normal.sexes==split.by.cat[1]))]

      hm450_ssets.normal.2 <-
        hm450_ssets_normal[names(which(normal.sexes==split.by.cat[2]))]

      ExperimentHub::setExperimentHubOption("CACHE", idat.files.dir)

      ExperimentHub::ExperimentHub()

      treatment_idat_prefixes <- sesame::searchIDATprefixes(idat.files.dir,
                                                            recursive=TRUE)
      sesameData::sesameDataCacheAll()

      hm450_sset.1 <- sesame::openSesame(treatment_idat_prefixes,
                                          mask = TRUE,
                                          sum.TypeI = TRUE,
                                          platform = hm450.sesame.platform,
                                          what="sigset")

      hm450_seg.1 <- foreach(i = 1:length(names(hm450_sset.1))) %do% {
        sesame::cnSegmentation(hm450_sset.1[[i]],
                               hm450_ssets_normal.1,
                               refversion = hm450.sesame.ref.version)
      }
      names(hm450_seg.1) <- names(hm450_sset.1)

      hm450_sset.2 <- sesame::openSesame(treatment_idat_prefixes,
                                          mask = TRUE,
                                          sum.TypeI = TRUE,
                                          platform = hm450.sesame.platform,
                                          what="sigset")

      hm450_seg.2 <- foreach(i = 1:length(names(hm450_sset.2))) %do% {
        sesame::cnSegmentation(hm450_sset.2[[i]],
                               hm450_ssets_normal.2,
                               refversion = hm450.sesame.ref.version)
      }
      names(hm450_seg.2) <- names(hm450_sset.2)

      seg <- list(hm450_seg.1, hm450_seg.2)

    }

  }else{

    if(is.null(hm450.split.by)){

    sample.sheet.df <- read.csv(hm450.sample.sheet.path,
                                header=TRUE,
                                stringsAsFactors = FALSE)

    subset.sample.sheet.df <- sample.sheet.df[sample.sheet.df %in%
                                                hm450.comparison,]

    subset.sample.sheet.df.normal <-
      sample.sheet.df[sample.sheet.df$Sample_Group %in%
                        hm450.comparison[2],]

    ExperimentHub::setExperimentHubOption("CACHE", idat.files.dir)

    ExperimentHub::ExperimentHub()

    idat_prefixes <- searchIDATprefixes(idat.files.dir, recursive=TRUE)

    sesameData::sesameDataCacheAll()

    treatment.names <- sesame.sample.sheet.df[
      sesame.sample.sheet.df$Sample_Group %in%
        sesame.comparison[1],"Sample_Name"]

    treatment_idat_prefixes <-  idat_prefixes[gsub(paste0(".*",
                                .Platform$file.sep),
                                "",
                                idat_prefixes) %in%
                                treatment.names]

    treatment.platform <- unique(sesame.sample.sheet.df[
      sesame.sample.sheet.df$Sample_Name %in%
        treatment.names, "Platform"])

    control_idat_prefixes <-  idat_prefixes[gsub(paste0(".*",
                                            .Platform$file.sep),
                                            "",
                                            idat_prefixes) %in%
                                            control.names]

    control.platform <- unique(sesame.sample.sheet.df[
      sesame.sample.sheet.df$Sample_Name %in%
        control.names, "Platform"])

    sesame_ssets_control <- openSesame(idat_prefixes.control,
                                   mask = TRUE,
                                   sum.TypeI = TRUE,
                                   platform = control.platform,
                                   what="sigset")

    sesame_ssets_treatment <- openSesame(idat_prefixes.treatment,
                                   mask = TRUE,
                                   sum.TypeI = TRUE,
                                   platform = treatment.platform,
                                   what="sigset")
    }else{

      split.by.cat <- unique(subset.sample.sheet.df[[hm450.split.by]])

      treatment.samples.1 <-
        sesame.sample.sheet.df[
          sesame.sample.sheet.df$Sample_Group==sesame.comparison[1] &
            sesame.sample.sheet.df[[hm450.split.by]]==split.by.cat[1],
            "Sample_Name"]

      control.samples.1 <-
        sesame.sample.sheet.df[
          sesame.sample.sheet.df$Sample_Group==sesame.comparison[2] &
            sesame.sample.sheet.df[[hm450.split.by]]==split.by.cat[1],
            "Sample_Name"]

      treatment.paths.1 <-
        sesame.sample.sheet.df[
          sesame.sample.sheet.df$Sample_Group==sesame.comparison[1] &
            sesame.sample.sheet.df[[hm450.split.by]]==split.by.cat[1],
            "Basename"]

      control.paths.1   <-
        sesame.sample.sheet.df[
          sesame.sample.sheet.df$Sample_Group==sesame.comparison[2] &
            sesame.sample.sheet.df[[hm450.split.by]]==split.by.cat[1],
            "Basename"]

      treatment.platform.1 <-
        unique(sesame.sample.sheet.df[
          sesame.sample.sheet.df$Sample_Group==sesame.comparison[1] &
            sesame.sample.sheet.df[[hm450.split.by]]==split.by.cat[2],
            "Platform"])

      control.platform.1 <-
        unique(sesame.sample.sheet.df[
          sesame.sample.sheet.df$Sample_Group==sesame.comparison[2] &
            sesame.sample.sheet.df[[hm450.split.by]]==split.by.cat[2],
            "Platform"])

      treatment.samples.2 <-
        sesame.sample.sheet.df[
          sesame.sample.sheet.df$Sample_Group==sesame.comparison[1] &
            sesame.sample.sheet.df[[hm450.split.by]]==split.by.cat[2],
            "Sample_Name"]

      control.samples.2 <-
        sesame.sample.sheet.df[
          sesame.sample.sheet.df$Sample_Group==sesame.comparison[2] &
            sesame.sample.sheet.df[[hm450.split.by]]==split.by.cat[2],
            "Sample_Name"]

      treatment.paths.2 <-
        sesame.sample.sheet.df[
          sesame.sample.sheet.df$Sample_Group==sesame.comparison[1] &
            sesame.sample.sheet.df[[hm450.split.by]]==split.by.cat[2],
            "Basename"]

      control.paths.2   <-
        sesame.sample.sheet.df[
          sesame.sample.sheet.df$Sample_Group==sesame.comparison[2] &
            sesame.sample.sheet.df[[hm450.split.by]]==split.by.cat[2],
            "Basename"]

      treatment.platform.2 <-
        unique(hm450.sample.sheet.df[
          hm450.sample.sheet.df$Sample_Group==hm450.comparison[1] &
            hm450.sample.sheet.df[[hm450.split.by]]==split.by.cat[2],
            "Platform"])

      control.platform.2   <-
        unique(hm450.sample.sheet.df[
          hm450.sample.sheet.df$Sample_Group==hm450.comparison[2] &
            hm450.sample.sheet.df[[hm450.split.by]]==split.by.cat[2],
            "Platform"])

      ExperimentHub::setExperimentHubOption("CACHE",
                                          hm450.idat.files.dir)

      ExperimentHub::ExperimentHub()

      idat_prefixes <- sesame::searchIDATprefixes(hm450.idat.files.dir,
                                                recursive=TRUE)

      idat_prefixes.treatment <-
        idat_prefixes[gsub(".*/","",idat_prefixes) %in%
                      gsub(paste0(".*",.Platform$file.sep),
                           "",treatment.paths)]

      idat_prefixes.control   <-
        idat_prefixes[gsub(".*/","",idat_prefixes) %in%
                      gsub(paste0(".*",.Platform$file.sep),
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

      sesame_rgset.treatment <-
        sesame::SigSetsToRGChannelSet(sesame_sset.treatment)

      sesame_rgset.control <-
        sesame::SigSetsToRGChannelSet(sesame_sset.control)

      sesame_cn_methylset.treatment <-
        minfi::getCN(minfi::preprocessRaw(sesame_rgset.treatment))

      sesame_cn_methylset.control <-
        minfi::getCN(minfi::preprocessRaw(sesame_rgset.control))

      if(hm450.workflow=="A"){

        proc_treatment_A.1  <-
          hm450_cn_methylset_treatment.1[,
                                      1:length(colnames(
                                      sesame_rgset_treatment.1))]

        proc_treatment_A.2    <-
          hm450_cn_methylset_treatment.2[,
                                      1:length(colnames(
                                      sesame_rgset_treatment.2))]

        proc_control_A.1 <-
          hm450_cn_methylset_control.1[,
                                    1:length(colnames(
                                    sesame_rgset_control.1))] ##4 Cols
        med_control_A.1  <- apply(proc_normal_female_A, 1, "median")

        proc_control_A.2   <-
          hm450_cn_methylset_control.2[,
                                    1:length(colnames(
                                    sesame_rgset_control.2))]
        med_control_A.2    <- apply(proc_control_A.2, 1, "median")

        candidates_data_treatment_A.1 <-
          findSegments2(proc_tumor_female_A[, , drop = FALSE],
                      med_normal_female_A,
                      proc_normal_female_A) ##Scaled with B but not A

        candidates_data_treatment_A.2 <-
          findSegments2(proc_tumor_male_A[, , drop = FALSE],
                      med_normal_male_A,
                      proc_normal_male_A)

      }else if(hm450.workflow=="B"){

      ## With z-Transformation, illumina
        proc_normal_female_B  <-
          hm450_cn_methylset_normal_female[,
                                         1:length(colnames(
                                           sesame_rgset_normal_female))]
        proc_normal_female_B[is.infinite(proc_normal_female_B)] <- NA
        proc_normal_female_B <- scale(proc_normal_female_B) ##Z transform
        med_normal_female_B <- apply(proc_normal_female_B, 1, "median")

        proc_normal_male_B <-
          hm450_cn_methylset_normal_male[,
                                       1:length(colnames(
                                         sesame_rgset_normal_male))]
        proc_normal_male_B[is.infinite(proc_normal_male_B)] <- NA
        proc_normal_male_B <- scale(proc_normal_male_B)
        med_normal_male_B <- apply(proc_normal_male_B, 1, "median")

        proc_tumor_female_B  <-
          hm450_cn_methylset_tumor_female[,
                                        1:length(colnames(
                                          sesame_rgset_tumor_female))]
        proc_tumor_female_B[is.infinite(proc_tumor_female_B)] <- NA
        proc_tumor_female_B <- scale(proc_tumor_female_B)
        med_tumor_female_B <- apply(proc_tumor_female_B, 1, "median")

        proc_tumor_male_B  <-
          hm450_cn_methylset_tumor_male[,
                                      1:length(colnames(
                                        sesame_rgset_tumor_male))]
        proc_tumor_male_B[is.infinite(proc_tumor_male_B)] <- NA
        proc_tumor_male_B <- scale(proc_tumor_male_B)
        med_tumor_male_B <- apply(proc_tumor_male_B, 1, "median")

        candidates_data_normal_female_B <-
          findSegments2(proc_tumor_female_B[, , drop = FALSE],
                      med_normal_female_B,
                      proc_normal_female_B) ##Scaled with B but not A

        candidates_data_normal_male_B <-
          findSegments2(proc_tumor_male_B[, , drop = FALSE],
                      med_normal_male_B,
                      proc_normal_male_B)

      }else if(hm450.workflow=="C"){

      ## Conumee-path, Illumina

        proc_normal_female_C <-
          hm450_cn_methylset_normal_female[,
                                         1:length(colnames(
                                           sesame_rgset_normal_female))]

        proc_normal_male_C   <-
          hm450_cn_methylset_normal_male[,
                                       1:length(colnames(
                                         sesame_rgset_normal_male))]

        proc_tumor_female_C <-
          hm450_cn_methylset_tumor_female[,
                                        1:length(colnames(
                                          sesame_rgset_tumor_female))]

        proc_tumor_male_C <-
          hm450_cn_methylset_tumor_male[,
                                      1:length(colnames(
                                        sesame_rgset_tumor_male))]

        proc_cord_female_C  <-
          hm450_cn_methylset_cord_female[,
                                       1:length(colnames(
                                         sesame_rgset_cord_female))]

        proc_cord_male_C    <-
          hm450_cn_methylset_cord_male[,
                                     1:length(colnames(
                                       sesame_rgset_cord_male))]

        proc_cord_female_present_C <-
          proc_cord_female_C[rownames(proc_cord_female_C) %in%
                             rownames(proc_tumor_female_C),]

        proc_cord_male_present_C <-
          proc_cord_male_C[rownames(proc_cord_male_C) %in%
                           rownames(proc_tumor_male_C),]

        proc_treatment_sorted_C <-
          proc_tumor_female_sorted_C[female_shared_names,]

        proc_control_sorted_C   <-
          proc_tumor_male_sorted_C[male_shared_names,]

      ####################### RUN CONUMEE ##############################

        candidates_data_C <-
          cnAnalysis450k::runConumee(proc_treatment_sorted_C,
                                   proc_control_sorted_C)

      }

    }

}











if(k450.workflow=="A"){

  ##sub routine A
  load(paste0(output.dir,
              file.sep,
              "candidates_data_normal_female_A.RData"))
  load(paste0(output.dir,
              file.sep,
              "candidates_data_normal_male_A.RData"))
  load(paste0(output.dir,
              file.sep,
              "candidates_data_cord_female_A.RData"))
  load(paste0(output.dir,
              file.sep,
              "candidates_data_cord_male_A.RData"))

  candidates_data_normal_female_A_sig <-
    candidates_data_normal_female_A[
      candidates_data_normal_female_A$p.val <= 0.05,]
  candidates_data_normal_male_A_sig   <-
    candidates_data_normal_male_A[
      candidates_data_normal_male_A$p.val <= 0.05,]
  candidates_data_cord_female_A_sig   <-
    candidates_data_cord_female_A[
      candidates_data_cord_female_A$p.val <= 0.05,]
  candidates_data_cord_male_A_sig   <-
    candidates_data_cord_male_A[
      candidates_data_cord_male_A$p.val <= 0.05,]

  candidates_data_normal_female_A_sig$num.mark <- NA
  candidates_data_normal_female_A_sig$bstat    <- NA
  candidates_data_normal_female_A_sig$treatment <- "tumor"
  candidates_data_normal_male_A_sig$num.mark <- NA
  candidates_data_normal_male_A_sig$bstat    <- NA
  candidates_data_normal_male_A_sig$treatment <- "tumor"
  candidates_data_cord_female_A_sig$num.mark <- NA
  candidates_data_cord_female_A_sig$bstat    <- NA
  candidates_data_cord_female_A_sig$treatment <- "tumor"
  candidates_data_cord_male_A_sig$num.mark <- NA
  candidates_data_cord_male_A_sig$bstat    <- NA
  candidates_data_cord_male_A_sig$treatment <- "tumor"

  preferred.columns <- c(7,1,2,3,12,13,8,5,4,9,10,11,14)

  candidates_data_normal_female_A_sig <-
    candidates_data_normal_female_A_sig[, preferred.columns]
  candidates_data_normal_male_A_sig <-
    candidates_data_normal_male_A_sig[, preferred.columns]
  candidates_data_cord_female_A_sig <-
    candidates_data_cord_female_A_sig[, preferred.columns]
  candidates_data_cord_male_A_sig <-
    candidates_data_cord_male_A_sig[, preferred.columns]

  candidates_data_normal_A_sig_colnames <- c("ID",
                                             "chrom",
                                             "loc.start",
                                             "loc.end",
                                             "num.mark",
                                             "bstat",
                                             "pval",
                                             "seg.mean",
                                             "seg.median",
                                             "karyotype",
                                             "sex_reported",
                                             "sex_inferred",
                                             "treatment")

  colnames(candidates_data_normal_female_A_sig) <-
    candidates_data_normal_A_sig_colnames
  colnames(candidates_data_normal_male_A_sig) <-
    candidates_data_normal_A_sig_colnames
  colnames(candidates_data_cord_female_A_sig) <-
    candidates_data_normal_A_sig_colnames
  colnames(candidates_data_cord_male_A_sig) <-
    candidates_data_normal_A_sig_colnames

  seg <- do.call(rbind,
                 list(candidates_data_normal_female_A_sig,
                      candidates_data_normal_male_A_sig,
                      candidates_data_cord_female_A_sig,
                      candidates_data_cord_male_A_sig)
  )

}else if(k450.workflow==B){

  load(paste0(output.dir,
              file.sep,
              "candidates_data_normal_female_B.RData"))
  load(paste0(output.dir,
              file.sep,
              "candidates_data_normal_male_B.RData"))
  load(paste0(output.dir,
              file.sep,
              "candidates_data_cord_female_B.RData"))
  load(paste0(output.dir,
              file.sep,
              "candidates_data_cord_male_B.RData"))

  candidates_data_normal_female_B_sig <-
    candidates_data_normal_female_B[
      candidates_data_normal_female_B$p.val <= 0.05,]
  candidates_data_normal_male_B_sig  <-
    candidates_data_normal_male_B[
      candidates_data_normal_male_B$p.val <= 0.05,]
  candidates_data_cord_female_B_sig  <-
    candidates_data_cord_female_B[
      candidates_data_cord_female_B$p.val <= 0.05,]
  candidates_data_cord_male_B_sig  <-
    candidates_data_cord_male_B[
      candidates_data_cord_male_B$p.val <= 0.05,]

  candidates_data_normal_female_B_sig$num.mark <- NA
  candidates_data_normal_female_B_sig$bstat    <- NA
  candidates_data_normal_female_B_sig$treatment <- "tumor"
  candidates_data_normal_male_B_sig$num.mark <- NA
  candidates_data_normal_male_B_sig$bstat    <- NA
  candidates_data_normal_male_B_sig$treatment <- "tumor"
  candidates_data_cord_female_B_sig$num.mark <- NA
  candidates_data_cord_female_B_sig$bstat    <- NA
  candidates_data_cord_female_B_sig$treatment <- "tumor"
  candidates_data_cord_male_B_sig$num.mark <- NA
  candidates_data_cord_male_B_sig$bstat    <- NA
  candidates_data_cord_male_B_sig$treatment <- "tumor"

  candidates_data_normal_female_B_sig <-
    candidates_data_normal_female_B_sig[, preferred.columns]
  candidates_data_normal_male_B_sig <-
    candidates_data_normal_male_B_sig[, preferred.columns]
  candidates_data_cord_female_B_sig <-
    candidates_data_cord_female_B_sig[, preferred.columns]
  candidates_data_cord_male_B_sig <-
    candidates_data_cord_male_B_sig[, preferred.columns]

  candidates_data_normal_B_sig_colnames <- c("ID",
                                             "chrom",
                                             "loc.start",
                                             "loc.end",
                                             "num.mark",
                                             "bstat",
                                             "pval",
                                             "seg.mean",
                                             "seg.median",
                                             "karyotype",
                                             "sex_reported",
                                             "sex_inferred",
                                             "treatment")

  colnames(candidates_data_normal_female_B_sig) <-
    candidates_data_normal_B_sig_colnames
  colnames(candidates_data_normal_male_B_sig) <-
    candidates_data_normal_B_sig_colnames
  colnames(candidates_data_cord_female_B_sig) <-
    candidates_data_normal_B_sig_colnames
  colnames(candidates_data_cord_male_B_sig) <-
    candidates_data_normal_B_sig_colnames

  seg <- do.call(rbind,
                 list(candidates_data_normal_female_B_sig,
                      candidates_data_normal_male_B_sig,
                      candidates_data_cord_female_B_sig,
                      candidates_data_cord_male_B_sig))

}else if(k450.workflow=="C"){

  load(paste0(output.dir,
              file.sep,
              "candidates_data_normal_female_C.RData"))
  load(paste0(output.dir,
              file.sep,
              "candidates_data_normal_male_C.RData"))
  load(paste0(output.dir,
              file.sep,
              "candidates_data_cord_female_C.RData"))
  load(paste0(output.dir,
              file.sep,
              "candidates_data_cord_male_C.RData"))

  candidates_data_normal_female_C_sig <-
    candidates_data_normal_female_C$data[
      candidates_data_normal_female_C$data$pval <= 0.05,]
  candidates_data_normal_male_C_sig    <-
    candidates_data_normal_male_C$data[
      candidates_data_normal_male_C$data$pval <= 0.05,]
  candidates_data_cord_female_C_sig    <-
    candidates_data_cord_female_C$data[
      candidates_data_cord_female_C$data$pval <= 0.05,]
  candidates_data_cord_male_C_sig    <-
    candidates_data_cord_male_C$data[
      candidates_data_cord_male_C$data$pval <= 0.05,]

  candidates_data_normal_female_C_sig$treatment <- "tumor"
  candidates_data_normal_male_C_sig$treatment <- "tumor"
  candidates_data_cord_female_C_sig$treatment <- "tumor"
  candidates_data_cord_male_C_sig$treatment <- "tumor"

  seg <- do.call(rbind,
                 list(candidates_data_normal_female_C_sig,
                      candidates_data_normal_male_C_sig,
                      candidates_data_cord_female_C_sig,
                      candidates_data_cord_male_C_sig))

}else{

  stop(paste0("Error: need to select a ",
              "proper sub workflow for 450k",
              " <k450.workflow>"))

}

##Note differences in 450k standard and conumee:
##colnames(candidates_data_cord_male_B_sig)
##"chr"
##"startCG"
##"endCG"
##"median"
##"mean"
##"sd"
##"smp"
##"p.val"
##"karyotype"
##"sex_reported"
##"sex_inferred"
##"num.mark"
##"bstat"

##candidates_data_cord_male_C_sig
##"ID"
##"chrom"
##"loc.start"
##"loc.end"
##"num.mark"
##"bstat"
##"pval"
##"seg.mean"
##"seg.median"
##"karyotype"
##"sex_reported"
##"sex_inferred"

##Desired names:
##colnames(seg)
##"ID"
##"chrom"
##"loc.start"
##"loc.end"
##"num.mark"
##"bstat"
##"pval"
##"seg.mean"
##"seg.median"
##karyotype
##sex_reported
##sex_inferred
##treatment

colnames(seg)[1] <- "Sample_ID"
seg <- seg[,c(1,2,3,4,5,8,10,11,12,13)]
seg$state <- round(2^seg$seg.mean * 2)
seg$state[seg$state > 4] <- 4
seg$method <- "k450"
seg$sub.method <- sub.workflow
row.names(seg) <- NULL
seg <- na.omit(seg) ##Workflow C ends up with some NA rows

return(seg.out)

}
