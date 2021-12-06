#!/usr/bin/env Rscript

#' @title My 450k
#' @description My version of 450k analyses
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
#' @import from sesameData sesamesesameDataGet
#' @import from sesameData sesamesesameDataCache
#' @import from sesame SigSetsToRGChannel
#' @import from sesame openSesame
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
    sesameData::sesameDataCache(hm450.sesame.data.cache)
    sesame_ssets_ref <- sesameData::sesameDataGet(hm450.sesame.data.normal)
  }else{
    sample.sheet.df <- read.csv(hm450.sample.sheet.path,
                                header=TRUE,
                                stringsAsFactors = FALSE)

    subset.sample.sheet.df <- sample.sheet.df[sample.sheet.df %in% comparison,]
    subset.sample.sheet.df.normal <-
      sample.sheet.df[sample.sheet.df$Sample_Group %in% comparisons[2],]
    setExperimentHubOption("CACHE", idat.files.dir)
    ExperimentHub()
    idat_prefixes <- searchIDATprefixes(idat.files.dir, recursive=TRUE)
    idat_prefixes.normal <- idat_prefixes[idat_prefixes %in%
                      subset.sample.sheet.df.normal$Sample_Name]
    idat_prefixes.treatment <- idat_prefixes[!idat_prefixes %in%
                            subset.sample.sheet.df.normal$Sample_Name]
    sesameDataCacheAll()
    sesame_ssets_ref <- openSesame(idat_prefixes.normal,
                                      mask = TRUE,
                                      sum.TypeI = TRUE,
                                      platform = sesame.platform,
                                      what="sigset")
    sesame_ssets_treatment <- openSesame(idat_prefixes.treatment,
                                   mask = TRUE,
                                   sum.TypeI = TRUE,
                                   platform = sesame.platform,
                                   what="sigset")
  }

  if(!is.null(split.by)){
    split.names <- unique(subset.sample.sheet.df[[split.by]])
      sesame_ssets_treatment.1 <- sesame_ssets_treatment[sesame_ssets_treatment ]
      sesame_ssets_treatment.2 <-
      sesame_ssets_control.1   <-
      sesame_ssets_control.2   <-
    idat_prefixes.treatment
  }
  setExperimentHubOption("CACHE", idat.files.dir)

  ExperimentHub()
  idat_prefixes <- searchIDATprefixes(idat.files.dir,recursive=TRUE)
  sesameDataCacheAll()

  sesame_sset <- sesame::openSesame(idat_prefixes,
                            mask = TRUE,
                            sum.TypeI = TRUE,
                            platform = sesame.platform,
                            what="sigset")

  sesame_rgset <- sesame::SigSetsToRGChannelSet(sesame_sset)

hm450_cn_methylset_normal_female <-
  minfi::getCN(minfi::preprocessRaw(sesame_rgset_normal_female))

hm450_cn_methylset_normal_male <-
  minfi::getCN(minfi::preprocessRaw(sesame_rgset_normal_male))

hm450_cn_methylset_tumor_female <-
  minfi::getCN(minfi::preprocessRaw(sesame_rgset_tumor_female))

hm450_cn_methylset_tumor_male <-
  minfi::getCN(minfi::preprocessRaw(sesame_rgset_tumor_male))

hm450_cn_methylset_cord_female <-
  minfi::getCN(minfi::preprocessRaw(sesame_rgset_cord_female))

hm450_cn_methylset_cord_male <-
  minfi::getCN(minfi::preprocessRaw(sesame_rgset_cord_male))

save(hm450_cn_methylset_normal_female,
     file=paste0(work.dir,file.sep,"hm450_cn_methylset_normal_female.RData"))

save(hm450_cn_methylset_normal_male,
     file=paste0(work.dir,file.sep,"hm450_cn_methylset_normal_male.RData"))

save(hm450_cn_methylset_tumor_female,
     file=paste0(work.dir,file.sep,"hm450_cn_methylset_tumor_female.RData"))

save(hm450_cn_methylset_tumor_male,
     file=paste0(work.dir,file.sep,"hm450_cn_methylset_tumor_male.RData"))

save(hm450_cn_methylset_cord_female,
     file=paste0(work.dir,file.sep,"hm450_cn_methylset_cord_female.RData"))

save(hm450_cn_methylset_cord_male,
     file=paste0(work.dir,file.sep,"hm450_cn_methylset_cord_male.RData"))

## Choose a workflow, options are:
## A, B (Z-transform), or C (conumee)

switch(hm450.workflow,

       A = {

         ## Without z-transformation, illumina
         ## For every tumor sample we are comparing to the median of
         ## the normals

         proc_tumor_female_A  <-
           hm450_cn_methylset_tumor_female[,
           1:length(colnames(
           sesame_rgset_tumor_female))]

         proc_tumor_male_A    <-
           hm450_cn_methylset_tumor_male[,
           1:length(colnames(
           sesame_rgset_tumor_male))]

         proc_normal_female_A <-
           hm450_cn_methylset_normal_female[,
           1:length(colnames(
           sesame_rgset_normal_female))] ##4 Cols
         med_normal_female_A  <- apply(proc_normal_female_A, 1, "median")

         proc_normal_male_A   <-
           hm450_cn_methylset_normal_male[,
           1:length(colnames(
           sesame_rgset_normal_male))]
         med_normal_male_A    <- apply(proc_normal_male_A, 1, "median")

         proc_cord_female_A   <-
           hm450_cn_methylset_cord_female[,
           1:length(colnames(
           sesame_rgset_cord_female))]
         med_cord_female_A    <- apply(proc_cord_female_A, 1, "median")

         proc_cord_male_A     <-
           hm450_cn_methylset_cord_male[,
           1:length(colnames(
           sesame_rgset_cord_male))]
         med_cord_male_A      <- apply(proc_cord_male_A, 1, "median")

         candidates_data_normal_female_A <-
           findSegments2(proc_tumor_female_A[, , drop = FALSE],
                         med_normal_female_A,
                         proc_normal_female_A) ##Scaled with B but not A
         candidates_data_normal_female_A$karyotype    <- ""
         candidates_data_normal_female_A$sex_reported <- ""
         candidates_data_normal_female_A$sex_inferred <- ""
         for(i in 1:length(unique(candidates_data_normal_female_A$smp))){
           sample.now <- unique(candidates_data_normal_female_A$smp)[i]
           if(sample.now %in% rownames(tumor)){
             candidates_data_normal_female_A[
             candidates_data_normal_female_A$smp==sample.now,"karyotype"] <-
              tumor[sample.now, "karyotype"]
             candidates_data_normal_female_A[
             candidates_data_normal_female_A$smp==sample.now,"sex_reported"] <-
               tumor[sample.now, "gender_reported"]
             candidates_data_normal_female_A[
             candidates_data_normal_female_A$smp==sample.now,"sex_inferred"] <-
               tumor[sample.now, "sex_inferred"]
           }
         }

         candidates_data_normal_male_A <-
           findSegments2(proc_tumor_male_A[, , drop = FALSE],
                         med_normal_male_A,
                         proc_normal_male_A)
         candidates_data_normal_male_A$karyotype    <- ""
         candidates_data_normal_male_A$sex_reported <- ""
         candidates_data_normal_male_A$sex_inferred <- ""
         for(i in 1:length(unique(candidates_data_normal_male_A$smp))){
           sample.now <- unique(candidates_data_normal_male_A$smp)[i]
           if(sample.now %in% rownames(tumor)){
             candidates_data_normal_male_A[
               candidates_data_normal_male_A$smp==sample.now,"karyotype"] <-
               tumor[sample.now, "karyotype"]
             candidates_data_normal_male_A[
               candidates_data_normal_male_A$smp==sample.now,"sex_reported"] <-
               tumor[sample.now, "gender_reported"]
             candidates_data_normal_male_A[
               candidates_data_normal_male_A$smp==sample.now,"sex_inferred"] <-
               tumor[sample.now, "sex_inferred"]
           }
         }

         candidates_data_cord_female_A <-
           findSegments2(proc_tumor_female_A[, , drop = FALSE],
                         med_cord_female_A,
                         proc_cord_female_A)
         candidates_data_cord_female_A$karyotype    <- ""
         candidates_data_cord_female_A$sex_reported <- ""
         candidates_data_cord_female_A$sex_inferred <- ""
         for(i in 1:length(unique(candidates_data_cord_female_A$smp))){
           sample.now <- unique(candidates_data_cord_female_A$smp)[i]
           if(sample.now %in% rownames(tumor)){
             candidates_data_cord_female_A[
               candidates_data_cord_female_A$smp==sample.now,"karyotype"] <-
               tumor[sample.now, "karyotype"]
             candidates_data_cord_female_A[
               candidates_data_cord_female_A$smp==sample.now,"sex_reported"] <-
               tumor[sample.now, "gender_reported"]
             candidates_data_cord_female_A[
               candidates_data_cord_female_A$smp==sample.now,"sex_inferred"] <-
               tumor[sample.now, "sex_inferred"]
           }
         }

         candidates_data_cord_male_A <-
           findSegments2(proc_tumor_male_A[, , drop = FALSE],
                         med_cord_male_A,
                         proc_cord_male_A)
         candidates_data_cord_male_A$karyotype    <- ""
         candidates_data_cord_male_A$sex_reported <- ""
         candidates_data_cord_male_A$sex_inferred <- ""
         for(i in 1:length(unique(candidates_data_cord_male_A$smp))){
           sample.now <- unique(candidates_data_cord_male_A$smp)[i]
           if(sample.now %in% rownames(tumor)){
             candidates_data_cord_male_A[
               candidates_data_cord_male_A$smp==sample.now,"karyotype"] <-
               tumor[sample.now, "karyotype"]
             candidates_data_cord_male_A[
               candidates_data_cord_male_A$smp==sample.now,"sex_reported"] <-
               tumor[sample.now, "gender_reported"]
             candidates_data_cord_male_A[
               candidates_data_cord_male_A$smp==sample.now,"sex_inferred"] <-
               tumor[sample.now, "sex_inferred"]
           }
         }

         save(candidates_data_normal_female_A,
              file=paste0(work.dir,
                          file.sep,
                          "candidates_data_normal_female_A.RData"))

         save(candidates_data_normal_male_A,
              file=paste0(work.dir,
                          file.sep,
                          "candidates_data_normal_male_A.RData"))

         save(candidates_data_cord_female_A,
              file=paste0(work.dir,
                          file.sep,
                          "candidates_data_cord_female_A.RData"))

         save(candidates_data_cord_male_A,
              file=paste0(work.dir,
                          file.sep,
                          "candidates_data_cord_male_A.RData"))

         rm(proc_normal_female_A)
         rm(proc_normal_male_A)
         rm(proc_tumor_female_A)
         rm(proc_tumor_male_A)
         rm(proc_cord_female_A)
         rm(proc_cord_male_A)
         rm(candidates_data_normal_female_A)
         rm(candidates_data_normal_male_A)
         rm(candidates_data_cord_female_A)
         rm(candidates_data_cord_male_A)

       },

       B = {

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

         proc_cord_female_B <-
           hm450_cn_methylset_cord_female[,
           1:length(colnames(
           sesame_rgset_cord_female))]
         proc_cord_female_B[is.infinite(proc_cord_female_B)] <- NA
         proc_cord_female_B <- scale(proc_cord_female_B)
         med_cord_female_B <- apply(proc_cord_female_B, 1, "median")

         proc_cord_male_B <-
           hm450_cn_methylset_cord_male[,
           1:length(colnames(
           sesame_rgset_cord_male))]
         proc_cord_male_B[is.infinite(proc_cord_male_B)] <- NA
         proc_cord_male_B <- scale(proc_cord_male_B)
         med_cord_male_B <- apply(proc_cord_male_B, 1, "median")

         candidates_data_normal_female_B <-
           findSegments2(proc_tumor_female_B[, , drop = FALSE],
                         med_normal_female_B,
                         proc_normal_female_B) ##Scaled with B but not A
         candidates_data_normal_female_B$karyotype    <- ""
         candidates_data_normal_female_B$sex_reported <- ""
         candidates_data_normal_female_B$sex_inferred <- ""
         for(i in 1:length(unique(candidates_data_normal_female_B$smp))){
           sample.now <- unique(candidates_data_normal_female_B$smp)[i]
           if(sample.now %in% rownames(tumor)){
             candidates_data_normal_female_B[
             candidates_data_normal_female_B$smp==sample.now,"karyotype"] <-
               tumor[sample.now, "karyotype"]
             candidates_data_normal_female_B[
             candidates_data_normal_female_B$smp==sample.now,"sex_reported"] <-
               tumor[sample.now, "gender_reported"]
             candidates_data_normal_female_B[
             candidates_data_normal_female_B$smp==sample.now,"sex_inferred"] <-
               tumor[sample.now, "sex_inferred"]
           }
         }

         candidates_data_normal_male_B <-
           findSegments2(proc_tumor_male_B[, , drop = FALSE],
                         med_normal_male_B,
                         proc_normal_male_B)
         candidates_data_normal_male_B$karyotype    <- ""
         candidates_data_normal_male_B$sex_reported <- ""
         candidates_data_normal_male_B$sex_inferred <- ""
         for(i in 1:length(unique(candidates_data_normal_male_B$smp))){
           sample.now <- unique(candidates_data_normal_male_B$smp)[i]
           if(sample.now %in% rownames(tumor)){
             candidates_data_normal_male_B[
               candidates_data_normal_male_B$smp==sample.now,"karyotype"] <-
               tumor[sample.now, "karyotype"]
             candidates_data_normal_male_B[
               candidates_data_normal_male_B$smp==sample.now,"sex_reported"] <-
               tumor[sample.now, "gender_reported"]
             candidates_data_normal_male_B[
               candidates_data_normal_male_B$smp==sample.now,"sex_inferred"] <-
               tumor[sample.now, "sex_inferred"]
           }
         }

         candidates_data_cord_female_B <-
           findSegments2(proc_tumor_female_B[, , drop = FALSE],
                         med_cord_female_B,
                         proc_cord_female_B)
         candidates_data_cord_female_B$karyotype    <- ""
         candidates_data_cord_female_B$sex_reported <- ""
         candidates_data_cord_female_B$sex_inferred <- ""
         for(i in 1:length(unique(candidates_data_cord_female_B$smp))){
           sample.now <- unique(candidates_data_cord_female_B$smp)[i]
           if(sample.now %in% rownames(tumor)){
             candidates_data_cord_female_B[
               candidates_data_cord_female_B$smp==sample.now,"karyotype"] <-
               tumor[sample.now, "karyotype"]
             candidates_data_cord_female_B[
               candidates_data_cord_female_B$smp==sample.now,"sex_reported"] <-
               tumor[sample.now, "gender_reported"]
             candidates_data_cord_female_B[
               candidates_data_cord_female_B$smp==sample.now,"sex_inferred"] <-
               tumor[sample.now, "sex_inferred"]
           }
         }

         candidates_data_cord_male_B <-
           findSegments2(proc_tumor_male_B[, , drop = FALSE],
                         med_cord_male_B,
                         proc_cord_male_B)
         candidates_data_cord_male_B$karyotype    <- ""
         candidates_data_cord_male_B$sex_reported <- ""
         candidates_data_cord_male_B$sex_inferred <- ""
         for(i in 1:length(unique(candidates_data_cord_male_B$smp))){
           sample.now <- unique(candidates_data_cord_male_B$smp)[i]
           if(sample.now %in% rownames(tumor)){
             candidates_data_cord_male_B[
               candidates_data_cord_male_B$smp==sample.now,"karyotype"] <-
               tumor[sample.now, "karyotype"]
             candidates_data_cord_male_B[
               candidates_data_cord_male_B$smp==sample.now,"sex_reported"] <-
               tumor[sample.now, "gender_reported"]
             candidates_data_cord_male_B[
               candidates_data_cord_male_B$smp==sample.now,"sex_inferred"] <-
               tumor[sample.now, "sex_inferred"]
           }
         }

         save(candidates_data_normal_female_B,
              file=paste0(work.dir,
                          file.sep,
                          "candidates_data_normal_female_B.RData"))

         save(candidates_data_normal_male_B,
              file=paste0(work.dir,
                          file.sep,
                          "candidates_data_normal_male_B.RData"))

         save(candidates_data_cord_female_B,
              file=paste0(work.dir,
                          file.sep,
                          "candidates_data_cord_female_B.RData"))

         save(candidates_data_cord_male_B,
              file=paste0(work.dir,
                          file.sep,
                          "candidates_data_cord_male_B.RData"))

         rm(proc_normal_female_B)
         rm(proc_normal_male_B)
         rm(proc_tumor_female_B)
         rm(proc_tumor_male_B)
         rm(proc_cord_female_B)
         rm(proc_cord_male_B)
         rm(candidates_data_normal_female_B)
         rm(candidates_data_normal_male_B)
         rm(candidates_data_cord_female_B)
         rm(candidates_data_cord_male_B)

       },

       C = {

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

         ## Conumee:
         ## What information is of interest?
         ## Add the genes of interest
         ## Can select which feature of interest we want to look at
         ## Calculate segments with conumee

         ## MM Note: from warnings() "All infinite values are set to NA!"
         ## MM Note: cord data is from epic and tumor data is 450k

         ##proc_tumor_female_C %>% head(n=10)
         ##proc_cord_female_C  %>% head(n=10)

         ##rownames(proc_tumor_female_C) %in% rownames(proc_cord_female_C)
         ## Quite a few probes are present, excellent.
         ## Let's subset out the 450k from the EPIC cord ref
         ## to use in the analysis.

         proc_cord_female_present_C <-
           proc_cord_female_C[rownames(proc_cord_female_C) %in%
                                rownames(proc_tumor_female_C),]

         proc_cord_male_present_C <-
           proc_cord_male_C[rownames(proc_cord_male_C) %in%
                              rownames(proc_tumor_male_C),]

         ############ Samples and ref also need to be in same order ##########
         ## "CpG probe IDs not in the same order in data and ctrl!"

         proc_tumor_female_sorted_C = proc_tumor_female_C[order(
           rownames(proc_tumor_female_C)),]

         proc_tumor_male_sorted_C   = proc_tumor_male_C[order(
           rownames(proc_tumor_male_C)),]

         proc_cord_female_sorted_C = proc_cord_female_C[order(
           rownames(proc_cord_female_C)),]

         proc_cord_male_sorted_C   = proc_cord_male_C[order(
           rownames(proc_cord_male_C)),]

         female_shared_names <- intersect(rownames(proc_cord_female_sorted_C),
                                          rownames(proc_tumor_female_sorted_C))

         male_shared_names   <- intersect(rownames(proc_cord_male_sorted_C),
                                          rownames(proc_tumor_male_sorted_C))

         proc_tumor_female_sorted_C <-
           proc_tumor_female_sorted_C[female_shared_names,]

         proc_tumor_male_sorted_C   <-
           proc_tumor_male_sorted_C[male_shared_names,]

         proc_cord_female_sorted_C  <-
           proc_cord_female_sorted_C[female_shared_names,]

         proc_cord_male_sorted_C    <-
           proc_cord_male_sorted_C[male_shared_names,]

         ## Check ordering now:
         ##proc_tumor_female_sorted_C  %>% head(n=10)
         ##proc_cord_female_sorted_C   %>% head(n=10)

         ####################### RUN CONUMEE ##############################

         candidates_data_normal_female_C <-
           cnAnalysis450k::runConumee(proc_tumor_female_C,
                                      proc_normal_female_C)
         candidates_data_normal_female_C$data$karyotype    <- ""
         candidates_data_normal_female_C$data$sex_reported <- ""
         candidates_data_normal_female_C$data$sex_inferred <- ""
         for(i in 1:length(unique(candidates_data_normal_female_C$data$ID))){
           sample.now <- unique(candidates_data_normal_female_C$data$ID)[i]
           if(sample.now %in% rownames(tumor)){
             candidates_data_normal_female_C$data[
               candidates_data_normal_female_C$data$ID==sample.now,
               "karyotype"] <- tumor[sample.now, "karyotype"]
             candidates_data_normal_female_C$data[
               candidates_data_normal_female_C$data$ID==sample.now,
               "sex_reported"] <- tumor[sample.now, "gender_reported"]
             candidates_data_normal_female_C$data[
               candidates_data_normal_female_C$data$ID==sample.now,
               "sex_inferred"] <- tumor[sample.now, "sex_inferred"]
           }
         }

         candidates_data_normal_male_C   <-
           cnAnalysis450k::runConumee(proc_tumor_male_C,
                                      proc_normal_male_C)
         candidates_data_normal_male_C$data$karyotype    <- ""
         candidates_data_normal_male_C$data$sex_reported <- ""
         candidates_data_normal_male_C$data$sex_inferred <- ""
         for(i in 1:length(unique(candidates_data_normal_male_C$data$ID))){
           sample.now <- unique(candidates_data_normal_male_C$data$ID)[i]
           if(sample.now %in% rownames(tumor)){
             candidates_data_normal_male_C$data[
               candidates_data_normal_male_C$data$ID==sample.now,
               "karyotype"] <- tumor[sample.now, "karyotype"]
             candidates_data_normal_male_C$data[
               candidates_data_normal_male_C$data$ID==sample.now,
               "sex_reported"] <- tumor[sample.now, "gender_reported"]
             candidates_data_normal_male_C$data[
               candidates_data_normal_male_C$data$ID==sample.now,
               "sex_inferred"] <- tumor[sample.now, "sex_inferred"]
           }
         }

         candidates_data_cord_female_C <-
           cnAnalysis450k::runConumee(proc_tumor_female_sorted_C,
                                      proc_cord_female_sorted_C)
         candidates_data_cord_female_C$data$karyotype    <- ""
         candidates_data_cord_female_C$data$sex_reported <- ""
         candidates_data_cord_female_C$data$sex_inferred <- ""
         for(i in 1:length(unique(candidates_data_cord_female_C$data$ID))){
           sample.now <- unique(candidates_data_cord_female_C$data$ID)[i]
           if(sample.now %in% rownames(tumor)){
             candidates_data_cord_female_C$data[
               candidates_data_cord_female_C$data$ID==sample.now,
               "karyotype"] <- tumor[sample.now, "karyotype"]
             candidates_data_cord_female_C$data[
               candidates_data_cord_female_C$data$ID==sample.now,
               "sex_reported"] <- tumor[sample.now, "gender_reported"]
             candidates_data_cord_female_C$data[
               candidates_data_cord_female_C$data$ID==sample.now,
               "sex_inferred"] <- tumor[sample.now, "sex_inferred"]
           }
         }

         candidates_data_cord_male_C   <-
           cnAnalysis450k::runConumee(proc_tumor_male_sorted_C,
                                      proc_cord_male_sorted_C)
         candidates_data_cord_male_C$data$karyotype    <- ""
         candidates_data_cord_male_C$data$sex_reported <- ""
         candidates_data_cord_male_C$data$sex_inferred <- ""
         for(i in 1:length(unique(candidates_data_cord_male_C$data$ID))){
           sample.now <- unique(candidates_data_cord_male_C$data$ID)[i]
           if(sample.now %in% rownames(tumor)){
             candidates_data_cord_male_C$data[
               candidates_data_cord_male_C$data$ID==sample.now,
               "karyotype"] <- tumor[sample.now, "karyotype"]
             candidates_data_cord_male_C$data[
               candidates_data_cord_male_C$data$ID==sample.now,
               "sex_reported"] <- tumor[sample.now, "gender_reported"]
             candidates_data_cord_male_C$data[
               candidates_data_cord_male_C$data$ID==sample.now,
               "sex_inferred"] <- tumor[sample.now, "sex_inferred"]
           }
         }

         save(candidates_data_normal_female_C,
              file=paste0(work.dir,
                          file.sep,
                          "candidates_data_normal_female_C.RData"))

         save(candidates_data_normal_male_C,
              file=paste0(work.dir,
                          file.sep,
                          "candidates_data_normal_male_C.RData"))

         save(candidates_data_cord_female_C,
              file=paste0(work.dir,
                          file.sep,
                          "candidates_data_cord_female_C.RData"))

         save(candidates_data_cord_male_C,
              file=paste0(work.dir,
                          file.sep,
                          "candidates_data_cord_male_C.RData"))

         rm(proc_normal_female_C)
         rm(proc_normal_male_C)
         rm(proc_tumor_female_C)
         rm(proc_tumor_male_C)
         rm(proc_cord_female_C)
         rm(proc_cord_male_C)
         rm(candidates_data_normal_female_C)
         rm(candidates_data_normal_male_C)
         rm(candidates_data_cord_female_C)
         rm(candidates_data_cord_male_C)

       }

)

rm(sesame_rgset_normal_female)
rm(sesame_rgset_normal_male)
rm(sesame_rgset_tumor_female)
rm(sesame_rgset_tumor_male)
rm(sesame_rgset_cord_female)
rm(sesame_rgset_cord_male)

rm(hm450_cn_methylset_normal_female)
rm(hm450_cn_methylset_normal_male)
rm(hm450_cn_methylset_tumor_female)
rm(hm450_cn_methylset_tumor_male)
rm(hm450_cn_methylset_cord_female)
rm(hm450_cn_methylset_cord_male)

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
