#!/usr/bin/env Rscript

#' @title methyl_master_formatting_sesame
#' @description Main function:
#' MultiMethylv1.0 : CNV calling from methylation data
#' Michael Mariani PhD Dartmouth College 2021
#' Possible routines:
#' "test"            ##Run a quick test
#' "process_sesame", ##Preprocess the TCGA and cord data in sesame format
#' "sesame",         ##Run Sesame  CNV calling (get segments)
#' "hm450" ,         ##Run 450K    CNV calling (get segments)
#' "champ" ,         ##Run ChAMP   CNV calling (get segments)
#' "epicopy" ,       ##Run EpiCopy CNV calling (get segments)
#' "compare"         ##Run algorithm comparison functionality
#' @param sesame.form.seg
#' @param sesame.form.output.dir
#' @param output.form.dir
#' @param sample.form.sample.sheet.path
#' @param sesame.form.reference
#' @param sesame.form.split.by
#' @param sesame.form.comparison
#' @param sesame.form.save.seg
#' @param sesame.form.plot.individual
#' @param ...
#' @import CNVRanger
#' @import matter
#' @importFrom magrittr %>%
#' @return #Formatted seg list object for visualizing
#' @export
methyl_master_formatting_sesame <- function(sesame.form.seg,
                                            sesame.form.output.dir,
                                            sesame.form.sample.sheet.path,
                                            sesame.form.reference,
                                            sesame.form.split.by,
                                            sesame.form.comparison,
                                            sesame.form.save.seg,
                                            sesame.form.plot.individual,
                                            ...){

##if(sesame.form.save.seg==TRUE){
##  save(sesame.form.seg, file =  paste0(sesame.form.output.dir,
##                                       .Platform$file.sep,
##                                       "seg.RData"))
##}

if(sesame.form.plot.individual==TRUE){

  for(i in 1:length(sesame.form.seg)){
    methyl_master_plot_individual(pi.seg = sesame.form.seg[[i]],
                                  pi.output.dir = sesame.form.output.dir,
                                  pi.name = as.character(i)
                                  )
  }

}

if(sesame.form.reference=="internal"){

  if(is.null(sesame.form.split.by)){

    sesame.seg <- sesame.form.seg[[1]]
    rm(sesame.form.seg)

    sesame_seg_treatment.df <- binding_frames_mm(
      sesame.seg)

    sesame_seg_treatment.df <- sesame_seg_treatment.df[
      sesame_seg_treatment.df$pval <= 0.05 &
        !is.na(sesame_seg_treatment.df$pval),]

    ##sesame_seg_treatment.df %>% colnames
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

    seg <- sesame_seg_treatment.df
    rm(sesame_seg_treatment.df)

    colnames(seg)[1] <- "Sample_ID"

    seg$chrom <-
      unlist(strsplit(seg$chrom,split="chr"))[c(FALSE,TRUE)]

    seg$bstat      <- NA
    seg$state <- round(2^seg$seg.mean * 2)
    seg$state[seg$state > 4] <- 4
    seg$treatment  <- sesame.form.comparison[1]
    seg$method <- "sesame"
    seg$sub.method <- NA
    row.names(seg) <- NULL

    ##seg <- na.omit(seg) ##Workflow C ends up with some NA rows
    ##Above will remove rows with ANY NAs present!

    preferred.column.names <- c("Sample_ID",
                                "chrom",
                                "loc.start",
                                "loc.end",
                                "num.mark",
                                "bstat",
                                "seg.mean",
                                "seg.median",
                                "pval",
                                "state",
                                "treatment",
                                "method",
                                "sub.method")

    seg <- seg[,preferred.column.names]

    seg.out <- list(seg)

    names(seg.out) <- sesame.form.comparison[1]

  }else{

    ##sesame_seg_treatment <- sesame.seg[names(sesame.seg) %in%
    ##sesame.sample.sheet.df[sesame.sample.sheet.df$Sample_Group %in%
    ##                         sesame.comparison[1], "Sample_Name"]]

    sesame_seg_treatment_1 <- sesame.form.seg[[1]]
    sesame_seg_treatment_2 <- sesame.form.seg[[2]]

    rm(sesame.form.seg)

    split.by.1 <- unique(sesame.sample.sheet.df[[sesame.form.split.by]])[1]
    split.by.2 <- unique(sesame.sample.sheet.df[[sesame.form.split.by]])[2]

    ##sesame_seg_treatment.1 <- sesame.seg.treatment[
    ##  names(sesame.seg.treatment) %in%
    ##    sesame.sample.sheet.df[sesame.sample.sheet.df$gender_reported %in%
    ##                             split.by.1, "Sample_Name"]]

    ##sesame_seg_treatment.2 <- sesame.seg[names(sesame.seg) %in%
    ##  sesame.sample.sheet.df[sesame.sample.sheet.df$gender_reported %in%
    ##                           split.by.2, "Sample_Name"]]

    sesame_seg_treatment.1.df <- binding_frames_mm(
      sesame_seg_treatment.1)

    sesame_seg_treatment.1.df <- sesame_seg_treatment.1.df[
      sesame_seg_treatment.1.df$pval <= 0.05,]

    sesame_seg_treatment.2.df <- binding_frames_mm(
      sesame_seg_treatment.2)

    sesame_seg_treatment.2.df <- sesame_seg_treatment.2.df[
      sesame_seg_treatment.2.df$pval <= 0.05,]

    colnames(sesame_seg_treatment.1.df)[1] <- "Sample_ID"
    seg.1 <- sesame_seg_treatment.1.df [,c(1,2,3,4,5,6,10)]
    seg.1$state <- round(2^seg.1$seg.mean * 2)
    seg.1$state[seg.1$state > 4] <- 4
    seg.1$method <- "sesame"
    row.names(seg.1) <- NULL
    seg.1 <- na.omit(seg.1)
    seg.1$chrom <-
      unlist(strsplit(seg.1$chrom,split="chr"))[c(FALSE,TRUE)]

    colnames(sesame_seg_treatment.2.df)[1] <- "Sample_ID"
    seg.2 <- sesame_seg_treatment.2.df [,c(1,2,3,4,5,6,10)]
    seg.2$state <- round(2^seg.2$seg.mean * 2)
    seg.2$state[seg.2$state > 4] <- 4
    seg.2$method <- "sesame"
    row.names(seg.2) <- NULL
    seg.2 <- na.omit(seg.2)
    seg.2$chrom <-
      unlist(strsplit(seg.2$chrom,split="chr"))[c(FALSE,TRUE)]

    seg.out <- list(seg.1,
                    seg.2)

    names(seg.out) <- c(paste0(sesame.form.comparison[1],"_1"),
                        paste0(sesame.form.comparison[1],"_2"))

    ##For sesame we want to omit the "chr"
    ##from chrom for intial plotting
    ##then add back in downstream for ggplot

  }

}else{

  if(is.null(sesame.form.split.by)){

    ##sesame_seg_treatment <- sesame.seg[names(sesame.seg) %in%
    ##  sesame.sample.sheet.df[sesame.sample.sheet.df$Sample_Group %in%
    ##                           sesame.comparison[1], "Sample_Name"]]

    sesame_seg_treatment <- sesame.form.seg
    rm(sesame.form.seg)

    sesame_seg_treatment.df <- binding_frames_mm(
      sesame_seg_treatment)

    sesame_seg_treatment.df <- sesame_seg_treatment.df[
      sesame_seg_treatment.df$pval <= 0.05,]

    colnames(sesame_seg_treatment.df)[1] <- "Sample_ID"
    seg.treatment <- sesame_seg_treatment.df [,c(1,2,3,4,5,6,10)]
    seg.treatment$state <- round(2^seg.treatment$seg.mean * 2)
    seg.treatment$state[seg.treatment$state > 4] <- 4
    seg.treatment$method <- "sesame"
    row.names(seg.treatment) <- NULL
    seg.treatment <- na.omit(seg.treatment)
    seg.treatment$chrom <-
      unlist(strsplit(seg.treatment$chrom,split="chr"))[c(FALSE,TRUE)]

    seg.out <- list(seg.treatment)

    names(seg.out) <- sesame.form.comparison[1]

    ##For sesame we want to omit the "chr"
    ##from chrom for intial plotting
    ##then add back in downstream for ggplot

  }else{

    sesame_seg_treatment_1 <- sesame.form.seg[1]

    sesame_seg_treatment.df.1 <- binding_frames_mm(
      sesame_seg_treatment_1)
    sesame_seg_treatment.df.1 <- sesame_seg_treatment.df.1[
      sesame_seg_treatment.df.1$pval <= 0.05,]

    sesame_seg_treatment_2 <- sesame.form.seg[2]
    sesame_seg_treatment.df.2 <- binding_frames_mm(
      sesame_seg_treatment_2)
    sesame_seg_treatment.df.2 <- sesame_seg_treatment.df.2[
      sesame_seg_treatment.df.2$pval <= 0.05,]

    rm(sesame.form.seg)

    ##sesame_seg_treatment <- sesame.seg[names(sesame.seg) %in%
    ##sesame.sample.sheet.df[sesame.sample.sheet.df$Sample_Group %in%
    ##                           sesame.comparison[1], "Sample_Name"]]

    ##sesame_seg_treatment.df <- binding_frames_mm(
    ##    sesame_seg_treatment)

    ##sesame_seg_treatment.df <- sesame_seg_treatment.df[
    ##    sesame_seg_treatment.df$pval <= 0.05,]

    split.categories <- unique(sesame.sample.sheet.df[[sesame.form.split.by]])

    ##sesame_seg_treatment.df.1 <- sesame_seg_treatment.df[
    ##    sesame.sample.sheet.df[[sesame.split.by]] %in% split.categories[1],]

    ##sesame_seg_treatment.df.2 <- sesame_seg_treatment.df[
    ##    sesame.sample.sheet.df[[sesame.split.by]] %in% split.categories[2],]

    colnames(sesame_seg_treatment.df.1)[1] <- "Sample_ID"
    seg.treatment.1 <- sesame_seg_treatment.df.1[,c(1,2,3,4,5,6,10)]
    seg.treatment.1$state <- round(2^seg.treatment.1$seg.mean * 2)
    seg.treatment.1$state[seg.treatment.1$state > 4] <- 4
    seg.treatment.1$method <- "sesame"
    row.names(seg.treatment.1) <- NULL
    seg.treatment.1 <- na.omit(seg.treatment.1)
    seg.treatment.1$chrom <-
      unlist(strsplit(seg.treatment.1$chrom,split="chr"))[c(FALSE,TRUE)]

    colnames(sesame_seg_treatment.df.2)[1] <- "Sample_ID"
    seg.treatment.2 <- sesame_seg_treatment.df.2[,c(1,2,3,4,5,6,10)]
    seg.treatment.2$state <- round(2^seg.treatment.2$seg.mean * 2)
    seg.treatment.2$state[seg.treatment.2$state > 4] <- 4
    seg.treatment.2$method <- "sesame"
    row.names(seg.treatment.2) <- NULL
    seg.treatment.2 <- na.omit(seg.treatment.2)
    seg.treatment.2$chrom <-
      unlist(strsplit(seg.treatment.2$chrom,split="chr"))[c(FALSE,TRUE)]

    ##For sesame we want to omit the "chr"
    ##from chrom for intial plotting
    ##then add back in downstream for ggplot
    seg.out <- list(seg.treatment.1,
                    seg.treatment.2)

    names(seg.out) <- c(paste0(sesame.form.comparison[1],"_1"),
                        paste0(sesame.form.comparison[1],"_2"))

  }

} ##End formating

if(sesame.form.save.seg==TRUE){
  save(seg.out,
       file=paste0(sesame.form.output.dir,
                   .Platform$file.sep,
                   "seg.out.RData"))
}

return(seg.out)

}
