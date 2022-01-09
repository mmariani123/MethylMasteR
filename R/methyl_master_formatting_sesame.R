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
#' @param sesame.form.auto.corrected
#' @param sesame.form.thresholds
#' @param sesame.form.add.meta
#' @param ...
#' @import CNVRanger
#' @import matter
#' @importFrom magrittr %>%
#' @return #Formatted seg list object for visualizing
#' @export
methyl_master_formatting_sesame <- function(sesame.form.seg=NULL,
                                            sesame.form.output.dir=getwd(),
                                            sesame.form.sample.sheet.path=NULL,
                                            sesame.form.reference=NULL,
                                            sesame.form.split.by=NULL,
                                            sesame.form.comparison=NULL,
                                            sesame.form.save.seg=FALSE,
                                            sesame.form.plot.individual=FALSE,
                                            sesame.form.auto.corrected=FALSE,
                                            sesame.form.thresholds=c(-0.3,0.3),
                                            sesame.form.add.meta=NULL,
                                            ...){

##if(sesame.form.save.seg==TRUE){
##  save(sesame.form.seg, file =  paste0(sesame.form.output.dir,
##                                       .Platform$file.sep,
##                                       "seg.RData"))
##}

sesame.form.sample.sheet.df <- read.csv(file = sesame.form.sample.sheet.path,
                                        header = TRUE,
                                        stringsAsFactors = FALSE)

if(sesame.form.plot.individual==TRUE){

  if(!is.null(sesame.form.split.by)){

    split.cats <- unique(sesame.form.sample.sheet.df[[sesame.form.split.by]])

    for(i in 1:length(sesame.form.seg)){
      methyl_master_plot_individual(pi.seg = sesame.form.seg[[i]],
                                  pi.output.dir = sesame.form.output.dir,
                                  pi.name = paste0("tumor_",split.cats[i]))
    }

    }else{
      methyl_master_plot_individual(pi.seg = sesame.form.seg[[1]],
                                  pi.output.dir = sesame.form.output.dir,
                                  pi.name = "tumor")
  }

}

if(sesame.form.reference=="internal"){

  if(is.null(sesame.form.split.by)){

    sesame.seg <- sesame.form.seg[[1]]
    rm(sesame.form.seg)

    sesame.form.col <- sesame.form.sample.sheet.df[[sesame.form.add.meta]]
    names(sesame.form.col) <- sesame.form.add.col
    sesame_seg_treatment.df <- binding_frames_mm(
      x=sesame.seg,
      auto.corrected=sesame.form.auto.corrected,
      add.col=sesame.form.col
      )

    if(sesame.form.auto.corrected==FALSE){
      sesame_seg_treatment.df <- sesame_seg_treatment.df[
        sesame_seg_treatment.df$pval <= 0.05 &
          !is.na(sesame_seg_treatment.df$pval),]
    }

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

    if(sesame.form.auto.corrected==TRUE){
      ##colnames(seg)
      ##"Sample"
      ##"Chromosome"
      ##"bp.Start"
      ##"bp.End"
      ##"Num.of.Markers"
      ##"Mean"
      ##"ID"
      seg<-seg[,c(2:7)]
      colnames(seg)[6] <- "Sample_ID"
      colnames(seg)[1] <- "chrom"
      colnames(seg)[2] <- "loc.start"
      colnames(seg)[3] <- "loc.end"
      seg$chrom <-
        unlist(strsplit(seg$chrom,split="chr"))[c(FALSE,TRUE)]
      colnames(seg)[4] <- "num.mark"
      colnames(seg)[5] <- "seg.mean"
      seg$bstat        <- NA
      seg$seg.median   <- NA
      seg$pval         <- 0.05
      seg$state        <- calc_seg_state(seg.means = seg$seg.mean,
                                         upper.thresh = 4,
                                         use.cutoff = TRUE,
                                         cutoff = sesame.form.thresholds
                                        )
      seg$treatment    <- sesame.form.comparison[1]
      seg$method       <- "sesame"
      seg$sub.method   <- NA
      row.names(seg)   <- NULL

    }else{

      colnames(seg)[1] <- "Sample_ID"

      seg$chrom <-
        unlist(strsplit(seg$chrom,split="chr"))[c(FALSE,TRUE)]

      seg$bstat      <- NA
      seg$state      <- calc_seg_state(seg.means = seg$seg.mean,
                                       upper.thresh = 4,
                                       use.cutoff = TRUE,
                                       cutoff = sesame.form.thresholds,
                                       )
      seg$treatment  <- sesame.form.comparison[1]
      seg$method <- "sesame"
      seg$sub.method <- NA
      row.names(seg) <- NULL

    }

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
    ##sesame.form.sample.sheet.df[sesame.form.sample.sheet.df$Sample_Group %in%
    ##                         sesame.comparison[1], "Sample_Name"]]

    sesame_seg_treatment_1 <- sesame.form.seg[[1]]
    sesame_seg_treatment_2 <- sesame.form.seg[[2]]

    rm(sesame.form.seg)

    split.by.1 <- unique(sesame.form.sample.sheet.df[[sesame.form.split.by]])[1]
    split.by.2 <- unique(sesame.form.sample.sheet.df[[sesame.form.split.by]])[2]

    sesame_seg_treatment.df.1 <- binding_frames_mm(
      sesame_seg_treatment_1)

    sesame_seg_treatment.df.2 <- binding_frames_mm(
      sesame_seg_treatment_2)

    sesame_seg_treatment.df.1 <- sesame_seg_treatment.df.1[
      sesame_seg_treatment.df.1$pval <= 0.05 &
        !is.na(sesame_seg_treatment.df.1$pval),]

    sesame_seg_treatment.df.2 <- sesame_seg_treatment.df.2[
      sesame_seg_treatment.df.2$pval <= 0.05 &
        !is.na(sesame_seg_treatment.df.2$pval),]

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

    seg.1 <- sesame_seg_treatment.df.1
    rm(sesame_seg_treatment.df.1)

    seg.2 <- sesame_seg_treatment.df.2
    rm(sesame_seg_treatment.df.2)

    colnames(seg.1)[1] <- "Sample_ID"
    colnames(seg.2)[1] <- "Sample_ID"

    seg.1$chrom <-
      unlist(strsplit(seg.1$chrom,split="chr"))[c(FALSE,TRUE)]

    seg.2$chrom <-
      unlist(strsplit(seg.2$chrom,split="chr"))[c(FALSE,TRUE)]

    seg.1$bstat      <- NA
    seg.1$state <- round(2^seg.1$seg.mean * 2)
    seg.1$state[seg.1$state > 4] <- 4
    seg.1$treatment  <- paste0(sesame.form.comparison[1],
                               "_",
                               split.by.1)
    seg.1$method <- "sesame"
    seg.1$sub.method <- NA
    row.names(seg.1) <- NULL

    seg.2$bstat      <- NA
    seg.2$state <- round(2^seg.2$seg.mean * 2)
    seg.2$state[seg.2$state > 4] <- 4
    seg.2$treatment  <- paste0(sesame.form.comparison[1],
                               "_",
                               split.by.2)
    seg.2$method <- "sesame"
    seg.2$sub.method <- NA
    row.names(seg.2) <- NULL

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

    seg.1 <- seg.1[,preferred.column.names]
    seg.2 <- seg.2[,preferred.column.names]

    seg.out <- list(seg.1,seg.2)

    names(seg.out) <- c(paste0(sesame.form.comparison[1],"_",split.by.1),
                        paste0(sesame.form.comparison[1],"_",split.by.2))


  }

}else{

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
    ##sesame.form.sample.sheet.df[sesame.form.sample.sheet.df$Sample_Group %in%
    ##                         sesame.comparison[1], "Sample_Name"]]

    sesame_seg_treatment_1 <- sesame.form.seg[[1]]
    sesame_seg_treatment_2 <- sesame.form.seg[[2]]

    rm(sesame.form.seg)

    split.by.1 <- unique(sesame.form.sample.sheet.df[[sesame.form.split.by]])[1]
    split.by.2 <- unique(sesame.form.sample.sheet.df[[sesame.form.split.by]])[2]

    sesame_seg_treatment.df.1 <- binding_frames_mm(
      sesame_seg_treatment_1)

    sesame_seg_treatment.df.2 <- binding_frames_mm(
      sesame_seg_treatment_2)

    sesame_seg_treatment.df.1 <- sesame_seg_treatment.df.1[
      sesame_seg_treatment.df.1$pval <= 0.05 &
        !is.na(sesame_seg_treatment.df.1$pval),]

    sesame_seg_treatment.df.2 <- sesame_seg_treatment.df.2[
      sesame_seg_treatment.df.2$pval <= 0.05 &
        !is.na(sesame_seg_treatment.df.2$pval),]

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

    seg.1 <- sesame_seg_treatment.df.1
    rm(sesame_seg_treatment.df.1)

    seg.2 <- sesame_seg_treatment.df.2
    rm(sesame_seg_treatment.df.2)

    colnames(seg.1)[1] <- "Sample_ID"
    colnames(seg.2)[1] <- "Sample_ID"

    seg.1$chrom <-
      unlist(strsplit(seg.1$chrom,split="chr"))[c(FALSE,TRUE)]

    seg.2$chrom <-
      unlist(strsplit(seg.2$chrom,split="chr"))[c(FALSE,TRUE)]

    seg.1$bstat      <- NA
    seg.1$state <- round(2^seg.1$seg.mean * 2)
    seg.1$state[seg.1$state > 4] <- 4
    seg.1$treatment  <- paste0(sesame.form.comparison[1],
                               "_",
                               split.by.1)
    seg.1$method <- "sesame"
    seg.1$sub.method <- NA
    row.names(seg.1) <- NULL

    seg.2$bstat      <- NA
    seg.2$state <- round(2^seg.2$seg.mean * 2)
    seg.2$state[seg.2$state > 4] <- 4
    seg.2$treatment  <- paste0(sesame.form.comparison[1],
                               "_",
                               split.by.2)
    seg.2$method <- "sesame"
    seg.2$sub.method <- NA
    row.names(seg.2) <- NULL

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

    seg.1 <- seg.1[,preferred.column.names]
    seg.2 <- seg.2[,preferred.column.names]

    seg.out <- list(seg.1,seg.2)

    names(seg.out) <- c(paste0(sesame.form.comparison[1],"_",split.by.1),
                        paste0(sesame.form.comparison[1],"_",split.by.2))

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
