#!/usr/bin/env Rscript

#' @title methyl_master_formatting_custom
#' @description Formatting the results of the custom functionality to prepare
#' for comparison and output
#' @param custom.form.seg The input sesame segmentation list for formatting
#' @param custom.form.output.dir The output dir for sesame formatting
#' @param custom.form.sample.sheet.path The path to the MethylMaster
#' sample sheet
#' @param custom.form.comparison The MethylMaster comparison vector used in the
#' custom routine analysis
#' @param custom.form.save.seg Whether or not to save the formatted seg output
#' @param custom.form.thresholds The thresholds used to determine the CNV
#' state, if NULL, the equation  seg.state <- round(2^seg.means * 2) is used
#' @param ... additional parameters to pass to methyl_master_formatting_sesame
#' @import CNVRanger
#' @import matter
#' @importFrom magrittr %>%
#' @return Formatted seg list object for visualizing etc.
#' @export
methyl_master_formatting_custom <- function(custom.form.seg=NULL,
                                            custom.form.output.dir=getwd(),
                                            custom.form.sample.sheet.path=NULL,
                                            custom.form.comparison=NULL,
                                            custom.form.save.seg=FALSE,
                                            custom.form.thresholds=NULL,
                                            ...){

    ##if(sesame.form.save.seg==TRUE){
    ##  save(sesame.form.seg, file =  paste0(sesame.form.output.dir,
    ##                                       .Platform$file.sep,
    ##                                       "seg.RData"))
    ##}

    custom.form.sample.sheet.df <-
      read.csv(file = custom.form.sample.sheet.path,
               header = TRUE,
               stringsAsFactors = FALSE)

    custom.seg <- custom.form.seg[[1]]
    rm(custom.form.seg)

    ##custom.form.col <- custom.form.sample.sheet.df[[custom.form.add.meta]]
    ##names(custom.form.col) <-
    ##    rep(custom.form.add.meta, times=length(custom.form.col))

    custom_seg_treatment.df <- binding_frames_mm(
      x=custom.seg,
      auto.corrected=TRUE,
      add.col=NULL
      )

    seg <- custom_seg_treatment.df
    rm(custom_seg_treatment.df)

    colnames(seg)[7] <- "Sample_ID"
    colnames(seg)[2] <- "chrom"
    seg$chrom <-
      unlist(strsplit(seg$chrom,split="chr"))[c(FALSE,TRUE)]
    colnames(seg)[3] <- "loc.start"
    colnames(seg)[4] <- "loc.end"
    colnames(seg)[5] <- "num.mark"
    colnames(seg)[6] <- "seg.mean"
    seg$bstat        <- NA
    seg$seg.median   <- NA
    seg$pval         <- 0.05
    seg$state        <- calc_seg_state(seg.means = seg$seg.mean,
                                       upper.thresh = 4,
                                       cutoff = custom.form.thresholds
    )
    seg$treatment  <- custom.form.comparison[1]
    seg$method <- "custom"
    seg$sub.method <- NA
    row.names(seg) <- NULL

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

    names(seg.out) <- custom.form.comparison[1]

    if(custom.form.save.seg==TRUE){
      save(seg.out,
         file=paste0(custom.form.output.dir,
                   .Platform$file.sep,
                   "seg.out.RData"))
    }

    return(seg.out)

}
