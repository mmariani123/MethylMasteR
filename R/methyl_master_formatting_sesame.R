#!/usr/bin/env Rscript

#' @title methyl_master_formatting_sesame
#' @description formatting the results of the sesame functionality to prepare
#' for comparison and output
#' @param sesame.form.seg The input sesame segmentation list for formatting
#' @param sesame.form.output.dir The output dir for sesame formatting
#' @param sample.form.sample.sheet.path The path to the MethylMaster
#' sample sheet
#' @param sesame.form.reference The sesame reference used in the analysis
#' @param sesame.form.split.by The split.by field used in the sesame routine
#' in the MethylMaster sample sheet
#' @param sesame.form.comparison The MethylMaster comparison vector used in the
#' Sesame routine analysis
#' @param sesame.form.save.seg Whether or not to save the formatted seg output
#' @param sesame.form.plot.individual Whether to plot the individual sesame
#' sample signal plots
#' @param sesame.form.thresholds The sesame thresholds used to determine the CNV
#' state, if NULL, the equation  seg.state <- round(2^seg.means * 2) is used
#' @param ... additional parameters to pass to methyl_master_formatting_sesame
#' @import CNVRanger
#' @import matter
#' @importFrom magrittr %>%
#' @return Formatted seg list object for visualizing etc.
#' @export
methyl_master_formatting_sesame <- function(sesame.form.seg=NULL,
                                            sesame.form.output.dir=getwd(),
                                            sesame.form.sample.sheet.path=NULL,
                                            sesame.form.reference=NULL,
                                            sesame.form.split.by=NULL,
                                            sesame.form.comparison=NULL,
                                            sesame.form.save.seg=FALSE,
                                            sesame.form.plot.individual=FALSE,
                                            sesame.form.thresholds=NULL,
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

methyl_master_plot_individual(pi.seg = sesame.form.seg[[1]],
                              pi.output.dir = sesame.form.output.dir,
                              pi.name = "tumor")
}

sesame.seg <- sesame.form.seg[[1]]
rm(sesame.form.seg)

sesame.form.col <- sesame.form.sample.sheet.df[[sesame.form.add.meta]]

names(sesame.form.col) <-
  rep(sesame.form.add.meta,times=length(sesame.form.col))

sesame_seg_treatment.df <- binding_frames_mm(
  x=sesame.seg,
  auto.corrected=sesame.form.auto.corrected,
  add.col=sesame.form.col
)

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
seg$state      <- calc_seg_state(seg.means = seg$seg.mean,
                                 upper.thresh = 4,
                                 cutoff = sesame.form.thresholds
)
seg$treatment  <- sesame.form.comparison[1]
seg$method <- "sesame"
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
                            "sub.method",
                            sesame.form.add.meta)

seg <- seg[,preferred.column.names]

seg.out <- list(seg)

names(seg.out) <- sesame.form.comparison[1]

if(sesame.form.save.seg==TRUE){
  save(seg.out,
       file=paste0(sesame.form.output.dir,
                   .Platform$file.sep,
                   "seg.out.RData"))
}

return(seg.out)

}
