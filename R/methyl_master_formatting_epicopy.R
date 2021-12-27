#!/usr/bin/env Rscript

#' @title methyl_master_formatting_epicopy
#' @description Formatting Epicopy output into MethylMaster seg object
#' for downstream analyses
#' @param epi.form.seg
#' @param epi.form.output.dir
#' @param epi.form.save.seg
#' @param epi.form.comparison
#' @return Formatted epicopy result as MethylMaster seg object
#' for downstream analysis
#' @export
methyl_master_formatting_epicopy <- function(epi.form.seg=NULL,
                                             epi.form.output.dir=getwd(),
                                             epi.form.save.seg=FALSE,
                                             epi.form.comparison=NULL
                                             ){

  epicopy_results <- epi.form.seg[[1]]
  rm(epi.form.seg)

  epicopy_results$output$ID <-
    gsub("^X",
         "",
         epicopy_results$output$ID,perl = TRUE)

  seg <- epicopy_results$output
  rm(epicopy_results)

  ##colnames(seg) %>% cat(sep="\n")
  ##ID
  ##chrom
  ##loc.start
  ##loc.end
  ##num.mark
  ##seg.mean

  ##intersect(rownames(tumor), rownames(normal))
  ##intersect(tumor$X, normal$X)
  ##save(epicopy_results,
  ##     file = paste0(work.dir,
  ##                   file.sep,
  ##                   "epicopy_example_results.RData"))
  colnames(seg)[1] <- "Sample_ID"
  seg$seg.median <- NA
  seg$bstat      <- NA
  seg$pval       <- 0.05
  seg$state <- round(2^seg$seg.mean * 2)
  seg$state[seg$state > 4] <- 4
  seg$treatment  <- epi.form.comparison[1]
  seg$method <- "epicopy"
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

  names(seg.out) <- epi.form.comparison[1]

  if(epi.form.save.seg==TRUE){
    save(seg.out,
         file=paste0(epi.form.output.dir,
                     .Platform$file.sep,
                     "seg.out.RData"))
  }

  return(seg.out)

}
