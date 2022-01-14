#!/usr/bin/env Rscript

#' @title methyl_master_formatting_epicopy
#' @description Formatting the Epicopy output into MethylMaster seg object
#' for downstream analyses
#' @param epi.form.seg The input Epicopy seg results
#' @param epi.form.output.dir The output dir for the formatted seg object
#' @param epi.form.save.seg Whether to save the formatted Epicopy seg results
#' @param epi.form.comparison The comparison vector used earlier in the anlysis
#' @param epi.form.thresholds The thresholds used to determine the CNV
#' state, if NULL, the equation  seg.state <- round(2^seg.means * 2) is used
#' @return Formatted epicopy result as a MethylMaster seg object
#' for downstream analysis
#' @export
methyl_master_formatting_epicopy <- function(epi.form.seg=NULL,
                                             epi.form.output.dir=getwd(),
                                             epi.form.save.seg=FALSE,
                                             epi.form.comparison=NULL,
                                             epi.form.thresholds=NULL
                                             ){

    ##If not using GISTIC output:

    ##epicopy_results <- epi.form.seg[[1]]
    ##rm(epi.form.seg)

    ##epicopy_results$output$ID <-
    ##  gsub("^X",
    ##      "",
    ##      epicopy_results$output$ID,perl = TRUE)

    ##seg <- epicopy_results$output
    ##rm(epicopy_results)

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
    ##colnames(seg)[1] <- "Sample_ID"
    ##seg$seg.median <- NA
    ##seg$bstat      <- NA
    ##seg$pval       <- 0.05
    ##seg$state      <- calc_seg_state(seg.means = seg$seg.mean,
    ##                                 upper.thresh = 4,
    ##                                 cutoff = epi.form.thresholds
    ##)
    ##seg$treatment  <- epi.form.comparison[1]
    ##seg$method <- "epicopy"
    ##seg$sub.method <- NA
    ##row.names(seg) <- NULL

#############################################################################

  ##GISTIC formatted colnames:
  ##colnames(epicopy_results)
  ##"ID"
  ##"Chromosome"
  ##"Start"
  ##"End"
  ##"Probe_count"
  ##"Segment_Mean"

  epicopy_results <- epi.form.seg[[1]]
  rm(epi.form.seg)

  epicopy_results$ID <-
    gsub("^X",
    "",
    epicopy_results$ID,
    perl = TRUE)

  seg <- epicopy_results
  rm(epicopy_results)

  colnames(seg)[colnames(seg)=="ID"]           <- "Sample_ID"
  colnames(seg)[colnames(seg)=="Chromosome"]   <- "chrom"
  colnames(seg)[colnames(seg)=="Start"]        <- "loc.start"
  colnames(seg)[colnames(seg)=="End"]          <- "loc.end"
  colnames(seg)[colnames(seg)=="Probe_count"]  <- "num.mark"
  colnames(seg)[colnames(seg)=="Segment_Mean"] <- "seg.mean"

  seg$seg.median <- NA
  seg$bstat      <- NA
  seg$pval       <- 0.05
  seg$state      <- calc_seg_state(seg.means = seg$seg.mean,
                                   upper.thresh = 4,
                                   cutoff = epi.form.thresholds
  )
  seg$treatment  <- epi.form.comparison[1]
  seg$method <- "epicopy"
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

    names(seg.out) <- epi.form.comparison[1]

  if(epi.form.save.seg==TRUE){
    save(seg.out,
          file=paste0(epi.form.output.dir,
                      .Platform$file.sep,
                      "seg.out.RData"))
  }

  return(seg.out)

}
