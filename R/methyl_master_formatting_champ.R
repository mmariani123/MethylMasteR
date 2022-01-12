#!/usr/bin/env Rscript

#' @title methyl_master_formatting_champ
#' @description formatting the champ seg results for
#' overlaps and visualization etc.a
#' Michael Mariani PhD Dartmouth College 2021
#' @param champ.form.seg The input ChAMP CNV segmentation results
#' @param champ.form.output.dir The output directory
#' @param champ.form.reference The reference that was used
#' @param champ.form.split.by The splt.by variable in the sample sheet that was
#' used
#' @param champ.form.save.seg Whether to save the formatted ChAMP results
#' @param champ.form.comparison The MethylMaster comparison vector
#' @param champ.form.padj The padj that was used in the ChAMP routine
#' @param ... Additional parameters passed to methyl_master_champ_formatting
#' @import CNVRanger
#' @import matter
#' @importFrom magrittr %>%
#' @return Formatted champ result as MethylMaster seg object
#' for downstream analysis
#' @export
methyl_master_formatting_champ <- function(champ.form.seg=NULL,
                                           champ.form.output.dir=getwd(),
                                           champ.form.reference=NULL,
                                           champ.form.split.by=NULL,
                                           champ.form.save.seg=FALSE,
                                           champ.form.comparison=NULL,
                                           champ.form.padj=NULL
                                           ){

    ##For testing:
    ##load(file=paste0("G:\\My Drive\\dartmouth",
    ##                 "\\salas_lab_working\\cnv",
    ##                 "\\cnv_testing\\mm_testing_old",
    ##                 "\\CHAMP_RESULT\\myCNA.rda"))
    ##myCNA$sampleResult
    ##champ.form.results <- list()
    ##champ.form.results$champ.CNA <- myCNA

    champ.form.results <- champ.form.seg[[1]]
    rm(champ.form.seg)

    for(i in 1:length(champ.form.results$champ.CNA$sampleResult)){

      champ.form.results$champ.CNA$sampleResult[[i]]$ID <-
        names(champ.form.results$champ.CNA$sampleResult[i])

      colnames(champ.form.results$champ.CNA$sampleResult[[i]])[1] <-
        "Sample_ID"
      ##colnames(myCNA$sampleResult[[i]])[1] ##All set

    }

    ##champ.betas <- lapply(champ.form.results[[1]],
    ##                      FUN = function(x){return(x$beta)})

    ##champ.intensities <- lapply(champ.form.results[[1]],
    ##                      FUN = function(x){return(x$intensity)})

    seg <- do.call(rbind,
                 champ.form.results$champ.CNA$sampleResult)

    rm(champ.form.results)

    ##load(file=paste0(output.dir,
    ##                 .Platform$file.sep,
    ##                 "champ_seg.RData"))
    ##seg <- champ_seg
    ##rm(champ_seg)

    ##seg$chrom <-
    ##  as.integer(unlist(strsplit(seg$chrom,split="chr"))[c(FALSE,TRUE)])
    ##seg$chrom <- as.integer(unlist(seg$chrom))
    ##Don't do as.integer() becasue may contain X,Y etc.
    seg$chrom <- unlist(seg$chrom)

    ##"Sample_ID"
    ##"chrom"
    ##"loc.start"
    ##"loc.end"   "
    ##"num.mark"
    ##"seg.mean"

    seg$seg.median <- NA
    seg$bstat      <- NA
    seg$pval       <- champ.form.padj
    seg$state      <- calc_seg_state(seg.means = seg$seg.mean,
                                     upper.thresh = 4,
                                     cutoff = champ.form.thresholds
                      )
    seg$treatment  <- champ.form.comparison[1]
    seg$method <- "champ"
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

    names(seg.out) <- champ.form.comparison[1]

  if(champ.form.save.seg==TRUE){
    save(seg.out,
         file=paste0(champ.form.output.dir,
                     .Platform$file.sep,
                     "seg.out.RData"))
  }

  return(seg.out)

}
