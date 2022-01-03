#!/usr/bin/env Rscript

#' @title methyl_master_formatting_champ
#' @description formatting the champ seg results for
#' overlaps and visualization etc.a
#' Michael Mariani PhD Dartmouth College 2021
#' @param champ.form.seg
#' @param champ.form.output.dir
#' @param champ.form.reference
#' @param champ.form.split.by
#' @param champ.form.save.seg
#' @param champ.form.comparison
#' @param champ.form.padj
#' @param ...
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
if(champ.form.reference=="internal"){

  if(is.null(chmap.form.split.by))
  {

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
    seg$state <- round(2^seg$seg.mean * 2)
    seg$state[seg$state > 4] <- 4
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

  }else{

    champ.form.results.1 <- champ.form.seg[[1]]
    champ.form.results.2 <- champ.form.seg[[2]]
    rm(champ.form.seg)

    for(i in 1:length(champ.form.results.1$champ.CNA$sampleResult)){

      champ.form.results.1$champ.CNA$sampleResult[[i]]$ID <-
        names(champ.form.results.1$champ.CNA$sampleResult[i])

      colnames(champ.form.results.1$champ.CNA$sampleResult[[i]])[1] <-
        "Sample_ID"

    }

    for(i in 1:length(champ.form.results.2$champ.CNA$sampleResult)){

      champ.form.results$champ.CNA$sampleResult[[i]]$ID <-
        names(champ.form.results.2$champ.CNA$sampleResult[i])

      colnames(champ.form.results.2$champ.CNA$sampleResult[[i]])[1] <-
        "Sample_ID"

    }

    seg.1 <- do.call(rbind,
                   champ.form.results.1$champ.CNA$sampleResult)

    seg.2 <- do.call(rbind,
                   champ.form.results.2$champ.CNA$sampleResult)

    rm(champ.form.results.1)
    rm(champ.form.results.2)

    seg.1$chrom <- unlist(seg.1$chrom)
    seg.2$chrom <- unlist(seg.2$chrom)

    seg.1$seg.median <- NA
    seg.1$bstat      <- NA
    seg.1$pval       <- champ.form.padj
    seg.1$state <- round(2^seg.1$seg.mean * 2)
    seg.1$state[seg.1$state > 4] <- 4
    seg.1$treatment  <- champ.form.comparison[1]
    seg.1$method <- "champ"
    seg.1$sub.method <- NA
    row.names(seg.1) <- NULL

    seg.2$seg.median <- NA
    seg.2$bstat      <- NA
    seg.2$pval       <- champ.form.padj
    seg.2$state <- round(2^seg.2$seg.mean * 2)
    seg.2$state[seg.2$state > 4] <- 4
    seg.2$treatment  <- champ.form.comparison[1]
    seg.2$method <- "champ"
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

    seg.out <- list(seg.1, seg.2)

    names(seg.out) <- c(paste0(champ.form.comparison[1],"_1"),
                        paste0(champ.form.comparison[1],"_2"))


  }

}else{

  if(is.null(champ.form.split.by)){

    champ.form.results <- champ.form.seg[[1]]
    rm(champ.form.seg)

    for(i in 1:length(champ.form.results$champ.CNA$sampleResult)){

      champ.form.results$champ.CNA$sampleResult[[i]]$ID <-
        names(champ.form.results$champ.CNA$sampleResult[i])

      colnames(champ.form.results$champ.CNA$sampleResult[[i]])[1] <-
        "Sample_ID"

    }

    seg <- do.call(rbind,
                   champ.form.results$champ.CNA$sampleResult)

    rm(champ.form.results)

    seg$chrom <- unlist(seg$chrom)

    seg$seg.median <- NA
    seg$bstat      <- NA
    seg$pval       <- champ.form.padj
    seg$state <- round(2^seg$seg.mean * 2)
    seg$state[seg$state > 4] <- 4
    seg$treatment  <- champ.form.comparison[1]
    seg$method <- "champ"
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

    names(seg.out) <- champ.form.comparison[1]

  }else{

    champ.form.results.1 <- champ.form.seg[[1]]
    champ.form.results.2 <- champ.form.seg[[2]]
    rm(champ.form.seg)

    for(i in 1:length(champ.form.results.1$champ.CNA$sampleResult)){

      champ.form.results.1$champ.CNA$sampleResult[[i]]$ID <-
        names(champ.form.results.1$champ.CNA$sampleResult[i])

      colnames(champ.form.results.1$champ.CNA$sampleResult[[i]])[1] <-
        "Sample_ID"

    }

    for(i in 1:length(champ.form.results.2$champ.CNA$sampleResult)){

      champ.form.results$champ.CNA$sampleResult[[i]]$ID <-
        names(champ.form.results.2$champ.CNA$sampleResult[i])

      colnames(champ.form.results.2$champ.CNA$sampleResult[[i]])[1] <-
        "Sample_ID"

    }

    seg.1 <- do.call(rbind,
                     champ.form.results.1$champ.CNA$sampleResult)

    seg.2 <- do.call(rbind,
                     champ.form.results.2$champ.CNA$sampleResult)

    rm(champ.form.results.1)
    rm(champ.form.results.2)

    seg.1$chrom <- unlist(seg.1$chrom)
    seg.2$chrom <- unlist(seg.2$chrom)

    seg.1$seg.median <- NA
    seg.1$bstat      <- NA
    seg.1$pval       <- champ.form.padj
    seg.1$state <- round(2^seg.1$seg.mean * 2)
    seg.1$state[seg.1$state > 4] <- 4
    seg.1$treatment  <- champ.form.comparison[1]
    seg.1$method <- "champ"
    seg.1$sub.method <- NA
    row.names(seg.1) <- NULL

    seg.2$seg.median <- NA
    seg.2$bstat      <- NA
    seg.2$pval       <- champ.form.padj
    seg.2$state <- round(2^seg.2$seg.mean * 2)
    seg.2$state[seg.2$state > 4] <- 4
    seg.2$treatment  <- champ.form.comparison[1]
    seg.2$method <- "champ"
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

    seg.out <- list(seg.1, seg.2)

    names(seg.out) <- c(paste0(champ.form.comparison[1],"_1"),
                        paste0(champ.form.comparison[1],"_2"))

  }

}

  if(champ.form.save.seg==TRUE){
    save(seg.out,
         file=paste0(champ.form.output.dir,
                     .Platform$file.sep,
                     "seg.out.RData"))
  }

  return(seg.out)

}
