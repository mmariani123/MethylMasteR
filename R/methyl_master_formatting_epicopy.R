#!/usr/bin/env Rscript

#' @title methyl_master_formatting_epicopy
#' @description Formatting Epicopy output into MethylMaster seg object
#' for downstream analyses
#' @param epi.form.seg
#' @param epi.form.output.dir
#' @param epi.form.reference
#' @param epi.form.split.by
#' @param epi.form.save.seg
#' @param epi.form.comparison
#' @return Formatted epicopy result as MethylMaster seg object
#' for downstream analysis
#' @export
methyl_master_formatting_epicopy <- function(epi.form.seg=NULL,
                                             epi.form.output.dir=getwd(),
                                             epi.form.reference=NULL,
                                             epi.form.split.by=NULL,
                                             epi.form.save.seg=FALSE,
                                             epi.form.comparison=NULL
                                             ){

if(epi.form.reference=="internal"){

  if(is.null(epi.split.by)){

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

  }else{

    epicopy_results.1 <- epi.form.seg[[1]]
    epicopy_results.2 <- epi.form.seg[[2]]
    rm(epi.form.seg)

    epicopy_results.1$output$ID <-
      gsub("^X",
           "",
           epicopy_results.1$output$ID,perl = TRUE)

    epicopy_results.2$output$ID <-
      gsub("^X",
           "",
           epicopy_results.2$output$ID,perl = TRUE)

    seg.1 <- epicopy_results.1$output
    seg.2 <- epicopy_results.2$output
    rm(epicopy_results.1)
    rm(epicopy_results.2)

    colnames(seg.1)[1] <- "Sample_ID"
    seg.1$seg.median <- NA
    seg.1$bstat      <- NA
    seg.1$pval       <- 0.05
    seg.1$state <- round(2^seg.1$seg.mean * 2)
    seg.1$state[seg.1$state > 4] <- 4
    seg.1$treatment  <- epi.form.comparison[1]
    seg.1$method <- "epicopy"
    seg.1$sub.method <- NA
    row.names(seg.1) <- NULL

    colnames(seg.2)[1] <- "Sample_ID"
    seg.2$seg.median <- NA
    seg.2$bstat      <- NA
    seg.2$pval       <- 0.05
    seg.2$state <- round(2^seg.2$seg.mean * 2)
    seg.2$state[seg.2$state > 4] <- 4
    seg.2$treatment <- epi.form.comparison[1]
    seg.2$method <- "epicopy"
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

    names(seg.out) <- c(paste0(epi.form.comparison[1],"_1",
                               epi.form.comparison[1],"_2"))
  }

}else{

  if(is.null(epi.form.split.by)){

    epicopy_results <- epi.form.seg[[1]]
    rm(epi.form.seg)

    epicopy_results$output$ID <-
      gsub("^X",
           "",
           epicopy_results$output$ID,perl = TRUE)

    seg <- epicopy_results$output
    rm(epicopy_results)

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

  }else{

    epicopy_results.1 <- epi.form.seg[[1]]
    epicopy_results.2 <- epi.form.seg[[2]]
    rm(epi.form.seg)

    epicopy_results.1$output$ID <-
      gsub("^X",
           "",
           epicopy_results.1$output$ID,perl = TRUE)

    epicopy_results.2$output$ID <-
      gsub("^X",
           "",
           epicopy_results.2$output$ID,perl = TRUE)

    seg.1 <- epicopy_results.1$output
    seg.2 <- epicopy_results.2$output
    rm(epicopy_results.1)
    rm(epicopy_results.2)

    colnames(seg.1)[1] <- "Sample_ID"
    seg.1$seg.median <- NA
    seg.1$bstat      <- NA
    seg.1$pval       <- 0.05
    seg.1$state <- round(2^seg.1$seg.mean * 2)
    seg.1$state[seg.1$state > 4] <- 4
    seg.1$treatment  <- epi.form.comparison[1]
    seg.1$method <- "epicopy"
    seg.1$sub.method <- NA
    row.names(seg.1) <- NULL

    colnames(seg.2)[1] <- "Sample_ID"
    seg.2$seg.median <- NA
    seg.2$bstat      <- NA
    seg.2$pval       <- 0.05
    seg.2$state <- round(2^seg.2$seg.mean * 2)
    seg.2$state[seg.2$state > 4] <- 4
    seg.2$treatment <- epi.form.comparison[1]
    seg.2$method <- "epicopy"
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

    names(seg.out) <- c(paste0(epi.form.comparison[1],"_1",
                               epi.form.comparison[1],"_2"))

  }

}

if(epi.form.save.seg==TRUE){
   save(seg.out,
        file=paste0(epi.form.output.dir,
                    .Platform$file.sep,
                    "seg.out.RData"))
}

return(seg.out)

}
