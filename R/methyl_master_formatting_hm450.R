#!/usr/bin/env Rscript

#' @title methyl_master_formatting_hm450
#' @description formatting the hm450 output to prepare for plotting etc.
#' Michael Mariani PhD Dartmouth College 2021
#' @param seg
#' @param hm450.form.output.dir
#' @param hm450.form.sample.sheet.path
#' @param hm450.form.reference
#' @param hm450.form.split.by
#' @param hm450.form.workflow
#' @param hm450.form.comparison
#' @param hm450.form.save.seg
#' @param ...
#' @import CNVRanger
#' @import matter
#' @importFrom magrittr %>%
#' @return #Formatted seg list object for visualizing
#' @export
methyl_master_formatting_hm450 <- function(seg=NULL,
                                    hm450.form.output.dir=getwd(),
                                    hm450.form.sample.sheet.path=NULL,
                                    hm450.form.reference="internal",
                                    hm450.form.split.by=NULL,
                                    hm450.form.workflow="B",
                                    hm450.form.comparison=NULL,
                                    hm450.form.save.seg=FALSE,
                                    ...
                                  ){

if(hm450.form.reference=="internal"){

  if(is.null(hm450.form.split.by)){

  }else{

  }

}else{

  if(is.null(hm450.form.split.by)){

  }else{

    if(hm450.workflow=="A"){

      ##sub routine A

      candidates_data_treatment_A_sig.1 <-
        candidates_data_treatment_A.1[
        candidates_data_treatment_A.1$p.val <= 0.05,]

      candidates_data_treatment_A_sig.2   <-
        candidates_data_treatment_A.2[
        candidates_data_treatment_A.2$p.val <= 0.05,]

      candidates_data_treatment_A_sig.1$num.mark    <- NA
        candidates_data_treatment_A_sig.1$bstat     <- NA
        candidates_data_treatment_A_sig.1$treatment <- hm450.form.comparison[1]

      candidates_data_treatment_A_sig.2$num.mark    <- NA
        candidates_data_treatment_A_sig.2$bstat     <- NA
        candidates_data_treatment_A_sig.2$treatment <- hm450.form.comparison[1]

      preferred.columns <- c(7,1,2,3,12,13,8,5,4,9,10,11,14)

      candidates_data_treatment_A_sig.1 <-
        candidates_data_treatment_A_sig.1[, preferred.columns]

      candidates_data_treatment_A_sig.2 <-
        candidates_data_treatment_A_sig.2[, preferred.columns]

      candidates_data_treatment_A_sig_colnames <- c("ID",
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

      colnames(candidates_data_treatment_A_sig.1) <-
        candidates_data_treatment_A_sig_colnames

      colnames(candidates_data_treatment_A_sig.2) <-
        candidates_data_treatment_A_sig_colnames

      colnames(seg.1)[1] <- "Sample_ID"
      seg.1 <- seg.1[,c(1,2,3,4,5,8,10,11,12,13)]
      seg.1$state <- round(2^seg.1$seg.mean * 2)
      seg.1$state[seg.1$state > 4] <- 4
      seg.1$method <- "hm450"
      seg.1$sub.method <- hm450.form.sub.workflow
      row.names(seg.1) <- NULL
      seg.1 <- na.omit(seg.1) ##Workflow C ends up with some NA rows

      colnames(seg.2)[1] <- "Sample_ID"
      seg.2 <- seg.2[,c(1,2,3,4,5,8,10,11,12,13)]
      seg.2$state <- round(2^seg.2$seg.mean * 2)
      seg.2$state[seg.2$state > 4] <- 4
      seg.2$method <- "hm450"
      seg.2$sub.method <- hm450.form.sub.workflow
      row.names(seg.2) <- NULL
      seg.2 <- na.omit(seg.2) ##Workflow C ends up with some NA rows

      seg.out <- list(candidates_data_treatment_A_sig.1,
                    candidates_data_treatment_A_sig.2)

      names(seg.out) <- c(paste0(hm450.form.comparison[1],"_1"),
                          paste0(hm450.form.comparison[1],"_2"))

    }else if(k450.workflow==B){

      candidates_data_treatment_B_sig.1 <-
      candidates_data_treatment_B.1[
      candidates_data_treatment_B.1$p.val <= 0.05,]

      candidates_data_treatment_B_sig.2  <-
      candidates_data_treatment_B.2[
      candidates_data_treatment_B.2$p.val <= 0.05,]

      candidates_data_treatment_B_sig.1$num.mark  <- NA
      candidates_data_treatment_B_sig.1$bstat     <- NA
      candidates_data_treatment_B_sig.1$treatment <- hm450.form.comparison[1]

      candidates_data_treatment_B_sig.2$num.mark  <- NA
      candidates_data_treatment_B_sig.2$bstat     <- NA
      candidates_data_treatment_B_sig.2$treatment <- hm450.form.comparison[1]

      preferred.columns <- c(7,1,2,3,12,13,8,5,4,9,10,11,14)

      candidates_data_treatment_B_sig.1 <-
      candidates_data_treatment_B_sig.1[, preferred.columns]

      candidates_data_treatment_B_sig.2 <-
      candidates_data_treatment_B_sig.2[, preferred.columns]

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

      colnames(candidates_data_treatment_B_sig.1) <-
        candidates_data_normal_B_sig_colnames

      colnames(candidates_data_treatment_B_sig.2) <-
        candidates_data_normal_B_sig_colnames

      seg.1 <- candidates_data_treatment_B_sig.1
      seg.2 <- candidates_data_treatment_B_sig.2

      colnames(seg.1)[1] <- "Sample_ID"
      seg.1 <- seg.1[,c(1,2,3,4,5,8,10,11,12,13)]
      seg.1$state <- round(2^seg.1$seg.mean * 2)
      seg.1$state[seg.1$state > 4] <- 4
      seg.1$method <- "hm450"
      seg.1$sub.method <- hm450.form.sub.workflow
      row.names(seg.1) <- NULL
      seg.1 <- na.omit(seg.1) ##Workflow C ends up with some NA rows

      colnames(seg.2)[1] <- "Sample_ID"
      seg.2 <- seg.2[,c(1,2,3,4,5,8,10,11,12,13)]
      seg.2$state <- round(2^seg.2$seg.mean * 2)
      seg.2$state[seg.2$state > 4] <- 4
      seg.2$method <- "hm450"
      seg.2$sub.method <- hm450.form.sub.workflow
      row.names(seg.2) <- NULL
      seg.2 <- na.omit(seg.2) ##Workflow C ends up with some NA rows

      seg.out <- list(seg.1,
                      seg.2)

      names(seg) <- c(paste0(hm450.form.comparison[1],"_1"),
                      paste0(hm450.form.comparison[1],"_2"))

    }else if(k450.workflow=="C"){

      candidates_data_treatment_C_sig.1 <-
      candidates_data_treatment_C.1$data[
      candidates_data_treatment_C.1$data$pval <= 0.05,]

      candidates_data_treatment_C_sig.2    <-
      candidates_data_treatment_C.2$data[
      candidates_data_treatment_C.2$data$pval <= 0.05,]

      seg.1 <- candidates_data_treatment_C_sig.1
      seg.2 <- candidates_data_treatment_C_sig.2

      seg.1$treatment <- hm450.form.comparison[1]
      seg.2$treatment <- hm450.form.comparison[1]

      colnames(seg.1)[1] <- "Sample_ID"
      seg.1 <- seg.1[,c(1,2,3,4,5,8,10,11,12,13)]
      seg.1$state <- round(2^seg.1$seg.mean * 2)
      seg.1$state[seg.1$state > 4] <- 4
      seg.1$method <- "hm450"
      seg.1$sub.method <- hm450.form.sub.workflow
      row.names(seg.1) <- NULL
      seg.1 <- na.omit(seg.1) ##Workflow C ends up with some NA rows

      colnames(seg.2)[1] <- "Sample_ID"
      seg.2 <- seg.2[,c(1,2,3,4,5,8,10,11,12,13)]
      seg.2$state <- round(2^seg.2$seg.mean * 2)
      seg.2$state[seg.2$state > 4] <- 4
      seg.2$method <- "hm450"
      seg.2$sub.method <- hm450.form.sub.workflow
      row.names(seg.2) <- NULL
      seg.2 <- na.omit(seg.2) ##Workflow C ends up with some NA rows

      seg.out <- list(candidates_data_treatment_C_sig.1,
                    candidates_data_treatment_C_sig.2)

      names(seg.out) <- c(paste0(hm450.form.comparison[1],"_1"),
                          paste0(hm450.form.comparison[1],"_2"))

    }else{

      stop(paste0("Error: need to select a ",
                "proper sub workflow for 450k",
              " <hm450.workflow>"))

    }

  }

}

return(seg.out)

}
