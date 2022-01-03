#!/usr/bin/env Rscript

#' @title methyl_master_formatting_hm450
#' @description formatting the hm450 output to prepare for plotting etc.
#' Michael Mariani PhD Dartmouth College 2021
#' @param hm450.form.seg
#' @param hm450.form.output.dir
#' @param hm450.form.sample.sheet.path
#' @param hm450.form.reference
#' @param hm450.form.split.by
#' @param hm450.form.workflow
#' @param hm450.form.comparison
#' @param hm450.form.save.seg
#' @param hm450.form.anno.file.path
#' @param ...
#' @import CNVRanger
#' @import matter
#' @importFrom magrittr %>%
#' @return #Formatted seg list object for visualizing
#' @export
methyl_master_formatting_hm450 <- function(hm450.form.seg=NULL,
                                    hm450.form.output.dir=getwd(),
                                    hm450.form.sample.sheet.path=NULL,
                                    hm450.form.reference="internal",
                                    hm450.form.split.by=NULL,
                                    hm450.form.workflow="B",
                                    hm450.form.comparison=NULL,
                                    hm450.form.save.seg=FALSE,
                                    hm450.form.anno.file.path,
                                    ...
                                  ){

##load(hm450.anno.file.path)
##load("G:\\My Drive\\dartmouth\\salas_lab_working\\cnv
##\\cnv_testing\\probe450kfemanno.rda")
##load("G:\\My Drive\\dartmouth\\salas_lab_working\\cnv
##\\cnv_testing\\hm450.manifest.hg38.rda")
##https://www.bioconductor.org/packages/release/BiocViews.html#___IlluminaChip
##annotation1 <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
##annotation1 <- as.data.frame(annotation1)

if(hm450.form.reference=="internal"){

  if(is.null(hm450.form.split.by)){

    if(hm450.form.workflow=="A"){

      ##sub routine A

      candidates_data_treatment_A_sig <-
        candidates_data_treatment_A[candidates_data_treatment_B$p.val <= 0.05,]

      rm(candidates_data_treatment_A)

      colnames(candidates_data_treatment_A_sig)[
        colnames(candidates_data_treatment_A_sig)=="smp"] <- "Sample_ID"

      colnames(candidates_data_treatment_A_sig)[
        colnames(candidates_data_treatment_A_sig)=="chr"] <- "chrom"

      candidates_data_treatment_A_sig$chrom <-
        unlist(strsplit(candidates_data_treatment_A_sig$chrom,
                        split="chr"))[c(FALSE,TRUE)]

      colnames(candidates_data_treatment_A_sig)[
        colnames(candidates_data_treatment_A_sig)=="startPos"] <- "loc.start"

      colnames(candidates_data_treatment_A_sig)[
        colnames(candidates_data_treatment_A_sig)=="endPos"] <- "loc.end"

      colnames(candidates_data_treatment_A_sig)[
        colnames(candidates_data_treatment_A_sig)=="mean"] <- "seg.mean"

      colnames(candidates_data_treatment_A_sig)[
        colnames(candidates_data_treatment_A_sig)=="median"] <- "seg.median"

      colnames(candidates_data_treatment_A_sig)[
        colnames(candidates_data_treatment_A_sig)=="p.val"] <- "pval"

      candidates_data_treatment_A_sig$num.mark  <- NA
      candidates_data_treatment_A_sig$bstat     <- NA
      candidates_data_treatment_A_sig$state <-
        round(2^candidates_data_treatment_A_sig$seg.mean * 2)
      candidates_data_treatment_A_sig$state[
        candidates_data_treatment_A_sig$state > 4] <- 4
      candidates_data_treatment_A_sig$treatment <- hm450.form.comparison[1]
      candidates_data_treatment_A_sig$method <- "hm450"
      candidates_data_treatment_A_sig$sub.method <- "A"
      row.names(candidates_data_treatment_B_sig) <- NULL
      ##seg <- na.omit(seg) ##Workflow C ends up with some NA rows

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

      candidates_data_treatment_A_sig <-
        candidates_data_treatment_A_sig[,preferred.column.names]

      seg.out <- list(candidates_data_treatment_A_sig)

      names(seg.out) <- hm450.form.comparison[1]

    }else if(hm450.form.workflow=="B"){

      candidates_data_treatment_B <- hm450.form.seg[[1]]
      rm(hm450.form.seg)

      ##colnames(candidates_data_treatment_B)
      ##"chr"
      ##"startPos"
      ##"endPos"
      ##"median"
      ##"mean"
      ##"sd"
      ##"smp"
      ##"p.val"

      ##any(is.na=(candidates_data_treatment_B$chr))
      ##No NA for chr field which is good

      candidates_data_treatment_B_sig <-
        candidates_data_treatment_B[candidates_data_treatment_B$p.val <= 0.05,]

      rm(candidates_data_treatment_B)

      colnames(candidates_data_treatment_B_sig)[
        colnames(candidates_data_treatment_B_sig)=="smp"] <- "Sample_ID"

      colnames(candidates_data_treatment_B_sig)[
        colnames(candidates_data_treatment_B_sig)=="chr"] <- "chrom"

      candidates_data_treatment_B_sig$chrom <-
        unlist(strsplit(candidates_data_treatment_B_sig$chrom,
                        split="chr"))[c(FALSE,TRUE)]

      colnames(candidates_data_treatment_B_sig)[
        colnames(candidates_data_treatment_B_sig)=="startPos"] <- "loc.start"

      colnames(candidates_data_treatment_B_sig)[
        colnames(candidates_data_treatment_B_sig)=="endPos"] <- "loc.end"

      colnames(candidates_data_treatment_B_sig)[
        colnames(candidates_data_treatment_B_sig)=="mean"] <- "seg.mean"

      colnames(candidates_data_treatment_B_sig)[
        colnames(candidates_data_treatment_B_sig)=="median"] <- "seg.median"

      colnames(candidates_data_treatment_B_sig)[
        colnames(candidates_data_treatment_B_sig)=="p.val"] <- "pval"

      candidates_data_treatment_B_sig$num.mark  <- NA
      candidates_data_treatment_B_sig$bstat     <- NA
      candidates_data_treatment_B_sig$state <-
        round(2^candidates_data_treatment_B_sig$seg.mean * 2)
      candidates_data_treatment_B_sig$state[
        candidates_data_treatment_B_sig$state > 4] <- 4
      candidates_data_treatment_B_sig$treatment <- hm450.form.comparison[1]
      candidates_data_treatment_B_sig$method <- "hm450"
      candidates_data_treatment_B_sig$sub.method <- "B"
      row.names(candidates_data_treatment_B_sig) <- NULL
      ##seg <- na.omit(seg) ##Workflow C ends up with some NA rows

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

      candidates_data_treatment_B_sig <-
        candidates_data_treatment_B_sig[,preferred.column.names]

      seg.out <- list(candidates_data_treatment_B_sig)

      names(seg.out) <- hm450.form.comparison[1]

      ##annotation1$probe <- rownames(annotation1)
      ##seg.1$loc.start
      ##hm450.manifest.hg38$addressA
      ##hm450.manifest.hg38$probeStart
      ##hm450.manifest.hg38$probeEnd

    }else if(hm450.form.workflow=="C"){

      candidates_data_treatment_C <- hm450.form.seg[[1]]
      rm(hm450.form.seg)

      ##colnames(candidates_data_treatment_B)
      ##"chr"
      ##"startPos"
      ##"endPos"
      ##"median"
      ##"mean"
      ##"sd"
      ##"smp"
      ##"p.val"

      ##any(is.na=(candidates_data_treatment_B$chr))
      ##No NA for chr field which is good

      candidates_data_treatment_C_sig <-
        candidates_data_treatment_C[candidates_data_treatment_C$p.val <= 0.05,]

      rm(candidates_data_treatment_C)

      colnames(candidates_data_treatment_C_sig)[
        colnames(candidates_data_treatment_C_sig)=="smp"] <- "Sample_ID"

      colnames(candidates_data_treatment_C_sig)[
        colnames(candidates_data_treatment_C_sig)=="chr"] <- "chrom"

      candidates_data_treatment_C_sig$chrom <-
        unlist(strsplit(candidates_data_treatment_C_sig$chrom,
                        split="chr"))[c(FALSE,TRUE)]

      colnames(candidates_data_treatment_C_sig)[
        colnames(candidates_data_treatment_C_sig)=="startPos"] <- "loc.start"

      colnames(candidates_data_treatment_C_sig)[
        colnames(candidates_data_treatment_C_sig)=="endPos"] <- "loc.end"

      colnames(candidates_data_treatment_C_sig)[
        colnames(candidates_data_treatment_C_sig)=="mean"] <- "seg.mean"

      colnames(candidates_data_treatment_C_sig)[
        colnames(candidates_data_treatment_C_sig)=="median"] <- "seg.median"

      colnames(candidates_data_treatment_C_sig)[
        colnames(candidates_data_treatment_C_sig)=="p.val"] <- "pval"

      candidates_data_treatment_C_sig$num.mark  <- NA
      candidates_data_treatment_C_sig$bstat     <- NA
      candidates_data_treatment_C_sig$state <-
        round(2^candidates_data_treatment_C_sig$seg.mean * 2)
      candidates_data_treatment_C_sig$state[
        candidates_data_treatment_C_sig$state > 4] <- 4
      candidates_data_treatment_C_sig$treatment <- hm450.form.comparison[1]
      candidates_data_treatment_C_sig$method <- "hm450"
      candidates_data_treatment_C_sig$sub.method <- "C"
      row.names(candidates_data_treatment_C_sig) <- NULL

      candidates_data_treatment_C_sig <-
      na.omit(candidates_data_treatment_C_sig)
      ##Workflow C ends up with some NA rows

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

      candidates_data_treatment_C_sig <-
        candidates_data_treatment_C_sig[,preferred.column.names]

      seg.out <- list(candidates_data_treatment_C_sig)

      names(seg.out) <- hm450.form.comparison[1]

    }else{

      stop(paste0("Error: need to select a ",
                  "proper sub workflow for 450k",
                  " <hm450.form.workflow>"))

    }

  }else{

    if(hm450.form.workflow=="A"){

      ##sub routine A
      candidates_data_treatment_A.1 <- hm450.form.seg[[1]]
      candidates_data_treatment_A.2 <- hm450.form.seg[[2]]
      rm(hm450.form.seg)

      candidates_data_treatment_A_sig.1 <-
        candidates_data_treatment_A_sig.1[
          candidates_data_treatment_A.1$p.val <= 0.05,]
      candidates_data_treatment_A_sig.2 <-
        candidates_data_treatment_A_sig.2[
          candidates_data_treatment_A.2$p.val <= 0.05,]

      rm(candidates_data_treatment_A.1)
      rm(candidates_data_treatment_A.2)

      colnames(candidates_data_treatment_A_sig.1)[
        colnames(candidates_data_treatment_A_sig.1)=="smp"] <- "Sample_ID"
      colnames(candidates_data_treatment_A_sig.2)[
        colnames(candidates_data_treatment_A_sig.2)=="smp"] <- "Sample_ID"

      colnames(candidates_data_treatment_A_sig.1)[
        colnames(candidates_data_treatment_A_sig.1)=="chr"] <- "chrom"
      colnames(candidates_data_treatment_A_sig.2)[
        colnames(candidates_data_treatment_A_sig.2)=="chr"] <- "chrom"

      candidates_data_treatment_A_sig.1$chrom <-
        unlist(strsplit(candidates_data_treatment_A_sig.1$chrom,
                        split="chr"))[c(FALSE,TRUE)]
      candidates_data_treatment_A_sig.2$chrom <-
        unlist(strsplit(candidates_data_treatment_A_sig.2$chrom,
                        split="chr"))[c(FALSE,TRUE)]

      colnames(candidates_data_treatment_A_sig.1)[
        colnames(candidates_data_treatment_A_sig.1)=="startPos"] <- "loc.start"
      colnames(candidates_data_treatment_A_sig.2)[
        colnames(candidates_data_treatment_A_sig.2)=="startPos"] <- "loc.start"

      colnames(candidates_data_treatment_A_sig.1)[
        colnames(candidates_data_treatment_A_sig.1)=="endPos"] <- "loc.end"
      colnames(candidates_data_treatment_A_sig.2)[
        colnames(candidates_data_treatment_A_sig.2)=="endPos"] <- "loc.end"

      colnames(candidates_data_treatment_A_sig.1)[
        colnames(candidates_data_treatment_A_sig.1)=="mean"] <- "seg.mean"
      colnames(candidates_data_treatment_A_sig.2)[
        colnames(candidates_data_treatment_A_sig.2)=="mean"] <- "seg.mean"

      colnames(candidates_data_treatment_A_sig.1)[
        colnames(candidates_data_treatment_A_sig.1)=="median"] <- "seg.median"
      colnames(candidates_data_treatment_A_sig.2)[
        colnames(candidates_data_treatment_A_sig.2)=="median"] <- "seg.median"

      colnames(candidates_data_treatment_A_sig.1)[
        colnames(candidates_data_treatment_A_sig.1)=="p.val"] <- "pval"
      colnames(candidates_data_treatment_A_sig.2)[
        colnames(candidates_data_treatment_A_sig.2)=="p.val"] <- "pval"

      candidates_data_treatment_A_sig.1$num.mark  <- NA
      candidates_data_treatment_A_sig.1$bstat     <- NA
      candidates_data_treatment_A_sig.1$state <-
        round(2^candidates_data_treatment_A_sig.1$seg.mean * 2)
      candidates_data_treatment_A_sig.1$state[
        candidates_data_treatment_A_sig.1$state > 4] <- 4
      candidates_data_treatment_A_sig.1$treatment <- hm450.form.comparison[1]
      candidates_data_treatment_A_sig.1$method <- "hm450"
      candidates_data_treatment_A_sig.1$sub.method <- "A"
      row.names(candidates_data_treatment_A_sig.1) <- NULL
      ##candidates_data_treatment_C_sig <-
      ##  na.omit(candidates_data_treatment_C_sig)
      ##Workflow C ends up with some NA rows

      candidates_data_treatment_A_sig.2$num.mark  <- NA
      candidates_data_treatment_A_sig.2$bstat     <- NA
      candidates_data_treatment_A_sig.2$state <-
        round(2^candidates_data_treatment_A_sig.2$seg.mean * 2)
      candidates_data_treatment_A_sig.2$state[
        candidates_data_treatment_A_sig.2$state > 4] <- 4
      candidates_data_treatment_A_sig.2$treatment <- hm450.form.comparison[1]
      candidates_data_treatment_A_sig.2$method <- "hm450"
      candidates_data_treatment_A_sig.2$sub.method <- "A"
      row.names(candidates_data_treatment_A_sig.2) <- NULL
      ##candidates_data_treatment_C_sig <-
      ##  na.omit(candidates_data_treatment_C_sig)
      ##Workflow C ends up with some NA rows


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

      candidates_data_treatment_A_sig.1 <-
        candidates_data_treatment_A_sig.1[,preferred.column.names]
      candidates_data_treatment_A_sig.2 <-
        candidates_data_treatment_A_sig.2[,preferred.column.names]

      seg.out <- list(candidates_data_treatment_A_sig.1,
                      candidates_data_treatment_A_sig.2)

      names(seg.out) <- c(paste0(hm450.form.comparison[1],"_1"),
                          paste0(hm450.form.comparison[1],"_2"))

    }else if(hm450.form.workflow=="B"){

      ##sub routine B
      candidates_data_treatment_B.1 <- hm450.form.seg[[1]]
      candidates_data_treatment_B.2 <- hm450.form.seg[[2]]
      rm(hm450.form.seg)

      candidates_data_treatment_B_sig.1 <-
        candidates_data_treatment_B_sig.1[
          candidates_data_treatment_B.1$p.val <= 0.05,]
      candidates_data_treatment_B_sig.2 <-
        candidates_data_treatment_B_sig.2[
          candidates_data_treatment_B.2$p.val <= 0.05,]

      rm(candidates_data_treatment_B.1)
      rm(candidates_data_treatment_B.2)

      colnames(candidates_data_treatment_B_sig.1)[
        colnames(candidates_data_treatment_B_sig.1)=="smp"] <- "Sample_ID"
      colnames(candidates_data_treatment_B_sig.2)[
        colnames(candidates_data_treatment_B_sig.2)=="smp"] <- "Sample_ID"

      colnames(candidates_data_treatment_B_sig.1)[
        colnames(candidates_data_treatment_B_sig.1)=="chr"] <- "chrom"
      colnames(candidates_data_treatment_B_sig.2)[
        colnames(candidates_data_treatment_B_sig.2)=="chr"] <- "chrom"

      candidates_data_treatment_B_sig.1$chrom <-
        unlist(strsplit(candidates_data_treatment_B_sig.1$chrom,
                        split="chr"))[c(FBLSE,TRUE)]
      candidates_data_treatment_B_sig.2$chrom <-
        unlist(strsplit(candidates_data_treatment_B_sig.2$chrom,
                        split="chr"))[c(FBLSE,TRUE)]

      colnames(candidates_data_treatment_B_sig.1)[
        colnames(candidates_data_treatment_B_sig.1)=="startPos"] <- "loc.start"
      colnames(candidates_data_treatment_B_sig.2)[
        colnames(candidates_data_treatment_B_sig.2)=="startPos"] <- "loc.start"

      colnames(candidates_data_treatment_B_sig.1)[
        colnames(candidates_data_treatment_B_sig.1)=="endPos"] <- "loc.end"
      colnames(candidates_data_treatment_B_sig.2)[
        colnames(candidates_data_treatment_B_sig.2)=="endPos"] <- "loc.end"

      colnames(candidates_data_treatment_B_sig.1)[
        colnames(candidates_data_treatment_B_sig.1)=="mean"] <- "seg.mean"
      colnames(candidates_data_treatment_B_sig.2)[
        colnames(candidates_data_treatment_B_sig.2)=="mean"] <- "seg.mean"

      colnames(candidates_data_treatment_B_sig.1)[
        colnames(candidates_data_treatment_B_sig.1)=="median"] <- "seg.median"
      colnames(candidates_data_treatment_B_sig.2)[
        colnames(candidates_data_treatment_B_sig.2)=="median"] <- "seg.median"

      colnames(candidates_data_treatment_B_sig.1)[
        colnames(candidates_data_treatment_B_sig.1)=="p.val"] <- "pval"
      colnames(candidates_data_treatment_B_sig.2)[
        colnames(candidates_data_treatment_B_sig.2)=="p.val"] <- "pval"

      candidates_data_treatment_B_sig.1$num.mark  <- NB
      candidates_data_treatment_B_sig.1$bstat     <- NB
      candidates_data_treatment_B_sig.1$state <-
        round(2^candidates_data_treatment_B_sig.1$seg.mean * 2)
      candidates_data_treatment_B_sig.1$state[
        candidates_data_treatment_B_sig.1$state > 4] <- 4
      candidates_data_treatment_B_sig.1$treatment <- hm450.form.comparison[1]
      candidates_data_treatment_B_sig.1$method <- "hm450"
      candidates_data_treatment_B_sig.1$sub.method <- "B"
      row.names(candidates_data_treatment_B_sig.1) <- NULL
      ##candidates_data_treatment_C_sig <-
      ##  na.omit(candidates_data_treatment_C_sig)
      ##Workflow C ends up with some NB rows

      candidates_data_treatment_B_sig.2$num.mark  <- NB
      candidates_data_treatment_B_sig.2$bstat     <- NB
      candidates_data_treatment_B_sig.2$state <-
        round(2^candidates_data_treatment_B_sig.2$seg.mean * 2)
      candidates_data_treatment_B_sig.2$state[
        candidates_data_treatment_B_sig.2$state > 4] <- 4
      candidates_data_treatment_B_sig.2$treatment <- hm450.form.comparison[1]
      candidates_data_treatment_B_sig.2$method <- "hm450"
      candidates_data_treatment_B_sig.2$sub.method <- "B"
      row.names(candidates_data_treatment_B_sig.2) <- NULL
      ##candidates_data_treatment_C_sig <-
      ##  na.omit(candidates_data_treatment_C_sig)
      ##Workflow C ends up with some NB rows


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

      candidates_data_treatment_B_sig.1 <-
        candidates_data_treatment_B_sig.1[,preferred.column.names]
      candidates_data_treatment_B_sig.2 <-
        candidates_data_treatment_B_sig.2[,preferred.column.names]

      seg.out <- list(candidates_data_treatment_B_sig.1,
                      candidates_data_treatment_B_sig.2)

      names(seg.out) <- c(paste0(hm450.form.comparison[1],"_1"),
                          paste0(hm450.form.comparison[1],"_2"))

    }else if(hm450.form.workflow=="C"){

      ##sub routine C
      candidates_data_treatment_C.1 <- hm450.form.seg[[1]]
      candidates_data_treatment_C.2 <- hm450.form.seg[[2]]
      rm(hm450.form.seg)

      candidates_data_treatment_C_sig.1 <-
        candidates_data_treatment_C_sig.1[
          candidates_data_treatment_C.1$p.val <= 0.05,]
      candidates_data_treatment_C_sig.2 <-
        candidates_data_treatment_C_sig.2[
          candidates_data_treatment_C.2$p.val <= 0.05,]

      rm(candidates_data_treatment_C.1)
      rm(candidates_data_treatment_C.2)

      colnames(candidates_data_treatment_C_sig.1)[
        colnames(candidates_data_treatment_C_sig.1)=="smp"] <- "Sample_ID"
      colnames(candidates_data_treatment_C_sig.2)[
        colnames(candidates_data_treatment_C_sig.2)=="smp"] <- "Sample_ID"

      colnames(candidates_data_treatment_C_sig.1)[
        colnames(candidates_data_treatment_C_sig.1)=="chr"] <- "chrom"
      colnames(candidates_data_treatment_C_sig.2)[
        colnames(candidates_data_treatment_C_sig.2)=="chr"] <- "chrom"

      candidates_data_treatment_C_sig.1$chrom <-
        unlist(strsplit(candidates_data_treatment_C_sig.1$chrom,
                        split="chr"))[c(FCLSE,TRUE)]
      candidates_data_treatment_C_sig.2$chrom <-
        unlist(strsplit(candidates_data_treatment_C_sig.2$chrom,
                        split="chr"))[c(FCLSE,TRUE)]

      colnames(candidates_data_treatment_C_sig.1)[
        colnames(candidates_data_treatment_C_sig.1)=="startPos"] <- "loc.start"
      colnames(candidates_data_treatment_C_sig.2)[
        colnames(candidates_data_treatment_C_sig.2)=="startPos"] <- "loc.start"

      colnames(candidates_data_treatment_C_sig.1)[
        colnames(candidates_data_treatment_C_sig.1)=="endPos"] <- "loc.end"
      colnames(candidates_data_treatment_C_sig.2)[
        colnames(candidates_data_treatment_C_sig.2)=="endPos"] <- "loc.end"

      colnames(candidates_data_treatment_C_sig.1)[
        colnames(candidates_data_treatment_C_sig.1)=="mean"] <- "seg.mean"
      colnames(candidates_data_treatment_C_sig.2)[
        colnames(candidates_data_treatment_C_sig.2)=="mean"] <- "seg.mean"

      colnames(candidates_data_treatment_C_sig.1)[
        colnames(candidates_data_treatment_C_sig.1)=="median"] <- "seg.median"
      colnames(candidates_data_treatment_C_sig.2)[
        colnames(candidates_data_treatment_C_sig.2)=="median"] <- "seg.median"

      colnames(candidates_data_treatment_C_sig.1)[
        colnames(candidates_data_treatment_C_sig.1)=="p.val"] <- "pval"
      colnames(candidates_data_treatment_C_sig.2)[
        colnames(candidates_data_treatment_C_sig.2)=="p.val"] <- "pval"

      candidates_data_treatment_C_sig.1$num.mark  <- NC
      candidates_data_treatment_C_sig.1$bstat     <- NC
      candidates_data_treatment_C_sig.1$state <-
        round(2^candidates_data_treatment_C_sig.1$seg.mean * 2)
      candidates_data_treatment_C_sig.1$state[
        candidates_data_treatment_C_sig.1$state > 4] <- 4
      candidates_data_treatment_C_sig.1$treatment <- hm450.form.comparison[1]
      candidates_data_treatment_C_sig.1$method <- "hm450"
      candidates_data_treatment_C_sig.1$sub.method <- "C"
      row.names(candidates_data_treatment_C_sig.1) <- NULL
      candidates_data_treatment_C_sig <-
        na.omit(candidates_data_treatment_C_sig)
      ##Workflow C ends up with some NC rows

      candidates_data_treatment_C_sig.2$num.mark  <- NC
      candidates_data_treatment_C_sig.2$bstat     <- NC
      candidates_data_treatment_C_sig.2$state <-
        round(2^candidates_data_treatment_C_sig.2$seg.mean * 2)
      candidates_data_treatment_C_sig.2$state[
        candidates_data_treatment_C_sig.2$state > 4] <- 4
      candidates_data_treatment_C_sig.2$treatment <- hm450.form.comparison[1]
      candidates_data_treatment_C_sig.2$method <- "hm450"
      candidates_data_treatment_C_sig.2$sub.method <- "C"
      row.names(candidates_data_treatment_C_sig.2) <- NULL
      candidates_data_treatment_C_sig <-
        na.omit(candidates_data_treatment_C_sig)
      ##Workflow C ends up with some NC rows

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

      candidates_data_treatment_C_sig.1 <-
        candidates_data_treatment_C_sig.1[,preferred.column.names]
      candidates_data_treatment_C_sig.2 <-
        candidates_data_treatment_C_sig.2[,preferred.column.names]

      seg.out <- list(candidates_data_treatment_C_sig.1,
                      candidates_data_treatment_C_sig.2)

      names(seg.out) <- c(paste0(hm450.form.comparison[1],"_1"),
                          paste0(hm450.form.comparison[1],"_2"))
  }

}else{

  if(is.null(hm450.form.split.by)){

    if(hm450.form.workflow=="A"){

      ##sub routine A

      candidates_data_treatment_A_sig <-
        candidates_data_treatment_A[candidates_data_treatment_B$p.val <= 0.05,]

      rm(candidates_data_treatment_A)

      colnames(candidates_data_treatment_A_sig)[
        colnames(candidates_data_treatment_A_sig)=="smp"] <- "Sample_ID"

      colnames(candidates_data_treatment_A_sig)[
        colnames(candidates_data_treatment_A_sig)=="chr"] <- "chrom"

      candidates_data_treatment_A_sig$chrom <-
        unlist(strsplit(candidates_data_treatment_A_sig$chrom,
                        split="chr"))[c(FALSE,TRUE)]

      colnames(candidates_data_treatment_A_sig)[
        colnames(candidates_data_treatment_A_sig)=="startPos"] <- "loc.start"

      colnames(candidates_data_treatment_A_sig)[
        colnames(candidates_data_treatment_A_sig)=="endPos"] <- "loc.end"

      colnames(candidates_data_treatment_A_sig)[
        colnames(candidates_data_treatment_A_sig)=="mean"] <- "seg.mean"

      colnames(candidates_data_treatment_A_sig)[
        colnames(candidates_data_treatment_A_sig)=="median"] <- "seg.median"

      colnames(candidates_data_treatment_A_sig)[
        colnames(candidates_data_treatment_A_sig)=="p.val"] <- "pval"

      candidates_data_treatment_A_sig$num.mark  <- NA
      candidates_data_treatment_A_sig$bstat     <- NA
      candidates_data_treatment_A_sig$state <-
        round(2^candidates_data_treatment_A_sig$seg.mean * 2)
      candidates_data_treatment_A_sig$state[
        candidates_data_treatment_A_sig$state > 4] <- 4
      candidates_data_treatment_A_sig$treatment <- hm450.form.comparison[1]
      candidates_data_treatment_A_sig$method <- "hm450"
      candidates_data_treatment_A_sig$sub.method <- "A"
      row.names(candidates_data_treatment_B_sig) <- NULL
      ##seg <- na.omit(seg) ##Workflow C ends up with some NA rows

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

      candidates_data_treatment_A_sig <-
        candidates_data_treatment_A_sig[,preferred.column.names]

      seg.out <- list(candidates_data_treatment_A_sig)

      names(seg.out) <- hm450.form.comparison[1]

    }else if(hm450.form.workflow=="B"){

      candidates_data_treatment_B <- hm450.form.seg[[1]]
      rm(hm450.form.seg)

      ##colnames(candidates_data_treatment_B)
      ##"chr"
      ##"startPos"
      ##"endPos"
      ##"median"
      ##"mean"
      ##"sd"
      ##"smp"
      ##"p.val"

      ##any(is.na=(candidates_data_treatment_B$chr))
      ##No NA for chr field which is good

      candidates_data_treatment_B_sig <-
        candidates_data_treatment_B[candidates_data_treatment_B$p.val <= 0.05,]

      rm(candidates_data_treatment_B)

      colnames(candidates_data_treatment_B_sig)[
        colnames(candidates_data_treatment_B_sig)=="smp"] <- "Sample_ID"

      colnames(candidates_data_treatment_B_sig)[
        colnames(candidates_data_treatment_B_sig)=="chr"] <- "chrom"

      candidates_data_treatment_B_sig$chrom <-
        unlist(strsplit(candidates_data_treatment_B_sig$chrom,
                        split="chr"))[c(FALSE,TRUE)]

      colnames(candidates_data_treatment_B_sig)[
        colnames(candidates_data_treatment_B_sig)=="startPos"] <- "loc.start"

      colnames(candidates_data_treatment_B_sig)[
        colnames(candidates_data_treatment_B_sig)=="endPos"] <- "loc.end"

      colnames(candidates_data_treatment_B_sig)[
        colnames(candidates_data_treatment_B_sig)=="mean"] <- "seg.mean"

      colnames(candidates_data_treatment_B_sig)[
        colnames(candidates_data_treatment_B_sig)=="median"] <- "seg.median"

      colnames(candidates_data_treatment_B_sig)[
        colnames(candidates_data_treatment_B_sig)=="p.val"] <- "pval"

      candidates_data_treatment_B_sig$num.mark  <- NA
      candidates_data_treatment_B_sig$bstat     <- NA
      candidates_data_treatment_B_sig$state <-
        round(2^candidates_data_treatment_B_sig$seg.mean * 2)
      candidates_data_treatment_B_sig$state[
        candidates_data_treatment_B_sig$state > 4] <- 4
      candidates_data_treatment_B_sig$treatment <- hm450.form.comparison[1]
      candidates_data_treatment_B_sig$method <- "hm450"
      candidates_data_treatment_B_sig$sub.method <- "B"
      row.names(candidates_data_treatment_B_sig) <- NULL
      ##seg <- na.omit(seg) ##Workflow C ends up with some NA rows

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

      candidates_data_treatment_B_sig <-
        candidates_data_treatment_B_sig[,preferred.column.names]

      seg.out <- list(candidates_data_treatment_B_sig)

      names(seg.out) <- hm450.form.comparison[1]

      ##annotation1$probe <- rownames(annotation1)
      ##seg.1$loc.start
      ##hm450.manifest.hg38$addressA
      ##hm450.manifest.hg38$probeStart
      ##hm450.manifest.hg38$probeEnd

    }else if(hm450.form.workflow=="C"){

      candidates_data_treatment_C <- hm450.form.seg[[1]]
      rm(hm450.form.seg)

      ##colnames(candidates_data_treatment_B)
      ##"chr"
      ##"startPos"
      ##"endPos"
      ##"median"
      ##"mean"
      ##"sd"
      ##"smp"
      ##"p.val"

      ##any(is.na=(candidates_data_treatment_B$chr))
      ##No NA for chr field which is good

      candidates_data_treatment_C_sig <-
        candidates_data_treatment_C[candidates_data_treatment_C$p.val <= 0.05,]

      rm(candidates_data_treatment_C)

      colnames(candidates_data_treatment_C_sig)[
        colnames(candidates_data_treatment_C_sig)=="smp"] <- "Sample_ID"

      colnames(candidates_data_treatment_C_sig)[
        colnames(candidates_data_treatment_C_sig)=="chr"] <- "chrom"

      candidates_data_treatment_C_sig$chrom <-
        unlist(strsplit(candidates_data_treatment_C_sig$chrom,
                        split="chr"))[c(FALSE,TRUE)]

      colnames(candidates_data_treatment_C_sig)[
        colnames(candidates_data_treatment_C_sig)=="startPos"] <- "loc.start"

      colnames(candidates_data_treatment_C_sig)[
        colnames(candidates_data_treatment_C_sig)=="endPos"] <- "loc.end"

      colnames(candidates_data_treatment_C_sig)[
        colnames(candidates_data_treatment_C_sig)=="mean"] <- "seg.mean"

      colnames(candidates_data_treatment_C_sig)[
        colnames(candidates_data_treatment_C_sig)=="median"] <- "seg.median"

      colnames(candidates_data_treatment_C_sig)[
        colnames(candidates_data_treatment_C_sig)=="p.val"] <- "pval"

      candidates_data_treatment_C_sig$num.mark  <- NA
      candidates_data_treatment_C_sig$bstat     <- NA
      candidates_data_treatment_C_sig$state <-
        round(2^candidates_data_treatment_C_sig$seg.mean * 2)
      candidates_data_treatment_C_sig$state[
        candidates_data_treatment_C_sig$state > 4] <- 4
      candidates_data_treatment_C_sig$treatment <- hm450.form.comparison[1]
      candidates_data_treatment_C_sig$method <- "hm450"
      candidates_data_treatment_C_sig$sub.method <- "B"
      row.names(candidates_data_treatment_C_sig) <- NULL

      candidates_data_treatment_C_sig <-
        na.omit(candidates_data_treatment_C_sig)
      ##Workflow C ends up with some NA rows

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

      candidates_data_treatment_C_sig <-
        candidates_data_treatment_C_sig[,preferred.column.names]

      seg.out <- list(candidates_data_treatment_C_sig)

      names(seg.out) <- hm450.form.comparison[1]

    }else{

      stop(paste0("Error: need to select a ",
                  "proper sub workflow for 450k",
                  " <hm450.form.workflow>"))

    }

  }else{

    if(hm450.form.workflow=="A"){

      ##sub routine A
      candidates_data_treatment_A.1 <- hm450.form.seg[[1]]
      candidates_data_treatment_A.2 <- hm450.form.seg[[2]]
      rm(hm450.form.seg)

      candidates_data_treatment_A_sig.1 <-
        candidates_data_treatment_A_sig.1[
          candidates_data_treatment_A.1$p.val <= 0.05,]
      candidates_data_treatment_A_sig.2 <-
        candidates_data_treatment_A_sig.2[
          candidates_data_treatment_A.2$p.val <= 0.05,]

      rm(candidates_data_treatment_A.1)
      rm(candidates_data_treatment_A.2)

      colnames(candidates_data_treatment_A_sig.1)[
        colnames(candidates_data_treatment_A_sig.1)=="smp"] <- "Sample_ID"
      colnames(candidates_data_treatment_A_sig.2)[
        colnames(candidates_data_treatment_A_sig.2)=="smp"] <- "Sample_ID"

      colnames(candidates_data_treatment_A_sig.1)[
        colnames(candidates_data_treatment_A_sig.1)=="chr"] <- "chrom"
      colnames(candidates_data_treatment_A_sig.2)[
        colnames(candidates_data_treatment_A_sig.2)=="chr"] <- "chrom"

      candidates_data_treatment_A_sig.1$chrom <-
        unlist(strsplit(candidates_data_treatment_A_sig.1$chrom,
                        split="chr"))[c(FALSE,TRUE)]
      candidates_data_treatment_A_sig.2$chrom <-
        unlist(strsplit(candidates_data_treatment_A_sig.2$chrom,
                        split="chr"))[c(FALSE,TRUE)]

      colnames(candidates_data_treatment_A_sig.1)[
        colnames(candidates_data_treatment_A_sig.1)=="startPos"] <- "loc.start"
      colnames(candidates_data_treatment_A_sig.2)[
        colnames(candidates_data_treatment_A_sig.2)=="startPos"] <- "loc.start"

      colnames(candidates_data_treatment_A_sig.1)[
        colnames(candidates_data_treatment_A_sig.1)=="endPos"] <- "loc.end"
      colnames(candidates_data_treatment_A_sig.2)[
        colnames(candidates_data_treatment_A_sig.2)=="endPos"] <- "loc.end"

      colnames(candidates_data_treatment_A_sig.1)[
        colnames(candidates_data_treatment_A_sig.1)=="mean"] <- "seg.mean"
      colnames(candidates_data_treatment_A_sig.2)[
        colnames(candidates_data_treatment_A_sig.2)=="mean"] <- "seg.mean"

      colnames(candidates_data_treatment_A_sig.1)[
        colnames(candidates_data_treatment_A_sig.1)=="median"] <- "seg.median"
      colnames(candidates_data_treatment_A_sig.2)[
        colnames(candidates_data_treatment_A_sig.2)=="median"] <- "seg.median"

      colnames(candidates_data_treatment_A_sig.1)[
        colnames(candidates_data_treatment_A_sig.1)=="p.val"] <- "pval"
      colnames(candidates_data_treatment_A_sig.2)[
        colnames(candidates_data_treatment_A_sig.2)=="p.val"] <- "pval"

      candidates_data_treatment_A_sig.1$num.mark  <- NA
      candidates_data_treatment_A_sig.1$bstat     <- NA
      candidates_data_treatment_A_sig.1$state <-
        round(2^candidates_data_treatment_A_sig.1$seg.mean * 2)
      candidates_data_treatment_A_sig.1$state[
        candidates_data_treatment_A_sig.1$state > 4] <- 4
      candidates_data_treatment_A_sig.1$treatment <- hm450.form.comparison[1]
      candidates_data_treatment_A_sig.1$method <- "hm450"
      candidates_data_treatment_A_sig.1$sub.method <- "A"
      row.names(candidates_data_treatment_A_sig.1) <- NULL
      ##candidates_data_treatment_C_sig <-
      ##  na.omit(candidates_data_treatment_C_sig)
      ##Workflow C ends up with some NA rows

      candidates_data_treatment_A_sig.2$num.mark  <- NA
      candidates_data_treatment_A_sig.2$bstat     <- NA
      candidates_data_treatment_A_sig.2$state <-
        round(2^candidates_data_treatment_A_sig.2$seg.mean * 2)
      candidates_data_treatment_A_sig.2$state[
        candidates_data_treatment_A_sig.2$state > 4] <- 4
      candidates_data_treatment_A_sig.2$treatment <- hm450.form.comparison[1]
      candidates_data_treatment_A_sig.2$method <- "hm450"
      candidates_data_treatment_A_sig.2$sub.method <- "A"
      row.names(candidates_data_treatment_A_sig.2) <- NULL
      ##candidates_data_treatment_C_sig <-
      ##  na.omit(candidates_data_treatment_C_sig)
      ##Workflow C ends up with some NA rows


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

      candidates_data_treatment_A_sig.1 <-
        candidates_data_treatment_A_sig.1[,preferred.column.names]
      candidates_data_treatment_A_sig.2 <-
        candidates_data_treatment_A_sig.2[,preferred.column.names]

      seg.out <- list(candidates_data_treatment_A_sig.1,
                      candidates_data_treatment_A_sig.2)

      names(seg.out) <- c(paste0(hm450.form.comparison[1],"_1"),
                          paste0(hm450.form.comparison[1],"_2"))

    }else if(hm450.form.workflow=="B"){

      ##sub routine B
      candidates_data_treatment_B.1 <- hm450.form.seg[[1]]
      candidates_data_treatment_B.2 <- hm450.form.seg[[2]]
      rm(hm450.form.seg)

      candidates_data_treatment_B_sig.1 <-
        candidates_data_treatment_B_sig.1[
          candidates_data_treatment_B.1$p.val <= 0.05,]
      candidates_data_treatment_B_sig.2 <-
        candidates_data_treatment_B_sig.2[
          candidates_data_treatment_B.2$p.val <= 0.05,]

      rm(candidates_data_treatment_B.1)
      rm(candidates_data_treatment_B.2)

      colnames(candidates_data_treatment_B_sig.1)[
        colnames(candidates_data_treatment_B_sig.1)=="smp"] <- "Sample_ID"
      colnames(candidates_data_treatment_B_sig.2)[
        colnames(candidates_data_treatment_B_sig.2)=="smp"] <- "Sample_ID"

      colnames(candidates_data_treatment_B_sig.1)[
        colnames(candidates_data_treatment_B_sig.1)=="chr"] <- "chrom"
      colnames(candidates_data_treatment_B_sig.2)[
        colnames(candidates_data_treatment_B_sig.2)=="chr"] <- "chrom"

      candidates_data_treatment_B_sig.1$chrom <-
        unlist(strsplit(candidates_data_treatment_B_sig.1$chrom,
                        split="chr"))[c(FBLSE,TRUE)]
      candidates_data_treatment_B_sig.2$chrom <-
        unlist(strsplit(candidates_data_treatment_B_sig.2$chrom,
                        split="chr"))[c(FBLSE,TRUE)]

      colnames(candidates_data_treatment_B_sig.1)[
        colnames(candidates_data_treatment_B_sig.1)=="startPos"] <- "loc.start"
      colnames(candidates_data_treatment_B_sig.2)[
        colnames(candidates_data_treatment_B_sig.2)=="startPos"] <- "loc.start"

      colnames(candidates_data_treatment_B_sig.1)[
        colnames(candidates_data_treatment_B_sig.1)=="endPos"] <- "loc.end"
      colnames(candidates_data_treatment_B_sig.2)[
        colnames(candidates_data_treatment_B_sig.2)=="endPos"] <- "loc.end"

      colnames(candidates_data_treatment_B_sig.1)[
        colnames(candidates_data_treatment_B_sig.1)=="mean"] <- "seg.mean"
      colnames(candidates_data_treatment_B_sig.2)[
        colnames(candidates_data_treatment_B_sig.2)=="mean"] <- "seg.mean"

      colnames(candidates_data_treatment_B_sig.1)[
        colnames(candidates_data_treatment_B_sig.1)=="median"] <- "seg.median"
      colnames(candidates_data_treatment_B_sig.2)[
        colnames(candidates_data_treatment_B_sig.2)=="median"] <- "seg.median"

      colnames(candidates_data_treatment_B_sig.1)[
        colnames(candidates_data_treatment_B_sig.1)=="p.val"] <- "pval"
      colnames(candidates_data_treatment_B_sig.2)[
        colnames(candidates_data_treatment_B_sig.2)=="p.val"] <- "pval"

      candidates_data_treatment_B_sig.1$num.mark  <- NB
      candidates_data_treatment_B_sig.1$bstat     <- NB
      candidates_data_treatment_B_sig.1$state <-
        round(2^candidates_data_treatment_B_sig.1$seg.mean * 2)
      candidates_data_treatment_B_sig.1$state[
        candidates_data_treatment_B_sig.1$state > 4] <- 4
      candidates_data_treatment_B_sig.1$treatment <- hm450.form.comparison[1]
      candidates_data_treatment_B_sig.1$method <- "hm450"
      candidates_data_treatment_B_sig.1$sub.method <- "B"
      row.names(candidates_data_treatment_B_sig.1) <- NULL
      ##candidates_data_treatment_C_sig <-
      ##  na.omit(candidates_data_treatment_C_sig)
      ##Workflow C ends up with some NB rows

      candidates_data_treatment_B_sig.2$num.mark  <- NB
      candidates_data_treatment_B_sig.2$bstat     <- NB
      candidates_data_treatment_B_sig.2$state <-
        round(2^candidates_data_treatment_B_sig.2$seg.mean * 2)
      candidates_data_treatment_B_sig.2$state[
        candidates_data_treatment_B_sig.2$state > 4] <- 4
      candidates_data_treatment_B_sig.2$treatment <- hm450.form.comparison[1]
      candidates_data_treatment_B_sig.2$method <- "hm450"
      candidates_data_treatment_B_sig.2$sub.method <- "B"
      row.names(candidates_data_treatment_B_sig.2) <- NULL
      ##candidates_data_treatment_C_sig <-
      ##  na.omit(candidates_data_treatment_C_sig)
      ##Workflow C ends up with some NB rows


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

      candidates_data_treatment_B_sig.1 <-
        candidates_data_treatment_B_sig.1[,preferred.column.names]
      candidates_data_treatment_B_sig.2 <-
        candidates_data_treatment_B_sig.2[,preferred.column.names]

      seg.out <- list(candidates_data_treatment_B_sig.1,
                      candidates_data_treatment_B_sig.2)

      names(seg.out) <- c(paste0(hm450.form.comparison[1],"_1"),
                          paste0(hm450.form.comparison[1],"_2"))

    }else if(hm450.form.workflow=="C"){

      ##sub routine C
      candidates_data_treatment_C.1 <- hm450.form.seg[[1]]
      candidates_data_treatment_C.2 <- hm450.form.seg[[2]]
      rm(hm450.form.seg)

      candidates_data_treatment_C_sig.1 <-
        candidates_data_treatment_C_sig.1[
          candidates_data_treatment_C.1$p.val <= 0.05,]
      candidates_data_treatment_C_sig.2 <-
        candidates_data_treatment_C_sig.2[
          candidates_data_treatment_C.2$p.val <= 0.05,]

      rm(candidates_data_treatment_C.1)
      rm(candidates_data_treatment_C.2)

      colnames(candidates_data_treatment_C_sig.1)[
        colnames(candidates_data_treatment_C_sig.1)=="smp"] <- "Sample_ID"
      colnames(candidates_data_treatment_C_sig.2)[
        colnames(candidates_data_treatment_C_sig.2)=="smp"] <- "Sample_ID"

      colnames(candidates_data_treatment_C_sig.1)[
        colnames(candidates_data_treatment_C_sig.1)=="chr"] <- "chrom"
      colnames(candidates_data_treatment_C_sig.2)[
        colnames(candidates_data_treatment_C_sig.2)=="chr"] <- "chrom"

      candidates_data_treatment_C_sig.1$chrom <-
        unlist(strsplit(candidates_data_treatment_C_sig.1$chrom,
                        split="chr"))[c(FCLSE,TRUE)]
      candidates_data_treatment_C_sig.2$chrom <-
        unlist(strsplit(candidates_data_treatment_C_sig.2$chrom,
                        split="chr"))[c(FCLSE,TRUE)]

      colnames(candidates_data_treatment_C_sig.1)[
        colnames(candidates_data_treatment_C_sig.1)=="startPos"] <- "loc.start"
      colnames(candidates_data_treatment_C_sig.2)[
        colnames(candidates_data_treatment_C_sig.2)=="startPos"] <- "loc.start"

      colnames(candidates_data_treatment_C_sig.1)[
        colnames(candidates_data_treatment_C_sig.1)=="endPos"] <- "loc.end"
      colnames(candidates_data_treatment_C_sig.2)[
        colnames(candidates_data_treatment_C_sig.2)=="endPos"] <- "loc.end"

      colnames(candidates_data_treatment_C_sig.1)[
        colnames(candidates_data_treatment_C_sig.1)=="mean"] <- "seg.mean"
      colnames(candidates_data_treatment_C_sig.2)[
        colnames(candidates_data_treatment_C_sig.2)=="mean"] <- "seg.mean"

      colnames(candidates_data_treatment_C_sig.1)[
        colnames(candidates_data_treatment_C_sig.1)=="median"] <- "seg.median"
      colnames(candidates_data_treatment_C_sig.2)[
        colnames(candidates_data_treatment_C_sig.2)=="median"] <- "seg.median"

      colnames(candidates_data_treatment_C_sig.1)[
        colnames(candidates_data_treatment_C_sig.1)=="p.val"] <- "pval"
      colnames(candidates_data_treatment_C_sig.2)[
        colnames(candidates_data_treatment_C_sig.2)=="p.val"] <- "pval"

      candidates_data_treatment_C_sig.1$num.mark  <- NC
      candidates_data_treatment_C_sig.1$bstat     <- NC
      candidates_data_treatment_C_sig.1$state <-
        round(2^candidates_data_treatment_C_sig.1$seg.mean * 2)
      candidates_data_treatment_C_sig.1$state[
        candidates_data_treatment_C_sig.1$state > 4] <- 4
      candidates_data_treatment_C_sig.1$treatment <- hm450.form.comparison[1]
      candidates_data_treatment_C_sig.1$method <- "hm450"
      candidates_data_treatment_C_sig.1$sub.method <- "C"
      row.names(candidates_data_treatment_C_sig.1) <- NULL
      candidates_data_treatment_C_sig <-
        na.omit(candidates_data_treatment_C_sig)
      ##Workflow C ends up with some NC rows

      candidates_data_treatment_C_sig.2$num.mark  <- NC
      candidates_data_treatment_C_sig.2$bstat     <- NC
      candidates_data_treatment_C_sig.2$state <-
        round(2^candidates_data_treatment_C_sig.2$seg.mean * 2)
      candidates_data_treatment_C_sig.2$state[
        candidates_data_treatment_C_sig.2$state > 4] <- 4
      candidates_data_treatment_C_sig.2$treatment <- hm450.form.comparison[1]
      candidates_data_treatment_C_sig.2$method <- "hm450"
      candidates_data_treatment_C_sig.2$sub.method <- "C"
      row.names(candidates_data_treatment_C_sig.2) <- NULL
      candidates_data_treatment_C_sig <-
        na.omit(candidates_data_treatment_C_sig)
      ##Workflow C ends up with some NC rows


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

      candidates_data_treatment_C_sig.1 <-
        candidates_data_treatment_C_sig.1[,preferred.column.names]
      candidates_data_treatment_C_sig.2 <-
        candidates_data_treatment_C_sig.2[,preferred.column.names]

      seg.out <- list(candidates_data_treatment_C_sig.1,
                      candidates_data_treatment_C_sig.2)

      names(seg.out) <- c(paste0(hm450.form.comparison[1],"_1"),
                          paste0(hm450.form.comparison[1],"_2"))

  }

}

if(hm450.form.save.seg==TRUE){
    save(seg.out,
         file=paste0(hm450.form.output.dir,
                     .Platform$file.sep,
                     "seg.out.RData"))
}

return(seg.out)

}
