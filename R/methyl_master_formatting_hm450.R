#!/usr/bin/env Rscript

#' @title methyl_master_formatting_hm450
#' @description Formatting the hm450 output to prepare for plotting etc.
#' @param hm450.form.seg The hm450 routine seg results as input to be formatted
#' @param hm450.form.output.dir The output dir for the formatting results
#' @param hm450.form.sample.sheet.path The MehtylMaster sample sheet path
#' @param hm450.form.thresholds The thresholds used to determine the CNV
#' state, if NULL, the equation  seg.state <- round(2^seg.means * 2) is used
#' @param hm450.form.workflow The specific HM450 workflow that was used
#' @param hm450.form.comparison The MethylMaster two-element comparison vector
#' @param hm450.form.save.seg Whether to save the formatted HM450 results
#' @param ... Additional parameters passed to methyl_master_formatting_hm450
#' @import CNVRanger
#' @import matter
#' @importFrom magrittr %>%
#' @return #Formatted seg list object for overlaps and visualizing
#' @export
methyl_master_formatting_hm450 <- function(hm450.form.seg=NULL,
                                    hm450.form.output.dir=getwd(),
                                    hm450.form.sample.sheet.path=NULL,
                                    hm450.form.reference="internal",
                                    hm450.form.thresholds=NULL,
                                    hm450.form.split.by=NULL,
                                    hm450.form.workflow="B",
                                    hm450.form.comparison=NULL,
                                    hm450.form.save.seg=FALSE,
                                    ...
                                  ){

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
    calc_seg_state(seg.means = candidates_data_treatment_B_sig$seg.mean,
                               upper.thresh = 4,
                               cutoff = hm450.form.thresholds)
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

if(hm450.form.save.seg==TRUE){
    save(seg.out,
         file=paste0(hm450.form.output.dir,
                     .Platform$file.sep,
                     "seg.out.RData"))
}

return(seg.out)

}
