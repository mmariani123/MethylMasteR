#!/usr/bin/env Rscript

#' @title methyl_master_prepare_sample_sheet
#' @description compare results from different routine outputs
#' For example, create a MethylMaster style sample sheet
#' from selected 3p-del (positive control) KIRC samples
#' from the firehose collection on CBioPortal
#' @param prep.output.path.name The output file destination
#' @param prep.file.sep The file separator to use
#' @param prep.idat.dir The .idat file directory
#' @param prep.sample.sheet The sample sheet to use
#' @param prep.primary.samples The specified primary sample IDs
#' @param prep.sample.group The specified sample groups
#' @param prep.platform The platform
#' @param prep.batch the batch variable
#' @return A sample sheet for use with MethylMaster
#' @export
methyl_master_prepare_sample_sheet <- function(prep.output.path.name=NULL,
                                               prep.file.sep=NULL,
                                               prep.idat.dir=NULL,
                                               prep.sample.sheet=NULL,
                                               prep.primary.samples=NULL,
                                               prep.sample.group=NULL,
                                               prep.platform=NULL,
                                               prep.batch=1
                                              ){

file.sep <- prep.file.sep

idat.dir.base <- list.dirs(idat.dir)

colm.names <- c("Sample_Name",
                "Sample_Group",
                "Sample_Plate",
                "Pool_ID",
                "Sample_Well",
                "Sentrix_Position",
                "Sentrix_ID",
                "Basename",
                "Batch",
                "gender_reported"
                )

##sample.sheet.df <- data.frame(matrix(data=0,
##                                     ncol = 9,
##                                     nrow = length(select.cbio.primaries)),
##                              stringsAsFactors = FALSE)

df.num.rows <- length(prep.primary.samples)
sample.sheet.df <- data.frame(Sample_Name     = character(length = df.num.rows),
                              Sample_Group    = character(length = df.num.rows),
                              Sample_Plate    = character(length = df.num.rows),
                              Pool_ID         = character(length = df.num.rows),
                              Sample_Well     = character(length = df.num.rows),
                              Sentrix_Position= character(length = df.num.rows),
                              Sentrix_ID      = character(length = df.num.rows),
                              Basename        = character(length = df.num.rows),
                              stringsAsFactors= FALSE)

sample.sheet.df$Sample_Name      <- select.cbio.primaries
sample.sheet.df$Sample_Group     <- "tumor"
sample.sheet.df$Sample_Plate     <- NA
sample.sheet.df$Pool_ID          <- NA
sample.sheet.df$Sample_Well      <- NA
sample.sheet.df$Sentrix_Position <- gsub(".*_","",prep.primary.samples)
sample.sheet.df$Sentrix_ID       <- gsub("_.*","",prep.primary.samples)
sample.sheet.df$Basename         <- paste0(idat.dir.base,
                                    prep.file.sep,
                                    prep.primary.samples)
sample.sheet.df$Platform <- prep.platform
sample.sheet.df$Batch <- prep.batch

write.csv(sample.sheet.df,
          file = output.path.name,
          row.names = FALSE,
          quote = FALSE)

}
