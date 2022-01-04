#!/usr/bin/env Rscript

#' @title methyl_master_get_cord_sex
#' @description get sex information for cord samples
#' Michael Mariani PhD Dartmouth College 2022
#' @param mmcs.idat.files.dir
#' @param mmcs.platform
#' @param mmcs.output.file.name
#' @import sesame
#' @import sesameData
#' @import ExperimentHub
#' @return #outputs a file of the sexes for select cord samples
#' @export
methyl_measter_get_cord_sex <- function(mmcs.idat.files.dir=getwd(),
                                        mmcs.platform='EPIC',
                                        mmcs.output.file.name=NULL
                                        ){

mmcs.idat.files.dir <- paste0("G:\\My Drive\\dartmouth",
                           "\\salas_lab_working\\cnv",
                           "\\blca_idat_files_testing")

mmcs.platform <- "EPIC"

mmcs.output.file.name <- paste0("G:\\My Drive\\dartmouth",
                                "\\salas_lab_working\\cnv",
                                "\\cord_sex_info.txt")

############################# PROCESS CORD ###################################

## Here we can use sesame to infer the sex and karyotype of samples
## such as the cord samples from GSE153668
## the foreach do paradigm is preserved in case we hook it up to
## future parallelization

setExperimentHubOption("CACHE",
                       mmcs.idat.files.dir)

ExperimentHub::ExperimentHub()

idat_prefixes <- sesame::searchIDATprefixes(mmcs.idat.files.dir,
                                            recursive=TRUE)

idat_prefixes <- idat_prefixes[grepl(idat_prefixes,pattern = "GSM")]

sesameData::sesameDataCacheAll()

sesame_sset <- sesame::openSesame(idat_prefixes,
                                  mask = TRUE,
                                  sum.TypeI = TRUE,
                                  platform = mmcs.platform,
                                  what="sigset")

sex_cord <- foreach(i = 1:length(names(sesame_sset))) %do% {
  sesame::inferSex(sesame_sset[[i]])
}

karyotype_cord <- foreach(i = 1:length(names(sesame_sset))) %do%{
  sesame::inferSexKaryotypes(sesame_sset[[i]])
}

length(sex_cord)
length(karyotype_cord)

sex.frame <- data.frame(sample=names(sesame_sset),
                        sex=unlist(sex_cord),
                        karyotype=unlist(karyotype_cord),
                        stringsAsFactors=FALSE)

write.table(sex.frame,
            file=mmcs.output.file.name,
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE,
            sep="\t")

}
