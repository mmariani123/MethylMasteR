#!/usr/bin/env Rscript

#' @title methyl_master_tcga_clin_data
#' @description get primary ids from sample names from tcga clinical file
#' @param clin.tcga.file.path The path to the tcga files
#' @param clin.sample.names The sample names
#' @param clin.copy.files Whether to copy the files to a different directory
#' @param clin.sample.sheet Whether to use a sample sheet specifying files
#' @param clin.idat.dir The idat dir
#' @param clin.sub.dir The sub directory
#' @importFrom readxl read_excel
#' @return outputs a vector or excel of
#' primary ids for input tcga smamples names
#' @export
methyl_master_tcga_clin_data <- function(clin.tcga.file.path=NULL,
                                         clin.sample.names=NULL,
                                         clin.copy.files=FALSE,
                                         clin.sample.sheet=TRUE,
                                         clin.idat.dir=getwd(),
                                         clin.sub.dir=NULL){

  clin.data <- read.table(file = clin.tcga.file.path,
                          header=TRUE,
                          sep=",",
                          stringsAsFactors = FALSE)

  if(clin.sample.sheet==TRUE){

    select.cbio.samples <-
      readxl::read_excel(path = clin.tcga.file.path,
                       col_names = TRUE,
                       sheet = 1)

  }else{

    select.cbio.samples <- clin.sample.names

  }

  select.cbio.primaries <-
    clin.data[clin.data$submitter_id %in%
                select.cbio.samples$sample,"primary"]

  kirc.sample.paths <- list.files(idat.dir,
                                  full.names = TRUE,
                                  recursive = TRUE)

  ##unique.primary.samples <-
  ##  unique(basename(gsub("_Grn.idat$|_Red.idat","",kirc.sample.paths)))
  ##length(unique.primary.samples) ##487
  ##kirc.sample.paths[grepl(select.cbio.primaries, select.cbio.primaries)]

  copy.samples <-
    unlist(lapply(select.cbio.primaries,
                  FUN = function(x,kirc.sample.paths){
                    grep(x,
                         kirc.sample.paths,
                         value=TRUE)
                    }
                  )
           )

  if(clin.copy.files==TRUE){

    file.copy(copy.samples, clin.sub.dir)

  }else{

    return(copy.samples)

  }

}
