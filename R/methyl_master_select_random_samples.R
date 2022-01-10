#!/usr/bin/env Rscript

#' @title methyl_master_load_data_objects
#' @description Possible function for loading data objects
#' Michael Mariani PhD Dartmouth College 2021
#' norm.path<-"G:\\My Drive\\dartmouth\\salas_lab_working\\cnv\\BLCA_normal.txt"
#' tum.path<-"G:\\My Drive\\dartmouth\\salas_lab_working\\cnv\\BLCA_normal.txt"
#' @param normal.path The normal path
#' @param tumor.path The tumor path
#' @return Return randomly sampled TCGA IDAT samples
#' @export
methyl_master_select_random_samples <- function(normal.path,
                                                tumor.path)
{

normal.df <-
  read.table(normal.path,
             sep="\t",
             header=TRUE,
             stringsAsFactors = FALSE
             )

tumor.df <-
  read.table(tumor.path,
             sep="\t",
             header=TRUE,
             stringsAsFactors = FALSE
             )

normal.male   <- normal.df[normal.df$gender_reported=="male","X"]
normal.female <- normal.df[normal.df$gender_reported=="female","X"]
tumor.male    <- tumor.df[normal.df$gender_reported=="male","X"]
tumor.female  <- tumor.df[normal.df$gender_reported=="female","X"]

set.seed(42)
normal.male   <- sample(normal.male,   size=10, replace=FALSE)
normal.female <- sample(normal.female, size=10, replace=FALSE)
tumor.male    <- sample(tumor.male,    size=10, replace=FALSE)
tumor.female  <- sample(tumor.female,  size=10, replace=FALSE)

cat(normal.male,sep="\n")
cat(normal.female,sep="\n")
cat(tumor.male,sep="\n")
cat(tumor.female,sep="\n")

}

