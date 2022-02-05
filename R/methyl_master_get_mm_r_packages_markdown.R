#---
#title: "Get MM R Packages"
#author: "Michael Mariani PhD"
#date: "2/5/2022"
#output: html_document
#---
#
#```{r setup, include=FALSE}
#
#knitr::opts_chunk$set(echo = TRUE)
#
#```
#
### Old packages for MethylMasteR
#
#Getting as many of the originally MethylMasteR packages as
#possible, focus now will be on getting evrything up to date.
#
#```{r installs, inlcude=TRUE}
#
#install.packages(pkgs="data.table",
#                 lib=paste0("C:\\Users\\Mike_2\\Desktop",
#                          "\\methyl_master_packages"),
#                 destdir = paste0("C:\\Users\\Mike_2\\Desktop",
#                          "\\methyl_master_packages"),
#                 update=FALSE,
#                 version="1.9.4")
#
#install.packages(pkgs="stats",
#                 lib=paste0("C:\\Users\\Mike_2\\Desktop",
#                          "\\methyl_master_packages"),
#                 destdir = paste0("C:\\Users\\Mike_2\\Desktop",
#                          "\\methyl_master_packages"),
#                 update=FALSE,
#                 version="1.16.0")
#
#devtools::install_version("stats",
#                          version = "1.16.0",
#                          upgrade="never",
#                          repos = "http://cran.us.r-project.org",
#                          type="source")
#
##Maybe just use newest version of champ? unless can find 2.20.1
##did find 2.20.0, I set min to 2.20.0
#BiocManager::install("ChAMP",
#                     lib=paste0("C:\\Users\\Mike_2\\Desktop",
#                                "\\methyl_master_packages"),
#                     destdir = paste0("C:\\Users\\Mike_2\\Desktop",
#                                "\\methyl_master_packages"),
#                     update=FALSE,
#                     version="2.20.1")
#
#download.file("https://www.bioconductor.org/packages/release/bioc/s#rc/contrib/ChAMP_2.24.0.tar.gz",
#              destdir=paste0("C:\\Users\\Mike_2\\Desktop",
#                                "\\methyl_master_packages"),
#              method="wget",
#              quiet = FALSE,
#              mode = "w",
#              cacheOK = TRUE,
#              extra = getOption("download.file.extra"),
#              headers = NULL)
#
#devtools::install_version("limma",
#                          version = "3.46.0",
#                          upgrade="never",
#                          repos = "http://cran.us.r-project.org",
#                          type="source")
#
##3.45.0 can't be foiund, jusr downloaded newest: 3.5.0
#download.file("https://www.bioconductor.org/packages/3.4/bioc/src/c#ontrib/limma_3.46.0.tar.gz",
#              destfile=paste0("C:\\Users\\Mike_2\\Desktop",
#                                "\\methyl_master_packages",
#                              "limma_3.46.0.tar.gz"),
#              method="auto",
#              quiet = FALSE,
#              mode = "w",
#              cacheOK = TRUE,
#              extra = options(repos = c(CRAN = #"https://cloud.r-project.org/")),
#              headers = NULL)
#
#download.file("https://cran.r-project.org/src/contrib/igraph_1.2.11#.tar.gz",
#              destfile=paste0("C:\\Users\\Mike_2\\Desktop",
#                                "\\methyl_master_packages",
#                              "\\igraph_1.2.6.tar.gz"),
#              method="curl",
#              quiet = FALSE,
#              mode = "w",
#              cacheOK = TRUE,
#              extra = getOption("download.file.extra"),
#              headers = NULL)
#
##I am not sure where to find 3.12.0
#download.file("https://bioconductor.org/packages/release/data/annot#ation/src/contrib/org.Hs.eg.db_3.14.0.tar.gz",
#              destfile=paste0("C:\\Users\\Mike_2\\Desktop",
#                                "\\methyl_master_packages",
#                              "\\org.Hs.eg.db_3.14.0.tar.gz"),
#              method="curl",
#              quiet = FALSE,
#              mode = "w",
#              cacheOK = TRUE,
#              extra = getOption("download.file.extra"),
#              headers = NULL)
#
#download.file("https://bioconductor.org/packages/release/data/annot#ation/src/contrib/org.Hs.eg.db_3.14.0.tar.gz",
#              destfile=paste0("C:\\Users\\Mike_2\\Desktop",
#                                "\\methyl_master_packages",
#                              "\\org.Hs.eg.db_3.14.0.tar.gz"),
#              method="curl",
#              quiet = FALSE,
#              mode = "w",
#              cacheOK = TRUE,
#              extra = getOption("download.file.extra"),
#              headers = NULL)
#
##marray (== 1.68.0)
##Could onlu find 1.7.1 and 1.7.2
#
##shape (== 1.4.6)
##Appears to be most recent on CRAN
#
##This will install all dependencies as well:
#withr::with_libpaths(
#  new = paste0("C:\\Users\\Mike_2\\Desktop",
#             "\\methyl_master_packages"),
#  devtools::install_github('sean-cho/Epicopy',
#  upgrade = "never"))
#
##Epicopy (== 0.99.0)
##I believe this output the newest .tar.gz
#download.file("https://github.com/sean-cho/Epicopy/branches/master.#tar.gz",
#              destfile=paste0("C:\\Users\\Mike_2\\Desktop",
#                              "\\methyl_master_packages",
#                              "\\Epicopy.tar.gz"),
#              method="curl",
#              quiet = FALSE,
#              mode = "w",
#              cacheOK = TRUE,
#              extra = getOption("download.file.extra"),
#              headers = NULL)

#foreach (== 1.5.1)
#Downloaded from CRAN archive

#ggplot2 (== 3.3.3)
#Downloaded from CRAN archive

#magrittr (== 2.0.1)
#Download from CRAN archive

#plyr (== 1.8.6)
#This is the newest one on the CRAN archive

#dplyr (== 1.0.7)
#This is the newest one on CRAN

#TCGAbiolinks (== 1.14.0)
#this is old version and don't know where to find
#so lets use the newest from Bioc 2.22.4

#ExperimentHub (== 1.16.1)
#Found on Bioc Archive

#sesame (== 1.8.10)
#I was able to find this

#sesameData (== 1.8.10)
#Lets use the newest 1.12 version

#minfi (== 1.18.0)
#I could find 1.18.6
#couldn't find so I downoaded the newest,
#1.4.0 from bioc

#GEOquery (== 2.62.1)
#Found in bioc archive

#DNAcopy (== 2.0)
#typo? 1.68.0 is the newest on Bioc

#CNVRanger (== 1.6.1)
#Found on Bioc Archive

#cnAnalysis450k (== 0.99.26)
#this is the version at mknoll/cnAnalysis450k
#download.file("https://github.com/mknoll/cnAnalysis450k/branches/ma#ster.tar.gz",
#              destfile=paste0("C:\\Users\\Mike_2\\Desktop",
#                              "\\methyl_master_packages",
#              "\\cnAnalysis450k_0.99.26.tar.gz"),
#              method="curl",
#              quiet = FALSE,
#              mode = "w",
#              cacheOK = TRUE,
#              extra = getOption("download.file.extra"),
#              headers = NULL)
#
##CNAclinic (== 1.0)
##this is the version at sdchandra/CNAclinic
#download.file("https://github.com/sdchandra/CNAclinic/branches/mast#er.tar.gz",
#              destfile=paste0("C:\\Users\\Mike_2\\Desktop",
#                              "\\methyl_master_packages",
#              "\\CNAclinic_1.0.tar.gz"),
#              method="curl",
#              quiet = FALSE,
#              mode = "w",
#              cacheOK = TRUE,
#              extra = getOption("download.file.extra"),
#              headers = NULL)

#IlluminaHumanMethylation450kanno.ilmn12.hg19 (== 0.6.0)
#Still newest on Bioc

#ENmix (== 1.26.10)
#Found on Bioc Archive

#profmem (== 0.6.0)
#current one on CRAN

#profvis (== 0.3.7)
#Current version on CRAN

#htmlwidgets (== 1.5.4)
#default version on CRAN

#ungeviz (== 0.1.0)
#this is the current github
#version at https://github.com/wilkelab/ungeviz

#download.file("https://github.com/wilkelab/ungeviz/blob/master.tar.#gz",
#              destfile=paste0("C:\\Users\\Mike_2\\Desktop",
#                              "\\methyl_master_packages",
#              "\\ungeviz_0.1.0.tar.gz"),
#              method="curl",
#              quiet = FALSE,
#              mode = "w",
#              cacheOK = TRUE,
#              extra = getOption("download.file.extra"),
#              headers = NULL)

#AnnotationDbi (== 1.52.0)
#cant find the aboe version so I just donwloaded
#the newest one 1.56.2

#Biobase (== 2.50.0)
#Found on Bioc Archive

#BiocGenerics (== 0.36.1)
#Found on Bioc Archive

#BiocParallel (== 1.24.1)
#Found on Bioc archive

#ChAMPdata (== 2.22.0)
#I can't find the version so I just use
#2.26.0 the newest from bioc

#CopyNumber450kCancer (== 1.0.6)
#https://github.com/cran/CopyNumber450kCancer
#This should be 1.0.4 I believe

#DelayedArray (== 0.16.3)
#found on Bioc archive

#GenomicRanges (== 1.42.0)
#Found on Bioc Archive

#Matrix (== 1.3-2)
#found on CRAN

#RaggedExperiment (== 1.14.2)
#Found on Bioc Archive

#S4Vectors (== 0.28.1)
#Found on Bioc Archive

#cowplot (== 1.1.1)
#This is the newest version on CRAN

#dendextend (== 1.15.1)
#Found on CRAN archive

#doParallel (== 1.0.16)
#Newest version on CRAN archive

#future (== 1.21.0)
#Found on CRAN archive

#gtools (== 3.8.2)
#Found on CRAn archive

#illuminaio (== 0.32.0)
#found on CRAN archive

#iterators (== 1.0.13)
#Found on CRAn archive

#matter (== 1.16.0)
#Found on BIOC archive

#pheatmap (== 1.0.12)
#newest version on CRAN

#ramify (== 0.3.3)
#newest version on CRAn

#readxl (== 1.3.1)
#Newest version on CRANs

#rlist (== 0.4.6.2)
#Newest version on CRAN

#scales (== 1.1.1)
#Newest version on CRAN

##```
