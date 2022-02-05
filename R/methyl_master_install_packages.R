#!/usr/bin/env Rscript

####install packages for MethylMasteR
##
##install.packages("data.table",
##                 version="1.9.4",
##                 lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                 update=FALSE,
##                 destdir="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib")
##
##install.packages("stats",
##                 version = "1.16.0",
##                 lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                 update=FALSE,
##                 destdir="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib")
####Warning in install.packages :
####  package 'stats' is a base package, and should not be updated
##
##BiocManager::install("ChAMP",
##                     version = "2.20.1",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
####'2.24.0'
##
####install.packages("limma",
####                     version = "3.46.0",
####                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
####                     update=FALSE,
####                     destdir="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib")
####'3.50.0'
##
####BiocManager::install("igraph",
####                     version = "1.2.6",
####                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
####                     update=FALSE)
##
####BiocManager::install("org.Hs.eg.db",
####                     version = "3.12.0",
####                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
####                     update=FALSE)
##
####BiocManager::install("marray",
####                     version = "1.68.0",
####                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
####                     update=FALSE)
##
##install.packages("shape",
##                 version = "1.4.6",
##                 lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                 update=FALSE,
##                 destdir="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib")
##
##BiocManager::install("Epicopy",
##                     version = "0.99.0",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("foreach",
##                     version = "1.5.1",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("ggplot2",
##                     version = "3.3.3",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("magrittr",
##                     version = "2.0.1",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("plyr",
##                     version = "1.8.6",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("dplyr",
##                     version = "1.0.7",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("TCGAbiolinks",
##                     version = "1.14.0",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("ExperimentHub",
##                     version = "1.16.1",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("sesame",
##                     version = "1.8.10",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("sesameData",
##                     version = "1.8.10",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("minfi",
##                     version = "1.18.0",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("GEOquery",
##                     version = "2.62.1",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("DNAcopy",
##                     version = "2.0",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("CNVRanger",
##                     version = "1.6.1",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("cnAnalysis450k",
##                     version = "0.99.26",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("CNAclinic",
##                     version = "1.0",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19",
##                     version = "0.6.0",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("ENmix",
##                     version = "1.26.10",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("profmem",
##                     version = "0.6.0",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("profvis",
##                     version = "0.3.7",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("htmlwidgets",
##                     version = "1.5.4",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("ungeviz",
##                     version = "0.1.0",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("AnnotationDbi",
##                     version = "1.52.0",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("Biobase",
##                     version = "2.50.0",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("BiocGenerics",
##                     version = "0.36.1",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("BiocParallel",
##                     version = "1.24.1",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("ChAMPdata",
##                     version = "2.22.0",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("CopyNumber450kCancer",
##                     version = "1.0.6",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("DelayedArray",
##                     version = "0.16.3",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("GenomicRanges",
##                     version = "1.42.0",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("Matrix",
##                     version = "1.3-2",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("RaggedExperiment",
##                     version = "1.14.2",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("S4Vectors",
##                     version = "0.28.1",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("cowplot",
##                     version = "1.1.1",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("dendextend",
##                     version = "1.15.1",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("doParallel",
##                     version = "1.0.16",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("future",
##                     version = "1.21.0",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("gtools",
##                     version = "3.8.2",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("illuminaio",
##                     version = "0.32.0",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("iterators",
##                     version = "1.0.13",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("matter",
##                     version = "1.16.0",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("pheatmap",
##                     version = "1.0.12",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("ramify",
##                     version = "0.3.3",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("readxl",
##                     version = "1.3.1",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("rlist",
##                     version = "0.4.6.2",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
##
##BiocManager::install("scales",
##                     version = "1.1.1",
##                     lib="C:\\Users\\Mike\\Desktop\\methyl_master_r_lib",
##                     update=FALSE)
