#!/usr/bin env Rscript

#' @title methyl_master_load_packages
#' @description A function that attempts to load all packages for MethylMaster
#' Note loading packages is borken up so that when the
## pipeline is run on windows it won't overload RStudio
## and spontaneously abort
#' @param routine The routine to be used
#' @return Loads all pertinent R packages
#' @export

methyl_master_load_packages <- function(routine
                                        ){

if(routine=="test"){

  library(magrittr)
  library(profvis)
  library(profmem)

}else if(routine=="champ"){

  ##Just load champ libraries when running champ
  ##so wont crash on windows
  ##For ChAMP:
  ##library(sesame)
  library(dplyr)
  library(magrittr)
  library(data.table)
  library(ggplot2)
  library(sesameData)
  library(stats) ##Model.matrix()
  library(ChAMP)
  library(Epicopy)
  library(limma) ##lmFit() , contrasts.fit(), eBayes(), topTable()
  library(igraph) ##graph.adjency(), clusters(), set.edge.attribute(), V(), E()
  library(org.Hs.eg.db) ##Contains org.Hs.egSYMBOL
  library(marray) ##mPalette()
  library(shape) ##colorlegend()
  library(CNVRanger)
  library(profmem)
  library(profvis)
  library(htmlwidgets)
  library(ungeviz)

}else if(routine=="epicopy"){

  library(dplyr)
  library(magrittr)
  library(data.table)
  library(ggplot2)
  library(minfi)
  library(epicopy)
  library(profmem)
  library(profvis)
  library(htmlwidgets)
  library(ungeviz)

}else{

  library(foreach)
  library(ggplot2)
  library(magrittr)
  library(plyr)
  library(dplyr)
  library(foreach)
  library(data.table)
  library(TCGAbiolinks)
  library(ExperimentHub)
  library(sesame)
  library(sesameData)
  library(minfi)
  library(GEOquery)
  library(DNAcopy)
  library(CNVRanger)
  library(cnAnalysis450k)
  library(CNAclinic)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  library(FlowSorted.CordBloodCombined.450k)
  ##library(ENmix)

  ##For ChAMPP
  library(stats) ##Model.matrix()
  library(ChAMP)
  library(Epicopy)
  library(limma) ##lmFit() , contrasts.fit(), eBayes(), topTable()
  library(igraph) ##graph.adjency(), clusters(), set.edge.attribute(), V(), E()
  library(org.Hs.eg.db) ##Contains org.Hs.egSYMBOL
  library(marray) ##mPalette()
  library(shape) ##colorlegend()

  library(profmem)
  library(profvis)
  library(htmlwidgets)
  library(ungeviz)

}

}
