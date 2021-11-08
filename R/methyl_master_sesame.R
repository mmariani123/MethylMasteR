#!/usr/bin/env Rscript

#' @title ##Methyl MasteR Sesame
#' @description MethylMasteR SeSAMe run
#'
#' First part:: Very similar to PROCESS SESAME
#' command to prepare for segementation
#' this is included so that the total
#' time and mem for running sesame
#' can be recorded
#'
#' @param sesame.idat.files.dir
#' @param sesame.samples
#' @import data.table
#' @import dplyr
#' @import sesame
#' @import sesameData
#' @export

methyl_master_sesame <- function(sesame.idat.files.dir=getwd(),
                                 sesame.samples=NULL,
                                 sesame.platform="HM450", ##"EPIC" for cord
                                 sesame.data.cache="EPIC",
                                 sesame.data.normal="EPIC.5.normal",
                                 sesame.ref.version="hg38",
                                 sesame.calc.betas=FALSE,
                                 sesame.infer.sex.karyotypes=FALSE,
                                 sesamle.ref=NULL,
                                 sesame.output.dir,
                                 sesame.reference="internal", ##"Sample_Group"
                                 sesame.reference.names=c("tumor","normal"),
                                 split.by=NULL,
                                 ...){

if(sesame.reference=="internal"){
  ##Get sesame normal samples:
  sesameData::sesameDataCache(sesame.data.cache)
  sesame_ssets_normal <- sesameData::sesameDataGet(sesame.data.normal)
  setExperimentHubOption("CACHE", idat.files.dir)
  ExperimentHub()
  idat_prefixes <- searchIDATprefixes(idat.files.dir,recursive=TRUE)
  sesameDataCacheAll()

  sesame_betas <- openSesame(idat_prefixes,
                             mask = TRUE,
                             sum.TypeI = TRUE,
                             platform = sesame.platform)

  sesame_sset <- openSesame(idat_prefixes,
                            mask = TRUE,
                            sum.TypeI = TRUE,
                            platform = sesame.platform,
                            what="sigset")

  sesame_rgset <- sesame::SigSetsToRGChannelSet(sesame_sset)

  sesame_karyotype <- foreach(i = 1:length(names(sesame_sset))) %do% {
    sesame::inferSexKaryotypes(sesame_sset[[i]])
  }
  names(sesame_karyotype) <- names(sesame_sset)
  sesame_karyotype <- as.data.frame(unlist(sesame_karyotype))
  colnames(sesame_karyotype) <- "karyotype"

  sesame_seg <- foreach(i = 1:length(names(sesame_sset))) %do% {
     sesame::cnSegmentation(sesame_sset[[i]],
                            sesame_ssets_normal,
                            refversion = sesame.ref.version)
  }
  names(sesame_seg) <- names(sesame_sset)

  sesame_qc <- foreach(i = 1:length(names(sesame_sset))) %do% {
    sesame::sesameQC(sesame_sset[[i]])
  }
  names(sesame_qc) <- names(sesame_sset)
  sesame_qc <- as.data.frame(data.table::rbindlist(sesame_qc))
  rownames(sesame_qc) <- names(sesame_sset)

}else if(sesame.reference!="internal"){
  ref.subs <- list()
  for(i in seq_along(sesame.reference.names)){
  ref.subs[[i]] <- data.table::subset(sesame.samples,
                               {{sesame_reference}}==sesame.reference.names[i])
  ref.subs[[i]]$Sample_ID
  setExperimentHubOption("CACHE",
                         idat.files.dir)
  ExperimentHub()
  idat_prefixes[[i]] <- searchIDATprefixes(idat.files.dir,
                                      recursive=TRUE)

  sesameDataCacheAll()
  if(sesame.calc.betas==TRUE){
    sesame_betas[[i]] <- openSesame(idat_prefixes[[i]],
                               mask = TRUE,
                               sum.TypeI = TRUE,
                               platform = sesame.platform)
  }
  sesame_sset[[i]] <- openSesame(idat_prefixes[[i]],
                            mask = TRUE,
                            sum.TypeI = TRUE,
                            platform = sesame.platform,
                            what="sigset")

  if(sesame.infer.sex.karyotypes==TRUE){
    sesame.karyotype <- foreach(i = 1:length(names(sesame_sset))) %do% {
      sesame::inferSexKaryotypes(sesame_sset[[i]])
    }
    names(sesame.karyotype) <- names(sesame_sset)
    sesame.karyotype <- as.data.frame(unlist(sesame_karyotype))
    colnames(sesame.karyotype) <- "karyotype"
    data.table::fwrite(sesame.karyotype,
                       file = sesame.output.dir,
                       col.names = TRUE,
                       row.names = FALSE,
                       quote = FALSE)
  }

  sesame_qc[[i]] <- foreach(i = 1:length(names(sesame_sset))) %do% {
    sesame::sesameQC(sesame_sset[[i]])
  }
  names(sesame_qc) <- names(sesame_sset)
  sesame_qc <- as.data.frame(data.table::rbindlist(sesame_qc))
  rownames(sesame_qc) <- names(sesame_sset)

}

  sesame_seg <- foreach(i = 1:length(names(sesame_sset))) %do% {
    sesame::cnSegmentation(sesame_sset[[i]],
                           ses.ref,
                           refversion = sesame.ref.version)
  }
  names(sesame_seg) <- names(sesame_sset)

########## separate baseline_seg into normal and tumor #################

sesame_seg  <- sesame_seg[
  names(sesame_seg) %in%
    tumor$X[tumor[[gender_reported]]=={{split.by}}]]
for(i in 1:length(sesame_seg_normal_tumor_female)){
  sample.now <- names(sesame_seg_normal_tumor_female)[i]
  sesame_seg[[i]]$seg.signals$karyotype    <- ""
  sesame_seg[[i]]$seg.signals$sex_reported <- ""
  sesame_seg[[i]]$seg.signals$sex_inferred <- ""
}

##Note the below are actual normal samples
save(sesame_seg_normal_normal,
     file=paste0(work.dir,
                 file.sep,
                 "sesame_seg_normal_normal.RData"))
save(sesame_seg_normal_tumor,
     file=paste0(work.dir,
                 file.sep,
                 "sesame_seg_normal_tumor.RData"))

#################### PROCESS CORD #############################

## Subset cord blood into female and male samples

         karyotype_cord <- NA

         karyotype_cord <- foreach(i = 1:length(names(sesame_sset_cord))) %do%
           {sesame::inferSexKaryotypes(sesame_sset_cord[[i]])}

         sesame_female_cord_id <- character(length=0)
         sesame_male_cord_id   <- character(length=0)
         ##sesame_female_cord_id
         ##sesame_male_cord_id

         ##length(karyotype.cord) ##20
         for(i in 1:length(karyotype_cord)){
           print(karyotype_cord[i])
           if(as.character(karyotype_cord[i]) == "XaXi"){
             sesame_female_cord_id <- append(sesame_female_cord_id,
                                             names(sesame_sset_cord)[i])
           }else{
             sesame_male_cord_id <- append(sesame_male_cord_id,
                                           names(sesame_sset_cord)[i])
           }
         }

}
