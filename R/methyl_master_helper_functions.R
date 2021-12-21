#!/usr/bin env Rscript

##Michael Mariani Dartmouth College 2021

##Helper functions for MultiMethyl

######################## FUNCTIONS ########################################

#' Harmonize RG data
#'
#'
#' @param rgset1 Path to the input file
#' @param rgset2
#' @return Harmonized data set
#' @export
harmonized_rg_sets <- function(rgset1,rgset2){
  rgset.harm <- minfi:::.harmonizeDataFrames(minfi:::.pDataFix(colData(rgset1)),
                                             minfi:::.pDataFix(colData(rgset2)))
  return(rgset.harm)
}

#' Use sesame to infer sex karyotype
#'
#'
#' @param sset Path to the input file
#' @param ref.sset
#' @return sset with inferred sex data
#' @export
infer_sex_karyotype <- function(sset, ref.sset){
  karyotype.inferred <- foreach(i = 1:length(names(ref.sset))) %do%
    {
      sesame::inferSexKaryotypes(sset[[i]])
    }
}

#' Weighted-mean function
#'
#'
#' @param scores  input scores
#' @param ranges  input Granges
#' @param qranges input Qranges
#' @return weighted mean values
#' @export
weightedmean <- function(scores, ranges, qranges){
    return(sum(scores * width(ranges)) / sum(width(ranges)))
}

#' Epicopy-style weighted mean funtion
#'
#'
#' @param scores  input scores
#' @param ranges  input Granges
#' @param qranges input Qranges
#' @return epicopy-style weighted mean values
#' @export
weightedmean.epicopy <- function(scores, ranges, qranges){
    isects <- pintersect(ranges, qranges)
    return(sum(scores * width(isects)) / sum(width(isects)))
}

#' Bind individual sesame frames together
#'
#'
#' @param x sesame-style data frame
#' @return segment dataframe output of multiple dfs bound together
#' @export
binding_frames_mm <- function(x){
  binding.list <- list()
  for(i in 1:length(names(x))){
    binding.list[[i]] <- x[[i]]$seg.signals
    binding.list[[i]]$ID <- names(x)[i]
  }
  seg.out <- do.call(rbind,binding.list)
  return(seg.out)
}

#' Bind individual sesame frames together
#'
#'
#' @param x sesame-style data frame
#' @param cutoff cutoff CN, any number above this will be assigned this value
#' @return seg frame with state assigned
#' @export
calc_seg_state <- function(seg, cutoff){
  seg$state <- round(2^seg$seg.mean * 2)
  seg$state[seg$state > cutoff] <- cutoff
  return(seg)
}

#' This is a faster version of binding 450K bins
#'
#'
#' @param proc.samples samples input
#' @param med.ref.values median reference values
#' @param proc.ref.samples refernce samples
#' @param bin.size bin size
#' @return 450k-style candidates data frame
#' @export
########################  ########################

fast_450k_binning <- function(proc.samples,
                              med.ref.values,
                              proc.ref.samples,
                              bin.size = 500000){

  ## calculate bins, binsize=50000
  ## CHANGE BINSIZE HERE

  ##candidates_data_normal_female <-
  ##  cnAnalysis450k::createBinsFast(proc_tumor_female[, 1:3, drop = FALSE],
  ##                                 med_normal_female,
  ##                                 proc_normal_female,
  ##                                 binsize = 500000)

  candidates.data <-
    cnAnalysis450k::createBinsFast(proc.samples[, 1:3, drop = FALSE],
                                   med.ref.values,
                                   proc.ref.samples,
                                   binsize = bin.size)

  return(candidates.data)
  ## Note the offset function in the above script

}
