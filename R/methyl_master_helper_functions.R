#!/usr/bin env Rscript

##Michael Mariani PhD Dartmouth College 2021-2022

##Helper functions for MethylMaster

######################## FUNCTIONS ########################################

#' @title harmonized_rg_sets
#' @description Harmonize RG data
#' Functionality by Lucas Salas PhD Dartmouth College 2021
#' @param rgset1 Path to the first RGset
#' @param rgset2 Path to the second RGset
#' @return Harmonized data set
#' @export
harmonized_rg_sets <- function(rgset1,rgset2){
  rgset.harm <- minfi:::.harmonizeDataFrames(minfi:::.pDataFix(colData(rgset1)),
                                             minfi:::.pDataFix(colData(rgset2)))
  return(rgset.harm)
}

#' @title infer_sex_karyotype
#' @description Use sesame to infer sex karyotype
#' functionality by Jennifer Chen Dartmouth College 2021
#' @param sset SeSAMe sset
#' @param ref.sset SeSAMe reference sset
#' @return sset with inferred sex data
#' @export
infer_sex_karyotype <- function(sset, ref.sset){
  karyotype.inferred <- foreach(i = 1:length(names(ref.sset))) %do%
    {
      sesame::inferSexKaryotypes(sset[[i]])
    }
}

#' @title weightedmean
#' @description Weighted-mean function
#' Original function from Epicopy by Cho et al., 2019
#' @param scores  input scores
#' @param ranges  input Granges
#' @param qranges input Qranges
#' @return weighted mean values
#' @export
weightedmean <- function(scores, ranges, qranges){
    return(sum(scores * width(ranges)) / sum(width(ranges)))
}

#' @title weightedmean.epicopy
#' @description Epicopy-style weighted mean function
#' Original function from Epicopy by Cho et al., 2019
#' @param scores  input scores
#' @param ranges  input Granges
#' @param qranges input Qranges
#' @return epicopy-style weighted mean values
#' @export
weightedmean.epicopy <- function(scores, ranges, qranges){
    isects <- pintersect(ranges, qranges)
    return(sum(scores * width(isects)) / sum(width(isects)))
}

#' @title binding_frames_mm
#' @description Bind individual sesame frames together
#' @param x Wesame-style data frame
#' @param auto.corrected Whether or not data went through AutoCorrecPeaks
#' @param add.col Addtional (numeric) column to add to output dataframe to be
#' used in downstream analysis/plotting
#' @return Segment dataframe output of multiple dfs bound together
#' @export
binding_frames_mm <- function(x,
                              auto.corrected,
                              add.col
                              ){
  binding.list <- list()
  for(i in 1:length(names(x))){
    if(auto.corrected==TRUE){
      binding.list[[i]] <- x[[i]]
      binding.list[[i]]$ID <- names(x)[i]
      if(!is.null(add.col)){
        binding.list[[i]][[names(add.col)[i]]] <- unname(add.col[i])
      }
      if(any(is.na(binding.list[[i]]$Chromosome))){
        stop(paste0("ERROR one or more entries in the 'chrom' ",
                    "field in the sesame seg output are NA"))
      }
    }else{
      binding.list[[i]] <- x[[i]]$seg.signals
      binding.list[[i]]$ID <- names(x)[i]
      if(!is.null(add.col)){
        binding.list[[i]][[names(add.col)[i]]] <- unname(add.col[i])
      }
      if(any(is.na(binding.list[[i]]$chrom))){
        stop(paste0("ERROR one or more entries in the 'chrom' ",
        "field in the sesame seg output are NA"))
      }
    }
  }
  seg.out <- do.call(rbind,binding.list)
  return(seg.out)
}

#' @title calc_seg_state
#' @description common function used to calculate seg state from seg mean
#' @param seg.means vector of seg.means
#' @param upper.thresh CN state threshold, any number above this will be
#' assigned this value
#' @param cutoff vector of cutoff threshold for loss and
#' gain respectively
#' @return vector of calculated CNV states
#' @export
calc_seg_state <- function(seg.means=NULL,
                           upper.thresh=4,
                           cutoff=NULL
                           ){
  if(is.null(cutoff)){
    seg.state <- round(2^seg.means * 2)
    seg.state[seg.state > upper.thresh] <- upper.thresh
  }else{
    seg.state <-
      ifelse(seg.means>=cutoff[2],3,
             ##seg.means>=cutoff[2],2+floor(seg.means/0.3),
             ifelse(seg.means<=cutoff[1],1,2))
    seg.state[seg.state > upper.thresh] <- upper.thresh
  }
  return(seg.state)
}

#' @title fast_450_k_binning
#' @description A faster version for binding 450K bins
#' Originally function by Knoll et al., 2017
#' @param proc.samples The samples used as input
#' @param med.ref.values The median reference values
#' @param proc.ref.samples the refernce samples
#' @param bin.size The bin size
#' @return A 450k-style candidates data frame
#' @export
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

#' @title AutoCorrectPeak.mm
#' @description MethylMaster version of CopyNumber450kCancer::AutoCorrectPeak
#' function. Original function by Marzouka et al. 2016
#' @param object The sesame CNV object to have its peaks adjusted/corrected
#' @param output.dir The specified output directory
#' @param cutoff The cutoff threshold
#' @param markers The number of markers
#' @param ... Additional parameters to AutoCorrectPeak.mm
#' @import CopyNumber450kCancer
#' @return Corrected regions plots and files and returns a corrected object
#' @export
AutoCorrectPeak.mm <- function(object,
                               output.dir=getwd(),
                               cutoff = 0.1,
                               markers = 20,
                               ...){

  if(exists(output.dir)){

    unlink(output.dir, recursive = TRUE)
    dir.create(output.dir)

  }else{

    dir.create(output.dir)

  }

  ##Set NAs to 0
  if(missing(markers)){
    markers <- c(0)
  }

  QC <- object$SL[, 1:2]

  QC[, 3:7] <- 0

  colnames(QC) <- c("Sample",
                    "Comment",
                    "peak.sharpness",
                    "number.of.regions",
                    "IQR",
                    "SD",
                    "MAPD")

  par(mfrow = c(2, 2),
      mar = c(5.1, 0, 4.1, 0),
      oma = c(2, 0, 0, 4))

  layout(matrix(c(1, 2, 3, 4), 2, 2,
                byrow = TRUE),
         widths = c(3, 21),
         heights = c(10, 10),
         TRUE)

  ##for(i in 1:length(object$SL[, "Sample"])){

    ##print(paste("Auto Correction.... Sample number ",i))
    print(paste("Auto Correction.... Sample number ",1))

    sam <- object$regions_auto[which(object$regions_auto$Sample %in%
                                  ##as.character(object$SL[i, "Sample"])), ]
                                    as.character(object$SL[1, "Sample"])), ]

    forDen <- sam[which(sam$Chromosome != "chrX" &
                          sam$Chromosome != "chrY"), c("Num.of.Markers",
                                                       "Mean")]

    sam.original <- sam

    d <- density(forDen$Mean,
                 weights = forDen$Num.of.Markers/sum(forDen$Num.of.Markers),
                 na.rm = TRUE,
                 kernel = c("gaussian"),
                 adjust = 0.15,
                 n = 512)

    max.peak.value <- d$x[which.max(d$y)]

    sam$Mean <- sam$Mean - max.peak.value

    point <- which(d$x == (max.peak.value))

    QC.peak.sharpness <-
      ((d$y[point + 20] + d$y[point - 20])/2)/d$y[which(d$x == max.peak.value)]

    ##QC[i, "peak.sharpness"] <- as.numeric(QC.peak.sharpness)
    QC[1, "peak.sharpness"] <- as.numeric(QC.peak.sharpness)

    ##QC[i, "number.of.regions"] <- length(sam[, 1])
    QC[1, "number.of.regions"] <- length(sam[, 1])

    ##QC[i, "IQR"] <- IQR(forDen[, "Mean"], na.rm = TRUE, type = 7)
    QC[1, "IQR"] <- IQR(forDen[, "Mean"], na.rm = TRUE, type = 7)

    ##QC[i, "SD"] <- sd(forDen[, "Mean"], na.rm = TRUE)
    QC[1, "SD"] <- sd(forDen[, "Mean"], na.rm = TRUE)

    ##QC[i, "MAPD"] <- median(abs(diff(forDen[, "Mean"],
      ##                               na.rm = TRUE)), na.rm = TRUE)
    QC[1, "MAPD"] <- median(abs(diff(forDen[, "Mean"],
                                   na.rm = TRUE)), na.rm = TRUE)

    ##object$regions_auto[which(object$regions_auto$Sample %in%
    ##              as.character(object$SL[i, "Sample"])), "Mean"] <- sam$Mean
    object$regions_auto[which(object$regions_auto$Sample %in%
                          as.character(object$SL[1, "Sample"])), "Mean"] <-
      sam$Mean

  ##object$mod_auto[which(object$mod_auto$Sample %in% as.character(object$SL[i,
  ##              "Sample"])), 2:4] <- c(round(max.peak.value,
  ##                                    3), round(-max.peak.value, 3), "Auto")
    object$mod_auto[which(object$mod_auto$Sample %in% as.character(object$SL[1,
                            "Sample"])), 2:4] <- c(round(max.peak.value,
                                  3), round(-max.peak.value, 3), "Auto")

    ##print(paste("Plotting.... Sample number ", i))
    print(paste("Plotting.... Sample number "))

    pdf(file = paste0(output.dir,
                      .Platform$file.sep,
                      ##i,
                      ##"_",
                      ##object$SL[i, "Sample"],
                      ##object$SL[1, "Sample"],
                      gsub(".*/(?!.*/)","",output.dir,perl = TRUE),
                      "_",
                      "auto_corrected_plot.pdf"),
                      width = 12,
                      height = 8)

    par(mfrow = c(2, 2), mar = c(5.1, 0, 4.1, 0), oma = c(2, 0, 0, 4))

    layout(matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE), widths = c(3, 21),
           heights = c(10, 10), TRUE)

    plot(d$y, d$x, ylim = c(-1, 1), type = "l", ylab = "",
         xlab = "", axes = FALSE, xlim = rev(range(d$y)))

    abline(h = c(0, -cutoff, cutoff), lty = 3)

    box()

    legend("bottomleft", legend = "Density", cex = 1)

    plotRegions(sam.original,
                cutoff = cutoff,
                markers = markers,
                main = paste("Sample::",
                             ##object$SL[i, "Sample"],
                             ##object$SL[1, "Sample"],
                             gsub(".*/(?!.*/)","",output.dir,perl = TRUE)
                             ##"       Info::",
                             ##object$SL[i, "Comment"]
                             ##object$SL[1, "Comment"]
                             ),
                ...)

    plot(d$y, d$x - max.peak.value, ylim = c(-1, 1), type = "l",
         ylab = "", xlab = "", axes = FALSE, xlim = rev(range(d$y)))

    abline(h = c(0, -cutoff, cutoff), lty = 3)

    box()

    legend("bottomleft", legend = "Density", cex = 1)

    plotRegions(sam,
                cutoff = cutoff,
                markers = markers,
                main = paste("Autocorrected plot"))

    dev.off()

  ##}

  object$QC <- QC

  object$regions <- object$regions_auto

  print(paste0("Saving the autocorrection files to ",
               output.dir,
               " ..."))

  write.csv(object$regions_auto,
            file = paste0(output.dir,
                          .Platform$file.sep,
                          "autocorrected_regions.csv"))

  write.csv(object$mod_auto,
            file = paste0(output.dir,
                          .Platform$file.sep,
                          "autocorrections.csv"))

  write.csv(QC,
            file = paste0(output.dir,
                          .Platform$file.sep,
                          "QC.csv"),
            row.names = FALSE)

  print("Auto Correction Done.")

  return(object$regions)

}
