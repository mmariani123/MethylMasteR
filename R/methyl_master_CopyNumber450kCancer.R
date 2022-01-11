#!/usr/bin/env Rscript

#################### CopyNumber450kCancer ###################################
############### Note the offset function ####################################
#############################################################################
#############################################################################
#############################################################################

#' @title ReadData
#' @description MethylMaster version of the ReadData()
#' function from CopyNumber450kCancer.R
#'
#' 'Function to read the data------
#' Functions for checking the headers:::
#' this function to check the header of the input file:
#' ("Sample" "Chromosome" "bp.Start"  "bp.End" "Num.of.Markers" "Mean")'
#'
#' @param regions_file The regions file
#' @param Sample_list The Sample_list file
#' @param copynumber450k The copynumber450 parameter
#' @return CopyNumber450kCancer object
#' @export
ReadData <- function(regions_file,
                     Sample_list,
                     copynumber450k=FALSE
                     ){

    checkHeaderRegions <- function(file){
        colnamesFile<-c("Sample",
                        "Chromosome",
                        "bp.Start",
                        "bp.End",
                        "Num.of.Markers",
                        "Mean")
        if(sum(colnames(file)==colnamesFile)==6){
            print("The header of the regions is ... OK")
        }
        else{
            print(paste0("The header of the file is not OK...should be: ",
                    "Sample/Chromosome/bp.Start/bp.End/Num.of.Markers/Mean"))
        }
        invisible(sum(colnames(file)==colnamesFile))
        ##if 6 then it is OK ekse there is problem in the header
    }

    ##this function to check the Sample info
    ##The file of the samples info shoud have the header ("Sample","Comment")
    checkHeaderSL <- function(file){
        colnamesSL<-c("Number","Sample","Comment")
        if(sum(colnames(file)==colnamesSL)==3){
            print("The header of the samples list is ... OK")
        }
        else{
          print(paste0("The header of the sample file is not OK...should be: ",
                         "Number / Sample / Comment"))
        }
        invisible(sum(colnames(file)==colnamesSL))
        ##if 3 then it is OK if else there is problem in the header
    }
    if(is.character(regions_file)){
        regions<-read.csv(regions_file,stringsAsFactors =FALSE)
        ##load the file that contains regions and means
    } else {
        regions <- regions_file
    }
    if (copynumber450k==TRUE){
        regions<-regions[,c("Sample",
                            "chrom",
                            "loc.start",
                            "loc.end",
                            "num.mark",
                            "seg.mean")]

        colnames(regions)<-c("Sample",
                             "Chromosome",
                             "bp.Start",
                             "bp.End",
                             "Num.of.Markers",
                             "Mean")

    }

    if(is.character(regions_file)){
        SL<-read.csv(Sample_list,stringsAsFactors =FALSE)
        ##Load the file that contains the names of samples and the comments
    } else {
        SL <- Sample_list
    }

    checkHeaderRegions(regions)
    checkHeaderSL(SL)

    ##regions[is.na(regions)] <- 0
    object <- list(
        mainDir = getwd(),
        regions = regions,
        regions_save = regions,
        regions_auto = regions,
        SL = SL)
    class(object) <- "CopyNumber450kCancer_data"
    mod<-as.data.frame(SL$Sample,stringsAsFactors =FALSE)
    ## copy to store the modification in it
    mod[,2:6]<-0
    mod[,6]<-"No"
    mod[is.na(mod)] <- 0
    colnames(mod)<-c("Sample","Lower_selected_level",
                     "Upper_selected_level",
                     "Mean_of_selected_regions",
                     "Shifting",
                     "Reviewed?")
    object$mod <- mod

  mod_auto<-mod[,c("Sample","Mean_of_selected_regions","Shifting","Reviewed?")]
    colnames(mod_auto)<-c("Sample",
                          "Auto_Maximum_Peak",
                          "Shifting",
                          "Auto_Corrected?")
    object$mod_auto <- mod_auto

    object
}

#' @title print.CopyNumber450kCancer_data
#' @description MethylMaster version of the print.CopyNumber450kCancer_data()
#' function from CopyNumber450kCancer.R
#'
#' "print.CopyNumber450kCancer_data
#'
#' Do the printing"
#'
#' @param x The x parameter
#' @param ... Additional parameters passed to print.CopyNumber450kCancer
#' @return
#' @export
print.CopyNumber450kCancer_data <- function(x, ...){
    cat(sprintf("CopyNumber450kCancer data with %i samples in regions and %i
              samples in sample list.\n",
                length(unique(x$regions$Sample)),
                length(unique(x$SL$Sample))
    ))
    cat("Contains:\n", paste("  ", names(x), collapse = "\n"), "\n", sep="")
}

#' @title CopyNumber450kCancer
#' @description The MethylMaster version of the CopyNumber450kCancer() function
#' from "CopyNumber450kCancer.R
#'
#' functions for plotting:::
#' ---to plot-------------
#' this uses the regions file: Chromosome should be in this format: "chr1"
#' similar to the original function, this one use only the cutoff
#' This function plots the chromosomal regions (segments) with colored segments
#' based on the cutoff. This function was built based on "plotSample" function
#' in "CopyNumber450k" package (http://www.bioconductor.org/packages/release/
#' bioc/html/CopyNumber450k.html), and uses a modified "minor.tick" function in
#' Hmisc" package to draw small tick in the plots."
#'
#' @param object The object parameter
#' @param chr The chr parameter
#' @param start The start parameter
#' @param end The end parameter
#' @param cutoff The cutoff parameter
#' @param markers The markers parameter
#' @param ... Additional parameters passed to opyNumber450kCancer
#' @return Return plot
#' @export
plotRegions <- function(object,
                        chr,
                        start,
                        end,
                        cutoff=0.1,
                        markers=20,
                        ...
                        ){
    sample_segments <- object

    if(hasArg(markers)){
        sample_segments$Mean[which(sample_segments$Num.of.Markers<=markers)]<-0
    }

    segment_values <- as.numeric(sample_segments[,"Mean"])
    segment_colors <- rep("black", nrow(sample_segments))

    if (missing(cutoff)) {
        cutoff<-(0.1)
    }

    segment_colors[as.numeric(segment_values) >= cutoff] <- "green"
    segment_colors[as.numeric(segment_values) <= -cutoff] <- "red"

    if (missing(chr)) {
        # Plotting the whole genome
        chromosomes <- unique(sample_segments[, "Chromosome"])
        site_per_chr <- cumsum(c(0,
                                 sapply(
                                     chromosomes, function(chr) max(as.numeric(
          sample_segments[sample_segments[,"Chromosome"] == chr, "bp.End"])))))
        offset <- site_per_chr - min(
            as.numeric(
      sample_segments[sample_segments[, "Chromosome"] == "chr1", "bp.Start"]))
        ## 1 instead of "chr1" #as.numeric(gsub("\\D", "", x))
        start <- 0
        end <- as.numeric(max(site_per_chr))
        x_axis_type <- "n"
    } else {
        # Plotting a region
        if (missing(start)) {
            start <- 0
        }

        if (missing(end)) {
            end <- as.numeric(
      max(sample_segments[sample_segments[, "Chromosome"] == chr, "bp.End"]))
        }

        chromosomes <- chr
        offset <- 0
        x_axis_type <- NULL
    }

    yMin <- (-1) ##min(c(-1, as.numeric(sample_segments[significant_segments,
    ##"Mean"])))
    yMax <- 1 ##max(c(1, as.numeric(sample_segments[significant_segments,
    ##"Mean"])))

    ##if (missing(ylab)) {
    # #ylab<-""
    ##}

    myPlot <- plot(range(start, end),
                   range(yMin, yMax),
                   type = "n",
                   axes=FALSE,
                   xaxt = x_axis_type,
                   xlab="",
                   ylab="", ...) ##ylab="Mean",

    ##---this function to plot the tick on the right side
    tick.tick<-function (nx = 2, ny = 2, tick.ratio = 0.5) {
        ax <- function(w, n, tick.ratio) {
            range <- par("usr")[if (w == "x")
                1:2
                else 3:4]
            tick.pos <- if (w == "x")
                par("xaxp")
            else par("yaxp")
            distance.between.minor <- (tick.pos[2] - tick.pos[1])/tick.pos[3]/n
            possible.minors <- tick.pos[1] - (0:100) * distance.between.minor
            low.minor <- min(possible.minors[possible.minors >= range[1]])
            if (is.na(low.minor))
                low.minor <- tick.pos[1]
            possible.minors <- tick.pos[2] + (0:100) * distance.between.minor
            hi.minor <- max(possible.minors[possible.minors <= range[2]])
            if (is.na(hi.minor))
                hi.minor <- tick.pos[2]
            axis(if (w == "x")
                1
                else 4, seq(low.minor, hi.minor, by = distance.between.minor),
                labels = FALSE, tcl = par("tcl") * tick.ratio)
        }
        if (nx > 1)
            ax("x", nx, tick.ratio = tick.ratio)
        if (ny > 1)
            ax("y", ny, tick.ratio = tick.ratio)
        invisible()
    }

    if (missing(chr)) {
        xlabs <- sapply(2:length(site_per_chr), function(j) {
            ((site_per_chr[j] - site_per_chr[j - 1])/2) + site_per_chr[j - 1]
        })

        axis(1, at = xlabs, labels = chromosomes, lty = 0, las = 2, ...)
        axis(4)
        tick.tick(nx=0,ny=2, tick.ratio=1.6)
        tick.tick(nx=0,ny=10, tick.ratio=0.6)
        mtext("L-value", side = 4, line = 2, cex = par("cex.lab"))
        box()
        abline(v = site_per_chr, lty = 3)
        abline(h = c(0,-cutoff,cutoff), lty = 3)
    }

    lapply(1:length(chromosomes), function(i) {
        used_segments <- sample_segments[, "Chromosome"] == chromosomes[i]
        colors <- segment_colors[used_segments]
  starts <- as.numeric(sample_segments[used_segments, "bp.Start"]) + offset[i]
        ends <- as.numeric(sample_segments[used_segments, "bp.End"]) + offset[i]
        y <- as.numeric(sample_segments[used_segments, "Mean"])
        graphics::segments(starts, y, ends, y, col = colors, lwd = 2, lty = 1)
    })

    return(myPlot)
}
