#!/usr/bin/env Rscript

#' @title methyl_master_plot_individual
#' @description output signal inensity plots for individual samples
#' ##Michael Mariani PhD Dartmouth College 2021
#' @import ggplot2
#' @param pi.seg
#' @param pi.output.dir
#' @param pi.name
#' @param ...
#' @return #NULL
#' @export
methyl_master_plot_individual <- function(pi.seg = NULL,
                           pi.output.dir = NULL,
                           pi.name = NULL,
                           ...
                           ){

  individual.plots.dir <- paste0(pi.output.dir,
                                 .Platform$file.sep,
                                 "indivdual_plots_",
                                  pi.name)

  if(!dir.exists(individual.plots.dir)){

    dir.create(individual.plots.dir)

  }else{

    print("Individual plots dir exists, overwriting")
    unlink(individual.plots.dir,
           recursive=TRUE,
           force=TRUE)
    dir.create(individual.plots.dir)

  }

  ##For others willl need to convert to sesame object or recode graph
  ##visualize individual sesame by sample:
  debug(visualizeSegments.mm)
  pi.seg <- readRDS(file = paste0("C:\\Users\\Mike\\Desktop\\pi.seg.RDS"))
  samples <- unique(pi.seg$Sample_ID)
  sample.now <- pi.seg[pi.seg$Sample_ID==samples[1],]
  visualizeSegments.mm(sample.now)

  for(i in 1:length(samples)){

    sample.now <- pi.seg[pi.seg$Sample_ID==samples[i],]
    ##ggplot2::ggsave({sesame::visualizeSegments(pi.seg[[i]]) +
    ggplot2::ggsave({visualizeSegments.mm(sample.now) +
        ggplot2::theme_bw() +
        ggplot2::ylab("Signal") +
        ggplot2::labs(color="Sesame seg. \nsignal") +
        ggplot2::ggtitle(paste0("Sample ",
                       sample.now,
                       "\n",
                       pi.name)) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90),
              plot.title = ggplot2::element_text(hjust=0.5))},
        filename = paste0(individual.plots.dir,
                          .Platform$file.sep,
                          sample.now,
                          "_",
                          pi.name,
                          "_sesame_plot.pdf"),
        width=12,
        height=8,
        device="pdf")

  }

}

#' @title visualizeSegments.mm
#' @description modified sesame::visualizeSegments function
#' Michael Mariani PhD Dartmouth College 2021, originally by SeSAMe team
#' @param seg
#' @param to.plot
#' @import ggplot2
#' @import scales
#' @import GenomicRanges
#' @return #A ggplot
#' @export
visualizeSegments.mm <- function(seg,
                                 to.plot = NULL
                                 ){

  ##stopifnot(is(seg, "CNSegment"))

  ##pkgload::pkgTest("ggplot2")
  ##pkgload::pkgTest("scales")
  ##pkgload::pkgTest("GenomicRanges")

  ##bin.coords  <- seg$bin.coords
  ##bin.seqinfo <- GenomicRanges::seqinfo(bin.coords)
  ##bin.signals <- seg$bin.signals
  ##seg.signals <- seg$seg.signals
  ##total.length <- sum(as.numeric(bin.seqinfo@seqlengths), na.rm = TRUE)
  ##if (is.null(to.plot)){
  ##  to.plot <- (bin.seqinfo@seqlengths > total.length * 0.01)
  ##}
  ##seqlen <- as.numeric(bin.seqinfo@seqlengths[to.plot])
  ##seq.names <- bin.seqinfo@seqnames[to.plot]
  ##total.length <- sum(seqlen, na.rm = TRUE)
  ##seqcumlen <- cumsum(seqlen)
  ##seqcumlen <- seqcumlen[-length(seqcumlen)]
  ##seqstart <- setNames(c(0, seqcumlen), seq.names)
  ##bin.coords <- bin.coords[as.vector(GenomicRanges::seqnames(bin.coords)) %in%
  ##                           seq.names]
  ##bin.signals <- bin.signals[names(bin.coords)]
  ##GenomicRanges::values(bin.coords)$bin.mids <-
  ##  (GenomicRanges::start(bin.coords) +
  ##  GenomicRanges::end(bin.coords))/2
  ##GenomicRanges::values(bin.coords)$bin.x <-
  ##  seqstart[as.character(GenomicRanges::seqnames(bin.coords))] +
  ##  bin.coords$bin.mids

  ##p <- ggplot2::qplot(bin.coords$bin.x/total.length,
  ##                    bin.signals,
  ##                    color = bin.signals,
  ##                    alpha = I(0.8))

  p <- ggplot2::ggplot(data=seg,
                       ggplot2::aes(x=floor((loc.start+loc.end)/2),
                                     y=seg.mean,
                                     color=seg.mean)
                      ) +
    ggplot2::geom_point() +
    ggplot2::geom_segment(ggplot2::aes(x    = loc.start,
                                       xend = loc.end,
                                       y    = seg.mean,
                                       yend = seg.mean),
                 size = 1.5,
                 color = "blue") +
    ggplot2::theme_bw() +
    ggplot2::scale_color_gradient2(limits = c(-0.3,0.3),
                                   low = "red",
                                   mid = "grey",
                                   high = "green",
                                   oob = scales::squish) +
      ggplot2::xlab("") +
      ggplot2::ylab("") +
      ##ggplot2::facet_wrap(~chrom) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::geom_vline(xintercept = seqstart[-1]/total.length,
                                 alpha = I(0.5)) +
      ggplot2::scale_x_continuous(labels = seq.names,
               breaks = (loc.start + loc.end/2)/total.length) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                        hjust = 0.5))

    return(p)

  ##seg.beg <-(seqstart[seg.signals$chrom] + seg.signals$loc.start)/total.length
  ##seg.end <-(seqstart[seg.signals$chrom] + seg.signals$loc.end)/total.length

  ##p <- p + ggplot2::geom_vline(xintercept = seqstart[-1]/total.length,
  ##                             alpha = I(0.5))
  ##p <- p + ggplot2::scale_x_continuous(labels = seq.names,
  ##  breaks = (seqstart + seqlen/2)/total.length) +
  ##  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
  ##                                                    hjust = 0.5))

}
