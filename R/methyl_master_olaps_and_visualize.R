#!/usr/bin/env Rscript

##Michael Mariani Dartmouth College 2021

##For Testing:
##library(magrittr)
##load("C:\\Users\\Mike\\Desktop\\cnv_testing\\sesame_cord\\seg.RData")
##routine    <- "sesame"
##output.dir <- "C:\\Users\\Mike\\Desktop\\cnv_testing\\sesame_cord"
##file.sep   <- .Platform$file.sep
##reference  <- "internal"
##split.by <- "gender_reported"
##comparison <- c("tumor", "normal")
##
##sample.sheet.path <- paste0("C:\\Users\\Mike",
##                            "\\Desktop\\cnv_testing",
##                            "\\Sample_Sheet.csv")
##
##sample.sheet.csv <- paste0(sample.sheet.path) %>%
##  read.csv(header = TRUE,
##           stringsAsFactors = FALSE)

#' @title methyl_master_olaps_and_visualize
#' @description perform overlaps analysis and visualization
#'
#' @param ov.seg
#' @param ov.split.field
#' @param ov.keep.extra.columns
#' @param ov.overlap.density
#' @param ov.estimate.recursion
#' @param ov.simplify.reduce
#' @param ...
#' @return #seg.out
#' @export
olaps_and_visualize <- function(ov.seg                = NULL,
                                ov.split.field        = "Sample_ID",
                                ov.keep.extra.columns = TRUE,
                                ov.overlap.density    = 0.1,
                                ov.estimate.recursion = TRUE,
                                ov.simplify.reduce    = weighted.mean,
                                ...
                                ){

  grl <- GenomicRanges::makeGRangesListFromDataFrame(seg,
                                                     split.field=split.field,
                                                     keep.extra.columns=keep.extra.columns)

  grl <- GenomicRanges::sort(grl)

  ##ra  <- RaggedExperiment::RaggedExperiment(grl)

  ##Excluding 997367 copy-number neutral regions (CN state = 2, diploid)
  ##Sign-in-by error below with usual function
  ##cnvrs <- CNVRanger::populationRanges(grl,
  ##                                     density=0.1,
  ##                                     est.recur=TRUE)

  ##RaggedExperiment::assay(ra[1:5,1:5])

  ##Trims low-density areas (usually <10%
  ##of the total contributing individual
  ##calls within a summarized region).

  cnvrs <- salas_mm_population_ranges(grl,
                                      ##density=0.1,
                                      ##We can also like the density
                                      ##Jenn likes .01
                                      density         = ov.overlap.density,
                                      est.recur       = ov.estimate.recursion,
                                      simplify.reduce = ov.weighted.mean
                                      )
  ##cnvrs$type %>% unique()##
  cnvrs.filt <- subset(cnvrs, pvalue < 0.05)

  ra  <- RaggedExperiment::RaggedExperiment(grl)

  ##RaggedExperiment::assay(ra[1:5,1:5])

  ##below the less filtered ragged experiment
  if(less.stringent.ra.setting==TRUE){
    cvns.matrix <- RaggedExperiment::assay(ra)
  }else{
    ##This will likley not work as well with
    ##fewer samples so can default to the
    ##above, both return overlaps.
    cvns.matrix <-
      RaggedExperiment::qreduceAssay(ra,
                                     cnvrs.filt,
                                     simplifyReduce = weightedmean)
  }

  cvns.matrix <- cvns.matrix[order(rownames(cvns.matrix)),]
  cvns.matrix[is.na(cvns.matrix)] <- 2
  cvns.matrix <- round(cvns.matrix, 0)

  write.table(cvns.matrix,
              file = paste0(output.dir,
                            file.sep,
                            routine,
                            "_",
                            names(seg),
                            "_overlaps.csv"),
              sep=",",
              col.names=TRUE,
              row.names=TRUE,
              quote=FALSE)

  pheatmap.out <- pheatmap::pheatmap(cvns.matrix,
                                     cluster_rows = T,
                                     cluster_cols = T,
                                     show_colnames = T)

  ggsave(pheatmap.out,
         file = paste0(output.dir,
                       file.sep,
                       routine,
                       "_",
                       names(seg),
                       "_pheatmap.pdf"),
         device="pdf",
         width=12,
         height=8)

  ######################### HEATMAPS ######################################

  seg$status <- ifelse(seg$state>2,
                       "gain",
                       ifelse(seg$state<2,
                              "loss",
                              "neutral"))

  seg$chrom  <- paste0("chr",seg$chrom)
  seg$chrom  <- factor(seg$chrom,
                       ##levels=gtools::mixedsort(
                       ##unique(ex.data.in$chrom))
                       levels=c("chr1",
                                "chr2",
                                "chr3",
                                "chr4",
                                "chr5",
                                "chr6",
                                "chr7",
                                "chr8",
                                "chr9",
                                "chr10",
                                "chr11",
                                "chr12",
                                "chr13",
                                "chr14",
                                "chr15",
                                "chr16",
                                "chr17",
                                "chr18",
                                "chr19",
                                "chr20",
                                "chr21",
                                "chr22",
                                "chrX",
                                "chrY"))

  seg$Sample_ID <-
    factor(seg$Sample_ID,
           levels=unique(gtools::mixedsort(seg$Sample_ID)))

  seg$status <- factor(seg$status,
                       levels = c("loss",
                                  "neutral",
                                  "gain"))

  ##floor(seg.vis.total$loc.start + seg.vis.total$loc.end) /2)
  seg.heatmap <- ggplot(seg,
                        aes(x=as.integer((floor((loc.start + loc.end)/2))),
                            ##x=chrom,
                            ##y=seg.mean,
                            y=0.5,
                            width=loc.start - loc.end,
                            fill=status)) +
    geom_tile() +
    ##geom_vline(xintercept=0) +
    facet_grid(rows=vars(Sample_ID),
               cols=vars(chrom),
               scales = "free",
               space = "free",
               switch = "y") +
    theme_bw() +
    theme(panel.spacing = unit(0,"cm"),
          axis.ticks.x  = element_blank(),
          axis.ticks.y  = element_blank(),
          axis.text.x   = element_blank(),
          axis.text.y   = element_blank(),
          axis.text.x.bottom  = element_blank(),
          axis.text.y.left = element_blank(),
          axis.title.y.left = element_blank(),
          axis.title.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),
          strip.text.y.left = element_text(angle = 0),
          strip.text.x.top = element_text(angle = 90)) +
    ylim(c(0,1)) +
    scale_fill_manual(values=c("orange",
                               "white",
                               "blue"))

  ggsave(seg.heatmap,
         filename = paste0(output.dir,
                           file.sep,
                           routine,
                           "_",
                           sub.workflow,
                           "_",
                           ref,
                           "_mmheatmap.pdf"),
         width=12,
         height=8,
         device="pdf")

}
