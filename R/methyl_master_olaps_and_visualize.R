#!/usr/bin/env Rscript

#' @title methyl_master_olaps_and_visualize
#' @description perform overlaps analysis and visualization
#' ##Michael Mariani Dartmouth College 2021
#' @param ov.seg
#' @param ov.name
#' @param ov.output.dir
#' @param ov.routine
#' @param ov.split.field
#' @param ov.keep.extra.columns
#' @param ov.overlap.density
#' @param ov.estimate.recurrence
#' @param ov.keep.extra.columns
#' @param ov.simplify.reduce
#' @param ov.less.stringent.ra.setting
#' @param ov.pvalue
#' @param ov.plot.individual
#' @param ...
#' @import pheatmap
#' @import ggplot2
#' @importFrom gtools mixedorder unbox
#' @importFrom GenomicRanges as.data.frame unbox
#' @importFrom GenomicRanges makeGRangesListFromDataFrame unbox
#' @importFrom ramify flatten unbox
#' @return #seg.out
#' @export
methyl_master_olaps_and_visualize <- function(ov.seg   = NULL,
                                              ov.name  = NULL,
                                ov.output.dir          = getwd(),
                                ov.routine             = NULL,
                                ov.split.field         = NULL,
                                ov.keep.extra.columns  = TRUE,
                                ov.overlap.density     = 0.1,
                                ov.estimate.recurrence = FALSE,
                                ov.simplify.reduce     = weightedmean,
                                ov.less.stringent.ra.setting = FALSE,
                                ov.pvalue              = 0.05,
                                ov.plot.individual     = FALSE,
                                ...
                                ){

  seg <- ov.seg

  ######################### MY HEATMAPS ######################################

  seg <- seg[seg$pval <= ov.pvalue,]

  if(ov.plot.individual==TRUE & ov.routine!="sesame"){

    print("Individual plots only works for sesmae routine currently")

  }else if(ov.plot.individual==TRUE){
    methyl_master_plot_individual(pi.seg = seg,
                                pi.output.dir = ov.output.dir,
                                pi.name = ov.routine)
  }

  seg$status <- ifelse(seg$state>2,
                       "gain",
                       ifelse(seg$state<2,
                              "loss",
                              "neutral"))

  if(!grepl("chr", seg$chrom, perl=TRUE)){
    seg$chrom  <- paste0("chr", seg$chrom)
  }

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

  seg$loc.start <- as.numeric(seg$loc.start)
  seg$loc.end   <- as.numeric(seg$loc.end)

  ##floor(seg.vis.total$loc.start + seg.vis.total$loc.end) /2)
  seg.heatmap <- ggplot2::ggplot(seg,
                                 ggplot2::aes(
                                 x=as.integer((floor((loc.start + loc.end)/2))),
                                              ##x=chrom,
                                              ##y=seg.mean,
                                              y=0.5,
                                              width=loc.start - loc.end,
                                              fill=status)) +
    ggplot2::geom_tile() +
    ##geom_vline(xintercept=0) +
    ggplot2::facet_grid(rows=vars(Sample_ID),
                        cols=vars(chrom),
                        scales = "free",
                        space = "free",
                        switch = "y") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.spacing = ggplot2::unit(0,"cm"),
                   axis.ticks.x        = ggplot2::element_blank(),
                   axis.ticks.y        = ggplot2::element_blank(),
                   axis.text.x         = ggplot2::element_blank(),
                   axis.text.y         = ggplot2::element_blank(),
                   axis.text.x.bottom  = ggplot2::element_blank(),
                   axis.text.y.left    = ggplot2::element_blank(),
                   axis.title.y.left   = ggplot2::element_blank(),
                   axis.title.x        = ggplot2::element_blank(),
                   panel.grid.minor.x  = ggplot2::element_blank(),
                   panel.grid.major.x  = ggplot2::element_blank(),
                   panel.grid.minor.y  = ggplot2::element_blank(),
                   panel.grid.major.y  = ggplot2::element_blank(),
                   strip.text.y.left   = ggplot2::element_text(angle = 0),
                   strip.text.x.top    = ggplot2::element_text(angle = 90)) +
    ggplot2::ylim(c(0,1)) +
    ggplot2::scale_fill_manual(values=c("orange",
                                        "white",
                                        "blue"))

  ggplot2::ggsave(seg.heatmap,
                  filename = paste0(ov.output.dir,
                                    .Platform$file.sep,
                                    ov.routine,
                                    "_",
                                    ov.name,
                                    "_mmheatmap.pdf"),
                  width=12,
                  height=8,
                  device="pdf")

############################### PHEATMAP #####################################

  cnvs.matrix <- reshape2::dcast(seg,
                         paste0(chrom,
                         ":",
                         as.integer(loc.start),
                         "-",
                         as.integer(loc.end)
                         )~Sample_ID, value.var = "state")
  cnvs.rownames <- cnvs.matrix[,1]
  cnvs.matrix <- cnvs.matrix[,-1]
  rownames(cnvs.matrix) <- cnvs.rownames
  cnvs.matrix[is.na(cnvs.matrix)] <- 2

  cluster.cols <- TRUE
  cluster.rows <- TRUE

  if(length(unique(ramify::flatten(cnvs.matrix,
                                   across = c("rows", "columns"))))==1){

    print("All the same cnv values in matrix, not creating pheatmap")

  }else{

    if(is.null(ncol(cnvs.matrix))){
      cluster.cols <- FALSE
    }
    if(is.null(nrow(cnvs.matrix))){
      cluster.rows <- FALSE
    }
    pheatmap.out <- pheatmap::pheatmap(cnvs.matrix,
                                       cluster_rows = cluster.cols,
                                       cluster_cols = cluster.rows,
                                       show_colnames = T,
                                       show_rownames = T)
    ggplot2::ggsave(pheatmap.out,
                    file = paste0(ov.output.dir,
                                  .Platform$file.sep,
                                  ov.routine,
                                  "_",
                                  ov.name,
                                  "_pheatmap.pdf"),
                    device="pdf",
                    width=12,
                    height=8)
  }

  cnvs.matrix <- cnvs.matrix[gtools::mixedorder(rownames(cnvs.matrix)),]

  write.table(GenomicRanges::as.data.frame(cnvs.matrix),
              file = paste0(ov.output.dir,
                            .Platform$file.sep,
                            ov.routine,
                            "_",
                            ov.name,
                            "_pheatmap_table.csv"),
              sep=",",
              col.names=NA,
              row.names=TRUE,
              quote=FALSE)


########################### OVERLAPS ##########################################

  seg <- ov.seg
  rm(ov.seg)

  grl <- GenomicRanges::makeGRangesListFromDataFrame(seg,
                              split.field=ov.split.field,
                 keep.extra.columns=ov.keep.extra.columns)

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

  cnvrs <- methyl_master_population_ranges(grl,
                                      ##density=0.1,
                                      ##We can also like the density
                                      ##Jenn likes .01
                                      density         = ov.overlap.density,
                                      est.recur       = ov.estimate.recurrence
                                      )
  ##cnvrs$type %>% unique()##

  cnvrs.filt <- subset(cnvrs, pvalue <= ov.pvalue)

  ##ra  <- RaggedExperiment::RaggedExperiment(grl)
  ##http://bioconductor.org/packages/release/bioc/vignettes/
  ##RaggedExperiment/inst/doc/RaggedExperiment.html

  ##The simplifyReduce argument in qreduceAssay allows
  ##the user to summarize overlapping regions with a custom
  ##method for the given “query” region of interest.
  ##We provide one for calculating a weighted average
  ##score per query range, where the weight is proportional
  ##to the overlap widths between overlapping ranges and a query range.

  ##Note that there are three arguments to this function.
  ##See the documentation for additional details.

  ##From Bioconductor manual it appears that RaggedExperiment::Assay()
  ##defaults to sparseAssay(): A matrix() with dimensions dim(x).
  ##Elements contain the assay value for the ith range and jth sample.
  ##Use ’sparse=TRUE’ to obtain a sparseMatrix assay representation.

  ##RaggedExperiment::assay(ra[1:5,1:5])

  ##below the less filtered ragged experiment
  ##if(ov.less.stringent.ra.setting==TRUE){
  ##  cnvs.matrix <- RaggedExperiment::assay(x=ra,i="state")
  ##}else{
  ##  ##This will likley not work as well with
  ##  ##fewer samples so can default to the
  ##  ##above, both return overlaps.
  ##  cnvs.matrix <-
  ##    methyl_master_qreduceassay(ra,
  ##                               cnvrs.filt,
  ##                               simplifyReduce = ov.simplify.reduce)
  ##}

  ##cnvs.matrix <- methyl_master_assay(ra, i="state")

  ##cnvs.matrix <- cnvs.matrix[order(rownames(cnvs.matrix)),]
  ##cnvs.matrix[is.na(cnvs.matrix)] <- 2
  ##cnvs.matrix <- round(cnvs.matrix, 0)

  write.table(GenomicRanges::as.data.frame(cnvrs),
              file = paste0(ov.output.dir,
                            .Platform$file.sep,
                            ov.routine,
                            "_",
                            ov.name,
                            "_overlaps.csv"),
              sep=",",
              col.names=NA,
              row.names=TRUE,
              quote=FALSE)

  write.table(GenomicRanges::as.data.frame(cnvrs.filt),
              file = paste0(ov.output.dir,
                            .Platform$file.sep,
                            ov.routine,
                            "_",
                            ov.name,
                            "_overlaps_filt.csv"),
              sep=",",
              col.names=NA,
              row.names=TRUE,
              quote=FALSE)

}
