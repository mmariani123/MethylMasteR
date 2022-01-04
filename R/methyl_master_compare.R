#!/usr/bin/env Rscript

#' @title methyl_master_compare
#' @description compare results from different routine outputs
#' Michael Mariani PhD Dartmouth College 2021
#' @param compare.list.files
#' @param compare.files.in
#' @param compare.input.dir
#' @param compare.output.dir
#' @param compare.names
#' @param compare.olaps.size
#' @param ...
#' @import ggplot2
#' @import foreach
#' @importFrom GenomicRanges GRangesList unbox
#' @importFrom GenomicRanges makeGRangesFromDataFrame unbox
#' @importFrom cowplot plotgrid
#' @importFrom S4Vectors mcols unbox
#' @importFrom seqsetvis ssvOverlapIntervalSets unbox
#' @importFrom seqsetvis ssvFeatureUpset unbox
#' @importFrom rlist list.append unbox
#' @return #outputs a comparison plot to <compare.output.dir>
#' @export
methyl_master_compare <- function(compare.list.files=FALSE,
                                  compare.files.in=NULL,
                                  compare.input.dir=getwd(),
                                  compare.output.dir=getwd(),
                                  compare.names=NULL,
                                  compare.olaps.size=1,
                                  ...
){

############################# Time-Mem ######################################

if(compare.list.files==TRUE){

files.in <- compare.files.in
files.in <- files.in[grepl("time_mem",files.in)]

}else{

files.in <- list.files(path=compare.input.dir,
                       pattern="time_mem",
                       full.names = TRUE,
                       recursive = TRUE)

compare.names <- gsub(".*/","",dirname(files.in))

}

items <- list()
foreach(i = files.in) %do% {
  con = file(i, "r")
  while(TRUE){
    lines <- readLines(con,n=1)
    if(length(lines)==0){break}
    if(grepl("Total time:",lines,perl=TRUE) |
        grepl("peakRAM::peakRAM\\(\\) max output \\(bytes\\):",lines)){
      lines <- readLines(con,n=2)
      items <- rlist::list.append(items,lines)
    }
  }
  close(con)
  items
}
items   <- unlist(items)[c(FALSE,TRUE)]
time    <- items[c(TRUE,FALSE)]
max.mem <- items[c(FALSE,TRUE)]

time.mem.df <- data.frame(sample=compare.names,
                          time=round(as.numeric(time),digits = 2),
                          max.mem=max.mem,
                          stringsAsFactors = FALSE)

time.plot <- ggplot2::ggplot(time.mem.df,
                             ggplot2::aes(x=sample,y=time,fill=sample)) +
         ggplot2::geom_col() +
         ggplot2::theme_bw() +
         ggplot2::ylab("Time (m)") +
         ggplot2::xlab("Sample") +
         ggplot2::geom_text(ggplot2::aes(y=time, label=time), vjust=-0.2) +
         ggplot2::geom_blank(ggplot2::aes(x = sample, y = time*1.1)) +
  theme(axis.text.x = element_text(angle=90))

mem.plot <- ggplot2::ggplot(time.mem.df,
                            ggplot2::aes(x=sample,y=as.numeric(max.mem),
                                         fill=sample)) +
         ggplot2::geom_col() +
         ggplot2::theme_bw() +
         ggplot2::ylab("Max mem. (MiB)") +
         ggplot2::xlab("Sample") +
  ggplot2::geom_text(ggplot2::aes(y=as.numeric(max.mem),
                                  label=max.mem), vjust=-0.2) +
  ggplot2::geom_blank(ggplot2::aes(x = sample, y = as.numeric(max.mem)*1.1)) +
  theme(axis.text.x = element_text(angle=90))

############################# Overlaps #######################################

if(compare.list.files==TRUE){

  olaps.in <- compare.files.in
  olaps.in <- olaps.in[grepl(".csv$",olaps.in)]

}else{

  olaps.in <- list.files(path=compare.input.dir,
                         pattern="overlaps_filt.csv$",
                         full.names = TRUE,
                         recursive = TRUE)

}

olaps.df.list <- foreach(i = olaps.in) %do% {
  olaps.df <- read.csv(i,
                       header = TRUE,
                       stringsAsFactors = FALSE)
   ##splitted <-
    ##unlist(strsplit(unlist(strsplit(rownames(olaps.df),split=":")),split="-"))
   ##olaps.df$chr   <- splitted[c(TRUE,FALSE,FALSE)]
   ##olaps.df$start <- splitted[c(FALSE,TRUE,FALSE)]
   ##olaps.df$end   <- splitted[c(FALSE,FALSE,TRUE)]
   ##olaps.gr <- GenomicRanges::makeGRangesFromDataFrame(olaps.df.list[[1]],
   ##                                      keep.extra.columns = TRUE)

  olaps.gr <- GenomicRanges::makeGRangesFromDataFrame(olaps.df,
                                                      keep.extra.columns = TRUE)

  if(nrow(olaps.gr!=0)){
    name(olaps.gr) <- compare.names[i]
    return(olaps.gr)
  }

}

##olaps.out <- findOverlaps(olaps.gr.list)
olaps.gr.list <- GenomicRanges::GRangesList(olaps.df.list)
ssv.out <- S4Vectors::mcols(
  seqsetvis::ssvOverlapIntervalSets(olaps.gr.list,
                                    minoverlap=compare.olaps.size)
  )
colnames(ssv.out) <- compare.names
olaps.plot <- seqsetvis::ssvFeatureUpset(ssv.out,
                                         sets.x.label="Overlap category",
                                         mainbar.y.label=paste0("Number of ",
                                         as.character(compare.olaps.size),
                                         " bp \nCNV overlaps"))

##findOverlapsOfPeaks(unlist(olaps.gr.list))
##findOverlaps(olaps.gr.list)
##UpSetR::upset(make_comb_mat(olaps.gr.list))
##x <- genomicElementUpSetR(overlaps)
##upset(x$plotData, nsets=13, nintersects=NA)

final.plot <- cowplot::plot_grid(time.plot,
                                 mem.plot,
                                 olaps.plot,
                                 nrow=3,
                                 align="vh",
                                 axis="lr",
                                 rel_widths = c(1,1,10),
                                 rel_heights = c(1,1,2))

ggplot2::ggsave(final.plot,
       filename = paste0(compare.output.dir,
                         .Platform$file.sep,
                         "sample.comparison.plots.pdf"),
       height=12,
       width=6,
       device="pdf")

}
