#!/usr/bin/env Rscript

#' @title methyl_master_compare
#' @description compare results from different routine outputs
#' Michael Mariani PhD Dartmouth College 2021
#' @param compare.list.files
#' @param compare.files.in
#' @param compare.input.dir
#' @param compare.output.dir
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
                                  compare.input.dir=NULL,
                                  compare.output.dir=NULL,
                                  ...
){

############################# Time-Mem ######################################

if(compare.list.files==TRUE){

files.in <- compare.files.in

}else{

files.in <- list.files(path=compare.input.dir,
                       pattern="time_mem",
                       full.names = TRUE,
                       recursive = TRUE)

}

items <- list()
foreach(i = files.in) %do% {
  con = file(i, "r")
  while(TRUE){
    lines <- readLines(con,n=1)
    if(length(lines)==0){break}
    if(grepl("Total time:",lines,perl=TRUE) |
        grepl("Profmem max output \\(bytes\\):",lines)){
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

time.mem.df <- data.frame(sample=c("1","2","3","4"),
                          time=round(as.numeric(time),digits = 2),
                          max.mem=max.mem,
                          stringsAsFactors = FALSE)

time.plot <- ggplot2::ggplot(time.mem.df,
                             ggplot2::aes(x=sample,y=time,fill=sample)) +
         ggplot2::geom_col() +
         ggplot2::theme_bw() +
         ggplot2::ylab("Max mem. (MB)") +
         ggplot2::xlab("Sample")

mem.plot <- ggplot2::ggplot(time.mem.df,
                            ggplot2::aes(x=sample,y=max.mem, fill=sample)) +
         ggplot2::geom_col() +
         ggplot2::theme_bw() +
         ggplot2::ylab("Max mem. (MB)") +
         ggplot2::xlab("Sample")

############################# Overlaps #######################################

if(compare.list.files==TRUE){

  olaps.in <- compare.files.in

}else{

  olaps.in <- list.files(path=compare.input.dir,
                         pattern="overlaps.csv$",
                         full.names = TRUE,
                         recursive = TRUE)

}

olaps.df.list <- foreach(i = olaps.in) %do% {
  olaps.df <- read.csv(i,
                       header = TRUE,
                       stringsAsFactors = FALSE)
   splitted <-
     unlist(strsplit(unlist(strsplit(rownames(olaps.df),split=":")),split="-"))
   olaps.df$chr   <- splitted[c(TRUE,FALSE,FALSE)]
   olaps.df$start <- splitted[c(FALSE,TRUE,FALSE)]
   olaps.df$end   <- splitted[c(FALSE,FALSE,TRUE)]
   olaps.gr <- GenomicRanges::makeGRangesFromDataFrame(olaps.df.list[[1]],
                                        keep.extra.columns = TRUE)
   return(olaps.gr)

}

names(olaps.df.list) <- c("sample 1", "sample 2", "sample 3" , "sample 4")
##olaps.out <- findOverlaps(olaps.gr.list)
olaps.gr.list <- GenomicRanges::GRangesList(olaps.df.list)
ssv.out <- S4Vectors::mcols(seqsetvis::ssvOverlapIntervalSets(olaps.gr.list))
colnames(ssv.out) <- c("sample 1", "sample 2", "sample 3" , "sample 4")
olaps.plot <- seqsetvis::ssvFeatureUpset(ssv.out)

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
                                 rel_widths = c(1,1,10))

ggplot2::ggsave(final.plot,
       filename = paste0(compare.output.dir,
                         .Platform$file.sep,
                         "sample.comparison.plots.pdf"),
       height=8,
       width=6,
       device="pdf")

}
