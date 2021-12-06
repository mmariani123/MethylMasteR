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

  ##visualize individual sesame by sample:
  for(i in 1:length(pi.seg)){

    sample.now <- names(pi.seg)[i]
    ggplot2::ggsave({sesame::visualizeSegments(pi.seg[[i]]) +
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
