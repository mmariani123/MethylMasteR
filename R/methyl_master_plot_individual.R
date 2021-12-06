#!/usr/bin/env Rscript

##Michael Mariani PhD Dartmouth College 2021

plot_indivdual <- function(output.dir,
                           seg){

  individual.plots.dir <- paste0(output.dir,
                                 .Platform$file.sep,
                                 "indivdual_plots_",
                                 names(seg))

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
  for(i in 1:length(seg)){

    sample.now <- names(seg)[i]
    ggsave({sesame::visualizeSegments(seg[[i]]) +
        theme_bw() +
        ylab("Signal") +
        labs(color="Sesame seg. \nsignal") +
        ggtitle(paste0("Sample ",
                       sample.now,
                       "\nRef = ",
                       ref)) +
        theme(axis.text.x = element_text(angle=90),
              plot.title = element_text(hjust=0.5))},
        filename = paste0(individual.plots.dir,
                          file.sep,
                          sample.now,
                          "_",
                          ref,
                          "_sesame_plot.pdf"),
        width=12,
        height=8,
        device="pdf")

  }

}
