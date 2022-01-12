# MethylMasteR
Michael Mariani PhD, Salas Lab, Dartmouth College 2022

#Install

require(devtools)\
devtools::install_github(mmariani123/MethylMasteR)
library(MethylMasteR)

## Run MethylMasteR

## Set your paths:

input.dir <- "path.to.idat.dir"\
output.dir <- "path.to.output.dir"\
sample.sheet.path <- "path.to.sample.sheet"\

## Select your routine :

1.) "sesame" for SeSAMEe\
2.) "hm450" for HM450\
3.) "champ" for ChAMP\
4.) "epicopy" for Epicopy\
5.) "compare" for comparison routine*

*For the compairson routine you will need to have output from two more
results folders from the other four routines

## Select the other parameters and run!

methyl_master(\
  routine                   = "sesame",          #The routine to run\
  input.dir                 = input.dir,         #The input (idat.files) directory\
  output.dir                = output.dir,        #The output directory\
  sample.sheet.path         = sample.sheet.path, The path to the MethylMasteR 
                                                 sample sheet\
  r.lib.path                = .libPaths()[1],    #The path to the R Library path\
  file.sep                  = "\\\\\\\\", #for windows or "/" for Linux\
  n.cores                   = 1, #Multicore does not work for all software on all platforms
  os.type                   = "windows", #or "linux"\
  proj                      = "TCGA-KIRC", #"TCGA-BLCA" etc.\
  visualize                 = TRUE, #output plots,
  visualize.individual      = FALSE, #Indvidual sample plots only works for routine sesame\
  reference                 = "internal", #"comparison" or 'internal"\
  reference.name            = NA, #For Epicopy use NA for median, "all" is not currently supported by Epicopy\
  comparison                = c("tumor","cord"), #Always required treatment first then control\
                              #can also be more specific when designing sample sheet and use values like\
                              such as c("tumor_male","cord_male etc") etc.\
                              #in routines where internal reference is used, second argument is ignored\
  form.thresholds           = NULL, # Used to calculate final CNV state. If NULL, equation is is used,\
                              otherwise specify threshold vector of lower and upper beta values such\
                              as c(-0.3,0.3),\
  overlap.density           = 0.1, #For combining final CNV calls for confidence\
  save.seg                  = TRUE, #Whether to save segmentation results\
  sesame.data.cache         = "EPIC", #The default sesame reference platform\
  sesame.data.normal        = 'EPIC.5.normal', #The default sesame and hm450 internal reference samples\
  sesame.ref.version        = "hg19", #Or can set to "hg38"\
  hm450.workflow            = "B", #The HM450 subworkflow to use 9Only B is running currently:\
                              #"A" no correct,\
                              #"B" median correction (default),\
                              #"C" run Conumee\
  champ.padj                = 0.05, #padj to filter champ results\
  champ.control             = FALSE, #run champ.control etc.\
  champ.run.combat          = FALSE, #run champ.run.combat etc.\
  champ.run.dmp             = FALSE, #If only one pheno var must = FALSE\
  champ.run.dmr             = FALSE, #If only one pheno var must = FALSE\
  champ.run.block           = FALSE, #If only one pheno var must = FALSE\
  champ.run.gsea            = FALSE, #Requires dmp and dmr results\
  champ.run.epimod          = FALSE, #If only one pheno var must = FALSE\
  epi.run.gistic            = TRUE, #Whether to Run GISTIC in Epicopy workflow\
  olaps.split.field         = "Sample_ID", #Split field to ise during overlaps\
                                      #Don't change unless you know what you are doing\
  estimate.recurrence       = TRUE, #Estimate recursion to produce p values when finding overlaps\
                               #with population ranges functions\
  ov.pvalue                 = 0.05, #pvalue threshold for overlaps identified\
  ov.keep.extra.columns     = TRUE, #Keep extra metadata columns when finding overlaps\
  simplify.reduce           = weightedmean, #Equation to use during reduction\
  create.dir                = TRUE, #create directory if does not exist?\
  compare.list.files        = FALSE, #Set to true to list files below, otherwise\
                              The "compare" routine will search the input.dir\
                              recursively for output to compare\
  compare.files.in          = compare.files.in, #vector list of individual files for the compare routine\
  compare.names             = compare.names, #Specify name of each comparison result\
  compare.olaps.size        = 1 #Overlap of one or more base pairs\
                              #considered as an overlap in the "compare" routine\
)
