# MethylMasteR

Michael Mariani PhD, Salas Lab, Dartmouth College 2022

## Install

```r

require(devtools)
devtools::install_github("mmariani123/methylmaster")
library(MethylMasteR)

```
## Run MethylMasteR

## Load the test data:

```r

input.dir <- system.file("extdata",
              package = 'MethylMasteR')
              
output.dir <- "path/to/output/dir"

sample.sheet.path <- system.file("extdata",
                   "Sample_Sheet_Test.csv",
                   package = 'MethylMasteR')

```

## Select your routine :

Run SeSAMe with the test data specified above

```r

routine.run <- "sesame"
#1.) "sesame" for SeSAMEe
#2.) "hm450" for HM450
#3.) "champ" for ChAMP
#4.) "epicopy" for Epicopy
#5.) "custom" for custom routine

```
## Select the other parameters and run!

```r

methyl_master(
  routine                   = routine.run, #The routine to run
  input.dir                 = input.dir, #The input (idat.files) directory
  output.dir                = output.dir, #The output directory
  sample.sheet.path         = sample.sheet.path, #The path to the MethylMasteR sample sheet
  r.lib.path                = .libPaths()[1], #The path to the R Library path
  file.sep                  = "\\\\", #For windows or "/" for Linux
  create.dir                = TRUE, #Whether to cretae directory if does not 
                                    #exist?
  save.seg                  = TRUE, #Whether to save segmentation results
  n.cores                   = 1, #Multicore does not work for all routines 
                                 #on all operating systems
  os.type                   = "windows", #Or "linux"
  proj                      = "TCGA-KIRC", #"TCGA-BLCA" etc.
  visualize                 = TRUE, #Whether to output plots,
  visualize.individual      = FALSE, #Whether to output indvidual sample plots,
                                     #only works for routine sesame
  reference                 = "internal", #"comparison" or 'internal"
  reference.name            = NA, #For Epicopy use NA for median, 
                                  #"all" is not currently supported 
                                  #by Epicopy
  comparison                = c("tumor","normal"), #Always required treatment 
                              #first then contro can also be more specific when 
                              #designing sample sheet and use values like
                              # c("tumor_male","cord_male etc") etc.
                              #Note: in routines where internal reference is 
                              #used, second argument is ignored
  form.thresholds           = NULL, #Used to calculate final CNV state. 
                                    #If NULL, equation is is used;
                                    #otherwise, specify threshold vector 
                                    #of lower and upper beta values such
                                    #as c(-0.3,0.3), we used c(-0.2,0.2) 
                                    #in the paper
  overlap.density           = 0.1, #For combining final CNV calls for confidence
  sesame.data.cache         = "EPIC", #The default sesame reference platform
  sesame.data.normal        = 'EPIC.5.normal', #The default sesame and hm450 
                                               #internal reference samples
  sesame.ref.version        = "hg19", #Or can set to "hg38"
  hm450.workflow            = "B", #The HM450 subworkflow to use - only B 
                                   #is running currently.
                              #"A" no correction,
                              #"B" median correction (default),
                              #"C" run Conumee
  champ.padj                = 0.05, #padj to filter champ results
  champ.control             = FALSE, #run champ.control etc.
  champ.run.combat          = FALSE, #run champ.run.combat etc.
  champ.run.dmp             = FALSE, #If only one pheno var must = FALSE
  champ.run.dmr             = FALSE, #If only one pheno var must = FALSE
  champ.run.block           = FALSE, #If only one pheno var must = FALSE
  champ.run.gsea            = FALSE, #Requires dmp and dmr results
  champ.run.epimod          = FALSE, #If only one pheno var must = FALSE
  epi.run.gistic            = TRUE, #Whether to Run GISTIC in Epicopy workflow
  olaps.split.field         = "Sample_ID", #Split field to ise during overlaps
                                           #Don't change unless you know what 
                                           #you are doing
  estimate.recurrence       = TRUE, #Estimate recursion to produce p values when 
                                    #finding overlaps with population_ranges 
                                    #functions
  ov.pvalue                 = 0.05, #pvalue threshold for overlaps identified
  ov.keep.extra.columns     = TRUE, #Keep extra metadata columns when finding 
                                    #overlaps
  simplify.reduce           = weightedmean #Equation to use during reduction
)

```
