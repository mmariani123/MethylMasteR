# MethylMasteR

Michael Mariani PhD, Salas Lab, Dartmouth College 2022

## Install 

Follow the commands below

1.) install Docker on your OS of choice <br>
2.) Then run "docker pull mmariani123/methylmaster:lastest" <br>
3.) Then start the container with local_path set to the desired output directory <br>
    e.g.: <br>
    user_name=<your_user_name> <br>
    local_path="C:\Users\$user_name\Desktop\rocker_test" <br>
    docker run --rm -v $local_path:/home/rstudio -p 127.0.0.1:8787:8787 -e DISABLE_AUTH=true methylmaster <br>
4.) Open your web browser of choice and navigate to http://127.0.0.1:8787/ <br>
5.) Follow the commands below in the Rocker RStudio session ... <br>

## Run MethylMasteR

## Load the test data:

```r

# If you are having SeSAMe caching issues the below update code from:
# https://bioconductor.org/packages/devel/bioc/vignettes/ExperimentHub/inst/doc/ExperimentHub.html#default-caching-location- update
# Should solve the issue

moveFiles<-function(package){
olddir <- path.expand(rappdirs::user_cache_dir(appname=package))
newdir <- tools::R_user_dir(package, which="cache")
dir.create(path=newdir, recursive=TRUE)
files <- list.files(olddir, full.names =TRUE)
moveres <- vapply(files,
    FUN=function(fl){
    filename = basename(fl)
    newname = file.path(newdir, filename)
    file.rename(fl, newname)
    },
    FUN.VALUE = logical(1))
if(all(moveres)) unlink(olddir, recursive=TRUE)

package="ExperimentHub"
moveFiles(package)
       
library(MethylMasteR)

input.dir <- system.file("extdata",
              package = 'MethylMasteR')

#Create an output dir for the SeSAMe test data:
output.dir <- paste0(getwd(),.Platform$file.sep,"sesame")
              
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

Note that it is recommended to set <form.thresholds> to c(-0.2,0.2)
when running the test data specified above.

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
