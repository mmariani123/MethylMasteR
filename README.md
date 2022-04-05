# MethylMasteR

Michael Mariani PhD, Salas Lab, Dartmouth College 2022

## Install Rstudio version:

1.) install Docker on your OS of choice <br>
<br>
2.) Then run: <br>
    docker pull mmariani123/methylmaster:latest <br>
<br>
3.) Then start the container with local_path set to the desired output directory <br>
    e.g.: <br>
    user_name=<your_user_name> <br>
    local_path="C:\Users\$user_name\Desktop\rocker_test" <br>
    docker run --rm -v $local_path:/home/rstudio -p 127.0.0.1:8787:8787 -e DISABLE_AUTH=true mmariani123/methylmaster <br>
<br>
4.) Open your web browser of choice and navigate to http://127.0.0.1:8787/ <br>
<br>
5.) Follow the commands below in the Rocker RStudio session ... <br>
    <b>(skip down to the "Run MethylMasteR" section)</b> <br>

## Install command line version for use on Slurm HPC:

<b>***If MethylMasteR singularity image (".sif" file") has already been created, skip to the next section </b> <br>

1.) Make sure that Singularity is installed on your cluster <br>
<br>
2.) Then run: <br>
    "singularity pull methylmaster_base_slim_script.sif docker://mmariani123/methylmaster_base_slim_script:latest"" <br>
    This will pull the Docker container and create a singularity container file (".sif") <br>
    (Note: you can change the name of "methylmaster_base_slim_script.sif" to anything that you want) <br>

## Run command line version for use on Slurm HPC:

1.) create a working directory on your HPC where you would like to run MethylMasteR <br>
    e.g. "mkdir /dartfs/rc/lab/S/SalasLab/programs/methylmaster_testing" <br>
    (Note: you can name the folder whatever you like) <br>
<br>
2.) Create an R script using the example R code below in the "Run MethylMasteR" section <br>
    You can adjust your own parameters as desired but MAKE SURE to name the R file: "run_methylmaster.R". <br>
    Place the "run_methylmaster.R" file in the working directory you created above. <br>
<br>
3.) Then, include the following command in a .bash file (name it whatever you like) <br>
    and submit this bash file on your Slurm HPC <br>
    /usr/bin/singularity run \ <br>
    -B /dartfs/rc/lab/S/SalasLab/programs/methylmaster_testing<b>:</b>/home/docker/work \ <br>
    /dartfs/rc/lab/S/SalasLab/programs/methylmaster_base_slim_script.sif <br>
<br>
    (Note: make sure that "/dartfs/rc/lab/S/SalasLab/programs/methylmaster_testing" is set <br>
    to your working directory (remember it does not have to be called "methylmaster_testing" <br>
    but DO NOT change "/home/docker/work" because this path is set inside the container. <br>
    Also, make sure that "/dartfs/rc/lab/S/SalasLab/programs/methylmaster_base_slim_script.sif" <br>
    is the path to the .sif file that was created) <br>
<br>
4.) Files will be output to the working directory that you created above <br>

# Run MethylMasteR

## Troubleshooting SeSAMe data caching

If you are having SeSAMe caching issues the below update code from: <br>
https://bioconductor.org/packages/devel/bioc/vignettes/ExperimentHub/inst/doc/ExperimentHub.html#default-caching-location- update <br>
Should solve the issue

```r

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
}
package="ExperimentHub"
moveFiles(package)

```

## If the above function doesn't work try caching the sesameData directly

```r

library(sesameData)
sesameDataCacheAll()

```
## Load the test data:

```r

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

Set <routine.run>="sesame" to try the sesame routine <br>
with the test data specified above.

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
  genome.version            = "hg19", #Or can set to "hg38"
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
