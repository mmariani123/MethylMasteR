#!/usr/bin/env Rscript

## Can test ChAMP with the test data: data(EPICSimData)

##The full pipeline can be run at once with one command:

##Testing:
##data(EPICSimData)
##champ.test.directory <- system.file("extdata", 
##                      package="ChAMPdata")
##myLoad.test <- champ.load(champ.test.directory,
##                     arraytype="450K")
##myLoad.test$pd$Sample_Name
##myLoad.test$pd$Sample_Group
##champ.results <- champ.process(directory=champ.test.directory)
##myLoad <- champ.load(idat.pooled.files.dir,
##                     arraytype="450K")
##myLoad$pd$Pool_ID
##myLoad$pd$Pool_ID
##myLoad$pd$Pool_ID
##myLoad$pd$Sample_Group

champ_results <- ChAMP::champ.process(directory = idat.pooled.files.dir,
                                      arraytype = champ.array.type,
                                      ##batchname = c("Slide"),
                                      batchname = champ.batch.name,
                                      adjPVal = champ.padj,
                                      cores = champ.ncores,
                                      control = champ.control,
                                      controlGroup = champ.contol.group,
                                      runimpute = champ.runimpute,
                                      runQC     = champ.runQC,
                                      runnorm   = champ.runnorm,
                                      runSVD    = champ.runSVD,
                                      runCombat = champ.runCombat,
                                      runDMP    = champ.runDMP,
                                      runDMR    = champ.runDMR,
                                      runBlock  = champ.runBlock,
                                      runGSEA   = champ.runGSEA)

##save(champ_results,
##     file=paste0(work.dir,
##                 file.sep,
##                 "champ_results.RData"))

##From champ SampleSheet.CSV:
##champ.sample.sheet.df <- list.files(idat.pooled.files.dir,
##                                    pattern = ".csv",
##                                    full.names = TRUE
##                                    ) %>% read.csv(header=TRUE,
##                                                   stringsAsFactors = FALSE)
##colnames(champ.sample.sheet.df)
##"ï..Sample_Name"
##"Sample_Group"
##"Sample_Plate"
##"Pool_ID"
##"Sample_Well"     
##"Sentrix_Position"
##"Sentrix_ID"
##"Basename"
##"Batch"   

##champ.sample.sheet.df %>% dplyr::pull(`ï..Sample_Name`)
##champ.sample.sheet.df %>% dplyr::pull(Sample_Group)

## 6164655052_R05C02, ##N
## 6164655053_R05C02, ##N
## 7796806075_R04C01, ##N
## 7796806101_R06C01, ##N
## 7796806096_R01C02, ##N
## 6285650044_R06C02, ##N
## 7796806095_R04C01, ##N
## 6264496107_R03C02, ##N
## 8795194084_R01C01, ##N
## 7796806101_R04C01, ##N
## 7786915032_R01C01, ##T
## 7796806101_R03C02, ##T
## 9422491111_R06C02, ##T
## 7786915031_R06C01, ##T
## 9806233125_R01C01, ##T
## 3999510104_R03C01, ##T
## 3999510104_R04C02, ##T
## 8795194085_R04C01, ##T
## 9630789147_R02C01, ##T
## 8795194084_R03C02) ##T

##Testing:
##load(file=paste0("G:\\My Drive\\dartmouth",
##                 "\\salas_lab_working\\cnv\\testing",
##                 "\\champ_kidney_results\\CHAMP_RESULT\\myCNA.rda"))
##length(myCNA$sampleResult) ##20
##names(myCNA$sampleResult) ## All ""
##names(myCNA$sampleResult) <- champ.sample.sheet.df %>% 
##                              dplyr::pull(`ï..Sample_Name`)
##names(myCNA$sampleResult) ##Looks good
##for(i in 1:length(myCNA$sampleResult)){
##  myCNA$sampleResult[[i]]$ID <- names(myCNA$sampleResult[i])
##  colnames(myCNA$sampleResult[[i]])[1] <- "Sample_ID"
##}
##colnames(myCNA$sampleResult[[i]])[1] ##All set

##load(paste0(work.dir,
##            file.sep,
##            "champ_example_result.RData"))
##results.data <- champ.results$champ.CNA$sampleResult
##rm(champ.results)
##colnames(champ.resultschamp.CNA$sampleResult[[1]]) 

##load(file=paste0("G:\\My Drive\\dartmouth",
##                 "\\salas_lab_working\\cnv",
##                 "\\testing\\champ_results.RData"))
##length(champ_results$champ.CNA$sampleResult) ##20
##names(champ_results$champ.CNA$sampleResult) ## All ""
##names(champ_results$champ.CNA$sampleResult) <- champ.sample.sheet.df %>% 
##  dplyr::pull(`ï..Sample_Name`)
##names(champ_results$champ.CNA$sampleResult) ##Looks good
for(i in 1:length(champ_results$champ.CNA$sampleResult)){
  champ_results$champ.CNA$sampleResult[[i]]$ID <- 
    names(champ_results$champ.CNA$sampleResult[i])
  colnames(champ_results$champ.CNA$sampleResult[[i]])[1] <- "Sample_ID"
  ##colnames(myCNA$sampleResult[[i]])[1] ##All set
}
champ_results <- do.call(rbind, champ_results$champ.CNA$sampleResult)

##Add sex info for samples:
##champ_results <- champ_seg
champ_results$karyotype    <- ""
champ_results$sex_reported <- ""
champ_results$sex_inferred <- ""
champ_results$treatment    <- ""
for(i in 1:length(unique(champ_results$Sample_ID))){
  print(i)
  sample.now <- unique(champ_results$Sample_ID)[i]
  if(sample.now %in% rownames(tumor)){
    champ_results[champ_results$Sample_ID==sample.now,]$karyotype    <- 
      tumor[tumor$X==sample.now, "karyotype"]
    champ_results[champ_results$Sample_ID==sample.now,]$sex_reported <- 
      tumor[tumor$X==sample.now, "gender_reported"]
    champ_results[champ_results$Sample_ID==sample.now,]$sex_inferred  <- 
      tumor[tumor$X==sample.now, "sex_inferred"]
    champ_results[champ_results$Sample_ID==sample.now,]$treatment <- "tumor"
  }else if(sample.now %in% rownames(normal)){
    champ_results[champ_results$Sample_ID==sample.now,]$karyotype    <- 
      normal[normal$X==sample.now, "karyotype"]
    champ_results[champ_results$Sample_ID==sample.now,]$sex_reported <- 
      normal[normal$X==sample.now, "gender_reported"]
    champ_results[champ_results$Sample_ID==sample.now,]$sex_inferred <- 
      normal[normal$X==sample.now, "sex_inferred"]
    champ_results[champ_results$Sample_ID==sample.now,]$treatment <- "normal"
  }else{
    stop(paste0("Error: champ processing:",
                "samples names must belong ",
                "to either tumor or normal ",
                "treatments"))
  }
}

save(champ_results, file=paste0(work.dir, file.sep, "champ_results.RData"))
rm(champ_results)
