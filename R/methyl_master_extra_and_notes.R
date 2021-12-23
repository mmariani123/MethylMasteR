#1/usr/bin/file Rscript

##Methyl Master Extra and Notes

##ChAMP

##<control>
##If champ.CNA() should be calculate copy number
##variance between case and control? (The other
##option for champ.CNA() is calculate copy number
##variance for each sample to the averaged value).
##(default = TRUE)

##<controlGroup>
##Which pheno should be treated as control
##group while running champ.CNA().(default = "champCtls")

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
##"?..Sample_Name"
##"Sample_Group"
##"Sample_Plate"
##"Pool_ID"
##"Sample_Well"
##"Sentrix_Position"
##"Sentrix_ID"
##"Basename"
##"Batch"

##champ.sample.sheet.df %>% dplyr::pull(`?..Sample_Name`)
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
##                              dplyr::pull(`?..Sample_Name`)
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
##  dplyr::pull(`?..Sample_Name`)
##names(champ_results$champ.CNA$sampleResult) ##Looks good
