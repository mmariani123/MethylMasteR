#!/usr/bin/env Rscript

######################## MINFI #############################################
############################################################################

##Add in the preprocessing from k450
##Also remember the di-bias correction issue 
##preproccessnoob

load(paste0(work.dir,file.sep,"sesame_rgset_normal_female.RData"))
load(paste0(work.dir,file.sep,"sesame_rgset_normal_male.RData"))
load(paste0(work.dir,file.sep,"sesame_rgset_tumor_female.RData"))
load(paste0(work.dir,file.sep,"sesame_rgset_tumor_male.RData"))
load(paste0(work.dir,file.sep,"sesame_rgset_cord_female.RData"))
load(paste0(work.dir,file.sep,"sesame_rgset_cord_male.RData"))

k450_cn_methylset_normal_female <-
  minfi::getCN(minfi::preprocessRaw(sesame_rgset_normal_female))

k450_cn_methylset_normal_male <-
  minfi::getCN(minfi::preprocessRaw(sesame_rgset_normal_male))

k450_cn_methylset_tumor_female <- 
  minfi::getCN(minfi::preprocessRaw(sesame_rgset_tumor_female))

k450_cn_methylset_tumor_male <- 
  minfi::getCN(minfi::preprocessRaw(sesame_rgset_tumor_male))

k450_cn_methylset_cord_female <-
  minfi::getCN(minfi::preprocessRaw(sesame_rgset_cord_female))

k450_cn_methylset_cord_male <-
  minfi::getCN(minfi::preprocessRaw(sesame_rgset_cord_male))

save(k450_cn_methylset_normal_female, 
     file=paste0(work.dir,file.sep,"k450_cn_methylset_normal_female.RData"))

save(k450_cn_methylset_normal_male, 
     file=paste0(work.dir,file.sep,"k450_cn_methylset_normal_male.RData"))

save(k450_cn_methylset_tumor_female, 
     file=paste0(work.dir,file.sep,"k450_cn_methylset_tumor_female.RData"))

save(k450_cn_methylset_tumor_male, 
     file=paste0(work.dir,file.sep,"k450_cn_methylset_tumor_male.RData"))

save(k450_cn_methylset_cord_female, 
     file=paste0(work.dir,file.sep,"k450_cn_methylset_cord_female.RData"))

save(k450_cn_methylset_cord_male, 
     file=paste0(work.dir,file.sep,"k450_cn_methylset_cord_male.RData"))

## Choose a workflow, options are:
## A, B (Z-transform), or C (conumee)

switch(k450.workflow,
       
       A = {
         
         ## Without z-transformation, illumina
         ## For every tumor sample we are comparing to the median of 
         ## the normals
         
         proc_tumor_female_A  <- k450_cn_methylset_tumor_female[,
                                                                1:length(colnames(
                                                                  sesame_rgset_tumor_female))]
         
         proc_tumor_male_A    <- k450_cn_methylset_tumor_male[,
                                                              1:length(colnames(
                                                                sesame_rgset_tumor_male))]
         
         proc_normal_female_A <- k450_cn_methylset_normal_female[,
                                                                 1:length(colnames(
                                                                   sesame_rgset_normal_female))] ##4 Cols
         med_normal_female_A  <- apply(proc_normal_female_A, 1, "median")
         
         proc_normal_male_A   <- k450_cn_methylset_normal_male[,
                                                               1:length(colnames(
                                                                 sesame_rgset_normal_male))]
         med_normal_male_A    <- apply(proc_normal_male_A, 1, "median")
         
         proc_cord_female_A   <- k450_cn_methylset_cord_female[,
                                                               1:length(colnames(
                                                                 sesame_rgset_cord_female))]
         med_cord_female_A    <- apply(proc_cord_female_A, 1, "median")
         
         proc_cord_male_A     <- k450_cn_methylset_cord_male[,
                                                             1:length(colnames(
                                                               sesame_rgset_cord_male))]
         med_cord_male_A      <- apply(proc_cord_male_A, 1, "median")
         
         candidates_data_normal_female_A <-
           findSegments2(proc_tumor_female_A[, , drop = FALSE], 
                         med_normal_female_A, 
                         proc_normal_female_A) ##Scaled with B but not A
         candidates_data_normal_female_A$karyotype    <- ""
         candidates_data_normal_female_A$sex_reported <- ""
         candidates_data_normal_female_A$sex_inferred <- ""
         for(i in 1:length(unique(candidates_data_normal_female_A$smp))){
           sample.now <- unique(candidates_data_normal_female_A$smp)[i]
           if(sample.now %in% rownames(tumor)){
             candidates_data_normal_female_A[
               candidates_data_normal_female_A$smp==sample.now,"karyotype"] <- 
               tumor[sample.now, "karyotype"]
             candidates_data_normal_female_A[
               candidates_data_normal_female_A$smp==sample.now,"sex_reported"] <- 
               tumor[sample.now, "gender_reported"]
             candidates_data_normal_female_A[
               candidates_data_normal_female_A$smp==sample.now,"sex_inferred"] <- 
               tumor[sample.now, "sex_inferred"]
           }
         }
         
         candidates_data_normal_male_A <-
           findSegments2(proc_tumor_male_A[, , drop = FALSE], 
                         med_normal_male_A, 
                         proc_normal_male_A)
         candidates_data_normal_male_A$karyotype    <- ""
         candidates_data_normal_male_A$sex_reported <- ""
         candidates_data_normal_male_A$sex_inferred <- ""
         for(i in 1:length(unique(candidates_data_normal_male_A$smp))){
           sample.now <- unique(candidates_data_normal_male_A$smp)[i]
           if(sample.now %in% rownames(tumor)){
             candidates_data_normal_male_A[
               candidates_data_normal_male_A$smp==sample.now,"karyotype"] <- 
               tumor[sample.now, "karyotype"]
             candidates_data_normal_male_A[
               candidates_data_normal_male_A$smp==sample.now,"sex_reported"] <- 
               tumor[sample.now, "gender_reported"]
             candidates_data_normal_male_A[
               candidates_data_normal_male_A$smp==sample.now,"sex_inferred"] <- 
               tumor[sample.now, "sex_inferred"]
           }
         }
         
         candidates_data_cord_female_A <-
           findSegments2(proc_tumor_female_A[, , drop = FALSE], 
                         med_cord_female_A, 
                         proc_cord_female_A)
         candidates_data_cord_female_A$karyotype    <- ""
         candidates_data_cord_female_A$sex_reported <- ""
         candidates_data_cord_female_A$sex_inferred <- ""
         for(i in 1:length(unique(candidates_data_cord_female_A$smp))){
           sample.now <- unique(candidates_data_cord_female_A$smp)[i]
           if(sample.now %in% rownames(tumor)){
             candidates_data_cord_female_A[
               candidates_data_cord_female_A$smp==sample.now,"karyotype"] <- 
               tumor[sample.now, "karyotype"]
             candidates_data_cord_female_A[
               candidates_data_cord_female_A$smp==sample.now,"sex_reported"] <- 
               tumor[sample.now, "gender_reported"]
             candidates_data_cord_female_A[
               candidates_data_cord_female_A$smp==sample.now,"sex_inferred"] <- 
               tumor[sample.now, "sex_inferred"]
           }
         }
         
         candidates_data_cord_male_A <-
           findSegments2(proc_tumor_male_A[, , drop = FALSE], 
                         med_cord_male_A, 
                         proc_cord_male_A)
         candidates_data_cord_male_A$karyotype    <- ""
         candidates_data_cord_male_A$sex_reported <- ""
         candidates_data_cord_male_A$sex_inferred <- ""
         for(i in 1:length(unique(candidates_data_cord_male_A$smp))){
           sample.now <- unique(candidates_data_cord_male_A$smp)[i]
           if(sample.now %in% rownames(tumor)){
             candidates_data_cord_male_A[
               candidates_data_cord_male_A$smp==sample.now,"karyotype"] <- 
               tumor[sample.now, "karyotype"]
             candidates_data_cord_male_A[
               candidates_data_cord_male_A$smp==sample.now,"sex_reported"] <- 
               tumor[sample.now, "gender_reported"]
             candidates_data_cord_male_A[
               candidates_data_cord_male_A$smp==sample.now,"sex_inferred"] <- 
               tumor[sample.now, "sex_inferred"]
           }
         }
         
         save(candidates_data_normal_female_A, 
              file=paste0(work.dir,
                          file.sep,
                          "candidates_data_normal_female_A.RData"))
         
         save(candidates_data_normal_male_A, 
              file=paste0(work.dir,
                          file.sep,
                          "candidates_data_normal_male_A.RData"))
         
         save(candidates_data_cord_female_A, 
              file=paste0(work.dir,
                          file.sep,
                          "candidates_data_cord_female_A.RData"))
         
         save(candidates_data_cord_male_A, 
              file=paste0(work.dir,
                          file.sep,
                          "candidates_data_cord_male_A.RData"))
         
         rm(proc_normal_female_A)
         rm(proc_normal_male_A)
         rm(proc_tumor_female_A)
         rm(proc_tumor_male_A)
         rm(proc_cord_female_A)
         rm(proc_cord_male_A)
         rm(candidates_data_normal_female_A)
         rm(candidates_data_normal_male_A)
         rm(candidates_data_cord_female_A)
         rm(candidates_data_cord_male_A)
         
       },
       
       B = {
         
         ## With z-Transformation, illumina
         proc_normal_female_B  <- k450_cn_methylset_normal_female[,
                                                                  1:length(colnames(
                                                                    sesame_rgset_normal_female))]
         proc_normal_female_B[is.infinite(proc_normal_female_B)] <- NA
         proc_normal_female_B <- scale(proc_normal_female_B) ##Z transform
         med_normal_female_B <- apply(proc_normal_female_B, 1, "median")
         
         proc_normal_male_B <- k450_cn_methylset_normal_male[,
                                                             1:length(colnames(
                                                               sesame_rgset_normal_male))]
         proc_normal_male_B[is.infinite(proc_normal_male_B)] <- NA
         proc_normal_male_B <- scale(proc_normal_male_B)
         med_normal_male_B <- apply(proc_normal_male_B, 1, "median")
         
         proc_tumor_female_B  <- k450_cn_methylset_tumor_female[,
                                                                1:length(colnames(
                                                                  sesame_rgset_tumor_female))]
         proc_tumor_female_B[is.infinite(proc_tumor_female_B)] <- NA
         proc_tumor_female_B <- scale(proc_tumor_female_B)
         med_tumor_female_B <- apply(proc_tumor_female_B, 1, "median")
         
         proc_tumor_male_B  <- k450_cn_methylset_tumor_male[,
                                                            1:length(colnames(
                                                              sesame_rgset_tumor_male))]
         proc_tumor_male_B[is.infinite(proc_tumor_male_B)] <- NA
         proc_tumor_male_B <- scale(proc_tumor_male_B)
         med_tumor_male_B <- apply(proc_tumor_male_B, 1, "median")
         
         proc_cord_female_B <- k450_cn_methylset_cord_female[,
                                                             1:length(colnames(
                                                               sesame_rgset_cord_female))]
         proc_cord_female_B[is.infinite(proc_cord_female_B)] <- NA
         proc_cord_female_B <- scale(proc_cord_female_B)
         med_cord_female_B <- apply(proc_cord_female_B, 1, "median")
         
         proc_cord_male_B <- k450_cn_methylset_cord_male[,
                                                         1:length(colnames(
                                                           sesame_rgset_cord_male))]
         proc_cord_male_B[is.infinite(proc_cord_male_B)] <- NA
         proc_cord_male_B <- scale(proc_cord_male_B)
         med_cord_male_B <- apply(proc_cord_male_B, 1, "median")
         
         candidates_data_normal_female_B <-
           findSegments2(proc_tumor_female_B[, , drop = FALSE], 
                         med_normal_female_B, 
                         proc_normal_female_B) ##Scaled with B but not A
         candidates_data_normal_female_B$karyotype    <- ""
         candidates_data_normal_female_B$sex_reported <- ""
         candidates_data_normal_female_B$sex_inferred <- ""
         for(i in 1:length(unique(candidates_data_normal_female_B$smp))){
           sample.now <- unique(candidates_data_normal_female_B$smp)[i]
           if(sample.now %in% rownames(tumor)){
             candidates_data_normal_female_B[
               candidates_data_normal_female_B$smp==sample.now,"karyotype"] <- 
               tumor[sample.now, "karyotype"]
             candidates_data_normal_female_B[
               candidates_data_normal_female_B$smp==sample.now,"sex_reported"] <- 
               tumor[sample.now, "gender_reported"]
             candidates_data_normal_female_B[
               candidates_data_normal_female_B$smp==sample.now,"sex_inferred"] <- 
               tumor[sample.now, "sex_inferred"]
           }
         }
         
         candidates_data_normal_male_B <-
           findSegments2(proc_tumor_male_B[, , drop = FALSE], 
                         med_normal_male_B, 
                         proc_normal_male_B)
         candidates_data_normal_male_B$karyotype    <- ""
         candidates_data_normal_male_B$sex_reported <- ""
         candidates_data_normal_male_B$sex_inferred <- ""
         for(i in 1:length(unique(candidates_data_normal_male_B$smp))){
           sample.now <- unique(candidates_data_normal_male_B$smp)[i]
           if(sample.now %in% rownames(tumor)){
             candidates_data_normal_male_B[
               candidates_data_normal_male_B$smp==sample.now,"karyotype"] <- 
               tumor[sample.now, "karyotype"]
             candidates_data_normal_male_B[
               candidates_data_normal_male_B$smp==sample.now,"sex_reported"] <- 
               tumor[sample.now, "gender_reported"]
             candidates_data_normal_male_B[
               candidates_data_normal_male_B$smp==sample.now,"sex_inferred"] <- 
               tumor[sample.now, "sex_inferred"]
           }
         }
         
         candidates_data_cord_female_B <-
           findSegments2(proc_tumor_female_B[, , drop = FALSE], 
                         med_cord_female_B, 
                         proc_cord_female_B)
         candidates_data_cord_female_B$karyotype    <- ""
         candidates_data_cord_female_B$sex_reported <- ""
         candidates_data_cord_female_B$sex_inferred <- ""
         for(i in 1:length(unique(candidates_data_cord_female_B$smp))){
           sample.now <- unique(candidates_data_cord_female_B$smp)[i]
           if(sample.now %in% rownames(tumor)){
             candidates_data_cord_female_B[
               candidates_data_cord_female_B$smp==sample.now,"karyotype"] <- 
               tumor[sample.now, "karyotype"]
             candidates_data_cord_female_B[
               candidates_data_cord_female_B$smp==sample.now,"sex_reported"] <- 
               tumor[sample.now, "gender_reported"]
             candidates_data_cord_female_B[
               candidates_data_cord_female_B$smp==sample.now,"sex_inferred"] <- 
               tumor[sample.now, "sex_inferred"]
           }
         }
         
         candidates_data_cord_male_B <-
           findSegments2(proc_tumor_male_B[, , drop = FALSE], 
                         med_cord_male_B, 
                         proc_cord_male_B)
         candidates_data_cord_male_B$karyotype    <- ""
         candidates_data_cord_male_B$sex_reported <- ""
         candidates_data_cord_male_B$sex_inferred <- ""
         for(i in 1:length(unique(candidates_data_cord_male_B$smp))){
           sample.now <- unique(candidates_data_cord_male_B$smp)[i]
           if(sample.now %in% rownames(tumor)){
             candidates_data_cord_male_B[
               candidates_data_cord_male_B$smp==sample.now,"karyotype"] <- 
               tumor[sample.now, "karyotype"]
             candidates_data_cord_male_B[
               candidates_data_cord_male_B$smp==sample.now,"sex_reported"] <- 
               tumor[sample.now, "gender_reported"]
             candidates_data_cord_male_B[
               candidates_data_cord_male_B$smp==sample.now,"sex_inferred"] <- 
               tumor[sample.now, "sex_inferred"]
           }
         }
         
         save(candidates_data_normal_female_B, 
              file=paste0(work.dir,
                          file.sep,
                          "candidates_data_normal_female_B.RData"))
         
         save(candidates_data_normal_male_B, 
              file=paste0(work.dir,
                          file.sep,
                          "candidates_data_normal_male_B.RData"))
         
         save(candidates_data_cord_female_B, 
              file=paste0(work.dir,
                          file.sep,
                          "candidates_data_cord_female_B.RData"))
         
         save(candidates_data_cord_male_B, 
              file=paste0(work.dir,
                          file.sep,
                          "candidates_data_cord_male_B.RData"))
         
         rm(proc_normal_female_B)
         rm(proc_normal_male_B)
         rm(proc_tumor_female_B)
         rm(proc_tumor_male_B)
         rm(proc_cord_female_B)
         rm(proc_cord_male_B)
         rm(candidates_data_normal_female_B)
         rm(candidates_data_normal_male_B)
         rm(candidates_data_cord_female_B)
         rm(candidates_data_cord_male_B)
         
       },
       
       C = {
         
         ## Conumee-path, Illumina
         
         proc_normal_female_C <- k450_cn_methylset_normal_female[,
                                                                 1:length(colnames(
                                                                   sesame_rgset_normal_female))]
         
         proc_normal_male_C   <- k450_cn_methylset_normal_male[,
                                                               1:length(colnames(
                                                                 sesame_rgset_normal_male))]
         
         proc_tumor_female_C <- k450_cn_methylset_tumor_female[,
                                                               1:length(colnames(
                                                                 sesame_rgset_tumor_female))]
         
         proc_tumor_male_C   <- k450_cn_methylset_tumor_male[,
                                                             1:length(colnames(
                                                               sesame_rgset_tumor_male))]
         
         proc_cord_female_C  <- k450_cn_methylset_cord_female[,
                                                              1:length(colnames(
                                                                sesame_rgset_cord_female))]
         
         proc_cord_male_C    <- k450_cn_methylset_cord_male[,
                                                            1:length(colnames(
                                                              sesame_rgset_cord_male))]
         
         ## Conumee:
         ## What information is of interest?
         ## Add the genes of interest
         ## Can select which feature of interest we want to look at
         ## Calculate segments with conumee
         
         ## MM Note: from warnings() "All infinite values are set to NA!"
         ## MM Note: cord data is from epic and tumor data is 450k
         
         ##proc_tumor_female_C %>% head(n=10)
         ##proc_cord_female_C  %>% head(n=10)
         
         ##rownames(proc_tumor_female_C) %in% rownames(proc_cord_female_C)
         ## Quite a few probes are present, excellent. 
         ## Let's subset out the 450k from the EPIC cord ref 
         ## to use in the analysis.
         
         proc_cord_female_present_C <- 
           proc_cord_female_C[rownames(proc_cord_female_C) %in% 
                                rownames(proc_tumor_female_C),]
         
         proc_cord_male_present_C <- 
           proc_cord_male_C[rownames(proc_cord_male_C) %in% 
                              rownames(proc_tumor_male_C),]
         
         ############ Samples and ref also need to be in same order ##########
         ## "CpG probe IDs not in the same order in data and ctrl!"
         
         proc_tumor_female_sorted_C = proc_tumor_female_C[order(
           rownames(proc_tumor_female_C)),]
         
         proc_tumor_male_sorted_C   = proc_tumor_male_C[order(
           rownames(proc_tumor_male_C)),]
         
         proc_cord_female_sorted_C = proc_cord_female_C[order(
           rownames(proc_cord_female_C)),]
         
         proc_cord_male_sorted_C   = proc_cord_male_C[order(
           rownames(proc_cord_male_C)),]
         
         female_shared_names <- intersect(rownames(proc_cord_female_sorted_C),
                                          rownames(proc_tumor_female_sorted_C))
         
         male_shared_names   <- intersect(rownames(proc_cord_male_sorted_C),
                                          rownames(proc_tumor_male_sorted_C))
         
         proc_tumor_female_sorted_C <- 
           proc_tumor_female_sorted_C[female_shared_names,]
         
         proc_tumor_male_sorted_C   <- 
           proc_tumor_male_sorted_C[male_shared_names,]
         
         proc_cord_female_sorted_C  <- 
           proc_cord_female_sorted_C[female_shared_names,]
         
         proc_cord_male_sorted_C    <- 
           proc_cord_male_sorted_C[male_shared_names,]
         
         ## Check ordering now:
         ##proc_tumor_female_sorted_C  %>% head(n=10)
         ##proc_cord_female_sorted_C   %>% head(n=10)
         
         ####################### RUN CONUMEE ##############################
         
         candidates_data_normal_female_C <- 
           cnAnalysis450k::runConumee(proc_tumor_female_C, 
                                      proc_normal_female_C)
         candidates_data_normal_female_C$data$karyotype    <- ""
         candidates_data_normal_female_C$data$sex_reported <- ""
         candidates_data_normal_female_C$data$sex_inferred <- ""
         for(i in 1:length(unique(candidates_data_normal_female_C$data$ID))){
           sample.now <- unique(candidates_data_normal_female_C$data$ID)[i]
           if(sample.now %in% rownames(tumor)){
             candidates_data_normal_female_C$data[
               candidates_data_normal_female_C$data$ID==sample.now,
               "karyotype"] <- tumor[sample.now, "karyotype"]
             candidates_data_normal_female_C$data[
               candidates_data_normal_female_C$data$ID==sample.now,
               "sex_reported"] <- tumor[sample.now, "gender_reported"]
             candidates_data_normal_female_C$data[
               candidates_data_normal_female_C$data$ID==sample.now,
               "sex_inferred"] <- tumor[sample.now, "sex_inferred"]
           }
         }
         
         candidates_data_normal_male_C   <- 
           cnAnalysis450k::runConumee(proc_tumor_male_C, 
                                      proc_normal_male_C)
         candidates_data_normal_male_C$data$karyotype    <- ""
         candidates_data_normal_male_C$data$sex_reported <- ""
         candidates_data_normal_male_C$data$sex_inferred <- ""
         for(i in 1:length(unique(candidates_data_normal_male_C$data$ID))){
           sample.now <- unique(candidates_data_normal_male_C$data$ID)[i]
           if(sample.now %in% rownames(tumor)){
             candidates_data_normal_male_C$data[
               candidates_data_normal_male_C$data$ID==sample.now,
               "karyotype"] <- tumor[sample.now, "karyotype"]
             candidates_data_normal_male_C$data[
               candidates_data_normal_male_C$data$ID==sample.now,
               "sex_reported"] <- tumor[sample.now, "gender_reported"]
             candidates_data_normal_male_C$data[
               candidates_data_normal_male_C$data$ID==sample.now,
               "sex_inferred"] <- tumor[sample.now, "sex_inferred"]
           }
         }
         
         candidates_data_cord_female_C <- 
           cnAnalysis450k::runConumee(proc_tumor_female_sorted_C, 
                                      proc_cord_female_sorted_C)
         candidates_data_cord_female_C$data$karyotype    <- ""
         candidates_data_cord_female_C$data$sex_reported <- ""
         candidates_data_cord_female_C$data$sex_inferred <- ""
         for(i in 1:length(unique(candidates_data_cord_female_C$data$ID))){
           sample.now <- unique(candidates_data_cord_female_C$data$ID)[i]
           if(sample.now %in% rownames(tumor)){
             candidates_data_cord_female_C$data[
               candidates_data_cord_female_C$data$ID==sample.now,
               "karyotype"] <- tumor[sample.now, "karyotype"]
             candidates_data_cord_female_C$data[
               candidates_data_cord_female_C$data$ID==sample.now,
               "sex_reported"] <- tumor[sample.now, "gender_reported"]
             candidates_data_cord_female_C$data[
               candidates_data_cord_female_C$data$ID==sample.now,
               "sex_inferred"] <- tumor[sample.now, "sex_inferred"]
           }
         }
         
         candidates_data_cord_male_C   <- 
           cnAnalysis450k::runConumee(proc_tumor_male_sorted_C, 
                                      proc_cord_male_sorted_C)
         candidates_data_cord_male_C$data$karyotype    <- ""
         candidates_data_cord_male_C$data$sex_reported <- ""
         candidates_data_cord_male_C$data$sex_inferred <- ""
         for(i in 1:length(unique(candidates_data_cord_male_C$data$ID))){
           sample.now <- unique(candidates_data_cord_male_C$data$ID)[i]
           if(sample.now %in% rownames(tumor)){
             candidates_data_cord_male_C$data[
               candidates_data_cord_male_C$data$ID==sample.now,
               "karyotype"] <- tumor[sample.now, "karyotype"]
             candidates_data_cord_male_C$data[
               candidates_data_cord_male_C$data$ID==sample.now,
               "sex_reported"] <- tumor[sample.now, "gender_reported"]
             candidates_data_cord_male_C$data[
               candidates_data_cord_male_C$data$ID==sample.now,
               "sex_inferred"] <- tumor[sample.now, "sex_inferred"]
           }
         }
         
         save(candidates_data_normal_female_C, 
              file=paste0(work.dir,
                          file.sep,
                          "candidates_data_normal_female_C.RData"))
         
         save(candidates_data_normal_male_C, 
              file=paste0(work.dir,
                          file.sep,
                          "candidates_data_normal_male_C.RData"))
         
         save(candidates_data_cord_female_C, 
              file=paste0(work.dir,
                          file.sep,
                          "candidates_data_cord_female_C.RData"))
         
         save(candidates_data_cord_male_C, 
              file=paste0(work.dir,
                          file.sep,
                          "candidates_data_cord_male_C.RData"))
         
         rm(proc_normal_female_C)
         rm(proc_normal_male_C)
         rm(proc_tumor_female_C)
         rm(proc_tumor_male_C)
         rm(proc_cord_female_C)
         rm(proc_cord_male_C)
         rm(candidates_data_normal_female_C)
         rm(candidates_data_normal_male_C)
         rm(candidates_data_cord_female_C)
         rm(candidates_data_cord_male_C)
         
       }
       
)

rm(sesame_rgset_normal_female)
rm(sesame_rgset_normal_male)
rm(sesame_rgset_tumor_female)
rm(sesame_rgset_tumor_male)
rm(sesame_rgset_cord_female)
rm(sesame_rgset_cord_male)

rm(k450_cn_methylset_normal_female)
rm(k450_cn_methylset_normal_male)
rm(k450_cn_methylset_tumor_female)
rm(k450_cn_methylset_tumor_male)
rm(k450_cn_methylset_cord_female)
rm(k450_cn_methylset_cord_male)
