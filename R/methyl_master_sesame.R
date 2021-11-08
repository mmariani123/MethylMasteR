#!/usr/bin/env Rscript

##Methyl Master Sesame

#################### INITIAL SESAME #########################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################

##Very similar to PROCESS SESAME command to prepare for segementation
##this is included so that the total time and mem for running sesame 
##can be recorded

methyl_master_sesame <- function(idat.pooled.files.dir,
                                 sample.sheet){

normal <- subset(sample.sheet.csv, Sample_Group=="normal")
tumor  <- subset(sample.sheet.csv, Sample_Group=="tumor")
  
rownames(normal) <- normal$Sample_Name
rownames(tumor)  <- tumor$Sample_Name
  
setExperimentHubOption("CACHE", 
                       idat.pooled.files.dir)

ExperimentHub()

idat_prefixes <- searchIDATprefixes(idat.pooled.files.dir, 
                                    recursive=TRUE)

sesameDataCacheAll()

##sesame_betas <- openSesame(idat_prefixes, 
##                           mask = TRUE,
##                           sum.TypeI = TRUE, 
##                           platform = sesame.platform)

sesame_sset <- openSesame(idat_prefixes, 
                          mask = TRUE,
                          sum.TypeI = TRUE, 
                          platform = sesame.platform,
                          what="sigset")

sesame_karyotype <- foreach(i = 1:length(names(sesame_sset))) %do% 
  {sesame::inferSexKaryotypes(sesame_sset[[i]])}

names(sesame_karyotype) <- names(sesame_sset)

sesame_karyotype <- as.data.frame(unlist(sesame_karyotype))

colnames(sesame_karyotype) <- "karyotype"

sesame_seg <- foreach(i = 1:length(names(sesame_sset))) %do% 
  {sesame::cnSegmentation(sesame_sset[[i]], 
                          sesame_ssets_normal, 
                          refversion = sesame.ref.version)}
names(sesame_seg) <- names(sesame_sset)

sesame_qc <- foreach(i = 1:length(names(sesame_sset))) %do% 
  {sesame::sesameQC(sesame_sset[[i]])}

names(sesame_qc) <- names(sesame_sset)

sesame_qc <- as.data.frame(data.table::rbindlist(sesame_qc))

rownames(sesame_qc) <- names(sesame_sset)

########## separate baseline_seg into normal and tumor #################

sesame_seg_normal_normal <- sesame_seg[names(sesame_seg) %in% normal$X]
sesame_seg_normal_tumor  <- sesame_seg[names(sesame_seg) %in% tumor$X]

rm(sesame_seg)

##Note the below are actual normal samples
save(sesame_seg_normal_normal, 
     file=paste0(work.dir,
                 file.sep,
                 "sesame_seg_normal_normal.RData"))
save(sesame_seg_normal_tumor, 
     file=paste0(work.dir,
                 file.sep,
                 "sesame_seg_normal_tumor.RData"))

##save(sesame_betas,  file=paste0(work.dir,file.sep,"sesame_betas.RData"))
##save(sesame_sset,   file=paste0(work.dir,file.sep,"sesame_sset.RData"))
##save(sesame_rgset,  file=paste0(work.dir,file.sep,"sesame_rgset.RData"))
##save(sesame_seg,    file=paste0(work.dir,file.sep,"sesame_seg.RData"))
##save(sesame_qc,     file=paste0(work.dir,file.sep,"sesame_qc.RData"))

##rm(sesame_betas)
##rm(sesame_sset)
##rm(sesame_rgset)
##rm(sesame_seg)
##rm(sesame_qc)

##save(sesame_sset_tumor, 
##     file=paste0(work.dir,file.sep,"sesame_sset_tumor.RData"))
##save(sesame_sset_normal,           
##     file=paste0(work.dir,file.sep,"sesame_sset_normal.RData"))
##save(sesame_sset_normal_female,
##     file=paste0(work.dir,file.sep,"sesame_sset_normal_female.RData"))
##save(sesame_sset_normal_male  ,
##     file=paste0(work.dir,file.sep,"sesame_sset_normal_male.RData"))
##save(sesame_sset_tumor_female,  
##     file=paste0(work.dir,file.sep,"sesame_sset_tumor_female.RData"))
##save(sesame_sset_tumor_male ,   
##     file=paste0(work.dir,file.sep,"sesame_sset_tumor_male.RData"))

##rm(sesame_sset_tumor)
##rm(sesame_sset_normal)
##rm(sesame_sset_tumor_male)
##rm(sesame_sset_tumor_female)
##rm(sesame_sset_normal_male)
##rm(sesame_sset_normal_female)

switch(sesame.ref,
       
       norm = {
         
         ############# segmentation: normal ref ##################################
         #########################################################################    
         
         ##load(paste0(work.dir, file.sep, "sesame_sset_normal_male.RData"  ))
         ##load(paste0(work.dir, file.sep, "sesame_sset_normal_female.RData"))
         ##load(paste0(work.dir, file.sep, "sesame_sset_tumor_male.RData"   ))
         ##load(paste0(work.dir, file.sep, "sesame_sset_tumor_female.RData" ))
         
         ##The below are tumor samples with the name "normal" because that 
         ##is the reference
         
         sesame_seg_normal_tumor_female  <- sesame_seg_normal_tumor[
           names(sesame_seg_normal_tumor) %in% 
             tumor$X[tumor$gender_reported=="female"]]
         for(i in 1:length(sesame_seg_normal_tumor_female)){
           sample.now <- names(sesame_seg_normal_tumor_female)[i]
           sesame_seg_normal_tumor_female[[i]]$seg.signals$karyotype    <- ""
           sesame_seg_normal_tumor_female[[i]]$seg.signals$sex_reported <- ""
           sesame_seg_normal_tumor_female[[i]]$seg.signals$sex_inferred <- ""
           if(sample.now %in% rownames(tumor)){
             sesame_seg_normal_tumor_female[[i]]$seg.signals$karyotype     <- 
               tumor[sample.now, "karyotype"]
             sesame_seg_normal_tumor_female[[i]]$seg.signals$sex_reported  <- 
               tumor[sample.now, "gender_reported"]
             sesame_seg_normal_tumor_female[[i]]$seg.signals$sex_inferred  <- 
               tumor[sample.now, "sex_inferred"]
           }
         }
         
         sesame_seg_normal_tumor_male <- sesame_seg_normal_tumor[
           names(sesame_seg_normal_tumor) %in% 
             tumor$X[tumor$gender_reported=="male"]]
         for(i in 1:length(sesame_seg_normal_tumor_male)){
           sample.now <- names(sesame_seg_normal_tumor_male)[i]
           sesame_seg_normal_tumor_male[[i]]$seg.signals$karyotype <- ""
           sesame_seg_normal_tumor_male[[i]]$seg.signals$sex_reported <- ""
           sesame_seg_normal_tumor_male[[i]]$seg.signals$sex_inferred <- ""
           if(sample.now %in% rownames(tumor)){
             sesame_seg_normal_tumor_male[[i]]$seg.signals$karyotype    <- 
               tumor[sample.now, "karyotype"]
             sesame_seg_normal_tumor_male[[i]]$seg.signals$sex_reported  <- 
               tumor[sample.now, "gender_reported"]
             sesame_seg_normal_tumor_male[[i]]$seg.signals$sex_inferred  <- 
               tumor[sample.now, "sex_inferred"]
           }
         }
         
         rm(sesame_seg_normal_tumor)
         
         save(sesame_seg_normal_tumor_female, 
              file=paste0(work.dir,
                          file.sep,
                          "sesame_seg_normal_tumor_female.RData"))
         
         save(sesame_seg_normal_tumor_male, 
              file=paste0(work.dir,
                          file.sep,
                          "sesame_seg_normal_tumor_male.RData"))
         
       },
       
       sp={
         
############# segmentation: tumor and sex paired  ref ###################
#########################################################################
         
         setExperimentHubOption("CACHE", 
                                idat.pooled.files.dir)
         
         ExperimentHub()
         
         idat_prefixes <- searchIDATprefixes(idat.pooled.files.dir, 
                                             recursive=TRUE)
         
         length(idat_prefixes)  ##30
         length(names(idat_prefixes)) ##30
         table(duplicated(idat_prefixes)) ##No duplicates
         table(names(idat_prefixes) %in% rownames(tumor))["TRUE"] ##10
         table(names(idat_prefixes) %in% rownames(normal))["TRUE"] ##10
         tumor.idat.samples  <- names(idat_prefixes)[names(idat_prefixes) %in% 
                                                       rownames(tumor)]
         normal.idat.samples <- names(idat_prefixes)[names(idat_prefixes) %in% 
                                                       rownames(normal)]
         
         intersect(tumor.idat.samples,normal.idat.samples) ##0, check
         
         idat_tumor_prefixes <- idat_prefixes[tumor.idat.samples]
         length(idat_tumor_prefixes) ##10
         
         sesameDataCacheAll()
         
         sesame_sset_tumor <- openSesame(idat_tumor_prefixes, 
                                         mask = TRUE,
                                         sum.TypeI = TRUE, 
                                         platform = sesame.platform,
                                         what="sigset")
         
         ##length(sesame_sset_tumor) ##10
         
         names(sesame_sset_tumor)  <- rownames(tumor)
         
         sesame_karyotype_tumor <- 
           foreach(i = 1:length(names(sesame_sset_tumor))) %do% 
           {sesame::inferSexKaryotypes(sesame_sset_tumor[[i]])}
         
         names(sesame_karyotype_tumor) <- names(sesame_sset_tumor)
         
         sesame_karyotype_tumor <- as.data.frame(unlist(sesame_karyotype_tumor))
         
         colnames(sesame_karyotype_tumor) <- "karyotype"
         
         sesame_qc_tumor <- foreach(i = 1:length(names(sesame_sset_tumor))) %do% 
           {sesame::sesameQC(sesame_sset_tumor[[i]])}
         
         names(sesame_qc_tumor) <- names(sesame_sset_tumor)
         
         sesame_qc_tumor <- as.data.frame(data.table::rbindlist(sesame_qc_tumor))
         
         rownames(sesame_qc_tumor) <- names(sesame_sset_tumor)
         
################################ NORMAL ######################################
         
         sesame_sset_normal <- openSesame(idat_normal_prefixes, 
                                          mask = TRUE,
                                          sum.TypeI = TRUE, 
                                          platform = sesame.platform,
                                          what="sigset")
         
         ##length(sesame_sset_normal) ##10
         
         names(sesame_sset_normal)  <- rownames(normal)
         
         sesame_karyotype_normal <- 
           foreach(i = 1:length(names(sesame_sset_normal))) %do% 
           {sesame::inferSexKaryotypes(sesame_sset_normal[[i]])}
         
         names(sesame_karyotype_normal) <- names(sesame_sset_normal)
         
         sesame_karyotype_normal <- as.data.frame(unlist(sesame_karyotype_normal))
         
         colnames(sesame_karyotype_normal) <- "karyotype"
         
         sesame_qc_normal <- foreach(i = 1:length(names(sesame_sset_normal))) %do% 
           {sesame::sesameQC(sesame_sset_normal[[i]])}
         
         names(sesame_qc_normal) <- names(sesame_sset_normal)
         
         sesame_qc_normal <- as.data.frame(data.table::rbindlist(sesame_qc_normal))
         
         rownames(sesame_qc_normal) <- names(sesame_sset_normal)
         
########################## Split into male and female #######################     
         
         ## subset sset into female and male for normal and tumor; 
         ## used to subset into sex using gender_reported, 
         ## need to be sex_inferred instead
         
         sesame_sset_tumor_female        <- sesame_sset_tumor[
           tumor[tumor$sex_inferred == "FEMALE",]]
         names(sesame_sset_tumor_female) <- 
           tumor[tumor$sex_inferred == "FEMALE",]
         
         sesame_sset_tumor_male          <- sesame_sset_tumor[
           tumor[tumor$sex_inferred == "MALE",]]
         names(sesame_sset_tumor_male)   <-  
           tumor[tumor$sex_inferred == "MALE",]
         
         sesame_sset_tumor_female        <- sesame_sset_tumor[
           tumor[tumor$sex_inferred == "FEMALE",]]
         names(sesame_sset_tumor_female) <- 
           tumor[tumor$sex_inferred == "FEMALE",]
         
         sesame_sset_tumor_male          <- sesame_sset_tumor[
           tumor[tumor$sex_inferred == "MALE",]]
         names(sesame_sset_tumor_male)   <-  
           tumor[tumor$sex_inferred == "MALE",]
         
         
######################### sesame segmentation ###########################
         
######################### female ########################################
         
         sesame_seg_paired_female <- foreach(i = 1:length(names(
           sesame_sset_tumor_female))) %do% 
           {sesame::cnSegmentation(sesame_sset_tumor_female[[i]], 
                                   sesame_sset_normal_female, 
                                   refversion = sesame.ref.version)}
         names(sesame_seg_paired_female) <- names(sesame_sset_tumor_female)
         for(i in 1:length(sesame_seg_paired_female)){
           sample.now <- names(sesame_seg_paired_female)[i]
           sesame_seg_paired_female[[i]]$seg.signals$karyotype <- ""
           sesame_seg_paired_female[[i]]$seg.signals$sex_reported <- ""
           sesame_seg_paired_female[[i]]$seg.signals$sex_inferred <- ""
           if(sample.now %in% rownames(tumor)){
             sesame_seg_paired_female[[i]]$seg.signals$karyotype    <- 
               tumor[sample.now, "karyotype"]
             sesame_seg_paired_female[[i]]$seg.signals$sex_reported  <- 
               tumor[sample.now, "gender_reported"]
             sesame_seg_paired_female[[i]]$seg.signals$sex_inferred  <- 
               tumor[sample.now, "sex_inferred"]
           }}
         
         save(sesame_seg_paired_female, 
              file=paste0(work.dir,
                          "\\",
                          "sesame_seg_paired_female.RData"))
         
         rm(sesame_seg_paired_female )
         rm(sesame_sset_normal_female)
         rm(sesame_sset_tumor_female )
         
######################### male #########################################
         
         sesame_seg_paired_male <- foreach(i = 1:length(
           names(sesame_sset_tumor_male))) %do% 
           {sesame::cnSegmentation(sesame_sset_tumor_male[[i]], 
                                   sesame_sset_normal_male, 
                                   refversion = sesame.ref.version)}
         names(sesame_seg_paired_male) <- names(sesame_sset_tumor_male)
         for(i in 1:length(sesame_seg_paired_male)){
           sample.now <- names(sesame_seg_paired_male)[i]
           sesame_seg_paired_male[[i]]$seg.signals$karyotype <- ""
           sesame_seg_paired_male[[i]]$seg.signals$sex_reported <- ""
           sesame_seg_paired_male[[i]]$seg.signals$sex_inferred <- ""
           if(sample.now %in% rownames(tumor)){
             sesame_seg_paired_male[[i]]$seg.signals$karyotype    <- 
               tumor[sample.now, "karyotype"]
             sesame_seg_paired_male[[i]]$seg.signals$sex_reported  <- 
               tumor[sample.now, "gender_reported"]
             sesame_seg_paired_male[[i]]$seg.signals$sex_inferred  <- 
               tumor[sample.now, "sex_inferred"]
           }}
         
         save(sesame_seg_paired_male, 
              file=paste0(work.dir,
                          file.sep,
                          "sesame_seg_paired_male.RData"))
         
         rm(sesame_seg_paired_male )
         rm(sesame_sset_normal_male)  
         rm(sesame_sset_tumor_male )
         
       },
       
       cord={
         
#################### PROCESS Tumor #############################
         
         setExperimentHubOption("CACHE", 
                                idat.pooled.files.dir)
         
         ExperimentHub()
         
         idat_prefixes <- searchIDATprefixes(idat.pooled.files.dir, 
                                             recursive=TRUE)
         
         length(idat_prefixes)  ##30
         length(names(idat_prefixes)) ##30
         table(duplicated(idat_prefixes)) ##No duplicates
         table(names(idat_prefixes) %in% rownames(tumor))["TRUE"] ##10
         table(names(idat_prefixes) %in% rownames(normal))["TRUE"] ##10
         tumor.idat.samples  <- names(idat_prefixes)[names(idat_prefixes) %in% 
                                                       rownames(tumor)]
         normal.idat.samples <- names(idat_prefixes)[names(idat_prefixes) %in% 
                                                       rownames(normal)]
         
         intersect(tumor.idat.samples,normal.idat.samples) ##0, check
         
         idat_tumor_prefixes <- idat_prefixes[tumor.idat.samples]
         length(idat_tumor_prefixes) ##10
         
         sesameDataCacheAll()
         
         sesame_sset_tumor <- openSesame(idat_tumor_prefixes, 
                                         mask = TRUE,
                                         sum.TypeI = TRUE, 
                                         platform = sesame.platform,
                                         what="sigset")
         
         ##length(sesame_sset_tumor) ##10
         
         names(sesame_sset_tumor)  <- rownames(tumor)
         
         sesame_karyotype_tumor <- 
           foreach(i = 1:length(names(sesame_sset_tumor))) %do% 
           {sesame::inferSexKaryotypes(sesame_sset_tumor[[i]])}
         
         names(sesame_karyotype_tumor) <- names(sesame_sset_tumor)
         
         sesame_karyotype_tumor <- as.data.frame(unlist(sesame_karyotype_tumor))
         
         colnames(sesame_karyotype_tumor) <- "karyotype"
         
         sesame_qc_tumor <- foreach(i = 1:length(names(sesame_sset_tumor))) %do% 
           {sesame::sesameQC(sesame_sset_tumor[[i]])}
         
         names(sesame_qc_tumor) <- names(sesame_sset_tumor)
         
         sesame_qc_tumor <- as.data.frame(data.table::rbindlist(sesame_qc_tumor))
         
         rownames(sesame_qc_tumor) <- names(sesame_sset_tumor)
         
         meta_tumor_female         <- tumor[tumor$sex_inferred == "FEMALE",]
         sesame_sset_tumor_female        <- sesame_sset_tumor[
           rownames(meta_tumor_female)]
         names(sesame_sset_tumor_female) <- rownames(meta_tumor_female)
         
         meta_tumor_male           <- tumor[tumor$sex_inferred == "MALE",]
         sesame_sset_tumor_male          <- sesame_sset_tumor[
           rownames(meta_tumor_male)]
         names(sesame_sset_tumor_male)   <- rownames(meta_tumor_male)
         
#################### PROCESS CORD #############################
         
         
         setExperimentHubOption("CACHE",
                                cord.files.path)
         
         ExperimentHub()
         
         cord_idat_prefixes <- searchIDATprefixes(cord.files.path,
                                                  recursive=TRUE)
         
         sesameDataCacheAll()  
         
         sesame_sset_cord <- openSesame(cord_idat_prefixes, 
                                        mask = TRUE,
                                        sum.TypeI = TRUE, 
                                        platform = sesame.cord.platform,
                                        what="sigset")
         
         ## Subset cord blood into female and male samples
         
         karyotype_cord <- NA
         
         karyotype_cord <- foreach(i = 1:length(names(sesame_sset_cord))) %do% 
           {sesame::inferSexKaryotypes(sesame_sset_cord[[i]])}
         
         sesame_female_cord_id <- character(length=0)
         sesame_male_cord_id   <- character(length=0)
         ##sesame_female_cord_id
         ##sesame_male_cord_id
         
         ##length(karyotype.cord) ##20
         for(i in 1:length(karyotype_cord)){
           print(karyotype_cord[i])
           if(as.character(karyotype_cord[i]) == "XaXi"){
             sesame_female_cord_id <- append(sesame_female_cord_id, 
                                             names(sesame_sset_cord)[i])
           }else{
             sesame_male_cord_id <- append(sesame_male_cord_id, 
                                           names(sesame_sset_cord)[i])
           }
         }
         
         ##print(sesame_female_cord_id) ##length = 4
         ##print(sesame_male_cord_id)   ##length = 6
         
         ##subset female sset:
         
         sesame_sset_cord_female <- NULL
         sesame_sset_cord_male   <- NULL
         
         for(i in 1:length(sesame_female_cord_id)){
           sesame_sset_cord_female <- append(sesame_sset_cord_female,
                                             sesame_sset_cord[[sesame_female_cord_id[i]]])
         }
         
         ##length(sesame_sset_cord_female) ##4
         names(sesame_sset_cord_female) <- sesame_female_cord_id
         
         for(i in 1:length(sesame_male_cord_id)){
           sesame_sset_cord_male <- append(sesame_sset_cord_male,
                                           sesame_sset_cord[[sesame_male_cord_id[i]]])
         }
         
         ##length(sesame_sset_cord_male) ##6
         names(sesame_sset_cord_male) <- sesame_male_cord_id
         
############################ female ##########################################
         
         sesame_seg_cord_female <- foreach(i = 1:length(names(
           sesame_sset_tumor_female))) %do% 
           {sesame::cnSegmentation(sesame_sset_tumor_female[[i]], 
                                   sesame_sset_cord_female, 
                                   refversion = sesame.ref.version)}
         names(sesame_seg_cord_female) <- names(sesame_sset_tumor_female)
         for(i in 1:length(sesame_seg_cord_female)){
           sample.now <- names(sesame_seg_cord_female)[i]
           sesame_seg_cord_female[[i]]$seg.signals$karyotype <- ""
           sesame_seg_cord_female[[i]]$seg.signals$sex_reported <- ""
           sesame_seg_cord_female[[i]]$seg.signals$sex_inferred <- ""
           if(sample.now %in% rownames(tumor)){
             sesame_seg_cord_female[[i]]$seg.signals$karyotype    <- 
               tumor[sample.now, "karyotype"]
             sesame_seg_cord_female[[i]]$seg.signals$sex_reported  <- 
               tumor[sample.now, "gender_reported"]
             sesame_seg_cord_female[[i]]$seg.signals$sex_inferred  <- 
               tumor[sample.now, "sex_inferred"]
           }}
         
         save(sesame_seg_cord_female, 
              file=paste0(work.dir,
                          file.sep,
                          "sesame_seg_cord_female.RData"))
         
         rm(sesame_sset_tumor_female)
         rm(sesame_sset_cord_female)
         rm(sesame_seg_cord_female)
         
########################### male #############################################
         
         sesame_seg_cord_male <- foreach(i = 1:length(names(
           sesame_sset_tumor_male))) %do% 
           {sesame::cnSegmentation(sesame_sset_tumor_male[[i]], 
                                   sesame_sset_cord_male, 
                                   refversion = sesame.ref.version)}
         names(sesame_seg_cord_male) <- names(sesame_sset_tumor_male)
         for(i in 1:length(sesame_seg_cord_male)){
           sample.now <- names(sesame_seg_cord_male)[i]
           sesame_seg_cord_male[[i]]$seg.signals$karyotype <- ""
           sesame_seg_cord_male[[i]]$seg.signals$sex_reported <- ""
           sesame_seg_cord_male[[i]]$seg.signals$sex_inferred <- ""
           if(sample.now %in% rownames(tumor)){
             sesame_seg_cord_male[[i]]$seg.signals$karyotype    <- 
               tumor[sample.now, "karyotype"]
             sesame_seg_cord_male[[i]]$seg.signals$sex_reported  <- 
               tumor[sample.now, "gender_reported"]
             sesame_seg_cord_male[[i]]$seg.signals$sex_inferred  <- 
               tumor[sample.now, "sex_inferred"]
           }}
         
         save(sesame_seg_cord_male, 
              file=paste0(work.dir,
                          file.sep,
                          "sesame_seg_cord_male.RData"))
         
         rm(sesame_sset_tumor_male)  
         rm(sesame_sset_cord_male)
         rm(sesame_seg_cord_male)
         
       }
)

}
