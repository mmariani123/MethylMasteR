#!/usr/bin/env Rscript

##Methyl Master Process Sesame 

#################### PROCESS SESAME #########################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################

setExperimentHubOption("CACHE", 
                       idat.pooled.files.dir)

ExperimentHub()

idat_prefixes <- searchIDATprefixes(idat.pooled.files.dir, 
                                    recursive=TRUE)

sesameDataCacheAll()

sesame_betas <- openSesame(idat_prefixes, 
                           mask = TRUE,
                           sum.TypeI = TRUE, 
                           platform = sesame.platform)

sesame_sset <- openSesame(idat_prefixes, 
                          mask = TRUE,
                          sum.TypeI = TRUE, 
                          platform = sesame.platform,
                          what="sigset")

sesame_rgset <- sesame::SigSetsToRGChannelSet(sesame_sset)

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

##rm(sset)

##length(idat_prefixes) ##30

##Below takes a long time to run (at least a couple hours on my desktop)
##Completed over night
##Warning message:
##  In `[.data.frame`(chrom.windows, min.window, "end") :
##  restarting interrupted promise evaluation
##Seems to have completed fine though

sesame_sset_tumor         <- sesame_sset[rownames(tumor)]
length(sesame_sset_tumor) ##10
names(sesame_sset_tumor)  <- rownames(tumor)

meta_tumor_female         <- tumor[tumor$sex_inferred == "FEMALE",]
sesame_sset_tumor_female        <- sesame_sset_tumor[
  rownames(meta_tumor_female)]
names(sesame_sset_tumor_female) <- rownames(meta_tumor_female)

meta_tumor_male           <- tumor[tumor$sex_inferred == "MALE",]
sesame_sset_tumor_male          <- sesame_sset_tumor[
  rownames(meta_tumor_male)]
names(sesame_sset_tumor_male)   <- rownames(meta_tumor_male)

sesame_sset_normal         <- sesame_sset[rownames(normal)]
length(sesame_sset_normal) ##10
names(sesame_sset_normal)  <- rownames(normal)

## subset sset into female and male for normal and tumor; 
## used to subset into sex using gender_reported, 
## need to be sex_inferred instead !!

meta_normal_female               <- normal[normal$sex_inferred == "FEMALE",]
sesame_sset_normal_female        <- sesame_sset_normal[
  rownames(meta_normal_female)]
names(sesame_sset_normal_female) <- rownames(meta_normal_female)

meta_normal_male          <- normal[normal$sex_inferred == "MALE",]
sesame_sset_normal_male   <- sesame_sset_normal[
  rownames(meta_normal_male)]
names(sesame_sset_normal_male)   <- rownames(meta_normal_male)

sesame_rgset_normal        <- SigSetsToRGChannelSet(sesame_sset_normal)
sesame_rgset_tumor         <- SigSetsToRGChannelSet(sesame_sset_tumor)
sesame_rgset_normal_female <- SigSetsToRGChannelSet(sesame_sset_normal_female)
sesame_rgset_normal_male   <- SigSetsToRGChannelSet(sesame_sset_normal_male)
sesame_rgset_tumor_female  <- SigSetsToRGChannelSet(sesame_sset_tumor_female)
sesame_rgset_tumor_male    <- SigSetsToRGChannelSet(sesame_sset_tumor_male)

save(sesame_betas,  file=paste0(work.dir,file.sep,"sesame_betas.RData"))
save(sesame_sset,   file=paste0(work.dir,file.sep,"sesame_sset.RData"))
save(sesame_rgset,  file=paste0(work.dir,file.sep,"sesame_rgset.RData"))
save(sesame_seg,    file=paste0(work.dir,file.sep,"sesame_seg.RData"))
save(sesame_qc,     file=paste0(work.dir,file.sep,"sesame_qc.RData"))

rm(sesame_betas)
rm(sesame_sset)
rm(sesame_rgset)
rm(sesame_seg)
rm(sesame_qc)

save(sesame_sset_tumor, 
     file=paste0(work.dir,file.sep,"sesame_sset_tumor.RData"))
save(sesame_sset_normal,           
     file=paste0(work.dir,file.sep,"sesame_sset_normal.RData"))
save(sesame_sset_normal_female,
     file=paste0(work.dir,file.sep,"sesame_sset_normal_female.RData"))
save(sesame_sset_normal_male  ,
     file=paste0(work.dir,file.sep,"sesame_sset_normal_male.RData"))
save(sesame_sset_tumor_female,  
     file=paste0(work.dir,file.sep,"sesame_sset_tumor_female.RData"))
save(sesame_sset_tumor_male ,   
     file=paste0(work.dir,file.sep,"sesame_sset_tumor_male.RData"))

save(sesame_rgset_normal, 
     file=paste0(work.dir,file.sep,"sesame_rgset_normal.RData"))
save(sesame_rgset_tumor,        
     file=paste0(work.dir,file.sep,"sesame_rgset_tumor.RData"))
save(sesame_rgset_normal_female,
     file=paste0(work.dir,file.sep,"sesame_rgset_normal_female.RData"))
save(sesame_rgset_normal_male  ,
     file=paste0(work.dir,file.sep,"sesame_rgset_normal_male.RData"))
save(sesame_rgset_tumor_female,
     file=paste0(work.dir,file.sep,"sesame_rgset_tumor_female.RData"))
save(sesame_rgset_tumor_male , 
     file=paste0(work.dir,file.sep,"sesame_rgset_tumor_male.RData"))

rm(sesame_sset_tumor)
rm(sesame_sset_normal)
rm(sesame_sset_tumor_male)
rm(sesame_sset_tumor_female)
rm(sesame_sset_normal_male)
rm(sesame_sset_normal_female)

rm(sesame_rgset_tumor)
rm(sesame_rgset_normal)
rm(sesame_rgset_tumor_male)
rm(sesame_rgset_tumor_female)
rm(sesame_rgset_normal_male)
rm(sesame_rgset_normal_female)


#################### PROCESS CORD #############################
###############################################################

setExperimentHubOption("CACHE", cord.files.path)

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

sesame_rgset_cord_female <- SigSetsToRGChannelSet(sesame_sset_cord_female)
sesame_rgset_cord_male   <- SigSetsToRGChannelSet(sesame_sset_cord_male)

save(sesame_sset_cord, 
     file=paste0(work.dir,
                 file.sep,
                 "sesame_sset_cord.RData"))
save(sesame_sset_cord_female, 
     file=paste0(work.dir,
                 file.sep,
                 "sesame_sset_cord_female.RData"))
save(sesame_sset_cord_male, 
     file=paste0(work.dir,
                 file.sep,
                 "sesame_sset_cord_male.RData"))
save(sesame_rgset_cord_female, 
     file=paste0(work.dir,
                 file.sep,
                 "sesame_rgset_cord_female.RData"))
save(sesame_rgset_cord_male, 
     file=paste0(work.dir,
                 file.sep,
                 "sesame_rgset_cord_male.RData"))

rm(sesame_sset_cord)
rm(sesame_sset_cord_female)
rm(sesame_sset_cord_male)
rm(sesame_rgset_cord_female)
rm(sesame_rgset_cord_male)

