#!/usr/bin/env Rscript

switch(routine,
       
       ## 0: homozygous deletion (2-copy loss)
       ## 1: heterozygous deletion (1-copy loss)
       ## 2: normal diploid state
       ## 3: 1-copy gain
       ## 4: amplification (>= 2-copy gain)
       
       ##Overall the pipeline results structure is
       ##a dataframe with the following fields:
       ##"Sample_ID"        
       ##"chrom"     
       ##"loc.start" 
       ##"loc.end"   
       ##"num.mark"  
       ##"seg.mean"
       ##"state"
       ##"method"
       
       ####################### Load SeSAMe data: ############################
       
       sesame = {
         
         load(paste0(work.dir,
                     file.sep,
                     "sesame_seg_normal_tumor_female.RData"))
         load(paste0(work.dir,
                     file.sep,
                     "sesame_seg_normal_tumor_male.RData"))
         load(paste0(work.dir,
                     file.sep,
                     "sesame_seg_paired_female.RData"))
         load(paste0(work.dir,
                     file.sep,
                     "sesame_seg_paired_male.RData"))
         load(paste0(work.dir,
                     file.sep,
                     "sesame_seg_cord_female.RData"))
         load(paste0(work.dir,
                     file.sep,
                     "sesame_seg_cord_male.RData"))
         
         sesame_seg_normal_tumor_female <- binding_frames_mm(
           sesame_seg_normal_tumor_female)
         sesame_seg_normal_tumor_female <- sesame_seg_normal_tumor_female[
           sesame_seg_normal_tumor_female$pval <= 0.05,]
         
         sesame_seg_normal_tumor_male   <- binding_frames_mm(
           sesame_seg_normal_tumor_male)
         sesame_seg_normal_tumor_male   <- sesame_seg_normal_tumor_female[
           sesame_seg_normal_tumor_male$pval <= 0.05,]
         
         sesame_paired_female <- binding_frames_mm(sesame_seg_paired_female)
         sesame_paired_female <- sesame_paired_female[
           sesame_paired_female$pval <= 0.05,]
         
         sesame_paired_male   <- binding_frames_mm(sesame_seg_paired_male)
         sesame_paired_male   <- sesame_paired_female[
           sesame_paired_male$pval <= 0.05,]
         
         sesame_cord_female   <- binding_frames_mm(sesame_seg_cord_female)
         sesame_cord_female   <- sesame_cord_female[
           sesame_cord_female$pval <= 0.05,]
         
         sesame_cord_male     <- binding_frames_mm(sesame_seg_cord_male)
         sesame_cord_male     <- sesame_cord_male[
           sesame_cord_male$pval <= 0.05,]
         
         colnames(sesame_cord_male)
         seg <- sesame_cord_male[,c(1,2,3,4,5,6,13,14,15)]
         colnames(seg)
         colnames(seg)[1] <- "Sample_ID"
         seg$state <- round(2^seg$seg.mean * 2)
         seg$state[seg$state > 4] <- 4
         seg$method <- "sesame"
         row.names(seg) <- NULL
         seg <- na.omit(seg)
         ##For sesame we want to omit the "chr" from chrom for intial plotting
         ##then add back in downstream for ggplot
         seg$chrom <- 
           unlist(strsplit(seg$chrom,split="chr"))[c(FALSE,TRUE)]
         
       },
       
       ####################### Load 450K data: ##############################
       
       k450 = {
         
         if(k450.workflow=="A"){
           
           ##sub routine A
           load(paste0(work.dir,file.sep,"candidates_data_normal_female_A.RData"))
           load(paste0(work.dir,file.sep,"candidates_data_normal_male_A.RData"))
           load(paste0(work.dir,file.sep,"candidates_data_cord_female_A.RData"))
           load(paste0(work.dir,file.sep,"candidates_data_cord_male_A.RData"))
           
           candidates_data_normal_female_A_sig <- 
             candidates_data_normal_female_A[
               candidates_data_normal_female_A$p.val <= 0.05,]
           candidates_data_normal_male_A_sig   <- 
             candidates_data_normal_male_A[
               candidates_data_normal_male_A$p.val <= 0.05,]
           candidates_data_cord_female_A_sig   <- 
             candidates_data_cord_female_A[
               candidates_data_cord_female_A$p.val <= 0.05,]
           candidates_data_cord_male_A_sig   <- 
             candidates_data_cord_male_A[
               candidates_data_cord_male_A$p.val <= 0.05,]
           
           candidates_data_normal_female_A_sig$num.mark <- NA
           candidates_data_normal_female_A_sig$bstat    <- NA
           candidates_data_normal_female_A_sig$treatment <- "tumor"
           candidates_data_normal_male_A_sig$num.mark <- NA
           candidates_data_normal_male_A_sig$bstat    <- NA
           candidates_data_normal_male_A_sig$treatment <- "tumor"
           candidates_data_cord_female_A_sig$num.mark <- NA
           candidates_data_cord_female_A_sig$bstat    <- NA
           candidates_data_cord_female_A_sig$treatment <- "tumor"
           candidates_data_cord_male_A_sig$num.mark <- NA
           candidates_data_cord_male_A_sig$bstat    <- NA
           candidates_data_cord_male_A_sig$treatment <- "tumor"
           
           preferred.columns <- c(7,1,2,3,12,13,8,5,4,9,10,11,14)
           
           candidates_data_normal_female_A_sig <- 
             candidates_data_normal_female_A_sig[, preferred.columns]
           candidates_data_normal_male_A_sig <- 
             candidates_data_normal_male_A_sig[, preferred.columns]
           candidates_data_cord_female_A_sig <- 
             candidates_data_cord_female_A_sig[, preferred.columns]
           candidates_data_cord_male_A_sig <- 
             candidates_data_cord_male_A_sig[, preferred.columns]
           colnames(candidates_data_normal_female_A_sig) <-  c("ID", "chrom", "loc.start", "loc.end", "num.mark", "bstat", "pval", "seg.mean", "seg.median", "karyotype", "sex_reported", "sex_inferred", "treatment")
           colnames(candidates_data_normal_male_A_sig) <-  c("ID", "chrom", "loc.start", "loc.end", "num.mark", "bstat", "pval", "seg.mean", "seg.median", "karyotype", "sex_reported", "sex_inferred", "treatment")
           colnames(candidates_data_cord_female_A_sig) <-  c("ID", "chrom", "loc.start", "loc.end", "num.mark", "bstat", "pval", "seg.mean", "seg.median", "karyotype", "sex_reported", "sex_inferred", "treatment")
           colnames(candidates_data_cord_male_A_sig) <-  c("ID", "chrom", "loc.start", "loc.end", "num.mark", "bstat", "pval", "seg.mean", "seg.median", "karyotype", "sex_reported", "sex_inferred", "treatment")
           
           seg <- do.call(rbind,
                          list(candidates_data_normal_female_A_sig,  
                               candidates_data_normal_male_A_sig,    
                               candidates_data_cord_female_A_sig,  
                               candidates_data_cord_male_A_sig))
           
         }else if(k450.workflow==B){
           
           load(paste0(work.dir,file.sep,"candidates_data_normal_female_B.RData"))
           load(paste0(work.dir,file.sep,"candidates_data_normal_male_B.RData"))
           load(paste0(work.dir,file.sep,"candidates_data_cord_female_B.RData"))
           load(paste0(work.dir,file.sep,"candidates_data_cord_male_B.RData"))
           
           candidates_data_normal_female_B_sig <- 
             candidates_data_normal_female_B[
               candidates_data_normal_female_B$p.val <= 0.05,]
           candidates_data_normal_male_B_sig  <- 
             candidates_data_normal_male_B[
               candidates_data_normal_male_B$p.val <= 0.05,]
           candidates_data_cord_female_B_sig  <- 
             candidates_data_cord_female_B[
               candidates_data_cord_female_B$p.val <= 0.05,]
           candidates_data_cord_male_B_sig  <- 
             candidates_data_cord_male_B[
               candidates_data_cord_male_B$p.val <= 0.05,]
           
           candidates_data_normal_female_B_sig$num.mark <- NA
           candidates_data_normal_female_B_sig$bstat    <- NA
           candidates_data_normal_female_B_sig$treatment <- "tumor"
           candidates_data_normal_male_B_sig$num.mark <- NA
           candidates_data_normal_male_B_sig$bstat    <- NA
           candidates_data_normal_male_B_sig$treatment <- "tumor"
           candidates_data_cord_female_B_sig$num.mark <- NA
           candidates_data_cord_female_B_sig$bstat    <- NA
           candidates_data_cord_female_B_sig$treatment <- "tumor"
           candidates_data_cord_male_B_sig$num.mark <- NA
           candidates_data_cord_male_B_sig$bstat    <- NA
           candidates_data_cord_male_B_sig$treatment <- "tumor"
           
           candidates_data_normal_female_B_sig <- 
             candidates_data_normal_female_B_sig[, preferred.columns]
           candidates_data_normal_male_B_sig <- 
             candidates_data_normal_male_B_sig[, preferred.columns]
           candidates_data_cord_female_B_sig <- 
             candidates_data_cord_female_B_sig[, preferred.columns]
           candidates_data_cord_male_B_sig <- 
             candidates_data_cord_male_B_sig[, preferred.columns]
           colnames(candidates_data_normal_female_B_sig) <-  c("ID", "chrom", "loc.start", "loc.end", "num.mark", "bstat", "pval", "seg.mean", "seg.median", "karyotype", "sex_reported", "sex_inferred", "treatment")
           colnames(candidates_data_normal_male_B_sig) <-  c("ID", "chrom", "loc.start", "loc.end", "num.mark", "bstat", "pval", "seg.mean", "seg.median", "karyotype", "sex_reported", "sex_inferred", "treatment")
           colnames(candidates_data_cord_female_B_sig) <-  c("ID", "chrom", "loc.start", "loc.end", "num.mark", "bstat", "pval", "seg.mean", "seg.median", "karyotype", "sex_reported", "sex_inferred", "treatment")
           colnames(candidates_data_cord_male_B_sig) <-  c("ID", "chrom", "loc.start", "loc.end", "num.mark", "bstat", "pval", "seg.mean", "seg.median", "karyotype", "sex_reported", "sex_inferred", "treatment")
           
           seg <- do.call(rbind,
                          list(candidates_data_normal_female_B_sig,  
                               candidates_data_normal_male_B_sig,    
                               candidates_data_cord_female_B_sig,   
                               candidates_data_cord_male_B_sig))
           
         }else if(k450.workflow=="C"){
           
           load(paste0(work.dir,file.sep,"candidates_data_normal_female_C.RData"))
           load(paste0(work.dir,file.sep,"candidates_data_normal_male_C.RData"))
           load(paste0(work.dir,file.sep,"candidates_data_cord_female_C.RData"))
           load(paste0(work.dir,file.sep,"candidates_data_cord_male_C.RData"))
           
           candidates_data_normal_female_C_sig <- 
             candidates_data_normal_female_C$data[
               candidates_data_normal_female_C$data$pval <= 0.05,]
           candidates_data_normal_male_C_sig    <- 
             candidates_data_normal_male_C$data[
               candidates_data_normal_male_C$data$pval <= 0.05,]
           candidates_data_cord_female_C_sig    <- 
             candidates_data_cord_female_C$data[
               candidates_data_cord_female_C$data$pval <= 0.05,]
           candidates_data_cord_male_C_sig    <- 
             candidates_data_cord_male_C$data[
               candidates_data_cord_male_C$data$pval <= 0.05,]
           
           candidates_data_normal_female_C_sig$treatment <- "tumor"
           candidates_data_normal_male_C_sig$treatment <- "tumor"
           candidates_data_cord_female_C_sig$treatment <- "tumor"
           candidates_data_cord_male_C_sig$treatment <- "tumor"
           
           seg <- do.call(rbind,
                          list(candidates_data_normal_female_C_sig,  
                               candidates_data_normal_male_C_sig,    
                               candidates_data_cord_female_C_sig,   
                               candidates_data_cord_male_C_sig))
           
         }else{
           
           stop("Error: need to select a proper sub workflow for 450k <k450.workflow>")
           
         }
         
         ##Note differences in 450k standard and conumee:
         ##colnames(candidates_data_cord_male_B_sig)
         ##"chr"          
         ##"startCG"      
         ##"endCG"        
         ##"median"       
         ##"mean"         
         ##"sd"          
         ##"smp"
         ##"p.val"
         ##"karyotype"
         ##"sex_reported"
         ##"sex_inferred"
         ##"num.mark"
         ##"bstat"
         
         ##candidates_data_cord_male_C_sig
         ##"ID"
         ##"chrom"
         ##"loc.start"
         ##"loc.end"
         ##"num.mark"
         ##"bstat"       
         ##"pval"
         ##"seg.mean"
         ##"seg.median"
         ##"karyotype"
         ##"sex_reported"
         ##"sex_inferred"
         
         ##Desired names:
         ##colnames(seg)
         ##"ID"         
         ##"chrom"
         ##"loc.start"
         ##"loc.end"
         ##"num.mark"
         ##"bstat"     
         ##"pval"
         ##"seg.mean"
         ##"seg.median"
         ##karyotype
         ##sex_reported
         ##sex_inferred
         ##treatment
         
         colnames(seg)[1] <- "Sample_ID"
         seg <- seg[,c(1,2,3,4,5,8,10,11,12,13)]
         seg$state <- round(2^seg$seg.mean * 2)
         seg$state[seg$state > 4] <- 4
         seg$method <- "k450"
         seg$sub.method <- sub.workflow
         row.names(seg) <- NULL
         seg <- na.omit(seg) ##Workflow C ends up with some NA rows
         
       },
       
       ####################### Load champ data: ##############################
       
       champ = {
         
         load(file=paste0(work.dir, file.sep, "champ_seg.RData"))
         seg <- champ_seg
         rm(champ_seg)
         ##colnames(seg)
         ##"Sample_ID"        
         ##"chrom"     
         ##"loc.start" 
         ##"loc.end"   
         ##"num.mark"  
         ##"seg.mean"
         seg$state <- round(2^seg$seg.mean * 2)
         seg$state[seg$state > 4] <- 4
         seg$method <- "champ"
         row.names(seg) <- NULL
         seg <- na.omit(seg)
         ##For k450 we want to omit the "chr" from chrom for intial plotting
         ##then add back in downstream for ggplot
         
         seg$chrom <- 
           as.integer(unlist(strsplit(seg$chrom,split="chr"))[c(FALSE,TRUE)])
         
       },
       
       ####################### Load epicopy data: ##############################
       
       epicopy = {
         
         load(paste0(work.dir,
                     file.sep,
                     "epicopy_results.rda"))
         
         seg <- epicopy_results$output
         ##rm(epicopy_results)
         ##colnames(seg)
         ##"Sample_ID"        
         ##"chrom"     
         ##"loc.start" 
         ##"loc.end"   
         ##"num.mark"  
         ##"seg.mean"
         seg$state <- round(2^seg$seg.mean * 2)
         seg$state[seg$state > 4] <- 4
         seg$method <- "epicopy"
         row.names(seg) <- NULL
         seg <- na.omit(seg)
         
         ##seg2 <- epi_seg$output
         ##seg3 <- seg2[abs(seg2$seg.mean) >= 0.3 & 
         ##               !seg2$chrom %in% c("chrX", "chrY") & 
         ##               seg2$num.mark > 4,]
         ##seg3$state<- round(2^seg3$seg.mean * 2)
         ##seg3$state[seg3$state > 4] <- 4
         ##table(seg3$state)
         ##0  1  2  3  4 
         ##1 18 11 41  2 
         ##seg3$Sample_ID <- seg3$ID
         ##seg3$method <- "epicopy"
         ##colnames(seg)
         ##seg <- seg3[, c("Sample_ID", 
         ##                "chrom", 
         ##                "loc.start",   
         ##                "loc.end", 
         ##                "num.mark",
         ##                "seg.mean",
         ##                "state",
         ##                "method")]
         
         table(seg$state)
         ##0    1    2    3    4 
         ##72 2954 7118 1803  325 
         
         colnames(seg)
         ##"Sample_ID"        
         ##"chrom"     
         ##"loc.start" 
         ##"loc.end"   
         ##"num.mark"  
         ##"seg.mean"
         ##"state"
         ##"method"
         
         ##seg %>% dplyr::pull(Sample_ID) %>% unique()
         
       })

################# OVERLAPS AND VISUALIZE ####################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################

##samples.list <- as.list(c(paste0("normal_",seq(1,10)),
##                            paste0("tumor_",seq(1,10))))

##seg$chrom <- as.numeric(seg$chrom)
##seg <- na.omit(seg)
##changing the chrom to number and removing NAs from X and Y doesn't seem
##to fix the population ranges problem. 

grl <- GenomicRanges::makeGRangesListFromDataFrame(seg, 
                                                   split.field="Sample_ID",
                                                   keep.extra.columns=TRUE)

grl <- GenomicRanges::sort(grl)

##ra  <- RaggedExperiment::RaggedExperiment(grl)

##Excluding 997367 copy-number neutral regions (CN state = 2, diploid)
##Sign-in-by error below with usual function
##cnvrs <- CNVRanger::populationRanges(grl, 
##                                     density=0.1, 
##                                     est.recur=TRUE)

##RaggedExperiment::assay(ra[1:5,1:5])

##Trims low-density areas (usually <10% 
##of the total contributing individual 
##calls within a summarized region).

cnvrs <- salas_mm_population_ranges(grl, 
                                    ##density=0.1, 
                                    ##We can also like the density
                                    ##Jenn likes .01
                                    density=.1,
                                    est.recur=TRUE)
##cnvrs$type %>% unique()##
cnvrs.filt <- subset(cnvrs, pvalue < 0.05)

ra  <- RaggedExperiment::RaggedExperiment(grl)

##RaggedExperiment::assay(ra[1:5,1:5])

##below the less filtered ragged experiment
if(less.stringent.ra.setting==TRUE){
  cvns.matrix <- RaggedExperiment::assay(ra)
}else{
  ##This will likley not work as well with 
  ##fewer samples so can default to the 
  ##above, both return overlaps. 
  cvns.matrix <- 
    RaggedExperiment::qreduceAssay(ra, 
                                   cnvrs.filt, 
                                   simplifyReduce = weightedmean)
}

cvns.matrix <- cvns.matrix[order(rownames(cvns.matrix)),]
cvns.matrix[is.na(cvns.matrix)] <- 2
cvns.matrix <- round(cvns.matrix, 0)

write.table(cvns.matrix,
            file = paste0(work.dir,
                          file.sep,
                          routine,
                          "_",
                          sub.workflow,
                          "_",
                          ref,
                          "_overlaps.csv"),
            sep=",",
            col.names=TRUE,
            row.names=TRUE,
            quote=FALSE)

pheatmap.out <- pheatmap::pheatmap(cvns.matrix, 
                                   cluster_rows = T, 
                                   cluster_cols = T, 
                                   show_colnames = T)

ggsave(pheatmap.out,
       file = paste0(work.dir,
                     file.sep,
                     routine,
                     "_",
                     sub.workflow,
                     "_",
                     ref,
                     "_pheatmap.pdf"),
       device="pdf",
       width=12,
       height=8)

##Also from epicopy:
##cvnsmatrix <- cvnsmatrix[order(rownames(cvnsmatrix)),]
##cvnsmatrix[1:5,1:5]
##write.csv(cvnsmatrix, file=paste0(work.dir,
##                                  file.sep,
##                                  "cvnsmatrix_epicopy.csv"))
##cvnsmatrix2 <- cvnsmatrix
##cvnsmatrix2[is.na(cvnsmatrix2)] <- 2
##cvnsmatrix2 <- round(cvnsmatrix2,0)
##identical(rownames(pheno), colnames(cvnsmatrix))
##write.csv(pheno, file=paste0(work.dir,
##                             "\\",
##                             "pheno_sesame.CSV"))
##pheatmap::pheatmap(cvns.matrix, 
##                   cluster_rows = F, 
##                   cluster_cols = F, 
##                   show_colnames = F)
##pheatmap::pheatmap(cvnsmatrix2, 
##                   show_colnames = F)

########################## Heatmaps #########################################

##this is ~0.33 feber et al calc. and used with tcga
##dont forget the glioblastoma paper and parameters:

##Below adds "status field" to seg frame for ggplot heatmap

##Seg state calculation:
##frame.in$state <- round(2^frame.in$seg.mean * 2)
##frame.in$state[frame.in$state > 4] <- 4
##frame.in$method <- ""
##frame.in$Sample_ID <- sample.in

seg$status <- ifelse(seg$state>2,
                     "gain",
                     ifelse(seg$state<2,
                            "loss",
                            "neutral"))

seg$chrom  <- paste0("chr",seg$chrom)
seg$chrom  <- factor(seg$chrom, 
                     ##levels=gtools::mixedsort(
                     ##unique(ex.data.in$chrom))
                     levels=c("chr1", 
                              "chr2", 
                              "chr3",  
                              "chr4",
                              "chr5", 
                              "chr6", 
                              "chr7",  
                              "chr8", 
                              "chr9", 
                              "chr10",
                              "chr11", 
                              "chr12", 
                              "chr13", 
                              "chr14", 
                              "chr15",
                              "chr16", 
                              "chr17",
                              "chr18", 
                              "chr19",
                              "chr20",
                              "chr21", 
                              "chr22",
                              "chrX",
                              "chrY"))

seg$Sample_ID <- 
  factor(seg$Sample_ID,
         levels=unique(gtools::mixedsort(seg$Sample_ID)))

seg$status <- factor(seg$status,
                     levels = c("loss",
                                "neutral",
                                "gain"))

##floor(seg.vis.total$loc.start + seg.vis.total$loc.end) /2)
seg.heatmap <- ggplot(seg,
                      aes(x=as.integer((floor((loc.start + loc.end)/2))),
                          ##x=chrom,
                          ##y=seg.mean, 
                          y=0.5,
                          width=loc.start - loc.end,
                          fill=status)) +
  geom_tile() +
  ##geom_vline(xintercept=0) +
  facet_grid(rows=vars(Sample_ID),
             cols=vars(chrom),
             scales = "free",
             space = "free",
             switch = "y") +
  theme_bw() +
  theme(panel.spacing = unit(0,"cm"),
        axis.ticks.x  = element_blank(),
        axis.ticks.y  = element_blank(),
        axis.text.x   = element_blank(),
        axis.text.y   = element_blank(),
        axis.text.x.bottom  = element_blank(),
        axis.text.y.left = element_blank(),
        axis.title.y.left = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.text.y.left = element_text(angle = 0),
        strip.text.x.top = element_text(angle = 90)) +
  ylim(c(0,1)) +
  scale_fill_manual(values=c("orange", 
                             "white", 
                             "blue"))

ggsave(seg.heatmap,
       filename = paste0(work.dir,
                         file.sep,
                         routine,
                         "_",
                         sub.workflow,
                         "_",
                         ref,
                         "_mmheatmap.pdf"),
       width=12,
       height=8,
       device="pdf")

#################### VISUALISE INDIVIDUAL ######################################

if(visualize.individual==TRUE){
  if(!dir.exists(sesame.plots.dir)){
    dir.create(sesame.plots.dir)
  }else{
    print("sesame plots dir exists, overwriting")
    unlink(sesame.plots.dir,
           recursive=TRUE,
           force=TRUE)
    dir.create(sesame.plots.dir)
  }
  ##visualize individual sesame by sample:
  for(i in 1:length(sesame_seg_cord_male)){
    sample.now <- names(sesame_seg_cord_male)[i]
    ggsave({sesame::visualizeSegments(sesame_seg_cord_male[[1]]) + 
        theme_bw() +
        ylab("Signal") +
        labs(color="Sesame seg. \nsignal") +
        ggtitle(paste0("Sample ",
                       sample.now,
                       "\nRef = ",
                       ref)) +
        theme(axis.text.x = element_text(angle=90),
              plot.title = element_text(hjust=0.5))},
        filename = paste0(sesame.plots.dir,
                          file.sep,
                          sample.now,
                          "_",
                          ref,
                          "_sesame_plot.pdf"),
        width=12,
        height=8,
        device="pdf")
  } 
  
} ##End visualize individual

######################### End visualize ######################################
