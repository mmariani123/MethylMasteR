#!/usr/bin/env Rscript

##Michael Mariani Dartmouth College 2021

################## Multi Methyl extra and notes ###############################

##Pipeline to analyze and compare the CNV calling 
##algorithms that use methylation data

##Convert column names so matches sesame dataframe column names
##candidates_data_normal_male_x$seg.mean  <- candidates_data_normal_male_x$mean
##candidates_data_cord_female_x$pval      <- candidates_data_cord_female_x$p.val
##candidates_data_cord_female_x$chrom     <- candidates_data_cord_female_x$chr
##candidates_data_norm_male_x$Sample_ID   <- candidates_data_normal_male_x$smp

##Note:
##https://github.com/Bioconductor/RaggedExperiment/issues/10

##Note:
##epicopy internal normals are: 'thyroid', 'breast', 'lung', or 'all',

##Not in use currently:
##For visualizing cuttoffs thresholds: Can be "man", "auto" or "skip"
##visualize.cutoffs="skip"
##Visualize K450_WORKFLOW, Can be "A" , "B" , or "C" 
##visualize.workflow="C"

## Which routine would you like to run?
## "test"            ##Run a quick test
## "download" ,      ##For downloading TCGA data
## "process_sesame", ##Preprocess the TCGA and cord data in sesame format
## "sesame",         ##Run Sesame  CNV calling (get segments)
## "k450" ,          ##Run 450K    CNV calling (get segments)
## "champ" ,         ##Run ChAMP   CNV calling (get segments)
## "epicopy" ,       ##Run EpiCopy CNV calling (get segments)

##normal.path <- paste0(work.dir,
##					            file.sep,
##                      "BLCA_normal_subsample_10.txt")

##tumor.path <- paste0(work.dir,
##					           file.sep,
##                      "BLCA_tumor_subsample_10.txt")

##normal = read.table(file = normal.path,
##                    header = TRUE,
##                    sep = "\t",
##                    stringsAsFactors = FALSE
##                    )	
##rownames(normal) <- normal$X

##tumor  = read.table(file = tumor.path,
##                    header = TRUE,
##                    sep="\t",
##                    stringsAsFactors = FALSE
##                    )
##rownames(tumor) <- tumor$X				

##Sample Sheet
##sample.sheet.csv <- paste0(idat.pooled.files.dir,
##                           file.sep,
##                           "Sample_Sheet.csv") %>% 
##  read.csv(header = TRUE,
##           stringsAsFactors = FALSE)

##subset(sample.sheet.csv, Sample_Name %in% tumor$X)  %>% nrow() ##10
##subset(sample.sheet.csv, Sample_Name %in% normal$X) %>% nrow() ##10

##length(intersect(sample.sheet.csv$Sample_Name, tumor$X))  ##10
##length(intersect(sample.sheet.csv$Sample_Name, normal$X)) ##10

config.in <- scan(file=paste0("C:\\Users\\Mike",
                              "\\Desktop\\github",
                              "\\scripts\\DARTMOUTH",
                              "\\multi_methylv1.0",
                              "\\multi_methyl.config"),
                  comment.char="#",
                  what="character")

config.frame <- 
  data.frame(param=unlist(strsplit(config.in,split="="))[c(TRUE,FALSE)],
             value=unlist(strsplit(config.in,split="="))[c(FALSE,TRUE)])

