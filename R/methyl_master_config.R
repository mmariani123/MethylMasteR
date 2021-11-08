#!/usr/bin/env Rscript

##General
os.type="linux" ##Set to <"windows"> or <"mac"> or <"linux">
sample.sheet="Sample_Sheet.csv"
r.lib.path="/dartfs/rc/lab/S/SalasLab/R/x86_64-pc-linux-gnu-library/4.0"
##Need to set on Discovery or else will default to home
proj="TCGA-BLCA"
visualize=FALSE ##output plots
visualize.individual=FALSE ##output sesame style individual plots
weighted.mean="normal" ##can be 'nomral' or 'epicopy'
## Possible routines:
## "test"            ##Run a quick test
## "download",       ##For downloading TCGA data
## "process_sesame", ##Preprocess the TCGA and cord data in sesame format
## "sesame",         ##Run Sesame  CNV calling (get segments)
## "k450" ,          ##Run 450K    CNV calling (get segments)
## "champ" ,         ##Run ChAMP   CNV calling (get segments)
## "epicopy" ,       ##Run EpiCopy CNV calling (get segments)
routine="sesame"
sex="gender_reported"

####################### For download ####################################
download.barcodes.file="BLCA_tumor_subsample_sample_names.txt"
download.project=proj
download.data.category="Raw microarray data"
download.data.type="Raw intensities"
download.exp.strategy="Methylation array"
download.legacy=TRUE
download.file.type=".idat"
download.platform="Illumina Human Methylation 450"
download.method="api"
download.files.per.chunk=20
download.method="client"

######################## For sesame ######################################
sesame.platform="HM450" ##For TCGA
sesame.cord.platform="EPIC" ##For cord reference
sesame.ref="norm" ##Sesame ref can be
                  ##"norm" - built-in,
                  ##"sp"   - make sure sex variable is in sample sheet
                  ##"cord" - built in
sesame.data.cache="EPIC"
sesame.data.normal='EPIC.5.normal'
sesame.ref.version="hg38"

####################### For k450 ##########################################
k450.ref="norm" ##"cord"
k450.workflow="C" ##which k450 workflow? Can be "A" , "B-Ztransform", "C" conumee or "none"

####################### For champ #########################################
champ.array.type="450K"
champ.batch.name=c("Batch")
champ.padj=0.05
champ.ncores=16
champ.control=TRUE
champ.control.group="normal" ##champ.contol.group="champCtls"
champ.runimpute=TRUE
champ.runQC=TRUE
champ.runnorm=TRUE
champ.runSVD=TRUE
champ.runCombat=TRUE
champ.runDMP=TRUE
champ.runDMR=TRUE
champ.runBlock=TRUE
champ.runGSEA=TRUE

##for epicopy
epi.platform="hm450" ##"hm450" or "epic" mode depending on input data
use.epi.normal=FALSE ##Use the epicopy built in ref?
epi.ref="median" ##or can be 'mode' (only used if <epi.platform>==TRUE)
less.stringent.ra=FALSE ##FALSE: use popranges to get more confident ranges
epi.output.dir=FALSE
epi.normals="Sample_Group" ##In sample sheet
##sampNames="Sample_Name"
epi.samp.names=NULL
epi.qn=FALSE
epi.mode.bw=0.1
epi.mode.method="naive"
epi.normal.cnv=TRUE
epi.mean.center=TRUE
epi.filter.probes=FALSE
epi.retained.probes=NULL
epi.keepfnobj=TRUE
epi.fn.output=NULL

################ For champ extra functions ###############################
champ.beta          = myLoad$beta
champ.pheno         = myLoad$pd$Sample_G
champ.mds.plot      = TRUE
champ.density.plot  = TRUE
champ.dendrogram    = TRUE
champ.pdf.plot      = TRUE
champ.rplot         = TRUE
champ.feature.sel   = "None"
champ.results.dir   = "./CHAMP_QCimages/"
