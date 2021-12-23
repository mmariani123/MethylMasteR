#!/usr/bin/env Rscript

##Michael Mariani PhD Dartmouth College 2021

#' @title champ.process.mm
#' @description champ.process function modified by
#' Michael Mariani PhD Dartmouth College 2021
#' @param runload
#' @param directory
#' @param resultsDir
#' @param sample.sheet
#' @param arraytype
#' @param filters
#' @param runimpute
#' @param imputemethod
#' @param runQC
#' @param QCplots
#' @param runnorm
#' @param normalizationmethod
#' @param runSVD
#' @param RGEffect
#' @param runCombat
#' @param batchname
#' @param runDMP
#' @param runDMR
#' @param DMRmethod
#' @param runBlock
#' @param runGSEA
#' @param runEpiMod
#' @param runCNA
#' @param control
#' @param controlGroup
#' @param runRefBase
#' @param compare.group
#' @param adjPVal
#' @param resultsDir
#' @param PDFplot
#' @param Rplot
#' @param cores
#' @param saveStepresults
#' @import ChAMP
#' @import ChAMPdata
#' @importFrom dendextend labels_colors unbox
#' @return Return a champ.result object
#' @export
champ.process.mm <- function(runload = TRUE,
                             directory = getwd(),
                             resultsDir = "./CHAMP_RESULT/",
                             sample.sheet = NULL,
                             arraytype = "450K",
                             controlGroup = "champCtls",
                             filters = c("XY",
                                         "DetP",
                                         "Beads",
                                         "NoCG",
                                         "SNP",
                                         "MultiHit"),
                             runimpute = TRUE,
                             imputemethod = "Combine",
                             runQC = TRUE,
                             QCplots = c("mdsPlot",
                                         "densityPlot",
                                         "dendrogram"),
                             runnorm = TRUE,
                             normalizationmethod = "BMIQ",
                             runSVD              = TRUE,
                             RGEffect            = FALSE,
                             runCombat           = TRUE,
                             batchname           = c("Slide"),
                             runDMP              = TRUE,
                             runDMR              = TRUE,
                             DMRmethod           = "Bumphunter",
                             runBlock            = TRUE,
                             runGSEA             = TRUE,
                             runEpiMod           = TRUE,
                             runCNA              = TRUE,
                             control             = TRUE,
                             runRefBase          = FALSE,
                             compare.group       = NULL,
                             adjPVal             = 0.05,
                             PDFplot             = TRUE,
                             Rplot               = TRUE,
                             cores               = 1,
                             saveStepresults     = TRUE
                             ){
  message("[===========================]")
  message("[<<< ChAMP.PROCESS START >>>]")
  message("-----------------------------")
  message("This is champ.process() function, ",
          "which would run ALL analysis ChAMP ",
          "package can be done on your package. ",
          "But sometimes it's not always very ",
          "smoothly to conduct all function ",
          "because of mistake in dataset or ",
          "unsignificant results. ",
          "This if you encounter problem in ",
          "champ.proces(), you shall try run ",
          "champ step by step.^_^")
  if (!file.exists(resultsDir)) {
    dir.create(resultsDir)
    message("Creating results directory. All Results will be saved in ",
            resultsDir)
  }
  CHAMP.RESULT <- list()
  if (runload) {
    myfilter <- rep(FALSE, 6)
    names(myfilter) <- c("XY", "DetP", "Beads", "NoCG",
                         "SNP", "MultiHit")
    myfilter[filters] <- TRUE
    message("\nRunning champ.load()...")
    myLoad <- ChAMP::champ.load(directory = directory,
                         sample.sheet = sample.sheet,
                         compare.group = compare.group,
                         control = FALSE,
                         methValue = "B",
                         autoimpute = TRUE,
                         filterXY = myfilter["XY"],
                         filterDetP = myfilter["DetP"],
                         ProbeCutoff = 0,
                         SampleCutoff = 0.1,
                         detPcut = 0.01,
                         filterBeads = myfilter["Beads"],
                         beadCutoff = 0.05,
                         filterNoCG = myfilter["NoCG"],
                         filterSNPs = myfilter["SNP"],
                         filterMultiHit = myfilter["MultiHit"],
                         arraytype = arraytype)
    if (saveStepresults) {
      save(myLoad, file = paste(resultsDir, "/myLoad.rda",
                                sep = ""))
      message("champ.load()'s result \"myLoad\" has been saved in ",
              resultsDir, " as \"myLoad.rda.\"")
    }
    message("Run champ.load() Over!\n")
    gc()
    CHAMP.RESULT[["champ.load"]] <- myLoad
  }
  tmpbeta <- myLoad$beta
  tmppd <- myLoad$pd
  if (runimpute) {
    message("\nRunning champ.impute()...")
    myImpute <- ChAMP::champ.impute(beta = myLoad$beta,
                             pd = myLoad$pd,
                             method = imputemethod,
                             k = 5,
                             ProbeCutoff = 0.2,
                             SampleCutoff = 0.1)
    if (saveStepresults) {
      save(myImpute, file = paste(resultsDir, "/myImpute.rda",
                                  sep = ""))
      message("champ.impute()'s result \"myImpute\" has been saved in ",
              resultsDir, " as \"myImpute.rda.\"")
    }
    gc()
    CHAMP.RESULT[["champ.impute"]] <- myImpute
    message("Run champ.impute() Over!\n")
  }
  tmppd <- myImpute$pd
  if (runQC) {
    message("\nRunning champ.QC()...")
    myQCplots <- rep(FALSE, 3)
    names(myQCplots) <- c("mdsPlot", "densityPlot", "dendrogram")
    myQCplots[QCplots] <- TRUE
    ChAMP::champ.QC(beta = tmpbeta, pheno = tmppd$Sample_Group,
             mdsPlot = myQCplots["mdsPlot"],
             densityPlot = myQCplots["densityPlot"],
             dendrogram = myQCplots["dendrogram"],
             PDFplot = PDFplot,
             Rplot = Rplot, Feature.sel = "None",
             resultsDir = paste(resultsDir,
                                "/CHAMP_QCimages/", sep = ""))
    if (PDFplot) {
      message("Plots of champ.QC() has been saved in ",
              paste(resultsDir, "/CHAMP_QCimages/", sep = ""))
    }
    gc()
    message("Run champ.QC() Over!\n")
  }
  if (runnorm) {
    message("\nRunning champ.norm()...")
    myNorm <- ChAMP::champ.norm(beta = tmpbeta, rgSet = myLoad$rgSet,
                         mset = myLoad$mset,
                         resultsDir = paste(resultsDir,
                          "/CHAMP_Normalization/", sep = ""),
                         method = normalizationmethod,
                         plotBMIQ = PDFplot,
                         arraytype = arraytype,
                         cores = cores)
    if (saveStepresults) {
      save(myNorm, file = paste(resultsDir, "/myNorm.rda",
                                sep = ""))
      message("champ.norm()'s result \"myNorm\" has been saved in ",
              resultsDir, " as \"myNorm.rda.\"")
      if (normalizationmethod == "BMIQ" & PDFplot == TRUE)
        message("Plots of champ.norm() has been saved in ",
                paste(resultsDir, "/CHAMP_Normalization/",
                      sep = ""))
    }
    gc()
    CHAMP.RESULT[["champ.norm"]] <- myNorm
    message("Run champ.norm() Over!\n")
  }
  if (runSVD) {
    message("\nRunning champ.SVD()...")
    ChAMP::champ.SVD(beta = myNorm, rgSet = myLoad$rgSet, pd = tmppd,
              RGEffect = RGEffect, PDFplot = PDFplot, Rplot = Rplot,
              resultsDir = paste(resultsDir, "/CHAMP_SVDimages/",
                                 sep = ""))
    if (PDFplot) {
      message("Plots of champ.SVD() has been saved in ",
              paste(resultsDir, "/CHAMP_SVDimages/", sep = ""))
    }
    gc()
    message("Run champ.SVD() Over!\n")
  }
  if (runCombat) {
    message("\nRunning champ.runCombat()...")
    myCombat <- ChAMP::champ.runCombat(beta = myNorm, pd = tmppd,
                                batchname = batchname, logitTrans = TRUE)
    if (saveStepresults) {
      save(myCombat, file = paste(resultsDir, "/myCombat.rda",
                                  sep = ""))
      message("champ.runCombat()'s result \"myCombat\" has been saved in ",
              resultsDir, " as \"myCombat.rda.\"")
    }
    gc()
    CHAMP.RESULT[["champ.runCombat"]] <- myCombat
    message("Run champ.runCombat() Over!\n")
  }
  if (runDMP) {
    message("\nRunning champ.DMP()...")
    myDMP <- ChAMP::champ.DMP(beta = myNorm,
                       pheno = tmppd$Sample_Group,
                       adjPVal = adjPVal,
                       adjust.method = "BH",
                       compare.group = compare.group,
                       arraytype = arraytype)
    if (saveStepresults) {
      save(myDMP, file = paste(resultsDir, "/myDMP.rda",
                               sep = ""))
      message("champ.DMP()'s result \"myDMP\" has been saved in ",
              resultsDir, " as \"myDMP.rda.\"")
    }
    gc()
    CHAMP.RESULT[["champ.DMP"]] <- myDMP
    message("Run champ.DMP() Over!\n")
  }
  if (runDMR) {
    message("\nRunning champ.DMR()...")
    myDMR <- ChAMP::champ.DMR(beta = myNorm,
                       pheno = tmppd$Sample_Group,
                       arraytype = arraytype,
                       method = DMRmethod,
                       minProbes = 7,
                       adjPvalDmr = adjPVal,
                       maxGap = 300,
                       cutoff = NULL,
                       pickCutoff = TRUE,
                       smooth = TRUE,
                       smoothFunction = loessByCluster,
                       useWeights = FALSE,
                       permutations = NULL,
                       B = 250,
                       nullMethod = "bootstrap",
                       cores = cores,
                       meanLassoRadius = 375,
                       minDmrSep = 1000,
                       minDmrSize = 50,
                       adjPvalProbe = adjPVal,
                       Rplot = Rplot,
                       PDFplot = PDFplot,
                       resultsDir = paste(resultsDir,
                        "/CHAMP_ProbeLasso/", sep = ""),
                        rmSNPCH = T,
                       dist = 2, mafcut = 0.05)
    if (saveStepresults) {
      save(myDMR, file = paste(resultsDir, "/myDMR.rda",
                               sep = ""))
      message("champ.DMR()'s result \"myDMR\" has been saved in ",
              resultsDir, " as \"myDMR.rda.\"")
      if (DMRmethod == "ProbeLasso" & PDFplot == TRUE)
        message("Plots of champ.DMR() have been saved in ",
                paste(resultsDir, "/CHAMP_ProbeLasso/", sep = ""))
    }
    gc()
    CHAMP.RESULT[["champ.DMR"]] <- myDMR
    message("Run champ.DMR() Over!\n")
  }
  if (runBlock) {
    message("\nRunning champ.Block()...")
    myBlock <- ChAMP::champ.Block(beta = myNorm,
                           pheno = tmppd$Sample_Group,
                           arraytype = arraytype,
                           maxClusterGap = 250000,
                           B = 500,
                           bpSpan = 250000,
                           minNum = 10,
                           cores = cores)
    if (saveStepresults) {
      save(myBlock, file = paste(resultsDir, "/myBlock.rda",
                                 sep = ""))
      message("champ.Block()'s result \"myBlock\" has been saved in ",
              resultsDir, " as \"myBlock.rda.\"")
    }
    gc()
    CHAMP.RESULT[["champ.Block"]] <- myBlock
    message("Run champ.Block() Over!\n")
  }
  if (runGSEA) {
    message("\nRunning champ.GSEA()...")
    myGSEA <- ChAMP::champ.GSEA(beta = myNorm, DMP = myDMP, DMR = myDMR,
                         CpGlist = NULL, Genelist = NULL, arraytype = arraytype,
                         adjPval = adjPVal)
    if (saveStepresults) {
      save(myGSEA, file = paste(resultsDir, "/myGSEA.rda",
                                sep = ""))
      message("champ.GSEA()'s result \"myGSEA\" has been saved in ",
              resultsDir, " as \"myGSEA.rda.\"")
    }
    gc()
    CHAMP.RESULT[["champ.GSEA"]] <- myGSEA
    message("Run champ.GSEA() Over!\n")
  }
  if (runEpiMod) {
    message("\nRunning champ.EpiMod()...")
    myEpiMod <- champ.EpiMod(beta = myNorm,
                             pheno = tmppd$Sample_Group,
                             nseeds = 100,
                             gamma = 0.5,
                             nMC = 1000,
                             sizeR.v = c(1, 100),
                             minsizeOUT = 10,
                             resultsDir = paste(resultsDir,
                                                                                                                          "/CHAMP_EpiMod/", sep = ""), PDFplot = PDFplot,
                             arraytype = arraytype)
    if (saveStepresults) {
      save(myEpiMod, file = paste(resultsDir, "/myEpiMod.rda",
                                  sep = ""))
      message("champ.EpiMod()'s result \"myEpiMod\" has been saved in ",
              resultsDir, " as \"myEpiMod.rda.\"")
      if (PDFplot == TRUE)
        message("Plots of champ.EpiMod() have been saved in ",
                paste(resultsDir, "/CHAMP_EpiMod/", sep = ""))
    }
    gc()
    CHAMP.RESULT[["champ.EpiMod"]] <- myEpiMod
    message("Run champ.EpiMod() Over!\n")
  }
  if (runCNA) {
    message("\nRunning champ.CNA()...")
    myCNA <- ChAMP::champ.CNA(intensity = myLoad$intensity,
                       pheno = myLoad$pd$Sample_Group,
                       resultsDir = "./CHAMP_CNA",
                       control = TRUE,
                       controlGroup = "champCtls",
                       sampleCNA = TRUE,
                       groupFreqPlots = TRUE,
                       Rplot = Rplot,
                       PDFplot = PDFplot,
                       freqThreshold = 0.3,
                       arraytype = arraytype)
    if (saveStepresults) {
      save(myCNA, file = paste(resultsDir, "/myCNA.rda",
                               sep = ""))
      message("champ.CNA()'s result \"myCNA\" has been saved in ",
              resultsDir, " as \"myCNA.rda.\"")
    }
    gc()
    CHAMP.RESULT[["champ.CNA"]] <- myCNA
    message("Run champ.CNA() Over!\n")
  }
  if (runRefBase) {
    message("\nRunning champ.refbase()...")
    myRefbase <- ChAMP::champ.refbase(beta = myNorm, arraytype = arraytype)
    if (saveStepresults) {
      save(myRefbase, file = paste(resultsDir, "/myRefbase.rda",
                                   sep = ""))
      message("champ.refbase()'s result \"myRefbase\" has been saved in ",
              resultsDir, " as \"myRefbase.rda.\"")
    }
    gc()
    CHAMP.RESULT[["champ.refbase"]] <- myRefbase
    message("Run champ.refbase() Over!\n")
  }
  message("[<<<< ChAMP.PROCESS END >>>>]")
  message("[===========================]")
  return(CHAMP.RESULT)
}

#' @title champ.import.mm
#' @description champ.import function modified by Michael Mariani
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kept.
#' @param infile Path to the input file
#' @param directory
#' @param sample.sheet
#' @param offset
#' @param control
#' @param compare.group
#' @param arraytype
#' @importFrom illuminaio readIDAT unbox
#' @return Return a list
#' @export
champ.import.mm <- function(directory = getwd(),
                            sample.sheet = NULL,
                            offset = 100,
                            control = TRUE,
                            compare.group = compare.group,
                            arraytype = "450K"
                         ){

  message("[===========================]")
  message("[<<<< ChAMP.IMPORT START >>>>>]")
  message("-----------------------------")
  message("\n[ Section 1: Read PD Files Start ]")
  if (!file.exists(directory) || is.na(file.info(directory)$isdir) ||
      file.info(directory)$isdir == FALSE) {
    stop(paste0("  Your 'directory' for loading does not exists, ",
    "please assign a correct directory."))
  }
  ##MM
  if(!is.null(sample.sheet)){
    csvfile <- sample.sheet
  }
  else{
    csvfile <- list.files(directory,
                        recursive = TRUE,
                        pattern = "csv$",
                        full.names = TRUE)
  }
  if (length(csvfile) == 0) {
    stop(paste("  champ.import can not find any csv file in ",
               directory, ". Please check your folder."))
  }
  else if (length(csvfile) >= 2) {
    stop(paste("  champ.import finds more than one csv file in ",
               directory, ". Please check your folder."))
  }
  message("  CSV Directory: ", csvfile)
  message("  Find CSV Success")
  message("  Reading CSV File")
  skipline <- which(substr(readLines(csvfile), 1, 6) == "[Data]")
  if (length(skipline) == 0)
    suppressWarnings(pd <- read.csv(csvfile, stringsAsFactor = FALSE,
                                    header = TRUE))
  else suppressWarnings(pd <- read.csv(csvfile, skip = skipline,
                                       stringsAsFactor = FALSE, header = TRUE))

  #######################################################################

  if(control==TRUE){

    pd <- pd[pd$Sample_Group %in% compare.group,]

  }else{

    pd <- pd[pd$Sample_Group %in% compare.group[1],]

  }

  #######################################################################

  if ("Sentrix_Position" %in% colnames(pd)) {
    colnames(pd)[which(colnames(pd) == "Sentrix_Position")] <- "Array"
    message("  Replace Sentrix_Position into Array")
  }
  else {
    message("  Your pd file contains NO Array(Sentrix_Position) information.")
  }
  if ("Sentrix_ID" %in% colnames(pd)) {
    colnames(pd)[which(colnames(pd) == "Sentrix_ID")] <- "Slide"
    message("  Replace Sentrix_ID into Slide")
  }
  else {
    message("  Your pd file contains NO Slide(Sentrix_ID) information.")
  }
  sapply(c("Pool_ID", "Sample_Plate", "Sample_Well"),
         function(x) if (x %in% colnames(pd))
           pd[, x] <- as.character(pd[, x])
         else message("  There is NO ", x, " in your pd file."))

  GrnPath <- unlist(sapply(paste(pd$Slide, pd$Array, "Grn.idat",
                                 sep = "_"),
                           function(x) grep(x,
                              list.files(directory,
                              recursive = T, full.names = TRUE), value = TRUE)))

  RedPath <- unlist(sapply(paste(pd$Slide, pd$Array, "Red.idat",
                                 sep = "_"),
                           function(x) grep(x,
                                            list.files(directory,
                              recursive = T, full.names = TRUE), value = TRUE)))

  if (!identical(names(GrnPath), paste(pd$Slide, pd$Array,
                                       "Grn.idat", sep = "_")))
    stop("  Error Match between pd file and Green Channel IDAT file.")
  if (!identical(names(RedPath), paste(pd$Slide, pd$Array,
                                       "Red.idat", sep = "_")))
    stop("  Error Match between pd file and Red Channel IDAT file.")
  message("[ Section 1: Read PD file Done ]")
  message("\n\n[ Section 2: Read IDAT files Start ]")
  count <- 1
  G.idats <- lapply(GrnPath, function(x) {
    message("  Loading:", x, " ---- (", which(GrnPath ==
                                                x), "/", length(GrnPath), ")")
    illuminaio::readIDAT(x)
  })
  count <- 1
  R.idats <- lapply(RedPath, function(x) {
    message("  Loading:", x, " ---- (", which(RedPath ==
                                                x), "/", length(RedPath), ")")
    illuminaio::readIDAT(x)
  })
  names(G.idats) <- pd$Sample_Name
  names(R.idats) <- pd$Sample_Name
  checkunique <- unique(c(sapply(G.idats, function(x) nrow(x$Quants)),
                          sapply(R.idats, function(x) nrow(x$Quants))))
  if (length(checkunique) > 1) {

    message("\n  !!! Important !!! ")

    message(paste0("  Seems your IDAT files not from one Array, ",
                   "because they have different numbers of probe."))

    message(paste0("  ChAMP wil continue analysis with only COMMON ",
                   "CpGs exist across all your IDAt files. ",
                   "However we still suggest you to check your ",
                   "source of data.\n"))
  }
  CombineIDAT <- append(G.idats, R.idats)
  commonAddresses <- as.character(Reduce("intersect",
                                         lapply(CombineIDAT,
                                          function(x) rownames(x$Quants))))
  message("\n  Extract Mean value for Green and Red Channel Success")
  GreenMean <- do.call(cbind, lapply(G.idats,
                                     function(xx) xx$Quants[commonAddresses,
                                                                     "Mean"]))
  RedMean <- do.call(cbind, lapply(R.idats,
                                   function(xx) xx$Quants[commonAddresses,
                                                                   "Mean"]))
  message("    Your Red Green Channel contains ", nrow(GreenMean),
          " probes.")
  G.Load <- do.call(cbind, lapply(G.idats,
                                  function(x) x$Quants[commonAddresses,
                                                                "Mean"]))
  R.Load <- do.call(cbind, lapply(R.idats,
                                  function(x) x$Quants[commonAddresses,
                                                                "Mean"]))
  message("[ Section 2: Read IDAT Files Done ]")
  message("\n\n[ Section 3: Use Annotation Start ]")
  message("\n  Reading ", arraytype, " Annotation >>")
  ##########################################################################
  if(arraytype == "EPIC"){
    data(AnnoEPIC)
  }else{
    data(Anno450K)
  }
  ##########################################################################
  message("\n  Fetching NEGATIVE ControlProbe.")
  control_probe <- rownames(Anno$ControlProbe)[
    which(Anno$ControlProbe[, 1] == "NEGATIVE")]

  message("    Totally, there are ", length(control_probe),
          " control probes in Annotation.")

  control_probe <- control_probe[control_probe %in% rownames(R.Load)]

  message("    Your data set contains ", length(control_probe),
          " control probes.")

  rMu <- matrixStats::colMedians(R.Load[control_probe, ])
  rSd <- matrixStats::colMads(R.Load[control_probe, ])
  gMu <- matrixStats::colMedians(G.Load[control_probe, ])
  gSd <- matrixStats::colMads(G.Load[control_probe, ])
  rownames(G.Load) <- paste("G", rownames(G.Load), sep = "-")
  rownames(R.Load) <- paste("R", rownames(R.Load), sep = "-")
  IDAT <- rbind(G.Load, R.Load)
  message("\n  Generating Meth and UnMeth Matrix")
  message("    Extracting Meth Matrix...")
  M.check <- Anno$Annotation[, "M.index"] %in% rownames(IDAT)
  message("      Totally there are ", nrow(Anno$Annotation),
          " Meth probes in ", arraytype, " Annotation.")
  message("      Your data set contains ", length(M.check),
          " Meth probes.")
  M <- IDAT[Anno$Annotation[, "M.index"][M.check], ]
  message("    Extracting UnMeth Matrix...")
  U.check <- Anno$Annotation[, "U.index"] %in% rownames(IDAT)
  message("      Totally there are ", nrow(Anno$Annotation),
          " UnMeth probes in ", arraytype, " Annotation.")
  message("      Your data set contains ", length(U.check),
          " UnMeth probes.")
  U <- IDAT[Anno$Annotation[, "U.index"][U.check], ]
  if (!identical(M.check, U.check)) {
    stop("  Meth Matrix and UnMeth Matrix seems not paried correctly.")
  }
  else {
    CpG.index <- Anno$Annotation[, "CpG"][M.check]
  }
  rownames(M) <- CpG.index
  rownames(U) <- CpG.index
  message("\n  Generating beta Matrix")
  BetaValue <- M/(M + U + offset)
  message("  Generating M Matrix")
  MValue <- log2(M/U)
  message("  Generating intensity Matrix")
  intensity <- M + U
  message("  Calculating Detect P value")
  detP <- matrix(NA, nrow = nrow(intensity), ncol = ncol(intensity))
  rownames(detP) <- rownames(intensity)
  colnames(detP) <- colnames(intensity)
  type_II <- rownames(Anno$Annotation)[Anno$Annotation[, "Channel"] ==
                                         "g+r"]
  type_II <- type_II[type_II %in% rownames(detP)]
  type_I.red <- rownames(Anno$Annotation)[Anno$Annotation[,
                                                          "Channel"] == "r"]
  type_I.red <- type_I.red[type_I.red %in% rownames(detP)]
  type_I.grn <- rownames(Anno$Annotation)[Anno$Annotation[,
                                                          "Channel"] == "g"]
  type_I.grn <- type_I.grn[type_I.grn %in% rownames(detP)]
  for (i in 1:ncol(detP)) {
    detP[type_II, i] <- 1 - pnorm(intensity[type_II, i],
                                  mean = rMu[i] + gMu[i], sd = rSd[i] + gSd[i])
    detP[type_I.red, i] <- 1 - pnorm(intensity[type_I.red,
                          i], mean = rMu[i] * 2, sd = rSd[i] * 2)
    detP[type_I.grn, i] <- 1 - pnorm(intensity[type_I.grn,
                          i], mean = gMu[i] * 2, sd = gSd[i] * 2)
  }
  if (sum(is.na(detP)))
    message("    !!! There are NA values in your detP matrix.\n")
  message("  Counting Beads")
  NBeads <- do.call(cbind, lapply(R.idats, function(x) x$Quants[commonAddresses,
                                                                "NBeads"]))
  Mbead <- NBeads[substr(Anno$Annotation$M.index[M.check],
                         3, 100), ]
  Ubead <- NBeads[substr(Anno$Annotation$U.index[U.check],
                         3, 100), ]
  Ubead[Ubead < 3 | Mbead < 3] <- NA
  rownames(Ubead) <- rownames(intensity)
  message("[ Section 3: Use Annotation Done ]")
  message("\n[<<<<< ChAMP.IMPORT END >>>>>>]")
  message("[===========================]")
  message("[You may want to process champ.filter() next.]\n")
  return(list(beta = BetaValue, M = MValue, pd = pd, intensity = intensity,
              detP = detP, beadcount = Ubead, Meth = M, UnMeth = U))
}

#' @title champ.load.mm
#' @description champ.load function modified by Michael Mariani
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kept.
#' @param directory
#' @param sample.sheet
#' @param compare.group
#' @param control
#' @param method
#' @param methValue
#' @param autoimpute
#' @param filterDetP
#' @param ProbeCutoff
#' @param SampleCutoff
#' @param detPcut
#' @param filterBeads
#' @param beadCutoff
#' @param filterNoCG
#' @param filterSNPs
#' @param population
#' @param filterMultiHit
#' @param filterXY
#' @param force
#' @param arraytype
#' @importFrom ChAMP champ.filter unbox
#' @return Return a list of loaded data
#' @export
champ.load.mm <- function(directory   = getwd(),
                       sample.sheet   = NULL,
                       compare.group  = compare.group,
                       control        = FALSE,
                       method         = "ChAMP",
                       methValue      = "B",
                       autoimpute     = TRUE,
                       filterDetP     = TRUE,
                       ProbeCutoff    = 0,
                       SampleCutoff   = 0.1,
                       detPcut        = 0.01,
                       filterBeads    = TRUE,
                       beadCutoff     = 0.05,
                       filterNoCG     = TRUE,
                       filterSNPs     = TRUE,
                       population     = NULL,
                       filterMultiHit = TRUE,
                       filterXY       = TRUE,
                       force          = FALSE,
                       arraytype      = "450K")
{
  message("[===========================]")
  message("[<<<< ChAMP.LOAD START >>>>>]")
  message("-----------------------------")
  mybeadcount <- function(x) {
    nb <- getNBeads(x)
    typeIadd <- getProbeInfo(x, type = "I")
    typeImatchA <- match(typeIadd$AddressA, rownames(nb))
    typeImatchB <- match(typeIadd$AddressB, rownames(nb))
    typeIIadd <- getProbeInfo(x, type = "II")
    typeIImatch <- match(typeIIadd$Address, rownames(nb))
    nbcg <- nb
    locusNames <- getManifestInfo(x, "locusNames")
    bc_temp <- matrix(NA_real_, ncol = ncol(x), nrow = length(locusNames),
                      dimnames = list(locusNames, sampleNames(x)))
    TypeII.Name <- getProbeInfo(x, type = "II")$Name
    bc_temp[TypeII.Name, ] <- nbcg[getProbeInfo(x, type = "II")$AddressA,
    ]
    TypeI <- getProbeInfo(x, type = "I")
    bcB <- bc_temp
    bcA <- bc_temp
    bcB[TypeI$Name, ] <- nbcg[TypeI$AddressB, ]
    bcA[TypeI$Name, ] <- nbcg[TypeI$AddressA, ]
    bcB3 <- which(bcB < 3)
    bcA3 <- which(bcA < 3)
    bcA2 <- bcA
    bcB2 <- bcB
    bcA2[bcA3] <- NA
    bcA2[bcB3] <- NA
    bc <- data.frame(bcA2)
    bc
  }
  if (method == "minfi") {
    message("\n[ Loading Data with Minfi Method ]")
    message("----------------------------------")
    message("Loading data from ", directory)
    myDir <- directory
    suppressWarnings(targets <- read.metharray.sheet(myDir))
    rgSet <- read.metharray.exp(targets = targets,
                                extended = TRUE,
                                force = force)
    if (arraytype == "EPIC")
      rgSet@annotation <- c(array = "IlluminaHumanMethylationEPIC",
                            annotation = "ilm10b4.hg19")
    sampleNames(rgSet) = rgSet[[1]]
    pd <- pData(rgSet)
    mset <- preprocessRaw(rgSet)
    detP <- detectionP(rgSet)
    message("<< Read DataSet Success. >>\n")
    if (methValue == "B")
      tmp = getBeta(mset, "Illumina")
    else tmp = getM(mset)
    tmp[detP >= detPcut] <- NA
    message(paste0("The fraction of failed positions per sample\n",
                   "\n(You may need to delete samples with high proportion",
                   "of failed probes\n): "))
    numfail <- matrix(colMeans(is.na(tmp)))
    rownames(numfail) <- colnames(detP)
    colnames(numfail) <- "Failed CpG Fraction."
    print(numfail)
    RemainSample <- which(numfail < SampleCutoff)
    if (any(numfail >= SampleCutoff))
      message("The detSamplecut parameter is : ",
              SampleCutoff, "\nSamples : ",
              paste(rownames(numfail)[which(numfail >=
                SampleCutoff)], collapse = ","), " will be deleted.\n",
              "There are ", length(RemainSample),
              " samples left for analysis.\n")
    rgSet <- rgSet[, RemainSample]
    detP <- detP[, RemainSample]
    mset <- mset[, RemainSample]
    pd <- pd[RemainSample, ]
    tmp <- tmp[, RemainSample]
    if (filterDetP) {
      mset.f = mset[rowSums(is.na(tmp)) <= ProbeCutoff *
                      ncol(detP), ]
      if (ProbeCutoff == 0) {
        message("Filtering probes with a detection p-value above ",
                detPcut, " in one or more samples has removed ",
                dim(mset)[1] - dim(mset.f)[1],
                paste0(" probes from the analysis. ",
                "If a large number of probes have been removed,",
                "ChAMP suggests you to identify potentially bad samples."))
      }
      else {
        message("Filtering probes with a detection p-value above ",
                detPcut, " in at least ", ProbeCutoff *
                  100, "% of samples has removed ", dim(mset)[1] -
                  dim(mset.f)[1], paste0(" probes from the analysis. ",
                  "If a large number of probes have been removed, ",
                  "ChAMP suggests you look at the failedSample file ",
                  "to identify potentially bad samples."))
      }
      mset = mset.f
      tmp <- tmp[rowSums(is.na(tmp)) <= ProbeCutoff * ncol(detP),
      ]
      message("<< Filter DetP Done. >>\n")
    }
    if (sum(is.na(tmp)) == 0) {
      message(paste0("\nThere is no NA values in your matrix, ",
                     "there is no need to imputation.\n"))
    }
    else {
      message("\nThere are ", sum(is.na(tmp)), paste0(" NA remain in filtered ",
      "Data Set. Impute can be done for remain NAs, but not suitable ",
      "for small number samples. For small Data Set (like only 20 samples), ",
      "we suggest you set parameter ProbeCutoff as 0 in champ.load() ",
      "here, which would remove all NA involved probe no matter how ",
      "many samples of those probes are NA.\n"))
    }
    if (autoimpute & sum(is.na(tmp)) > 0) {
      message("Impute will be conducted here for remain ",
              sum(is.na(tmp)), paste0("  NAs. Note that if you don't do this, ",
              "NA values would be kept in your data set. ",
              "You may use champ.impute() function to do more ",
              "complex imputation as well."))
      message("\nImpute function is working now, it may need couple minutes...")
      zz <- file("ImputeMessage.Rout", open = "wt")
      sink(zz)
      sink(zz, type = "message")
      tmp <- impute.knn(tmp, k = 5)$data
      sink(type = "message")
      sink()
      message("<< Imputation Done. >>\n")
    }
    if (filterBeads) {
      bc = mybeadcount(rgSet)
      bc2 = bc[rowSums(is.na(bc)) < beadCutoff * (ncol(bc)),
      ]
      mset.f2 = mset[featureNames(mset) %in% row.names(bc2),
      ]
      tmp <- tmp[rownames(tmp) %in% row.names(bc2), ]
      message("Filtering probes with a beadcount <3 in at least ",
              beadCutoff * 100, "% of samples, has removed ",
              dim(mset)[1] - dim(mset.f2)[1], " from the analysis.")
      mset = mset.f2
      message("<< Filter Beads Done. >>\n")
    }
    if (filterNoCG) {
      mset.f2 = dropMethylationLoci(mset, dropCH = T)
      tmp <- tmp[rownames(tmp) %in% featureNames(mset.f2),
      ]
      message("Filtering non-cg probes, has removed ",
              dim(mset)[1] - dim(mset.f2)[1], " from the analysis.")
      mset <- mset.f2
      message("<< Filter NoCG Done. >>\n")
    }
    if (filterSNPs) {
      if (arraytype == "450K") {
        if (is.null(population)) {
          message("Using general 450K SNP list for filtering.")
          data(hm450.manifest.hg19)
          maskname <- rownames(hm450.manifest.hg19)[
            which(hm450.manifest.hg19$MASK_general == TRUE)]
        }
        else if (!population %in% c("AFR", "EAS",
                                    "EUR", "SAS", "AMR", "GWD",
                                    "YRI", "TSI", "IBS", "CHS",
                                    "PUR", "JPT", "GIH", "CHB",
                                    "STU", "ITU", "LWK", "KHV",
                                    "FIN", "ESN", "CEU", "PJL",
                                    "ACB", "CLM", "CDX", "GBR",
                                    "BEB", "PEL", "MSL", "MXL",
                                    "ASW")) {
          message("Using general 450K SNP list for filtering.")
          data(hm450.manifest.hg19)
          maskname <- rownames(hm450.manifest.hg19)[
            which(hm450.manifest.hg19$MASK_general == TRUE)]
        }
        else {
          message("Using ",
                  population,
                  " specific 450K SNP list for filtering.")
          data(hm450.manifest.pop.hg19)
          maskname <- rownames(hm450.manifest.pop.hg19)[
            which(hm450.manifest.pop.hg19[,
                                    paste("MASK_general_", population,
                                    sep = "")] == TRUE)]
        }
      }
      else {
        if (is.null(population)) {
          message("Using general EPIC SNP list for filtering.")
          data(EPIC.manifest.hg19)
          maskname <- rownames(EPIC.manifest.hg19)[
            which(EPIC.manifest.hg19$MASK_general == TRUE)]
        }
        else if (!population %in% c("AFR", "EAS",
                                    "EUR", "SAS", "AMR", "GWD",
                                    "YRI", "TSI", "IBS", "CHS",
                                    "PUR", "JPT", "GIH", "CHB",
                                    "STU", "ITU", "LWK", "KHV",
                                    "FIN", "ESN", "CEU", "PJL",
                                    "ACB", "CLM", "CDX", "GBR",
                                    "BEB", "PEL", "MSL", "MXL",
                                    "ASW")) {
          message("Using general EPIC SNP list for filtering.")
          data(EPIC.manifest.hg19)
          maskname <- rownames(EPIC.manifest.hg19)[
            which(EPIC.manifest.hg19$MASK_general == TRUE)]
        }
        else {
          message("Using ",
                  population,
                  " specific EPIC SNP list for filtering.")
          data(EPIC.manifest.pop.hg19)
          maskname <- rownames(EPIC.manifest.pop.hg19)[
            which(EPIC.manifest.pop.hg19[,
              paste("MASK_general_", population, sep = "")] == TRUE)]
        }
      }
      mset.f2 = mset[!featureNames(mset) %in% maskname,
      ]
      tmp <- tmp[!rownames(tmp) %in% maskname, ]
      message(paste0("Filtering probes with SNPs as identified in ",
      "Zhou's Nucleic Acids Research Paper, 2016, has removed "),
              dim(mset)[1] - dim(mset.f2)[1], " from the analysis.")
      mset = mset.f2
      message("<< Filter SNP Done. >>\n")
    }
    if (filterMultiHit) {
      data(multi.hit)
      mset.f2 = mset[!featureNames(mset) %in% multi.hit$TargetID,
      ]
      tmp <- tmp[!rownames(tmp) %in% multi.hit$TargetID,
      ]
      message(paste0("Filtering probes that align to multiple ",
      "locations as identified in Nordlund et al, has removed "),
              dim(mset)[1] - dim(mset.f2)[1], " from the analysis.")
      mset = mset.f2
      message("<< Filter MultiHit Done. >>\n")
    }
    if (filterXY) {
      if (arraytype == "EPIC")
        data(probe.features.epic)
      else data(probe.features)
      autosomes = probe.features[!probe.features$CHR %in%
                                   c("X", "Y"), ]
      mset.f2 = mset[featureNames(mset) %in% row.names(autosomes),
      ]
      tmp <- tmp[rownames(tmp) %in% row.names(autosomes),
      ]
      message("Filtering probes on the X or Y chromosome has removed ",
              dim(mset)[1] - dim(mset.f2)[1], " from the analysis.")
      mset = mset.f2
      message("<< Filter XY chromosome Done. >>\n")
    }
    message(paste(if (methValue == "B")
      "[Beta"
      else "[M", "value is selected as output.]\n"))
    beta.raw <- tmp
    intensity <- minfi::getMeth(mset) + minfi::getUnmeth(mset)
    detP <- detP[which(row.names(detP) %in% row.names(beta.raw)),
    ]
    if (min(beta.raw, na.rm = TRUE) <= 0)
      beta.raw[beta.raw <= 0] <- min(beta.raw[beta.raw >
                                                0])
    message(paste0("Zeros in your dataset have been replaced ",
      "with smallest positive value.\n"))
    if (max(beta.raw, na.rm = TRUE) >= 0)
      beta.raw[beta.raw >= 1] <- max(beta.raw[beta.raw <
                                                1])
    message(paste0("One in your dataset have been replaced ",
      "with largest value below 1.\n"))
    message("The analysis will be proceed with ", dim(beta.raw)[1],
            " probes and ", dim(beta.raw)[2], " samples.\n")
    message("Current Data Set contains ", sum(is.na(beta.raw)),
            " NA in ", if (methValue == "B")
              "[Beta]"
            else "[M]", " Matrix.\n")
    message("[<<<<< ChAMP.LOAD END >>>>>>]")
    message("[===========================]")
    message("[You may want to process champ.QC() next.]\n")
    return(list(mset = mset,
                rgSet = rgSet,
                pd = pd,
                intensity = intensity,
                beta = beta.raw,
                detP = detP))
  }
  else {
    message("\n[ Loading Data with ChAMP Method ]")
    message("----------------------------------")
    message(paste0("Note that ChAMP method will NOT return rgSet or mset, ",
    "they object defined by minfi. Which means, ",
    "if you use ChAMP method to load data, ",
    "you can not use SWAN or FunctionNormliazation method ",
    "in champ.norm() (you can use BMIQ or PBC still). ",
    "But All other function should not be influenced.\n"))

    myImport <- ChAMP::champ.import(directory,
                             arraytype = arraytype,
                             compare.group = compare.group,
                             control = FALSE,
                             sample.sheet=sample.sheet
                             )

    if (methValue == "B")
      myLoad <- ChAMP::champ.filter(beta = myImport$beta,
                             M = NULL,
                      pd = myImport$pd,
                      intensity = myImport$intensity,
                      Meth = NULL,
                      UnMeth = NULL,
                      detP = myImport$detP,
                      beadcount = myImport$beadcount,
                      autoimpute = autoimpute,
                      filterDetP = filterDetP,
                      ProbeCutoff = ProbeCutoff,
                      SampleCutoff = SampleCutoff,
                      detPcut = detPcut,
                      filterBeads = filterBeads,
                      beadCutoff = beadCutoff,
                      filterNoCG = filterNoCG,
                      filterSNPs = filterSNPs,
                      population = population,
                      filterMultiHit = filterMultiHit,
                      filterXY = filterXY,
                      arraytype = arraytype)

    else myLoad <- ChAMP::champ.filter(beta = NULL,
                                M = myImport$M,
                      pd = myImport$pd,
                      intensity = myImport$intensity,
                      Meth = NULL,
                      UnMeth = NULL,
                      detP = myImport$detP,
                      beadcount = myImport$beadcount,
                      autoimpute = autoimpute,
                      filterDetP = filterDetP,
                      ProbeCutoff = ProbeCutoff,
                      SampleCutoff = SampleCutoff,
                      detPcut = detPcut,
                      filterBeads = filterBeads,
                      beadCutoff = beadCutoff,
                      filterNoCG = filterNoCG,
                      filterSNPs = filterSNPs,
                      population = population,
                      filterMultiHit = filterMultiHit,
                      filterXY = filterXY,
                      arraytype = arraytype)

    message("[<<<<< ChAMP.LOAD END >>>>>>]")
    message("[===========================]")
    message("[You may want to process champ.QC() next.]\n")
    return(myLoad)

  }
}

#' # Extra ChAMP functions
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param infile Path to the input file
#' @param champ.beta
#' @param champ.pheno
#' @param champ.mds.plot
#' @param champ.density.plot
#' @param champ.dendrogram
#' @param champ.pdf.plot
#' @param champ.rplot
#' @param champ.feature.sel
#' @param champ.results.dir
#' @import dendextend
#' @return A matrix of the infile
#' @export
champ.QC <- function(beta        = NULL,
                     pheno       = NULL,
                     mdsPlot     = NULL,
                     densityPlot = NULL,
                     dendrogram  = NULL,
                     PDFplot     = NULL,
                     Rplot       = NULL,
                     Feature.sel = NULL,
                     resultsDir  = NULL
                     ){
    par(mar=c(1,1,1,1))
    message("[===========================]")
    message("[<<<<< ChAMP.QC START >>>>>>]")
    message("-----------------------------")

    if(!file.exists(resultsDir)){dir.create(resultsDir)}

    message("champ.QC Results will be saved in ", resultsDir)

    message("[QC plots will be proceed with ", dim(beta)[1],
            " probes and ", dim(beta)[2], " samples.]\n")

    if(min(beta, na.rm = TRUE) == 0){
        beta[beta == 0] <- 1e-06
        message("[", length(which(beta == 0)),
                paste0(" Zeros dectect in your ",
                "dataset, will be replaced with 0.000001]\n"))
    }

    if(ncol(beta) != length(pheno)){
        stop(paste0("Dimension of DataSet Samples, ",
                    "pheno and name must be the same. ",
                    "Please check your input."))
        message("<< Prepare Data Over. >>")
    }

    if(mdsPlot){
        if(Rplot){
            mdsPlot(beta,
                    numPositions = 1000,
                    sampGroups = pheno,
                    colnames(beta)
            )
            if(PDFplot){
                pdf(paste(resultsDir,
                          "raw_mdsPlot.pdf",
                          sep = "/"),
                    width = 6,
                    height = 4)
                mdsPlot(beta,
                        numPositions = 1000,
                        sampGroups = pheno,
                        colnames(beta)
                )
                dev.off()
            }
            message("<< plot mdsPlot Done. >>\n")
        }
    }
    if(densityPlot){
        if(Rplot)
            densityPlot(beta,
                        sampGroups = pheno,
                        main = paste("Density plot of raw data (",
                                     nrow(beta),
                                     " probes)",
                                     sep = ""),
                        xlab = "Beta")
        if(PDFplot){
            pdf(paste(resultsDir,
                      "raw_densityPlot.pdf",
                      sep = "/"),
                width = 6,
                height = 4)
            densityPlot(beta,
                        sampGroups = pheno,
                        main = paste("Density plot of raw data (",
                                     nrow(beta), " probes)", sep = ""),
                        xlab = "Beta")
            dev.off()
        }
        message("<< Plot densityPlot Done. >>\n")
    }
    if(dendrogram){
        if(Feature.sel == "None"){
            message(paste0("< Dendrogram Plot Feature Selection Method >: ",
                           "No Selection, directly use all CpGs to calculate",
                           " distance matrix."))
            hc <- hclust(dist(t(beta)))
        }
        else if (Feature.sel == "SVD") {
            message(paste0("< Dendrogram Feature Selection Method >: ",
                           "Use top SVD CpGs to calculate distance matrix."))
            SVD <- svd(beta)
            rmt.o <- EstDimRMT(beta - rowMeans(beta))
            M <- SVD$v[, 1:rmt.o$dim]
            rownames(M) <- colnames(beta)
            colnames(M) <- paste("Component", c(1:rmt.o$dim))
            hc <- hclust(dist(M))
        }
        dend <- as.dendrogram(hc)
        MyColor <- rainbow(length(table(pheno)))
        names(MyColor) <- names(table(pheno))
        dendextend::labels_colors(dend) <-
          MyColor[pheno[order.dendrogram(dend)]]
        dend <- dend %>% dendextend::set("labels_cex", 0.8)
        ##MM change below
        dend <- dend %>%
          dendextend::set("leaves_pch", 19) %>%
          dendextend::set("leaves_cex", 0.6) %>%
          dendextend::set("leaves_col",
                MyColor[pheno[order.dendrogram(dend)]][1])
        if(Rplot){
            plot(dend,
                 center = TRUE,
                 main = paste("All samples before normalization (",
                              nrow(beta), " probes)", sep = ""))
            legend("topright", fill = MyColor, legend = names(MyColor))
        }
        if(PDFplot){
            pdf(paste(resultsDir,
                      "raw_SampleCluster.pdf",
                      sep = "/"),
                width = floor(log(ncol(beta)) * 3),
                height = floor(log(ncol(beta)) * 2)
            )
            plot(dend,
                 center = TRUE,
                 main = paste("All samples before normalization (",
                              nrow(beta), " probes)", sep = "")
            )
            legend("topright", fill = MyColor, legend = names(MyColor))
            dev.off()
        }
        message("<< Plot dendrogram Done. >>\n")
    }
    message("[<<<<<< ChAMP.QC END >>>>>>>]")
    message("[===========================]")
    message("[You may want to process champ.norm() next.]\n")
}

#' GenStatM.R
#'
#' https://rdrr.io/bioc/FEM/src/R/GenStatM.R
#'
#'
#' @param dnaM.m
#' @param pheno.v
#' @param chiptype
#' @importFrom limma lmFit unbox
#' @importFrom limma contrasts.fit unbox
#' @importFrom limma eBayes unbox
#' @importFrom limma  topTable unbox
#' @return A matrix of the infile
#' @export
GenStatM <- function(dnaM.m,
                     pheno.v,
                     chiptype="450k"
                     ){

    if (chiptype == "450k"){
        data("probe450kfemanno")
        probefemanno <- probe450kfemanno
    }
    else if (chiptype == "EPIC" ){
        data("probeEPICfemanno")
        probefemanno <- probeEPICfemanno
    }
    else{
        print("ERROR: Please indicate correct data type!")
        break
    }

    extractFn <- function(tmp.v, ext.idx) {
        return(tmp.v[ext.idx])
    }
    map.idx <- match(rownames(dnaM.m), probefemanno$probeID);
    probeInfo.lv <- lapply(probefemanno, extractFn, map.idx)
    beta.lm <- list()
    for (g in 1:6) {
        group.idx <- which(probeInfo.lv[[3]] == g)
        tmp.m <- dnaM.m[group.idx, ]
        rownames(tmp.m) <- probeInfo.lv$eid[group.idx];
        sel.idx <- which(is.na(rownames(tmp.m)) == FALSE);
        tmp.m <- tmp.m[sel.idx,];
        nL <- length(factor(rownames(tmp.m)));
        nspg.v <- summary(factor(rownames(tmp.m)),maxsum=nL);
        beta.lm[[g]] <- rowsum(tmp.m,group=rownames(tmp.m))/nspg.v;
        print(paste("Done for regional gene group ", g, sep = ""))
    }
    unqEID.v <- unique(c(rownames(beta.lm[[2]]), rownames(beta.lm[[4]]),
                         rownames(beta.lm[[1]])))
    avbeta.m <- matrix(nrow = length(unqEID.v), ncol = ncol(dnaM.m))
    colnames(avbeta.m) <- colnames(dnaM.m)
    rownames(avbeta.m) <- unqEID.v
    for (gr in c(1, 4, 2)) {
        avbeta.m[match(rownames(beta.lm[[gr]]), rownames(avbeta.m)),
        ] <- beta.lm[[gr]]
    }
    data.m <- avbeta.m
    sampletype.f <- as.factor(pheno.v)
    design.sample <- model.matrix(~0 + sampletype.f)
    colnames(design.sample) <- levels(sampletype.f)
    sampletypes.v <- levels(sampletype.f)
    lmf.o <- limma::lmFit(data.m, design.sample)
    ntypes <- length(levels(sampletype.f))
    ncomp <- 0.5 * ntypes * (ntypes - 1)
    cont.m <- matrix(0, nrow = ncol(design.sample), ncol = ncomp)
    tmp.v <- vector()
    c <- 1
    for (i1 in 1:(ntypes - 1)) {
        for (i2 in (i1 + 1):ntypes) {
            cont.m[i1, c] <- -1
            cont.m[i2, c] <- 1
            tmp.v[c] <- paste(sampletypes.v[i2], "--", sampletypes.v[i1],
                              sep = "")
            c <- c + 1
        }
    }
    rownames(cont.m) <- sampletypes.v
    colnames(cont.m) <- tmp.v
    lmf2.o <- limma::contrasts.fit(lmf.o, cont.m)
    bay.o  <- limma::eBayes(lmf2.o)
    top.lm <- list()
    for (c in 1:ncol(cont.m)) {
        top.lm[[c]] <- limma::topTable(bay.o,
                                coef = c,
                                adjust.method = "fdr",
                                number = nrow(data.m))
    }

    return(list(top = top.lm, cont = cont.m, avbeta = avbeta.m))

}

#' DoIntEpi450k
#'
#' ##https://rdrr.io/bioc/FEM/src/R/DoIntEpi450k.R
#'
#'
#' @param statM.o
#' @param adj.m
#' @param c
#' @param dmaMode
#' @importFrom igraph graph.adjacency unbox
#' @importFrom igraph get.edgelist unbox
#' @importFrom igraph clusters unbox
#' @return A matrix of the infile
#' @export
DoIntEpi450k <-
    function(statM.o,
             adj.m,
             c,
             dmaMode="avbeta"
             ){

        if(length(grep("[a-zA-Z]",rownames(adj.m)))!=0){
            print("ERROR: The rownames of adj.m should be EntrezID");
            break
        }

        if (dmaMode == "avbeta"){
            avbeta.m <- statM.o$avbeta;
            commonEID.v <- intersect(rownames(adj.m),rownames(avbeta.m));
            mapA.idx <- match(commonEID.v,rownames(adj.m));
            tmpA.m <- adj.m[mapA.idx,mapA.idx];

            mapM.idx <- match(commonEID.v,rownames(avbeta.m));
            tmpM.m <- avbeta.m[mapM.idx,];

            gr.o <- igraph::graph.adjacency(tmpA.m,mode="undirected");
            comp.l <- igraph::clusters(gr.o);
            ngpc.v <- summary(factor(comp.l$member));
            maxCLid <- as.numeric(names(ngpc.v)[which.max(ngpc.v)]);
            maxc.idx <- which(comp.l$member==maxCLid);

            ##get the max connected network
            tmpA.m <- tmpA.m[maxc.idx,maxc.idx];
            gr.o <-  igraph::graph.adjacency(tmpA.m,mode="undirected");
            tmpE.m <- igraph::get.edgelist(gr.o);
            tmpM.m <- tmpM.m[maxc.idx,];

            #### now extract statistics
            selcol.idx <- match(c("t","P.Value"),colnames(statM.o$top[[c]]));
            ordrow.idx <- match(rownames(tmpM.m),rownames(statM.o$top[[c]]));
            if(length(intersect(c("ID"),colnames(statM.o$top[[c]])))==1){
                ordrow.idx <- match(rownames(tmpM.m),statM.o$top[[c]][,1]);
            }
            statM.m <- statM.o$top[[c]][ordrow.idx,selcol.idx];
            rownames(statM.m) <- rownames(tmpM.m);

            return(list(statM=statM.m,adj=tmpA.m,avbeta=avbeta.m));

        }
        else if(dmaMode == "singleProbe"){
            probeID.v <- statM.o$probeID[[c]];
            commonEID.v <- intersect(rownames(adj.m), names(probeID.v));
            mapA.idx <- match(commonEID.v,rownames(adj.m));
            tmpA.m <- adj.m[mapA.idx,mapA.idx];

            mapM.idx <- match(commonEID.v,names(probeID.v));
            tmpM.v <- probeID.v[mapM.idx];

            gr.o <- igraph::graph.adjacency(tmpA.m,mode="undirected");
            comp.l <- igraph::clusters(gr.o);
            ngpc.v <- summary(factor(comp.l$member));
            maxCLid <- as.numeric(names(ngpc.v)[which.max(ngpc.v)]);
            maxc.idx <- which(comp.l$member==maxCLid);

            #get the max connected network
            tmpA.m <- tmpA.m[maxc.idx,maxc.idx];
            gr.o <-  igraph::graph.adjacency(tmpA.m,mode="undirected");
            tmpE.m <- igraph::get.edgelist(gr.o);
            tmpM.v <- tmpM.v[maxc.idx];

            #### now extract statistics
            selcol.idx <- match(c("t","P.Value"),colnames(statM.o$top[[c]]));
            ordrow.idx <- match(names(tmpM.v),rownames(statM.o$top[[c]]));
            if(length(intersect(c("ID"),colnames(statM.o$top[[c]])))==1){
                ordrow.idx <- match(names(tmpM.v),statM.o$top[[c]][,1]);
            }
            statM.m <- statM.o$top[[c]][ordrow.idx,selcol.idx];
            rownames(statM.m) <- names(tmpM.v);

            return(list(statM=statM.m,adj=tmpA.m,probeID=probeID.v));

        }
        else{print("Please indicate correct mode for statM.o object!"); break}
    }

#' DoEpiMod
#'
#' https://rdrr.io/bioc/FEM/src/R/DoEpiMod.R
#'
#'
#' @param intEpi.o
#' @param nseeds
#' @param gamma
#' @param nMC=1000
#' @param sizeR.v
#' @param minsizeOUT
#' @param writeOUT
#' @param nameSTUDY
#' @param ew.v
#' @import igraph
#' @importFrom AnnotationDbi mappedkeys unbox
#' @importFrom org.Hs.eg.db org.Hs.egSYMBOL unbox
#' @return A matrix of the infile
#' @export
DoEpiMod <-
    function(intEpi.o,
             nseeds=100,
             gamma=0.5,
             nMC=1000,
             sizeR.v=c(1,100),
             minsizeOUT=10,
             writeOUT=TRUE,
             nameSTUDY="X",
             ew.v=NULL){

        PasteVector <- function(v){
            vt <- v[1];
            if(length(v) > 1){
                for(g in 2:length(v)){
                    vt <- paste(vt,v[g],sep=" ")

                }
            }
            vt <- paste(vt," EnD",sep="");
            out.v <- sub(" EnD","",vt);
            out.v <- sub("NA , ","",out.v);
            out.v <- sub(" , NA","",out.v);
            out.v <- sub(" , NA , "," , ",out.v);
            return(out.v);
        }

        Heaviside <- function(v){
            out.v <- v;
            out.v[which(v>=0)] <- 1;
            out.v[which(v<0)] <- 0;
            return(out.v);
        }

        WriteOutPval <- function(pv.v,round.min=3,round.max=50){
            round.v <- round.min:round.max
            th.v <- 10^(-round.v);
            outpv.v <- vector(length=length(pv.v));
            done.idx <- which(pv.v >= th.v[1]);
            outpv.v[done.idx] <- round(pv.v[done.idx],round.min);
            todo.idx <- setdiff(1:length(pv.v),done.idx);
            for(i in todo.idx){
                if(length(which(th.v <= pv.v[i]))>0){
                    outpv.v[i] <-
                      round(pv.v[i],round.v[min(which(th.v <= pv.v[i]))]);
                }
                else{
                    outpv.v[i] <- 0;
                }
            }
            return(outpv.v);
        }

        ##require(igraph);
        ##require(org.Hs.eg.db);

        statM.m <- intEpi.o$statM;
        adj.m <- intEpi.o$adj;
        statM.v <- statM.m[,1];
        nameSTUDY <- paste("Epi-",nameSTUDY,sep="");
        statI.v <- abs(statM.v);
        names(statI.v) <- rownames(statM.m);

        x <- org.Hs.eg.db::org.Hs.egSYMBOL;
        mapped_genes <- AnnotationDbi::mappedkeys(x)
        xx <- as.list(x[mapped_genes])
        mapEIDtoSYM.v <- unlist(xx);
        map.idx <- match(rownames(adj.m), names(mapEIDtoSYM.v));
        anno.m <- cbind(rownames(adj.m), mapEIDtoSYM.v[map.idx]);
        colnames(anno.m) <- c("EntrezID","Symbol");

        ### find subnetworks around seeds
        ntop <- nseeds;

        ### use greedy approach: rank nodes to define seeds
        rankN.s <- sort(statI.v,decreasing=TRUE,index.return=TRUE);
        seedsN.v <- names(statI.v)[rankN.s$ix[1:ntop]];

        ### now define weights on network
        print("Constructing weighted network");
        tmpA.m <- adj.m;
        gr.o <-  igraph::graph.adjacency(tmpA.m,mode="undirected");
        tmpE.m <- igraph::get.edgelist(gr.o);
        if(is.null(ew.v)){
            tmpW.v <- vector(length=nrow(tmpE.m));
            for(e in 1:nrow(tmpE.m)){
                match(tmpE.m[e,],rownames(tmpA.m)) -> map.idx;
                tmpW.v[e] <- mean(statI.v[map.idx]);
                print(e);
            }
        }
        else{
            tmpW.v <- ew.v
        }

        ### a number of edges might have a weight of zero which would
        ### later alter the topology of network. this is not desired,
        ### hence we replace 0 by the minimum non-zero value.
        minval <- min(setdiff(tmpW.v,0))
        tmpW.v[tmpW.v==0] <- minval;

        ### defined weighted graph and adjacency matrix
        grW.o <- igraph::set.edge.attribute(gr.o,"weight",value=tmpW.v);
        adjW.m <- igraph::get.adjacency(grW.o,attr="weight")

        ### Run Spin-Glass algorithm
        print("Running Spin-Glass algorithm");
        sizeN.v <- vector();
        sgcN.lo <- list();
        for(v in 1:ntop){
            sgcN.o <- igraph::spinglass.community(gr.o,
                                          weights=tmpW.v,
                                          spins=25,
                                          start.temp=1,
                                          stop.temp=0.1,
                                          cool.fact=0.99,
                                          update.rule=c("config"),
                                          gamma=gamma,
                                          vertex=rankN.s$ix[v]);
            sizeN.v[v] <- length(sgcN.o$comm);
            sgcN.lo[[v]] <- sgcN.o;
            print(paste("Done for seed ",v,sep=""));
        }
        names(sizeN.v) <- seedsN.v;
        print("Module sizes=");
        print(sizeN.v);
        ### compute modularities
        modN.v <- vector();
        for(v in 1:ntop){
            subgr.o <- igraph::induced.subgraph(grW.o,sgcN.lo[[v]]$comm);
            modN.v[v] <- mean(igraph::get.edge.attribute(subgr.o,name="weight"))
        }
        names(modN.v) <- seedsN.v;
        print("Modularity values=");
        print(modN.v);

        ### now determine significance against randomisation of profiles
        print("Starting Monte Carlo Runs");
        modNmc.m <- matrix(nrow=ntop,ncol=nMC);
        for(m in 1:ntop){
            subgr.o <- igraph::induced.subgraph(gr.o,sgcN.lo[[m]]$comm);
            nN <- sizeN.v[m];
            if( (nN> sizeR.v[1]) && (nN< sizeR.v[2])){
                tmpEL.m <- igraph::get.edgelist(subgr.o);
                for(run in 1:nMC){
                    permN.idx <-
                      sample(1:nrow(tmpA.m),nrow(tmpA.m),replace=FALSE);
                    tmpEW.v <- vector();
                    for(e in 1:nrow(tmpEL.m)){
                        match(tmpEL.m[e,],rownames(tmpA.m)[permN.idx]) ->
                          map.idx;
                        tmpEW.v[e] <- mean(statI.v[map.idx]);
                    }
                    subgrW.o <-
                      igraph::set.edge.attribute(subgr.o,"weight",value=tmpEW.v)
                    modNmc.m[m,run] <-
                      mean(igraph::get.edge.attribute(subgrW.o,name="weight"));
                }
            }
            print(paste("Done for seed/module ",m,sep=""));
        }

        modNpv.v <- rep(1,ntop);
        for(v in 1:ntop){
            if( (sizeN.v[v] > sizeR.v[1]) && (sizeN.v[v]< sizeR.v[2])){
                modNpv.v[v] <- length(which(modNmc.m[v,] > modN.v[v]))/nMC;
            }
        }
        names(modNpv.v) <- seedsN.v;
        print(modNpv.v);

        ### summarize hits
        print("Summarising and generating output");
        selpvN.idx <- which(modNpv.v < 0.05);
        selSize.idx <- which(sizeN.v >= minsizeOUT);
        selMod.idx <- intersect(selpvN.idx,selSize.idx);
        print(selMod.idx);
        print(seedsN.v);
        topmodN.m <- matrix(nrow=length(selMod.idx),ncol=6);
        match(seedsN.v[selMod.idx],anno.m[,1]) -> map.idx;
        seedsSYM.v <- anno.m[map.idx,2];

        topmodN.m[,1] <- seedsN.v[selMod.idx];
        topmodN.m[,2] <- seedsSYM.v;
        topmodN.m[,3:5] <- cbind(sizeN.v[selMod.idx],
                                 modN.v[selMod.idx],
                                 modNpv.v[selMod.idx]);
        mi <- 1;
        for(m in selMod.idx){
            tmpEID.v <- rownames(tmpA.m)[sgcN.lo[[m]]$comm];
            genes.v <- anno.m[match(tmpEID.v,anno.m[,1]),2];
            topmodN.m[mi,6] <- PasteVector(genes.v);
            mi <- mi+1;
        }
        colnames(topmodN.m) <- c("EntrezID(Seed)",
                                 "Symbol(Seed)",
                                 "Size",
                                 "Mod",
                                 "P",
                                 "Genes");

        if(writeOUT){
            write.table(topmodN.m,file=paste("topEPI-",
                                             nameSTUDY,
                                             ".txt",
                                             sep=""),
                        quote=FALSE,
                        sep="\t",
                        row.names=FALSE);
        }

        seltopmodN.lm <- list();
        for(m in 1:length(selMod.idx)){
            tmpEID.v <- rownames(tmpA.m)[sgcN.lo[[selMod.idx[m]]]$comm]

            match(tmpEID.v,anno.m[,1]) -> map.idx;

            match(tmpEID.v,rownames(tmpA.m)) -> map1.idx;

            seltopmodN.m <- cbind(anno.m[map.idx,1:2],statM.m[map1.idx,],
                                  statI.v[map1.idx]);

            seltopmodN.lm[[m]] <- seltopmodN.m;

            colnames(seltopmodN.lm[[m]]) <- c("EntrezID",
                                              "Symbol",
                                              "stat(DNAm)",
                                              "P(DNAm)",
                                              "stat(Int)");
        }

        names(seltopmodN.lm) <- seedsSYM.v

        if(writeOUT){

            for(m in 1:length(selMod.idx)){

                out.m <- seltopmodN.lm[[m]];

                out.m[,3] <- round(as.numeric(out.m[,3]),2);

                out.m[,4] <- WriteOutPval(as.numeric(out.m[,4]),round.max=100);

                out.m[,5] <- round(as.numeric(out.m[,5]),2);

                write(paste(seedsSYM.v[m],
                            " (",nrow(seltopmodN.lm[[m]]),
                            " genes)",sep=""),
                      file=paste("topEpiModLists-",
                                 nameSTUDY,".txt",
                                 sep=""),
                      ncolumns=1,
                      append=TRUE);

                write.table(out.m,
                            file=paste("topEpiModLists-",
                                       nameSTUDY,
                                       ".txt",
                                       sep=""),
                            quote=FALSE,
                            sep="\t",
                            row.names=FALSE,
                            append=TRUE);
            }

        }

        return(list(size=sizeN.v,
                    mod=modN.v,
                    pv=modNpv.v,
                    selmod=selMod.idx,
                    fem=topmodN.m,
                    topmod=seltopmodN.lm,
                    sgc=sgcN.lo,
                    ew=tmpW.v,
                    adj=intEpi.o$adj));

    }

#' FemModShow
#'
#' https://rdrr.io/bioc/FEM/src/R/FemModShow.R
#'
#'
#' @param mod
#' @param name
#' @param fem.o
#' @param mode
#' @import igraph
#' @importFrom shape colorlegend unbox
#' @importFrom marray maPalette unbox
#' @return A matrix of the infile
#' @export
FemModShow <-
    function(mod,
             name,
             fem.o,
             mode= "Integration"
             ){

        edgeweight=fem.o$ew
        adjacency=fem.o$adj

        ###############add shape to the vertex of igraph

        mycircle <- function(coords, v=NULL, params) {
            vertex.color <- S4Vectors::params("vertex", "color")
            if (length(vertex.color) != 1 && !is.null(v)) {
                vertex.color <- vertex.color[v]
            }
            vertex.size  <- 1/200 * S4Vectors::params("vertex", "size")
            if (length(vertex.size) != 1 && !is.null(v)) {
                vertex.size <- vertex.size[v]
            }
            vertex.frame.color <- S4Vectors::params("vertex", "frame.color")
            if (length(vertex.frame.color) != 1 && !is.null(v)) {
                vertex.frame.color <- vertex.frame.color[v]
            }
            vertex.frame.width <- S4Vectors::params("vertex", "frame.width")
            if (length(vertex.frame.width) != 1 && !is.null(v)) {
                vertex.frame.width <- vertex.frame.width[v]
            }

            mapply(coords[,1],
                   coords[,2],
                   vertex.color,
                   vertex.frame.color,
                   vertex.size,
                   vertex.frame.width,
                   FUN=function(x, y, bg, fg, size, lwd) {
                       symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd,
                               circles=size, add=TRUE, inches=FALSE)
                   })
        }

        #mod is from fem result such as HAND2,
        #hand2<-fembi.o$topmod$HAND2.
        #overall ideas: give the mtval, rtval,
        #vmcolor, vrcolor, label to vertex;
        #give weight, edgewidth to the edge of
        #realgraph. then subgraph
        #the HAND2 mod

        realgraph = igraph::graph.adjacency(adjacency,mode="undirected")

        #add the values
        igraph::E(realgraph)$weight = edgeweight #give the weight to the edges

        ############################################################
        # 1 edge width size

        edge.width.v = igraph::E(realgraph)$weight

        idxbw02=which(edge.width.v<=2 && edge.width.v>0)
        idxbw25=which(edge.width.v>2 && edge.width.v<5)
        idxbw510=which(edge.width.v>=5 && edge.width.v<10)
        idxgt10=which(edge.width.v>=10)

        edge.width.v[idxbw02]=1/4*edge.width.v[idxbw02]
        edge.width.v[idxbw25]=1/2*edge.width.v[idxbw25]
        edge.width.v[idxbw510]=3/4*edge.width.v[idxbw510]
        edge.width.v[idxgt10]=10#the max edgewidth is fixed as 10

        idxlt025=which(edge.width.v<0.25)
        edge.width.v[idxlt025]=0.25
        # if the edgewidth is less than 0.25,
        # it's too narrow to see. So fix the them with 0.25

        igraph::E(realgraph)$edgewidth=edge.width.v

        ##################################################################
        # 2 and 3  transform methylation value to node color, rna expression
        mod=as.data.frame(mod);
        if(mode=="Epi"){
            mod[,"stat(mRNA)"]=mod[,"stat(DNAm)"];
        }else if(mode=="Exp"){
            mod[,"stat(DNAm)"]=mod[,"stat(mRNA)"];
        }
        #if  there is mode is epi we add the stat(mRNA) colum and set them as 0;
        #sub graph mod and inhebited the edgewidth add the mval and rval
        mod.graph=igraph::induced.subgraph(realgraph,v=as.vector(mod[,1]))
        print(mod[,1])
        # the fembi.o$topmod$HAND2 is a dataframe and it has "Symbol",
        # "stat(DNAm)"
        # "stat(mRNA)"
        # which is useful
        mtval.v=vector()
        for(i in igraph::V(mod.graph)$name){
            mtval.v=c(mtval.v,(as.vector(mod[i,"stat(DNAm)"])))}
        # add the mtval to the mod graph one by one
        # according to the squence of mod graph name
        igraph::V(mod.graph)$mtval=mtval.v;

        rtval.v=vector()
        for(i in igraph::V(mod.graph)$name){
            rtval.v=c(rtval.v,(as.vector(mod[i,"stat(mRNA)"])))}
        # add the
        igraph::V(mod.graph)$rtval=rtval.v;
        print(rtval.v)

        # add the vm.color, vr.color
        vm.color=rep(0,length(igraph::V(mod.graph)));
        vr.color=rep(0,length(igraph::V(mod.graph)));
        # color scheme generate 100 colors
        tmcolor.scheme<-marray::maPalette(low = "yellow",
                                  high ="blue",
                                  mid="grey",
                                  k =100);
        trcolor.scheme<-marray::maPalette(low = "green",
                                  high ="red",
                                  mid="grey",
                                  k =100);

        tmcolor.scheme[13:88]="#BEBEBE";
        #( floor(-1.5/0.04)+51,floor(1.5/0.04)+51)
        #which is [13,88] is grey "#BEBEBE". 1.5 is
        #the thresh hold for the t values
        tmcolor.scheme[1:12]=marray::maPalette(low = "yellow",
                                       high="lightyellow",
                                       k =12);
        tmcolor.scheme[89:100]=marray::maPalette(low = "lightblue",
                                         high="blue",
                                         k =12);
        trcolor.scheme[13:88]="#BEBEBE";
        #[ floor(-1.5/0.04)+51,floor(1.5/0.04)+51]
        #which is [13,88] is grey "#BEBEBE". 1.5 is
        #the thresh hold for the t values
        trcolor.scheme[1:12]=marray::maPalette(low = "green",
                                       high="lightgreen",
                                       k =12);
        trcolor.scheme[89:100]=marray::maPalette(low = "#DC6868",
                                         high="red",
                                         k =12);
        # give the color according the mtval.
        # (-2,2),floor get the integer mod on 0.04 + 51
        # then get the according color.
        tmcolor.position=floor(as.numeric(igraph::V(mod.graph)$mtval)/0.04)+51;
        tmcolor.position[which(tmcolor.position<1)]<-1;
        tmcolor.position[which(tmcolor.position>100)]<-100;
        igraph::V(mod.graph)$vmcolor<-vm.color
        ##add the vm.color to the vertex value
        vm.color=tmcolor.scheme[tmcolor.position];
        print(vm.color)

        # add the frame color idea: get the tr color position then get
        # the color from trcolor.scheme
        trcolor.position=floor(as.numeric(igraph::V(mod.graph)$rtval)/0.04)+51;
        print(trcolor.position)
        trcolor.position[which(trcolor.position<1)]<-1;
        trcolor.position[which(trcolor.position>100)]<-100;
        vr.color=trcolor.scheme[trcolor.position];

        igraph::V(mod.graph)$vrcolor<-vr.color  ## the rna expression color

        if(mode=="Exp"){
            igraph::V(mod.graph)$color<-igraph::V(mod.graph)$vrcolor;
        }else{
            igraph::V(mod.graph)$color<-igraph::V(mod.graph)$vmcolor
            #use the vmcolor as the vertex color but if the mod is Exp
            #them the vertex color is vrcolor
        }
        print(vr.color)
        ####################################################################
        #add the mod label value

        label.v=vector()
        for(i in igraph::V(mod.graph)$name){
            label.v=c(label.v,(as.vector(mod[i,"Symbol"])))}
        #add the V(mod.graph)$name's
        #labels one by one from mod["$name","Symbol"]

        #####################################################################
        #create subgraph label.cex value
        igraph::V(mod.graph)$label.cex=rep(0.5,
                                      length(as.vector(igraph::V(mod.graph))));

        #all the cex first set as 0.7
        igraph::V(mod.graph)$label.cex[which(
            as.vector(igraph::V(mod.graph)$name)==as.vector(mod[1,1]))]=0.8
        #only the firt mod name was set as 1

        ####################################################################
        #generate the plot
        #when you want to plot the vertex shape,
        #and its frame width first you should load the api script api bellow
        igraph::add.vertex.shape("fcircle",
                         clip=igraph::igraph.shape.noclip,
                         plot=mycircle,
                         parameters=list(vertex.frame.color=1,
                                         vertex.frame.width=1))

        pdf(paste(name,".mod.pdf",sep=""))
        if(mode =="Integration"){

            plot(mod.graph,
                 layout=igraph::layout.fruchterman.reingold,
                 vertex.shape="fcircle",
                 vertex.frame.color=igraph::V(mod.graph)$vrcolor,
                 vertex.frame.width=4,
                 vertex.size=10,
                 vertex.label=label.v,
                 vertex.label.dist=0.6,
                 vertex.label.cex=igraph::V(mod.graph)$label.cex,
                 vertex.label.font=3,
                 edge.color="grey",
                 edge.width=igraph::E(mod.graph)$edgewidth)

            shape::colorlegend(trcolor.scheme,
                        seq(-2,2,0.5),
                        ratio.colbar=0.3,
                        xlim=c(-1.55,-1.4),
                        ylim=c(-0.5,0),
                        align="r",
                        cex=0.5)

            shape::colorlegend(tmcolor.scheme,
                        seq(-2,2,0.5),
                        ratio.colbar=0.3,
                        xlim=c(-1.55,-1.4),
                        ylim=c(0.5,1),
                        align="r",
                        cex=0.5)

            text(-1.50, 0.43, c("t(DNAm)\nCore"), cex=0.6)
            text(-1.50, -0.57,c("t(mRNA)\nBorder"), cex=0.6)

        }
        else if(mode =="Epi"){
            # if the mode is Epi the frame need not to show
            plot(mod.graph,
                 layout=igraph::layout.fruchterman.reingold,
                 vertex.frame.color=NA,
                 vertex.size=10,
                 vertex.label=label.v,
                 vertex.label.dist=0.6,
                 vertex.label.cex=igraph::V(mod.graph)$label.cex,
                 vertex.label.font=3,
                 edge.color="grey",
                 edge.width=igraph::E(mod.graph)$edgewidth)

            shape::colorlegend(tmcolor.scheme,
                        seq(-2,2,0.5),
                        ratio.colbar=0.3,
                        xlim=c(-1.55,-1.4),
                        ylim=c(0.5,1),
                        align="r",
                        cex=0.5)
            text(-1.50, 0.43, c("t(DNAm)"),cex=0.6)

        }
        else if(mode =="Exp"){

            plot(mod.graph,
                 layout=igraph::layout.fruchterman.reingold,
                 vertex.frame.color=NA,
                 vertex.size=10,
                 vertex.label=label.v,
                 vertex.label.dist=0.6,
                 vertex.label.cex=igraph::V(mod.graph)$label.cex,
                 vertex.label.font=3,
                 edge.color="grey",
                 edge.width=igraph::E(mod.graph)$edgewidth)

            shape::colorlegend(trcolor.scheme,
                        seq(-2,2,0.5),
                        ratio.colbar=0.3,
                        xlim=c(-1.55,-1.4),
                        ylim=c(-0.5,0),
                        align="r",
                        cex=0.5)

            text(-1.50, -0.57,c("t(mRNA)"),cex=0.6)

        }
        dev.off()
        return(igraph::igraph.to.graphNEL(mod.graph));

}

#' @title champ.EpiMod
#' @description champ.EpiMod functionality

#' BiocManager::install("FEM")
#' https://rdrr.io/github/gaberosser/ChAMP/src/R/champ.EpiMod.R
#'
#'
#' @param beta
#' @param pheno
#' @param nseeds
#' @param gamma
#' @param nMC
#' @param sizeR.v
#' @param minsizeOUT
#' @param resultsDir
#' @param PDFplot=TRUE
#' @param arraytype
#' @return A matrix of the infile
#' @export
champ.EpiMod <- function(beta=myNorm,
                         pheno=myLoad$pd$Sample_Group,
                         nseeds=100,
                         gamma=0.5,
                         nMC=1000,
                         sizeR.v=c(1,100),
                         minsizeOUT=10,
                         resultsDir="./CHAMP_EpiMod/",
                         PDFplot=TRUE,
                         arraytype="450K")
{

    if(getRversion() >= "3.1.0") utils::globalVariables(c("myNorm",
                                                        "myLoad",
                                                        "hprdAsigH.m")
    )
    message("[===========================]")
    message("[<<< ChAMP.EpiMod START >>>>]")
    message("-----------------------------")

    ### Prepare Checking ###
    if (!file.exists(resultsDir)) dir.create(resultsDir)
    message("champ.EpiMod Results will be saved in ",resultsDir)

    data(hprdAsigH)
    message("<< Load PPI network hprdAsigH >>")

    if(arraytype=="EPIC")
        statM.o <- GenStatM(beta,pheno,arraytype)
    else
        statM.o <- GenStatM(beta,pheno,"450k")
    message("<< Generate statM.o >>")

    intEpi.o=DoIntEpi450k(statM.o,hprdAsigH.m,c=1)

    message("<< Calculate EpiMod.o >>")
    EpiMod.o=DoEpiMod(intEpi.o,
                      nseeds=nseeds,
                      gamma=gamma,
                      nMC=nMC,
                      sizeR.v=sizeR.v,
                      minsizeOUT=minsizeOUT,
                      writeOUT=TRUE,
                      ew.v=NULL);

    if(PDFplot)
    {
        message("<< Draw All top significant module plot in PDF >>")
        tmpdir <- getwd()
        setwd(resultsDir)
        for(i in names(EpiMod.o$topmod)) FemModShow(EpiMod.o$topmod[[i]],
                                                    name=i,
                                                    EpiMod.o,
                                                    mode="Epi")
        setwd(tmpdir)
    }

    message("[<<<< ChAMP.EpiMod END >>>>>]")
    message("[===========================]")
    return(EpiMod.o)
}
