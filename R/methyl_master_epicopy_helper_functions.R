#!/usr/bin/env Rscripts

#' @title .coerce.pData.mm
#' @description MethylMaster .coerce.pData() function
#' This will fix the issue with LRR2
#' @param pdat Path to the input file
#' @return Compatible pData data frame
#' @export
.coerce.pData.mm <- function(pdat){
  pDat <- apply(pdat, 2, as.character)
  dimnames(pDat) <- dimnames(pdat)
  pDat <- as(pDat,
             "DataFrame"##,
             ##stringsAsFactors = FALSE ##Will cause error
  )
  return(pDat)
}

#' @title LRRtoCNA.mm
#' @description MethylMaster version of the Epicopy LRRtoCNA function
#' Epicopy convert log r ratio data to cna calls:
#' Similar to Lucas Salas' PhD modified LRRtoCNA2 function
#' @param LRR The LRR parameter
#' @param sampID The sampID parameter
#' @param ncores The ncores parameter, note that multiple cores may not work
#' on every system
#' @param seg.undo.splits The seg.undo.splits parameter
#' @param seg.undo.SD The seg.undo.SD parameter
#' @param epic.manifest.path The epic.manifest.path parameter
#' @param platform The platform parameter
#' @import minfi
#' @return Returns CNA segments
#' @export
LRRtoCNA.mm <- function(LRR,
                      sampID = NULL,
                      ncores = 1L,
                      seg.undo.splits = "sdundo",
                      seg.undo.SD = 2,
                      epic.manifest.path = getwd(),
                      platform=c("hm450","EPIC"),
                      ...) {
  if(platform=="EPIC"){
    load(file=epic.manifest.path) ##"EPIC.manifest.hg38.rda"
    EPIC.manifest.hg38$probeTarget <-
      GenomicRanges::start(EPIC.manifest.hg38)
    hm450 <- EPIC.manifest.hg38
  }else{
    hm450 <- Epicopy:::hm450
  }

  if(!all(rownames(LRR) %in% names(hm450))){
    stop(paste0("Some probe names in LRR ",
                "do not exist in 450K marker file."))
  }

  lrr.map = hm450[rownames(LRR)]
  lrr.chrom = as.numeric(gsub("chr",
                              "",
                              seqnames(lrr.map)))
  lrr.loc = lrr.map$probeTarget
  if(is.null(sampID)){
    sampID = colnames(LRR)
  }
  cat("Getting CNA object.\n")
  cna = DNAcopy::CNA(genomdat = LRR,
                     chrom = lrr.chrom,
                     maploc = lrr.loc,
                     data.type = "logratio",
                     sampleid = sampID)
  cat("Smoothing CNA object.\n")
  scna = DNAcopy::smooth.CNA(cna)
  if (!inherits(ncores, c("numeric", "integer"))) {
    stop("ncores must be numeric/integer")
  }
  ncores = as.integer(ncores)
  if (ncores == 1) {
    cat("Running segmentation with 1 core.\n")
    seg = DNAcopy::segment(scna,
                           undo.splits = seg.undo.splits,
                           undo.SD = seg.undo.SD,
                           ...)
  }
  else if(ncores > 1){
    cat("Running segmentation with",
        ncores,
        "cores.\n")
    seg = Epicopy:::parSegment(CNAobj = scna,
                               undo.splits = seg.undo.splits,
                               undo.SD = seg.undo.SD,
                               njobs = ncores, ...)
  }
  else{
    stop("ncores is not specified correctly.")
  }
  return(seg)
}

#' @title .funnorm.mm
#' @description MethylMaster version of the Epicopy .funnorm() function
#' @param rgSet The rgSet parameter
#' @param nPCs The nPCs parameter
#' @param sex The sex parameter
#' @param verbose the verbose parameter
#' @import minfi
#' @import BiocGenerics
#' @return Return a final gm set
#' @export
.funnorm.mm <- function(rgSet,
                        nPCs = 2,
                        sex = NULL,
                        verbose = TRUE){
  minfi:::.isRGOrStop(rgSet)
  ##rgSet <- BiocGenerics::updateObject(rgSet)
  if(verbose)
    cat("[preprocessFunnorm] Mapping to genome\n")

  gmSet <- minfi::mapToGenome(rgSet)
  subverbose <- max(as.integer(verbose) - 1L, 0)

  if (verbose)
    cat("[preprocessFunnorm] Quantile extraction\n")

  extractedData <- minfi:::.extractFromRGSet450k(rgSet)

  if (is.null(sex)) {
    gmSet <- minfi::addSex(gmSet, getSex(gmSet, cutoff = -3))
    sex <- rep(1L, length(gmSet$predictedSex))
    sex[gmSet$predictedSex == "F"] <- 2L
  }
  rm(rgSet)
  if (verbose)
    cat("[preprocessFunnorm] Normalization\n")
  CN <- minfi::getCN(gmSet)
  gmSet <- minfi:::.normalizeFunnorm450k(object = gmSet,
                                         extractedData = extractedData,
                                         sex = sex,
                                         nPCs = nPCs,
                                         verbose = subverbose)
  return(gmSet)

}

#' @title getLRR.mm
#' @description MethylMaster version of Epicopy getLRR function
#' @param rgSet The rgSet parameter
#' @param Normals The Normals parameter
#' @param sampNames The sampNames parameter
#' @param QN The QN parameter
#' @param Ref The Ref parameter
#' @param mode.bw The mode.bw parameter
#' @param mode.method The mode.method parameter
#' @param normal.cnv The normal.cnv parameter
#' @param mean.center The mean.center parameter
#' @param filter.autosomal.probes The filter.autosomal.probes parameter
#' @param autosomal The autosomal parameter
#' @param filterProbes The filterProbes parameter
#' @param ... Additional parameter passed to getLRR.mm
#' @param import EpicopyData
#' @import stats
#' @import minfi
#' @return Return final LRR values
#' @export
getLRR.mm <- function(rgSet = rgset,
                      Normals = Normals,
                      sampNames = sampNames,
                      QN = FALSE,
                      Ref = "median",
                      mode.bw = 0.1,
                      mode.method = "naive",
                      normal.cnv = TRUE,
                      mean.center = TRUE,
                      filter.autosomal.probes = FALSE,
                      autosomal = NULL,
                      filterProbes = FALSE,
                      keep_fnobj = FALSE,
                     ...
                     ){
##debug: {
  if (!inherits(rgSet, "RGChannelSet")) {
    stop("rgset has to be RGChannelSet.")
  }
  if (!inherits(Normals, c("numeric", "logical",
                           "RGChannelSet", "NULL", "character"))) {
    stop("Normals in wrong format.")
  }
  if (!Ref %in% c("mode", "median")) {
    stop("Ref has to be 'mode' or 'median'.")
  }
  pData(rgSet) <- .coerce.pData.mm(pData(rgSet))
  if (!is.null(Normals) & length(Normals) == 1)
    if (is.na(Normals)) {
      cat(paste0("Normals specified as NA. ",
                 "Reference will calculated ",
                 "using all user provided samples.\n"))
      Normals <- rep(TRUE, ncol(rgSet))
      normal.cnv <- TRUE
    }
  if (inherits(Normals, c("numeric", "logical"))) {
    ##if (nrow(rgSet) != nrow(Normals)){
    ##  stop("Features in normal and tumor set do not match.")
    ##}else{
    rgset <- rgSet
    Normals <- Normals
    multi.platform <- FALSE
    ##}
  }
  if (inherits(Normals, "RGChannelSet")) {
    if (nrow(rgSet) != nrow(Normals)){
      ##if(all.equal(colnames(rgSet),colnames(Normals))==TRUE){
      ##  stop("Error: treatment and normal samples are the same samples!")
      ##}
      stop("Features in normal and tumor set do not match.")
      ##multi.platform <- TRUE
      ##combined.colnames <- c(colnames(rgSet),colnames(Normals))
      #norm.rgSet <- Normals
      ##Normals <- combined.colnames  %in% colnames(Normals)
    }else{
      rgset <- BiocGenerics::combine(rgSet, Normals)
      Normals <- colnames(rgset) %in% colnames(Normals)
      multi.platform <- FALSE
    }
  }
  if (is.null(Normals)) {
    if ("EpicopyData" %in% installed.packages()) {
      library(EpicopyData)
      cat("Using all epicopy normals as reference set.\n")
      if (!"all.normals" %in% ls(globalenv()))
        data("allNormals")
      rgset <- BiocGenerics::combine(rgSet, all.normals)
      Normals <- colnames(rgset) %in% colnames(all.normals)
    }
    else {
      cat(paste0("EpicopyData not installed. ",
                 "Reference will calculated ",
                 "using all user provided samples.\n"))
      Normals <- rep(TRUE, ncol(rgSet))
      normal.cnv <- TRUE
      rgset <- rgSet
    }
  }
  if (inherits(Normals, "character")) {
    if (!Normals %in% c("all", "thyroid", "breast",
                        "lung")) {
      stop(paste0("Character normals have to be one of ",
                  "'all', 'thyroid', 'breast', or 'lung'."))
    }
    .Normals <- Normals
    if (Normals == "all") {
      library(EpicopyData)
      if (!"all.normals" %in% ls(globalenv())){
      ##data("allNormals")
      all.normals <- updateObject(data("allNormals"))
      Normals <- all.normals
      }
    }
    else if (Normals == "thyroid") {
      library(EpicopyData)
      if (!"all.normals" %in% ls(globalenv()))
        data("allNormals")
      Normals <- all.normals[, pData(all.normals)$Tissue_Type %in%
                               "thyroid"]
    }
    else if (Normals == "breast") {
      library(EpicopyData)
      if (!"all.normals" %in% ls(globalenv()))
        data("allNormals")
      Normals <- all.normals[, pData(all.normals)$Tissue_Type %in%
                               "breast"]
    }
    else if (Normals == "lung") {
      library(EpicopyData)
      if (!"all.normals" %in% ls(globalenv()))
        data("allNormals")
      Normals <- all.normals[, pData(all.normals)$Tissue_Type %in%
                               "lung"]
    }
    cat("Using", .Normals, "normals as reference set.\n")
    pData(Normals) <- .coerce.pData.mm(pData(Normals))
    rgset <- BiocGenerics::combine(rgSet, Normals)
    Normals <- colnames(rgset) %in% colnames(Normals)
  }

  Normals <- as.logical(Normals)
  if (any(is.na(Normals))) {
    stop("Error in classifying normals. Please review input.")
  }

  if(multi.platform==TRUE){

    if (is.null(sampNames)) {
      colnames(rgSet) <- pData(rgSet)$Sample_Name
      colnames(norm.rgSet) <- pData(norm.rgSet)$Sample_Name
    }
    else {
      if (ncol(rgSet) == length(sampNames))
        stop("sampNames is not equal to number of samples in rgset.\n")
      colnames(rgSet) <- sampNames
    }
    cat("Performing functional normalization.\n")
    fn.set <- .funnorm.mm(rgSet=rgSet)
    fn.intensity <- minfi::getUnmeth(fn.set) + minfi::getMeth(fn.set)
    fn.intensity[fn.intensity <= 0] <- 0.01
    fn.intensity <- log2(fn.intensity)

    fn.set.norm <- .funnorm.mm(rgSet=norm.rgSet)
    fn.intensity.norm <- minfi::getUnmeth(fn.set.norm) +
      minfi::getMeth(fn.set.norm)
    fn.intensity.norm[fn.intensity.norm <= 0] <- 0.01
    fn.intensity.norm <- log2(fn.intensity.norm)
    normal.fn <- fn.intensity.norm

    if (keep_fnobj) {
      cat("Saving funnorm objects regular and normals to specified directory.")
      if (is.null(fn_output))
        fn_output <- "."
      fn_name <- paste0(fn_output, "/epicopy_fn.rda")
      save(fn.set, file = fn_name)
      fn_name.norm <- paste0(fn_output, "/epicopy_fn_norm.rda")
      save(fn.set.norm, file = fn_name_norm)
    }
    if (QN) {
      cat("Quantile normalization of normals.\n")
      normal.t1g <- normal.fn[t1g.probes, ]
      normal.t1r <- normal.fn[t1r.probes, ]
      normal.t2 <- normal.fn[t2.probes, ]
      stopifnot(stats::complete.cases(normal.t1g),
                stats::complete.cases(normal.t1r),
                stats::complete.cases(normal.t2))
      normal.t1g.qn <- normalize.quantiles(normal.t1g)
      dimnames(normal.t1g.qn) <- dimnames(normal.t1g)
      normal.t1r.qn <- normalize.quantiles(normal.t1r)
      dimnames(normal.t1r.qn) <- dimnames(normal.t1r)
      normal.t2.qn <- normalize.quantiles(normal.t2)
      dimnames(normal.t2.qn) <- dimnames(normal.t2)
      normal.fn <- rbind(normal.t1g.qn, normal.t1r.qn, normal.t2.qn)
    }
    if(filter.autosomal.probes==TRUE){
      if(is.null(autosomal)){
        stop("Error: need to set <autosomal> parameter")
      }
      if (!setequal(autosomal, rownames(normal.fn))){
        stop("Normalization returns discordant number of probes.")
      }
      cat("Filtering for autosomal probes.\n")
      normal.fn <- normal.fn[autosomal, ]
      tumor.fn <- fn.intensity[autosomal, !Normals]
    }

    if (filterProbes) {
      if (is.null(retainedProbes)) {
        cat("Filtering for TCGA probeset.")
        normal.fn <- normal.fn[tcga.probeset, ]
        tumor.fn <- tumor.fn[tcga.probeset, ]
      }
      else {
        if (!inherits(retainedProbes, "character"))
          stop("RetainedProbes has to be of\n                                                      class character.\n")
        if (!all(retainedProbes %in% rownames(tumor.fn)))
          stop("Some retainedProbes are\n                                                            not found in dataset.\n")
        cat("Filtering for user selected probeset.")
        normal.fn <- normal.fn[retainedProbes, ]
        tumor.fn <- tumor.fn[retainedProbes, ]
      }
    }
    cat("Calculating reference intensity using", Ref, ".\n")
    if (Ref == "mode") {
      normal.ref <- apply(normal.fn, 1, function(x) {
        modeest::mlv(x, bw = mode.bw, method = mode.method,
                     na.rm = TRUE)$M
      })
    }
    if (Ref == "median") {
      normal.ref <- apply(normal.fn, 1, median, na.rm = TRUE)
    }
    cat("Obtaining LRR.\n")
    tumor.lrr <- sweep(tumor.fn, 1, normal.ref, "-")
    if (!normal.cnv) {
      final.lrr <- tumor.lrr
    }
    else {
      normal.lrr <- sweep(normal.fn, 1, normal.ref, "-")
      final.lrr <- cbind(tumor.lrr, normal.lrr)
    }
    if (mean.center) {
      cat("Mean centering data.\n")
      lrr.mean <- apply(final.lrr, 2, mean)
      final.lrr <- sweep(final.lrr, 2, lrr.mean, "-")
    }
    return(final.lrr)

  }else{
    if (is.null(sampNames)) {
      colnames(rgset) <- pData(rgset)$Sample_Name
    }
    else {
      if (ncol(rgset) == length(sampNames))
        stop("sampNames is not equal to number of samples in rgset.\n")
      colnames(rgset) <- sampNames
    }
    cat("Performing functional normalization.\n")
    fn.set <- .funnorm.mm(rgSet=rgset, ...)
    fn.intensity <- minfi::getUnmeth(fn.set) + minfi::getMeth(fn.set)
    fn.intensity[fn.intensity <= 0] <- 0.01
    fn.intensity <- log2(fn.intensity)
    normal.fn <- fn.intensity[, Normals]
    if (keep_fnobj) {
      cat("Saving funnorm object with normals to specified directory.")
      if (is.null(fn_output))
        fn_output <- "."
      fn_name <- paste0(fn_output, "/epicopy_fn.rda")
      save(fn.set, file = fn_name)
    }
    if (QN) {
      cat("Quantile normalization of normals.\n")
      normal.t1g <- normal.fn[t1g.probes, ]
      normal.t1r <- normal.fn[t1r.probes, ]
      normal.t2 <- normal.fn[t2.probes, ]
      stopifnot(stats::complete.cases(normal.t1g),
              stats::complete.cases(normal.t1r),
              stats::complete.cases(normal.t2))
      normal.t1g.qn <- normalize.quantiles(normal.t1g)
      dimnames(normal.t1g.qn) <- dimnames(normal.t1g)
      normal.t1r.qn <- normalize.quantiles(normal.t1r)
      dimnames(normal.t1r.qn) <- dimnames(normal.t1r)
      normal.t2.qn <- normalize.quantiles(normal.t2)
      dimnames(normal.t2.qn) <- dimnames(normal.t2)
      normal.fn <- rbind(normal.t1g.qn, normal.t1r.qn, normal.t2.qn)
      if (!setequal(autosomal, rownames(normal.fn)))
        stop("Normalization returns discordant number of probes.")
    }
    if(filter.autosomal.probes==TRUE){
      if(is.null(autosomal)){
        stop("Error: you must select <autosomal> probes")
      }
      cat("Filtering for autosomal probes.\n")
      normal.fn <- normal.fn[autosomal, ]
      tumor.fn <- fn.intensity[autosomal, !Normals]
    }else{
      tumor.fn <- fn.intensity[, !Normals]
    }
    if (filterProbes) {
      if (is.null(retainedProbes)) {
        cat("Filtering for TCGA probeset.")
        normal.fn <- normal.fn[tcga.probeset, ]
        tumor.fn <- tumor.fn[tcga.probeset, ]
      }
      else {
        if (!inherits(retainedProbes, "character"))
          stop("RetainedProbes has to be of\n                                                      class character.\n")
        if (!all(retainedProbes %in% rownames(tumor.fn)))
          stop("Some retainedProbes are\n                                                            not found in dataset.\n")
        cat("Filtering for user selected probeset.")
        normal.fn <- normal.fn[retainedProbes, ]
        tumor.fn <- tumor.fn[retainedProbes, ]
      }
    }
    cat("Calculating reference intensity using", Ref, ".\n")
    if (Ref == "mode") {
      normal.ref <- apply(normal.fn, 1, function(x) {
        modeest::mlv(x, bw = mode.bw, method = mode.method,
                   na.rm = TRUE)$M
      })
    }
    if (Ref == "median") {
      normal.ref <- apply(normal.fn, 1, median, na.rm = TRUE)
    }
    cat("Obtaining LRR.\n")
    tumor.lrr <- sweep(tumor.fn, 1, normal.ref, "-")
    if (!normal.cnv) {
      final.lrr <- tumor.lrr
    }
    else {
      normal.lrr <- sweep(normal.fn, 1, normal.ref, "-")
      final.lrr <- cbind(tumor.lrr, normal.lrr)
    }
    if (mean.center) {
      cat("Mean centering data.\n")
      lrr.mean <- apply(final.lrr, 2, mean)
      final.lrr <- sweep(final.lrr, 2, lrr.mean, "-")
    }
    return(final.lrr)
  }
}

#' @title readmetharray.mm
#' @description MethylMaster version of Minfi read.metharray() function
#' @param basenames The basenames parameter
#' @param extended The extended parameter
#' @param verbose The verbose parameter
#' @param force The force parameter
#' @param diff.array.type The diff.array.type parameter
#' @import DelayedArray
#' @import BiocParallel
#' @import illuminaio
#' @import ENmix
#' @import minfi
#' @return Methylation array data
#' @export
read.metharray.mm <- function(basenames,
                              extended = FALSE,
                              verbose = FALSE,
                              force = FALSE){
  BACKEND <- DelayedArray::getAutoRealizationBackend()
  BPREDO <- list()
  BPPARAM <- BiocParallel::SerialParam()
  basenames <- sub("_Grn\\.idat.*", "", basenames)
  basenames <- sub("_Red\\.idat.*", "", basenames)
  stopifnot(!anyDuplicated(basenames))
  G.files <- paste(basenames, "_Grn.idat", sep = "")
  names(G.files) <- basename(basenames)
  these.dont.exists <- !file.exists(G.files)
  if (any(these.dont.exists)) {
    G.files[these.dont.exists] <- paste0(G.files[these.dont.exists],
                                         ".gz")
    ##G.files[these.dont.exists] <- G.files[these.dont.exists]
  }
  if (!all(file.exists(G.files))) {
    noexits <- sub("\\.gz", "", G.files[!file.exists(G.files)])
    ##noexits <- G.files[!file.exists(G.files)]
    stop("The following specified files do not exist:",
         paste(noexits, collapse = ", "))
  }
  R.files <- paste(basenames, "_Red.idat", sep = "")
  names(R.files) <- basename(basenames)
  these.dont.exists <- !file.exists(R.files)
  if (any(these.dont.exists)) {
    R.files[these.dont.exists] <- paste0(R.files[these.dont.exists],
                                         ".gz")
  }
  if (!all(file.exists(R.files))) {
    noexits <- R.files[!file.exists(G.files)]
    stop("The following specified files do not exist:",
         paste(noexits, collapse = ", "))
  }
  stime <- system.time({
    G.Quants <- BiocParallel::bplapply(G.files, function(xx) {
      if (verbose)
        message("[read.metharray] Reading ", basename(xx))
      Quants <- illuminaio::readIDAT(xx)[["Quants"]]
      if (!extended) {
        Quants <- Quants[, "Mean", drop = FALSE]
      }
      if (!is.null(BACKEND)) {
        Quants <- realize(Quants, BACKEND = BACKEND)
      }
      Quants
    }, BPREDO = BPREDO, BPPARAM = BPPARAM)
    R.Quants <- BiocParallel::bplapply(R.files, function(xx) {
      if (verbose)
        message("[read.metharray] Reading ", basename(xx))
      Quants <- illuminaio::readIDAT(xx)[["Quants"]]
      if (!extended) {
        Quants <- Quants[, "Mean", drop = FALSE]
      }
      if (!is.null(BACKEND)) {
        Quants <- realize(Quants, BACKEND = BACKEND)
      }
      Quants
    }, BPREDO = BPREDO, BPPARAM = BPPARAM)
  })[3]
  if (verbose) {
    message(sprintf("[read.metharray] Read idat files in %.1f seconds",
                    stime))
  }
  if (verbose) {
    message("[read.metharray] Creating data matrices ... ",
            appendLF = FALSE)
  }
  ptime1 <- proc.time()
  allNProbes <- vapply(G.Quants, nrow, integer(1L))
  arrayTypes <- cbind(do.call(rbind, lapply(allNProbes,
                                            ENmix:::.guessArrayTypes)),
                      size = allNProbes)
  sameLength <- (length(unique(arrayTypes[, "size"])) == 1)
  sameArray <- (length(unique(arrayTypes[, "array"])) == 1)
    if(!sameLength && !sameArray){
    cat("[read.metharray] Trying to parse IDAT files from different arrays.\n")
      cat("  Inferred Array sizes and types:\n")
      print(arrayTypes[, c("array", "size")])
      stop("[read.metharray] Trying to parse different IDAT files, of ",
          "different size and type.")
    }
    if(!sameLength && sameArray && !force){
      stop("[read.metharray] Trying to parse IDAT files with different ",
          "array size but seemingly all of the same type.\n  You can force ",
          "this by 'force=TRUE', see the man page ?read.metharray")
    }
  commonAddresses <- as.character(Reduce("intersect", lapply(G.Quants,
                                                             rownames)))
  GreenMean <- do.call(cbind, lapply(G.Quants, function(xx) xx[commonAddresses,
                                                               "Mean"]))
  colnames(GreenMean) <- names(G.Quants)
  RedMean <- do.call(cbind, lapply(R.Quants, function(xx) xx[commonAddresses,
                                                             "Mean"]))
  colnames(RedMean) <- names(R.Quants)
  if (extended) {
    GreenSD <- do.call(cbind, lapply(G.Quants, function(xx) xx[commonAddresses,
                                                               "SD"]))
    colnames(GreenSD) <- names(G.Quants)
    RedSD <- do.call(cbind, lapply(R.Quants, function(xx) xx[commonAddresses,
                                                             "SD"]))
    colnames(RedSD) <- names(R.Quants)
    NBeads <- do.call(cbind, lapply(G.Quants, function(xx) xx[commonAddresses,
                                                              "NBeads"]))
    colnames(NBeads) <- names(G.Quants)
  }
  ptime2 <- proc.time()
  stime <- (ptime2 - ptime1)[3]
  if (verbose) {
    message(sprintf("done in %.1f seconds", stime))
  }
  if (verbose) {
    message("[read.metharray] Instantiating final object ... ",
            appendLF = FALSE)
  }
  ptime1 <- proc.time()
  if (extended) {
    out <- minfi::RGChannelSetExtended(Red = RedMean, Green = GreenMean,
                          RedSD = RedSD, GreenSD = GreenSD, NBeads = NBeads)
  }
  else {
    out <- minfi::RGChannelSet(Red = RedMean, Green = GreenMean)
  }
  rownames(out) <- commonAddresses
  out@annotation <- c(array = arrayTypes[1, 1], annotation = arrayTypes[1,2])
  ptime2 <- proc.time()
  stime <- (ptime2 - ptime1)[3]
  if (verbose) {
    message(sprintf("done in %.1f seconds", stime))
  }
  return(out)
}

#' @title readSheet.mm
#' @description MethylMaster version of readSheet() from Minfi
#' @param file The file parameter
#' @param recursive The recursive parameter
#' @return Formatted csv as dataframe
#' @export
readSheet.mm <- function(file,
                      recursive=TRUE
                      ){
  dataheader <- grep("^\\[DATA\\]", readLines(file),
                     ignore.case = TRUE)
  if (length(dataheader) == 0)
    dataheader <- 0
  df <- read.csv(file, stringsAsFactor = FALSE, skip = dataheader)
  nam <- grep(pattern = "Sentrix_Position", x = names(df),
              ignore.case = TRUE, value = TRUE)
  if (length(nam) == 1) {
    df$Array <- as.character(df[, nam])
    df[, nam] <- NULL
  }
  nam <- grep(pattern = "Array[\\._]ID", x = names(df),
              ignore.case = TRUE, value = TRUE)
  if (length(nam) == 1) {
    df$Array <- as.character(df[, nam])
    df[, nam] <- NULL
  }
  if (!"Array" %in% names(df)) {
    warning(sprintf("Could not infer array name for file: %s",
                    file))
  }
  nam <- grep("Sentrix_ID", names(df), ignore.case = TRUE,
              value = TRUE)
  if (length(nam) == 1) {
    df$Slide <- as.character(df[, nam])
    df[, nam] <- NULL
  }
  nam <- grep(pattern = "Slide[\\._]ID", x = names(df),
              ignore.case = TRUE, value = TRUE)
  if (length(nam) == 1) {
    df$Slide <- as.character(df[, nam])
    df[, nam] <- NULL
  }
  if (!"Slide" %in% names(df)) {
    warning(sprintf("Could not infer slide name for file: %s",
                    file))
  }
  else {
    df[, "Slide"] <- as.character(df[, "Slide"])
  }
  nam <- grep(pattern = "Plate[\\._]ID", x = names(df),
              ignore.case = TRUE, value = TRUE)
  if (length(nam) == 1) {
    df$Plate <- as.character(df[, nam])
    df[, nam] <- NULL
  }
  for (nam in c("Pool_ID", "Sample_Plate",
                "Sample_Well")) {
    if (nam %in% names(df)) {
      df[[nam]] <- as.character(df[[nam]])
    }
  }
  if (!is.null(df$Array)) {
    patterns <- sprintf("%s_%s_Grn.idat", df$Slide,
                        df$Array)
    allfiles <- list.files(##path = dirname(file),
                           path = unique(dirname(df$Basename)),
                           recursive = recursive,
                           full.names = TRUE)
    basenames <- sapply(X = patterns, FUN = function(xx) grep(xx,
                              allfiles, value = TRUE))
    names(basenames) <- NULL
    basenames <- sub("_Grn\\.idat.*", "",
                     basenames,
                     ignore.case = TRUE)
    df$Basename <- basenames
  }
  df
}

#' @title read.metharray.sheet.mm
#' @description MethylMaster version of read.metharray.sheet() from minfi
#' @param base The base parameter
#' @param pattern The pattern parameter
#' @param ignore.case The ignore.case parameter
#' @param recursive The recursive parameter
#' @param single.file The single.file parameter
#' @param single.file.path The single.file.path parameter
#' @param verbose The verbose parameter
#' @param ... Additional arguments passed to read.metharray.mm
#' @return A formatted csv as a dataframe
#' @export
read.metharray.sheet.mm <- function(base,
                                    pattern = "csv$",
                                    ignore.case = TRUE,
                                    recursive = TRUE,
                                    single.file = TRUE,
                                    single.file.path = NULL,
                                    verbose = TRUE,
                                    ...)
{
  if(single.file==TRUE){
    if(is.null (single.file.path) | !file.exists(single.file.path))
      stop("'single.file.path' does not exist or is set to NULL")
    csvfiles <- single.file.path
  }else if (!all(file.exists(base)))
    stop("'base' does not exists")
  info <- file.info(base)
  if (!all(info$isdir) && !all(!info$isdir)) {
    stop("'base needs to be either directories or files")
  }
  if(!exists("csvfiles")){
    if (all(info$isdir)) {
      csvfiles <- list.files(path = base,
                           recursive = recursive,
                           pattern = pattern,
                           ignore.case = ignore.case,
                           full.names = TRUE)
      if (verbose) {
        message("[read.metharray.sheet] Found the following CSV files:")
        print(csvfiles)
      }
    }else {
      csvfiles <- list.files(base, full.names = TRUE)
    }
  }
  dfs <- lapply(csvfiles, readSheet.mm)
  namesUnion <- Reduce(union, lapply(dfs, names))
  df <- do.call(what = rbind, args = lapply(dfs, function(df) {
    newnames <- setdiff(namesUnion, names(df))
    newdf <- matrix(data = NA, ncol = length(newnames), nrow = nrow(df),
                    dimnames = list(NULL, newnames))
    cbind(df, as.data.frame(newdf))
  }))
  return(df)
}

#' @title epicopy.mm
#' @description MethylMaster version of the main Epicopy::epicopy() function
#' from Epicopy
#' @param target_dir The input directory
#' @param output_dir The output directory
#' @param single.file Specify a single file (logical)
#' @param single.file.path The single file path
#' @param epic.manifest.path The path to the epic.manifest file
#' @param comparison The MehtylMaster comparison vector
#' @param reference_group The reference group to specify
#' @param project_name The project name
#' @param Normals The specified Normals see documentaion
#' @param sampNames Customize the sample names here
#' @param QN The QN parameter to use
#' @param Ref The Ref parameter to use
#' @param mode.bw The mode.bw parameter to use
#' @param mode.method The mode.method parameter to use
#' @param normal.cnv The normal.cnv parameter to use
#' @param mean.center Whether to use mean centering (default=TRUE)
#' @param filterProbes The filterProbes parameter
#' @param retainedProbes The retained Probes parameter
#' @param ncores The number of cores to use, note that multiple cores may not
#' work on every system
#' @param seg.undo.splits The seg.undo.splits parameter
#' @param seg.undo.SD The seg.undo.SD parameter
#' @param filterbycount The filterbycount parameter
#' @param min_probes The min_probes parameter
#' @param verbose The verbose parameter
#' @param keep_fnobj The keep_fnobj parameter
#' @param fn_output The fn_output parameter
#' @param ... Additional parameters to pass to epicopy.mm()
#' @import parallel
#' @import foreach
#' @import iterators
#' @import doParallel
#' @import minfi
#' @return Return a cna list containing CNV frame(s)
#' @export
epicopy.mm <- function(target_dir,
                       output_dir = NULL,
                       single.file = TRUE,
                       single.file.path = NULL,
                       epic.manifest.path = getwd(),
                       comparison = NULL,
                       reference_group = "internal", ##comparison
                       project_name = NULL,
                       Normals = NULL,
                       sampNames = NULL,
                       QN = FALSE,
                       Ref = "mode",
                       mode.bw = 0.1,
                       mode.method = "naive",
                       normal.cnv = NULL,
                       mean.center = TRUE,
                       filterProbes = FALSE,
                       retainedProbes = NULL,
                       ncores = 1L,
                       seg.undo.splits = "sdundo",
                       seg.undo.SD = 2,
                       filterbycount = TRUE,
                       min_probes = 50,
                       verbose = TRUE,
                       keep_fnobj = FALSE,
                       fn_output = NULL,
                       ...)
{
  if (is.null(output_dir)) {
    cat("Using current directory as output directory.\n")
    output_dir <- "."
  }
  target_sheet <- read.metharray.sheet.mm(target_dir,
                                          single.file = single.file,
                                          single.file.path = single.file.path)

  if(!is.null(split.by)){
    split.groups <- unique(target_sheet[[split.by]])
  }

  #########################################################################

  if(reference_group=="internal"){
    target_sheet <- target_sheet[target_sheet$Sample_Group %in% comparison[1],]
  }else{
    target_sheet <- target_sheet[target_sheet$Sample_Group %in% comparison,]
  }

  if (!is.null(Normals) & !is.na(Normals)){
    if(!Normals %in% c("breast", "lung", "thyroid",
                        "all")){
      if(!inherits(Normals, "character")){
        stop(paste0("Normals input not recognized.
                    Please review and see ?epicopy.\n"))
      }
      cat("Looking for normal annotation under",
          Normals, "column in samplesheet.\n")
      normal_index <- tolower(target_sheet[, Normals])
      ##if (!any(normal_index %in% "normal"))
      if(!any(normal_index %in% comparison[2])){
        stop("No normal sample found in annotation.\n")
      }
      ##nnormals <- table(normal_index %in% "normal")["TRUE"]
      nnormals <- table(normal_index %in% comparison[2])["TRUE"]
      cat("Found", nnormals, "normal samples.\n")
      ##Normals <- normal_index %in% "normal"
      Normals <- normal_index %in% comparison[2]
      if(is.null(normal.cnv)){
        normal.cnv <- TRUE
      }
    }
    else{
      if(is.null(normal.cnv)){
        normal.cnv <- FALSE
      }
    }
  }
  else{
    if(is.null(normal.cnv)){
      normal.cnv <- FALSE
    }
  }

  platforms <- unique(target_sheet$Platform)

  if(length(unique(platforms))>1){

  treatment.platform <- unique(target_sheet[target_sheet$Sample_Group %in%
                                                comparison[1],"Platform"])

  control.platform   <- unique(target_sheet[target_sheet$Sample_Group %in%
                                                comparison[2],"Platform"])

  rgsets <- list()
  foreach(i=1:length(comparison),.combine=rbind, .packages=c("foreach",
                                                              "minfi")) %do% {
    sub.target.sheet <- target_sheet[target_sheet$Sample_Group==comparison[i],]
    rgset <-
      read.metharray.exp(targets = sub.target.sheet,
                           verbose = verbose)
    rgsets[i] <- rgset
  }
    ##"IlluminaHumanMethylation450k
    ##"IlluminaHumanMethylationEPIC
    rgset <- minfi::combineArrays(rgsets[[1]],
                                  rgsets[[2]],
                            outType = ifelse(control.platform=="EPIC",
                                      "IlluminaHumanMethylationEPIC",
                                      "IlluminaHumanMethylation450k"),
                                  verbose=TRUE)

    lrr <- getLRR.mm(rgSet = rgset[,rgset$Sample_Group %in% comparison[1]],
                     ifelse(is.na(Normals),
                            Normals = NA,
                            Normals = rgset[,rgset$Sample_Group %in%
                                              comparison[2]]),
                     sampNames = sampNames,
                     QN = QN,
                     Ref = Ref,
                     mode.bw = mode.bw,
                     mode.method = mode.method,
                     normal.cnv = normal.cnv,
                     mean.center = mean.center,
                     filterProbes = filterProbes,
                     ...
    )

    cna <- LRRtoCNA.mm(lrr,
                       ncores = ncores,
                       seg.undo.splits = seg.undo.splits,
                       seg.undo.SD = seg.undo.SD,
                       epic.manifest.path = epic.manifest.path,
                       platform = control.platform
    )

    ##return(cna)

  }else{

  control.platform <- unique(target_sheet$Platform)

  rgset <- read.metharray.exp(targets = target_sheet,
                              verbose = verbose)

  lrr <- getLRR.mm(rgSet = rgset[,rgset$Sample_Group %in% comparison[1]],
                   Normals = ifelse(is.na(Normals),
                                    NA,
                                    rgset[,rgset$Sample_Group %in%
                                            comparison[2]]),
                sampNames = sampNames,
                QN = QN,
                Ref = Ref,
                mode.bw = mode.bw,
                mode.method = mode.method,
                normal.cnv = normal.cnv,
                mean.center = mean.center,
                filterProbes = filterProbes,
                ...
                )

  cna <- LRRtoCNA.mm(lrr,
                     ncores = ncores,
                     seg.undo.splits = seg.undo.splits,
                     seg.undo.SD = seg.undo.SD,
                     epic.manifest.path = epic.manifest.path,
                     platform = control.platform)

  ##return(cna)

  }

  cnas.out <- list(cna)


  return(cnas.out)

}
