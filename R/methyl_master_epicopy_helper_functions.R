#!/usr/bin/env Rscripts

#' This will fix the issue with LRR2
#'
#'
#' @param pdat Path to the input file
#' @return compatable pData data frame
#' @export
.coerce.pData.mm <- function(pdat){
  pDat <- apply(pdat, 2, as.character)
  dimnames(pDat) <- dimnames(pdat)
  pDat <- as(pDat,
             "DataFrame",
             ##stringsAsFactors = FALSE ##Will cause error
  )
  return(pDat)
}

#' Epciopy modified convert log r ratio data to cna calls:
#'
#' Same as LucasLRRtoCNA2 function
#'
#' @param LRR
#' @param sampID
#' @param ncores
#' @param seg.undo.splits
#' @param seg.undo.SD
#' @param platform
#' @import minfi
#' @return returns segments
#' @export
LRRtoCNA.mm <- function(LRR,
                      sampID = NULL,
                      ncores = 1L,
                      seg.undo.splits = "sdundo",
                      seg.undo.SD = 2,
                      platform=c("hm450","EPIC"),
                      ...) {
  if(platform=="EPIC"){
    load(file=paste0(work.dir,
                     file.sep,
                     "EPIC.manifest.hg38.rda"))
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

#' Epciopy modified .funnorm function
#'
#'
#' @param rgSet
#' @param nPCs
#' @param sex
#' @param verbose
#' @import minfi
#' @import BiocGenerics
#' @return return final gm set
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

#' Epciopy modified getLRR function
#'
#'
#' @param rgSet
#' @param Normals
#' @param sampNames
#' @param QN
#' @param Ref
#' @param mode.bw
#' @param mode.method
#' @param normal.cnv
#' @param mean.center
#' @param filterProbes
#' @param ...
#' @import stats
#' @import minfi
#' @return return final LRR values
#' @export
getLRR.mm <- function(rgSet = rgset,
                      Normals = Normals,
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
debug: {
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
    rgset <- rgSet
    Normals <- Normals
  }
  if (inherits(Normals, "RGChannelSet")) {
    if (nrow(rgSet) != nrow(Normals)){
      ##  stop("Features in normal and tumor set do not match.")
      multi.platform <- TRUE
      combined.colnames <- c(colnames(rgset),colnames(Normals))
      norm.rgset <- Normals
      Normals <- combined.colnames  %in% colnames(Normals)
    }else{
      rgset <- BiocGenerics::combine(rgSet, Normals)
      Normals <- colnames(rgset) %in% colnames(Normals)
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
      if (!"all.normals" %in% ls(globalenv()))
        data("allNormals")
      Normals <- all.normals
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
      colnames(rgset) <- pData(rgset)$Sample_Name
      colnames(norm.rgset) <- pData(norm.rgset)$Sample_Name
    }
    else {
      if (ncol(rgset) == length(sampNames))
        stop("sampNames is not equal to number of samples in rgset.\n")
      colnames(rgset) <- sampNames
    }
    cat("Performing functional normalization.\n")
    fn.set <- .funnorm.mm(rgSet=rgset)
    fn.intensity <- minfi::getUnmeth(fn.set) + minfi::getMeth(fn.set)
    fn.intensity[fn.intensity <= 0] <- 0.01
    fn.intensity <- log2(fn.intensity)

    fn.set.norm <- .funnorm.mm(rgSet=norm.rgset)
    fn.intensity.norm <- minfi::getUnmeth(fn.set.norm) +
      minfi::getMeth(fn.set.norm)
    fn.intensity.norm[fn.intensity.norm <= 0] <- 0.01
    fn.intensity.norm <- log2(fn.intensity.norm)
    normal.fn <- fn.intensity.norm[, Normals]
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
      if (!setequal(autosomal, rownames(normal.fn)))
        stop("Normalization returns discordant number of probes.")
    }
    cat("Filtering for autosomal probes.\n")
    normal.fn <- normal.fn[autosomal, ]
    tumor.fn <- fn.intensity[autosomal, !Normals]
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
    cat("Filtering for autosomal probes.\n")
    normal.fn <- normal.fn[autosomal, ]
    tumor.fn <- fn.intensity[autosomal, !Normals]
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

#' Modified read.metharray from minfi
#'
#'
#' @param basenames
#' @param extended
#' @param verbose
#' @param force
#' @param diff.array.type
#' @import DelayedArray
#' @import BiocParallel
#' @import illuminaio
#' @import ENmix
#' @import minfi
#' @return methylation array data
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
  }
  if (!all(file.exists(G.files))) {
    noexits <- sub("\\.gz", "", G.files[!file.exists(G.files)])
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

#' Modified read.metharray.sheet from minfi
#'
#'
#' @param base
#' @param pattern
#' @param ignore.case
#' @param recursive
#' @param single.file
#' @param single.file.path
#' @param verbose
#' @param ... Additional arguments
#' @return formatted csv as df
#' @export
read.metharray.sheet.mm <- function(base,
                                    pattern = "csv$",
                                    ignore.case = TRUE,
                                    recursive = TRUE,
                                    single.file = TRUE,
                                    single.file.path = NULL,
                                    verbose = TRUE)
{
  readSheet <- function(file) {
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
      allfiles <- list.files(path = dirname(file), recursive = recursive,
                             full.names = TRUE)
      basenames <- sapply(X = patterns, FUN = function(xx) grep(xx,
                                                                allfiles, value = TRUE))
      names(basenames) <- NULL
      basenames <- sub("_Grn\\.idat.*", "",
                       basenames, ignore.case = TRUE)
      df$Basename <- basenames
    }
    df
  }
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
  if (all(info$isdir)) {
    csvfiles <- list.files(path = base, recursive = recursive,
                           pattern = pattern, ignore.case = ignore.case, full.names = TRUE)
    if (verbose) {
      message("[read.metharray.sheet] Found the following CSV files:")
      print(csvfiles)
    }
  }else {
    csvfiles <- list.files(base, full.names = TRUE)
  }
  dfs <- lapply(csvfiles, readSheet)
  namesUnion <- Reduce(union, lapply(dfs, names))
  df <- do.call(what = rbind, args = lapply(dfs, function(df) {
    newnames <- setdiff(namesUnion, names(df))
    newdf <- matrix(data = NA, ncol = length(newnames), nrow = nrow(df),
                    dimnames = list(NULL, newnames))
    cbind(df, as.data.frame(newdf))
  }))
  return(df)
}

#' Modified epicopy function from Epicopy
#'
#' I updated the read.metharray.sheet function
#' and incorporated it here
#'
#' @param target_dir
#' @param output_dir
#' @param project_name
#' @param Normals
#' @param sampNames
#' @param QN
#' @param Ref
#' @param mode.bw
#' @param mode.method
#' @param normal.cnv
#' @param mean.center
#' @param filterProbes
#' @param retainedProbes
#' @param ncores
#' @param seg.undo.splits
#' @param seg.undo.SD
#' @param filterbycount
#' @param min_probes
#' @param verbose
#' @param ...
#' @return return a CNV frame
#' @export
epicopy.mm <- function(target_dir,
                    output_dir = NULL,
                    single.file = TRUE,
                    single.file.path = NULL,
                    comparisons = NULL,
                    reference_group = "normal",
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
                    ...)
{
  if (is.null(output_dir)) {
    cat("Using current directory as output directory.\n")
    output_dir <- "."
  }
  target_sheet <- read.metharray.sheet.mm(target_dir,
                                          single.file = single.file,
                                          single.file.path = single.file.path)
  target_sheet <- target_sheet[target_sheet$Sample_Group %in% comparisons,]
  if (!is.null(Normals)) {
    if (!Normals %in% c("breast", "lung", "thyroid",
                        "all")) {
      if (!inherits(Normals, "character"))
        stop("Normals input not recognized. Please review and see ?epicopy.\n")
      cat("Looking for normal annotation under",
          Normals, "column in samplesheet.\n")
      normal_index <- tolower(target_sheet[, Normals])
      if (!any(normal_index %in% "normal"))
        stop("No normal sample found in annotation.\n")
      nnormals <- table(normal_index %in% "normal")["TRUE"]
      cat("Found", nnormals, "normal samples.\n")
      Normals <- normal_index %in% "normal"
      if (is.null(normal.cnv))
        normal.cnv <- TRUE
    }
    else {
      if (is.null(normal.cnv))
        normal.cnv <- FALSE
    }
  }
  else {
    if (is.null(normal.cnv))
      normal.cnv <- FALSE
  }

  rgset <- list()
  platforms <- unique(target_sheet$Platform)
  if(length(platforms)>1){
    rgsets <- list()
    foreach(i=1:length(platforms)) %do% {
      sub.target.sheet <- target_sheet[target_sheet$Platform==platforms[i],]
      rgset <-
        read.metharray.exp(targets = sub.target.sheet,
                           verbose = verbose)
      rgsets[[i]] <- rgset
      if(sub.target.sheet$Sample_Group=="normal"){
        names(rgsets[[i]]) <- "normal"
      }else{
        names(rgsets[[i]]) <- sub.target.sheet$Sample_Group
      }

    }
    lrr <- getLRR.mm(rgSet = rgset[names(rgset)!="normal"],
                     Normals = rgset[names(rgset)=="normal"],
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
                       seg.undo.SD = seg.undo.SD
    )
    if (!is.null(output_dir)) {
      if (output_dir == FALSE) {
        return(cna)
      }
    }
    export_gistic(input = cna, output_dir = output_dir, filterbycount = filterbycount,
                  min_probes = min_probes, segfile_name = project_name,
                  markerfile_name = project_name)
    return(cna)
  }else{

  rgset <- read.metharray.exp(targets = target_sheet,
                              verbose = verbose)
  lrr <- getLRR.mm(rgSet = rgset,
                Normals = Normals,
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
                     seg.undo.SD = seg.undo.SD
                     )
  if (!is.null(output_dir)) {
    if (output_dir == FALSE) {
      return(cna)
    }
  }
  export_gistic(input = cna, output_dir = output_dir, filterbycount = filterbycount,
                min_probes = min_probes, segfile_name = project_name,
                markerfile_name = project_name)
  return(cna)
  }
}
