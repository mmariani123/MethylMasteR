#!/usr/bin/env Rscript

#' @title methyl_master_population_ranges
#' @description MwthylMaster version of the populatonRanges function
#' originally by Ludwig Geistlinger 2017
#' @param data The data parameter
#' @param grl The grl parameter
#' @param mode The mode parameter
#' @param density The density parameter
#' @param ro.thresh The ro.thresh parameter
#' @param multi.assign The multi.assign parameter
#' @param verbose The verbose parameter
#' @param min.size The min.size parameter
#' @param classify.ranges The classify.ranges parameter
#' @param type.thresh The type.thresh parameter
#' @param est.recur The est.recur parameter
#' @import CNVRanger
#' @import S4Vectors
#' @return A populationRanges CNV overlaps object
#' @export
methyl_master_population_ranges <- function(grl,
                                       mode=c("density", "RO"),
                                       density=0.1,
                                       ro.thresh=0.5,
                                       multi.assign=FALSE,
                                       verbose=FALSE,
                                       min.size=2,
                                       classify.ranges=TRUE,
                                       type.thresh=0.1,
                                       est.recur=FALSE)
{
  mode <- match.arg(mode)

  grl.check <- grl

  if(classify.ranges) grl.check <- .excludeNeutralRanges(grl.check)

  if(isEmpty(grl.check)){message(paste0("All copy-number neutral ",
                                  "regions are diploid, ",
                               "rerunning with classify.ranges==FALSE"))
  }else{
                          grl <- grl.check
                        }

  if(mode == "density") pranges <- .densityPopRanges(grl, density)
  else pranges <- .roPopRanges(grl, ro.thresh, multi.assign, verbose)

  pranges <- pranges[BiocGenerics::width(pranges) >= min.size]
  if(classify.ranges) pranges <- .classifyRegs.mm(pranges, unlist(grl),
                                                  type.thresh)
  if(est.recur) pranges <- .estimateRecurrence(pranges, grl)
  return(pranges)
}

.excludeNeutralRanges <- function(grl)
{
  # exclude neutral regions (CN = 2, diploid)
  calls <- IRanges::stack(grl, index.var = "sample")
  if(!("state" %in% colnames(S4Vectors::mcols(calls))))
  stop("Required column \'state\' storing integer copy number state not found")

  is.neutral <- calls$state == "2"
  sum.neutral <- sum(is.neutral)
  if(sum.neutral)
  {
    message(paste("Excluding", sum.neutral, "copy-number neutral regions",
                  "(CN state = 2, diploid)"))
    calls <- calls[!is.neutral]
    sample <- calls$sample
    ind <- colnames(S4Vectors::mcols(calls)) != "sample"
    S4Vectors::mcols(calls) <- S4Vectors::mcols(calls)[,ind, drop = FALSE]
    grl <- split(calls, sample)
  }
  grl
}

#' @title plotRecurrentRegions
#' @description MethylMasteR version of plotReccurentRegions function
#' "Plot recurrent CNV regions
#' Illustrates summarized CNV regions along a chromosome."
#' @param regs The regs Parameter
#' @param genome The genome parameter
#' @param chr The chr parameter
#' @param pthresh The pthresh parameter
#' @export
plotRecurrentRegions <- function(regs,
                                 genome,
                                 chr,
                                 pthresh=0.05){
  if (!requireNamespace("Gviz", quietly = TRUE))
    stop(paste("Required package \'Gviz\' not found.",
               "Use \'BiocManager::install(\"Gviz\") to install it."))

  colM <- "gray24"
  colB <- "gray94"
  tlevels <- c("gain", "loss", "both")

  chr.regs <- subset(regs, seqnames == chr)
  if(!length(chr.regs)) "No CNV regions to plot on the specified chromosome"

  itrack <- Gviz::IdeogramTrack(genome=genome, chr=chr, fontsize=15)
  gtrack <- Gviz::GenomeAxisTrack(littleTicks=TRUE, fontsize=15)

  gain.regs <- subset(chr.regs, type=="gain")
  tlist <- list()
  if(length(gain.regs))
  {
    gain.track <- Gviz::DataTrack(gain.regs, data=gain.regs$freq,
                                type="h", groups=factor("gain", levels=tlevels),
                                  name="#samples", cex.title=1, cex.axis=1,
                                  font.axis=2, col.title=colM, col.axis=colM,
                                  background.title=colB, legend=TRUE)
    tlist <- c(tlist, gain.track)
  }

  loss.regs <- subset(chr.regs, type=="loss")
  if(length(loss.regs))
  {
    loss.track <- Gviz::DataTrack(loss.regs, data=loss.regs$freq,
                                type="h", groups=factor("loss", levels=tlevels),
                                  name="#samples", cex.title=1, cex.axis=1,
                                  font.axis=2, col.title=colM, col.axis=colM,
                                  background.title=colB, legend=TRUE)
    tlist <- c(tlist, loss.track)
  }

  both.regs <- subset(chr.regs, type=="both")
  if(length(both.regs))
  {
    both.track <- Gviz::DataTrack(both.regs, data=both.regs$freq,
                                type="h", groups=factor("both", levels=tlevels),
                                  name="#samples", cex.title=1, cex.axis=1,
                                  font.axis=2, col.title=colM, col.axis=colM,
                                  background.title=colB, legend=TRUE)
    tlist <- c(tlist, both.track)
  }

  otrack <- Gviz::OverlayTrack(trackList = tlist, background.title=colB)
  ylim <- vapply(tlist, function(tr) range(S4Vectors::values(tr)), numeric(2))
  ylim <- extendrange(range(as.vector(ylim)))

  tracklist <- list(itrack, gtrack, otrack)

  # significant regions
  sig.regs <- subset(chr.regs, pvalue < pthresh)
  if(length(sig.regs))
  {
    atrack <- Gviz::AnnotationTrack(sig.regs, name="recur", cex.title=1,
                                    cex.axis=1, font.axis=2, col.title=colM,
                                    col.axis=colM, background.title=colB)
    tracklist <- c(tracklist, atrack)
  }

  Gviz::plotTracks(tracklist, ylim=ylim)
}

#' @title cnvOncoPrint
#' @description The MethylMaster version of the OncoPrint() function
#' to:
#' "plot for CNV regions
#' Illustrates overlaps between CNV calls and genomic features across a
#' sample population."
#' @param calls The calls parameter
#' @param features The features parameter
#' @param multi.calls The multi.calls parameter
#' @param top.features The top.features parameter
#' @param top.samples The top.samples parameter
#' @param ... Additional parameters passed to cnvOncoPrint
#' @export
cnvOncoPrint <- function(calls,
                         features,
                         multi.calls=.largest,
                         top.features = 25,
                         top.samples = 100,
                         ...){
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE))
    stop(paste("Required package \'ComplexHeatmap\' not found.",
               "Use \'BiocManager::install(\"ComplexHeatmap\") ",
               "to install it."))

  # summarize CNV calls in feature regions
  calls <- as(calls, "RaggedExperiment")
  qassay <- RaggedExperiment::qreduceAssay(calls,
                                           features,
                                           simplifyReduce=multi.calls,
                                           background=2)

  if("symbol" %in% colnames(S4Vectors::mcols(features))){
    rownames(qassay) <- features$symbol
  }else{
    rownames(qassay) <- names(features)
  }

  qassay <- qassay[rownames(qassay) != "",]
  qassay <- qassay[!duplicated(rownames(qassay)),]

  # select most frequently altered features
  indr <- order(rowSums(qassay != 2), decreasing=TRUE)
  has.top.features <- top.features > 0 && top.features < nrow(qassay)
  if(has.top.features) indr <- indr[seq_len(top.features)]
  qassay <- qassay[indr,]

  # select most frequently altered samples
  indc <- order(colSums(qassay != 2), decreasing=TRUE)
  has.top.samples <- top.samples > 0 && top.samples < ncol(qassay)
  if(has.top.samples) indc <- indc[seq_len(top.samples)]
  qassay <- qassay[,indc]

  # recode genotypes
  qassay[qassay == 2] <- " "
  qassay[qassay == "0"] <- "HOMDEL"
  qassay[qassay == "1"] <- "HETDEL"
  qassay[qassay == "3"] <- "GAIN1"
  qassay[qassay == "4"] <- "GAIN2+"

  # assign colors to genotypes
  cb.red <- "#D55E00"
  cb.blue <- "#0072B2"
  cb.lightblue <- "#56B4E9"
  cb.orange <- "#E69F00"

  colors = c( "HOMDEL" = cb.blue,
              "HETDEL" = cb.lightblue,
              "GAIN1" = cb.orange,
              "GAIN2+" = cb.red)

  # map functions to CNV types
  .mf <- function(couleur)
  {
    args <- alist(x =, y =, w =, h =)
    args <- as.pairlist(args)
    body <- substitute({
      grid::grid.rect(x, y,
                      w - grid::unit(0.5, "mm"),
                      h - grid::unit(0.5, "mm"),
                      gp = grid::gpar(fill = z, col = NA))
    }, list(z = couleur))
    eval(call("function", args, body))
  }
  mutfuns <- lapply(colors, .mf)

  background <- function(x, y, w, h)
    grid::grid.rect(x, y,
                    w - grid::unit(0.5, "mm"),
                    h - grid::unit(0.5, "mm"),
                    gp = grid::gpar(fill = "#FFFFFF", col = "#FFFFFF"))
  mutfuns2 <- c(background = background, mutfuns)

  # plot
  heatmap_legend_param <- list(title = "CNV type",
                               at = c("HOMDEL", "HETDEL", "GAIN1", "GAIN2"),
                               labels = c("2-copy loss",
                                          "1-copy loss",
                                          "1-copy gain",
                                          ">= 2-copy gain"))

  suppressMessages(
    ComplexHeatmap::oncoPrint(
      qassay, alter_fun = mutfuns2, col = colors, show_pct = FALSE,
      heatmap_legend_param = heatmap_legend_param, ...)
  )
}

#' @title .densityPopRanges
#' @description The MethylMaster version of .densityPopRanges() function
#' @param grl param grl
#' @param density param density
#' @import GenomicRanges
#' @importFrom S4Vectors subjectHits
#' @return ranges
#' @export
.densityPopRanges <- function(grl,
                              density)
{
  gr <- unlist(GenomicRanges::GRangesList(grl))
  cover <- GenomicRanges::reduce(gr)
  disjoint <- GenomicRanges::disjoin(gr)

  olaps <- GenomicRanges::findOverlaps(disjoint, cover)
  covered.hits <- S4Vectors::subjectHits(olaps)
  cover.support <- GenomicRanges::countOverlaps(cover, gr)[covered.hits]
  ppn <- GenomicRanges::countOverlaps(disjoint, gr) / cover.support

  pranges <- GenomicRanges::reduce(disjoint[ppn >= density])
  return(pranges)
}

#' @title .classifyRegs.mm
#' @description The MethylMaster version of .classifyRegs() function
#' @param regs param regs
#' @param calls param calls
#' @param type.thresh param type.thresh
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors queryHits
#' @importFrom S4Vectors subjectHits
#' @return regs
#' @export
.classifyRegs.mm <- function(regs,
                             calls,
                             type.thresh=0.1
                             ){
  olaps <- GenomicRanges::findOverlaps(regs, calls)
  qh <- S4Vectors::queryHits(olaps)
  sh <- S4Vectors::subjectHits(olaps)

  # number of samples
  sampL      <- split(names(calls)[sh], qh)
  .nrSamples <- function(x) length(unique(x))
  samples    <- vapply(sampL, .nrSamples, numeric(1), USE.NAMES=FALSE)
  regs$freq  <- samples

  # type: gain, loss, both
  stateL     <- split(calls$state[sh], qh)
  types      <- vapply(stateL, .getType.mm, character(1),
                  type.thresh=type.thresh, USE.NAMES=FALSE)
  regs$type  <- types

  return(regs)
}

#' @title .getType.mm
#' @description The MethylMaster version of .getType() function
#' @param states integer vector
#' @param state.thresh numeric
#' @return character
#' @export
.getType.mm <- function(states, type.thresh=0.1)
{
  ##if(any(states!=2)){
  fract.amp <- mean(states > 2)
  fract.del <- mean(states < 2)
  type <- c(fract.amp, fract.del) > type.thresh
  type <- c("gain", "loss")[type]
  if(length(type) > 1){
    type <- "both"
  }
  ##}else if(!any(states!=2)){
  ##  type <- "neutral"
  ##}
  return(type)
}

#' @title .roPopRanges
#' @description  MethylMaster version of the .roPopRanges function
#' "(2) RO approach"
#' @param grl The grl parameter
#' @param ro.tresh The ro.thresh parameter
#' @param multi.assign The multu.assign parameter
#' @param verbose The verbose parameter
#' @return GenomicRanges::GRanges
#' @export
.roPopRanges <- function(grl,
                         ro.thresh=0.5,
                         multi.assign=FALSE,
                         verbose=FALSE)
{
  gr <- unlist(grl)
  S4Vectors::mcols(gr) <- NULL

  # build initial clusters
  init.clusters <- GenomicRanges::reduce(gr)

  if(verbose) message(paste("TODO:", length(init.clusters)))

  # cluster within each initial cluster
  olaps <- GenomicRanges::findOverlaps(init.clusters, gr)
  qh <- S4Vectors::queryHits(olaps)
  sh <- S4Vectors::subjectHits(olaps)
  cl.per.iclust <- lapply(seq_along(init.clusters),
                          function(i)
                          {
                            if(verbose) message(i)
                            # get calls of cluster
                            ind <- sh[qh==i]
                            ccalls <- gr[ind]
                            if(length(ccalls) < 2) return(ccalls)
                            clusters <- .clusterCalls(ccalls,
                                                      ro.thresh,
                                                      multi.assign)
                            if(is.list(clusters))
                              clusters <- IRanges::extractList(ccalls, clusters)
                            clusters <- range(clusters)
                            if(is(clusters, "GRangesList")){
                              clusters <- sort(unlist(clusters))
                            }
                            return(clusters)
                          })
  ro.ranges <- unname(unlist(GenomicRanges::GRangesList(cl.per.iclust)))
  return(ro.ranges)
}

#' @title .clusterCalls
#' @description  MethylMaster version of the .clusterCalls function
#' "the clustering itself then goes sequentially through the identified RO hits,
#' touching each hit once, and checks whether this hit could be merged to
#' already existing clusters"
#' @param calls GenomicRanges::GRanges
#' @param ro.thresh numeric
#' @param multi.assign logical
#' @return list of integer vectors
#' @export
.clusterCalls <- function(calls, ro.thresh=0.5, multi.assign=FALSE)
{
  hits <- .getROHits(calls, ro.thresh)

  # exit here if not 2 or more hits
  if(length(hits) < 2) return(calls)

  # worst case: there are as many clusters as hits
  cid <- seq_along(hits)
  qh <- S4Vectors::queryHits(hits)
  sh <- S4Vectors::subjectHits(hits)

  # touch each hit once and check whether ...
  # ... it could be merged to a previous cluster
  for(i in 2:length(hits))
  {
    # has this hit already been merged?
    if(cid[i] != i) next
    curr.hit <- hits[i]

    # check each previous cluster
    for(j in seq_len(i-1))
    {
      # has this hit already been merged?
      if(cid[j] != j) next

      # if not, check it
      prev.cluster <- hits[cid == j]
      mergeIndex <- .getMergeIndex(curr.hit, prev.cluster, hits)

      if(!is.null(mergeIndex))
      {
        if(length(mergeIndex) == 1) cid[i] <- j
        else cid[mergeIndex] <- j
        break
      }
    }
  }

  # compile hit clusters
  hit.clusters <- unname(S4Vectors::splitAsList(hits, cid))

  # extract call clusters
  call.clusters <- lapply(hit.clusters,
    function(h) union(S4Vectors::queryHits(h), S4Vectors::subjectHits(h)))

  # can calls be assigned to more than one cluster?
  if(!multi.assign) call.clusters <- .pruneMultiAssign(call.clusters)

  return(call.clusters)
}

#' @title .getROHits
#' @desscription The MethylMaster version of the .getROHits function
#'' "given a set individual calls, returns overlaps (hits) between them
#' that satisfy the RO threshold"
#' @param calls GenomicRanges::GRanges
#' @param ro.thresh numeric
#' @return S4Vectors::Hits
.getROHits <- function(calls, ro.thresh=0.5)
{
  # calculate pairwise ro
  hits <- GenomicRanges::findOverlaps(calls,
                                      drop.self=TRUE,
                                      drop.redundant=TRUE)

  x <- calls[S4Vectors::queryHits(hits)]
  y <- calls[S4Vectors::subjectHits(hits)]
  pint <- GenomicRanges::pintersect(x, y)
  rovlp1 <- BiocGenerics::width(pint) / BiocGenerics::width(x)
  rovlp2 <- BiocGenerics::width(pint) / BiocGenerics::width(y)

  # keep only hits with ro > threshold
  ind <- rovlp1 > ro.thresh & rovlp2 > ro.thresh
  hits <- hits[ind]

  # exit here if not 2 or more hits
  if(length(hits) < 2) return(hits)

  rovlp1 <- rovlp1[ind]
  rovlp2 <- rovlp2[ind]

  # order hits by RO
  pmins <- pmin(rovlp1, rovlp2)
  ind <- order(pmins, decreasing=TRUE)

  qh <- S4Vectors::queryHits(hits)
  sh <- S4Vectors::subjectHits(hits)
  hits <- S4Vectors::Hits(qh[ind], sh[ind],
                          S4Vectors::queryLength(hits),
                          S4Vectors::subjectLength(hits))
  S4Vectors::mcols(hits)$RO1 <- rovlp1[ind]
  S4Vectors::mcols(hits)$RO2 <- rovlp2[ind]

  return(hits)
}

# ' @ttitle .getMergeIndex
#' @description the MethylMaster version of the .getMergeIndex function
#' "decides whether a given hit can be merged to an already existing cluster
#' mergeability requires that all cluster members satisfy the pairwise RO
#' threshold"
#' @param hit S4Vectors::Hits
#' @param cluster S4Vectors::Hits
#' @param hits S4Vectors::Hits
#' @return integer vector
#' @export
.getMergeIndex <- function(hit, cluster, hits)
{
  # (1) check whether query / S4Vectors::subject of hit is part of cluster
  curr.qh <- S4Vectors::queryHits(hit)
  curr.sh <- S4Vectors::subjectHits(hit)

  prev.qh <- S4Vectors::queryHits(cluster)
  prev.sh <- S4Vectors::subjectHits(cluster)
  prev.members <- union(prev.qh, prev.sh)
  is.part <- c(curr.qh, curr.sh) %in% prev.members

  # (2) can it be merged?
  mergeIndex <- NULL

  # (2a) query *and* subject of hit are part of cluster
  if(all(is.part)) mergeIndex <- 1

  # (2b) query *or* subject of hit are part of cluster
  else if(any(is.part))
  {
    # check whether the call which is not part of the cluster
    # has sufficient RO with all others in the cluster
    npart <- c(curr.qh, curr.sh)[!is.part]
    len <- length(prev.members)
    req.hits <- S4Vectors::Hits(   rep.int(npart, len),
                                   prev.members,
                                   S4Vectors::queryLength(hits),
                                   S4Vectors::subjectLength(hits))
    is.gr <- S4Vectors::queryHits(req.hits) > S4Vectors::subjectHits(req.hits)
    req.hits[is.gr] <- t(req.hits[is.gr])
    mergeIndex <- match(req.hits, hits)
    if(any(is.na(mergeIndex))) mergeIndex <- NULL
  }
  return(mergeIndex)
}

#' @title .pruneMultiAssign
#' @description the MethylMaster version of the .pruneMultiAssign function
# "as the outlined procedure can assign a call to multiple clusters
# (in the most basic case a call A that has sufficient RO with a call B and
# a call C, but B and C do not have sufficient RO), this allows to optionally
# strip away such multi-assignments"
#' @param clusters list of integer vectors
#' @return list of integer vectors
#' @export
.pruneMultiAssign <- function(clusters)
{
  cid <- seq_along(clusters)
  times <- lengths(clusters)
  cid <- rep.int(cid, times)
  ind <- unlist(clusters)
  ndup <- !duplicated(ind)
  ind <- ind[ndup]
  cid <- cid[ndup]

  pruned.clusters <- split(ind, cid)
  return(pruned.clusters)
}

.estimateRecurrence <- function(regs, calls, mode=c("approx", "perm"),
                                perm=1000, genome="hg19")
{
  mode <- match.arg(mode)
  calls <- unlist(calls)
  obs.scores <- lapply(seq_along(regs),
                       function(i) .scoreRegion(regs[i], calls))

  if(mode == "approx") pvals <- .approxP.mm(obs.scores, regs$type)
  else pvals <- .permP(obs.scores, regs, calls, perm, genome)
  regs$pvalue <- pvals
  return(regs)
}

.approxP.mm <- function(obs.scores, types)
{
  #In some cases there are no "both" there are only "loss" and "gain"
  #so need to modify function;
  both.scores <- obs.scores[types == "both"]
  both.scores <- unlist(both.scores)
  gain.scores <- obs.scores[types == "gain"]
  gain.scores <- unlist(both.scores)
  loss.scores <- obs.scores[types == "loss"]
  loss.scores <- unlist(both.scores)
  if(is.null(both.scores) & is.null(loss.scores)){
    gain.ind    <- seq(1, length(gain.scores))
    gain.scores <- obs.scores[types == "gain"]
    gain.scores <- c(unlist(gain.scores), unlist(obs.scores)[gain.ind])
    gain.ecdf <- ecdf(gain.scores)
    loss.ecdf <- ecdf(2)
  }else if(is.null(both.scores) & is.null(gain.scores)){
    loss.ind <- seq(1,length(obs.scores),by=1)
    loss.scores <- obs.scores[types == "loss"]
    loss.scores <- c(unlist(loss.scores), unlist(obs.scores)[loss.ind])
    loss.ecdf <- ecdf(loss.scores)
    gain.ecdf <- ecdf(2)
  }else{
    gain.ind <- seq(1, length(both.scores) - 1, by=2)
    loss.ind <- seq(2, length(both.scores), by=2)
    gain.scores <- obs.scores[types == "gain"]
    gain.scores <- c(unlist(gain.scores), both.scores[gain.ind])
    loss.scores <-  obs.scores[types == "loss"]
    loss.scores <- c(unlist(loss.scores), both.scores[loss.ind])
    gain.ecdf <- ecdf(gain.scores)
    loss.ecdf <- ecdf(loss.scores)
  }

  .getP <- function(i##,
                    ##gain.ecdf(),
                    ##loss.ecdf()
                    )
  {
    s <- obs.scores[[i]]
    ty <- types[i]
    if(ty == "both")
    {
      m <- which.max(s)
      p <- ifelse(m == 1, gain.ecdf(s), loss.ecdf(s))
    }
    else p <- ifelse(ty == "gain",
                     gain.ecdf(s),
                     loss.ecdf(s))
    return(1 - p)
  }

  pvals <- vapply(seq_along(obs.scores),
                  .getP,
                  ##gain.ecdf(),
                  ##loss.ecdf(),
                  FUN.VALUE=numeric(1)
                  )

  return(pvals)
}

.permP <- function(obs.scores, regs, calls, perm, genome)
{
  # genome specified?
  g <- GenomeInfoDb::genome(regs)
  ind <- !is.na(g)
  if(any(ind)) genome <- g[ind]
  rregs <- regioneR::randomizeRegions(calls, genome=genome, mask=NA,
                                      allow.overlaps=TRUE, per.chromosome=TRUE)
  # TODO:
  pvals <- vapply(seq_along(regs), function(i) i, numeric(1))
  return(pvals)
}

#' @title .scoreRegion
#' @description The MethylMasteR version of the .scoreRegion function
#' @param r GenomicRanges::GRanges
#' @param calls GenomicRanges::GRanges
#' @param fract.thresh numeric
#' @return numeric
.scoreRegion <- function(r, calls, fract.thresh=0.1)
{
  chr <- GenomicRanges::seqnames(r)
  chr <- as.character(chr)

  ind <- GenomicRanges::seqnames(calls) == chr
  rcalls <- GenomicRanges::restrict(calls[ind],
                                    start=GenomicRanges::start(r),
                                    end=GenomicRanges::end(r))

  type <- .getType.mm(rcalls$state, fract.thresh)
  ##if(type == "neutral"){
    ##score=.getFreq.mm(rcalls, type)
    ##}else if(type != "both" & type != "neutral") score <- .getFreq.mm(rcalls,
  ##                                                                    type)
  if(type != "both") score <- .getFreq.mm(rcalls, type)
  else score <- vapply(c("gain", "loss"),
                       function(type) .getFreq.mm(rcalls, type), numeric(1))
  return(score)
}

#' @title .getfreq.mm
#' @description MethylMaster version of the .getFreq() function
#' @param calls GenomicRanges::GRanges
#' @param type character
#' @return numeric
#' @export
##.getFreq.mm <- function(calls, type=c("gain", "loss", "neutral"))
.getFreq.mm <- function(calls, type=c("gain", "loss"))
{
  type <- match.arg(type)
  is.amp <- type == "gain"

  chr <- GenomicRanges::seqnames(calls)[1]
  chr <- as.character(chr)

  ##if(type=="neutral"){
  ##  ind <- calls$state == 2
  ##  scalls <- calls[ind]
  ##}else{
    ind <- if(is.amp) calls$state > 2 else calls$state < 2
    scalls <- calls[ind]
  ##}

  ##if(type=="neutral"){
  ##  w <- scalls$state
  ##  ##w <- 0
  ##}else{
    w <- if(is.amp) scalls$state - 2 else 2 - scalls$state
  ##}

  x <- GenomicRanges::coverage(scalls, weight=w)[[chr]]
  cs <- cumsum(S4Vectors::runLength(x))
  len <- length(cs)
  if(len == 1) return(0)

  grid <- seq_len(len - 1)
  starts <-  cs[grid] + 1
  ends <- cs[grid + 1]
  cov <- S4Vectors::runValue(x)[grid + 1]

  max.score <- max(cov)
  return(max.score)
}
