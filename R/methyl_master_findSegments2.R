#!/usr/bin/env Rscript

#' Modified findSegments function
#'
#' Find segements
#'
#' @param data
#' @param ctrl
#' @param ctrlAll
#' @param statistic
#' @param plot
#' @param delta
#' @param output
#' @param ylim
#' @param arrayType
#' @import plyr
#' @return
#' @export
findSegments2 <-
  function(data,
           ctrl,
           ctrlAll,
           statistic = "wilcoxon",
           plot = FALSE,
           delta = 500,
           output = "diff",
           ylim = NULL,
           arrayType="auto"){
    ##check if fast version can be used


    print("### Find Segments in CN Data ...")

    ##get annotation
    ##cnAnalysis450k::getAnnoData()
    ##can be "450k" or "EPIC"
    ##cnAnalysis450k::determineArrayType()
    ##Tries to determine the array type (EPIC, 450k)

    if (arrayType=="auto") {
      anno <- cnAnalysis450k::getAnnoData(
        cnAnalysis450k::determineArrayType(data))
    } else {
      anno <- cnAnalysis450k::getAnnoData(arrayType)
    }
    annoSorted <- anno[order(anno$chr, anno$pos), ]
    annoSorted <- annoSorted[which(rownames(annoSorted) %in%
                                     rownames(ctrlAll)),]

    # Get cg Borders of Chromosomes
    chrs <- paste("chr", c(1:22, "X", "Y"), sep = "")
    chBorder <- NULL
    for (ch in chrs) {
      subCh <- data.frame(annoSorted[which(annoSorted$chr == ch), ])
      subCh <- subCh[order(subCh$pos), ]
      vec <-
        data.frame(
          chr = ch,
          start = rownames(subCh)[1],
          end = rownames(subCh)[length(subCh[, 1])]
        )
      chBorder <- rbind.fill(chBorder, vec)
    }
    rownames(chBorder) <- chBorder[, 1]

    ##Find segments
    ct <- ctrl[match(rownames(annoSorted), names(ctrl))]
    ctAll <- ctrlAll[match(rownames(annoSorted), rownames(ctrlAll)),
                     ,drop=FALSE]
    smp <- data[match(rownames(annoSorted), rownames(data)),
                ,drop=FALSE]
    patData <- NULL

    for (j in 1:length(data[1, ])) {
      #different patients
      da <- smp[, j]
      smpName <- colnames(data)[j]
      print(paste("Processing", smpName, "..."))

      sampleSegments <- NULL
      ##splitpos
      for (ch in levels(factor(chBorder$chr))) {
        from <- which(names(da) == chBorder[ch, "start"])
        to <- which(names(da) == chBorder[ch, "end"])

        ##ratio or difference?
        if (output == "ratio") {
          rat <- da[from:to] / ct[from:to]
        } else if (output == "diff") {
          rat <- da[from:to] - ct[from:to]
        }

        namesRat <- names(da)[from:to]
        sel <- !is.na(rat) & !is.infinite(rat)
        rat <- rat[sel]
        namesRat <- namesRat[sel]

        ##find changepoints
        res <- changepoint::cpt.var(rat, method = "BinSeg")

        ##calculate segment data
        p.val <- c()
        if (length(res@cpts) > 1) {
          ends <- res@cpts
          starts <- c(0, ends[1:(length(ends) - 1)]) + 1
          median <- c()
          mean <- c()
          sd <- c()
          for (pos in 1:length(starts)) {
            median <- c(median,
                        median(rat[starts[pos]:ends[pos]]))
            mean <- c(mean,
                      mean(rat[starts[pos]:ends[pos]]))
            sd <- c(sd, sd(rat[starts[pos]:ends[pos]]))
            ##statistics
            if (statistic == "t.test") {
              #standard t-test
              p.val <- c(p.val,
                         t.test(da[starts[pos]:ends[pos]],
                                ctAll[starts[pos]:ends[pos]])$p.value)
            } else if (statistic == "wilcoxon") {
              #mann-whitney-wilcoxon test
              p.val <- c(p.val,
                         wilcox.test(da[starts[pos]:ends[pos]],
                                     ctAll[starts[pos]:ends[pos]])$p.value)
            } else {
              stop("Unknown statistic!")
            }
          }
          segDF <- NULL
        } else {
          starts <- 1
          ends <- res@cpts
          median <- median(rat[starts:ends])
          mean <- mean(rat[starts:ends])
          sd <- sd(rat[starts:ends])
          ##statistics
          if (statistic == "t.test") {
            #standard t-test
            p.val <-
              t.test(da[starts:ends],
                     ctAll[starts:ends])$p.value
          } else if (statistic == "wilcoxon") {
            #mann-whitney-wilcoxon test
            p.val <-
              wilcox.test(da[starts:ends],
                          ctAll[starts:ends])$p.value
          } else {
            stop("Unknown statistic!")
          }
        }

        ## plot segment borders
        if (plot) {
          changepoints <- res@cpts
          par(mfrow = c(1, 1))
          if (!is.null(ylim)) {
            plot(rat, main = ch, ylim = ylim)
          } else {
            plot(rat, main = ch)#, ylim=c(0.9,1.1))
          }
          abline(v = changepoints,
                 col = 2,
                 lwd = 2)
          par(mfrow = c(3, 3))
          for (q in 1:length(changepoints)) {
            cg <- changepoints[q]
            ##get borders
            cgMin <- ifelse(cg - delta < 0,
                            0,
                            ifelse(
                              q == 1,
                              cg - delta,
                              ifelse(
                                cg - delta < changepoints[q - 1] - delta,
                                changepoints[q - 1],
                                cg - delta
                              )
                            ))
            segMedian1 <- median(rat[(cgMin):cg])
            sd1 <- sd(rat[(cgMin):cg])

            cgMax <- ifelse(
              cg + delta > length(rat),
              length(rat),
              ifelse(
                q == length(changepoints),
                cg + delta,
                ifelse(
                  cg + delta > changepoints[q + 1],
                  changepoints[q + 1],
                  cg + delta
                )
              )
            )
            segMedian2 <- median(rat[cg:cgMax])
            sd2 <- sd(rat[cg:cgMax])

            ##plot
            if (!is.null(ylim)) {
              plot(
                rat,
                xlim = c(cgMin, cgMax),
                ylim = ylim,
                main = paste(
                  ch,
                  ": ",
                  round(segMedian1, 3),
                  ",",
                  round(sd1, 3),
                  " / ",
                  round(segMedian2, 3),
                  ",",
                  round(sd2, 3)
                )
              )
            } else {
              plot(
                rat,
                xlim = c(cgMin, cgMax),
                #, ylim=c(0.9, 1.1),
                main = paste(
                  ch,
                  ": ",
                  round(segMedian1, 3),
                  ",",
                  round(sd1, 3),
                  " / ",
                  round(segMedian2, 3),
                  ",",
                  round(sd2, 3)
                )
              )
            }
            abline(v = cg,
                   col = 4,
                   lwd = 3)
            abline(h = segMedian1,
                   col = 3,
                   lwd = 2)
            abline(h = segMedian2,
                   col = 2,
                   lwd = 2)
          }
          abline(v = changepoints, col = 2)
        }

        segDF <- data.frame(
          chr = ch,
          startCG = namesRat[starts],
          #start=starts,
          endCG = namesRat[ends],
                    #end=ends,
                    median = median,
                    mean = mean,
                    sd = sd,
                    smp = smpName,
                    p.val = p.val
                )

                sampleSegments <- rbind.fill(sampleSegments, segDF)
            }

            #plot(sampleSegments$mean, col=sampleSegments$chr)
            patData <- rbind.fill(patData, sampleSegments)
        }
        return (patData)
  }
