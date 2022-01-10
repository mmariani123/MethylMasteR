#usr/bin/env Rscript

#' @title methyl_master_qreduceassay
#' @description MethylMaster version of the qreduceAssay function
#' Originally from the RaggedExperiment Package by Martin Morgan and
#' Marcel Ramos 2021.
#' @param x The x parameter
#' @param query The query parameter
#' @param simplifyReduce The simplifyReduce parameter
#' @param i The i parameter
#' @param withDimnames The withDimnames parameter
#' @param background the background parameter
#' @return A reduced ReaggedExperiment assay object
#' @export
methyl_master_qreduceassay <- function(x,
                                       query,
                                       simplifyReduce,
                                       i = 1,
                                       withDimnames = TRUE,
                                       background = NA_integer_){
  if (missing(i) && ncol(.mcols(x)) == 0)
    return(matrix(NA, 0, 0))
  stopifnot_simplify_ok(simplifyReduce, 3L)
  i <- .assay_i(x, i)
  mcol <- .mcols(x)[[i]][.rowidx(x)]
  dim <- .dim(x)
  subject <- unname(rowRanges(x))
  query <- granges(query)
  olap <- findOverlaps(query, subject)
  sidx <- subjectHits(olap)
  row <- queryHits(olap)
  col <- rep(seq_len(dim[[2]]), lengths(.assays(x)))[.rowidx(x)][sidx]
  score <- mcol[sidx]
  subject <- subject[sidx]
  qranges <- query[row]
  ranges <- restrict(subject, start = pmax(start(qranges),
                                           start(subject)), end = pmin(end(qranges), end(subject)))
  group <- (row - 1L) * max(col, 0) + col
  group <- match(group, unique(group))
  ugroup <- !duplicated(group)
  result <- simplifyReduce(unname(splitAsList(score, group)),
                           unname(splitAsList(ranges, group)), unname(qranges)[ugroup])
  group <- ugroup
  na <- as(background, class(result))
  dimnames <- list(NULL, NULL)
  if (withDimnames)
    dimnames <- list(as.character(query), .dimnames(x)[[2]])
  m <- matrix(na, nrow = length(query), ncol = dim[[2]], dimnames = dimnames)
  idx <- cbind(row = row[group], col = col[group])
  m[idx] <- result
  m[, .colidx(x), drop = FALSE]
}
