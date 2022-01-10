#usr/bin/env Rscript

#' @title methyl_master_assay
#' @description MethylMasteR Assay function
#' Originally from sparseAssay() from the
#' RaggedExperiment Package by Martin Morgan and Marcel Ramos 2021
#' I have modified my own version
#' here for use in methyl master
#' @param sesame.idat.files.dir The input dir with the .idat files
#' @param x A ragged experiment object
#' @param i Integer or character name of assay to be transformed
#' @param withDimnames Whether to include the dimnames on the retunred matrix
#' @param background A value for the returned matrix
#' @param sparse Whether to return a sparse matrix or not
#' @import RaggedExperiment
#' @import Matrix
#' @return A reduced ReaggedExperiment assay object
#' @export
methyl_master_assay <- function(x,
          i = 1,
          withDimnames = TRUE,
          background = NA_integer_,
          sparse = FALSE)
{
  i <- RaggedExperiment:::.assay_i(x, i)
  mcol <- RaggedExperiment:::.mcols(x)[[i]]
  dim <- RaggedExperiment:::.dim(x)
  if (withDimnames) {
    dimnames <- RaggedExperiment:::.dimnames(x)
    if (is.null(dimnames[[1]]))
      dimnames[[1]] <- as.character(RaggedExperiment:::.rowRanges(x))
  }
  else {
    dimnames <- NULL
  }

  ##idx <- cbind(row = seq_len(dim[[1L]]), col = rep(seq_len(dim[[2]]),
  ##                                    lengths(RaggedExperiment:::.assays(x))))
  rows.num <- dim[[1L]]/dim[[2L]]
  idx <- cbind(row = seq_len(rows.num),
               col = rep(1,times=rows.num))

  if (sparse) {
    M <- Matrix::sparseMatrix(i = idx[, 1], j = idx[, 2],
                              x = mcol, dims = dim)
    dimnames(M) <- dimnames
  }
  else {
    na <- as(background, class(mcol))
    M <- matrix(na, nrow = dim[[1]]/dim[[2]], ncol = 1)
    rownames(M) <- dimnames[[1]][1:(dim[[1]]/dim[[2]])]
    colnames(M) <- "State"
    M[idx] <- mcol[1:rows.num]
    ##M <- M[RaggedExperiment:::.rowidx(x),
    ##       RaggedExperiment:::.colidx(x),
    ##       drop = FALSE]
  }
  M
}
