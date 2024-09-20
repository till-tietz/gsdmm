#' @aliases gsdmm-package
#' @import Rcpp
#' @useDynLib gsdmm, .registration = TRUE
"_PACKAGE"

#' Fit Dirichlet Mixture Model
#'
#' @param texts a list of character vectors
#' @param n_iter integer number of iterations to run gibbs sampler for
#' @param n_clust integer upper bound on number of clusters.
#' The returned number of cluster will be smaller or equal to n_clust
#' @param alpha double governing the probability of assigning a text
#' to a currently empty cluster (larger alpha means higher probability).
#' @param beta double governing the tradeoff between cluster size and fit.
#' Smaller betas make clustering more sensitive to congruence
#' between cluster-word and document-word distributions
#' while larger betas make it more sensitive to cluster size.
#' @param progress logical indicating whether to print progress bar.
#' @return a list that contains an integer vector of clusters
#' and a word-cluster matrix.
#' @export
gsdmm <- function(texts,
                  n_iter = 30L,
                  n_clust = 8L,
                  alpha = 0.1,
                  beta = 0.1,
                  progress = TRUE) {
  tokens <- unlist(texts, use.names = FALSE)
  vocab <- unique(tokens)

  # turn texts into sequences of integers.
  # equivalent to `lapply(texts, \(i) match(i, vocab) - 1)`
  d <- vctrs::vec_chop(
    fastmatch::fmatch(tokens, vocab) - 1,
    sizes = vctrs::list_sizes(texts)
  )

  clusters <- gsdmm_gibbs(
    d = d,
    I = n_iter,
    K = n_clust,
    alpha = alpha,
    beta = beta,
    V = length(vocab), # num words in vocabulary
    progress = progress
  )

  # cast triplet data to a sparse matrix
  clusters$distribution <- with(
    clusters$distribution,
    Matrix::sparseMatrix(
      i = fastmatch::fmatch(doc, unique(doc)),
      j = fastmatch::fmatch(tok, unique(tok)),
      x = value,
      dimnames = list(unique(doc), unique(tok))
    )
  )
  rownames(clusters$distribution) <- vocab

  return(clusters)
}
