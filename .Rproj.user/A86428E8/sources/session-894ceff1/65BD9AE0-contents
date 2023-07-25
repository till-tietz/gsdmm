#' fit Dirichlet Mixture Model
#' 
#' @param texts a List of character vectors 
#' @param n_iter integer number of iterations to run gibbs sampler for 
#' @param n_clust integer upper bound on number of clusters. The returned number of cluster will be smaller or equal to n_clust
#' @param alpha double governing the probability of assigning a text to a currently empty cluster (larger alpha means higher probability).
#' @param beta double governing the tradeoff between cluster size and fit. Smaller betas make clustering more sensitive to congruence between cluster-word and document-word distributions while larger betas make it more sensitive to cluster size.
#' @return an integer vector of clusters


gsdmm <- function(texts, n_iter, n_clust, alpha, beta) {
  
  vocab <- unique(unlist(texts)) # vocabulary
  vocab <- stats::setNames(1:length(vocab), vocab) # add integer indices to vocabulary
  
  V <- length(vocab) # num words in vocabulary
  d <- lapply(texts, function(i) match(i, names(vocab))) # turn texts into sequences of integers
  d <- lapply(d, function(i) i - 1) # index from 0
  
  clusters <- gsdmm_gibbs(
    d = d,
    I = n_iter,
    K = n_clust,
    alpha = alpha,
    beta = beta,
    V = V
  )
  
  return(clusters)
}