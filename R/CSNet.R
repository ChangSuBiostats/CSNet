#' CSNet main function
#'
#' This function implements CSNet (estimation of Cell-type-Specific gene 
#' co-expression Networks using bulk gene expression data) from the paper
#' DOI:10.1080/01621459.2023.2297467
#'
#' @param P matrix of cell type proportions. n samples by K cell types.
#' @param X matrix of bulk gene expression. n samples by p genes.
#' @param ls_methods least squares method to estimate variance and covariance. 
#'  Default to list(var='nnls', covar = 'wls'), i.e. non-negative least squares
#'  for estimating variance and weighted least squares for estimating covariance.
#' @param weight_th the smallest value for weights in wls. Default to NULL.
#' @param plot_cv whether to plot error curves in cross validation.
#' @return A list of K cell-type-specific co-expression matrices.
#' @export
CSNet <- function(X, P,
                  ls_methods = list(var='nnls', covar = 'wls'),
                  weight_th = 0.01,
                  plot_cv = F){
  # run moment based regressions to obtain
  # dense CSNet estimator (prior to thresholding)
  dCSNet_cov_est <- mom_ls(P, X, ls_methods, weight_th)
  dCSNet_est <- lapply(dCSNet_cov_est, function(x) get_cor_from_cov(x))
  # run cross validation
  cv_res <- cross_validation(list(X = X, P = P),
                             list(methods_settings = ls_methods, th = weight_th),
                             n_splits = 1,
                             to_plot = plot_cv) 
  # apply thresholding and obtain CSNet estimator
  CSNet_est <- lapply(1:length(dCSNet_est), function(i){
    generalized_th(dCSNet_est[[i]], cv_res$th[i], 'SCAD', F) })
  return(CSNet_est)
}