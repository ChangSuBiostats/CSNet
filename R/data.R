#' Simulated bulk gene expression and cell type proportion data
#'
#' Simulated data to demonstrate the use of CSNet, following the simulation 
#' examples of DOI:10.1080/01621459.2023.2297467
#'
#' @format ## `who`
#' A list with three elements:
#' \describe{
#'   \item{X}{matrix of bulk gene expression, n samples by p genes}
#'   \item{P}{matrix of cell type proportions, n samples by K cell types}
#'   \item{ct_X}{list of two cell-type-specific gene expression matrices, 
#'   each with n samples by p genes}
#'   ...
#' }
#' @source <https://github.com/ChangSuBiostats/CSNet_analysis>
"sim_data_list"