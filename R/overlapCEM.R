###### Main Function ######

#' @title overlapCEM - main function
#'
#' @description Identifies overlapping co-expressed gene modules
#'
#' @docType package
#'
#' @author Quang-Huy Nguyen
#'
#' @usage overlapCEM(data = NULL, clinical = NULL, ncomp = NULL,
#'        standardize = T, method = c("ICA-FDR", "ICA-Zscore", "IPCA-FDR"),
#'        cex.text = 0.7)
#'
#' @param data a data frame or matrix. \code{data} has its rows are samples and its columns are genes.
#'
#' @param clinical a data frame or matrix. Input data serve to perform Pearson's orrelations between
#' each identified module and each clinical features. It includes its rows are samples, and its columns
#' are clinical features of your choice.
#'
#' @param ncomp positive integer. The optimal number of principal components. It should be >= 2.
#'
#' @param standardize logical. If your \code{data} are not standardized, just feed \code{T} or \code{TRUE}
#' to this paramerter. Default value is \code{T}.
#'
#' @param method Post-processing methods. Allowed values are "ICA-FDR", "ICA-Zscore",
#' or "IPCA-FDR".
#'
#' @param cex.text numeric. Change the font size of texts in cells of the heatmap showing correlations between
#' each identified module and each clinical features. Default value is 0.7.
#'
#' @return NULL
#'
#' @examples cem=overlapCEM(data = exp, clinical = cli, ncomp = 18,
#' method = 'ICA-Zscore')
#'
#' @export

overlapCEM = function(data = NULL, clinical = NULL, ncomp = NULL, standardize = T,
                      method = c("ICA-FDR", "ICA-Zscore", "IPCA-FDR"), cex.text = 0.7){

  #Errors
  if(missing(data)){
    stop("Error: gene expression data are missing. \n")
  }

  if(missing(clinical)){
    stop("Error: clinical data are missing. \n")
  }

  if(missing(ncomp)){
    stop("Error: You do not define the optimal number of principal components. \n")
  }

  if(missing(method)){
    stop("Error: Which post-processing method you want to choose? \n")
  }

  #main function
  if(method == 'ICA-FDR'){
    return(ICAfdr(data = data, clinical = clinical, standardize = standardize,
                  ncomp = ncomp, cex.text = cex.text))
  }

  if(method == 'ICA-Zscore'){
    return(ICAzscore(data = data, clinical = clinical, standardize = standardize,
                     ncomp = ncomp, cex.text = cex.text))
  }

  if(method == 'IPCA-FDR'){
    return(IPCAfdr(data = data, clinical = clinical, standardize = standardize,
                   ncomp = ncomp, cex.text = cex.text))
  }
}



