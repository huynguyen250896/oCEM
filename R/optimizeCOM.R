#' @title pcaCharts
#'
#' @description Visualizes the screeplot showing the optimal number of principal components
#'
#' @docType package
#'
#' @author Quang-Huy Nguyen
#'
#' @usage function(x,y)
#'
#' @param x %explained variance and cumulative %explained variance of real data.
#'
#' @param y average %explained variance of random matrices.
# NOT EXPORTED

pcaCharts <- function(x,y) {

  #%explained variance and cumulative %explained variance of real data
  x.var <- x$d^2/sum(x$d^2) * 100
  x.pvar <- cumsum(x.var)

  #average %explained variance of random matrices
  x.var1 <- y

  #plot
  par(mar = c(5,5,2,5))
  par(xaxs='i',yaxs='i')
  plot(x.var, xlab="k",  ylab="%explained by each component", axes = F, xlim=c(1,length(x.var)), ylim=c(round(min(x.var),1), round(max(x.var),1)), type='l',col='blue')
  lines(x.var1,col="purple")
  x.ticks = seq(1,length(x.var))
  y.ticks = seq(round(min(x.var),1), round(max(x.var),1))
  axis(1, at = x.ticks, cex.axis=0.85)
  axis(2, at = y.ticks, labels = paste0(y.ticks, "%"), las = 1, mgp = c(3, 0.75, 0), cex.axis=0.75)

  #Find point on intersection of blue and purple curves
  # Find points where x1 is above x2.
  above <- x.var > x.var1
  # Points always intersect when above=TRUE, then FALSE or reverse
  intersect.points <- which(diff(above) != 0)
  # Find the slopes for the line segment.
  x1.slopes <- x.var[intersect.points+1] - x.var[intersect.points]
  x2.slopes <- x.var1[intersect.points+1] - x.var1[intersect.points]
  # Find the intersection for the segment.
  x.points <- intersect.points + ((x.var1[intersect.points] - x.var[intersect.points]) / (x1.slopes-x2.slopes))
  y.points <- x.var[intersect.points] + (x1.slopes*(x.points-intersect.points))
  abline(v=x.points, h=round(y.points,1), col='red',lwd=2, lty=2)

  #Add a curve for cumulative explained variance
  par(new = T)
  par(xaxs='i',yaxs='i')
  plot(x.pvar, xlab=NA, ylab= NA, axes = F, ylim=c(0,100), type='l', col='green')
  y.ticks1 = seq(0,100,10)
  axis(side = 4, at = y.ticks1, labels = paste0(y.ticks1, "%"), las = 1, cex.axis=0.75)
  mtext("Cumulative explained variance", side = 4, line = 3)
}

#' @title optimizeCOM
#'
#' @description Determines the optimal number of principal components
#'
#' @docType package
#'
#' @author Quang-Huy Nguyen
#'
#' @usage overlapCEM(data = NULL, P=1000,
#'        standardize = T, verbose = T)
#'
#' @param data a data frame or matrix. \code{data} has its rows are samples and its columns are genes.
#'
#' @param P positive integer. The number of permutations. Default value is 1000.
#'
#' @param standardize logical. If your \code{data} are not standardized, just feed \code{T{} or \code{TRUE}
#' to this paramerter. Default value is \code{T}.
#'
#' @param cores positive integer. The number of cores used to parallely perform the function. Default value is 1. 
#'
#' @param verbose Default value is \code{TRUE}. A logical specifying whether to print details of analysis processes.
#'
#' @return NULL
#'
#' @examples optimizeCOM(data = exp)
#'
#' @export

optimizeCOM=function(data = NULL, P=1000, standardize = T, cores = 1, verbose = T){
  #Errors
  if(missing(data)){
    stop("Error: gene expression data are missing. \n")
  }

  #system time
  now = Sys.time()

  # defined log function
  mlog <- if(!verbose) function(...){} else function(...){
    message(...)
    flush.console()
  }

  # Determine the optimal number of indepedent components
  #singular value decomposition
  if(standardize == T | standardize == TRUE){
      svd=svd(scale(data))
  } else{
      svd=svd(data)
  }
  
  #%explained variance and cumulative %explained variance of the original data
  x.var = svd$d^2/sum(svd$d^2) * 100
  x.pvar <- cumsum(x.var)

  #P random matrices
  ranMat = parallel::mclapply(1:P, mc.cores = cores, FUN = function(p){
      #generate P random matrices gained by independently permuting samples
      permutated_matrix <- apply(data, 1, sample) %>% t()

      cc = list()
      if(standardize == T | TRUE){
          cc[[p]] <- svd(scale(permutated_matrix))$d
      } else {
          cc[[p]] <- svd(permutated_matrix)$d
      }

      #%explained variance of a random matrix
      x.var1 = list()
      x.var1[[p]] <- cc[[p]]^2/sum(cc[[p]]^2) * 100
  }) %>% do.call(what=rbind) %>% t()

  average = rowMeans(ranMat)
  #find point on intersection of blue and purple curves
  # Find points where x1 is above x2.
  # Points always intersect when above=TRUE
  intersect.points = which(diff(x.var > average) != 0)

  # Find the slopes for the line segment.
  x1.slopes <- x.var[intersect.points+1] - x.var[intersect.points]
  x2.slopes <- ranMat[intersect.points+1] - ranMat[intersect.points]

  # Find the intersection for the segment.
  optimal <- intersect.points + ((ranMat[intersect.points] - x.var[intersect.points]) / (x1.slopes-x2.slopes))
  optimal1 <- abs(min(ceiling(intersect.points + ((ranMat[intersect.points] - x.var[intersect.points]) / (x1.slopes-x2.slopes))))) #the optimal number of components

  #Scree plot
  if(verbose == T | verbose == TRUE){
      pcaCharts(svd, average)
  }

  # Run ICA
  if(standardize == T | standardize == TRUE){
      data = scale(data)
  }

  ICA = fastICA(t(data), n.comp = optimal1, fun = "logcosh", alpha = 1,
                  method = "C", row.norm = FALSE, maxit = 1000,
                  tol = 0.0001)
  rownames(ICA$S) = colnames(data)

  # Remove signatures whose kurtosis <= 3
  kurt_ICA=ICA$S #define signatures of genes

  # Run IPCA  
  IPCA = mixOmics::ipca(data, ncomp = optimal1, fun = "logcosh", mode="deflation")

  # Remove component whose kurtosis <= 3
  kurt_IPCA=IPCA[["loadings"]][["X"]] #define signatures of genes

  #extract kurtosis
  cc2=list(My_names_is = paste("Huy", 1:ncol(kurt_ICA)), ICA = NA, IPCA = NA)
  for (i in 1:ncol(kurt_IPCA)){
      cc2$ICA[i] = kurtosis(kurt_ICA[,i])
      cc2$IPCA[i] = kurtosis(kurt_IPCA[,i])
  }

  #messenger
  #TH1
  if( (all(cc2$ICA < 3) == TRUE) & (any(cc2$IPCA > 3) == TRUE) ){
      cat(">> oCEM suggests choosing the optimal number of components is:", optimal1, "\n")
      cat(">> oCEM also suggests using IPCA-FDR for your case", "\n")
  }

  #TH2
  if( (all(cc2$ICA < 3) == TRUE) & (all(cc2$IPCA < 3) == TRUE) ){
      cat(">> oCEM suggests choosing the optimal number of components is:", optimal1, "\n")
      cat(">> oCEM cannot suggest which method should be selected. Please use a more stringent approach to make the best decision for your case.", "\n")
  }

  #TH3
  if( (any(cc2$ICA > 3) == TRUE) & (all(cc2$IPCA < 3) == TRUE) ){
      cat(">> oCEM suggests choosing the optimal number of components is:", optimal1, "\n")
      cat(">> oCEM also suggests using ICA for your case.", "\n")
  }

  #TH4
  if( (any(cc2$ICA > 3) == TRUE) & (any(cc2$IPCA > 3) == TRUE) ){
      cat(">> oCEM suggests choosing the optimal number of components is:", optimal1, "\n")
      cat(">> Both ICA and IPCA-FDR are appropriate for your case. Please use a more stringent approach to make the best decision.", "\n")
  }

  #time difference
  timediff = Sys.time() - now;
  mlog("Done in ", timediff, " ", units(timediff), ".\n")

  return(optimal1)
}