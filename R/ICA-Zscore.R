#' @title ICA-Zscore
#'
#' @description decomposition method is ICA and post-processing is Z-score.
#'
#' @docType package
#'
#' @author Quang-Huy Nguyen
#'
#' @usage overlapCEM(data = NULL, clinical = NULL, ncomp = NULL,
#'        method = 'ICA-Zscore', standardize = T, cex.text = 0.7)
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
#' @param cex.text numeric. Change the font size of texts in cells of the heatmap showing correlations between
#' each identified module and each clinical features. Default value is 0.7.
#'
#' @return NULL
#'
#' @export

ICAzscore=function(data = NULL, clinical = NULL, ncomp = NULL, standardize = T,
                   cex.text = 0.7){

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

  #library
  library(dplyr)
  library(moments)
  library(tidyr)
  library(WGCNA)

  now = Sys.time()

  # Run ICA
  set.seed(25896)
  if(standardize == T | standardize == TRUE){
    x  = scale(data)
  } else{
    x = data
  }

  ICA = fastICA(t(x), n.comp = ncomp, fun = "logcosh", alpha = 1,
                method = "C", row.norm = FALSE, maxit = 1000,
                tol = 0.0001)
  rownames(ICA$S) = colnames(data)

  # Remove signatures whose kurtosis < 3
  kurt=ICA$S #define signatures of genes
  pattern = ICA$A #define pattern of samples

  #extract kurtosis
  cc2=c()
  for (i in 1:ncol(kurt)){
    cc2[i] = kurtosis(kurt[,i])
  }
  if(all(cc2 < 3) == TRUE){
    stop('>> ICA-Zscore is inappropriate for your case due to all kurtosis < 3. Please use one of the other methods provided by oCEM')
  }

  #extract index of columns of signatures whose kurtosis < 3
  ind=data.frame(My_name_is=paste0("Huy", 1:length(cc2)), index = NA)
  for (j in 1:length(cc2)){
    if(cc2[j] < 3){
      ind$index[j] = j
    }
  }; ind=ind[!(is.na(ind$index)),]

  #messenger
  if(length(ind$index) >= 1){
    if(length(ind$index) == 1){
      cat(">> oCEM detected only one component whose kurtosis < 3", "\n")
      cat(">> oCEM excluded it in order to server to find out modules of genes, leaving", dim(kurt)[2] - 1, "components for further analyses", "\n")
      kurt = kurt[,-ind$index] #remove that signatures
      pattern = pattern[-ind$index,] #remove that pattern
    } else{
      cat(">> oCEM detected", length(ind$index), "components whose kurtosis < 3", "\n")
      cat(">> oCEM excluded them in order to serve to find out modules of genes, leaving", dim(kurt)[2] - length(ind$index), "components for further analyses", "\n")
      kurt = kurt[,-ind$index] #remove those signatures
      pattern = pattern[-ind$index,] #remove those patterns
    }
  } else{
    cat(">> oCEM did not detect any component whose kurtosis < 3", "\n")
    kurt = kurt
    pattern = pattern
  }

  # Identification of gene modules
  scaled_expression = scale(t(kurt)) %>% #do Z-score transformation on components (columns)
    t() %>%
    data.frame()
  scaled_exp = scaled_expression
  names(scaled_expression) = gsub("X", "module ", names(scaled_exp))
  View(scaled_expression)
  cat(">> As shown in Z-score table related to expression levels of each gene within each module, please choose a standard deviation threshold to let oCEM be able to assign genes to modules specifically.", "\n")
  cat(">> *NOTE that your selection will affect the minimum number of genes within a certain module.", "\n")
  threshold = as.numeric(readline('>> The standard deviation threshold you want to choose is?'))
  for (p in 1:nrow(scaled_exp)){
    scaled_exp[p,] = ifelse(abs(scaled_exp[p,]) > threshold, 1, 0)
  }
  remove_comp = scaled_exp[,(which(colSums(scaled_exp) == 0))] #signatures whose all elements in the columns equal 0
  #empty signatures will be removed
  if( (ncol(remove_comp) != 0) == TRUE ){
    ind1 = which(colSums(scaled_exp) == 0) #index of signatures whose all elements equal 0
    scaled_exp = scaled_exp[,-(ind1)] # remove columns whose all elements equal 0
    pattern = pattern[-(ind1),] # remove columns whose all elements equal 0
    kurt1 = kurt[,-ind1]
  } else{
    scaled_exp = scaled_exp
    pattern = pattern
    kurt1 = kurt
  }

  #check whether oCEM does create the grey module or not
  if(any(as.logical(rowSums(scaled_exp) == 0) == TRUE)){ #oCEM does create the grey module :(

    #genes not assigned to any modules will go to grey module
    ind3 = c()
    for (k in 1:nrow(scaled_exp)){
      if(rowSums(scaled_exp)[k] == 0){
        ind3[k] = k #position of row vector of 0's
        ind3 <- ind3[!is.na(ind3)]
        zeroes <- rep(0, nrow(scaled_exp))
        scaled_exp = data.frame(scaled_exp) %>%
          dplyr::mutate(gray = zeroes)
        scaled_exp$gray[ind3] = rowSums(scaled_exp)[ind3] + 1
      }
    }
    #convert numeric lables into colors
    names(scaled_exp) = gsub('X', '', names(scaled_exp), perl = T)
    names(scaled_exp) = gsub('gray', '0', names(scaled_exp), perl = T)
    moduleColors = WGCNA::labels2colors(as.numeric(names(scaled_exp)))

    #determine the exact number of genes distributed to a certain module
    cc3 = data.frame(My_name_is=paste0("Huy", 1:ncol(scaled_exp)), module = NA, num = NA)
    for (n in 1:ncol(scaled_exp)){
      cc3$module = moduleColors
      cc3$num[n] = length(scaled_exp[which(scaled_exp[,n] == 1), n])
      rownames(cc3) = cc3$module
    }
    cc3 = cc3 %>% dplyr::select(-c(My_name_is, module)) %>% t()

    #mesenger
    if(ncol(remove_comp) >0){ #remove_comp exists
      if(ncol(remove_comp) >1){
        cat(">> oCEM excluded", ncol(remove_comp), "components as they are empty (i.e., no genes are assigned to them)", "\n")
        cat(">> The exact number of genes assigned to each module is: ", "\n"); print(cc3)
      } else{
        cat(">> oCEM excluded", ncol(remove_comp), "component as it is empty (i.e., no genes are assigned to it)", "\n")
        cat(">> The exact number of genes assigned to each module is: ", "\n"); print(cc3)
      }
    } else{
      cat(">> The exact number of genes assigned to each module is: ", "\n"); print(cc3)
    }

    # Clinical feature-gene modules Assocation
    #define numbers of genes and samples
    nGenes = ncol(data);
    nSamples = nrow(data);
    #re-define patterns with color labels
    pattern1 = data.frame(t(pattern))
    names(pattern1) = moduleColors[1:(length(moduleColors)-1)]

    moduleTraitCor = cor(pattern1, clinical, use = "p");
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

    sizeGrWindow(20,6)
    # Will display correlations and their p-values
    textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                       signif(moduleTraitPvalue, 1), ")", sep = "");
    dim(textMatrix) = dim(moduleTraitCor)
    par(mar = c(6, 8.5, 3, 3));
    # Display the correlation values within a heatmap plot
    WGCNA::labeledHeatmap(Matrix = moduleTraitCor,
                          xLabels = names(clinical),
                          yLabels = moduleColors[1:(length(moduleColors)-1)],
                          ySymbols = names(pattern1),
                          colorLabels = FALSE,
                          colors = blueWhiteRed(50),
                          textMatrix = textMatrix,
                          setStdMargins = FALSE,
                          cex.text = cex.text,
                          zlim = c(-1,1),
                          main = paste("Module-clinical feature relationships"))

    # Identify hub-genes in each module: most significant genes ~ most extreme in the signature distribution
    cat("\n","- Starting to detect hub-genes within each discovered module...")
    #expression level of genes outside a certain module equal NA values
    ind2 = list()
    for(n in 1:(ncol(scaled_exp)-1)){
      ind2[[n]] = which(scaled_exp[,n] == 1)
      kurt1[-ind2[[n]],n] = NA_real_
    }

    #message
    topgenename = list()
    bottomgenename = list()
    for(m in 1:(ncol(cc3) - 1) ){
      if(isTRUE((cc3 < 20)[m])){
        if(isTRUE((cc3 < 10)[m])){
          cat("\n",">>>> Top 2 hub genes identified in the",unique(moduleColors)[[m]],"module are:",paste(names(which.max(kurt1[,m])), names(which.min(kurt1[,m])),sep = " and "), "\n")
        } else{
          if(max(kurt1[,m], na.rm = T) > abs(min(kurt1[,m], na.rm = T))){
            topgenename[[m]]<-rownames(dplyr::top_n(data.frame(kurt1[,m]),3, wt = data.frame(kurt1[,m])))
            bottomgenename[[m]]<-rownames(dplyr::top_n(data.frame(kurt1[,m]),-2, wt = data.frame(kurt1[,m])))
            cat("\n",">>>> Top 5 hub genes identified in the",unique(moduleColors)[[m]],"module are:", topgenename[[m]],"",bottomgenename[[m]], "\n")
          } else{
            bottomgenename[[m]]<-rownames(dplyr::top_n(data.frame(kurt1[,m]),3, wt = data.frame(kurt1[,m])))
            topgenename[[m]]<-rownames(dplyr::top_n(data.frame(kurt1[,m]),-2, wt = data.frame(kurt1[,m])))
            cat("\n",">>>> Top 5 hub genes identified in the",unique(moduleColors)[[m]],"module are:", bottomgenename[[m]],"",topgenename[[m]], "\n")
          }
        }
      } else{
        topgenename[[m]]<-rownames(dplyr::top_n(data.frame(kurt1[,m]),5, wt = data.frame(kurt1[,m])))
        bottomgenename[[m]]<-rownames(dplyr::top_n(data.frame(kurt1[,m]),-5, wt = data.frame(kurt1[,m])))
        cat("\n",">>>> Top 10 hub genes identified in the",unique(moduleColors)[[m]],"module are:",topgenename[[m]],"",bottomgenename[[m]], "\n")
      }
    }

    #Design results
    results1 = list()
    names(scaled_exp) = moduleColors
    results1 = c(results1,list(pattern1, scaled_exp))
    rownames(results1[[1]]) = rownames(data)
    names(results1) = c('patterns', 'signatures')
    return(results1)

  } else{ #IDEALLY, oCEM does not create the grey module

    #genes not assigned to any modules will go to grey module
    #convert numeric lables into colors
    names(scaled_exp) = gsub('X', '', names(scaled_exp), perl = T)
    moduleColors = WGCNA::labels2colors(as.numeric(names(scaled_exp)))

    #determine the exact number of genes distributed to a certain module
    cc3 = data.frame(My_name_is=paste0("Huy", 1:ncol(scaled_exp)), module = NA, num = NA)
    for (n in 1:ncol(scaled_exp)){
      cc3$module = moduleColors
      cc3$num[n] = length(scaled_exp[which(scaled_exp[,n] == 1), n])
      rownames(cc3) = cc3$module
    }
    cc3 = cc3 %>% dplyr::select(-c(My_name_is, module)) %>% t()

    #message
    if(ncol(remove_comp) >0){ #remove_comp exists
      if(ncol(remove_comp) >1){
        cat(">> oCEM excluded", ncol(remove_comp), "components as they are empty (i.e., no genes are assigned to them)", "\n")
        cat(">> The exact number of genes assigned to each module is: ", "\n"); print(cc3)
      } else{
        cat(">> oCEM excluded", ncol(remove_comp), "component as it is empty (i.e., no genes are assigned to it)", "\n")
        cat(">> The exact number of genes assigned to each module is: ", "\n"); print(cc3)
      }
    } else{
      cat(">> The exact number of genes assigned to each module is: ", "\n"); print(cc3)
    }

    # Clinical feature-gene modules Assocation
    #define numbers of genes and samples
    nGenes = ncol(data);
    nSamples = nrow(data);
    #re-define patterns with color labels
    pattern1 = data.frame(t(pattern))
    names(pattern1) = moduleColors[1:(length(moduleColors))]

    moduleTraitCor = cor(pattern1, clinical, use = "p");
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

    sizeGrWindow(20,6)
    # Will display correlations and their p-values
    textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                       signif(moduleTraitPvalue, 1), ")", sep = "");
    dim(textMatrix) = dim(moduleTraitCor)
    par(mar = c(6, 8.5, 3, 3));
    # Display the correlation values within a heatmap plot
    WGCNA::labeledHeatmap(Matrix = moduleTraitCor,
                          xLabels = names(clinical),
                          yLabels = moduleColors[1:(length(moduleColors))],
                          ySymbols = names(pattern1),
                          colorLabels = FALSE,
                          colors = blueWhiteRed(50),
                          textMatrix = textMatrix,
                          setStdMargins = FALSE,
                          cex.text = cex.text,
                          zlim = c(-1,1),
                          main = paste("Module-clinical feature relationships"))

    # Identify hub-genes in each module: most significant genes ~ most extreme in the signature distribution
    cat("\n","- Starting to detect hub-genes within each discovered module...")
    #expression level of genes outside a certain module equal NA values
    ind2 = list()
    for(n in 1:(ncol(scaled_exp))){
      ind2[[n]] = which(scaled_exp[,n] == 1)
      kurt1[-ind2[[n]],n] = NA_real_
    }

    #messenger
    topgenename = list()
    bottomgenename = list()
    for(m in 1:ncol(cc3) ){
      if(isTRUE((cc3 < 20)[m])){
        if(isTRUE((cc3 < 10)[m])){
          cat("\n",">>>> Top 2 hub genes identified in the",unique(moduleColors)[[m]],"module are:",paste(names(which.max(kurt1[,m])), names(which.min(kurt1[,m])),sep = " and "), "\n")
        } else{
          if(max(kurt1[,m], na.rm = T) > abs(min(kurt1[,m], na.rm = T))){
            topgenename[[m]]<-rownames(dplyr::top_n(data.frame(kurt1[,m]),3, wt = data.frame(kurt1[,m])))
            bottomgenename[[m]]<-rownames(dplyr::top_n(data.frame(kurt1[,m]),-2, wt = data.frame(kurt1[,m])))
            cat("\n",">>>> Top 5 hub genes identified in the",unique(moduleColors)[[m]],"module are:", topgenename[[m]],"",bottomgenename[[m]], "\n")
          } else{
            bottomgenename[[m]]<-rownames(dplyr::top_n(data.frame(kurt1[,m]),3, wt = data.frame(kurt1[,m])))
            topgenename[[m]]<-rownames(dplyr::top_n(data.frame(kurt1[,m]),-2, wt = data.frame(kurt1[,m])))
            cat("\n",">>>> Top 5 hub genes identified in the",unique(moduleColors)[[m]],"module are:", bottomgenename[[m]],"",topgenename[[m]], "\n")
          }
        }
      } else{
        topgenename[[m]]<-rownames(dplyr::top_n(data.frame(kurt1[,m]),5, wt = data.frame(kurt1[,m])))
        bottomgenename[[m]]<-rownames(dplyr::top_n(data.frame(kurt1[,m]),-5, wt = data.frame(kurt1[,m])))
        cat("\n",">>>> Top 10 hub genes identified in the",unique(moduleColors)[[m]],"module are:",topgenename[[m]],"",bottomgenename[[m]], "\n")
      }
    }

    #Design results
    results1 = list()
    names(scaled_exp) = moduleColors
    results1 = c(results1,list(pattern1, scaled_exp))
    rownames(results1[[1]]) = rownames(data)
    names(results1) = c('patterns', 'signatures')
    return(results1)

  }
}
