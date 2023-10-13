#' @title ICA-FDR
#'
#' @description decomposition method is ICA and post-processing is FDR.
#'
#' @docType package
#'
#' @author Quang-Huy Nguyen
#'
#' @usage overlapCEM(data = NULL, clinical = NULL, ncomp = NULL,
#'        method = 'ICA-FDR', standardize = T, cex.text = 0.7, verbose = T)
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
#' @param verbose logical. Show the running time of this method. Default value is T.
#'
#' @return NULL
#'
#' @export

ICAfdr=function(data = NULL, clinical=NULL, ncomp = NULL, standardize = T,
                cex.text = 0.7, verbose = T){ 
  #system time
  now = Sys.time()

  # defined log function
  mlog <- if(!verbose) function(...){} else function(...){
    message(...)
    flush.console()
  }
  
  # Run ICA
  set.seed(25896)
  if(standardize == T | standardize == TRUE){
    data  = scale(data)
  }

  ICA = fastICA(t(data), n.comp = ncomp, fun = "logcosh", alpha = 1,
                method = "C", row.norm = FALSE, maxit = 1000,
                tol = 0.0001)
  rownames(ICA$S) = colnames(data)

  # Remove signatures whose kurtosis < 3
  kurt=ICA$S #define signatures of genes
  pattern = ICA$A #define pattern of samples

  #extract kurtosis
  kurtosis_res <- apply(kurt, 2, moments::kurtosis)
  if(all(kurtosis_res < 3) == TRUE){
    stop('>> ICA-FDR is inappropriate for your case due to all kurtosis < 3. Please use one of the other methods provided by oCEM')
  }

  #detect index of columns of signatures whose kurtosis < 3
  ind <- which( kurtosis_res <3 )
  if(length(ind) >= 1){
    if(length(ind) == 1){
      cat(">> oCEM detected only 1 component whose kurtosis < 3", "\n")
      cat(">> oCEM excluded it in order to server to find out modules of genes, leaving", dim(kurt)[2] - 1, "components for further analyses", "\n")
      kurt = kurt[,-ind] #remove that signatures
      pattern = pattern[-ind,] #remove that pattern
    } else{
      cat(">> oCEM detected", length(ind), "components whose kurtosis < 3", "\n")
      cat(">> oCEM excluded them in order to serve to find out modules of genes, leaving", dim(kurt)[2] - length(ind), "components for further analyses", "\n")
      kurt = kurt[,-ind] #remove those signatures
      pattern = pattern[-ind,] #remove those patterns
    }
  } else{
    cat(">> oCEM did not detect any component whose kurtosis < 3", "\n")
  }

  # Identification of gene modules
  fdr = apply(kurt,2,function(x) fdrtool::fdrtool(x, plot = F, verbose = F, color.figure = F)$pval) %>% as.data.frame()
  names(fdr) = gsub("V", "module", names(fdr))
  cat(">> Please choose a probability threshold to let oCEM be able to assign genes to modules specifically.", "\n")
  cat(">> *NOTE that your selection will affect the minimum number of genes within a certain module.", "\n")
  threshold = as.numeric(readline('>> The probability threshold you want to choose is? (e.g., 0.1, 0.01,...)'))
  fdr <- ifelse(fdr < threshold, 1, 0) %>% as.data.frame()

  #signatures whose all elements in the columns equal 0
  remove_comp <- which(colSums(fdr) == 0) 

  #In case the user chooses a too small threshold, all signatures will have their all entries in every column that equal 0
  test_fdr <- fdr[,-remove_comp]
  if( (length(remove_comp) > 0) && (ncol(test_fdr) == 0) ){
    stop("Threshold of your choice is too small, please choose another bigger one again!")
  }

  store<-fdr
  fdr<-store
  #empty signatures will be removed
  if( (length(remove_comp) > 0) ){ #have colSums = 0
    fdr <- fdr %>% dplyr::select(-all_of(names(remove_comp))) # remove columns whose all elements equal 0
    pattern = pattern[-remove_comp,] # remove columns whose all elements equal 0
    kurt1 = kurt[,-remove_comp]
  } else{ #have no colSums = 0
    kurt1 = kurt
  }

  #genes not assigned to any modules will go to grey module
  if( any(rowSums(fdr) == 0) ){
    ind3 <- which(rowSums(fdr) == 0)
    grey_module <- rep(0, nrow(fdr))
    for (i in ind3){
      grey_module[i] <- grey_module[[i]] + 1
    }
    fdr = data.frame(fdr) %>% dplyr::mutate("0" = grey_module)
  }
  colnames(fdr) <- gsub("module", "", colnames(fdr))

  #convert numeric lables into colors
  names(fdr) = gsub('X', '', names(fdr), perl = T)
  moduleColors <- WGCNA::labels2colors(as.numeric(colnames(fdr)))
  colnames(fdr) <-  moduleColors

  #determine the exact number of genes distributed to a certain module
  component_summary <- colSums(fdr)

  #mesenger
  if( length(remove_comp) > 0 ){ #remove_comp exists
    if(length(remove_comp) > 1){
      cat(">> oCEM excluded", length(remove_comp), "components as they are empty (i.e., no genes are assigned to them)", "\n")
      cat(">> The exact number of genes assigned to each module is: ", "\n"); print(component_summary)
    } else{
      cat(">> oCEM excluded", length(remove_comp), "component as it is empty (i.e., no genes are assigned to it)", "\n")
      cat(">> The exact number of genes assigned to each module is: ", "\n"); print(component_summary)
    }
  } else{
    cat(">> The exact number of genes assigned to each module is: ", "\n"); print(component_summary)
  }

  # Clinical feature-gene modules Assocation
  if( !is.null(clinical) ){
    #define numbers of genes and samples
    nGenes = ncol(data);
    nSamples = nrow(data);
    #re-define patterns with color labels
    pattern1 = data.frame(t(pattern))
    if( "grey" %in% colnames(fdr) ){
      names(pattern1) <- moduleColors[1:(length(moduleColors)-1)]
    } else{
      names(pattern1) <- moduleColors
    }
    
    moduleTraitCor = cor(pattern1, clinical, use = "p");
    moduleTraitPvalue = WGCNA::corPvalueStudent(moduleTraitCor, nSamples);

    WGCNA::sizeGrWindow(20,6)
    # Will display correlations and their p-values
    textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                      signif(moduleTraitPvalue, 1), ")", sep = "");
    dim(textMatrix) = dim(moduleTraitCor)
    par(mar = c(6, 8.5, 3, 3));
    # Display the correlation values within a heatmap plot
    WGCNA::labeledHeatmap(Matrix = moduleTraitCor,
                          xLabels = names(clinical),
                          yLabels = names(pattern1),
                          ySymbols = names(pattern1),
                          colorLabels = FALSE,
                          colors = blueWhiteRed(50),
                          textMatrix = textMatrix,
                          setStdMargins = FALSE,
                          cex.text = cex.text,
                          zlim = c(-1,1),
                          main = paste("Module-clinical feature relationships"))
  } else{
    moduleTraitCor = NULL
    moduleTraitPvalue = NULL
  }

  # Identify hub-genes in each module: most significant genes ~ most extreme in the signature distribution
  cat("\n","- Starting to detect hub-genes within each discovered module...")
  #make values of every gene sitting outside certain modules NA
  if( "grey" %in% colnames(fdr) ){
    fdr_without_grey <- fdr[,( 1:(ncol(fdr) -1) )]
    kurt1 <- kurt1 * fdr_without_grey %>% as.data.frame() %>% mutate_all(., function(x) as.numeric(as.character(x)))
    kurt1[kurt1 == 0] <- NA_integer_
  } else{
    kurt1 <- kurt1 * fdr %>% as.data.frame() %>% mutate_all(., function(x) as.numeric(as.character(x)))
    kurt1[kurt1 == 0] <- NA_integer_
  }

  #messenger
  topgenename = list()
  bottomgenename = list()    
  hub_genes_in_each_module <- lapply(1:ncol(kurt1), FUN=function(m){
    if( table(is.na(kurt1[,m]))[1] < 11 ){
      if( table(is.na(kurt1[,m]))[1] < 6 ){
        topgenename[[m]] <- rownames(kurt1[order(-kurt1[,m]), ])[1]
        bottomgenename[[m]] <- rownames(kurt1[order(kurt1[,m]), ])[1]   
        union_hub_genes <- union(topgenename[[m]], bottomgenename[[m]])
        cat("\n",">>>> Top hub genes identified in the",unique(moduleColors)[[m]],"module are:", union_hub_genes, "\n")
      } else{
        if(max(kurt1[,m], na.rm = T) > abs(min(kurt1[,m], na.rm = T))){
          topgenename[[m]] <- rownames(kurt1[order(-kurt1[,m]), ])[1:3]
          bottomgenename[[m]] <- rownames(kurt1[order(kurt1[,m]), ])[1:2]
          union_hub_genes <- union(topgenename[[m]], bottomgenename[[m]]) 
          cat("\n",">>>> Top 5 hub genes identified in the",unique(moduleColors)[[m]],"module are:", union_hub_genes, "\n")
        } else{
          topgenename[[m]] <- rownames(kurt1[order(-kurt1[,m]), ])[1:2]
          bottomgenename[[m]] <- rownames(kurt1[order(kurt1[,m]), ])[1:3]
          union_hub_genes <- union(topgenename[[m]], bottomgenename[[m]]) 
          cat("\n",">>>> Top 5 hub genes identified in the",unique(moduleColors)[[m]],"module are:", union_hub_genes, "\n")
        }
      }
    } else{
      topgenename[[m]] <- rownames(kurt1[order(-kurt1[,m]), ])[1:5]
      bottomgenename[[m]] <- rownames(kurt1[order(kurt1[,m]), ])[1:5]
      union_hub_genes <- union(topgenename[[m]], bottomgenename[[m]])
      cat("\n",">>>> Top 10 hub genes identified in the",unique(moduleColors)[[m]],"module are:", union_hub_genes, "\n")
    }
    res <- list(
      module = unique(moduleColors)[[m]],
      hubgenes = union_hub_genes
    )
  }) %>% do.call(what=rbind)

  #time difference
  timediff = Sys.time() - now;
  mlog("Done in ", timediff, " ", units(timediff), ".\n")

  res <- list(modules = as.data.frame(fdr),
              moduleTraitCor = moduleTraitCor,
              moduleTraitPvalue = moduleTraitPvalue,
              hubgenes = hub_genes_in_each_module)
  return(res)
}
