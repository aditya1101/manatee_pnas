#' Forward Search Function
#' @description
#' Forward search is useful for reducing the size of the gene set in your
#' \code{filterObject}. In general, forward search identifies a small set of genes with
#' maximum ability to distinguish cases from controls. \cr
#' 
#' Forward search is a method of optimizing a given set of significant genes to
#' maximize some metric of diagnostic performance. The function works by taking
#' a given set of genes (presumably a set that has been filtered for statistical
#' significance), and iteratively adding one gene at a time, until the stopping
#' threshold is reached. At each round, the gene whose addition contributes the
#' greatest increase in the designated performance metric is added.
#' @param pooledDataObject A MANATEE pooledDataObject
#' @param filterObject A MANATEE filterObject - if you just have a list of up/down genes, then you can use the upgenes/downgenes arguements
#' @param upgenes A vector of upgenes (if filterObject is provided, then it will override this)
#' @param downgenes A vector of downgenes (if filterObject is provided, then it will override this)
#' @param perf.meas Which diagnostic performance metric to use. Options are "AUC", "AUPRC", "AvgAUC", and "Specificity" (Default: "AUC")
#' @param forwardThresh Stopping threshold. As long as the performance of the best N genes is better than the performance
#' of the best N-1 genes plus the forwardThresh, then the forward search will continue. For example, if forwardThresh is 0.01,
#' then the performance of the best N genes must increase by at least 0.01 each time in order to continue (Default: 0)
#' @param weights \emph{Only for perf.meas=AvgAUC:} If left blank, then the simple average of class-specific AUCs will be used to
#' calculate the AvgAUC. However, if you want to take a weighted average instead, then pass in a vector of weights. Note that the
#' order of the weights will be matched to the order of the names in the \code{otherNames} vector.
#' @param caseNames \emph{Only for perf.meas=AvgAUC:} name of the class(es) that you're considering to be your case.
#' If left blank, then \code{caseNames} will be the unique set of names in \code{$pheno$group} that correspond to a class label of 1
#' @param otherNames \emph{Only for perf.meas=AvgAUC:} name of the class(es) that you're considering to be your controls.
#' If left blank, then \code{otherNames} will be the unique set of names in \code{$pheno$group} that correspond to a class label of 0.
#' The order of otherNames will be set by running the \code{sort} function, so if you are using weights, make sure that your \code{weights}
#' vector and \code{otherNames} vector match up correctly
#' @param group.colname \emph{Only for perf.meas=AvgAUC:} this designates the column name of the group vector within the pooledDataObject$pheno
#' dataframe (Default: "group")
#' @param max.genes If you don't want to look beyond a certain number of genes, then set \code{max.genes}
#' as the max number of genes to include in the resulting signature
#' @param init.pos If you want to start off with some upgenes already locked in, then indicate those genes with \code{init.pos}
#' @param init.neg If you want to start off with some downgenes genes already locked in, then indicate those genes with \code{init.pos}
#' @param force.posneg If you want to make sure that there are both positive and negative genes in the signature, you can set this to TRUE
#' and it will make sure that the first two genes selected are in different directions. This will not work if either init.pos or init.neg is being used.
#' @param replace.genes If this is set to TRUE, then genes will be replaced and can thus multiple copies of the same gene can be chosen.\
#' @param sensVec If perf.meas="Specificity" then this provides a vector of sensitivities at which to assess the specificity
#' @param numCores Number of cores to use for parallelization. If left blank, smart core usage will be used
#' @details
#' The runForwardSearch and runBackwardSearch functions are designed to assist
#' in selection of gene sets optimized for discriminatory power. The selection
#' of an optimized set is a non-convex problem, and hence both functions will
#' yield gene sets that are only locally optimized (ie, they are not global
#' optima). If the global optimum is desired, then a best subset selection
#' should be run instead, using the \code{bestSubsetSelection}
#' function. Both the runForwardSearch and runBackwardSearch functions follow a
#' greedy algorithm, either adding (or removing) genes that contribute the most
#' (or the least) to the designated performance metric in the discovery data. \cr
#' 
#' Both search functions allow a user to set a stopping threshold; the
#' fundamental tradeoff here will be sparsity of the returned gene set vs.
#' overall discriminatory power. The default threshold is 0, such the functions
#' will return the set of genes at which no gene could be added or removed for
#' the forward or backward functions, respectively, and improve the designated
#' performance metric.
#' @note As of now, the accepted diagnostic performance metrics are AUC (area under the Reciever Operating Curve),
#' AUPRC (area under the Precision Recall Curve), and AvgAUC
#' (the average AUC across all comparisons between your cases and each of your "control" categories)
#' @return
#'   \item{upgenes}{Vector of the upgenes that were chosen by the forward search}
#'   \item{downgenes}{Vector of the downgenes that were chosen by the forward search}
#' @examples
#'
#' @export
#' @import parallel zoo pracma MetaIntegrator ROCR
#' @author Aditya Rao
runForwardSearch <- function(pooledDataObject, filterObject=NULL, upgenes=NULL, downgenes=NULL, perf.meas="AUC", forwardThresh=0, weights=NULL,
                             caseNames = NULL, otherNames = NULL, group.colname = "group", max.genes=NULL, init.pos=NULL, init.neg=NULL,
                             force.posneg=TRUE, replace.genes=FALSE, sensVec = c(seq(0.9,0.99,0.01)), numCores=NULL){
  if(!checkManateeObject(pooledDataObject,"pooledDataObject")){stop("Invalid pooledDataObject provided")}
  if(!is.null(filterObject)){
    if(!checkManateeObject(filterObject,"filterObject")){
      if(is.null(filterObject$upGeneNames) && is.null(filterObject$downGeneNames)){
        stop("Invalid filterObject provided. If you are building a filterObject from scratch, look at ?runForwardSearch for tips.")
      }
    }
  }else if(!is.null(upgenes) || !is.null(downgenes)){
    filterObject = list(upGeneNames = upgenes, downGeneNames = downgenes)
  }else{
    stop("No filterObject, upgenes, or downgenes provided")
  }
  
  if(any(duplicated(filterObject$upGeneNames)) || any(duplicated(filterObject$downGeneNames))){
    stop("There are duplicate gene names provided, this unfortunately causes a bug with the function. Use the replace.genes argument instead.")
  }
  
  if(any(sensVec > 1) || any(sensVec < 0)){stop("Sensitivities must be within 0 and 1")}
  
  perf.meas = tolower(perf.meas)
  perf.meas = switch(perf.meas, "auc"="AUC", "auprc"="AUPRC", "avgauc"="AvgAUC", "specificity"="Specificity")
  if(!(perf.meas %in% c("AUC","AUPRC","AvgAUC","Specificity"))){stop("perf.meas must be \"AUC\", \"AUPRC\", \"AvgAUC\", or \"Specificity\"")}
  
  if(perf.meas=="AvgAUC"){
    forwardResults = .forwardSearchHelper(pos.genes=filterObject$upGeneNames, neg.genes=filterObject$downGeneNames, perf.meas = perf.meas, forwardThresh = forwardThresh, 
                                          geneMtx = pooledDataObject$genes, class = pooledDataObject$class, labelVec=pooledDataObject$pheno[,group.colname], weights=weights,
                                          caseNames=caseNames, otherNames=otherNames, max.genes=max.genes, init.pos=init.pos, init.neg=init.neg, force.posneg=force.posneg, replace.genes=replace.genes, numCores=numCores)
  }else{
    forwardResults = .forwardSearchHelper(pos.genes=filterObject$upGeneNames, neg.genes=filterObject$downGeneNames, perf.meas = perf.meas, forwardThresh = forwardThresh,
                                          geneMtx = pooledDataObject$genes, class = pooledDataObject$class, max.genes=max.genes, init.pos=init.pos, init.neg=init.neg, force.posneg=force.posneg, 
                                          replace.genes=replace.genes, sensVec=sensVec, numCores=numCores)
  }
  
  return(forwardResults)
}



#' Backward Search Function
#' @description
#' Backward search is useful for reducing the size of the gene set in your
#' \code{filterObject}. In general, forward search identifies a small set of genes with
#' maximum ability to distinguish cases from controls. \cr
#' 
#' Backward search is a method of optimizing a given set of significant genes to
#' maximize some metric of diagnostic performance. The function works by taking
#' a given set of genes (presumably a set that has been filtered for statistical
#' significance), and iteratively removing one gene at a time, until the stopping
#' threshold is reached. At each round, the gene whose removal contributes the
#' greatest increase in the designated performance metric is added.
#' @param pooledDataObject A MANATEE pooledDataObject
#' @param filterObject A MANATEE filterObject - if you just have a list of up/down genes, then you can use the upgenes/downgenes arguements
#' @param upgenes A vector of upgenes (if filterObject is provided, then it will override this)
#' @param downgenes A vector of downgenes (if filterObject is provided, then it will override this)
#' @param perf.meas Which diagnostic performance metric to use. Options are "AUC", "AUPRC", and "AvgAUC" (Default: "AUC")
#' @param backwardThresh Stopping threshold. As long as the performance of the best N genes is better than the performance of the best N+1 genes
#' plus the backwardThresh, then the backward search will continue. For example, if backwardThresh is -0.01, 
#' then the performance of the best N genes can only decrease by 0.01 or less each time in order to continue (Default: 0)
#' @param weights \emph{Only for perf.meas=AvgAUC:} If left blank, then the simple average of class-specific AUCs will be used to
#' calculate the AvgAUC. However, if you want to take a weighted average instead, then pass in a vector of weights. Note that the
#' order of the weights will be matched to the order of the names in the \code{otherNames} vector.
#' @param caseNames \emph{Only for perf.meas=AvgAUC:} name of the class(es) that you're considering to be your case.
#' If left blank, then \code{caseNames} will be the unique set of names in \code{$pheno$group} that correspond to a class label of 1
#' @param otherNames \emph{Only for perf.meas=AvgAUC:} name of the class(es) that you're considering to be your controls.
#' If left blank, then \code{otherNames} will be the unique set of names in \code{$pheno$group} that correspond to a class label of 0.
#' The order of otherNames will be set by running the \code{sort} function, so if you are using weights, make sure that your \code{weights}
#' vector and \code{otherNames} vector match up correctly
#' @param group.colname \emph{Only for perf.meas=AvgAUC:} this designates the column name of the group vector within the pooledDataObject$pheno
#' dataframe (Default: "group")
#' @param min.genes If you don't want to look beyond a certain number of genes, then set \code{min.genes}
#' as the minimum number of genes to include in the resulting signature
#' @param numCores Number of cores to use for parallelization. If left blank, smart core usage will be used
#' @details
#' The runForwardSearch and runBackwardSearch functions are designed to assist
#' in selection of gene sets optimized for discriminatory power. The selection
#' of an optimized set is a non-convex problem, and hence both functions will
#' yield gene sets that are only locally optimized (ie, they are not global
#' optima). If the global optimum is desired, then a best subset selection
#' should be run instead, using the \code{bestSubsetSelection}
#' function. Both the runForwardSearch and runBackwardSearch functions follow a
#' greedy algorithm, either adding (or removing) genes that contribute the most
#' (or the least) to the designated performance metric in the discovery data. \cr
#' 
#' Both search functions allow a user to set a stopping threshold; the
#' fundamental tradeoff here will be sparsity of the returned gene set vs.
#' overall discriminatory power. The default threshold is 0, such the functions
#' will return the set of genes at which no gene could be added or removed for
#' the forward or backward functions, respectively, and improve the designated
#' performance metric.
#' @note As of now, the accepted diagnostic performance metrics are AUC (area under the Reciever Operating Curve),
#' AUPRC (area under the Precision Recall Curve), and AvgAUC
#' (the average AUC across all comparisons between your cases and each of your "control" categories)
#' @return
#'   \item{upgenes}{Vector of the upgenes that were chosen by the backward search}
#'   \item{downgenes}{Vector of the downgenes that were chosen by the backward search}
#' @examples
#'
#' @export
#' @import parallel zoo pracma MetaIntegrator ROCR
#' @author Aditya Rao
runBackwardSearch <- function(pooledDataObject, filterObject=NULL, upgenes=NULL, downgenes=NULL, perf.meas="AUC", backwardThresh=0, weights=NULL,
                              caseNames = NULL, otherNames = NULL, group.colname = "group", min.genes=NULL, numCores=NULL){
  if(!checkManateeObject(pooledDataObject,"pooledDataObject")){stop("Invalid pooledDataObject provided")}
  if(!is.null(filterObject)){
    if(!checkManateeObject(filterObject,"filterObject")){
      if(is.null(filterObject$upGeneNames) && is.null(filterObject$downGeneNames)){
        stop("Invalid filterObject provided. If you are building a filterObject from scratch, look at ?runForwardSearch for tips.")
      }
    }
  }else if(!is.null(upgenes) || !is.null(downgenes)){
    filterObject = list(upGeneNames = upgenes, downGeneNames = downgenes)
  }else{
    stop("No filterObject, upgenes, or downgenes provided")
  }
  
  if(any(duplicated(filterObject$upGeneNames)) || any(duplicated(filterObject$downGeneNames))){
    stop("There are duplicate gene names provided, this unfortunately causes a bug with the function. Use the replace.genes argument instead.")
  }
  
  perf.meas = tolower(perf.meas)
  perf.meas = switch(perf.meas, "auc"="AUC", "auprc"="AUPRC", "avgauc"="AvgAUC")
  if(!(perf.meas %in% c("AUC","AUPRC","AvgAUC"))){stop("perf.meas must be \"AUC\", \"AUPRC\", or \"AvgAUC\"")}
  
  if(perf.meas=="AvgAUC"){
    backwardResults = .backwardSearchHelper(pos.genes=filterObject$upGeneNames, neg.genes=filterObject$downGeneNames, perf.meas = perf.meas, backwardThresh = backwardThresh, 
                                            geneMtx = pooledDataObject$genes, class = pooledDataObject$class, labelVec=pooledDataObject$pheno[,group.colname], weights=weights,
                                            caseNames=caseNames, otherNames=otherNames, min.genes=min.genes, numCores=numCores)
  }else{
    backwardResults = .backwardSearchHelper(pos.genes=filterObject$upGeneNames, neg.genes=filterObject$downGeneNames, perf.meas = perf.meas, backwardThresh = backwardThresh,
                                            geneMtx = pooledDataObject$genes, class = pooledDataObject$class, min.genes=min.genes, numCores=numCores)
  }
  
  return(backwardResults)
  
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#--------------------------------------------------------Helper Functions--------------------------------------------------------#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' Helper Function for Forward Search
#' @param pos.genes Vector of upgenes
#' @param neg.genes Vector of downgenes
#' @param perf.meas Which diagnostic performance metric to use. Options are "AUC", "AUPRC", "AvgAUC", and "Specificity" (Default: "AUC")
#' @param forwardThresh Stopping threshold. As long as the performance of the best N genes is better than the performance
#' of the best N-1 genes plus the forwardThresh, then the forward search will continue. For example, if forwardThresh is 0.01,
#' then the performance of the best N genes must increase by at least 0.01 each time in order to continue (Default: 0)
#' @param geneMtx The gene matrix to use for calculating signature scores (rows are genes and columns are samples)
#' @param class Vector of class labels
#' @param labelVec \emph{Only for perf.meas=AvgAUC:} Vector of sample category labels
#' @param weights \emph{Only for perf.meas=AvgAUC:} If left blank, then the simple average of class-specific AUCs will be used to
#' calculate the AvgAUC. However, if you want to take a weighted average instead, then pass in a vector of weights. Note that the
#' order of the weights will be matched to the order of the names in the \code{otherNames} vector.
#' @param caseNames \emph{Only for perf.meas=AvgAUC:} name of the class(es) that you're considering to be your case.
#' If left blank, then \code{caseNames} will be the unique set of names in \code{$pheno$group} that correspond to a class label of 1
#' @param otherNames \emph{Only for perf.meas=AvgAUC:} name of the class(es) that you're considering to be your controls.
#' If left blank, then \code{otherNames} will be the unique set of names in \code{$pheno$group} that correspond to a class label of 0.
#' The order of otherNames will be set by running the \code{sort} function, so if you are using weights, make sure that your \code{weights}
#' vector and \code{otherNames} vector match up correctly
#' @param max.genes If you don't want to look beyond a certain number of genes, then set \code{max.genes}
#' as the max number of genes to include in the resulting signature
#' @param init.pos If you want to start off with some upgenes already locked in, then indicate those genes with \code{init.pos}
#' @param init.neg If you want to start off with some downgenes genes already locked in, then indicate those genes with \code{init.pos}
#' @param force.posneg If you want to make sure that there are both positive and negative genes in the signature, you can set this to TRUE
#' and it will make sure that the first two genes selected are in different directions. This will not work if either init.pos or init.neg is being used.
#' @param replace.genes If this is set to TRUE, then genes will be replaced and can thus multiple copies of the same gene can be chosen.
#' @param sensVec If perf.meas="Specificity" then this provides a vector of sensitivities at which to assess the specificity
#' @param numCores Number of cores to use for parallelization. If left blank, smart core usage will be used
#' @note As of now, the accepted diagnostic performance metrics are AUC (area under the Reciever Operating Curve),
#' AUPRC (area under the Precision Recall Curve), and AvgAUC
#' (the average AUC across all comparisons between your cases and each of your "control" categories)
#' @return
#'   \item{upgenes}{Vector of the upgenes that were chosen by the forward search}
#'   \item{downgenes}{Vector of the downgenes that were chosen by the forward search}
#' @import parallel zoo pracma MetaIntegrator ROCR
#' @author Aditya Rao
.forwardSearchHelper <- function(pos.genes, neg.genes, perf.meas = "AUC", forwardThresh = 0, geneMtx, class, labelVec, weights=NULL, caseNames=NULL, otherNames=NULL,
                                 max.genes=NULL, init.pos=NULL, init.neg=NULL, force.posneg=TRUE, replace.genes=FALSE, sensVec = c(seq(0.9,0.99,0.01)), numCores=NULL){
  perf.meas = tolower(perf.meas)
  perf.meas = switch(perf.meas, "auc"="AUC", "auprc"="AUPRC", "avgauc"="AvgAUC", "specificity"="Specificity")
  if(!(perf.meas %in% c("AUC","AUPRC","AvgAUC","Specificity"))){stop("perf.meas must be \"AUC\", \"AUPRC\", \"AvgAUC\", or \"Specificity\"")}
  if(is.null(numCores)){
    #smart core usage = use 50% of resources unless you need less
    numCores <- round(detectCores()*5/10)
    if(length(c(pos.genes,neg.genes))<numCores){
      numCores <- length(c(pos.genes,neg.genes))
    }
  }
  if(any(!c(pos.genes,neg.genes) %in% rownames(geneMtx))){stop("Some of the designated genes are not present in the provided gene matrix")}
  if(any(!init.pos %in% pos.genes)){stop("There are some genes in init.pos that do not appear in pos.genes")}
  if(any(!init.neg %in% neg.genes)){stop("There are some genes in init.neg that do not appear in neg.genes")}
  if(is.null(init.pos)){
    pos=NULL
  }else{
    pos=init.pos
    if(!replace.genes){
      pos.genes = pos.genes[!pos.genes %in% init.pos]
    }
  }
  if(is.null(init.neg)){
    neg=NULL
  }else{
    neg=init.neg
    if(!replace.genes){
      neg.genes = neg.genes[!neg.genes %in% init.neg]
    }
  }
  if(is.null(pos) && is.null(neg)){
    count=0
  }else{
    count = length(pos) + length(neg)
  }
  if(force.posneg){
    if(!is.null(init.pos) || !is.null(init.neg)){
      warning("force.posneg cannot be used if init.pos or init.neg is set, so it will be ignored")
      force.posneg = FALSE
    }else if(length(pos.genes) == 0 || length(neg.genes) == 0){
      warning("force.posneg cannot be used if pos.genes or neg.genes is empty, so it will be ignored")
      force.posneg = FALSE
    }
  }
  
  perf.max = 0
  continue=TRUE
  genelist=c(pos.genes,neg.genes)
  bestgenelist = list()
  if(is.null(max.genes)){
    max.genes = length(pos.genes) + length(neg.genes)
    if(replace.genes){max.genes = 999999999}
  }else if(max.genes > (length(pos.genes) + length(neg.genes))){
    warning("The provided max.genes is greater than the total number of genes and will be ignored")
    max.genes = length(pos.genes) + length(neg.genes)
  }else if(max.genes < count){
    stop("max.genes cannot be less than the number of initial genes")
  }
  
  if(perf.meas=="AUC"){
    while(continue && (count < max.genes)){
      auclist=mclapply(mc.cores=numCores, genelist, function(gene){
        currpos=pos
        currneg=neg
        if(gene %in% pos.genes){
          currpos=append(currpos,gene)
        }else{
          currneg=append(currneg,gene)
        }
        scores = getGeneScores(geneMtx,currpos,currneg)
        auc = as.numeric(calculateROC(as.numeric(as.character(class)), as.numeric(scores))$auc)
        
        #hopefully this works
        if(force.posneg && count == 1){
          if(is.null(pos) && gene %in% neg.genes){auc=0}
          if(is.null(neg) && gene %in% pos.genes){auc=0}
        }
        return(auc)
      })
      #check for NULL values
      if(any(sapply(auclist,function(x) is.null(x)))){
        warning("Some values of auclist are NULL - this is due to an mclapply error. Reduce numCores and try again.")
      }
      auclist = unlist(auclist)
      curr.max = max(auclist,na.rm=T)
      perf.diff = curr.max-perf.max
      cat(sprintf("next best: %s\n",perf.diff))
      count = count+1
      if(perf.diff>forwardThresh){
        perf.max = curr.max
        topindex = which.max(auclist) #not 100% sure how this will behave if there is a tie
        if(length(topindex)>1){topindex=topindex[1]} #just make sure topindex is only one index
        topgene=genelist[topindex]
        if(topgene %in% pos.genes){
          cat(sprintf("Adding %s (up)\n",topgene))
          pos=append(pos,topgene)
        }else{
          cat(sprintf("Adding %s (down)\n",topgene))
          neg=append(neg,topgene)
        }
        if(!replace.genes){
          genelist=genelist[-topindex]
        }
        best.scores = getGeneScores(geneMtx,pos,neg)
        best.perf = calculateROC(as.numeric(as.character(class)), as.numeric(best.scores))
        best.auc = as.numeric(best.perf$auc)
        blurb = sprintf("%s Genes, %s=%s (95%% CI %s-%s): %s up and %s down",count,perf.meas,round(best.auc,3),round(best.perf$auc.CI[1],3),
                        round(best.perf$auc.CI[2],3),paste(pos,collapse=", "),paste(neg,collapse=", "))
        bestgenelist[[as.character(count)]] = list(numGenes = count, upgenes = pos, downgenes = neg, perf.meas = perf.meas, perf = best.auc,
                                                   perf.ci.lower = best.perf$auc.CI[1], perf.ci.upper = best.perf$auc.CI[2], blurb=blurb)
      }else{
        cat(sprintf("\nFinal %s=%s\n%s upgenes and %s downgenes chosen\n",perf.meas,curr.max,length(pos),length(neg)))
        continue = FALSE
      }
      if(continue && count >= max.genes){
        cat(sprintf("\nFinal %s=%s\n%s upgenes and %s downgenes chosen\n",perf.meas,curr.max,length(pos),length(neg)))
      }
    }
  }
  
  if(perf.meas=="AUPRC"){
    while(continue && (count < max.genes)){
      auprclist=mclapply(mc.cores=numCores, genelist, function(gene){
        currpos=pos
        currneg=neg
        if(gene %in% pos.genes){
          currpos=append(currpos,gene)
        }else{
          currneg=append(currneg,gene)
        }
        scores = getGeneScores(geneMtx,currpos,currneg)
        pred_PRC = prediction(as.numeric(scores), as.numeric(as.character(class)))
        perf = performance(pred_PRC, "prec", "rec")
        auprc=.calcauprc(perf@x.values[[1]],perf@y.values[[1]],class)$auprc
        
        #hopefully this works
        if(force.posneg && count == 1){
          if(is.null(pos) && gene %in% neg.genes){auprc=0}
          if(is.null(neg) && gene %in% pos.genes){auprc=0}
        }
        return(auprc)
      })
      #check for NULL values
      if(any(sapply(auprclist,function(x) is.null(x)))){
        warning("Some values of auprclist are NULL - this is due to an mclapply error. Reduce numCores and try again.")
      }
      auprclist = unlist(auprclist)
      curr.max = max(auprclist,na.rm=T)
      perf.diff = curr.max-perf.max
      cat(sprintf("next best: %s\n",perf.diff))
      count = count+1
      if(perf.diff>forwardThresh){
        perf.max = curr.max
        topindex = which.max(auprclist) #not 100% sure how this will behave if there is a tie
        if(length(topindex)>1){topindex=topindex[1]} #just make sure topindex is only one index
        topgene=genelist[topindex]
        if(topgene %in% pos.genes){
          cat(sprintf("Adding %s (up)\n",topgene))
          pos=append(pos,topgene)
        }else{
          cat(sprintf("Adding %s (down)\n",topgene))
          neg=append(neg,topgene)
        }
        if(!replace.genes){
          genelist=genelist[-topindex]
        }
        best.scores = getGeneScores(geneMtx,pos,neg)
        best.pred_PRC = prediction(as.numeric(best.scores), as.numeric(as.character(class)))
        best.perf_PRC = performance(best.pred_PRC, "prec", "rec")
        best.perf = .calcauprc(best.perf_PRC@x.values[[1]],best.perf_PRC@y.values[[1]],class)
        blurb = sprintf("%s Genes, %s=%s (95%% CI %s-%s): %s up and %s down",count,perf.meas,round(best.perf$auprc,3),round(best.perf$auprc.CI[1],3),
                        round(best.perf$auprc.CI[2],3),paste(pos,collapse=", "),paste(neg,collapse=", "))
        bestgenelist[[as.character(count)]] = list(numGenes = count, upgenes = pos, downgenes = neg, perf.meas = perf.meas, perf = best.perf$auprc,
                                                   perf.ci.lower = best.perf$auprc.CI[1], perf.ci.upper = best.perf$auprc.CI[2], blurb=blurb)
      }else{
        cat(sprintf("\nFinal %s=%s\n%s upgenes and %s downgenes chosen\n",perf.meas,curr.max,length(pos),length(neg)))
        continue = FALSE
      }
      if(continue && count >= max.genes){
        cat(sprintf("\nFinal %s=%s\n%s upgenes and %s downgenes chosen\n",perf.meas,curr.max,length(pos),length(neg)))
      }
    }
  }
  
  if(perf.meas=="AvgAUC"){
    labelVec = factor(labelVec)
    labelVec=droplevels(labelVec)
    if(is.null(caseNames)){
      caseNames=as.character(unique(labelVec[class==1]))
    }
    if(is.null(otherNames)){
      otherNames=sort(levels(labelVec)[-c(which(levels(labelVec) %in% caseNames))])
    }
    classIndex = lapply(otherNames, function(other){
      return(append(which(labelVec==other),which(labelVec %in% caseNames)))
    })
    if(is.null(weights)){
      weights = rep(1,length(classIndex))
    }
    if(length(weights) != length(classIndex)){stop("weights must be the same length as the number of class comparisons")}
    weights = weights/sum(weights,na.rm=TRUE) #normalize to one
    
    while(continue && (count < max.genes)){
      auclist=mclapply(mc.cores=numCores, genelist, function(gene){
        currpos=pos
        currneg=neg
        if(gene %in% pos.genes){
          currpos=append(currpos,gene)
        }else{
          currneg=append(currneg,gene)
        }
        scores = getGeneScores(geneMtx,currpos,currneg)
        aucVec=rep(0,length(classIndex))
        for(i in 1:length(classIndex)){
          aucVec[i] = as.numeric(calculateROC(as.numeric(as.character(class[classIndex[[i]] ])), as.numeric(scores[classIndex[[i]] ]))$auc)
        }
        avgauc=sum(aucVec*weights,na.rm = TRUE)
        
        #hopefully this works
        if(force.posneg && count == 1){
          if(is.null(pos) && gene %in% neg.genes){avgauc=0}
          if(is.null(neg) && gene %in% pos.genes){avgauc=0}
        }
        return(avgauc)
      })
      #check for NULL values
      if(any(sapply(auclist,function(x) is.null(x)))){
        warning("Some values of auclist are NULL - this is due to an mclapply error. Reduce numCores and try again.")
      }
      auclist = unlist(auclist)
      curr.max = max(auclist,na.rm=T)
      perf.diff = curr.max-perf.max
      cat(sprintf("next best: %s\n",perf.diff))
      count = count+1
      if(perf.diff>forwardThresh){
        perf.max = curr.max
        topindex = which.max(auclist) #not 100% sure how this will behave if there is a tie
        if(length(topindex)>1){topindex=topindex[1]} #just make sure topindex is only one index
        topgene=genelist[topindex]
        if(topgene %in% pos.genes){
          cat(sprintf("Adding %s (up)\n",topgene))
          pos=append(pos,topgene)
        }else{
          cat(sprintf("Adding %s (down)\n",topgene))
          neg=append(neg,topgene)
        }
        if(!replace.genes){
          genelist=genelist[-topindex]
        }
        best.scores = getGeneScores(geneMtx,pos,neg)
        best.aucVec=best.aucVecLo=best.aucVecHi=rep(0,length(classIndex))
        for(i in 1:length(classIndex)){
          best.perf = calculateROC(as.numeric(as.character(class[classIndex[[i]] ])), as.numeric(best.scores[classIndex[[i]] ]))
          best.aucVec[i] = as.numeric(best.perf$auc)
          best.aucVecLo[i] = best.perf$auc.CI[1]
          best.aucVecHi[i] = best.perf$auc.CI[2]
        }
        best.avgauc=sum(best.aucVec*weights,na.rm = TRUE)
        best.avgauclo=sum(best.aucVecLo*weights,na.rm = TRUE)
        best.avgauchi=sum(best.aucVecHi*weights,na.rm = TRUE)
        blurb = sprintf("%s Genes, %s=%s (95%% CI %s-%s): %s up and %s down",count,perf.meas,round(best.avgauc,3),round(best.avgauclo,3),
                        round(best.avgauchi,3),paste(pos,collapse=", "),paste(neg,collapse=", "))
        allAUCs=numeric(length(best.aucVec)*3)
        allAUCs[seq(1,by=3,length.out=length(best.aucVec))]=best.aucVec
        allAUCs[seq(2,by=3,length.out=length(best.aucVecLo))]=best.aucVecLo
        allAUCs[seq(3,by=3,length.out=length(best.aucVecHi))]=best.aucVecHi
        bestgenelist[[as.character(count)]] = list(numGenes = count, upgenes = pos, downgenes = neg, perf.meas = perf.meas, perf = best.avgauc,
                                                   perf.ci.lower = best.avgauclo, perf.ci.upper = best.avgauchi, blurb=blurb, allAUCs=allAUCs)
      }else{
        cat(sprintf("\nFinal %s=%s\n%s upgenes and %s downgenes chosen\n",perf.meas,curr.max,length(pos),length(neg)))
        continue = FALSE
      }
      if(continue && count >= max.genes){
        cat(sprintf("\nFinal %s=%s\n%s upgenes and %s downgenes chosen\n",perf.meas,curr.max,length(pos),length(neg)))
      }
    }
  }
  
  if(perf.meas=="Specificity"){
    while(continue && (count < max.genes)){
      speclist=mclapply(mc.cores=numCores, genelist, function(gene){
        currpos=pos
        currneg=neg
        if(gene %in% pos.genes){
          currpos=append(currpos,gene)
        }else{
          currneg=append(currneg,gene)
        }
        scores = getGeneScores(geneMtx,currpos,currneg)
        specTab = calcSpecs(scores, class, sensVec)
        specMean = mean(specTab$Spec)
        
        #hopefully this works
        if(force.posneg && count == 1){
          if(is.null(pos) && gene %in% neg.genes){specMean=0}
          if(is.null(neg) && gene %in% pos.genes){specMean=0}
        }
        return(specMean)
      })
      #check for NULL values
      if(any(sapply(speclist,function(x) is.null(x)))){
        warning("Some values of speclist are NULL - this is due to an mclapply error. Reduce numCores and try again.")
      }
      speclist = unlist(speclist)
      curr.max = max(speclist,na.rm=T)
      perf.diff = curr.max-perf.max
      cat(sprintf("next best: %s\n",perf.diff))
      count = count+1
      if(perf.diff>forwardThresh){
        perf.max = curr.max
        topindex = which.max(speclist) #not 100% sure how this will behave if there is a tie
        if(length(topindex)>1){topindex=topindex[1]} #just make sure topindex is only one index
        topgene=genelist[topindex]
        if(topgene %in% pos.genes){
          cat(sprintf("Adding %s (up)\n",topgene))
          pos=append(pos,topgene)
        }else{
          cat(sprintf("Adding %s (down)\n",topgene))
          neg=append(neg,topgene)
        }
        if(!replace.genes){
          genelist=genelist[-topindex]
        }
        best.scores = getGeneScores(geneMtx,pos,neg)
        best.specMean = mean(calcSpecs(best.scores,class,sensVec)$Spec)
        blurb = sprintf("%s Genes, Average %s=%s: %s up and %s down",count,perf.meas,round(best.specMean,3),paste(pos,collapse=", "),paste(neg,collapse=", "))
        bestgenelist[[as.character(count)]] = list(numGenes = count, upgenes = pos, downgenes = neg, perf.meas = perf.meas, perf = best.specMean,
                                                   perf.ci.lower = best.specMean, perf.ci.upper = best.specMean, blurb=blurb)
      }else{
        cat(sprintf("\nFinal Average %s=%s\n%s upgenes and %s downgenes chosen\n",perf.meas,curr.max,length(pos),length(neg)))
        continue = FALSE
      }
      if(continue && count >= max.genes){
        cat(sprintf("\nFinal Average %s=%s\n%s upgenes and %s downgenes chosen\n",perf.meas,curr.max,length(pos),length(neg)))
      }
    }
  }
  
  
  bestgeneDT = data.table(do.call(rbind,lapply(bestgenelist,function(x){
    data.frame(numGenes = x$numGenes, perf = x$perf, perf.ci.lower = x$perf.ci.lower, perf.ci.upper = x$perf.ci.upper,
               upgenes = paste(x$upgenes,collapse = " / "), downgenes = paste(x$downgenes,collapse = " / "), blurb=x$blurb,
               stringsAsFactors = FALSE)
  })))
  names(bestgeneDT)[2] = perf.meas
  names(bestgeneDT)[3] = paste0(perf.meas,".ci.lower")
  names(bestgeneDT)[4] = paste0(perf.meas,".ci.upper")
  
  if(perf.meas=="AvgAUC"){
    allAUCDT = data.table(do.call(rbind,lapply(bestgenelist,function(x){
      x$allAUCs
    })))
    otherLabs=rep("",length(classIndex)*3)
    for(i in 1:length(classIndex)){
      otherLabs[(3*(i-1))+1]=paste(otherNames[i],"auc",sep=".")
      otherLabs[(3*(i-1))+2]=paste(otherNames[i],"ci.lower",sep=".")
      otherLabs[(3*(i-1))+3]=paste(otherNames[i],"ci.upper",sep=".")
    }
    names(allAUCDT)=otherLabs
    bestgeneDT = cbind(bestgeneDT,allAUCDT)
  }
  if(perf.meas=="Specificity"){
    bestgeneDT$Specificity.ci.lower = NULL
    bestgeneDT$Specificity.ci.upper = NULL
    names(bestgeneDT)[2] = "AverageSpecificity"
    #could maybe adjust this to report the individual specificities as well? but that's low priority tbh
  }
  
  return(list(upgenes=pos, downgenes=neg, bestgenelist=bestgenelist, bestgeneDT=bestgeneDT))
}



#' Helper Function for Backward Search
#' @param pos.genes Vector of upgenes
#' @param neg.genes Vector of downgenes
#' @param perf.meas Which diagnostic performance metric to use. Options are "AUC", "AUPRC", and "AvgAUC" (Default: "AUC")
#' @param backwardThresh Stopping threshold. As long as the performance of the best N genes is better than the performance of the best N+1 genes
#' plus the backwardThresh, then the backward search will continue. For example, if backwardThresh is -0.01, 
#' then the performance of the best N genes can only decrease by 0.01 or less each time in order to continue (Default: 0)
#' @param geneMtx The gene matrix to use for calculating signature scores (rows are genes and columns are samples)
#' @param class Vector of class labels
#' @param labelVec \emph{Only for perf.meas=AvgAUC:} Vector of sample category labels
#' @param weights \emph{Only for perf.meas=AvgAUC:} If left blank, then the simple average of class-specific AUCs will be used to
#' calculate the AvgAUC. However, if you want to take a weighted average instead, then pass in a vector of weights. Note that the
#' order of the weights will be matched to the order of the names in the \code{otherNames} vector.
#' @param caseNames \emph{Only for perf.meas=AvgAUC:} name of the class(es) that you're considering to be your case.
#' If left blank, then \code{caseNames} will be the unique set of names in \code{$pheno$group} that correspond to a class label of 1
#' @param otherNames \emph{Only for perf.meas=AvgAUC:} name of the class(es) that you're considering to be your controls.
#' If left blank, then \code{otherNames} will be the unique set of names in \code{$pheno$group} that correspond to a class label of 0.
#' The order of otherNames will be set by running the \code{sort} function, so if you are using weights, make sure that your \code{weights}
#' vector and \code{otherNames} vector match up correctly
#' @param min.genes If you don't want to look beyond a certain number of genes, then set \code{min.genes}
#' as the minimum number of genes to include in the resulting signature
#' @param numCores Number of cores to use for parallelization. If left blank, smart core usage will be used
#' @note As of now, the accepted diagnostic performance metrics are AUC (area under the Reciever Operating Curve),
#' AUPRC (area under the Precision Recall Curve), and AvgAUC
#' (the average AUC across all comparisons between your cases and each of your "control" categories)
#' @return
#'   \item{upgenes}{Vector of the upgenes that were chosen by the backward search}
#'   \item{downgenes}{Vector of the downgenes that were chosen by the backward search}
#' @import parallel zoo pracma MetaIntegrator ROCR
#' @author Aditya Rao
.backwardSearchHelper <- function(pos.genes, neg.genes, perf.meas = "AUC", backwardThresh = 0, geneMtx, class, labelVec, weights=NULL,
                                  caseNames=NULL, otherNames=NULL, min.genes=NULL, numCores=NULL){
  perf.meas = tolower(perf.meas)
  perf.meas = switch(perf.meas, "auc"="AUC", "auprc"="AUPRC", "avgauc"="AvgAUC")
  if(!(perf.meas %in% c("AUC","AUPRC","AvgAUC"))){stop("perf.meas must be \"AUC\", \"AUPRC\", or \"AvgAUC\"")}
  if(is.null(numCores)){
    #smart core usage = use 50% of resources unless you need less
    numCores <- round(detectCores()*5/10)
    if(length(c(pos.genes,neg.genes))<numCores){
      numCores <- length(c(pos.genes,neg.genes))
    }
  }
  if(any(!c(pos.genes,neg.genes) %in% rownames(geneMtx))){stop("Some of the designated genes are not present in the provided gene matrix")}
  pos=pos.genes
  neg=neg.genes
  count=length(c(pos,neg))
  perf.max = 0
  continue=TRUE
  genelist=c(pos.genes,neg.genes)
  bestgenelist = list()
  if(is.null(min.genes)){
    min.genes = 1
  }else if(min.genes >= (length(pos.genes) + length(neg.genes))){
    warning("The provided min.genes is greater than or equal to the total number of genes and will be ignored")
    min.genes = 1
  }
  
  if(perf.meas=="AUC"){
    #get bestgenelist entry for set of all genes
    scores.all = getGeneScores(geneMtx,pos,neg)
    perf.all = calculateROC(as.numeric(as.character(class)), as.numeric(scores.all))
    auc.all = as.numeric(perf.all$auc)
    blurb.all = sprintf("%s Genes, %s=%s (95%% CI %s-%s): %s up and %s down",count,perf.meas,round(auc.all,3),round(perf.all$auc.CI[1],3),
                        round(perf.all$auc.CI[2],3),paste(pos,collapse=", "),paste(neg,collapse=", "))
    bestgenelist[[as.character(count)]] = list(numGenes = count, upgenes = pos, downgenes = neg, perf.meas = perf.meas, perf = auc.all,
                                               perf.ci.lower = perf.all$auc.CI[1], perf.ci.upper = perf.all$auc.CI[2], blurb=blurb.all)
    
    while(continue && (count > min.genes)){
      auclist=mclapply(mc.cores=numCores, genelist, function(gene){
        currpos=pos
        currneg=neg
        if(gene %in% pos.genes){
          currpos=currpos[-which(currpos==gene)]
        }else{
          currneg=currneg[-which(currneg==gene)]
        }
        scores = getGeneScores(geneMtx,currpos,currneg)
        auc = as.numeric(calculateROC(as.numeric(as.character(class)), as.numeric(scores))$auc)
        return(auc)
      })
      #check for NULL values
      if(any(sapply(auclist,function(x) is.null(x)))){
        warning("Some values of auclist are NULL - this is due to an mclapply error. Reducing numCores and try again.")
      }
      auclist = unlist(auclist)
      curr.max = max(auclist,na.rm=T)
      perf.diff = curr.max-perf.max
      cat(sprintf("next best: %s\n",perf.diff))
      count = count-1
      if(perf.diff>backwardThresh){
        perf.max = curr.max
        botindex = which.max(auclist) #not 100% sure how this will behave if there is a tie
        if(length(botindex)>1){botindex=botindex[1]} #just make sure topindex is only one index
        botgene=genelist[botindex]
        if(botgene %in% pos.genes){
          cat(sprintf("Removing %s (up)\n",botgene))
          pos=pos[-which(pos==botgene)]
        }else{
          cat(sprintf("Removing %s (down)\n",botgene))
          neg=neg[-which(neg==botgene)]
        }
        genelist=genelist[-botindex]
        best.scores = getGeneScores(geneMtx,pos,neg)
        best.perf = calculateROC(as.numeric(as.character(class)), as.numeric(best.scores))
        best.auc = as.numeric(best.perf$auc)
        blurb = sprintf("%s Genes, %s=%s (95%% CI %s-%s): %s up and %s down",count,perf.meas,round(best.auc,3),round(best.perf$auc.CI[1],3),
                        round(best.perf$auc.CI[2],3),paste(pos,collapse=", "),paste(neg,collapse=", "))
        bestgenelist[[as.character(count)]] = list(numGenes = count, upgenes = pos, downgenes = neg, perf.meas = perf.meas, perf = best.auc,
                                                   perf.ci.lower = best.perf$auc.CI[1], perf.ci.upper = best.perf$auc.CI[2], blurb=blurb)
      }else{
        cat(sprintf("\nFinal %s=%s\n%s upgenes and %s downgenes chosen\n",perf.meas,curr.max,length(pos),length(neg)))
        continue = FALSE
      }
      if(continue && count <= min.genes){
        cat(sprintf("\nFinal %s=%s\n%s upgenes and %s downgenes chosen\n",perf.meas,curr.max,length(pos),length(neg)))
      }
    }
  }
  
  if(perf.meas=="AUPRC"){
    #get bestgenelist entry for set of all genes
    scores.all = getGeneScores(geneMtx,pos,neg)
    pred_PRC.all = prediction(as.numeric(scores.all), as.numeric(as.character(class)))
    perf_PRC.all = performance(pred_PRC.all, "prec", "rec")
    perf.all = .calcauprc(perf_PRC.all@x.values[[1]],perf_PRC.all@y.values[[1]],class)
    blurb.all = sprintf("%s Genes, %s=%s (95%% CI %s-%s): %s up and %s down",count,perf.meas,round(perf.all$auprc,3),round(perf.all$auprc.CI[1],3),
                    round(perf.all$auprc.CI[2],3),paste(pos,collapse=", "),paste(neg,collapse=", "))
    bestgenelist[[as.character(count)]] = list(numGenes = count, upgenes = pos, downgenes = neg, perf.meas = perf.meas, perf = perf.all$auprc,
                                               perf.ci.lower = perf.all$auprc.CI[1], perf.ci.upper = perf.all$auprc.CI[2], blurb=blurb.all)
    
    while(continue && (count > min.genes)){
      auprclist=mclapply(mc.cores=numCores, genelist, function(gene){
        currpos=pos
        currneg=neg
        if(gene %in% pos.genes){
          currpos=currpos[-which(currpos==gene)]
        }else{
          currneg=currneg[-which(currneg==gene)]
        }
        scores = getGeneScores(geneMtx,currpos,currneg)
        pred_PRC = prediction(as.numeric(scores), as.numeric(as.character(class)))
        perf = performance(pred_PRC, "prec", "rec")
        auprc=.calcauprc(perf@x.values[[1]],perf@y.values[[1]],class)$auprc
        return(auprc)
      })
      #check for NULL values
      if(any(sapply(auprclist,function(x) is.null(x)))){
        warning("Some values of auprclist are NULL - this is due to an mclapply error. Reducing numCores and try again.")
      }
      auprclist = unlist(auprclist)
      curr.max = max(auprclist,na.rm=T)
      perf.diff = curr.max-perf.max
      cat(sprintf("next best: %s\n",perf.diff))
      count = count-1
      if(perf.diff>backwardThresh){
        perf.max = curr.max
        botindex = which.max(auprclist) #not 100% sure how this will behave if there is a tie
        if(length(botindex)>1){botindex=botindex[1]} #just make sure topindex is only one index
        botgene=genelist[botindex]
        if(botgene %in% pos.genes){
          cat(sprintf("Removing %s (up)\n",botgene))
          pos=pos[-which(pos==botgene)]
        }else{
          cat(sprintf("Removing %s (down)\n",botgene))
          neg=neg[-which(neg==botgene)]
        }
        genelist=genelist[-botindex]
        best.scores = getGeneScores(geneMtx,pos,neg)
        best.pred_PRC = prediction(as.numeric(best.scores), as.numeric(as.character(class)))
        best.perf_PRC = performance(best.pred_PRC, "prec", "rec")
        best.perf = .calcauprc(best.perf_PRC@x.values[[1]],best.perf_PRC@y.values[[1]],class)
        blurb = sprintf("%s Genes, %s=%s (95%% CI %s-%s): %s up and %s down",count,perf.meas,round(best.perf$auprc,3),round(best.perf$auprc.CI[1],3),
                        round(best.perf$auprc.CI[2],3),paste(pos,collapse=", "),paste(neg,collapse=", "))
        bestgenelist[[as.character(count)]] = list(numGenes = count, upgenes = pos, downgenes = neg, perf.meas = perf.meas, perf = best.perf$auprc,
                                                   perf.ci.lower = best.perf$auprc.CI[1], perf.ci.upper = best.perf$auprc.CI[2], blurb=blurb)
      }else{
        cat(sprintf("\nFinal %s=%s\n%s upgenes and %s downgenes chosen\n",perf.meas,curr.max,length(pos),length(neg)))
        continue = FALSE
      }
      if(continue && count <= min.genes){
        cat(sprintf("\nFinal %s=%s\n%s upgenes and %s downgenes chosen\n",perf.meas,curr.max,length(pos),length(neg)))
      }
    }
  }
  
  if(perf.meas=="AvgAUC"){
    labelVec = factor(labelVec)
    labelVec=droplevels(labelVec)
    if(is.null(caseNames)){
      caseNames=as.character(unique(labelVec[class==1]))
    }
    if(is.null(otherNames)){
      otherNames=levels(labelVec)[-c(which(levels(labelVec) %in% caseNames))]
    }
    classIndex = lapply(otherNames, function(other){
      return(append(which(labelVec==other),which(labelVec %in% caseNames)))
    })
    if(is.null(weights)){
      weights = rep(1,length(classIndex))
    }
    weights = weights/sum(weights,na.rm=TRUE) #normalize to one
    
    #get bestgenelist entry for set of all genes
    scores.all = getGeneScores(geneMtx,pos,neg)
    aucVec.all=aucVecLo.all=aucVecHi.all=rep(0,length(classIndex))
    for(i in 1:length(classIndex)){
      perf.all = calculateROC(as.numeric(as.character(class[classIndex[[i]] ])), as.numeric(scores.all[classIndex[[i]] ]))
      aucVec.all[i] = as.numeric(perf.all$auc)
      aucVecLo.all[i] = perf.all$auc.CI[1]
      aucVecHi.all[i] = perf.all$auc.CI[2]
    }
    avgauc.all=sum(aucVec.all*weights,na.rm = TRUE)
    avgauclo.all=sum(aucVecLo.all*weights,na.rm = TRUE)
    avgauchi.all=sum(aucVecHi.all*weights,na.rm = TRUE)
    blurb.all = sprintf("%s Genes, %s=%s (95%% CI %s-%s): %s up and %s down",count,perf.meas,round(avgauc.all,3),round(avgauclo.all,3),
                        round(avgauchi.all,3),paste(pos,collapse=", "),paste(neg,collapse=", "))
    allAUCs.all=numeric(length(aucVec.all)*3)
    allAUCs.all[seq(1,by=3,length.out=length(aucVec.all))]=aucVec.all
    allAUCs.all[seq(2,by=3,length.out=length(aucVecLo.all))]=aucVecLo.all
    allAUCs.all[seq(3,by=3,length.out=length(aucVecHi.all))]=aucVecHi.all
    bestgenelist[[as.character(count)]] = list(numGenes = count, upgenes = pos, downgenes = neg, perf.meas = perf.meas, perf = avgauc.all,
                                               perf.ci.lower = avgauclo.all, perf.ci.upper = avgauchi.all, blurb=blurb.all, allAUCs=allAUCs.all)
    
    while(continue && (count > min.genes)){
      auclist=mclapply(mc.cores=numCores, genelist, function(gene){
        currpos=pos
        currneg=neg
        if(gene %in% pos.genes){
          currpos=currpos[-which(currpos==gene)]
        }else{
          currneg=currneg[-which(currneg==gene)]
        }
        scores = getGeneScores(geneMtx,currpos,currneg)
        aucVec=rep(0,length(classIndex))
        for(i in 1:length(classIndex)){
          aucVec[i] = as.numeric(calculateROC(as.numeric(as.character(class[classIndex[[i]] ])), as.numeric(scores[classIndex[[i]] ]))$auc)
        }
        avgauc=sum(aucVec*weights,na.rm = TRUE)
        return(avgauc)
      })
      #check for NULL values
      if(any(sapply(auclist,function(x) is.null(x)))){
        warning("Some values of auclist are NULL - this is due to an mclapply error. Reducing numCores and try again.")
      }
      auclist = unlist(auclist)
      curr.max = max(auclist,na.rm=T)
      perf.diff = curr.max-perf.max
      cat(sprintf("next best: %s\n",perf.diff))
      count = count-1
      if(perf.diff>backwardThresh){
        perf.max = curr.max
        botindex = which.max(auclist) #not 100% sure how this will behave if there is a tie
        if(length(botindex)>1){botindex=botindex[1]} #just make sure topindex is only one index
        botgene=genelist[botindex]
        if(botgene %in% pos.genes){
          cat(sprintf("Removing %s (up)\n",botgene))
          pos=pos[-which(pos==botgene)]
        }else{
          cat(sprintf("Removing %s (down)\n",botgene))
          neg=neg[-which(neg==botgene)]
        }
        genelist=genelist[-botindex]
        best.scores = getGeneScores(geneMtx,pos,neg)
        best.aucVec=best.aucVecLo=best.aucVecHi=rep(0,length(classIndex))
        for(i in 1:length(classIndex)){
          best.perf = calculateROC(as.numeric(as.character(class[classIndex[[i]] ])), as.numeric(best.scores[classIndex[[i]] ]))
          best.aucVec[i] = as.numeric(best.perf$auc)
          best.aucVecLo[i] = best.perf$auc.CI[1]
          best.aucVecHi[i] = best.perf$auc.CI[2]
        }
        best.avgauc=sum(best.aucVec*weights,na.rm = TRUE)
        best.avgauclo=sum(best.aucVecLo*weights,na.rm = TRUE)
        best.avgauchi=sum(best.aucVecHi*weights,na.rm = TRUE)
        blurb = sprintf("%s Genes, %s=%s (95%% CI %s-%s): %s up and %s down",count,perf.meas,round(best.avgauc,3),round(best.avgauclo,3),
                        round(best.avgauchi,3),paste(pos,collapse=", "),paste(neg,collapse=", "))
        allAUCs=numeric(length(best.aucVec)*3)
        allAUCs[seq(1,by=3,length.out=length(best.aucVec))]=best.aucVec
        allAUCs[seq(2,by=3,length.out=length(best.aucVecLo))]=best.aucVecLo
        allAUCs[seq(3,by=3,length.out=length(best.aucVecHi))]=best.aucVecHi
        bestgenelist[[as.character(count)]] = list(numGenes = count, upgenes = pos, downgenes = neg, perf.meas = perf.meas, perf = best.avgauc,
                                                   perf.ci.lower = best.avgauclo, perf.ci.upper = best.avgauchi, blurb=blurb, allAUCs=allAUCs)
      }else{
        cat(sprintf("\nFinal %s=%s\n%s upgenes and %s downgenes chosen\n",perf.meas,curr.max,length(pos),length(neg)))
        continue = FALSE
      }
      if(continue && count <= min.genes){
        cat(sprintf("\nFinal %s=%s\n%s upgenes and %s downgenes chosen\n",perf.meas,curr.max,length(pos),length(neg)))
      }
    }
  }
  
  bestgeneDT = data.table(do.call(rbind,lapply(bestgenelist,function(x){
    data.frame(numGenes = x$numGenes, perf = x$perf, perf.ci.lower = x$perf.ci.lower, perf.ci.upper = x$perf.ci.upper,
               upgenes = paste(x$upgenes,collapse = " / "), downgenes = paste(x$downgenes,collapse = " / "), blurb=x$blurb,
               stringsAsFactors = FALSE)
  })))
  names(bestgeneDT)[2] = perf.meas
  names(bestgeneDT)[3] = paste0(perf.meas,".ci.lower")
  names(bestgeneDT)[4] = paste0(perf.meas,".ci.upper")
  
  if(perf.meas=="AvgAUC"){
    allAUCDT = data.table(do.call(rbind,lapply(bestgenelist,function(x){
      x$allAUCs
    })))
    otherLabs=rep("",length(classIndex)*3)
    for(i in 1:length(classIndex)){
      otherLabs[(3*(i-1))+1]=paste(otherNames[i],"auc",sep=".")
      otherLabs[(3*(i-1))+2]=paste(otherNames[i],"ci.lower",sep=".")
      otherLabs[(3*(i-1))+3]=paste(otherNames[i],"ci.upper",sep=".")
    }
    names(allAUCDT)=otherLabs
    bestgeneDT = cbind(bestgeneDT,allAUCDT)
  }
  
  return(list(upgenes=pos, downgenes=neg, bestgenelist=bestgenelist, bestgeneDT=bestgeneDT))
}









#NOT ADDING COMMENTS FOR THIS FUNCTION CUZ I'M VERY SLEEPY
calcSpecs <- function(scores, class, sensitivities = c(0.9,0.925,0.95,0.975,0.98,0.99,0.995,0.999), add_thresh=FALSE){
  df = data.frame(score = scores, class = class)
  optList = lapply(sensitivities, function(sens){
    return(optimal.cutpoints(X="score", status="class", tag.healthy=0,
                             data=df,methods=c("MinValueSe"),
                             control=control.cutpoints(valueSe=sens)))
  })
  
  resultsTab = data.frame(t(sapply(optList, function(opt){
    c(opt$MinValueSe$Global$optimal.cutoff$Se,opt$MinValueSe$Global$optimal.cutoff$Sp)
  })))
  
  if(add_thresh){
    resultsTab$Thresh = sapply(optList, function(opt) opt$MinValueSe$Global$optimal.cutoff$cutoff)
  }
  
  colnames(resultsTab)[1:2] = c("Sens","Spec")
  rownames(resultsTab) = paste0("Sens>",sensitivities*100,"%")
  return(resultsTab)
}





