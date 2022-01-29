#' Best Subset Selection Function
#' @description
#' Conducts best subset selection on a set of genes and returns a ranked list of
#' all gene combinations, separated by model size (i.e. number of genes). Ranking can be
#' done by AUC, AUPRC, or Average AUC. \cr
#' 
#' Best subset selection is a method of optimizing a given set of significant
#' genes to maximize some metric of diagnostic performance. The function works
#' by taking a given set of genes (presumably a set that has been filtered for
#' statistical significance), and testing the performance of every single
#' possible combination of those genes. At the end, the best performing gene
#' combinations are selected.
#' 
#' \strong{WARNING:} the time this function takes increases exponentially, so be very
#' careful with it. If running more than ~15 genes, it will likely require a multi-core server.
#' @param pooledDataObject A MANATEE pooledDataObject
#' @param filterObject A MANATEE filterObject
#' @param perf.meas Which diagnostic performance metric to use. Options are "AUC", "AUPRC", "AvgAUC", or "All" (Default: "AUC")
#' @param upgenes If you want to manually input genes instead of providing a filterObject, this is a vector of the upgenes in your signature
#' @param downgenes If you want to manually input genes instead of providing a filterObject, this is a vector of the downgenes in your signature
#' @param bss.min Minimum model size. For example, if \code{bss.min} is 1 and \code{bss.max} is 5, then all of the combinations of 1, 2, 3, 4, and 5 genes will be tested (Default: 1)
#' @param bss.max Maximum model size. For example, if \code{bss.min} is 1 and \code{bss.max} is 5, then all of the combinations of 1, 2, 3, 4, and 5 genes will be tested.
#' If left blank, then it will be set as the combined length of your upgenes and downgenes.
#' @param weights \emph{Only for perf.meas=AvgAUC or perf.meas=All:} If left blank, then the simple average of class-specific AUCs will be used to
#' calculate the AvgAUC. However, if you want to take a weighted average instead, then pass in a vector of weights. Note that the
#' order of the weights will be matched to the order of the names in the \code{otherNames} vector.
#' @param caseNames \emph{Only for perf.meas=AvgAUC or perf.meas=All:} name of the class(es) that you're considering to be your case.
#' If left blank, then \code{caseNames} will be the unique set of names in \code{$pheno$group} that correspond to a class label of 1
#' @param otherNames \emph{Only for perf.meas=AvgAUC or perf.meas=All:} name of the class(es) that you're considering to be your controls.
#' If left blank, then \code{otherNames} will be the unique set of names in \code{$pheno$group} that correspond to a class label of 0.
#' The order of otherNames will be set by running the \code{sort} function, so if you are using weights, make sure that your \code{weights}
#' vector and \code{otherNames} vector match up correctly
#' @param group.colname \emph{Only for perf.meas=AvgAUC or perf.meas=All:} this designates the column name of the group vector within the pooledDataObject$pheno
#' dataframe (Default: "group")
#' @param init.pos If you want to start off with some upgenes already locked in, then indicate those genes with \code{init.pos}
#' @param init.neg If you want to start off with some downgenes genes already locked in, then indicate those genes with \code{init.pos}
#' @param sort.by \emph{Only for perf.meas=All:} Which performance measure should be used to rank the genes (Default: "AUC")
#' @param numCores Number of cores to use for parallelization. If left blank, smart core usage will be used
#' @details
#' The bestSubsetSelection function is designed to assist in selection of gene
#' sets optimized for discriminatory power. Unlike a greedy search (like forward
#' search or backward search), best subset selection will yield gene sets that
#' are only globally optimized (ie, they are not merely local optima). \cr
#' @note As of now, the accepted diagnostic performance metrics are AUC (area under the Reciever Operating Curve),
#' AUPRC (area under the Precision Recall Curve), AvgAUC
#' (the average AUC across all comparisons between your cases and each of your "control" categories),
#' and All (all three metrics are calculated and then sorted according to the \code{sort.by} value)
#' @return Returns a modified version of the \code{pooledDataObject} with the best subset selection results
#' stored in \code{pooledDataObject$bssResults} under a unique identifier
#'   \item{\emph{X}.Genes}{The ranked results for each model size is stored in a separate data frame that contains information about
#'   performance, the included genes, and a short blurb for each gene combination tested}
#'   \item{bssDT}{Data table containing the combined results of the entire best subset selection}
#'   \item{bssDescription}{List with information that describes the input parameters for the best subset selection}
#' @examples
#' 
#' @export
#' @import parallel zoo pracma MetaIntegrator ROCR data.table pbapply
#' @importFrom pbmcapply pbmclapply
#' @author Aditya Rao
bestSubsetSelection <- function(pooledDataObject, filterObject=NULL, perf.meas="AUC", upgenes, downgenes, bss.min=1, bss.max=NULL,
                                weights=NULL, caseNames=NULL, otherNames=NULL, group.colname="group", init.pos=NULL, init.neg=NULL, sort.by="AUC", numCores=NULL){
  progress.bar.type="perGeneNum"
  if(!checkManateeObject(pooledDataObject,"pooledDataObject")){stop("Invalid pooledDataObject provided")}
  if(!is.null(filterObject)){
    upgenes = filterObject$upGeneNames
    downgenes = filterObject$downGeneNames
  }
  perf.meas = tolower(perf.meas)
  perf.meas = switch(perf.meas, "auc"="AUC", "auprc"="AUPRC", "avgauc"="AvgAUC", "all"="All")
  if(!(perf.meas %in% c("AUC","AUPRC","AvgAUC","All"))){stop("perf.meas must be \"AUC\", \"AUPRC\", \"AvgAUC\", or \"All\"")}
  
  BestSubsetList = .bestSubsetSelectionHelper(perf.meas=perf.meas, geneMtx=pooledDataObject$genes, class=pooledDataObject$class, pos.genes=upgenes, neg.genes=downgenes, bss.min=bss.min, bss.max=bss.max,
                                              labelVec=pooledDataObject$pheno[,group.colname], weights=weights, caseNames=caseNames, otherNames=otherNames, init.pos=init.pos, init.neg=init.neg, sort.by=sort.by, numCores=numCores, progress.bar.type=progress.bar.type)
  
  bss.name = sprintf("%s.%sUp%sDown",perf.meas,length(upgenes),length(downgenes))
  if("bssResults" %in% names(pooledDataObject)){
    if(bss.name %in% names(pooledDataObject$bssResults)){
      continue=TRUE
      iter = 1
      while(continue){
        temp.name = paste0(bss.name,".",iter)
        if(temp.name %in% names(pooledDataObject$bssResults)){
          iter = iter+1
        }else{
          bss.name=temp.name
          continue=FALSE
        }
      }
    }
    pooledDataObject$bssResults[[bss.name]]=BestSubsetList
  }else{
    pooledDataObject$bssResults = list()
    pooledDataObject$bssResults[[bss.name]]=BestSubsetList
  }
  return(pooledDataObject)
}



#' Forward Abridged Best Subset Selection
#' @description
#' Conducts a forward abridged best subset selection on a set of genes and returns a ranked list of
#' all gene combinations as well as a shortlist of the best gene combination for each model size
#' (i.e. number of genes). Ranking can be done by AUC, AUPRC, or Average AUC. \cr
#' 
#' Forward abridged best subset selection is a modified version of best subset selection where
#' BSS is first run to generate the best 1-gene through N-gene models, and then forward search
#' is subsequently used to determine the best (N+1)-gene and larger models.
#' 
#' \strong{WARNING:} This function is not as computationally intensive as \code{runBestSubsetSelection},
#' but it can still take a long time, so be careful with using it on very large gene sets.
#' @param pooledDataObject A MANATEE pooledDataObject
#' @param filterObject A MANATEE filterObject
#' @param perf.meas Which diagnostic performance metric to use. Options are "AUC", "AUPRC", or "AvgAUC" (Default: "AUC")
#' @param upgenes If you want to manually input genes instead of providing a filterObject, this is a vector of the upgenes in your signature
#' @param downgenes If you want to manually input genes instead of providing a filterObject, this is a vector of the downgenes in your signature
#' @param bss.min Initial/minimum model size for the best subset selection. For example, if \code{bss.min} is 1 then the best subset selection will start by testing all 1-gene models.
#' @param bss.max Final/maximum model size for the best subset selection. For example, if \code{bss.max} is 5 then the best subset selection will stop after testing all 5-gene models,
#' and forward search will be used for the rest of the model sizes.
#' @param forward.max Maximum model size for the forward search. For example, if \code{forward.max} is 10, then the forward search will stop after selecting a 10-gene model.
#' If left blank, then it will be set as the combined length of your upgenes and downgenes.
#' @param weights \emph{Only for perf.meas=AvgAUC:} If left blank, then the simple average of class-specific AUCs will be used to
#' calculate the AvgAUC. However, if you want to take a weighted average instead, then pass in a vector of weights. Note that the
#' order of the weights will be matched to the order of the names in the \code{otherNames} vector.
#' @param caseNames \emph{Only for perf.meas=AvgAUC:} name of the class(es) that you're considering to be your case.
#' If left blank, then \code{caseNames} will be the unique set of names in \code{$pheno$group} that correspond to a class label of 1
#' @param otherNames \emph{Only for perf.meas=AvgAUC:} name of the class(es) that you're considering to be your controls.
#' If left blank, then \code{otherNames} will be the unique set of names in \code{$pheno$group} that correspond to a class label of 0.
#' The order of otherNames will be set by running the \code{sort} function, so if you are using weights, make sure that your \code{weights}
#' vector and \code{otherNames} vector match up correctly
#' @param group.colname \emph{Only for perf.meas=AvgAUC or perf.meas=All:} this designates the column name of the group vector within the pooledDataObject$pheno
#' dataframe (Default: "group")
#' @param init.pos If you want to start off with some upgenes already locked in, then indicate those genes with \code{init.pos}
#' @param init.neg If you want to start off with some downgenes genes already locked in, then indicate those genes with \code{init.pos}
#' @param numCores Number of cores to use for parallelization. If left blank, smart core usage will be used
#' @details
#' Because best subset selection is computationally intensive, it is often unfeasible to run
#' it without first reducing the number of initial genes by a substantial amount. This function
#' is an alternative approach that allows a BSS-like method to be applied to larger gene sets.
#' Unlike a true BSS, it will not always select globally optimal gene combinations, but it will
#' likely select more optimal gene combinations than a greedy search like forward or backward search.
#' 
#' Other alternative BSS-like methods include \code{backwardAbridgedBSS}, \code{forwardGroupedBSS},
#' and \code{backwardGroupedBSS}
#' @note As of now, the accepted diagnostic performance metrics are AUC (area under the Reciever Operating Curve),
#' AUPRC (area under the Precision Recall Curve), and AvgAUC
#' (the average AUC across all comparisons between your cases and each of your "control" categories). \cr
#' 
#' Use \code{bss.min}, \code{bss.max}, and \code{forward.max} to control how the abridged BSS is run. For example,
#' if \code{bss.min}=2, \code{bss.max}=10, and \code{forward.max}=30, then BSS will be used to generate all 2-gene
#' through 10-gene models, and then forward search will be used to generate all 11-gene through 30-gene models. \cr
#' 
#' Note that if any inital genes are indicated with either \code{init.pos} or \code{init.neg}, then \code{bss.min}
#' will be disregarded.
#' @return
#'   \item{allResults}{Data table containing the combined results of the entire best subset selection and forward search}
#'   \item{bestResults}{Data table containing the best results for each tested model size}
#'   \item{abridgedBSSdescription}{List with information that describes the input parameters for the abridged best subset selection}
#' @examples
#' 
#' @export
#' @import parallel zoo pracma MetaIntegrator ROCR data.table pbapply
#' @importFrom pbmcapply pbmclapply
#' @author Aditya Rao
forwardAbridgedBSS <- function(pooledDataObject, filterObject=NULL, perf.meas="AUC", upgenes, downgenes, bss.min=1, bss.max, forward.max=NULL,
                               weights=NULL, caseNames=NULL, otherNames=NULL, group.colname="group", init.pos=NULL, init.neg=NULL, numCores=NULL){
  progress.bar.type="perGeneNum"
  if(!checkManateeObject(pooledDataObject,"pooledDataObject")){stop("Invalid pooledDataObject provided")}
  if(!is.null(filterObject)){
    upgenes = filterObject$upGeneNames
    downgenes = filterObject$downGeneNames
  }
  perf.meas = tolower(perf.meas)
  perf.meas = switch(perf.meas, "auc"="AUC", "auprc"="AUPRC", "avgauc"="AvgAUC")
  if(!(perf.meas %in% c("AUC","AUPRC","AvgAUC"))){stop("perf.meas must be \"AUC\", \"AUPRC\", or \"AvgAUC\"")}
  if(any(!init.pos %in% upgenes)){stop("There are some genes in init.pos that do not appear in upgenes")}
  if(any(!init.neg %in% downgenes)){stop("There are some genes in init.neg that do not appear in downgenes")}
  if(is.null(forward.max)){
    forward.max = length(c(upgenes,downgenes))
  }
  if(forward.max > length(c(upgenes,downgenes))){
    warning("forward.max is greater than the total number of genes and will be ignored")
    forward.max = length(c(upgenes,downgenes))
  }
  if(bss.max > forward.max){stop("bss.max cannot be greater than forward.max")}
  if(bss.max == forward.max){stop("bss.max is equal to forward.max - this is equivalent to just running a best subset selection")}
  if(length(c(init.pos,init.neg)) > 0){
    if(length(c(init.pos,init.neg)) >= bss.max){stop("bss.max must be greater than the number of initial genes, otherwise just a forward search would be run")}
    bss.max = bss.max-length(c(init.pos,init.neg))
    bss.min = 1
  }
  
  BestSubsetList = .bestSubsetSelectionHelper(perf.meas=perf.meas, geneMtx=pooledDataObject$genes, class=pooledDataObject$class, pos.genes=upgenes, neg.genes=downgenes,
                                              bss.min=bss.min, bss.max=bss.max, labelVec=pooledDataObject$pheno[,group.colname], weights=weights, caseNames=caseNames, otherNames=otherNames,
                                              init.pos=init.pos, init.neg=init.neg, numCores=numCores, progress.bar.type=progress.bar.type)
  
  max.subset = BestSubsetList[[paste0(bss.max,".Genes")]][1,]
  init.pos2 = unlist(strsplit(max.subset$upgenes,split=" / "))
  init.neg2 = unlist(strsplit(max.subset$downgenes,split=" / "))
  if(length(init.pos2)==0){
    init.pos2=NULL
  }
  if(length(init.neg2)==0){
    init.neg2=NULL
  }
  
  ForwardSearchList = .forwardSearchHelper(pos.genes=upgenes, neg.genes=downgenes, perf.meas=perf.meas, forwardThresh=-1, geneMtx=pooledDataObject$genes, class=pooledDataObject$class, labelVec=pooledDataObject$pheno[,group.colname],
                                           weights=weights, caseNames=caseNames, otherNames=otherNames, init.pos=init.pos2, init.neg=init.neg2, max.genes=forward.max)
  
  forwardDT = ForwardSearchList$bestgeneDT[order(-rank(numGenes))]
  abridgedDT = rbind(forwardDT,BestSubsetList$bssDT)
  abridgedDT.best = abridgedDT[,.SD[1,],by=numGenes]
  
  if(perf.meas=="AvgAUC"){
    my.caseNames = BestSubsetList$bssDescription$caseNames
    my.otherNames = BestSubsetList$bssDescription$otherNames
  }else{
    my.caseNames = NULL
    my.otherNames = NULL
  }
  description = list(type="Forward Abridged Best Subset Selection", all.upgenes=upgenes, all.downgenes=downgenes,
                     perf.meas=perf.meas, bss.min=bss.min, bss.max=bss.max, forward.max=forward.max, weights=weights,
                     caseNames=my.caseNames, otherNames=my.otherNames, init.pos=init.pos, init.neg=init.neg)
  
  return(list(allResults=abridgedDT, bestResults=abridgedDT.best, abridgedBSSdescription=description))
}



#' Backward Abridged Best Subset Selection
#' @description
#' Conducts a backward abridged best subset selection on a set of genes and returns a ranked list of
#' all gene combinations as well as a shortlist of the best gene combination for each model size
#' (i.e. number of genes). Ranking can be done by AUC, AUPRC, or Average AUC. \cr
#' 
#' Backward abridged best subset selection is a modified version of best subset selection where
#' a backward search is first run until the best (N+1)-gene model is determined, and then BSS is run
#' to generate the best N-gene through 1-gene models.
#' 
#' \strong{WARNING:} This function is not as computationally intensive as \code{runBestSubsetSelection},
#' but it can still take a long time, so be careful with using it on very large gene sets.
#' @param pooledDataObject A MANATEE pooledDataObject
#' @param filterObject A MANATEE filterObject
#' @param perf.meas Which diagnostic performance metric to use. Options are "AUC", "AUPRC", or "AvgAUC" (Default: "AUC")
#' @param upgenes If you want to manually input genes instead of providing a filterObject, this is a vector of the upgenes in your signature
#' @param downgenes If you want to manually input genes instead of providing a filterObject, this is a vector of the downgenes in your signature
#' @param bss.min Final/minimum model size for the best subset selection. For example, if \code{bss.min} is 1 then the best subset selection will stop after testing all 1-gene models.
#' @param bss.max Initial/maximum model size for the best subset selection. For example, if \code{bss.max} is 5 then the backward search will stop after finding the best 6-gene model,
#' and then best subset selection will start by testing all 5-gene models.
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
#' @param numCores Number of cores to use for parallelization. If left blank, smart core usage will be used
#' @details
#' Because best subset selection is computationally intensive, it is often unfeasible to run
#' it without first reducing the number of initial genes by a substantial amount. This function
#' is an alternative approach that allows a BSS-like method to be applied to larger gene sets.
#' Unlike a true BSS, it will not always select globally optimal gene combinations, but it will
#' likely select more optimal gene combinations than a greedy search like forward or backward search.
#' 
#' Other alternative BSS-like methods include \code{forwardAbridgedBSS}, \code{forwardGroupedBSS},
#' and \code{backwardGroupedBSS}
#' @note As of now, the accepted diagnostic performance metrics are AUC (area under the Reciever Operating Curve),
#' AUPRC (area under the Precision Recall Curve), and AvgAUC
#' (the average AUC across all comparisons between your cases and each of your "control" categories). \cr
#' 
#' Use \code{bss.min} and \code{bss.max} to control how the abridged BSS is run. For example, if \code{bss.min}=2,
#' \code{bss.max}=10, and there are 30 total genes, then backward search will be used to generate all 30-gene
#' through 11-gene models, and then BSS will be used to generate all 10-gene through 2-gene models.
#' @return
#'   \item{allResults}{Data table containing the combined results of the entire backward search and best subset selection}
#'   \item{bestResults}{Data table containing the best results for each tested model size}
#'   \item{abridgedBSSdescription}{List with information that describes the input parameters for the abridged best subset selection}
#' @examples
#' 
#' @export
#' @import parallel zoo pracma MetaIntegrator ROCR data.table pbapply
#' @importFrom pbmcapply pbmclapply
#' @author Aditya Rao
backwardAbridgedBSS <- function(pooledDataObject, filterObject=NULL, perf.meas="AUC", upgenes, downgenes, bss.min=1, bss.max,
                                weights=NULL, caseNames=NULL, otherNames=NULL, group.colname="group", numCores=NULL){
  progress.bar.type="perGeneNum"
  if(!checkManateeObject(pooledDataObject,"pooledDataObject")){stop("Invalid pooledDataObject provided")}
  if(!is.null(filterObject)){
    upgenes = filterObject$upGeneNames
    downgenes = filterObject$downGeneNames
  }
  perf.meas = tolower(perf.meas)
  perf.meas = switch(perf.meas, "auc"="AUC", "auprc"="AUPRC", "avgauc"="AvgAUC")
  if(!(perf.meas %in% c("AUC","AUPRC","AvgAUC"))){stop("perf.meas must be \"AUC\", \"AUPRC\", or \"AvgAUC\"")}
  if(bss.max > length(c(upgenes,downgenes))){
    stop("bss.max is greater than the combined length of your upgenes and downgenes")
  }
  if(bss.max == length(c(upgenes,downgenes))){
    stop("bss.max is equal to the combined length of your upgenes and downgenes - this is equivalent to just running a best subset selection")
  }
  
  BackwardSearchList = .backwardSearchHelper(pos.genes=upgenes, neg.genes=downgenes, perf.meas=perf.meas, backwardThresh=-1, geneMtx=pooledDataObject$genes, class=pooledDataObject$class, 
                                             labelVec=pooledDataObject$pheno[,group.colname], weights=weights, caseNames=caseNames, otherNames=otherNames, min.genes=bss.max)
  
  init.pos = unlist(strsplit(BackwardSearchList$bestgeneDT[numGenes==bss.max,upgenes],split=" / "))
  init.neg = unlist(strsplit(BackwardSearchList$bestgeneDT[numGenes==bss.max,downgenes],split=" / "))
  if(length(init.pos)==0){
    init.pos=NULL
  }
  if(length(init.neg)==0){
    init.neg=NULL
  }
  
  BestSubsetList = .bestSubsetSelectionHelper(perf.meas=perf.meas, geneMtx=pooledDataObject$genes, class=pooledDataObject$class, pos.genes=init.pos, neg.genes=init.neg, bss.min=bss.min, bss.max=bss.max,
                                              labelVec=pooledDataObject$pheno[,group.colname], weights=weights, caseNames=caseNames, otherNames=otherNames, numCores=numCores, progress.bar.type=progress.bar.type)
  
  backwardDT = BackwardSearchList$bestgeneDT
  abridgedDT = rbind(backwardDT,BestSubsetList$bssDT)
  abridgedDT.best = abridgedDT[,.SD[1,],by=numGenes]
  
  if(perf.meas=="AvgAUC"){
    my.caseNames = BestSubsetList$bssDescription$caseNames
    my.otherNames = BestSubsetList$bssDescription$otherNames
  }else{
    my.caseNames = NULL
    my.otherNames = NULL
  }
  description = list(type="Backward Abridged Best Subset Selection", all.upgenes=upgenes, all.downgenes=downgenes, perf.meas=perf.meas,
                     bss.min=bss.min, bss.max=bss.max, weights=weights, caseNames=my.caseNames, otherNames=my.otherNames)
  
  return(list(allResults=abridgedDT, bestResults=abridgedDT.best, abridgedBSSdescription=description))
}



#' Forward Grouped Best Subset Selection
#' @description
#' Conducts a forward grouped best subset selection on a set of genes and returns a ranked list of
#' all gene combinations as well as a shortlist of the best gene combination for each model size
#' (i.e. number of genes). Ranking can be done by AUC, AUPRC, or Average AUC. \cr
#' 
#' Forward grouped best subset selection is a modified version of best subset selection where
#' BSS is iteratively run to get N best subsets at a time. More specifically, a BSS will first
#' be run to generate the best 1-gene through N-gene combinations, after which the genes from the
#' best N-gene signature are locked in. Next, a second BSS is run to generate the best (N+1)-gene
#' through 2N-gene signatures, and the genes from the best 2N-gene signature are locked in. This
#' process continues until the best signatures for all desired model sizes are generated.
#' 
#' \strong{WARNING:} This function is not as computationally intensive as \code{runBestSubsetSelection},
#' but it can still take a long time, so be careful with using it on very large gene sets.
#' @param pooledDataObject A MANATEE pooledDataObject
#' @param filterObject A MANATEE filterObject
#' @param perf.meas Which diagnostic performance metric to use. Options are "AUC", "AUPRC", or "AvgAUC" (Default: "AUC")
#' @param upgenes If you want to manually input genes instead of providing a filterObject, this is a vector of the upgenes in your signature
#' @param downgenes If you want to manually input genes instead of providing a filterObject, this is a vector of the downgenes in your signature
#' @param group.size The number of different model sizes to generate with each BSS
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
#' @param numCores Number of cores to use for parallelization. If left blank, smart core usage will be used
#' @details
#' Because best subset selection is computationally intensive, it is often unfeasible to run
#' it without first reducing the number of initial genes by a substantial amount. This function
#' is an alternative approach that allows a BSS-like method to be applied to larger gene sets.
#' Unlike a true BSS, it will not always select globally optimal gene combinations, but it will
#' likely select more optimal gene combinations than a greedy search like forward or backward search.
#' 
#' Other alternative BSS-like methods include \code{forwardAbridgedBSS}, \code{backwardAbridgedBSS},
#' and \code{backwardGroupedBSS}
#' @note As of now, the accepted diagnostic performance metrics are AUC (area under the Reciever Operating Curve),
#' AUPRC (area under the Precision Recall Curve), and AvgAUC
#' (the average AUC across all comparisons between your cases and each of your "control" categories). \cr
#' 
#' Use \code{group.size} to control how the grouped BSS is run. For example, if \code{group.size}=5 and there are
#' 13 genes total, then the first BSS will generate all 1-gene through 5-gene models, the second BSS will generate
#' all 6-gene through 10-gene models, and the last BSS will generate all 11-gene through 13-gene models.
#' @return
#'   \item{allResults}{Data table containing the combined results of the entire grouped best subset selection}
#'   \item{bestResults}{Data table containing the best results for each tested model size}
#'   \item{groupedBSSdescription}{List with information that describes the input parameters for the grouped best subset selection}
#' @examples
#' 
#' @export
#' @import parallel zoo pracma MetaIntegrator ROCR data.table pbapply
#' @importFrom pbmcapply pbmclapply
#' @author Aditya Rao
forwardGroupedBSS <- function(pooledDataObject, filterObject=NULL, perf.meas="AUC", upgenes, downgenes, group.size, 
                              weights=NULL, caseNames=NULL, otherNames=NULL, group.colname="group", max.genes=NULL, init.pos=NULL, init.neg=NULL, numCores=NULL){
  progress.bar.type="perGeneNum"
  if(!checkManateeObject(pooledDataObject,"pooledDataObject")){stop("Invalid pooledDataObject provided")}
  if(!is.null(filterObject)){
    upgenes = filterObject$upGeneNames
    downgenes = filterObject$downGeneNames
  }
  perf.meas = tolower(perf.meas)
  perf.meas = switch(perf.meas, "auc"="AUC", "auprc"="AUPRC", "avgauc"="AvgAUC")
  if(!(perf.meas %in% c("AUC","AUPRC","AvgAUC"))){stop("perf.meas must be \"AUC\", \"AUPRC\", or \"AvgAUC\"")}
  if(any(!init.pos %in% upgenes)){stop("There are some genes in init.pos that do not appear in upgenes")}
  if(any(!init.neg %in% downgenes)){stop("There are some genes in init.neg that do not appear in downgenes")}
  if(is.null(init.pos)){
    pos=NULL
  }else{
    pos=init.pos
  }
  if(is.null(init.neg)){
    neg=NULL
  }else{
    neg=init.neg
  }
  if(is.null(max.genes)){
    max.genes = length(upgenes) + length(downgenes)
  }else if(max.genes > (length(upgenes) + length(downgenes))){
    warning("The provided max.genes is greater than the total number of genes and will be ignored")
    max.genes = length(upgenes) + length(downgenes)
  }else if(max.genes < length(c(pos,neg))){
    stop("max.genes cannot be less than the number of initial genes")
  }
  if(group.size >= max.genes){stop("group.size cannot be greater than or equal to max.genes")}
  
  #calculate number of total gene combinations DO THIS LATER ASFKGAFJGAGJNAGJNAGJAGJNAEJNAEJNDBNADFKBNADKBNADFKJBNADFJBNADKBNAKDBN
  combo.count = 0
  all.length = length(c(upgenes,downgenes))
  total.iter = ceil(all.length/group.size)
  for(i in 1:total.iter){
    combos = c(1:min(all.length,group.size))
    for(j in combos){
      combo.count=combo.count+choose(all.length,j)
    }
    all.length = all.length-group.size
  }
  cat(sprintf("Across the entire Grouped BSS, %s gene combinations will be tested\n\n",combo.count))
  
  continue=TRUE
  groupedDT = data.table()
  curr.bss.min = length(c(pos,neg))
  curr.bss.max = curr.bss.min+group.size
  if(curr.bss.min == 0){curr.bss.min=1}
  while(continue){
    BestSubsetList = .bestSubsetSelectionHelper(perf.meas=perf.meas, geneMtx=pooledDataObject$genes, class=pooledDataObject$class, pos.genes=upgenes,
                                                neg.genes=downgenes, bss.min=1, bss.max=curr.bss.max-curr.bss.min+1, labelVec=pooledDataObject$pheno[,group.colname],
                                                weights=weights, caseNames=caseNames, otherNames=otherNames, init.pos=pos, init.neg=neg, numCores=numCores, progress.bar.type=progress.bar.type)
    bssDT = BestSubsetList$bssDT
    groupedDT = rbind(bssDT,groupedDT)
    pos = unlist(strsplit(bssDT[numGenes==curr.bss.max,upgenes[1]],split=" / "))
    neg = unlist(strsplit(bssDT[numGenes==curr.bss.max,downgenes[1]],split=" / "))
    curr.bss.min = length(c(pos,neg))+1
    curr.bss.max = curr.bss.min+group.size-1
    if(curr.bss.max >= max.genes){
      curr.bss.max = max.genes
      continue = FALSE
    }
  }
  BestSubsetList = .bestSubsetSelectionHelper(perf.meas=perf.meas, geneMtx=pooledDataObject$genes, class=pooledDataObject$class, pos.genes=upgenes,
                                              neg.genes=downgenes, bss.min=1, bss.max=curr.bss.max-curr.bss.min+1, labelVec=pooledDataObject$pheno[,group.colname],
                                              weights=weights, caseNames=caseNames, otherNames=otherNames, init.pos=pos, init.neg=neg, numCores=numCores, progress.bar.type=progress.bar.type)
  groupedDT = rbind(BestSubsetList$bssDT,groupedDT)
  groupedDT.best = groupedDT[,.SD[1,],by=numGenes]
  
  if(perf.meas=="AvgAUC"){
    my.caseNames = BestSubsetList$bssDescription$caseNames
    my.otherNames = BestSubsetList$bssDescription$otherNames
  }else{
    my.caseNames = NULL
    my.otherNames = NULL
  }
  description = list(type="Forward Grouped Best Subset Selection", all.upgenes=upgenes, all.downgenes=downgenes, perf.meas=perf.meas,
                     group.size=group.size, max.genes=max.genes, init.pos=init.pos, init.neg=init.neg, weights=weights, caseNames=my.caseNames, otherNames=my.otherNames)
  
  return(list(allResults=groupedDT, bestResults=groupedDT.best, groupedBSSdescription=description))
}



#' Backward Grouped Best Subset Selection
#' @description
#' Conducts a backward grouped best subset selection on a set of genes and returns a ranked list of
#' all gene combinations as well as a shortlist of the best gene combination for each model size
#' (i.e. number of genes). Ranking can be done by AUC, AUPRC, or Average AUC. \cr
#' 
#' Backward grouped best subset selection is a modified version of best subset selection where
#' BSS is iteratively run to get N best subsets at a time. If there are K total genes, a BSS will first
#' be run to generate the best K-gene through (K-N)-gene combinations, after which the genes from the
#' best (K-N)-gene signature are locked in. Next, a second BSS is run to generate the best (K-N-1)-gene
#' through (K-2N)-gene signatures, and the genes from the best (K-2N)-gene signature are locked in. This
#' process continues until the best signatures for all desired model sizes are generated.
#' 
#' \strong{WARNING:} This function is not as computationally intensive as \code{runBestSubsetSelection},
#' but it can still take a long time, so be careful with using it on very large gene sets.
#' @param pooledDataObject A MANATEE pooledDataObject
#' @param filterObject A MANATEE filterObject
#' @param perf.meas Which diagnostic performance metric to use. Options are "AUC", "AUPRC", or "AvgAUC" (Default: "AUC")
#' @param upgenes If you want to manually input genes instead of providing a filterObject, this is a vector of the upgenes in your signature
#' @param downgenes If you want to manually input genes instead of providing a filterObject, this is a vector of the downgenes in your signature
#' @param group.size The number of different model sizes to generate with each BSS
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
#' Because best subset selection is computationally intensive, it is often unfeasible to run
#' it without first reducing the number of initial genes by a substantial amount. This function
#' is an alternative approach that allows a BSS-like method to be applied to larger gene sets.
#' Unlike a true BSS, it will not always select globally optimal gene combinations, but it will
#' likely select more optimal gene combinations than a greedy search like forward or backward search.
#' 
#' Other alternative BSS-like methods include \code{forwardAbridgedBSS}, \code{backwardAbridgedBSS},
#' and \code{forwardGroupedBSS}
#' @note As of now, the accepted diagnostic performance metrics are AUC (area under the Reciever Operating Curve),
#' AUPRC (area under the Precision Recall Curve), and AvgAUC
#' (the average AUC across all comparisons between your cases and each of your "control" categories). \cr
#' 
#' Use \code{group.size} to control how the grouped BSS is run. For example, if \code{group.size}=5 and there are
#' 14 genes total, then the first BSS will generate all 13-gene through 9-gene models, the second BSS will generate
#' all 8-gene through 4-gene models, and the last BSS will generate all 3-gene through 1-gene models.
#' @return
#'   \item{allResults}{Data table containing the combined results of the entire grouped best subset selection}
#'   \item{bestResults}{Data table containing the best results for each tested model size}
#'   \item{groupedBSSdescription}{List with information that describes the input parameters for the grouped best subset selection}
#' @examples
#' 
#' @export
#' @import parallel zoo pracma MetaIntegrator ROCR data.table pbapply
#' @importFrom pbmcapply pbmclapply
#' @author Aditya Rao
backwardGroupedBSS <- function(pooledDataObject, filterObject=NULL, perf.meas="AUC", upgenes, downgenes, group.size, 
                               weights=NULL, caseNames=NULL, otherNames=NULL, group.colname="group", min.genes=NULL, numCores=NULL){
  progress.bar.type="perGeneNum"
  if(!checkManateeObject(pooledDataObject,"pooledDataObject")){stop("Invalid pooledDataObject provided")}
  if(!is.null(filterObject)){
    upgenes = filterObject$upGeneNames
    downgenes = filterObject$downGeneNames
  }
  perf.meas = tolower(perf.meas)
  perf.meas = switch(perf.meas, "auc"="AUC", "auprc"="AUPRC", "avgauc"="AvgAUC")
  if(!(perf.meas %in% c("AUC","AUPRC","AvgAUC"))){stop("perf.meas must be \"AUC\", \"AUPRC\", or \"AvgAUC\"")}
  if(is.null(min.genes)){
    min.genes = 1
  }else if(min.genes >= length(c(upgenes,downgenes))){
    warning("The provided min.genes is greater than or equal to the total number of genes and will be ignored")
    min.genes = 1
  }
  if(group.size>=(length(c(upgenes,downgenes))-min.genes)){stop("group.size is too large")}
  
  #calculate number of total gene combinations
  combo.count = 0
  all.length = length(c(upgenes,downgenes))
  total.times = all.length-min.genes-1
  total.iter = ceil(total.times/group.size)
  for(i in 1:total.iter){
    combos = c((all.length-1):max(min.genes,all.length-group.size))
    for(j in combos){
      combo.count=combo.count+choose(all.length,j)
    }
    all.length = all.length-group.size
  }
  combo.count = combo.count+1
  if(all.length-1 %% group.size == 1){combo.count = combo.count+2}
  cat(sprintf("Across the entire Grouped BSS, %s gene combinations will be tested\n\n",combo.count))
  
  continue=TRUE
  pos=upgenes
  neg=downgenes
  groupedDT = data.table()
  curr.bss.max = length(c(pos,neg))
  curr.bss.min = curr.bss.max-group.size
  while(continue){
    BestSubsetList = .bestSubsetSelectionHelper(perf.meas=perf.meas, geneMtx=pooledDataObject$genes, class=pooledDataObject$class, pos.genes=pos, neg.genes=neg, bss.min=curr.bss.min, bss.max=curr.bss.max,
                                                labelVec=pooledDataObject$pheno[,group.colname], weights=weights, caseNames=caseNames, otherNames=otherNames, numCores=numCores, progress.bar.type=progress.bar.type)
    bssDT = BestSubsetList$bssDT
    groupedDT = rbind(groupedDT,bssDT)
    pos = unlist(strsplit(bssDT[numGenes==curr.bss.min,upgenes[1]],split=" / "))
    neg = unlist(strsplit(bssDT[numGenes==curr.bss.min,downgenes[1]],split=" / "))
    curr.bss.max = length(c(pos,neg))-1
    curr.bss.min = curr.bss.max-group.size+1
    if(curr.bss.min <= min.genes){
      curr.bss.min = min.genes
      continue = FALSE
    }
  }
  BestSubsetList = .bestSubsetSelectionHelper(perf.meas=perf.meas, geneMtx=pooledDataObject$genes, class=pooledDataObject$class, pos.genes=pos, neg.genes=neg, bss.min=curr.bss.min, bss.max=curr.bss.max,
                                              labelVec=pooledDataObject$pheno[,group.colname], weights=weights, caseNames=caseNames, otherNames=otherNames, numCores=numCores, progress.bar.type=progress.bar.type)
  groupedDT = rbind(groupedDT,BestSubsetList$bssDT)
  groupedDT.best = groupedDT[,.SD[1,],by=numGenes]
  
  if(perf.meas=="AvgAUC"){
    my.caseNames = BestSubsetList$bssDescription$caseNames
    my.otherNames = BestSubsetList$bssDescription$otherNames
  }else{
    my.caseNames = NULL
    my.otherNames = NULL
  }
  description = list(type="Backward Grouped Best Subset Selection", all.upgenes=upgenes, all.downgenes=downgenes, perf.meas=perf.meas,
                     group.size=group.size, min.genes=min.genes, weights=weights, caseNames=my.caseNames, otherNames=my.otherNames)
  
  return(list(allResults=groupedDT, bestResults=groupedDT.best, groupedBSSdescription=description))
}



###-###-###-###-###-###-##-#
###   plotBSSResults()   ###
###-###-###-###-###-###-##-#

#DESCRIPTION
#This function plots the results of a best subset selection by displaying the
#performance of the best subset at each gene number

#PARAMETERS
#bssResults - the output of running bestSubsetSelection on a MANATEE pooledDataObject
#perf.meas - what measure of diagnostic performance should be plotted? Options are: "AUC", "AUPRC", "AvgAUC", and "All"
#gene.nums - the sequence of gene numbers to be plotted
#breaks - the graphical breaks on the y-axis
#blurb - if TRUE, the function will return a "blurb" about the best subset at each gene number

#RETURN VALUE
#if blurb = TRUE, a vector containing a short descriptive sentence about the best subset at each gene number will be output

#REQUIRED PACKAGES: ggplot2

#code NOTE THIS IS NOT GENERALIZABLE YET, NEED TO DO MORE TESTING
plotBSSResults <- function(bssResults, perf.meas = "AUC", gene.nums = NULL, breaks = NULL, blurb = TRUE){
  if(!checkManateeObject(bssResults,"bssResults")){stop("Invalid bssResults provided")}
  if(!(perf.meas %in% c("AUC","AUPRC","AvgAUC"))){stop("perf.meas must be \"AUC\", \"AUPRC\", or \"AvgAUC\"")}
  
  if(is.null(gene.nums)){
    gene.nums = c(bssResults[[1]]$numGenes[1]:bssResults[[length(bssResults)-2]]$numGenes[1])
  }
  n=length(bssResults)-2
  best.meas = rep(0,n)
  best.measlo = rep(0,n)
  best.measup = rep(0,n)
  best.upgenes = rep("",n)
  best.downgenes = rep("",n)
  best.blurb = rep("",n)
  for(i in 1:(n)){
    best.meas[i]=bssResults[[i]][,perf.meas][1] #might need to change this
    best.measlo[i]=bssResults[[i]][,paste0(perf.meas,".ci.lower")][1] #might need to change this
    best.measup[i]=bssResults[[i]][,paste0(perf.meas,".ci.upper")][1] #might need to change this
    
    best.upgenes[i]=bssResults[[i]]$upgenes[1]
    best.downgenes[i]=bssResults[[i]]$downgenes[1]
    best.blurb[i]=sprintf("%s Genes, %s=%s (95%% CI %s-%s): upgenes are %s; downgenes are %s",bssResults[[i]]$numGenes[1],perf.meas,round(best.meas[i],3),
                          round(best.measlo[i],3),round(best.measup[i],3),paste(unlist(strsplit(best.upgenes[i],split =" / ")),collapse=", "),
                          paste(unlist(strsplit(best.downgenes[i],split =" / ")),collapse=", "))
  }
  #best AUCs line graph
  best_plot = data.frame(cbind(best.meas,best.measlo,best.measup,gene.nums))
  if(is.null(breaks)){
    plot = ggplot(best_plot,aes(y=best.meas,x=gene.nums))+geom_line()+geom_point()+scale_x_reverse(breaks=c(gene.nums))+
      geom_errorbar(aes(ymax=best.measup,ymin=best.measlo),width=0.5)+
      labs(title=paste("Best Subset for Each Gene Number by",perf.meas),y=perf.meas,x="Number of Genes")
  } else{
    plot = ggplot(best_plot,aes(y=best.meas,x=gene.nums))+geom_line()+geom_point()+scale_x_reverse(breaks=c(gene.nums))+
      scale_y_continuous(breaks=breaks,limits=c(min(breaks),max(breaks)))+geom_errorbar(aes(ymax=best.measup,ymin=best.measlo),width=0.5)+
      labs(title=paste("Best Subset for Each Gene Number by",perf.meas),y=perf.meas,x="Number of Genes")
  }
  
  if(blurb){cat(paste(best.blurb,collapse = "\n"))}
  plot
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#--------------------------------------------------------Helper Functions--------------------------------------------------------#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' Helper Function for Best Subset Selection
#' @description
#' Conducts best subset selection on a set of genes and returns a ranked list of
#' all gene combinations, separated by model size (i.e. number of genes). Ranking can be
#' done by AUC, AUPRC, or Average AUC. \cr
#' 
#' \strong{WARNING:} the time this function takes increases exponentially, so be very
#' careful with it. If running more than ~15 genes, it will likely require a multi-core server.
#' @param perf.meas Which diagnostic performance metric to use. Options are "AUC", "AUPRC", "AvgAUC", or "All" (Default: "AUC")
#' @param geneMtx The gene matrix to use for calculating signature scores (rows are genes and columns are samples)
#' @param class Vector of class labels
#' @param pos.genes Vector of upgenes
#' @param neg.genes Vector of downgenes
#' @param bss.min Minimum model size. For example, if \code{bss.min} is 1 and \code{bss.max} is 5, then all of the combinations of 1, 2, 3, 4, and 5 genes will be tested (Default: 1)
#' @param bss.max Maximum model size. For example, if \code{bss.min} is 1 and \code{bss.max} is 5, then all of the combinations of 1, 2, 3, 4, and 5 genes will be tested.
#' If left blank, then it will be set as the combined length of your upgenes and downgenes.
#' @param labelVec \emph{Only for perf.meas=AvgAUC or perf.meas=All:} Vector of sample category labels
#' @param weights \emph{Only for perf.meas=AvgAUC or perf.meas=All:} If left blank, then the simple average of class-specific AUCs will be used to
#' calculate the AvgAUC. However, if you want to take a weighted average instead, then pass in a vector of weights. Note that the
#' order of the weights will be matched to the order of the names in the \code{otherNames} vector.
#' @param caseNames \emph{Only for perf.meas=AvgAUC or perf.meas=All:} name of the class(es) that you're considering to be your case.
#' If left blank, then \code{caseNames} will be the unique set of names in \code{$pheno$group} that correspond to a class label of 1
#' @param otherNames \emph{Only for perf.meas=AvgAUC or perf.meas=All:} name of the class(es) that you're considering to be your controls.
#' If left blank, then \code{otherNames} will be the unique set of names in \code{$pheno$group} that correspond to a class label of 0.
#' The order of otherNames will be set by running the \code{sort} function, so if you are using weights, make sure that your \code{weights}
#' vector and \code{otherNames} vector match up correctly
#' @param init.pos If you want to start off with some upgenes already locked in, then indicate those genes with \code{init.pos}
#' @param init.neg If you want to start off with some downgenes genes already locked in, then indicate those genes with \code{init.pos}
#' @param sort.by \emph{Only for perf.meas=All:} Which performance measure should be used to rank the genes (Default: "AUC")
#' @param numCores Number of cores to use for parallelization. If left blank, smart core usage will be used
#' @param progress.bar.type If \code{progress.bar.type}="perGeneNum", then a progress bar is shown for each model size (i.e. each number of genes).
#' If \code{progress.bar.type}="overall", then a single progress bar is used for the entire best subset selection (note that this doesn't work with parallelization).
#' @return
#'   \item{\emph{X}.Genes}{The ranked results for each model size is stored in a separate data frame that contains information about
#'   performance, the included genes, and a short blurb for each gene combination tested}
#'   \item{bssDT}{Data table containing the combined results of the entire best subset selection}
#'   \item{bssDescription}{List with information that describes the input parameters for the best subset selection}
#' @import parallel zoo pracma MetaIntegrator ROCR data.table pbapply
#' @importFrom pbmcapply pbmclapply
#' @author Aditya Rao
.bestSubsetSelectionHelper <- function(perf.meas="AUC", geneMtx, class, pos.genes, neg.genes, bss.min=1, bss.max=NULL, labelVec, weights=NULL,
                                       caseNames=NULL, otherNames=NULL, init.pos=NULL, init.neg=NULL, sort.by="AUC", numCores=NULL, progress.bar.type="perGeneNum"){
  if(!is.null(caseNames) && any(!c(caseNames %in% unique(labelVec[class==1])))){
    stop("If caseNames is set by the user, then all of the provided names must be labels that are used for case samples (i.e. samples that are labeled with 1 in the provided class vector)")
  }
  if(!is.null(otherNames) && any(!c(otherNames %in% unique(labelVec[class==0])))){
    stop("If otherNames is set by the user, then all of the provided names must be labels that are used for control samples (i.e. samples that are labeled with 0 in the provided class vector)")
  }
  perf.meas = tolower(perf.meas)
  perf.meas = switch(perf.meas, "auc"="AUC", "auprc"="AUPRC", "avgauc"="AvgAUC", "all"="All")
  if(!(perf.meas %in% c("AUC","AUPRC","AvgAUC","All"))){stop("perf.meas must be \"AUC\", \"AUPRC\", \"AvgAUC\", or \"All\"")}
  if(any(!c(pos.genes,neg.genes) %in% rownames(geneMtx))){stop("Some of the designated genes are not present in the provided gene matrix")}
  if(any(!c(init.pos,init.neg) %in% rownames(geneMtx))){stop("Some of the designated initial genes are not present in the provided gene matrix")}
  pos.genes = pos.genes[!pos.genes %in% init.pos]
  neg.genes = neg.genes[!neg.genes %in% init.neg]
  
  #set gene.combos
  if(is.null(bss.max)){
    bss.max=length(c(pos.genes,neg.genes))
  }
  if(bss.max > length(c(pos.genes,neg.genes))){
    warning("bss.max is greater than the total number of genes and will be ignored")
    bss.max = length(c(pos.genes,neg.genes))
  }
  if(bss.min > bss.max){stop("bss.min cannot be greater than bss.max")}
  if(bss.min < 0){stop("bss.min cannot be less than zero")}
  gene.combos = c(bss.max:bss.min)
  #figure out how many total combinations
  combo.count = 0
  all.length = length(c(pos.genes,neg.genes))
  for(i in gene.combos){combo.count=combo.count+choose(all.length,i)}
  cat(sprintf("Testing this many gene combinations: %s\nBe wary of values over 1 million; depending on the number of samples, this can take days to complete on a multi-core machine.\n",combo.count))
  
  require(ROCR)
  require(MetaIntegrator)
  require(parallel)
  require(data.table)
  progress.bar=TRUE
  if(progress.bar){
    if(progress.bar.type=="perGeneNum"){
      require(pbmcapply)
      bss.apply = match.fun(pbmclapply)
    }
    else if(progress.bar.type=="overall"){
      require(pbapply)
      bss.apply = match.fun(mclapply)
      total.num = length(pos.genes) + length(neg.genes)
      count = 0
      for(i in gene.combos){
        count = count + choose(total.num,i)
      }
      pb <- startpb(min=0,max=count)
      best_subset_pb_count = 1
    }
  }else{
    bss.apply = match.fun(mclapply)
  }
  
  if(perf.meas=="AvgAUC" || perf.meas=="All"){
    if(is.factor(labelVec)){labelVec=droplevels(labelVec)}
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
  }
  
  names(gene.combos)=paste(gene.combos,"Genes",sep=".")
  genesIndex = which(rownames(geneMtx) %in% c(pos.genes,neg.genes))
  BestSubsetList = lapply(gene.combos, function(numGenes){
    comb = combn(genesIndex,numGenes)
    comblist = split(comb, rep(1:ncol(comb), each = nrow(comb)))
    if(is.null(numCores)){
      #smart core usage = use 75% of resources unless you need less
      maxCores <- round(detectCores()*7.5/10)
      if(ncol(comb)<maxCores){
        maxCores <- ncol(comb)
      }
    }else{
      maxCores = numCores
    }
    geneCombinations = do.call(rbind,bss.apply(mc.cores=maxCores, comblist, function(chosenGenes){
      genes.curr = geneMtx[chosenGenes,,drop=F]
      pos = pos.genes[pos.genes %in% rownames(genes.curr)]
      if(!is.null(init.pos)){pos=c(pos,init.pos)}
      neg = neg.genes[neg.genes %in% rownames(genes.curr)]
      if(!is.null(init.neg)){neg=c(neg,init.neg)}
      init.length = length(c(init.pos,init.neg))
      scores = getGeneScores(genes.curr,pos,neg)
      
      if(progress.bar && progress.bar.type=="overall"){
        setpb(pb, best_subset_pb_count)
        best_subset_pb_count <<- best_subset_pb_count+1
      }
      
      if(perf.meas=="AUC"){
        MI_ROC = calculateROC(as.numeric(as.character(class)), as.numeric(scores))
        auc = MI_ROC$auc
        auc.lo = MI_ROC$auc.CI[1]
        auc.hi = min(MI_ROC$auc.CI[2],1)
        blurb=sprintf("%s Genes, AUC=%s (95%% CI %s-%s): %s up and %s down",numGenes+init.length,round(auc,3),round(auc.lo,3),round(auc.hi,3),paste(pos,collapse=", "),paste(neg,collapse=", "))
        return(c(numGenes+init.length,auc,auc.lo,auc.hi,paste(pos,collapse = " / "),paste(neg,collapse = " / "),blurb))
        
      }else if(perf.meas=="AUPRC"){
        pred_ROC = prediction(scores, class)
        perf = performance(pred_ROC, "prec", "rec")
        #add extra points at (1,0) and (0,1) to complete the plot
        perf@x.values[[1]] = c(0,perf@x.values[[1]],1)
        perf@y.values[[1]] = c(1,perf@y.values[[1]],0)
        auprc.list=.calcauprc(perf@x.values[[1]],perf@y.values[[1]],class)
        auprc=auprc.list$auprc
        auprc.lo = auprc.list$auprc.CI[1]
        auprc.hi = min(auprc.list$auprc.CI[2],1)
        blurb=sprintf("%s Genes, AUPRC=%s (95%% CI %s-%s): %s up and %s down",numGenes+init.length,round(auprc,3),round(auprc.lo,3),round(auprc.hi,3),paste(pos,collapse=", "),paste(neg,collapse=", "))
        return(c(numGenes+init.length,auprc,auprc.lo,auprc.hi,paste(pos,collapse = " / "),paste(neg,collapse = " / "),blurb))
        
      }else if(perf.meas=="AvgAUC"){
        aucVec=loVec=upVec=rep(0,length(classIndex))
        for(i in 1:length(classIndex)){
          MI_ROC = calculateROC(as.numeric(as.character(class[classIndex[[i]] ])), as.numeric(scores[classIndex[[i]] ]))
          aucVec[i] = MI_ROC$auc
          loVec[i] = max(MI_ROC$auc.CI[1],0)
          upVec[i] = min(MI_ROC$auc.CI[2],1)
        }
        avgauc=sum(aucVec*weights,na.rm = TRUE)
        avgauc.lo=sum(loVec*weights,na.rm = TRUE)
        avgauc.up=sum(upVec*weights,na.rm = TRUE)
        blurb=sprintf("%s Genes, Average AUC=%s (95%% CI %s-%s): %s up and %s down",numGenes+init.length,round(avgauc,3),round(avgauc.lo,3),round(avgauc.up,3),paste(pos,collapse=", "),paste(neg,collapse=", "))
        allAUCs=numeric(length(aucVec)*3)
        allAUCs[seq(1,by=3,length.out=length(aucVec))]=aucVec
        allAUCs[seq(2,by=3,length.out=length(loVec))]=loVec
        allAUCs[seq(3,by=3,length.out=length(upVec))]=upVec
        return(c(numGenes+init.length,avgauc,avgauc.lo,avgauc.up,paste(pos,collapse = " / "),paste(neg,collapse = " / "),blurb,allAUCs))
        
      }else if(perf.meas=="All"){
        MI_ROC = calculateROC(as.numeric(as.character(class)), as.numeric(scores))
        auc = MI_ROC$auc
        auc.lo = MI_ROC$auc.CI[1]
        auc.hi = min(MI_ROC$auc.CI[2],1)
        pred_ROC = prediction(scores, class)
        perf = performance(pred_ROC, "prec", "rec")
        perf@x.values[[1]] = c(0,perf@x.values[[1]],1)
        perf@y.values[[1]] = c(1,perf@y.values[[1]],0)
        auprc.list=.calcauprc(perf@x.values[[1]],perf@y.values[[1]],class)
        auprc=auprc.list$auprc
        auprc.lo = auprc.list$auprc.CI[1]
        auprc.hi = min(auprc.list$auprc.CI[2],1)
        aucVec=loVec=upVec=rep(0,length(classIndex))
        for(i in 1:length(classIndex)){
          MI_ROC = calculateROC(as.numeric(as.character(class[classIndex[[i]] ])), as.numeric(scores[classIndex[[i]] ]))
          aucVec[i] = MI_ROC$auc
          loVec[i] = max(MI_ROC$auc.CI[1],0)
          upVec[i] = min(MI_ROC$auc.CI[2],1)
        }
        avgauc=sum(aucVec*weights,na.rm = TRUE)
        avgauc.lo=sum(loVec*weights,na.rm = TRUE)
        avgauc.up=sum(upVec*weights,na.rm = TRUE)
        allAUCs=numeric(length(aucVec)*3)
        allAUCs[seq(1,by=3,length.out=length(aucVec))]=aucVec
        allAUCs[seq(2,by=3,length.out=length(loVec))]=loVec
        allAUCs[seq(3,by=3,length.out=length(upVec))]=upVec
        return(c(numGenes+init.length,auc,auc.lo,auc.hi,auprc,auprc.lo,auprc.hi,avgauc,avgauc.lo,avgauc.up,paste(pos,collapse = " / "),paste(neg,collapse = " / "),allAUCs))
      }else{
        stop("Invalid perf.meas")
      }
    }))
    
    if(perf.meas=="AUC" || perf.meas=="AUPRC" || perf.meas=="AvgAUC"){
      geneCombinations=data.frame(geneCombinations,stringsAsFactors=FALSE)
      colnames(geneCombinations)=c("numGenes",perf.meas,paste0(perf.meas,".ci.lower"),paste0(perf.meas,".ci.upper"),"upgenes","downgenes","blurb")
      geneCombinations$numGenes = as.numeric(geneCombinations$numGenes);geneCombinations[,2] = as.numeric(geneCombinations[,2])
      geneCombinations[,paste0(perf.meas,".ci.lower")] = as.numeric(geneCombinations[,paste0(perf.meas,".ci.lower")])
      geneCombinations[,paste0(perf.meas,".ci.upper")] = as.numeric(geneCombinations[,paste0(perf.meas,".ci.upper")])
      sortedCombinations = geneCombinations[order(geneCombinations[,2],decreasing=TRUE),]
      rownames(sortedCombinations) = c(1:dim(sortedCombinations)[1])
      if(perf.meas=="AvgAUC"){
        otherLabs=rep("",length(classIndex)*3)
        for(i in 1:length(classIndex)){
          otherLabs[(3*(i-1))+1]=paste(otherNames[i],"auc",sep=".")
          otherLabs[(3*(i-1))+2]=paste(otherNames[i],"ci.lower",sep=".")
          otherLabs[(3*(i-1))+3]=paste(otherNames[i],"ci.upper",sep=".")
          sortedCombinations[,5+(3*i)]=as.numeric(sortedCombinations[,5+(3*i)])
          sortedCombinations[,6+(3*i)]=as.numeric(sortedCombinations[,6+(3*i)])
          sortedCombinations[,7+(3*i)]=as.numeric(sortedCombinations[,7+(3*i)])
        }
        colnames(sortedCombinations)=append(c("numGenes","AvgAUC","AvgAUC.ci.lower","AvgAUC.ci.upper","upgenes","downgenes","blurb"),otherLabs)
      }
      cat(sprintf("Done with all %s-gene combinations%s\n",numGenes,
                  ifelse(length(c(init.pos,init.neg))>0,sprintf(" (with %s genes constant)",length(c(init.pos,init.neg))),"")))
      return(sortedCombinations)
      
    }else if(perf.meas=="All"){
      geneCombinations=data.frame(geneCombinations,stringsAsFactors=FALSE)
      numericCols=c(1:10, 13:(12+3*length(classIndex)))
      for(i in numericCols){
        geneCombinations[,i] = as.numeric(geneCombinations[,i])
      }
      if(sort.by=="AUC"){
        sortIndex=order(geneCombinations[,2],decreasing=TRUE)
      }else if(sort.by=="AUPRC"){
        sortIndex=order(geneCombinations[,5],decreasing=TRUE)
      }else if(sort.by=="AvgAUC"){
        sortIndex=order(geneCombinations[,8],decreasing=TRUE)
      }else{
        stop("Invalid sort.by")
      }
      sortedCombinations = geneCombinations[sortIndex,]
      rownames(sortedCombinations) = c(1:dim(sortedCombinations)[1])
      otherLabs=rep("",length(classIndex)*3)
      for(i in 1:length(classIndex)){
        otherLabs[(3*(i-1))+1]=paste(otherNames[i],"auc",sep=".")
        otherLabs[(3*(i-1))+2]=paste(otherNames[i],"ci.lower",sep=".")
        otherLabs[(3*(i-1))+3]=paste(otherNames[i],"ci.upper",sep=".")
      }
      colnames(sortedCombinations)=c("numGenes","AUC","AUC.ci.lower","AUC.ci.upper","AUPRC","AUPRC.ci.lower","AUPRC.ci.upper",
                                     "AvgAUC","AvgAUC.ci.lower","AvgAUC.ci.upper","upgenes","downgenes",otherLabs)
      
      if(!progress.bar || progress.bar.type!="overall"){
        cat(sprintf("Done with all %s-gene combinations%s\n",numGenes,
                    ifelse(length(c(init.pos,init.neg))>0,sprintf(" (with %s genes constant)",length(c(init.pos,init.neg))),"")))
      }
      return(sortedCombinations)
    }
  })
  
  bssDescription = list(all.upgenes=c(pos.genes,init.pos), all.downgenes=c(neg.genes,init.neg), perf.meas=perf.meas, bss.min=bss.min, bss.max=bss.max)
  if(!is.null(init.pos)){bssDescription$init.upgenes = init.pos}
  if(!is.null(init.neg)){bssDescription$init.downgenes = init.neg}
  if(perf.meas=="AvgAUC" || perf.meas=="All"){
    bssDescription$caseNames = caseNames
    bssDescription$otherNames = otherNames
    if(length(unique(weights))>1){
      bssDescription$weights = weights
    }
  }
  bssDT = data.table(do.call(rbind,BestSubsetList))
  BestSubsetList$bssDT = bssDT
  BestSubsetList$bssDescription = bssDescription
  
  if(progress.bar && progress.bar.type=="overall"){
    closepb(pb)
  }
  return(BestSubsetList)
}
