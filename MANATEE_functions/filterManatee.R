#' Filter out significant genes from MANATEE results
#' @description
#' After MANATEE has been run, and the results have been written to \code{pooledDataObject$manateeResults},
#' the genes can then be filtered using a variety of different statistical criteria. The user can set these criteria with the
#' various thresholds provided. Note that these thresholds vary slightly depending on the type of MANATEE analysis that was run.
#' @param pooledDataObject A MANATEE \code{pooledDataObject} with the results of a MANATEE analysis stored in \code{$manateeResults}
#' @param isLeaveOneOut If TRUE, then the designated thresholds will also be applied to the results of the leave-one-out analysis (Default: TRUE)
#' @param EffectSizeThresh A gene is selected if the absolute value of its Hedges' g effect size is greater than or equal to this threshold
#' @param ttestFDRThresh A gene is selected if the one-sided t-test FDR (using the t-test up for upgenes and the t-test down for downgenes) is less than or equal to this threshold
#' @param SAMscoreThresh A gene is selected if the absolute value of its SAM score is greater than or equal to this threshold
#' @param FoldChangeThresh A gene is selected if the absolute value of its fold change is greater than or equal to this threshold
#' @param NumStudiesThresh A gene is selected if the number of studies it is measured in is greater than or equal to this threshold
#' @param PropSamplesThresh A gene is selected if the proportion of samples it is measured in is greater than or equal to this threshold
#' @param NumTopGenesThresh Use this to select the top N genes across both the upgenes and downgenes, based on the specified \code{NumTopGenesMetric}
#' @param NumTopGenesMetric If \code{NumTopGenesThresh != 0}, then this specifies the metric used to select the top N genes
#' The options are "SAMscore", "EffectSize", "FoldChange", or "ttestFDR". If Bootstrapped MANATEE was run, "SAMscoreFDR", "EffectSizeFDR", or "FoldChangeFDR" can be used as well (Default: "EffectSize")
#' @param NumTopGenesStage Specifies when to select the top N genes. If \code{NumTopGenesStage}="pre-filter", then the top N genes are chosen before applying any of the other filters
#' If \code{NumTopGenesStage}="post-filter", then the top N genes are chosen after applying all of the other filters (Default: "post-filter")
#' @param consistencyFilter If TRUE, all genes with conflicting effect size, SAM scores, and/or fold change (i.e. some are negative and others are positive) will be removed.
#' Note that the majority of the conflicting genes are those with low effect sizes/SAM scores/fold changes anyway.
#' In addition, genes that have 0 for either numStudies1 or numStudies0 and genes with an effect size of 0 will be removed. (Default: TRUE)
#' @param returnType If \code{returnType}="pooledDataObject", then the function will return a modified version of the \code{pooledDataObject} with results stored in \code{pooledDataObject$filterResults}.
#' If \code{returnType}="filterObject", then the function will just return a \code{filterObject} that contains the results of the filtering (Default: "pooledDataObject")
#' @param SAM.FDRThresh \emph{Basic MANATEE:} A gene is selected if SAM returned an FDR for it that is less than or equal to this threshold
#' @param SAM.FDR.type \emph{Basic MANATEE:} Either "local" or "q-value". This determines whether to use local FDRs (which are computed for each gene using a local neighborhood of other genes) or q-values
#' (for which FDRs are calculated for all of the significant genes at each delta, and the q-value is the lowest FDR at which a gene was called significant) for the FDRThresh parameter (Default: "local")
#' @param topGeneProtection \emph{Basic MANATEE:} local FDRs are likely overestimated for the top genes, since the neighborhood of other genes used
#' to calculate FDR consists mostly of genes that are less significant than the top genes (see \url{https://statweb.stanford.edu/~tibs/SAM/sam.pdf}).
#' In order to account for this, if topGeneProtection=TRUE, then any genes with a q-value of 0 will be included, regardless of what their local FDR is (Default: TRUE)
#' @param EffectSizeFDRThresh \emph{Bootstrapped MANATEE:} A gene is selected if the one-sided effect size FDR
#' (using the effect size FDR up for upgenes and the effect size FDR down for downgenes) is less than or equal to this threshold
#' @param SAMscoreFDRThresh \emph{Bootstrapped MANATEE:} A gene is selected if the one-sided SAM score FDR
#' (using the SAM score FDR up for upgenes and the SAM score FDR down for downgenes) is less than or equal to this threshold
#' @param FoldChangeFDRThresh \emph{Bootstrapped MANATEE:} A gene is selected if the one-sided fold change FDR
#' (using the fold change FDR up for upgenes and the fold change FDR down for downgenes) is less than or equal to this threshold
#' @note The returned tables of filtered genes are sorted by effect size.
#' 
#' NumStudiesThresh and PropSamplesThresh are applied at the very beginning, even before the pre-filter NumTopGenes is applied
#' 
#' Note that LOO analysis uses the non-LOO results for NumTopGenes and consistencyFilter
#' @return The function either returns a modified version of the \code{pooledDataObject} with the results of the filtering stored
#' in \code{$filterResults} under a unique identifier, or it returns a \code{filterObject} containing the results of the filtering.
#' \item{upGeneNames}{character vector of the genes in the \code{upGeneTable}}
#' \item{downGeneNames}{character vector of the genes in the \code{downGeneTable}}
#' \item{upGeneTable}{dataframe containing the filtered upgenes, with gene names in rows and relevant statistics in columns}
#' \item{downGeneTable}{dataframe containing the filtered downgenes, with gene names in rows and relevant statistics in columns}
#' \item{filterDescription}{list with information that describes the results of the filtering and that catalogs the filtering options that were chosen}
#' @seealso	\code{\link{runManatee}}
#' @examples
#'
#' @export
#' @import data.table
#' @author Aditya Rao
filterManatee <- function(pooledDataObject, isLeaveOneOut=TRUE, EffectSizeThresh=0, ttestFDRThresh=0, SAMscoreThresh=0, FoldChangeThresh=0,
                          NumStudiesThresh=0, PropSamplesThresh=0, NumTopGenesThresh=0, NumTopGenesMetric="EffectSize", NumTopGenesStage="post-filter",
                          consistencyFilter=TRUE, returnType="pooledDataObject", SAM.FDRThresh=0, SAM.FDR.type="local", topGeneProtection=TRUE,
                          EffectSizeFDRThresh=0, SAMscoreFDRThresh=0, FoldChangeFDRThresh=0){
  if(is.null(pooledDataObject$manateeResults)){stop("manateeResults is missing from the provided pooledDataObject")}
  if(!checkManateeObject(pooledDataObject$manateeResults,"manateeResults")){stop("Invalid manateeResults provided")} #not sure i need this
  if(isLeaveOneOut && is.null(pooledDataObject$leaveOneOutAnalysis)){stop("leaveOneOutAnalysis is missing from the provided pooledDataObject")}
  if(returnType != "pooledDataObject" && returnType != "filterObject"){stop("returnType must be either \"pooledDataObject\" or \"filterObject\"")}
  if(NumTopGenesThresh != 0 && NumTopGenesStage != "pre-filter" && NumTopGenesStage != "post-filter"){
    stop("NumTopGenesStage must be either \"pre-filter\" or \"post-filter\"")
  }
  if(NumTopGenesStage=="pre-filter" && NumTopGenesThresh>(nrow(pooledDataObject$manateeResults$upgenes)+nrow(pooledDataObject$manateeResults$downgenes))){
    warning("NumTopGenesThresh is greater than the total number of genes and will be ignored")
    NumTopGenesThresh = 0
  }
  if(NumStudiesThresh > max(c(pooledDataObject$manateeResults$upgenes$numStudies,pooledDataObject$manateeResults$downgenes$numStudies))){
    stop("NumStudiesThresh is higher than the maximum numStudies value across all genes, and would thus eliminate all genes")
  }
  if(PropSamplesThresh > max(c(pooledDataObject$manateeResults$upgenes$propSamples,pooledDataObject$manateeResults$downgenes$propSamples))){
    stop("PropSamplesThresh is higher than the maximum propSamples value across all genes, and would thus eliminate all genes")
  }
  if(EffectSizeThresh < 0){stop("EffectSizeThresh cannot be negative")}
  if(ttestFDRThresh < 0){stop("ttestFDRThresh cannot be negative")}
  if(SAMscoreThresh < 0){stop("SAMscoreThresh cannot be negative")}
  if(FoldChangeThresh < 0){stop("FoldChangeThresh cannot be negative")}
  if(NumTopGenesThresh < 0){stop("NumTopGenesThresh cannot be negative")}
  if(NumStudiesThresh < 0){stop("NumStudiesThresh cannot be negative")}
  if(PropSamplesThresh < 0){stop("PropSamplesThresh cannot be negative")}
  
  require(data.table)
  
  if(pooledDataObject$manateeResults$type == "Basic MANATEE"){
    filterObject = .filterManateeBasic(pooledDataObject, isLeaveOneOut, EffectSizeThresh, ttestFDRThresh, SAMscoreThresh, FoldChangeThresh,
                                       NumStudiesThresh, PropSamplesThresh, NumTopGenesThresh, NumTopGenesMetric, NumTopGenesStage, consistencyFilter, returnType,
                                       SAM.FDRThresh=SAM.FDRThresh, SAM.FDR.type=SAM.FDR.type, topGeneProtection=topGeneProtection)
  }else if(pooledDataObject$manateeResults$type == "Bootstrapped MANATEE"){
    filterObject = .filterManateeBoot(pooledDataObject, isLeaveOneOut, EffectSizeThresh, ttestFDRThresh, SAMscoreThresh, FoldChangeThresh,
                                      NumStudiesThresh, PropSamplesThresh, NumTopGenesThresh, NumTopGenesMetric, NumTopGenesStage, consistencyFilter, returnType,
                                      EffectSizeFDRThresh=EffectSizeFDRThresh, SAMscoreFDRThresh=SAMscoreFDRThresh, FoldChangeFDRThresh=FoldChangeFDRThresh)
  }else if(pooledDataObject$manateeResults$type == "Pairwise MANATEE"){
    filterObject = .filterManateePairwise(pooledDataObject, isLeaveOneOut, EffectSizeThresh, ttestFDRThresh, SAMscoreThresh, FoldChangeThresh,
                                          NumStudiesThresh, PropSamplesThresh, NumTopGenesThresh, NumTopGenesMetric, NumTopGenesStage, consistencyFilter, returnType)
  }
  
  if(nrow(filterObject$upGeneTable)==0){
    warning("After filtering, there are no upgenes remaining")
  }
  if(nrow(filterObject$downGeneTable)==0){
    warning("After filtering, there are no downgenes remaining")
  }
  
  if(returnType == "filterObject"){
    return(filterObject)
  }else{ #if(returnType == "pooledDataObject")
    if("filterResults" %in% names(pooledDataObject)){
      pooledDataObject$filterResults[[filterObject$filterDescription$filterName]]=filterObject
    }else{
      pooledDataObject$filterResults = list()
      pooledDataObject$filterResults[[filterObject$filterDescription$filterName]]=filterObject
    }
    return(pooledDataObject)
  }
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#--------------------------------------------------------Basic MANATEE--------------------------------------------------------#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

###-###-###-###-###-###-###-###-#
###   .filterManateeBasic()   ###
###-###-###-###-###-###-###-###-#

.filterManateeBasic <- function(pooledDataObject, isLeaveOneOut=TRUE, EffectSizeThresh=0, ttestFDRThresh=0, SAMscoreThresh=0, FoldChangeThresh=0,
                                NumStudiesThresh=0, PropSamplesThresh=0, NumTopGenesThresh=0, NumTopGenesMetric="EffectSize", NumTopGenesStage="post-filter",
                                consistencyFilter=TRUE, returnType="pooledDataObject", SAM.FDRThresh=0, SAM.FDR.type="local", topGeneProtection=TRUE){
  upgenes = data.table(pooledDataObject$manateeResults$upgenes, keep.rownames = T)
  downgenes = data.table(pooledDataObject$manateeResults$downgenes, keep.rownames = T)
  
  if(SAM.FDR.type != "local" && SAM.FDR.type != "q-value"){stop("SAM.FDR.type must be either \"local\" or \"q-value\"")}
  if(SAM.FDR.type == "local" && !("SAM.localFDR" %in% c(colnames(upgenes),colnames(downgenes)))){
    stop("If SAM.FDR.type=\"local\", then local FDR must be calculated in your manateeResults")
  }
  if(NumTopGenesThresh != 0 && !NumTopGenesMetric %in% c("SAMscore", "EffectSize", "FoldChange", "ttestFDR")){
    stop("For Basic MANATEE, NumTopGenesMetric must be \"SAMscore\", \"EffectSize\", \"FoldChange\", or \"ttestFDR\"")
  }
  if(SAM.FDRThresh < 0){stop("SAM.FDRThresh cannot be negative")}
  
  if(NumStudiesThresh != 0){
    upgenes = upgenes[numStudies >= NumStudiesThresh]
    downgenes = downgenes[numStudies >= NumStudiesThresh]
  }
  
  if(PropSamplesThresh != 0){
    upgenes = upgenes[propSamples >= PropSamplesThresh]
    downgenes = downgenes[propSamples >= PropSamplesThresh]
  }
  
  if(NumTopGenesThresh != 0 && NumTopGenesStage=="pre-filter"){
    topgenes = .getNumTopGenesBasic(upgenes, downgenes, NumTopGenesThresh, NumTopGenesMetric)
    upgenes = topgenes$upgenes
    downgenes = topgenes$downgenes
  }
  
  if(isLeaveOneOut){
    updown.loo = lapply(pooledDataObject$leaveOneOutAnalysis,function(data){
      upgenes.loo = data.table(data$upgenes, keep.rownames = T)
      downgenes.loo = data.table(data$downgenes, keep.rownames = T)
      return(.filterManateeBasicCore(upgenes.loo, downgenes.loo, EffectSizeThresh=EffectSizeThresh, ttestFDRThresh=ttestFDRThresh,
                                     SAMscoreThresh=SAMscoreThresh, FoldChangeThresh=FoldChangeThresh, SAM.FDRThresh=SAM.FDRThresh,
                                     SAM.FDR.type=SAM.FDR.type, topGeneProtection=topGeneProtection))
    })
    up.names = upgenes$rn
    down.names = downgenes$rn
    for(i in 1:length(updown.loo)){
      up.names = intersect(up.names, updown.loo[[i]]$upgenes$rn)
      down.names = intersect(down.names, updown.loo[[i]]$downgenes$rn)
    }
    
    upgenes = upgenes[rn %in% up.names]
    downgenes = downgenes[rn %in% down.names]
  }else{
    updown = .filterManateeBasicCore(upgenes, downgenes, EffectSizeThresh=EffectSizeThresh, ttestFDRThresh=ttestFDRThresh,
                                     SAMscoreThresh=SAMscoreThresh, FoldChangeThresh=FoldChangeThresh, SAM.FDRThresh=SAM.FDRThresh,
                                     SAM.FDR.type=SAM.FDR.type, topGeneProtection=topGeneProtection)
    upgenes = updown$upgenes
    downgenes = updown$downgenes
  }
  
  if(consistencyFilter){
    upgenes = upgenes[(SAM.score*effectSize)>=0]
    downgenes = downgenes[(SAM.score*effectSize)>=0]
    upgenes = upgenes[(SAM.score*foldChange)>=0]
    downgenes = downgenes[(SAM.score*foldChange)>=0]
    upgenes = upgenes[effectSize!=0]
    downgenes = downgenes[effectSize!=0]
    upgenes = upgenes[numStudies1!=0]
    downgenes = downgenes[numStudies1!=0]
    upgenes = upgenes[numStudies0!=0]
    downgenes = downgenes[numStudies0!=0]
  }
  
  if(NumTopGenesStage=="post-filter" && NumTopGenesThresh>(nrow(upgenes)+nrow(downgenes))){
    warning("After filtering, NumTopGenesThresh is now greater than the total number of genes and will be ignored")
    NumTopGenesThresh = 0
  }
  if(NumTopGenesThresh != 0 && NumTopGenesStage=="post-filter"){
    topgenes = .getNumTopGenesBasic(upgenes, downgenes, NumTopGenesThresh, NumTopGenesMetric)
    upgenes = topgenes$upgenes
    downgenes = topgenes$downgenes
  }
  
  #sort tables
  upgenes = upgenes[order(effectSize,decreasing=T)]
  downgenes = downgenes[order(effectSize,decreasing=F)]
  
  upGeneNames = upgenes$rn
  downGeneNames = downgenes$rn
  upGeneTable = data.frame(upgenes)
  rownames(upGeneTable) = upGeneTable$rn
  upGeneTable$rn = NULL
  downGeneTable = data.frame(downgenes)
  rownames(downGeneTable) = downGeneTable$rn
  downGeneTable$rn = NULL
  
  #make filterName
  fname0 = ifelse(isLeaveOneOut,"loo_","")
  fname1 = paste(sep = "_", ifelse(EffectSizeThresh != 0, paste0("ES",EffectSizeThresh), ""),
                 ifelse(ttestFDRThresh != 0, paste0("tFDR",ttestFDRThresh), ""),
                 ifelse(SAMscoreThresh != 0, paste0("SAM",SAMscoreThresh), ""),
                 ifelse(FoldChangeThresh != 0, paste0("FC",FoldChangeThresh), ""),
                 ifelse(SAM.FDRThresh != 0, paste0("sFDR",ifelse(SAM.FDR.type=="local",".l",".q"),SAM.FDRThresh), ""),
                 ifelse(NumStudiesThresh != 0, paste0("nS",NumStudiesThresh), ""),
                 ifelse(PropSamplesThresh != 0, paste0("pS",PropSamplesThresh), ""),"")
  if(NumTopGenesThresh != 0){
    metric = switch(NumTopGenesMetric,"SAMscore"="s","EffectSize"="e","FoldChange"="f","ttestFDR"="t")
    fname2 = sprintf("nTop.%s.%s%s_",metric,ifelse(NumTopGenesStage=="post-filter","po","pr"),NumTopGenesThresh)
  }else{
    fname2 = ""
  }
  fname3 = paste(ifelse(topGeneProtection,"tgp",""),ifelse(consistencyFilter,"cf",""),sep="_")
  filterName = paste0(fname0,fname1,fname2,fname3)
  
  #get rid of extra underscores
  filterName = gsub("__+","_",filterName)
  filterName = gsub("^_","",filterName)
  filterName = gsub("_$","",filterName)
  
  #format filterObject
  filterDescription = list(n.upgenes = nrow(upgenes),n.downgenes = nrow(downgenes),
                           EffectSizeThresh = EffectSizeThresh,ttestFDRThresh = ttestFDRThresh,
                           SAMscoreThresh = SAMscoreThresh,FoldChangeThresh = FoldChangeThresh,
                           NumStudiesThresh=NumStudiesThresh, PropSamplesThresh=PropSamplesThresh,
                           NumTopGenesThresh = NumTopGenesThresh,SAM.FDRThresh = SAM.FDRThresh,
                           SAM.FDR.type = SAM.FDR.type,NumTopGenesMetric = NumTopGenesMetric,
                           NumTopGenesStage = NumTopGenesStage,topGeneProtection = topGeneProtection,
                           consistencyFilter = consistencyFilter,filterName = filterName,
                           version = "MANATEE 1.0",timestamp = Sys.time())
  
  filterObject = list(upGeneNames = upGeneNames,downGeneNames = downGeneNames,
                      upGeneTable = upGeneTable,downGeneTable = downGeneTable,
                      filterDescription = filterDescription)
  
  return(filterObject)
}



###-###-###-###-###-###-###-###-###-#
###   .filterManateeBasicCore()   ###
###-###-###-###-###-###-###-###-###-#

.filterManateeBasicCore <- function(upgenes, downgenes, EffectSizeThresh=0, ttestFDRThresh=0, SAMscoreThresh=0,
                                    FoldChangeThresh=0,SAM.FDRThresh=0, SAM.FDR.type="local", topGeneProtection=TRUE){
  if(EffectSizeThresh != 0){
    upgenes = upgenes[effectSize >= EffectSizeThresh]
    downgenes = downgenes[effectSize <= -EffectSizeThresh]
  }
  
  if(ttestFDRThresh != 0){
    upgenes = upgenes[ttestFDRUp <= ttestFDRThresh]
    downgenes = downgenes[ttestFDRDown <= ttestFDRThresh]
  }
  
  if(SAMscoreThresh != 0){
    upgenes = upgenes[SAM.score >= SAMscoreThresh]
    downgenes = downgenes[SAM.score <= -SAMscoreThresh]
  }
  
  if(FoldChangeThresh != 0){
    upgenes = upgenes[foldChange >= FoldChangeThresh]
    downgenes = downgenes[foldChange <= -FoldChangeThresh]
  }
  
  if(SAM.FDRThresh != 0){
    if(SAM.FDR.type == "local"){
      if(topGeneProtection){
        upgenes = upgenes[SAM.localFDR <= SAM.FDRThresh | SAM.qValue == 0]
        downgenes = downgenes[SAM.localFDR <= SAM.FDRThresh | SAM.qValue == 0]
      }else{
        upgenes = upgenes[SAM.localFDR <= SAM.FDRThresh]
        downgenes = downgenes[SAM.localFDR <= SAM.FDRThresh]
      }
    }else{ #if SAM.FDR.type == "q-value"
      upgenes = upgenes[SAM.qValue <= SAM.FDRThresh]
      downgenes = downgenes[SAM.qValue <= SAM.FDRThresh]
    }
  }
  
  return(list(upgenes=upgenes,downgenes=downgenes))
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#--------------------------------------------------------Bootstrapped MANATEE--------------------------------------------------------#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

###-###-###-###-###-###-###-##-#
###   .filterManateeBoot()   ###
###-###-###-###-###-###-###-##-#

.filterManateeBoot <- function(pooledDataObject, isLeaveOneOut=TRUE, EffectSizeThresh=0, ttestFDRThresh=0, SAMscoreThresh=0, FoldChangeThresh=0,
                               NumStudiesThresh=0, PropSamplesThresh=0, NumTopGenesThresh=0,NumTopGenesMetric="EffectSize",NumTopGenesStage="post-filter",consistencyFilter=TRUE,
                               returnType="pooledDataObject", EffectSizeFDRThresh=0, SAMscoreFDRThresh=0, FoldChangeFDRThresh=0){
  if(NumTopGenesThresh != 0 && !NumTopGenesMetric %in% c("SAMscore", "EffectSize", "FoldChange", "ttestFDR", "SAMscoreFDR", "EffectSizeFDR", "FoldChangeFDR")){
    stop("For Bootstrapped MANATEE, NumTopGenesMetric must be \"SAMscore\", \"EffectSize\",
         \"FoldChange\", \"ttestFDR\", \"SAMscoreFDR\", \"EffectSizeFDR\", or \"FoldChangeFDR\"")
  }
  if(EffectSizeFDRThresh < 0){stop("EffectSizeFDRThresh cannot be negative")}
  if(SAMscoreFDRThresh < 0){stop("SAMscoreFDRThresh cannot be negative")}
  if(FoldChangeFDRThresh < 0){stop("FoldChangeFDRThresh cannot be negative")}
  
  upgenes = data.table(pooledDataObject$manateeResults$upgenes, keep.rownames = T)
  downgenes = data.table(pooledDataObject$manateeResults$downgenes, keep.rownames = T)
  
  if(NumStudiesThresh != 0){
    upgenes = upgenes[numStudies >= NumStudiesThresh]
    downgenes = downgenes[numStudies >= NumStudiesThresh]
  }
  
  if(PropSamplesThresh != 0){
    upgenes = upgenes[propSamples >= PropSamplesThresh]
    downgenes = downgenes[propSamples >= PropSamplesThresh]
  }
  
  if(NumTopGenesThresh != 0 && NumTopGenesStage=="pre-filter"){
    topgenes = .getNumTopGenesBoot(upgenes, downgenes, NumTopGenesThresh, NumTopGenesMetric)
    upgenes = topgenes$upgenes
    downgenes = topgenes$downgenes
  }
  
  if(isLeaveOneOut){
    updown.loo = lapply(pooledDataObject$leaveOneOutAnalysis,function(data){
      upgenes.loo = data.table(data$upgenes, keep.rownames = T)
      downgenes.loo = data.table(data$downgenes, keep.rownames = T)
      return(.filterManateeBootCore(upgenes.loo, downgenes.loo, EffectSizeThresh=EffectSizeThresh, ttestFDRThresh=ttestFDRThresh,
                                    SAMscoreThresh=SAMscoreThresh, FoldChangeThresh=FoldChangeThresh, EffectSizeFDRThresh=EffectSizeFDRThresh,
                                    SAMscoreFDRThresh=SAMscoreFDRThresh, FoldChangeFDRThresh=FoldChangeFDRThresh))
    })
    up.names = upgenes$rn
    down.names = downgenes$rn
    for(i in 1:length(updown.loo)){
      up.names = intersect(up.names, updown.loo[[i]]$upgenes$rn)
      down.names = intersect(down.names, updown.loo[[i]]$downgenes$rn)
    }
    
    upgenes = upgenes[rn %in% up.names]
    downgenes = downgenes[rn %in% down.names]
  }else{
    updown = .filterManateeBootCore(upgenes, downgenes, EffectSizeThresh=EffectSizeThresh, ttestFDRThresh=ttestFDRThresh,
                                    SAMscoreThresh=SAMscoreThresh, FoldChangeThresh=FoldChangeThresh, EffectSizeFDRThresh=EffectSizeFDRThresh,
                                    SAMscoreFDRThresh=SAMscoreFDRThresh, FoldChangeFDRThresh=FoldChangeFDRThresh)
    upgenes = updown$upgenes
    downgenes = updown$downgenes
  }
  
  if(consistencyFilter){
    upgenes = upgenes[(boot.SAM.score*boot.effectSize)>=0]
    downgenes = downgenes[(boot.SAM.score*boot.effectSize)>=0]
    upgenes = upgenes[(boot.SAM.score*boot.foldChange)>=0]
    downgenes = downgenes[(boot.SAM.score*boot.foldChange)>=0]
    upgenes = upgenes[boot.effectSize!=0]
    downgenes = downgenes[boot.effectSize!=0]
    upgenes = upgenes[numStudies1!=0]
    downgenes = downgenes[numStudies1!=0]
    upgenes = upgenes[numStudies0!=0]
    downgenes = downgenes[numStudies0!=0]
  }
  
  if(NumTopGenesStage=="post-filter" && NumTopGenesThresh>(nrow(upgenes)+nrow(downgenes))){
    warning("After filtering, NumTopGenesThresh is now greater than the total number of genes and will be ignored")
    NumTopGenesThresh = 0
  }
  if(NumTopGenesThresh != 0 && NumTopGenesStage=="post-filter"){
    topgenes = .getNumTopGenesBoot(upgenes, downgenes, NumTopGenesThresh, NumTopGenesMetric)
    upgenes = topgenes$upgenes
    downgenes = topgenes$downgenes
  }
  
  #sort tables
  upgenes = upgenes[order(boot.effectSize,decreasing=T)]
  downgenes = downgenes[order(boot.effectSize,decreasing=F)]
  
  upGeneNames = upgenes$rn
  downGeneNames = downgenes$rn
  upGeneTable = data.frame(upgenes)
  rownames(upGeneTable) = upGeneTable$rn
  upGeneTable$rn = NULL
  downGeneTable = data.frame(downgenes)
  rownames(downGeneTable) = downGeneTable$rn
  downGeneTable$rn = NULL

  #make filterName
  fname0 = ifelse(isLeaveOneOut,"loo_","")
  fname1 = paste(sep = "_", ifelse(EffectSizeThresh != 0, paste0("ES",EffectSizeThresh), ""),
                 ifelse(ttestFDRThresh != 0, paste0("tFDR",ttestFDRThresh), ""),
                 ifelse(SAMscoreThresh != 0, paste0("SAM",SAMscoreThresh), ""),
                 ifelse(FoldChangeThresh != 0, paste0("FC",FoldChangeThresh), ""),
                 ifelse(EffectSizeFDRThresh != 0, paste0("eFDR",EffectSizeFDRThresh), ""),
                 ifelse(SAMscoreFDRThresh != 0, paste0("sFDR",SAMscoreFDRThresh), ""),
                 ifelse(FoldChangeFDRThresh != 0, paste0("fFDR",FoldChangeFDRThresh), ""),
                 ifelse(NumStudiesThresh != 0, paste0("nS",NumStudiesThresh), ""),
                 ifelse(PropSamplesThresh != 0, paste0("pS",PropSamplesThresh), ""), "")
  if(NumTopGenesThresh != 0){
    metric = switch(NumTopGenesMetric,"SAMscore"="s","EffectSize"="e","FoldChange"="f","ttestFDR"="t",
                    "SAMscoreFDR"="sf", "EffectSizeFDR"="sf", "FoldChangeFDR"="ff")
    fname2 = sprintf("nTop.%s.%s%s_",metric,ifelse(NumTopGenesStage=="post-filter","po","pr"),NumTopGenesThresh)
  }else{
    fname2 = ""
  }
  fname3 = ifelse(consistencyFilter,"cf","")
  filterName = paste0(fname0,fname1,fname2,fname3)

  #get rid of extra underscores
  filterName = gsub("__+","_",filterName)
  filterName = gsub("^_","",filterName)
  filterName = gsub("_$","",filterName)
  
  #format filterObject
  filterDescription = list(n.upgenes = nrow(upgenes),n.downgenes = nrow(downgenes),
                           EffectSizeThresh = EffectSizeThresh,ttestFDRThresh = ttestFDRThresh,
                           SAMscoreThresh = SAMscoreThresh,FoldChangeThresh = FoldChangeThresh,
                           NumStudiesThresh=NumStudiesThresh, PropSamplesThresh=PropSamplesThresh,
                           EffectSizeFDRThresh = EffectSizeFDRThresh,SAMscoreFDRThresh = SAMscoreFDRThresh,
                           FoldChangeFDRThresh = EffectSizeFDRThresh,NumTopGenesThresh = NumTopGenesThresh,
                           NumTopGenesMetric = NumTopGenesMetric,NumTopGenesStage = NumTopGenesStage,
                           consistencyFilter = consistencyFilter,filterName = filterName,
                           version = "MANATEE 1.0",timestamp = Sys.time())
  
  filterObject = list(upGeneNames = upGeneNames,downGeneNames = downGeneNames,
                      upGeneTable = upGeneTable,downGeneTable = downGeneTable,
                      filterDescription = filterDescription)
  
  return(filterObject)
}



###-###-###-###-###-###-###-###-##-#
###   .filterManateeBootCore()   ###
###-###-###-###-###-###-###-###-##-#

.filterManateeBootCore <- function(upgenes, downgenes, EffectSizeThresh=0, ttestFDRThresh=0, SAMscoreThresh=0,
                                   FoldChangeThresh=0, EffectSizeFDRThresh=0, SAMscoreFDRThresh=0, FoldChangeFDRThresh=0){
  if(EffectSizeThresh != 0){
    upgenes = upgenes[boot.effectSize >= EffectSizeThresh]
    downgenes = downgenes[boot.effectSize <= -EffectSizeThresh]
  }
  
  if(ttestFDRThresh != 0){
    upgenes = upgenes[original.ttestFDRUp <= ttestFDRThresh]
    downgenes = downgenes[original.ttestFDRDown <= ttestFDRThresh]
  }
  
  if(SAMscoreThresh != 0){
    upgenes = upgenes[boot.SAM.score >= SAMscoreThresh]
    downgenes = downgenes[boot.SAM.score <= -SAMscoreThresh]
  }
  
  if(FoldChangeThresh != 0){
    upgenes = upgenes[boot.foldChange >= FoldChangeThresh]
    downgenes = downgenes[boot.foldChange <= -FoldChangeThresh]
  }
  
  if(EffectSizeFDRThresh != 0){
    upgenes = upgenes[boot.effectSizeFDRUp <= EffectSizeFDRThresh]
    downgenes = downgenes[boot.effectSizeFDRDown <= EffectSizeFDRThresh]
  }
  
  if(SAMscoreFDRThresh != 0){
    upgenes = upgenes[boot.SAM.scoreFDRUp <= SAMscoreFDRThresh]
    downgenes = downgenes[boot.SAM.scoreFDRDown <= SAMscoreFDRThresh]
  }
  
  if(FoldChangeFDRThresh != 0){
    upgenes = upgenes[boot.foldChangeFDRUp <= FoldChangeFDRThresh]
    downgenes = downgenes[boot.foldChangeFDRDown <= FoldChangeFDRThresh]
  }
  
  return(list(upgenes=upgenes,downgenes=downgenes))
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#--------------------------------------------------------Pairwise MANATEE--------------------------------------------------------#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

###-###-###-###-###-###-###-###-##-#
###   .filterManateePairwise()   ###
###-###-###-###-###-###-###-###-##-#

.filterManateePairwise <- function(pooledDataObject, isLeaveOneOut=TRUE, EffectSizeThresh=0, ttestFDRThresh=0, SAMscoreThresh=0, FoldChangeThresh=0,
                                   NumStudiesThresh=0, PropSamplesThresh=0, NumTopGenesThresh=0, NumTopGenesMetric="EffectSize", NumTopGenesStage="post-filter", consistencyFilter=TRUE, returnType="pooledDataObject"){
  if(NumTopGenesThresh != 0 && !NumTopGenesMetric %in% c("SAMscore", "EffectSize", "FoldChange", "ttestFDR")){
    stop("For Pairwise MANATEE, NumTopGenesMetric must be \"SAMscore\", \"EffectSize\", \"FoldChange\", or \"ttestFDR\"")
  }
  
  upgenes = data.table(pooledDataObject$manateeResults$upgenes, keep.rownames = T)
  downgenes = data.table(pooledDataObject$manateeResults$downgenes, keep.rownames = T)
  
  if(NumStudiesThresh != 0){
    upgenes = upgenes[numStudies >= NumStudiesThresh]
    downgenes = downgenes[numStudies >= NumStudiesThresh]
  }
  
  if(PropSamplesThresh != 0){
    upgenes = upgenes[propSamples >= PropSamplesThresh]
    downgenes = downgenes[propSamples >= PropSamplesThresh]
  }
  
  if(NumTopGenesThresh != 0 && NumTopGenesStage=="pre-filter"){
    topgenes = .getNumTopGenesPairwise(upgenes, downgenes, NumTopGenesThresh, NumTopGenesMetric)
    upgenes = topgenes$upgenes
    downgenes = topgenes$downgenes
  }
  
  if(isLeaveOneOut){
    updown.loo = lapply(pooledDataObject$leaveOneOutAnalysis,function(data){
      upgenes.loo = data.table(data$upgenes, keep.rownames = T)
      downgenes.loo = data.table(data$downgenes, keep.rownames = T)
      return(.filterManateePairwiseCore(upgenes.loo, downgenes.loo, EffectSizeThresh=EffectSizeThresh, ttestFDRThresh=ttestFDRThresh,
                                        SAMscoreThresh=SAMscoreThresh, FoldChangeThresh=FoldChangeThresh))
    })
    up.names = upgenes$rn
    down.names = downgenes$rn
    for(i in 1:length(updown.loo)){
      up.names = intersect(up.names, updown.loo[[i]]$upgenes$rn)
      down.names = intersect(down.names, updown.loo[[i]]$downgenes$rn)
    }
    
    upgenes = upgenes[rn %in% up.names]
    downgenes = downgenes[rn %in% down.names]
  }else{
    updown = .filterManateePairwiseCore(upgenes, downgenes, EffectSizeThresh=EffectSizeThresh, ttestFDRThresh=ttestFDRThresh,
                                        SAMscoreThresh=SAMscoreThresh, FoldChangeThresh=FoldChangeThresh)
    upgenes = updown$upgenes
    downgenes = updown$downgenes
  }
  
  if(consistencyFilter){
    upgenes = upgenes[(pair.SAM.score*pair.effectSize)>=0]
    downgenes = downgenes[(pair.SAM.score*pair.effectSize)>=0]
    upgenes = upgenes[(pair.SAM.score*pair.foldChange)>=0]
    downgenes = downgenes[(pair.SAM.score*pair.foldChange)>=0]
    upgenes = upgenes[pair.effectSize!=0]
    downgenes = downgenes[pair.effectSize!=0]
    upgenes = upgenes[numStudies1!=0]
    downgenes = downgenes[numStudies1!=0]
    upgenes = upgenes[numStudies0!=0]
    downgenes = downgenes[numStudies0!=0]
  }
  
  if(NumTopGenesStage=="post-filter" && NumTopGenesThresh>(nrow(upgenes)+nrow(downgenes))){
    warning("After filtering, NumTopGenesThresh is now greater than the total number of genes and will be ignored")
    NumTopGenesThresh = 0
  }
  if(NumTopGenesThresh != 0 && NumTopGenesStage=="post-filter"){
    topgenes = .getNumTopGenesPairwise(upgenes, downgenes, NumTopGenesThresh, NumTopGenesMetric)
    upgenes = topgenes$upgenes
    downgenes = topgenes$downgenes
  }
  
  #sort tables
  upgenes = upgenes[order(pair.effectSize,decreasing=T)]
  downgenes = downgenes[order(pair.effectSize,decreasing=F)]
  
  upGeneNames = upgenes$rn
  downGeneNames = downgenes$rn
  upGeneTable = data.frame(upgenes)
  rownames(upGeneTable) = upGeneTable$rn
  upGeneTable$rn = NULL
  downGeneTable = data.frame(downgenes)
  rownames(downGeneTable) = downGeneTable$rn
  downGeneTable$rn = NULL
  
  #make filterName
  fname0 = ifelse(isLeaveOneOut,"loo_","")
  fname1 = paste(sep = "_", ifelse(EffectSizeThresh != 0, paste0("ES",EffectSizeThresh), ""),
                 ifelse(ttestFDRThresh != 0, paste0("tFDR",ttestFDRThresh), ""),
                 ifelse(SAMscoreThresh != 0, paste0("SAM",SAMscoreThresh), ""),
                 ifelse(FoldChangeThresh != 0, paste0("FC",FoldChangeThresh), ""),
                 ifelse(NumStudiesThresh != 0, paste0("nS",NumStudiesThresh), ""),
                 ifelse(PropSamplesThresh != 0, paste0("pS",PropSamplesThresh), ""), "")
  if(NumTopGenesThresh != 0){
    metric = switch(NumTopGenesMetric,"SAMscore"="s","EffectSize"="e","FoldChange"="f","ttestFDR"="t")
    fname2 = sprintf("nTop.%s.%s%s_",metric,ifelse(NumTopGenesStage=="post-filter","po","pr"),NumTopGenesThresh)
  }else{
    fname2 = ""
  }
  fname3 = ifelse(consistencyFilter,"cf","")
  filterName = paste0(fname0,fname1,fname2,fname3)

  #get rid of extra underscores
  filterName = gsub("__+","_",filterName)
  filterName = gsub("^_","",filterName)
  filterName = gsub("_$","",filterName)
  
  #format filterObject
  filterDescription = list(n.upgenes = nrow(upgenes),n.downgenes = nrow(downgenes),
                           EffectSizeThresh = EffectSizeThresh,ttestFDRThresh = ttestFDRThresh,
                           SAMscoreThresh = SAMscoreThresh,FoldChangeThresh = FoldChangeThresh,
                           NumStudiesThresh=NumStudiesThresh, PropSamplesThresh=PropSamplesThresh,
                           NumTopGenesThresh = NumTopGenesThresh,NumTopGenesMetric = NumTopGenesMetric,
                           NumTopGenesStage = NumTopGenesStage,consistencyFilter = consistencyFilter,
                           filterName = filterName,version = "MANATEE 1.0",timestamp = Sys.time())
  
  filterObject = list(upGeneNames = upGeneNames,downGeneNames = downGeneNames,
                      upGeneTable = upGeneTable,downGeneTable = downGeneTable,
                      filterDescription = filterDescription)
  
  return(filterObject)
}



###-###-###-###-###-###-###-###-###-##-#
###   .filterManateePairwiseCore()   ###
###-###-###-###-###-###-###-###-###-##-#

.filterManateePairwiseCore <- function(upgenes, downgenes, EffectSizeThresh=0, ttestFDRThresh=0, SAMscoreThresh=0, FoldChangeThresh=0){
  if(EffectSizeThresh != 0){
    upgenes = upgenes[pair.effectSize >= EffectSizeThresh]
    downgenes = downgenes[pair.effectSize <= -EffectSizeThresh]
  }
  
  if(ttestFDRThresh != 0){
    upgenes = upgenes[pair.ttestFDRUp <= ttestFDRThresh]
    downgenes = downgenes[pair.ttestFDRDown <= ttestFDRThresh]
  }
  
  if(SAMscoreThresh != 0){
    upgenes = upgenes[pair.SAM.score >= SAMscoreThresh]
    downgenes = downgenes[pair.SAM.score <= -SAMscoreThresh]
  }
  
  if(FoldChangeThresh != 0){
    upgenes = upgenes[pair.foldChange >= FoldChangeThresh]
    downgenes = downgenes[pair.foldChange <= -FoldChangeThresh]
  }
  
  return(list(upgenes=upgenes,downgenes=downgenes))
}



#################################################################################################################
#################################################################################################################
#################################################################################################################



###-###-###-###-###-###-###-###-##
###   .getNumTopGenesBasic()   ###
###-###-###-###-###-###-###-###-##
#gets the top N genes (across both up and downgenes) based on the metric the user chooses
#this is the version for use with Basic MANATEE

.getNumTopGenesBasic <- function(upgenes, downgenes, NumTopGenesThresh, NumTopGenesMetric){
  if(NumTopGenesMetric == "SAMscore"){
    all.up = abs(upgenes$SAM.score)
    all.down = abs(downgenes$SAM.score)
  }
  if(NumTopGenesMetric == "EffectSize"){
    all.up = abs(upgenes$effectSize)
    all.down = abs(downgenes$effectSize)
  }
  if(NumTopGenesMetric == "FoldChange"){
    all.up = abs(upgenes$foldChange)
    all.down = abs(downgenes$foldChange)
  }
  if(NumTopGenesMetric == "ttestFDR"){
    all.up = upgenes$ttestFDRUp
    all.down = downgenes$ttestFDRDown
  }
  names(all.up) = upgenes$rn
  names(all.down) = downgenes$rn
  all = c(all.up,all.down)
  if(NumTopGenesMetric %in% c("SAMscore","EffectSize","FoldChange")){
    top.genes = names(sort(all,decreasing=T))[1:NumTopGenesThresh]
  }else{ #if NumTopGenesMetric=="ttestFDR"
    top.genes = names(sort(all,decreasing=F))[1:NumTopGenesThresh]
  }
  upgenes = upgenes[rn %in% top.genes]
  downgenes = downgenes[rn %in% top.genes]
  return(list(upgenes = upgenes,downgenes = downgenes))
}



###-###-###-###-###-###-###-###-#
###   .getNumTopGenesBoot()   ###
###-###-###-###-###-###-###-###-#
#gets the top N genes (across both up and downgenes) based on the metric the user chooses
#this is the version for use with Bootstrapped MANATEE

.getNumTopGenesBoot <- function(upgenes, downgenes, NumTopGenesThresh, NumTopGenesMetric){
  if(NumTopGenesMetric == "SAMscore"){
    all.up = abs(upgenes$boot.SAM.score)
    all.down = abs(downgenes$boot.SAM.score)
  }
  if(NumTopGenesMetric == "EffectSize"){
    all.up = abs(upgenes$boot.effectSize)
    all.down = abs(downgenes$boot.effectSize)
  }
  if(NumTopGenesMetric == "FoldChange"){
    all.up = abs(upgenes$boot.foldChange)
    all.down = abs(downgenes$boot.foldChange)
  }
  if(NumTopGenesMetric == "ttestFDR"){
    all.up = upgenes$original.ttestFDRUp
    all.down = downgenes$original.ttestFDRDown
  }
  if(NumTopGenesMetric == "SAMscoreFDR"){
    all.up = upgenes$boot.SAM.scoreFDRUp
    all.down = downgenes$boot.SAM.scoreFDRDown
  }
  if(NumTopGenesMetric == "EffectSizeFDR"){
    all.up = upgenes$boot.effectSizeFDRUp
    all.down = downgenes$boot.effectSizeFDRDown
  }
  if(NumTopGenesMetric == "FoldChangeFDR"){
    all.up = upgenes$boot.foldChangeFDRUp
    all.down = downgenes$boot.foldChangeFDRDown
  }
  names(all.up) = upgenes$rn
  names(all.down) = downgenes$rn
  all = c(all.up,all.down)
  if(NumTopGenesMetric %in% c("SAMscore","EffectSize","FoldChange")){
    top.genes = names(sort(all,decreasing=T))[1:NumTopGenesThresh]
  }else{ #if NumTopGenesMetric is some FDR measurement
    top.genes = names(sort(all,decreasing=F))[1:NumTopGenesThresh]
  }
  upgenes = upgenes[rn %in% top.genes]
  downgenes = downgenes[rn %in% top.genes]
  return(list(upgenes = upgenes,downgenes = downgenes))
}



###-###-###-###-###-###-###-###-###-#
###   .getNumTopGenesPairwise()   ###
###-###-###-###-###-###-###-###-###-#
#gets the top N genes (across both up and downgenes) based on the metric the user chooses
#this is the version for use with Pairwise MANATEE

.getNumTopGenesPairwise <- function(upgenes, downgenes, NumTopGenesThresh, NumTopGenesMetric){
  if(NumTopGenesMetric == "SAMscore"){
    all.up = abs(upgenes$pair.SAM.score)
    all.down = abs(downgenes$pair.SAM.score)
  }
  if(NumTopGenesMetric == "EffectSize"){
    all.up = abs(upgenes$pair.effectSize)
    all.down = abs(downgenes$pair.effectSize)
  }
  if(NumTopGenesMetric == "FoldChange"){
    all.up = abs(upgenes$pair.foldChange)
    all.down = abs(downgenes$pair.foldChange)
  }
  if(NumTopGenesMetric == "ttestFDR"){
    all.up = upgenes$pair.ttestFDRUp
    all.down = downgenes$pair.ttestFDRDown
  }
  names(all.up) = upgenes$rn
  names(all.down) = downgenes$rn
  all = c(all.up,all.down)
  if(NumTopGenesMetric %in% c("SAMscore","EffectSize","FoldChange")){
    top.genes = names(sort(all,decreasing=T))[1:NumTopGenesThresh]
  }else{ #if NumTopGenesMetric=="ttestFDR"
    top.genes = names(sort(all,decreasing=F))[1:NumTopGenesThresh]
  }
  upgenes = upgenes[rn %in% top.genes]
  downgenes = downgenes[rn %in% top.genes]
  return(list(upgenes = upgenes,downgenes = downgenes))
}


