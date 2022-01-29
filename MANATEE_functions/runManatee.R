#' Run the MANATEE framework
#' @description
#' This function runs the specified version of MANATEE on a \code{pooledDataObject}
#' @param pooledDataObject A MANATEE \code{pooledDataObject}. For Bootstrapped MANATEE, sample weights are stored in \code{pooledDataObject$sample.weights} if weighted resampling is desired.
#' For Pairwise MANATEE, class weights are stored in \code{pooledDataObject$class.weights} if automated weighting is not desired
#' @param manatee.type The type of MANATEE analysis that is desired. The current options are Basic, Bootstrapped, and Pairwise (Default: "Basic")
#' @param runLeaveOneOutAnalysis If TRUE, then a leave-one-out analysis will also be performed on any datasets that account for >10% of the total datapoints
#' (or whatever LOOthresh is set to) and that do not account for all of the cases or controls in the pooled data.
#' @param LOOthresh \emph{If runLeaveOneOutAnalysis is TRUE:} Any datasets that account for greater than the \code{LOOthresh} of total datapoints will be included in the leave-one-out analysis
#' @param seed The random seed to use (Default: 1337)
#' @param SAM.perms \emph{Basic MANATEE:} The number of SAM permutations to perform - minimum is 50 (Default: 100)
#' @param compute.localfdr \emph{Basic MANATEE:} If TRUE, local FDRs will be computed for SAM - note that this can take some time (Default: TRUE)
#' @param boot.reps \emph{Bootstrapped MANATEE:} The number of bootstrap repetitions to perform (Default: 10)
#' @param weighted.resampling \emph{Bootstrapped MANATEE:} If TRUE, balanced importance resampling is performed using the weights in \code{pooledDataObject$sample.weights}.
#' If FALSE, then balanced resampling will be performed. If not specified, then weighting will automatically be performed if \code{pooledDataObject$sample.weights} is populated,
#' but weighting will not be performed if \code{pooledDataObject$sample.weights} is NULL
#' @param assume.normal \emph{Bootstrapped MANATEE:} If TRUE, p-values will be calculated assuming a normal distribution for each measured statistic.
#' If FALSE, then p-values will be calculated from the empirical distribution, but note that using the empirical distribution requires far more bootstrap reps to be accurate (Default: TRUE)
#' @param run.perm.test \emph{Bootstrapped MANATEE:} If TRUE, then a permutation test will be used to calculate t-test p-values. Note that the permutation test is likely
#' superior to the bootstrapped p-values for low sample sizes \strong{THIS IS NOT CURRENTLY IMPLEMENTED} (Default: FALSE)
#' @param numCores \emph{Bootstrapped MANATEE:} Number of CPUs to use if parallel computing is desired (Default: 1)
#' @param common.groups \emph{Pairwise MANATEE:} The group(s) that will be common to each comparison.
#' If not specified, common.groups will be the vector of unique names for all samples labeled 1 in pooledDataObject$class
#' @param comparison.groups \emph{Pairwise MANATEE:} The group(s) that will be changed for each comparison.
#' If not specified, comparison.groups will be the vector of unique names for all samples labeled 0 in pooledDataObject$class
#' @param group.colname \emph{Pairwise MANATEE:} this designates the column name of the group vector within the pooledDataObject$pheno
#' dataframe (Default: "group")
#' @param weighting.type \emph{Pairwise MANATEE:} If \code{pooledDataObject$class.weights} is NULL (i.e. the user has not defined their own class weights),
#' then another weighting will be implemented for pooling the statistics across all pairwise comparisons.
#' The current options are "inverse variance" (i.e. inverse of the effect size sample variance), "root sample size", or "none". (Default: "none")
#' @details
#' \strong{GENERIC OVERALL MANATEE DESCRIPTION - add this after writing the paper}
#'
#' \strong{Basic MANATEE}: \cr
#' Runs SAM on the pooled data and extracts the SAM score, q-value, and the local FDR (if specified). The Hedges' g effect size and the fold change
#' are then calculated for each gene. Finally, t-tests are performed for each gene (a one-sided t-test in either direction), and the p-values are FDR corrected.
#' The genes are split into upgenes and downgenes based on the effect size.
#'
#' \strong{Bootstrapped MANATEE}: \cr
#' The data is bootstrapped (with weighted resampling if sample importance weights are provided) to estimate the mean and standard error of the
#' Hedges' g effect size, the fold change, and the SAM score for each gene. The p-value for each of these statistics is then estimated and FDR corrected.
#' Finally, t-tests are performed for each gene (a one-sided t-test in either direction), and the p-values are FDR corrected.
#' The genes are split into upgenes and downgenes based on the effect size.
#'
#' \strong{Pairwise MANATEE}: \cr
#' The data is split into multiple pairwise comparisons (between a common control group and multiple comparison groups). For each comparison,
#' the Hedges' g effect size, the fold change, the SAM score, and the t-test p-values (from a one-sided t-test in either direction) are computed.
#' These statistics are then pooled in a weighted fashion, using either user-defined weights or automatically computed weights
#' (see .runManateePairwise description for specifics on the pooling and weighting). The p-values are FDR corrected,
#' and then the genes are split into upgenes and downgenes based on the effect size.
#'
#' \itemize{\emph{Notes:}
#'   \item throughout this description, "upgenes" refers to genes that were elevated in samples labeled with "1" in \code{$class},
#'         and "downgenes" refers to genes that were elevated in samples labeled with "0" in \code{$class}
#'   \item only for use with logged2 data
#'   \item ignore the "value out of range in 'gammafn'" error for Basic MANATEE
#' }
#' @return Returns a modified version of the \code{pooledDataObject} with the MANATEE results stored in \code{pooledDataObject$manateeResults}
#' and the leave-one-out analysis results stored in \code{pooledDataObject$leaveOneOutAnalysis}
#'   \item{upgenes}{dataframe of the genes that were up in the analysis, with genes in rows and the relevant statistics in columns}
#'   \item{downgenes}{dataframe of the genes that were down in the analysis, with genes in rows and the relevant statistics in columns}
#'   \item{pairwiseResults}{For Paired MANATEE, this is a list of the statistics that were calculated for each pairwise comparison}
#'   \item{type}{string describing the type of MANATEE that was run}
#'   \item{call}{function call}
#' @seealso	\code{\link{filterManatee}}
#' @examples
#'
#' @export
#' @import samr boot pbapply coin parallel
#' @importFrom pbmcapply pbmclapply
#' @author Aditya Rao
runManatee <- function(pooledDataObject, manatee.type="Basic", runLeaveOneOutAnalysis=TRUE, LOOthresh = 0.05, seed=1337, SAM.perms=100, compute.localfdr=TRUE, boot.reps=10, weighted.resampling=NULL,
                       assume.normal=TRUE, run.perm.test=FALSE, numCores=1, common.groups=NULL, comparison.groups=NULL, group.colname="group", weighting.type="none"){
  call = sys.call(which = 0)
  if(!checkManateeObject(pooledDataObject,"pooledDataObject")){stop("Invalid pooledDataObject provided")}
  if(!"class" %in% names(pooledDataObject)){
    stop("Must have a valid class vector stored in pooledDataObject$class")
  }
  if(sum(pooledDataObject$class==0) < 2 || sum(pooledDataObject$class==1) < 2){
    stop("Must have 2 or more cases and 2 or more controls for appropriate calculation of Hedges' g")
  }
  if(!(manatee.type %in% c("Basic","Bootstrapped","Pairwise"))){
    stop("manatee.type must be \"Basic\", \"Bootstrapped\", or \"Pairwise\" - other options may be added later")
  }
  if(runLeaveOneOutAnalysis){library(pbmcapply)} #TEMPORARY
  
  if(manatee.type == "Basic"){
    pooledDataObject = .runManateeBasic(pooledDataObject, seed, SAM.perms, compute.localfdr)
    if(runLeaveOneOutAnalysis){
      pooledDataObject$leaveOneOutAnalysis = .runManateeBasicLOO(pooledDataObject, LOOthresh, seed, SAM.perms, compute.localfdr)
    }
  }else if(manatee.type == "Bootstrapped"){
    pooledDataObject = .runManateeBoot(pooledDataObject, seed, boot.reps, weighted.resampling, assume.normal, run.perm.test, numCores)
    if(runLeaveOneOutAnalysis){
      pooledDataObject$leaveOneOutAnalysis = .runManateeBootLOO(pooledDataObject, LOOthresh, seed, boot.reps, weighted.resampling, assume.normal, run.perm.test, numCores)
    }
  }else if(manatee.type == "Pairwise"){
    pooledDataObject = .runManateePairwise(pooledDataObject, seed, common.groups,comparison.groups, group.colname, weighting.type)
    if(runLeaveOneOutAnalysis){
      pooledDataObject$leaveOneOutAnalysis = .runManateePairwiseLOO(pooledDataObject, LOOthresh, seed, common.groups,comparison.groups, group.colname, weighting.type)
    }
  }
  
  #TO-DO (after paper)
  #add support for .runManateePairwise for multiple common groups
  #try to estimate smoothed version of empirical probability distribution for .runManateeBoot
  #add permutation test to .runManateeBoot
  #maybe add .runManateeEmpirical, if I decide it's worthwhile
  
  #check if either up or down genes are empty
  if(nrow(pooledDataObject$manateeResults$upgenes)==0){
    warning("There are 0 genes in the upgenes table")
  }
  if(nrow(pooledDataObject$manateeResults$downgenes)==0){
    warning("There are 0 genes in the downgenes table")
  }
  
  pooledDataObject$manateeResults$call = call
  return(pooledDataObject)
}



###-###-###-###-###-###-###-##
###   .runManateeBasic()   ###
###-###-###-###-###-###-###-##

#just run a SAM and then calculate vanilla ES and foldchange and then do a t-test and correct to fdr

.runManateeBasic <- function(pooledDataObject, seed=1337, SAM.perms=100, compute.localfdr=TRUE){
  require(samr)
  cat("Running Basic MANATEE\n")
  
  #SAM
  data_Sam = list(x = na.omit(pooledDataObject$genes), y = as.numeric(pooledDataObject$class)+1, genenames = rownames(pooledDataObject$genes), logged2 = TRUE)
  Sam_obj = SAM(data.matrix(na.omit(pooledDataObject$genes)), as.numeric(pooledDataObject$class)+1, resp.type="Two class unpaired", genenames = rownames(pooledDataObject$genes),
                fdr.output = 1, nperms = SAM.perms, logged2=TRUE, random.seed = seed)
  Sam_Siggenes = samr.compute.siggenes.table(Sam_obj$samr.obj, del = Sam_obj$del, data_Sam, Sam_obj$delta.table, min.foldchange = 0,compute.localfdr = compute.localfdr)
  upgenes.SAM = Sam_Siggenes$genes.up[,colnames(Sam_Siggenes$genes.up) %in% c("Score(d)","q-value(%)","localfdr(%)")]
  downgenes.SAM = Sam_Siggenes$genes.lo[,colnames(Sam_Siggenes$genes.lo) %in% c("Score(d)","q-value(%)","localfdr(%)")]
  up.empty = nrow(upgenes.SAM)==0
  down.empty = nrow(downgenes.SAM)==0
  if(up.empty && down.empty){
    stop("SAM returned 0 upgenes and 0 downgenes")
  }
  upgenes.SAM = data.frame(apply(upgenes.SAM,2,as.numeric))
  downgenes.SAM = data.frame(apply(downgenes.SAM,2,as.numeric))
  upgenes.SAM$gene = Sam_Siggenes$genes.up[,2]
  downgenes.SAM$gene = Sam_Siggenes$genes.lo[,2]
  if(compute.localfdr){
    if(!up.empty){colnames(upgenes.SAM) = c("SAM.score","SAM.qValue","SAM.localFDR","gene")}
    if(!down.empty){colnames(downgenes.SAM) = c("SAM.score","SAM.qValue","SAM.localFDR","gene")}
    upgenes.SAM$SAM.localFDR = upgenes.SAM$SAM.localFDR / 100 #adjusting since SAM outputs these as percents
    downgenes.SAM$SAM.localFDR = downgenes.SAM$SAM.localFDR / 100 #adjusting since SAM outputs these as percents
  }else{
    if(!up.empty){colnames(upgenes.SAM) = c("SAM.score","SAM.qValue","gene")}
    if(!down.empty){colnames(downgenes.SAM) = c("SAM.score","SAM.qValue","gene")}
  }
  upgenes.SAM$SAM.qValue = upgenes.SAM$SAM.qValue / 100 #adjusting since SAM outputs these as percents
  downgenes.SAM$SAM.qValue = downgenes.SAM$SAM.qValue / 100 #adjusting since SAM outputs these as percents
  rownames(upgenes.SAM) = upgenes.SAM$gene
  rownames(downgenes.SAM) = downgenes.SAM$gene
  SAMgenes = rbind(upgenes.SAM,downgenes.SAM)
  
  #effect size and foldchange
  y = as.numeric(pooledDataObject$class)
  es.all = apply(pooledDataObject$genes,1,function(x){.getES(x,y)[9:10]})
  es = es.all[1,]
  es.se = es.all[2,]
  fc = apply(pooledDataObject$genes,1,function(x){.getFC(x,y)})
  
  #t-test p-value and FDR
  pval.up = apply(pooledDataObject$genes,1,function(x){
    if(length(na.omit(x[y==1])) < 2 || length(na.omit(x[y==0])) < 2){
      NA
    }else{
      t.test(x[y==1],x[y==0],alternative="greater")$p.value
    }
  })
  pval.down = apply(pooledDataObject$genes,1,function(x){
    if(length(na.omit(x[y==1])) < 2 || length(na.omit(x[y==0])) < 2){
      NA
    }else{
      t.test(x[y==1],x[y==0],alternative="less")$p.value
    }
  })
  fdr.up = p.adjust(pval.up,method="fdr")
  fdr.down = p.adjust(pval.down,method="fdr")
  
  #numStudies info
  numStudies.all = getNumStudies(pooledDataObject)
  
  gene.stats = data.frame(gene=rownames(pooledDataObject$genes), effectSize=es, effectSizeStandardError = es.se, foldChange=fc, ttestPvalUp=pval.up, ttestFDRUp=fdr.up, 
                          ttestPvalDown=pval.down, ttestFDRDown=fdr.down, numStudies=numStudies.all$numStudies, numStudies1=numStudies.all$numStudies1,
                          numStudies0=numStudies.all$numStudies0, propSamples=numStudies.all$propSamples, stringsAsFactors = FALSE)
  
  #adjust the genes that have NA values (i.e. the genes that have fewer than 1 case or control across all datasets)
  #basically make them all 0 (and make all p-values/FDRs 1)
  gene.stats = .NAreplaceVals(gene.stats,zero.string = "effectSize|foldChange", one.string = "FDR|Pval")
  
  upgenes.init = gene.stats[gene.stats$effectSize>=0,]
  downgenes.init = gene.stats[gene.stats$effectSize<0,]
  
  allgenes = c(upgenes.init$gene,downgenes.init$gene)
  SAMgenes = SAMgenes[allgenes,]
  SAMgenes$gene = allgenes
  
  #adjust all NAs for SAMgenes - make NA scores 0, and NA q-values/FDRs 1
  SAMgenes = .NAreplaceVals(SAMgenes,zero.string = "score", one.string = "qValue|localFDR")
  
  #merge results
  upgenes = merge(SAMgenes,upgenes.init,by = "gene")
  downgenes = merge(SAMgenes,downgenes.init,by = "gene")
  rownames(upgenes) = upgenes$gene
  upgenes$gene = NULL
  rownames(downgenes) = downgenes$gene
  downgenes$gene = NULL
  
  pooledDataObject$manateeResults = list(upgenes = upgenes, downgenes = downgenes, type = "Basic MANATEE")
  return(pooledDataObject)
}



###-###-###-###-###-###-###-#
###   .runManateeBoot()   ###
###-###-###-###-###-###-###-#

.runManateeBoot <- function(pooledDataObject, seed=1337, boot.reps=10, weighted.resampling=NULL, assume.normal=TRUE, run.perm.test=FALSE, numCores=1){
  cat("Running Bootstrapped MANATEE\n")
  progress.bar=TRUE #recommended, but can add some overhead (although the overhead is super small - around 2 extra seconds for 1000 parallel iterations)
  if(length(pooledDataObject$class) <= 170){ #factorial(n) for n>170 returns Inf
    max.perms = factorial(length(pooledDataObject$class))
    if(max.perms < boot.reps){
      warning("Number of requested bootstrap reps is greater than the maximum possible permutations of the data and will be replaced by the max permutations")
      boot.reps = max.perms
    }
  }
  if(is.null(pooledDataObject$sample.weights) && !is.null(weighted.resampling) && weighted.resampling){
    stop("If weighted.resampling is TRUE, then sample weights must be provided in pooledDataObject$sample.weights")
  }
  
  if(is.null(pooledDataObject$sample.weights) && is.null(weighted.resampling)){weighted.resampling=FALSE}
  if(!is.null(pooledDataObject$sample.weights) && is.null(weighted.resampling)){weighted.resampling=TRUE}
  
  require(boot)
  require(samr)
  require(pbapply)
  if(run.perm.test){require(coin)}
  if(numCores>1){
    require(parallel)
    require(pbmcapply)
  }
  
  #bootstrapping
  bs <- function(data, indices){
    d = data[indices,]
    xstar = t(d[,2:ncol(d)])
    ystar = as.numeric(d[,1,drop=T])
    #effect size
    n1 <- sum(ystar==0); n2 <- sum(ystar==1)
    if(n1 < 2 || n2 < 2){return(rep(NA,nrow(xstar)*3))}
    m1 <- rowMeans(xstar[, ystar == 1, drop = F],na.rm = TRUE)
    m2 <- rowMeans(xstar[, ystar == 0, drop = F],na.rm = TRUE)
    diff <- m1 - m2
    sd1 <- apply(xstar[, ystar == 1, drop = F],1,sd,na.rm=TRUE) #try to make this faster
    sd2 <- apply(xstar[, ystar == 0, drop = F],1,sd,na.rm=TRUE) #try to make this faster
    sp   <- sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2 )/(n1 + n2 - 2))
    J   <- 1 - 3/(4*(n1 + n2) - 9) #bias correction factor
    g    <- J * diff/sp
    #fold change (calculated with the assumption of log2 normalized data and adjusted to linear scale)
    fc = 2^(m1 - m2)
    fc.pos = which(fc>=1)
    fc.neg = which(fc<1)
    fc[fc.pos] = fc[fc.pos]-1
    fc[fc.neg] = -1*(1/fc[fc.neg]-1)
    #SAM
    init.fit = ttest.func(xstar, ystar+1, sd = NULL)
    if(nrow(xstar) < 500){
      s0 = quantile(init.fit$sd, 0.05, na.rm = TRUE)
    }else{
      s0 = SAM.est.s0(init.fit$tt, init.fit$sd)$s0.hat
    }
    SAM.sd = sqrt(((n2 - 1) * samr:::varr(xstar[, ystar == 1], meanx = m1) +  (n1 - 1) * samr:::varr(xstar[, ystar == 0], meanx = m2)) * (1/n1 + 1/n2)/(n1 + n2 - 2))
    sam = diff/(SAM.sd+s0)
    return(c(g,fc,sam))
  }
  
  boot.data = cbind(pooledDataObject$class,t(pooledDataObject$genes))
  
  n = ncol(boot.data)-1
  curr.pboptions = pboptions()
  pboptions(style = 1,type = "timer")
  set.seed(seed)
  if(weighted.resampling){
    if(progress.bar){
      cat("Bootstrapping progress:\n")
      boot.out = .boot.pb(data = boot.data, statistic = bs, R=boot.reps, sim = "balanced", weights = pooledDataObject$sample.weights, parallel = "multicore", ncpus=numCores)
    }else{
      boot.out = boot(data = boot.data, statistic = bs, R=boot.reps, sim = "balanced", weights = pooledDataObject$sample.weights, parallel = "multicore", ncpus=numCores)
    }
  }else{
    if(progress.bar){
      cat("Bootstrapping progress:\n")
      boot.out = .boot.pb(data = boot.data, statistic = bs, R=boot.reps, sim = "balanced", parallel = "multicore", ncpus=numCores)
    }else{
      boot.out = boot(data = boot.data, statistic = bs, R=boot.reps, sim = "balanced", parallel = "multicore", ncpus=numCores)
    }
  }
  pboptions(curr.pboptions)
  
  means = apply(boot.out$t,2,mean,na.rm=TRUE)
  es.mean = means[1:n]
  fc.mean = means[(n+1):(2*n)]
  sam.mean = means[(2*n+1):(3*n)]
  se = apply(boot.out$t,2,sd) #sd of the bootstrap is equivalent to se
  es.se = se[1:n]
  fc.se = se[(n+1):(2*n)]
  sam.se = se[(2*n+1):(3*n)]
  es.orig = boot.out$t0[1:n]
  fc.orig = boot.out$t0[(n+1):(2*n)]
  sam.orig = boot.out$t0[(2*n+1):(3*n)]
  boot.results = data.frame(boot.effectSize = es.mean, boot.effectSizeStandardError = es.se,
                            boot.SAM.score = sam.mean, boot.SAM.scoreStandardError = sam.se,
                            boot.foldChange = fc.mean, boot.foldChangeStandardError = fc.se,
                            original.effectSize = es.orig, original.SAM.score = sam.orig,
                            original.foldChange = fc.orig)
  
  #calculate p-values (for two-sided test) and FDRs for each statistic
  if(assume.normal){
    stat.pvals = data.frame(t(apply(data.matrix(boot.results),1,function(gene){
      es.pval.up = pnorm(0,mean=gene[1],sd=gene[2])
      es.pval.down = pnorm(0,mean=gene[1],sd=gene[2],lower.tail = FALSE)
      sam.pval.up = pnorm(0,mean=gene[3],sd=gene[4])
      sam.pval.down = pnorm(0,mean=gene[3],sd=gene[4],lower.tail = FALSE)
      fc.pval.up = pnorm(0,mean=gene[5],sd=gene[6])
      fc.pval.down = pnorm(0,mean=gene[5],sd=gene[6],lower.tail = FALSE)
      return(c(es.pval.up,es.pval.down,sam.pval.up,sam.pval.down,fc.pval.up,fc.pval.down))
    })))
  }else{ #getting weirdness here where the p-vals are the same for all three stats - have only tested up to boot.reps=100 tho
    p.bool.up = t(boot.out$t)>=0
    #pval.up = apply(p.bool.up,1,sum)/boot.reps #might be biased
    pval.up = (1+apply(p.bool.up,1,sum))/(boot.reps+1) #potential bias removal (Bootstrap Methods and their Application, p. 141)
    pval.down = 1-pval.up
    stat.pvals = data.frame(pval.up[1:n],pval.down[1:n],pval.up[(2*n+1):(3*n)],
                            pval.down[(2*n+1):(3*n)],pval.up[(n+1):(2*n)],pval.down[(n+1):(2*n)])
  }
  
  #get numStudies info
  numStudies.all = getNumStudies(pooledDataObject)
  boot.results$numStudies = numStudies.all$numStudies
  boot.results$numStudies1 = numStudies.all$numStudies1
  boot.results$numStudies0 = numStudies.all$numStudies0
  boot.results$propSamples = numStudies.all$propSamples
  
  boot.results$boot.effectSizePvalUp = stat.pvals[,1]
  boot.results$boot.effectSizeFDRUp = p.adjust(boot.results$boot.effectSizePvalUp,method="fdr")
  boot.results$boot.effectSizePvalDown = stat.pvals[,2]
  boot.results$boot.effectSizeFDRDown = p.adjust(boot.results$boot.effectSizePvalDown,method="fdr")
  boot.results$boot.SAM.scorePvalUp = stat.pvals[,3]
  boot.results$boot.SAM.scoreFDRUp = p.adjust(boot.results$boot.SAM.scorePvalUp,method="fdr")
  boot.results$boot.SAM.scorePvalDown = stat.pvals[,4]
  boot.results$boot.SAM.scoreFDRDown = p.adjust(boot.results$boot.SAM.scorePvalDown,method="fdr")
  boot.results$boot.foldChangePvalUp = stat.pvals[,5]
  boot.results$boot.foldChangeFDRUp = p.adjust(boot.results$boot.foldChangePvalUp,method="fdr")
  boot.results$boot.foldChangePvalDown = stat.pvals[,6]
  boot.results$boot.foldChangeFDRDown = p.adjust(boot.results$boot.foldChangePvalDown,method="fdr")
  boot.results = boot.results[,c(1:2,14:17,3:4,18:21,5:6,22:25,10:13,7:9)]
  
  #t-test p-value and FDR - MAYBE REMOVE THIS IF I ADD THE PERMUTATION TEST
  y = as.numeric(pooledDataObject$class)
  boot.results$original.ttestPvalUp = apply(pooledDataObject$genes,1,function(x){
    if(length(na.omit(x[y==1])) < 2 || length(na.omit(x[y==0])) < 2){
      NA
    }else{
      t.test(x[y==1],x[y==0],alternative="greater")$p.value
    }
  })
  boot.results$original.ttestPvalDown = apply(pooledDataObject$genes,1,function(x){
    if(length(na.omit(x[y==1])) < 2 || length(na.omit(x[y==0])) < 2){
      NA
    }else{
      t.test(x[y==1],x[y==0],alternative="less")$p.value
    }
  })
  boot.results$original.ttestFDRUp = p.adjust(pval.up,method="fdr")
  boot.results$original.ttestFDRDown = p.adjust(pval.down,method="fdr")
  
  #permutation test
  if(run.perm.test){
    cat("Beginning Permutation Testing:\n")
    warning("Permutation Test is not currently supported and will not be run, but will be added later")
    #use the package coin - go to the section on LocationTests
  }
  
  #adjust the genes that have NA values (i.e. the genes that have fewer than 1 case or control across all datasets)
  boot.results = .NAreplaceVals(boot.results,zero.string = "effectSize$|foldChange$|score$|Error", one.string = "FDR|Pval")
  
  #split into upgenes and downgenes
  up.index = which(boot.results$boot.effectSize>=0)
  down.index = which(boot.results$boot.effectSize<0)
  upgenes = boot.results[up.index,]
  downgenes = boot.results[down.index,]
  
  pooledDataObject$manateeResults = list(upgenes = upgenes, downgenes = downgenes, type = "Bootstrapped MANATEE")
  return(pooledDataObject)
}


###-###-###-###-###-###-###-###-#
###   .runManateePairwise()   ###
###-###-###-###-###-###-###-###-#

#FOR NOW, WRITE PAIRWISE MANATEE TO ONLY WORK WITH ONE COMMON GROUP, BUT ADD SUPPORT FOR IT LATER

#definitely need to do testing with user defined weights because i'm not testing it rn

#WEIGHTS
#if user-defined weights are NULL, either use the inverse variance of the effect size, or the square root of the sample size

#DESCRIPTION
#use the statistics from borenstein to get the pooled mean and variance (weighted version for mean is easy, but need to figure out how to do weighted version for variance)
#use one sided t-tests to get the p-values (up and down) for each gene
#use weighted Z-test to combine p-values from each comparison (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3135688/)
#use the modification for dependent p-values using the equation to calculate correlation given a common "control" group (https://www.tandfonline.com/doi/abs/10.1080/01621459.1955.10501294)
#use p.adjust to get q-values

.runManateePairwise <- function(pooledDataObject, seed=1337, common.groups=NULL,comparison.groups=NULL,group.colname="group",weighting.type="none"){
  cat("Running Pairwise MANATEE\n")
  if(!(weighting.type %in% c("root sample size","inverse variance","none"))){
    stop("weighting.type must be either \"root sample size\", \"inverse variance\", or \"none\" - other options may be added later")
  }
  
  require(samr)
  
  #separate into groups
  if(is.null(common.groups)){
    common.groups = as.character(unique(pooledDataObject$pheno[,group.colname][pooledDataObject$class==1]))
    cat(paste0("Found ",length(common.groups)," common group",ifelse(length(common.groups)>1,"s",""),"\n"))
  }
  if(length(common.groups)>1){stop("Running Pairwise MANATEE with multiple common groups is not currently supported, but will be added later")} #update
  if(is.null(comparison.groups)){
    comparison.groups = as.character(unique(pooledDataObject$pheno[,group.colname][pooledDataObject$class==0]))
    names(comparison.groups) = comparison.groups
    cat(paste0("Found ",length(comparison.groups)," comparison group",ifelse(length(comparison.groups)>1,"s",""),"\n"))
  }
  if(length(comparison.groups)<=1){stop("Pairwise MANATEE requires at least 2 comparison groups")}
  
  my.weights = NULL
  if(!is.null(pooledDataObject$class.weights)){
    if(any(!(comparison.groups %in% names(pooledDataObject$class.weights)))){
      stop("The names of pooledDataObject$class.weights don't include all of the comparison group names")
    }
    my.weights = pooledDataObject$class.weights[comparison.groups]
    my.weights = my.weights/sum(my.weights,na.rm=TRUE) #normalize to one
  }
  
  #get statistics for each comparison
  pairwiseResults = lapply(comparison.groups, function(comparison){
    cindex = which(pooledDataObject$pheno[,group.colname] %in% c(common.groups,comparison)) #update
    xstar = pooledDataObject$genes[,cindex]
    ystar = as.numeric(pooledDataObject$class[cindex])
    stats = apply(xstar,1,function(gene){
      es.all = .getES(gene,ystar)
      es1 = es.all[9]
      es.se1 = es.all[10]
      fc1 = .getFC(gene,ystar)
      if(length(na.omit(gene[ystar==1])) < 2 || length(na.omit(gene[ystar==0])) < 2){
        pval.up1 = NA
        pval.down1 = NA
      }else{
        pval.up1 = t.test(gene[ystar==1],gene[ystar==0],alternative="greater")$p.value
        pval.down1 = t.test(gene[ystar==1],gene[ystar==0],alternative="less")$p.value
      }
      n0.1 = es.all[1]
      n1.1 = es.all[4]
      return(c(es1,es.se1,fc1,pval.up1,pval.down1,n0.1,n1.1))
    })
    es = stats[1,]
    es.se = stats[2,]
    fc = stats[3,]
    pval.up = stats[4,]
    pval.down = stats[5,]
    n0 = stats[6,]
    n1 = stats[7,]
    #SAM
    xstar = as.matrix(xstar) # fixing this error: requires numeric/complex matrix/vector arguments
    init.fit = ttest.func(xstar, ystar+1, sd = NULL)
    if(nrow(xstar) < 500){
      s0 = quantile(init.fit$sd, 0.05, na.rm=TRUE)
    }else{
      s0 = SAM.est.s0(init.fit$tt, init.fit$sd)$s0.hat
    }
    sam = ttest.func(xstar, ystar+1, s0 = s0, sd = NULL)$tt
    cat(paste0("Done with comparison group: ",comparison,"\n"))
    result = data.frame(es=es, es.se=es.se, sam=sam, fc=fc, pval.up=pval.up, pval.down=pval.down, n0=n0, n1=n1)
    rownames(result) = names(es)
    return(result)
  })
  
  #right now i'm assuming that genes are equivalent between comparisons, but this may not be true in later updates
  all.genes = c()
  for(i in 1:length(pairwiseResults)){
    all.genes = union(all.genes,rownames(pairwiseResults[[i]]))
  }
  
  #combine all the statistics together using the given weighting
  pooledResults = data.frame(t(sapply(all.genes, function(gene){
    gene.stats = data.frame(do.call(rbind,lapply(pairwiseResults,function(x){x[gene,]})))
    if(!is.null(my.weights)){
      weights = my.weights
      gene.stats = gene.stats[names(my.weights),] #reorder by names of my.weights
    }else{
      if(weighting.type == "inverse variance"){
        weights = 1/(gene.stats$es.se^2) #inverse variance weighting
      }else if(weighting.type == "root sample size"){
        weights = sqrt(gene.stats$n0+gene.stats$n1) #square root of sample size weighting
      }else{
        weights = rep(1,length(gene.stats$es.se)) #no weighting
      }
      weights = weights/sum(weights,na.rm=TRUE) #normalize to one
    }
    es.pool = sum(gene.stats$es*weights,na.rm=TRUE)
    sam.pool = sum(gene.stats$sam*weights,na.rm=TRUE)
    fc.pool = sum(gene.stats$fc*weights,na.rm=TRUE)
    
    es.var.pool1 = sum((gene.stats$es.se^2)*weights,na.rm=TRUE)
    es.var.pool2 = 0
    sum.cor.Z = 0 #this is used for the weighted Z-test (sum of weights*correlation)
    combos = combn(1:nrow(gene.stats),2)
    for(i in 1:ncol(combos)){
      c = combos[,1]
      n0 = gene.stats$n0[1]
      n1 = gene.stats$n1[c[1]]
      n2 = gene.stats$n1[c[2]]
      corCtrl = sqrt((1/(1+n0/n1))*((1/(1+n0/n2))))
      incr.es = 2*prod(weights[c],na.rm=TRUE)*corCtrl*prod(gene.stats$es.se[c],na.rm=TRUE)
      incr.Z = 2*prod(weights[c],na.rm=TRUE)*corCtrl
      es.var.pool2 = es.var.pool2 + incr.es
      sum.cor.Z = sum.cor.Z + incr.Z
    }
    es.se.pool = sqrt(es.var.pool1 + es.var.pool2)
    
    pval.up.pool = .weightedZtest.ctrl(gene.stats$pval.up,weights,sum.cor.Z)
    pval.down.pool = .weightedZtest.ctrl(gene.stats$pval.down,weights,sum.cor.Z)
    
    all.pool = c(es.pool,es.se.pool,sam.pool,fc.pool,pval.up.pool,pval.down.pool)
    names(all.pool) = c("pair.effectSize","pair.effectSizeStandardError","pair.SAM.score","pair.foldChange","pair.ttestPvalUp","pair.ttestPvalDown")
    return(all.pool)
  })))
  
  #get FDRs
  pooledResults$pair.ttestFDRUp = p.adjust(pooledResults$pair.ttestPvalUp,method="fdr")
  pooledResults$pair.ttestFDRDown = p.adjust(pooledResults$pair.ttestPvalDown,method="fdr")
  
  #get numStudies info
  numStudies.all = getNumStudies(pooledDataObject)
  pooledResults$numStudies = numStudies.all$numStudies
  pooledResults$numStudies1 = numStudies.all$numStudies1
  pooledResults$numStudies0 = numStudies.all$numStudies0
  pooledResults$propSamples = numStudies.all$propSamples
  
  #get original stats
  y = as.numeric(pooledDataObject$class)
  orig.stats = apply(pooledDataObject$genes,1,function(gene){
    es.orig1 = .getES(gene,y)[9]
    fc.orig1 = .getFC(gene,y)
    return(c(es.orig1,fc.orig1))
  })
  init.fit = ttest.func(as.matrix(pooledDataObject$genes), y+1, sd = NULL)
  if(nrow(pooledDataObject$genes) < 500){
    s0 = quantile(init.fit$sd, 0.05, na.rm = TRUE)
  }else{
    s0 = SAM.est.s0(init.fit$tt, init.fit$sd)$s0.hat
  }
  orig.sam = ttest.func(as.matrix(pooledDataObject$genes), y+1, s0 = s0, sd = NULL)$tt
  pooledResults$original.effectSize = orig.stats[1,]
  pooledResults$original.SAM.score = orig.sam
  pooledResults$original.foldChange = orig.stats[2,]
  
  #adjust the genes that have NA values (i.e. the genes that have fewer than 1 case or control across all datasets)
  pooledResults = .NAreplaceVals(pooledResults,zero.string = "effectSize$|foldChange$|score$|Error", one.string = "FDR|Pval")
  
  #reformat - right now, using effect size to split into up and down genes
  pooledResults = pooledResults[,c(1:5,7,6,8:15)]
  up.index = which(pooledResults$pair.effectSize>=0)
  down.index = which(pooledResults$pair.effectSize<0)
  upgenes = pooledResults[up.index,]
  downgenes = pooledResults[down.index,]
  
  pooledDataObject$manateeResults = list(upgenes = upgenes, downgenes = downgenes, pairwiseResults = pairwiseResults, type = "Pairwise MANATEE")
  return(pooledDataObject)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#--------------------------------------------------------LOO Wrappers--------------------------------------------------------#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

###-###-###-###-###-###-###-###-#
###   .runManateeBasicLOO()   ###
###-###-###-###-###-###-###-###-#
.runManateeBasicLOO <- function(pooledDataObject, LOOthresh = 0.05, seed=1337, SAM.perms=100, compute.localfdr=TRUE){
  LOOnames = .getLOOdatasetNames(pooledDataObject, LOOthresh)
  
  if(length(LOOnames)>0){
    #smart core usage = use 75% of resources unless you need less
    maxCores <- round(detectCores()*7.5/10)
    if(length(LOOnames)<maxCores){
      maxCores <- length(LOOnames)
    }
    
    cat("Beginning Leave-one-out Analysis:\n")
    LOOresults = pbmclapply(mc.cores = maxCores, LOOnames, function(gsename){
      pooledDataObject.remove = removeSamples(pooledDataObject, which(pooledDataObject$pheno$name==gsename),expr=FALSE)
      manateeResults.remove = .runManateeBasic(pooledDataObject.remove, seed, SAM.perms, compute.localfdr)
      manateeResultsLOO.remove = manateeResults.remove$manateeResults
      rm(manateeResults.remove)
      return(manateeResultsLOO.remove)
    })
    names(LOOresults) = paste0("removed_",LOOnames)
    return(LOOresults)
  }else{
    warning(paste("No datasets qualified for LOO analysis with the given LOOthresh of ",LOOthresh))
    return("No datasets qualified for LOO analysis")
  }
}



###-###-###-###-###-###-###-##-#
###   .runManateeBootLOO()   ###
###-###-###-###-###-###-###-##-#
.runManateeBootLOO <- function(pooledDataObject, LOOthresh = 0.05, seed=1337, boot.reps=10, weighted.resampling=NULL, assume.normal=TRUE, run.perm.test=FALSE, numCores=1){
  LOOnames = .getLOOdatasetNames(pooledDataObject, LOOthresh)
  
  if(length(LOOnames)>0){
    #smart core usage = use 75% of resources unless you need less
    maxCores <- round(detectCores()*7.5/10)
    if(length(LOOnames)<maxCores){
      maxCores <- length(LOOnames)
    }
    
    cat("Beginning Leave-one-out Analysis:\n")
    LOOresults = pbmclapply(mc.cores = maxCores, LOOnames, function(gsename){
      pooledDataObject.remove = removeSamples(pooledDataObject, which(pooledDataObject$pheno$name==gsename),expr=FALSE)
      manateeResults.remove = .runManateeBoot(pooledDataObject.remove, seed, boot.reps, weighted.resampling, assume.normal, run.perm.test, numCores)
      manateeResultsLOO.remove = manateeResults.remove$manateeResults
      rm(manateeResults.remove)
      return(manateeResultsLOO.remove)
    })
    names(LOOresults) = paste0("removed_",LOOnames)
    return(LOOresults)
  }else{
    warning(paste("No datasets qualified for LOO analysis with the given LOOthresh of ",LOOthresh))
    return("No datasets qualified for LOO analysis")
  }
}



###-###-###-###-###-###-###-###-##-#
###   .runManateePairwiseLOO()   ###
###-###-###-###-###-###-###-###-##-#
.runManateePairwiseLOO <- function(pooledDataObject, LOOthresh = 0.05, seed=1337, common.groups=NULL, comparison.groups=NULL, group.colname="group", weighting.type="none"){
  LOOnames = .getLOOdatasetNames(pooledDataObject, LOOthresh)
  
  if(length(LOOnames)>0){
    #smart core usage = use 75% of resources unless you need less
    maxCores <- round(detectCores()*7.5/10)
    if(length(LOOnames)<maxCores){
      maxCores <- length(LOOnames)
    }
    
    cat("Beginning Leave-one-out Analysis:\n")
    LOOresults = pbmclapply(mc.cores = maxCores, LOOnames, function(gsename){
      pooledDataObject.remove = removeSamples(pooledDataObject, which(pooledDataObject$pheno$name==gsename),expr=FALSE)
      manateeResults.remove = .runManateePairwise(pooledDataObject.remove, seed, common.groups,comparison.groups, group.colname, weighting.type)
      manateeResultsLOO.remove = manateeResults.remove$manateeResults
      rm(manateeResults.remove)
      return(manateeResultsLOO.remove)
    })
    names(LOOresults) = paste0("removed_",LOOnames)
    return(LOOresults)
  }else{
    warning(paste("No datasets qualified for LOO analysis with the given LOOthresh of ",LOOthresh))
    return("No datasets qualified for LOO analysis")
  }
}



###-###-###-###-###-###-###-###-#
###   .getLOOdatasetNames()   ###
###-###-###-###-###-###-###-###-#
#get the names of all the datasets that qualify for the LOO analysis
.getLOOdatasetNames <- function(pooledDataObject, LOOthresh = 0.05){
  #get the names of all of the datasets that account for >10% (or whatever the LOOthresh is) of the datapoints
  #but make sure that removing the dataset doesn't remove all of one class from the analysis
  LOOdatasetNames = vector("character")
  for(name in unique(pooledDataObject$pheno$name)){
    if((length(which(pooledDataObject$pheno$name == name))/length(pooledDataObject$pheno$name)) >= LOOthresh){
      class.in = pooledDataObject$class[which(pooledDataObject$pheno$name != name)]
      if((0 %in% class.in) && (1 %in% class.in)){ #note that right now, this only accounts for having two classes
        LOOdatasetNames = c(LOOdatasetNames,name)
      }
    }
  }
  return(LOOdatasetNames)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#--------------------------------------------------------Other Supporting Functions--------------------------------------------------------#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



###-###-###-###-###-###-###-###-#
###   .weightedZtest.ctrl()   ###
###-###-###-###-###-###-###-###-#

#DESCRIPTION
#Performs the weighted Z-test to combine p-values for dependent tests done on multiple comparisons with a common "control" group

#NOTE: this is very precise for low p-values (up to 1e-323), but is not very precise for p values very close to 1 (only up to 1-1e-16)
#Any p-values that are below or above these precision limits, respectively, are adjusted to be equal to the precision limit

#PARAMETERS
#pvals - one sided p-values that you want to combine
#weights - weights for each p-value
#sum.cor.Z - the sum of the correlations*weights between all the comparisons - can be pre-calculated and provided to the function
#n0 - if sum.cor.Z=NULL, the number of samples in the "control" group
#n1.vec - if sum.cor.Z=NULL, numeric vector of the number of samples in each non-"control" group

#RETURN VALUE
#returns the pooled p-value

#REQUIRED PACKAGES: None

.weightedZtest.ctrl <- function(pvals,weights,sum.cor.Z=NULL,n0=NULL,n1.vec=NULL){
  if(is.null(sum.cor.Z)){
    sum.cor.Z = 0
    combos = combn(1:length(n1.vec),2)
    for(i in 1:ncol(combos)){
      c = combos[,1]
      n1 = n1.vec[c[1]]
      n2 = n1.vec[c[2]]
      corCtrl = sqrt((1/(1+n0/n1))*((1/(1+n0/n2))))
      incr.Z = 2*prod(weights[c],na.rm=TRUE)*corCtrl
      sum.cor.Z = sum.cor.Z + incr.Z
    }
  }
  
  for(i in 1:length(pvals)){
    if(!is.na(pvals[i])){
      if(pvals[i] < 1e-323){
        pvals[i] = 1e-323
      }
      if(pvals[i] > (1-1e-16)){
        pvals[i] = (1-1e-16)
      }
    }
  }
  
  Z = qnorm(1-pvals)
  if(any(is.infinite(Z))){ #can't really help for -Inf, but if p-value is close to 1 I don't really care anyway
    Z = qnorm(pvals, lower.tail = FALSE)
    num = sum(weights*Z,na.rm=TRUE)
    dnm = sqrt(sum(weights^2,na.rm=TRUE)+sum.cor.Z)
    pval.pool = pnorm(num/dnm, lower.tail = FALSE)
  }else{
    num = sum(weights*Z,na.rm=TRUE)
    dnm = sqrt(sum(weights^2,na.rm=TRUE)+sum.cor.Z)
    pval.pool = 1 - pnorm(num/dnm)
  }
  return(pval.pool)
}



###-###-###-###-###-##-#
###   SAM.est.s0()   ###
###-###-###-###-###-##-#
#version of the est.s0 function from samr that works with NAs
SAM.est.s0 <- function(tt, sd, s0.perc = seq(0, 1, by = 0.05)){
  br = unique(quantile(sd, seq(0, 1, len = 101),na.rm=TRUE))
  nbr = length(br)
  a <- cut(sd, br, labels = F)
  a[is.na(a)] <- 1
  cv.sd <- rep(0, length(s0.perc))
  for (j in 1:length(s0.perc)) {
    w <- quantile(sd, s0.perc[j], na.rm=TRUE)
    w[j == 1] <- 0
    tt2 <- tt * sd/(sd + w)
    tt2[tt2 == Inf] = NA
    sds <- rep(0, nbr - 1)
    for (i in 1:(nbr - 1)) {
      sds[i] <- mad(tt2[a == i], na.rm = TRUE)
    }
    cv.sd[j] <- sqrt(var(sds))/mean(sds)
  }
  o = (1:length(s0.perc))[cv.sd == min(cv.sd)]
  s0.hat = quantile(sd[sd != 0], s0.perc[o], na.rm=TRUE)
  return(list(s0.perc = s0.perc, cv.sd = cv.sd, s0.hat = s0.hat))
}



###-###-###-###-###-##
###   .boot.pb()   ###
###-###-###-###-###-##
#jerry-rigged version of boot that adds a progress bar for both sequential and parallel bootstrapping (with multicore only)
.boot.pb <- function(data, statistic, R, sim = "ordinary", stype = c("i", "f", "w"), strata = rep(1, n), L = NULL,
                     m = 0, weights = NULL, ran.gen = function(d, p) d, mle = NULL, simple = FALSE, ...,
                     parallel = c("no", "multicore", "snow"), ncpus = getOption("boot.ncpus", 1L), cl = NULL){
  loadNamespace("boot")
  call <- match.call()
  stype <- match.arg(stype)
  if (missing(parallel)) 
    parallel <- getOption("boot.parallel", "no")
  parallel <- match.arg(parallel)
  have_mc <- have_snow <- FALSE
  if (parallel != "no" && ncpus > 1L) {
    if (parallel == "multicore") 
      have_mc <- .Platform$OS.type != "windows"
    else if (parallel == "snow") 
      have_snow <- TRUE
    if (!have_mc && !have_snow) 
      ncpus <- 1L
    loadNamespace("parallel")
  }
  if (simple && (sim != "ordinary" || stype != "i" || sum(m))) {
    warning("'simple=TRUE' is only valid for 'sim=\"ordinary\", stype=\"i\", n=0', so ignored")
    simple <- FALSE
  }
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
    runif(1)
  seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  n <- NROW(data)
  if ((n == 0) || is.null(n)) 
    stop("no data in call to 'boot'")
  temp.str <- strata
  strata <- tapply(seq_len(n), as.numeric(strata))
  t0 <- if (sim != "parametric") {
    if ((sim == "antithetic") && is.null(L)) 
      L <- empinf(data = data, statistic = statistic, stype = stype, strata = strata, ...)
    if (sim != "ordinary") 
      m <- 0
    else if (any(m < 0)) 
      stop("negative value of 'm' supplied")
    if ((length(m) != 1L) && (length(m) != length(table(strata)))) 
      stop("length of 'm' incompatible with 'strata'")
    if ((sim == "ordinary") || (sim == "balanced")) {
      if (boot:::isMatrix(weights) && (nrow(weights) != length(R))) 
        stop("dimensions of 'R' and 'weights' do not match")
    }
    else weights <- NULL
    if (!is.null(weights)) 
      #weights <- t(apply(matrix(weights, n, length(R), byrow = TRUE), 2L, normalize, strata))
      weights = weights/sum(weights)
    if (!simple) 
      i <- boot:::index.array(n, R, sim, strata, m, L, weights)
    original <- if (stype == "f") 
      rep(1, n)
    else if (stype == "w") {
      ns <- tabulate(strata)[strata]
      1/ns
    }
    else seq_len(n)
    t0 <- if (sum(m) > 0L) 
      statistic(data, original, rep(1, sum(m)), ...)
    else statistic(data, original, ...)
    rm(original)
    t0
  }
  else statistic(data, ...)
  pred.i <- NULL
  fn <- if (sim == "parametric") {
    ran.gen
    data
    mle
    function(r) {
      dd <- ran.gen(data, mle)
      statistic(dd, ...)
    }
  }
  else {
    if (!simple && ncol(i) > n) {
      pred.i <- as.matrix(i[, (n + 1L):ncol(i)])
      i <- i[, seq_len(n)]
    }
    if (stype %in% c("f", "w")) {
      f <- freq.array(i)
      rm(i)
      if (stype == "w") 
        f <- f/ns
      if (sum(m) == 0L) 
        function(r) statistic(data, f[r, ], ...)
      else function(r) statistic(data, f[r, ], pred.i[r, ], ...)
    }
    else if (sum(m) > 0L) 
      function(r) statistic(data, i[r, ], pred.i[r, ], ...)
    else if (simple) 
      function(r) statistic(data, boot:::index.array(n, 1, sim, strata, m, L, weights), ...)
    else function(r) statistic(data, i[r, ], ...)
  }
  RR <- sum(R)
  res <- if (ncpus > 1L && (have_mc || have_snow)) {
    if (have_mc) {
      pbmcapply::pbmclapply(seq_len(RR), fn, mc.cores = ncpus)
    }
    else if (have_snow) {
      list(...)
      if (is.null(cl)) {
        cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
        if (RNGkind()[1L] == "L'Ecuyer-CMRG") 
          parallel::clusterSetRNGStream(cl)
        res <- parallel::parLapply(cl, seq_len(RR), fn)
        parallel::stopCluster(cl)
        res
      }
      else parallel::parLapply(cl, seq_len(RR), fn)
    }
  }
  else pbapply::pblapply(seq_len(RR), fn)
  t.star <- matrix(NA, RR, length(t0)) #the first arg here was blank so I changed it to NA
  for (r in seq_len(RR)) t.star[r, ] <- res[[r]]
  if (is.null(weights)) 
    weights <- 1/tabulate(strata)[strata]
  boot:::boot.return(sim, t0, t.star, temp.str, R, data, statistic, 
                     stype, call, seed, L, m, pred.i, weights, ran.gen, mle)
}
