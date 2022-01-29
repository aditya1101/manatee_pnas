# New MANATEE functions
# Author: Aditya Rao
# Date: 1/25/18
# Contact: adityamr@stanford.edu

# New functions for MANATEE that are aimed at streamlining the framework

# NOTE: This will eventually be combined with the MANATEE_functions file, but I
# first want to get all the functions organized here so that I can eventually
# organize the MANATEE_functions into "base" functions and "MANATEE-specific"
# functions


###-###-###-###-###-###-###-###
###   formatManateeData()   ###
###-###-###-###-###-###-###-###

#DESCRIPTION
#This function formats your data so that it can be used with the MANATEE pipeline

#PARAMETERS
#GSElist - named list of datasets, formatted as either MANATEE Dataset Objects or MetaIntegrator Dataset Objects
#MetaToManatee - make TRUE if your datasets are MetaIntegrator Dataset Objects, which will be converted into MANATEE Dataset Objects (Default: FALSE)
#AddDatasetName - if MetaToManatee is FALSE, and if AddDatasetName is TRUE, then a "name" column in $pheno and an all-uppercase string in $formattedName will be added
#               - to each of your datasets based on the names of GSElist. If either the "name" column or $formattedName already exist, then they will not be changed (Default: FALSE)
#AddControl0Class - make TRUE if your datasets don't already have the $control.0.class column in $pheno in order to generate this column for each dataset (Default: TRUE)
#ControlNames - if AddControl0Class=TRUE, then this specifies the name(s) of your "control" class for creation of the $pheno$control.0.class column.
#               If set as NULL, then each dataset's $class will be copied into $pheno$control.0.class (Default: "Healthy")
#knn.neighbors - What to set "k" equal to for knn imputation (Default: 10)
#NAcutoff - Any genes with a higher proportion of NAs than the NAcutoff will be automatically removed (Default: 0.3)
#ExcelDateFix - make TRUE if you want to correct for Excel changing some gene names to dates (Default: TRUE)

#RETURN VALUE
#A modified version of the input GSElist that is a valid MANATEE DatasetList Object

#REQUIRED PACKAGES: None

formatManateeData <- function(GSElist, MetaToManatee=FALSE, AddDatasetName=TRUE, AddControl0Class=TRUE, ControlNames="Healthy", ImputeMissing = TRUE, knn.neighbors = 10, NAcutoff = 0.3, ExcelDateFix=TRUE){
  if(MetaToManatee){
    GSElist = lapply(GSElist,function(gse){
      return(convertMetaToManatee(gse,ControlNames=ControlNames))
    })
  }else{
    #make sure class vector is present in gse$class and gse$pheno$control.0.class, and if not, add a dummy class vector
    GSElist = lapply(GSElist, function(gse){
      if(!("class" %in% names(gse))){
        gse$class = rep(0,nrow(gse$pheno))
      }
      if(AddControl0Class){
        if(is.null(ControlNames)){
          gse$pheno$control.0.class = gse$class
        }else{
          gse$pheno$control.0.class = makeControl0Class(gse$pheno$group,ControlNames)
        }
      }
      return(gse)
    })
    
    if(AddDatasetName){ #MAKE SURE THIS WORKS DUDE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      GSEnames = names(GSElist)
      GSElist = lapply(GSEnames,function(gsename){
        gse = GSElist[[gsename]]
        if(is.null(gse$pheno$name)){
          gse$pheno$name = gsename
        }
        if(is.null(gse$formattedName)){
          gse$formattedName = toupper(gsename)
        }
        return(gse)
      })
      names(GSElist) = GSEnames
    }
    #add a check for preformatted Manatee DatasetList
  }
  
  GSElist = lapply(GSElist, function(gse){
    if(!(any(gse$pheno$control.0.class == 0) && any(gse$pheno$control.0.class == 1)) || any(!gse$pheno$control.0.class %in% c(0,1))){
      stop(paste0("Error in Dataset ",gse$formattedName,": Dataset$pheno$control.0.class must contain both 0s (to indicate healthy samples) and 1s (to indicate non-healthy samples)",
                  " and cannot contain any other digits. This error is usually caused when none of the ControlNames are present in Dataset$pheno$group"))
    }
    return(gse)
  })

  if(!checkManateeObject(GSElist,"DatasetList")){stop("This is not a valid MANATEE DatasetList Object")}

  #combine duplicate genes and impute NAs with knn
  GSElist = lapply(GSElist,function(gse){
    if(is.null(gse$genes)){
      if(is.null(gse$expr) || is.null(gse$keys)){
        stop("If $genes is not present, then $expr and $keys must be available")
      }
      gse$genes = data.matrix(gse$expr)
      rownames(gse$genes)=gse$keys
    }
    gse$genes = combineDuplicateGenes(gse$genes)
    if (sum(is.na(gse$genes)) > 0) {
      gse$genes = apply(gse$genes,1,function(gene){
        if(sum(is.na(gene))/length(gene) <= NAcutoff){
          return(gene)
        }
      })
      if(class(gse$genes) == "matrix"){
        gse$genes = t(gse$genes)
      }else{
        gse$genes = do.call(rbind,gse$genes)
      }
      if(sum(is.na(gse$genes)) > 0){
        warning("Some datasets are missing genes, if ImputeMissing=TRUE then KNN imputation will be performed")
        if(ImputeMissing){
          require(impute)
          gse$genes = impute.knn(gse$genes, k = knn.neighbors,rng.seed = 1337)$data
        }
        # if (!is.matrix(gse$genes)) {
        #   gse$genes = gse$genes$data
        # }
      }

    }
    return(gse)
  })

  if(ExcelDateFix){
    GSElist = lapply(GSElist,function(gse){
      if("keys" %in% names(gse)){
        corrected_names = excelDateCorrection(rownames(gse$genes),names(gse$keys),print=FALSE)
      }else{
        corrected_names = excelDateCorrection(rownames(gse$genes),print=FALSE)
      }
      if(any(duplicated(corrected_names))){ #if the date or HGNC correction has introduced duplicate gene names
        gse$genes = data.matrix(gse$genes)
        rownames(gse$genes) = corrected_names
        gse$genes = combineDuplicateGenes(gse$genes)
      }else{
        rownames(gse$genes) = corrected_names
      }
      return(gse)
    })
  }
  
  return(GSElist)
}



###-###-###-###-###-###-###-###-##
###   convertMetaToManatee()   ###
###-###-###-###-###-###-###-###-##

#DESCRIPTION
#This function takes a MetaIntegrator Dataset Object and converts it into a MANATEE Dataset Object

#Note that unlike in MetaIntegrator, all datasets must have $class populated

#PARAMETERS
#dataset - the dataset to be converted
#ControlNames - this specifies the name(s) of your "control" class for creation of the $pheno$control.0.class column.
#               If set as NULL, then the dataset's $class will be copied into $pheno$control.0.class (Default: "Healthy")

#RETURN VALUE
#A modified version of the input dataset that is a valid MANATEE Dataset Object

#REQUIRED PACKAGES: MetaIntegrator

#NOTE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#currently, if the $class vector is NULL, checkDataObject returns TRUE, but if
#this is fixed in the MetaIntegrator package then I won't do anything about it
#here

convertMetaToManatee <- function(dataset, ControlNames="Healthy"){
  require(MetaIntegrator)

  if(is.null(dataset$class)){
    if("formattedName" %in% names(dataset)){
      stop(paste0(dataset$formattedName," is missing $class"))
    } else{
      stop("This dataset is missing $class, and does not have a \"formattedName\" field so the problematic dataset's name cannot be identified")
    }
  }

  #TEMPORARY - add names to class
  if(is.null(names(dataset$class)) && length(dataset$class) == nrow(dataset$pheno)){
    names(dataset$class) = rownames(dataset$pheno)
  }

  if(!checkDataObject(dataset, "Dataset")){
    if("formattedName" %in% names(dataset)){
      stop(paste0(dataset$formattedName," is not a valid MetaIntegrator Dataset Object"))
    } else{
      stop("This dataset is not a valid MetaIntegrator Dataset Object, and does not have a \"formattedName\" field so the problematic dataset's name cannot be identified")
    }
  }
  if(!("group" %in% colnames(dataset$pheno))){
    stop(c("\"group\" column is required for all datasets, but is missing from ",dataset$formattedName))
  }

  if(!("genes" %in% names(dataset))){
    dataset$genes = data.matrix(dataset$expr)
    rownames(dataset$genes)=dataset$keys
    dataset$genes <- combineDuplicateGenes(dataset$genes)
  }
  if(!("name" %in% colnames(dataset$pheno))){
    dataset$pheno$name = dataset$formattedName
  }
  if(!("control.0.class" %in% colnames(dataset$pheno))){
    if(is.null(ControlNames)){
      dataset$pheno$control.0.class = as.numeric(dataset$class) #potential to break down here because of MetaIntegrator bug; see above
    }else{
      dataset$pheno$control.0.class = makeControl0Class(dataset$pheno$group,ControlNames)
    }
  }
  if(!(any(dataset$pheno$control.0.class == 0) && any(dataset$pheno$control.0.class == 1)) || any(!dataset$pheno$control.0.class %in% c(0,1))){
    stop(paste0("Error in Dataset ",dataset$formattedName,": Dataset$pheno$control.0.class must contain both 0s (to indicate healthy samples) and 1s (to indicate non-healthy samples)",
                " and cannot contain any other digits. This error is usually caused when none of the ControlNames are present in Dataset$pheno$group"))
  }

  if(!checkManateeObject(dataset,"Dataset")){stop("The function has failed to produce a valid MANATEE Dataset Object")}

  return(dataset)
}



###-###-###-###-###-###-##-#
###   makeDiscoValid()   ###
###-###-###-###-###-###-##-#

#DESCRIPTION
#This randomly splits the given list of datasets into discovery and validation
#(healthy and non-healthy samples are split separately) and then COCONUT
#normalizes the discovery and validation separately. Next, any designated
#samples are removed, and a "class" is created for both discovery and
#validation

#PARAMETERS
#GSElist - A MANATEE DatasetList Object
#conorm.type - The type of conormalization to perform. Options are COCONUT, COCOImputation, ComBat, or ComBatImputation (not case sensitive) (Default: "COCONUT")
#Impute.method - Which method to use for COCOImputation. Options are "simple", "quantile", and "reverse quantile" (Default: "simple")
#Impute.NAcutoff - Any genes with a higher proportion of NAs than the NAcutoff will be automatically removed. Note that if keepGenes is used, then the NAcutoff will be ignored (Default: 0.8)
#Impute.keepGenes - If you only want to perform COCOImputation for a certain set of genes, pass in a vector of those genes with this argument.
#                   If this argument is used, then only the genes that are common to all datasets and the genes in keepGenes will be included,
#                   and all other genes will be discarded (regardless of the \code{NAcutoff}).
#use.pheno.class - if TRUE, then the "class" column in $pheno will be used to generate the class vector for your discovery and validation. If FALSE, then makeClassVector will be used instead. (Default: FALSE)
#caseNames - name of the class(es) that you're considering to be your case (e.g. "Malaria")
#remove.samples - character vector of the labels that correspond to any samples you want to remove after running COCONUT. This is often done for control/healthy samples
#seed - random seed
#split.min.healthy - The minimum number of healthy samples needed for a dataset to be split. If there are fewer healthies, then all of the dataset will go in discovery.
#split.prop.discovery - The proportion of data to put in discovery. Default is 70%
#healthy.names - if your control/healthy samples are not labeled as "Healthy" in the group column of $pheno, you can indicate the alternative "Healthy" name(s) here
#forced.discovery - Any datasets that you want to put only in discovery (using the dataset name stored in the $pheno$name vector)
#forced.validation - Any datasets that you want to put only in validation (using the dataset name stored in the $pheno$name vector)
#remove.names - The labels used to designate a sample that is not going to be included in the analysis
#dummy.cols - if you want to preserve a column that isn't present in the $pheno of all datasets, you can add it here and dummy columns will be created in all datasets in order to preserve that column post-COCONUT

#RETURN VALUE
#A list containing the newly split Discovery data in one field and the Validation data in the other field

#REQUIRED PACKAGES: COCONUT (add a bunch from cocoimputation stuff)

#NOTE: If I ever get this error again (Error in X$genes[common, ] : incorrect number of dimensions) it's because of a drop error in COCONUT. I'm currently using debugCOCONUT to fix this issue.

makeDiscoValid <- function(GSElist, conorm.type = "COCONUT", Impute.method="simple", Impute.NAcutoff=0.8, Impute.keepGenes=NULL, use.pheno.class = FALSE, caseNames = NULL, remove.samples = NULL, seed = 1337, split.min.healthy = 10,
                           split.prop.discovery = 0.7, healthy.names = "Healthy", forced.discovery = NULL, forced.validation=NULL, remove.names = "Remove", dummy.cols = NULL){
  if(!checkManateeObject(GSElist,"DatasetList")){stop("Invalid DatasetList provided")}
  conorm.type = tolower(conorm.type)
  if(!(conorm.type %in% c("coconut","cocoimputation","combat","combatimputation"))){
    stop("method must be \"COCONUT\", \"COCOImputation\", \"ComBat\", or\"ComBatImputation\"")
  }
  if(!use.pheno.class && is.null(caseNames)){
    stop("If use.pheno.class is FALSE, then caseNames cannot be NULL")
  }
  #check for and remove all datasets that only have healthy samples
  #THIS HASN'T BEEN VERIFIED TO WORK YET
  display.rmdata = TRUE
  for(i in length(GSElist):1){
    if(all(as.character(GSElist[[i]]$pheno$group) %in% c(healthy.names,remove.names))){
      if(display.rmdata){cat(paste0("Removing ",names(GSElist)[i]," (only healthy)\n"))}
      GSElist[[i]] = NULL
    }else if(!any(as.character(GSElist[[i]]$pheno$group) %in% healthy.names)){
      if(display.rmdata){cat(paste0("Removing ",names(GSElist)[i]," (no healthy)\n"))}
      GSElist[[i]] = NULL
    }
  }

  if(!is.null(dummy.cols)){
    GSElist = lapply(GSElist, function(gse){
      for(col in dummy.cols){
        if(!col %in% colnames(gse$pheno)){
          gse$pheno[[col]] = NA
        }
      }
      return(gse)
    })
  }
  
  split.datasets = splitData(GSElist, split.min.healthy, split.prop.discovery, healthy.names, forced.discovery, forced.validation, seed)

  if(conorm.type == "coconut"){
    cat("Discovery\n")
    COCO.disco.out = debugCOCONUT(GSEs = split.datasets$Discovery, control.0.col = "control.0.class")
    Discovery = combineCOCOoutput(COCO.disco.out)
    cat("\nValidation\n")
    COCO.valid.out = debugCOCONUT(GSEs = split.datasets$Validation, control.0.col = "control.0.class")
    Validation = combineCOCOoutput(COCO.valid.out)
  }else if(conorm.type == "cocoimputation"){
    cat("Discovery\n")
    Discovery = runCOCOImputation(GSElist = split.datasets$Discovery, seed=seed, method=Impute.method, NAcutoff=Impute.NAcutoff, keepGenes=Impute.keepGenes)
    cat("\nValidation\n")
    Validation = runCOCOImputation(GSElist = split.datasets$Validation, seed=seed, method=Impute.method, NAcutoff=Impute.NAcutoff, keepGenes=Impute.keepGenes)
  }else if(conorm.type == "combat"){
    cat("Discovery\n")
    Discovery = .runComBat(GSElist = split.datasets$Discovery)
    cat("\nValidation\n")
    Validation = .runComBat(GSElist = split.datasets$Validation)
  }else if(conorm.type == "combatimputation"){
    cat("Discovery\n")
    Discovery = runCOCOImputation(GSElist = split.datasets$Discovery, seed=seed, use.ComBat = TRUE, method=Impute.method, NAcutoff=Impute.NAcutoff, keepGenes=Impute.keepGenes)
    cat("\nValidation\n")
    Validation = runCOCOImputation(GSElist = split.datasets$Validation, seed=seed, use.ComBat = TRUE, method=Impute.method, NAcutoff=Impute.NAcutoff, keepGenes=Impute.keepGenes)
  }

  if(!is.null(remove.samples)){
    for(name in remove.samples){
      remove.d = which(Discovery$pheno$group == name)
      if(length(remove.d) > 0){
        Discovery$pheno = Discovery$pheno[-remove.d,]
        Discovery$genes = Discovery$genes[,-remove.d]
      }
      remove.v = which(Validation$pheno$group == name)
      if(length(remove.v) > 0){
        Validation$pheno = Validation$pheno[-remove.v,]
        Validation$genes = Validation$genes[,-remove.v]
      }
    }
  }

  if(use.pheno.class){
    if(is.null(Discovery$pheno$control.0.class)){
      warning("There is no \"control.0.class\" column in $pheno in your normalized results, which probably means that the \"control.0.class\" column was missing in one of your input datasets. No class vector will be assigned.")
    }else{
      Discovery$class = as.numeric(as.character(Discovery$pheno$control.0.class))
      Validation$class = as.numeric(as.character(Validation$pheno$control.0.class))
    }
  }else{
    Discovery$class = makeClassVector(Discovery$pheno,caseNames=caseNames)
    Validation$class = makeClassVector(Validation$pheno,caseNames=caseNames)
  }
  
  if(class(Discovery$genes) == "data.frame"){
    Discovery$genes = data.matrix(Discovery$genes)
  }
  if(class(Validation$genes) == "data.frame"){
    Validation$genes = data.matrix(Validation$genes)
  }

  return(list(Discovery=Discovery,Validation=Validation))
}



###-###-###-###-###-###-#
###   runLasso.cv()   ###
###-###-###-###-###-###-#

#DESCRIPTION
#This runs a cross-validated lasso on the given pooledDataObject using the genes
#in the filterObject and returns the resulting lasso object. It also outputs an
#AUC vs. lambda plot describing the results of the lasso

#PARAMETERS
#pooledDataObject - A MANATEE pooledDataObject. If weighting is desired, then sample weights must be stored in pooledDataObject$sample.weights
#filterObject - a MANATEE filterObject
#weighted - if TRUE, weighting is performed using the weights in pooledDataObject$sample.weights. If not specified, then weighting will automatically be performed if \code{pooledDataObject$sample.weights} is populated,
#           but weighting will not be performed if \code{pooledDataObject$sample.weights} is NULL
#nfolds - number of folds; default is 10. Although nfolds can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets. Smallest value allowable is nfolds=3
#alpha - The elasticnet mixing parameter, with 0 ≤ a ≤ 1. alpha=1 is the lasso penalty, and alpha=0 the ridge penalty
#nlambda - The number of lambda values; default is 100
#title - title of the AUC vs. lambda plot
#plot - if TRUE, the AUC vs. lambda plot is made
#... - other arguments for cv.glmnet()

#RETURN VALUE
#Returns a modified version of the \code{pooledDataObject} with the lasso results stored in \code{pooledDataObject$lassoResults}
#lasso.obj - the output of lasso, modified to include extra fields: "call" contains the initial function call and "topGenesMtx" contains the gene matrix used for training the lasso

#REQUIRED PACKAGES: glmnet

runLasso.cv <- function(pooledDataObject, filterObject, weighted = NULL, nfolds = 10, alpha = 1, nlambda = 100, title=NULL, plot=TRUE, ...){
  call = sys.call(which = 0)
  if(!checkManateeObject(pooledDataObject,"pooledDataObject")){stop("Invalid pooledDataObject provided")}
  if(!checkManateeObject(filterObject,"filterObject")){stop("Invalid filterObject provided")}
  topGeneNames = c(filterObject$upGeneNames,filterObject$downGeneNames)
  if(length(topGeneNames)==0){stop("Both upGeneNames and downGeneNames are empty in the provided filterObject")}
  if(length(topGeneNames)==1){stop("There must be more than one gene in the filterObject in order to run lasso")}

  if(is.null(pooledDataObject$sample.weights) && is.null(weighted)){weighted=FALSE}
  if(!is.null(pooledDataObject$sample.weights) && is.null(weighted)){weighted=TRUE}

  topGenesMtx = pooledDataObject$genes[topGeneNames,,drop=F]

  require(glmnet)
  if(weighted){
    lasso.obj = cv.glmnet(t(topGenesMtx), pooledDataObject$class, family = "binomial", type.measure="auc",
                          weights=pooledDataObject$sample.weights, nfolds=nfolds, alpha=alpha, nlambda=nlambda, ...)
  } else{
    lasso.obj = cv.glmnet(t(topGenesMtx), pooledDataObject$class, family = "binomial", type.measure="auc",
                          nfolds=nfolds, alpha=alpha, nlambda=nlambda, ...)
  }

  #AUC vs genes plot
  if(plot){
    plot(lasso.obj)
    title(title, line = 2.5)
    mtext("Number of Genes")
  }

  lasso.obj$call = call
  lasso.obj$topGenesMtx = topGenesMtx
  pooledDataObject$lassoResults = lasso.obj
  return(pooledDataObject)
}



###-###-###-###-###-###-###-##
###   getTopLassoGenes()   ###
###-###-###-###-###-###-###-##

#DESCRIPTION
#once you've decided how many genes you want to get out of LASSO, you can run
#this function to get the optimal set of genes, as well as the name and the
#direction (i.e. upregulated or downregulated) of each gene. This will also plot
#the ROC and PRC curves for those genes, with the LASSO model as well as with
#geometric mean score

#PARAMETERS
#pooledDataObject - A MANATEE pooledDataObject
#numGenes - the number of genes that you want to get
#returnType - If \code{returnType}="pooledDataObject", then the function will return a modified version of the \code{pooledDataObject} with the top lasso genes stored in
#             \code{pooledDataObject$lassoResults$topLassoGenes} under a unique identifier.
#             If \code{returnType}="topLassoGenes", then the function will just return a list that describes the top lasso genes (Default: "pooledDataObject")
#plot - if TRUE, then the ROC and PRC curves are plotted
#cex.legend - sometimes the legend text gets cut off, so use this to decrease the legend text size if that happens

#RETURN VALUE
#upgenes - character vector of the upregulated genes
#downgenes - character vector of the downregulated genes
#nonzeroGenes - named vector of the coefficients of the non-zero genes
#lambda - the lambda that was used to get these genes - this can be used for plotting

#REQUIRED PACKAGES: ROCR, MetaIntegrator

getTopLassoGenes <- function(pooledDataObject, numGenes, returnType = "pooledDataObject", plot = TRUE, cex.legend=1){
  if(!checkManateeObject(pooledDataObject,"pooledDataObject")){stop("Invalid pooledDataObject provided")}
  if(is.null(pooledDataObject$lassoResults)){stop("pooledDataObject$lassoResults is empty")}
  if(!(returnType %in% c("pooledDataObject","topLassoGenes"))){
    stop("returnType must either be \"pooledDataObject\" or \"topLassoGenes\"")
  }

  gene_nums = getGeneNums.lasso(lasso.obj=pooledDataObject$lassoResults,lambda.vals=pooledDataObject$lassoResults$lambda)
  cut = which(gene_nums>=numGenes)
  cut = which(gene_nums == min(gene_nums[cut]))
  cut = max(cut) #note: taking the max here is only optimal for lower number of genes
  if(gene_nums[cut] != numGenes){
    warning(sprintf("There was no lambda that corresponded to %s genes; chose the best %s gene set instead.",numGenes,gene_nums[cut]))
  }
  topLassoGenes = .getNonzeroGenes.lasso(lasso.obj=pooledDataObject$lassoResults,lambda=pooledDataObject$lassoResults$lambda[cut])
  topLassoGenes$lambda = pooledDataObject$lassoResults$lambda[cut]

  if(plot){
    #make plots
    par.mfrow = par("mfrow")
    par.mar = par("mar")
    par(mfrow=c(2,2))
    par(mar=c(4.3,4.1,2.3,1.1))
    makeROCplot.lasso(pooledDataObject,lambda.vals=pooledDataObject$lassoResults$lambda[cut],title="ROC plot using CV lasso model - Discovery",cex.legend=cex.legend)
    makePRCplot.lasso(pooledDataObject,lambda.vals=pooledDataObject$lassoResults$lambda[cut],title="PRC plot using CV lasso model - Discovery",cex.legend=cex.legend)
    
    scores = getGeneScores(pooledDataObject$genes,topLassoGenes$upgenes,topLassoGenes$downgenes)
    .makeROCplotHelper(prediction.type="geometric.smd",class=pooledDataObject$class,gene.nums = gene_nums[cut],
                       scoreList=list(scores),title="ROC plot using geometric mean score - Discovery",cex.legend=cex.legend)
    .makePRCplotHelper(prediction.type="geometric.smd",class=pooledDataObject$class,gene.nums = gene_nums[cut],
                       scoreList=list(scores),title="PRC plot using geometric mean score - Discovery",cex.legend=cex.legend)
    par(mfrow = par.mfrow)
    par(mar = par.mar)
  }

  if(returnType == "pooledDataObject"){
    if("topLassoGenes" %in% names(pooledDataObject$lassoResults)){
      pooledDataObject$lassoResults$topLassoGenes[[as.character(gene_nums[cut])]]=topLassoGenes
    }else{
      pooledDataObject$lassoResults$topLassoGenes = list()
      pooledDataObject$lassoResults$topLassoGenes[[as.character(gene_nums[cut])]]=topLassoGenes
    }
    return(pooledDataObject)
  }else if(returnType=="topLassoGenes"){
    return(topLassoGenes)
  }
}



###-###-###-###-###-###-##
###   getSigScores()   ###
###-###-###-###-###-###-##

#DESCRIPTION
#Given a pooledDataObject and a filterObject, this function calculates the
#signature scores (i.e. geometric mean scores) for each sample.

#PARAMETERS
#pooledDataObject - A MANATEE \code{pooledDataObject}
#filterObject - A MANATEE \code{filterObject}
#makePos - If TRUE, then if there are any negative values in your $genes, then the minimum value will be subtracted from all values to make every value positive
#out.missing - if TRUE, the names of the missing genes will be output

#RETURN VALUE
#A vector of signature scores for each sample

#REQUIRED PACKAGES: None

getSigScores <- function(pooledDataObject, filterObject, makePos = TRUE, remove.missing=TRUE, out.missing=TRUE){
  args <- sys.call()
  pos = filterObject$upGeneNames
  neg = filterObject$downGeneNames
  geneMtx = pooledDataObject$genes
  if(out.missing){
    missingpos = pos[!(pos %in% rownames(geneMtx))]
    missingneg = neg[!(neg %in% rownames(geneMtx))]
    missing = c(missingpos,missingneg)
    if(length(missing)>0){
      cat(paste0("Missing these genes: ",missing," (from ",as.character(args[-1][2]),")\n"))
    }
  }
  pos = pos[pos %in% rownames(geneMtx)]
  neg = neg[neg %in% rownames(geneMtx)]

  if(makePos){ #If I need to make everything positive
    if(any(geneMtx < 0, na.rm=T)){
      geneMtx <- geneMtx - min(geneMtx, na.rm=T)
    }
  }

  if(length(neg)>=1 && length(pos)>=1){
    scores <- apply(geneMtx[pos, , drop=F], 2, geomMean, na.rm=T) - apply(geneMtx[neg, , drop=F], 2, geomMean, na.rm=T)
  }else if(length(pos)>=1){
    scores <- apply(geneMtx[pos, , drop=F], 2, geomMean, na.rm=T)
  }else if(length(neg)>=1){
    scores <- -1*apply(geneMtx[neg, , drop=F], 2, geomMean, na.rm=T)
  }
  return(scores)
}



###-###-###-###-###-###-###-##
###   getSigScores.Tim()   ###
###-###-###-###-###-###-###-##

#DESCRIPTION
#Given a pooledDataObject and a filterObject, this function calculates the
#signature scores (i.e. geometric mean scores) for each sample. This uses Tim's score calculation method.

#PARAMETERS
#pooledDataObject - A MANATEE \code{pooledDataObject}
#filterObject - A MANATEE \code{filterObject}
#makePos - If TRUE, then if there are any negative values in your $genes, then the minimum value will be subtracted from all values to make every value positive
#out.missing - if TRUE, the names of the missing genes will be output

#RETURN VALUE
#A vector of signature scores for each sample

#REQUIRED PACKAGES: None

getSigScores.Tim <- function(pooledDataObject, filterObject, makePos = TRUE, remove.missing=TRUE, out.missing=TRUE){
  args <- sys.call()
  pos = filterObject$upGeneNames
  neg = filterObject$downGeneNames
  geneMtx = pooledDataObject$genes
  if(out.missing){
    missingpos = pos[!(pos %in% rownames(geneMtx))]
    missingneg = neg[!(neg %in% rownames(geneMtx))]
    missing = c(missingpos,missingneg)
    if(length(missing)>0){
      cat(paste0("Missing these genes: ",missing," (from ",as.character(args[-1][2]),")\n"))
    }
  }
  pos = pos[pos %in% rownames(geneMtx)]
  neg = neg[neg %in% rownames(geneMtx)]
  
  if(makePos){ #If I need to make everything positive
    if(any(geneMtx < 0, na.rm=T)){
      geneMtx <- geneMtx - min(geneMtx, na.rm=T)
    }
  }
  
  if(length(neg)>=1 && length(pos)>=1){
    scores <- apply(geneMtx[pos, , drop=F], 2, geomMean, na.rm=T) - ratio * apply(geneMtx[neg, , drop=F], 2, geomMean, na.rm=T)
  }else if(length(pos)>=1){
    scores <- apply(geneMtx[pos, , drop=F], 2, geomMean, na.rm=T)
  }else if(length(neg)>=1){
    scores <- -1*apply(geneMtx[neg, , drop=F], 2, geomMean, na.rm=T)
  }
  return(scores)
}



###-###-###-###-###-###-###-###
###   makeROCplot.lasso()   ###
###-###-###-###-###-###-###-###

#DESCRIPTION
#Wrapper function to create a ROC plot from a pooledDataObject using a lasso
#object (stored in pooledDataObject$lassoResults) to make predictions

#PARAMETERS
#pooledDataObject - A MANATEE \code{pooledDataObject}
#lambda.vals - the lambda values at which you want to make predictions - if left blank, then lassoResults$lambda.min and lassoResults$lambda.1se will be used by default
#cross.validation - set as TRUE if cross-validated lasso was run (using cv.glmnet) and FALSE if a regular lasso was run

#GRAPHICAL PARAMETERS
#title - title of the figure
#subtitle - subtitle of the figure
#textSize - use this to easily increase or decrease the size of all the text in the plot
#rounding - how many digits to round the AUC and CI to
#colors - vector of colors for the lines
#legend - if true, a legend will be included

#HELPER FUNCTION PARAMETERS
#ROC.lty - line type
#ROC.lwd - line width
#background.color - background color of the plot
#grid.marks - increment between grid lines
#grid.color - grid line color
#grid.lty - grid line type
#grid.lwd - grid line width
#diag.lty - diagonal line style
#legend.lty - legend style (0 is no box, 1 is boxed legend)
#legend.color - color of the legend box
#cex.main - title size
#cex.subtitle - subtitle size
#cex.legend - legend text size

#RETURN VALUE
#A ROC plot showing diagnostic performance

#REQUIRED PACKAGES: ROCR, MetaIntegrator, glmnet

makeROCplot.lasso <- function(pooledDataObject, lambda.vals=NULL, cross.validation=TRUE, title=NULL, subtitle=NULL, textSize=NULL,
                              rounding=3, colors=c(2,4,6,5,3,7:10), legend=TRUE, ...){
  if(!checkManateeObject(pooledDataObject,"pooledDataObject")){stop("Invalid pooledDataObject provided")}
  if(is.null(pooledDataObject$lassoResults)){stop("pooledDataObject$lassoResults is empty")}

  if(is.null(lambda.vals)){
    lambda.vals=c(pooledDataObject$lassoResults$lambda.min, pooledDataObject$lassoResults$lambda.1se)
  }
  if(cross.validation){
    if(is.null(pooledDataObject$lassoResults$cvm)){stop("This is not a cross-validated lasso object")}
    .makeROCplotHelper(prediction.type="cv.lasso", class=pooledDataObject$class, lasso.obj=pooledDataObject$lassoResults, geneMtx=pooledDataObject$lassoResults$topGenesMtx, lambda.vals=lambda.vals, scoreList=NULL,
                       gene.nums=NULL, title=title, subtitle=subtitle, textSize=textSize, rounding=rounding, colors=colors, legend=legend, ...)
  }else{
    .makeROCplotHelper(prediction.type="lasso", class=pooledDataObject$class, lasso.obj=pooledDataObject$lassoResults, geneMtx=pooledDataObject$lassoResults$topGenesMtx, lambda.vals=lambda.vals, scoreList=NULL,
                       gene.nums=NULL, title=title, subtitle=subtitle, textSize=textSize, rounding=rounding, colors=colors, legend=legend, ...)
  }
}



###-###-###-###-###-###-###-###
###   makePRCplot.lasso()   ###
###-###-###-###-###-###-###-###

#DESCRIPTION
#Wrapper function to create a PRC plot from a pooledDataObject using a lasso
#object (stored in pooledDataObject$lassoResults) to make predictions

#PARAMETERS
#pooledDataObject - A MANATEE \code{pooledDataObject}
#lambda.vals - the lambda values at which you want to make predictions - if left blank, then lassoResults$lambda.min and lassoResults$lambda.1se will be used by default
#cross.validation - set as TRUE if cross-validated lasso was run (using cv.glmnet) and FALSE if a regular lasso was run

#GRAPHICAL PARAMETERS
#title - title of the figure
#subtitle - subtitle of the figure
#textSize - use this to easily increase or decrease the size of all the text in the plot
#rounding - how many digits to round the AUC and CI to
#colors - vector of colors for the lines
#legend - if true, a legend will be included

#HELPER FUNCTION PARAMETERS
#PRC.lty - line type
#PRC.lwd - line width
#background.color - background color of the plot
#grid.marks - increment between grid lines
#grid.color - grid line color
#grid.lty - grid line type
#grid.lwd - grid line width
#legend.lty - legend style (0 is no box, 1 is boxed legend)
#cex.main - title size
#cex.subtitle - subtitle size
#cex.legend - legend text size

#RETURN VALUE
#A PRC plot showing diagnostic performance

#REQUIRED PACKAGES: ROCR, MetaIntegrator, glmnet

makePRCplot.lasso <- function(pooledDataObject, lambda.vals=NULL, cross.validation=TRUE, title=NULL, subtitle=NULL, textSize=NULL,
                              rounding=3, colors=c(2,4,6,5,3,7:10), legend=TRUE, ...){
  if(!checkManateeObject(pooledDataObject,"pooledDataObject")){stop("Invalid pooledDataObject provided")}
  if(is.null(pooledDataObject$lassoResults)){stop("pooledDataObject$lassoResults is empty")}

  if(is.null(lambda.vals)){
    lambda.vals=c(pooledDataObject$lassoResults$lambda.min, pooledDataObject$lassoResults$lambda.1se)
  }
  if(cross.validation){
    if(is.null(pooledDataObject$lassoResults$cvm)){stop("This is not a cross-validated lasso object")}
    .makePRCplotHelper(prediction.type="cv.lasso", class=pooledDataObject$class, lasso.obj=pooledDataObject$lassoResults, geneMtx=pooledDataObject$lassoResults$topGenesMtx, lambda.vals=lambda.vals, scoreList=NULL,
                       gene.nums=NULL, title=title, subtitle=subtitle, textSize=textSize, rounding=rounding, colors=colors, legend=legend, ...)
  }else{
    .makePRCplotHelper(prediction.type="lasso", class=pooledDataObject$class, lasso.obj=pooledDataObject$lassoResults, geneMtx=pooledDataObject$lassoResults$topGenesMtx, lambda.vals=lambda.vals, scoreList=NULL,
                       gene.nums=NULL, title=title, subtitle=subtitle, textSize=textSize, rounding=rounding, colors=colors, legend=legend, ...)
  }
}



###-###-###-###-###-###-###-###-##-#
###   makeClassROCplot.lasso()   ###
###-###-###-###-###-###-###-###-##-#

#DESCRIPTION
#Wrapper function to create a Class ROC plot (a plot with multiple ROC curves
#comparing your positive class (i.e. your cases) against each of your other
#classes (i.e. your controls)) from a pooledDataObject using a lasso object
#(stored in pooledDataObject$lassoResults) to make predictions

#PARAMETERS
#pooledDataObject - A MANATEE \code{pooledDataObject}
#lambda.vals - the lambda values at which you want to make predictions - if left blank, then lassoResults$lambda.1se will be used by default
#cross.validation - set as TRUE if cross-validated lasso was run (using cv.glmnet) and FALSE if a regular lasso was run
#caseNames - if left blank, this will be a character vector of the class(es) that correspond to your cases (i.e. samples labeled with a 1 in pooledDataObject$class).
#            If the user wants to exclude any of these case classes, then this can be set as the names of all of the case classes that should be included
#otherNames - if left blank, this will be a character vector of the class(es) that correspond to your controls (i.e. samples labeled with a 0 in pooledDataObject$class).
#             If the user wants to exclude any of these control classes, then this can be set as the names of all of the control classes that should be included

#GRAPHICAL PARAMETERS
#title - title of the figure
#textSize - use this to easily increase or decrease the size of all the text in the plot
#full.labels - if TRUE, then the full names of the case categories will be used in the legend, if FALSE then a generic "Cases" will be used instead
#rounding - how many digits to round the AUC and CI to
#colors - vector of colors for the lines
#plot.unweighted - should the unweighted ROC curve be plotted (i.e. the normal ROC curve without class divisions)

#HELPER FUNCTION PARAMETERS
#lwd - vector of line widths. if NULL, the lines will all have a width of 1
#unweighted.lwd - width of the unweighted curve
#unweighted.col - color of the unweighted curve
#background.color - background color of the plot
#grid.marks - increment between grid lines
#grid.color - grid line color
#grid.lty - grid line type
#grid.lwd - grid line width
#diag.lty - diagonal line style
#legend.lty - legend style (0 is no box, 1 is boxed legend)
#legend.color - color of the legend box
#cex.main - title size
#cex.subtitle - subtitle size
#cex.legend - legend text size

#RETURN VALUE
#A ROC plot showing diagnostic performance across multiple designated classes

#REQUIRED PACKAGES: ROCR, MetaIntegrator, glmnet

makeClassROCplot.lasso <- function(pooledDataObject, lambda.vals=NULL, cross.validation=TRUE, caseNames=NULL, otherNames=NULL, title=NULL, textSize=NULL, full.labels=TRUE,
                                   rounding=3, colors=c(2,4,6,5,3,7:10), plot.unweighted=TRUE, ...){
  if(!checkManateeObject(pooledDataObject,"pooledDataObject")){stop("Invalid pooledDataObject provided")}
  if(is.null(pooledDataObject$lassoResults)){stop("pooledDataObject$lassoResults is empty")}

  if(is.null(lambda.vals)){
    lambda.vals=pooledDataObject$lassoResults$lambda.1se
  }
  if(cross.validation){
    if(is.null(pooledDataObject$lassoResults$cvm)){stop("This is not a cross-validated lasso object")}
    .makeClassROCplotHelper(prediction.type="cv.lasso", class=pooledDataObject$class, lasso.obj=pooledDataObject$lassoResults, geneMtx=pooledDataObject$lassoResults$topGenesMtx, lambda.vals=lambda.vals,
                            pos.genes=NULL, neg.genes=NULL, labelVec=pooledDataObject$pheno$group, caseNames=caseNames, otherNames=otherNames, title=title, textSize=textSize, full.labels = full.labels,
                            rounding=rounding, colors=colors, plot.unweighted=plot.unweighted, ...)
  }else{
    .makeClassROCplotHelper(prediction.type="lasso", class=pooledDataObject$class, lasso.obj=pooledDataObject$lassoResults, geneMtx=pooledDataObject$lassoResults$topGenesMtx, lambda.vals=lambda.vals,
                            pos.genes=NULL, neg.genes=NULL, labelVec=pooledDataObject$pheno$group, caseNames=caseNames, otherNames=otherNames, title=title, textSize=textSize, full.labels = full.labels,
                            rounding=rounding, colors=colors, plot.unweighted=plot.unweighted, ...)
  }
}



###-###-###-###-###-###-###-###-##-#
###   makeClassPRCplot.lasso()   ###
###-###-###-###-###-###-###-###-##-#

#DESCRIPTION
#Wrapper function to create a Class PRC plot (a plot with multiple PRC curves
#comparing your positive class (i.e. your cases) against each of your other
#classes (i.e. your controls)) from a pooledDataObject using a lasso object
#(stored in pooledDataObject$lassoResults) to make predictions

#PARAMETERS
#pooledDataObject - A MANATEE \code{pooledDataObject}
#lambda.vals - the lambda values at which you want to make predictions - if left blank, then lassoResults$lambda.1se will be used by default
#cross.validation - set as TRUE if cross-validated lasso was run (using cv.glmnet) and FALSE if a regular lasso was run
#caseNames - if left blank, this will be a character vector of the class(es) that correspond to your cases (i.e. samples labeled with a 1 in pooledDataObject$class).
#            If the user wants to exclude any of these case classes, then this can be set as the names of all of the case classes that should be included
#otherNames - if left blank, this will be a character vector of the class(es) that correspond to your controls (i.e. samples labeled with a 0 in pooledDataObject$class).
#             If the user wants to exclude any of these control classes, then this can be set as the names of all of the control classes that should be included

#GRAPHICAL PARAMETERS
#title - title of the figure
#textSize - use this to easily increase or decrease the size of all the text in the plot
#full.labels - if TRUE, then the full names of the case categories will be used in the legend, if FALSE then a generic "Cases" will be used instead
#rounding - how many digits to round the AUC and CI to
#colors - vector of colors for the lines
#plot.unweighted - should the unweighted PRC curve be plotted (i.e. the normal PRC curve without class divisions)

#HELPER FUNCTION PARAMETERS
#lwd - vector of line widths. if NULL, the lines will all have a width of 1
#unweighted.lwd - width of the unweighted curve
#unweighted.col - color of the unweighted curve
#background.color - background color of the plot
#grid.marks - increment between grid lines
#grid.color - grid line color
#grid.lty - grid line type
#grid.lwd - grid line width
#legend.lty - legend style (0 is no box, 1 is boxed legend)
#legend.color - color of the legend box
#cex.main - title size
#cex.subtitle - subtitle size
#cex.legend - legend text size

#RETURN VALUE
#A PRC plot showing diagnostic performance across multiple designated classes

#REQUIRED PACKAGES: ROCR, MetaIntegrator, glmnet, zoo, pracma

makeClassPRCplot.lasso <- function(pooledDataObject, lambda.vals=NULL, cross.validation=TRUE, caseNames=NULL, otherNames=NULL, title=NULL, textSize=NULL, full.labels=TRUE,
                                   rounding=3, colors=c(2,4,6,5,3,7:10), plot.unweighted=TRUE, ...){
  if(!checkManateeObject(pooledDataObject,"pooledDataObject")){stop("Invalid pooledDataObject provided")}
  if(is.null(pooledDataObject$lassoResults)){stop("pooledDataObject$lassoResults is empty")}

  if(is.null(lambda.vals)){
    lambda.vals=pooledDataObject$lassoResults$lambda.1se
  }
  if(cross.validation){
    if(is.null(pooledDataObject$lassoResults$cvm)){stop("This is not a cross-validated lasso object")}
    .makeClassPRCplotHelper(prediction.type="cv.lasso", class=pooledDataObject$class, lasso.obj=pooledDataObject$lassoResults, geneMtx=pooledDataObject$lassoResults$topGenesMtx, lambda.vals=lambda.vals,
                            pos.genes=NULL, neg.genes=NULL, labelVec=pooledDataObject$pheno$group, caseNames=caseNames, otherNames=otherNames, title=title, textSize=textSize, full.labels = full.labels,
                            rounding=rounding, colors=colors, plot.unweighted=plot.unweighted, ...)
  }else{
    .makeClassPRCplotHelper(prediction.type="lasso", class=pooledDataObject$class, lasso.obj=pooledDataObject$lassoResults, geneMtx=pooledDataObject$lassoResults$topGenesMtx, lambda.vals=lambda.vals,
                            pos.genes=NULL, neg.genes=NULL, labelVec=pooledDataObject$pheno$group, caseNames=caseNames, otherNames=otherNames, title=title, textSize=textSize, full.labels = full.labels,
                            rounding=rounding, colors=colors, plot.unweighted=plot.unweighted, ...)
  }
}



###-###-###-###-###-###-#
###   makeROCplot()   ###
###-###-###-###-###-###-#

#DESCRIPTION
#Wrapper function to create a ROC plot from a pooledDataObject using the genes
#from a filterObject (or multiple filterObjects) to make predictions with the
#geometric mean score model

#PARAMETERS
#pooledDataObject - A MANATEE \code{pooledDataObject}
#filterObjectList - Either a single MANATEE \code{filterObject}, or a list of \code{filterObject}s
#upgenes - if desired, a character vector of the upgenes can be provided instead of filterObject(s). This can either be a single character vector, or a list of character vectors
#downgenes - if desired, a character vector of the downgenes can be provided instead of filterObject(s). This can either be a single character vector, or a list of character vectors

#GRAPHICAL PARAMETERS
#title - title of the figure
#subtitle - subtitle of the figure
#legend.names - if you don't want the default legend labels (i.e. the number of genes in each filterObject) then you can provide alternative names as a character vector
#textSize - use this to easily increase or decrease the size of all the text in the plot
#rounding - how many digits to round the AUC and CI to
#colors - vector of colors for the lines
#legend - if true, a legend will be included

#HELPER FUNCTION PARAMETERS
#ROC.lty - line type
#ROC.lwd - line width
#background.color - background color of the plot
#grid.marks - increment between grid lines
#grid.color - grid line color
#grid.lty - grid line type
#grid.lwd - grid line width
#diag.lty - diagonal line style
#legend.lty - legend style (0 is no box, 1 is boxed legend)
#legend.color - color of the legend box
#cex.main - title size
#cex.subtitle - subtitle size
#cex.legend - legend text size

#RETURN VALUE
#A ROC plot showing diagnostic performance

#REQUIRED PACKAGES: ROCR, MetaIntegrator, glmnet

makeROCplot <- function(pooledDataObject, filterObjectList=NULL, upgenes=NULL, downgenes=NULL, title=NULL, subtitle=NULL, legend.names=NULL, textSize=NULL,
                        rounding=3, colors=c(2,4,6,5,3,7:10), legend=TRUE, ...){
  if(!checkManateeObject(pooledDataObject,"pooledDataObject")){stop("Invalid pooledDataObject provided")}

  if(is.null(filterObjectList)){
    if(is.null(upgenes) && is.null(downgenes)){stop("If filterObjectList is null, then upgenes and downgenes must be provided")}
    if((class(upgenes)=="character" && class(downgenes)=="character") || (is.null(upgenes) && class(downgenes)=="character") ||
       (class(upgenes)=="character" && is.null(downgenes))){
      scoreList = list(getGeneScores(pooledDataObject$genes,upgenes,downgenes))
      gene.nums = length(upgenes) + length(downgenes)
    }else if(class(upgenes)=="list" && class(downgenes)=="list"){
      if(length(upgenes) != length(downgenes)){stop("If upgenes and downgenes are lists of character vectors, they must both be the same size")}
      scoreList = list()
      gene.nums = rep(0,length(upgenes))
      for(i in 1:length(upgenes)){
        scoreList[[i]] = getGeneScores(pooledDataObject$genes,upgenes[[i]],downgenes[[i]])
        gene.nums[i] = length(upgenes[[i]])+length(downgenes[[i]])
      }
    }else{
      stop("upgenes and downgenes must either both be character vectors or they must both be equally sized lists of character vectors")
    }
  }else{
    if(!is.null(filterObjectList$upGeneNames)){ #i.e. a single filterObject
      if(!checkManateeObject(filterObjectList,"filterObject")){stop("Invalid filterObject provided")}
      filterObjectList = list(filterObjectList)
    }
    scoreList = list()
    gene.nums = rep(0,length(filterObjectList))
    for(i in 1:length(filterObjectList)){
      if(!checkManateeObject(filterObjectList[[i]],"filterObject")){stop(paste("filterObject",i,"in filterObjectList is invalid"))}
      scoreList[[i]] = getSigScores(pooledDataObject,filterObjectList[[i]])
      gene.nums[i] = length(filterObjectList[[i]]$upGeneNames)+length(filterObjectList[[i]]$downGeneNames)
    }
  }



  .makeROCplotHelper(prediction.type="geometric.smd", class=pooledDataObject$class, lasso.obj=pooledDataObject$lassoResults, geneMtx=pooledDataObject$genes, lambda.vals=NULL, scoreList=scoreList,
                     gene.nums=gene.nums, title=title, subtitle=subtitle, legend.names=legend.names, textSize=textSize, rounding=rounding, colors=colors, legend=legend, ...)
}



###-###-###-###-###-###-#
###   makePRCplot()   ###
###-###-###-###-###-###-#

#DESCRIPTION
#Wrapper function to create a PRC plot from a pooledDataObject using the genes
#from a filterObject (or multiple filterObjects) to make predictions with the
#geometric mean score model

#PARAMETERS
#pooledDataObject - A MANATEE \code{pooledDataObject}
#filterObjectList - Either a single MANATEE \code{filterObject}, or a list of \code{filterObject}s
#upgenes - if desired, a character vector of the upgenes can be provided instead of filterObject(s). This can either be a single character vector, or a list of character vectors
#downgenes - if desired, a character vector of the downgenes can be provided instead of filterObject(s). This can either be a single character vector, or a list of character vectors

#GRAPHICAL PARAMETERS
#title - title of the figure
#subtitle - subtitle of the figure
#legend.names - if you don't want the default legend labels (i.e. the number of genes in each filterObject) then you can provide alternative names as a character vector
#textSize - use this to easily increase or decrease the size of all the text in the plot
#rounding - how many digits to round the AUC and CI to
#colors - vector of colors for the lines
#legend - if true, a legend will be included

#HELPER FUNCTION PARAMETERS
#PRC.lty - line type
#PRC.lwd - line width
#background.color - background color of the plot
#grid.marks - increment between grid lines
#grid.color - grid line color
#grid.lty - grid line type
#grid.lwd - grid line width
#legend.lty - legend style (0 is no box, 1 is boxed legend)
#legend.color - color of the legend box
#cex.main - title size
#cex.subtitle - subtitle size
#cex.legend - legend text size

#RETURN VALUE
#A PRC plot showing diagnostic performance

#REQUIRED PACKAGES: ROCR, MetaIntegrator, glmnet

makePRCplot <- function(pooledDataObject, filterObjectList=NULL, upgenes=NULL, downgenes=NULL, title=NULL, subtitle=NULL, legend.names=NULL, textSize=NULL,
                        rounding=3, colors=c(2,4,6,5,3,7:10), legend=TRUE, ...){
  if(!checkManateeObject(pooledDataObject,"pooledDataObject")){stop("Invalid pooledDataObject provided")}


  if(is.null(filterObjectList)){
    if(is.null(upgenes) && is.null(downgenes)){stop("If filterObjectList is null, then upgenes and downgenes must be provided")}
    if(class(upgenes)=="character" && class(downgenes)=="character" || (is.null(upgenes) && class(downgenes)=="character") ||
       (class(upgenes)=="character" && is.null(downgenes))){
      scoreList = list(getGeneScores(pooledDataObject$genes,upgenes,downgenes))
      gene.nums = length(upgenes) + length(downgenes)
    }else if(class(upgenes)=="list" && class(downgenes)=="list"){
      if(length(upgenes) != length(downgenes)){stop("If upgenes and downgenes are lists of character vectors, they must both be the same size")}
      scoreList = list()
      gene.nums = rep(0,length(upgenes))
      for(i in 1:length(upgenes)){
        scoreList[[i]] = getGeneScores(pooledDataObject$genes,upgenes[[i]],downgenes[[i]])
        gene.nums[i] = length(upgenes[[i]])+length(downgenes[[i]])
      }
    }else{
      stop("upgenes and downgenes must either both be character vectors or they must both be equally sized lists of character vectors")
    }
  }else{
    if(!is.null(filterObjectList$upGeneNames)){ #i.e. a single filterObject
      if(!checkManateeObject(filterObjectList,"filterObject")){stop("Invalid filterObject provided")}
      filterObjectList = list(filterObjectList)
    }
    scoreList = list()
    gene.nums = rep(0,length(filterObjectList))
    for(i in 1:length(filterObjectList)){
      if(!checkManateeObject(filterObjectList[[i]],"filterObject")){stop(paste("filterObject",i,"in filterObjectList is invalid"))}
      scoreList[[i]] = getSigScores(pooledDataObject,filterObjectList[[i]])
      gene.nums[i] = length(filterObjectList[[i]]$upGeneNames)+length(filterObjectList[[i]]$downGeneNames)
    }
  }

  .makePRCplotHelper(prediction.type="geometric.smd", class=pooledDataObject$class, lasso.obj=pooledDataObject$lassoResults, geneMtx=pooledDataObject$genes, lambda.vals=NULL, scoreList=scoreList,
                     gene.nums=gene.nums, title=title, subtitle=subtitle, legend.names=legend.names, textSize=textSize, rounding=rounding, colors=colors, legend=legend, ...)
}



###-###-###-###-###-###-###-##
###   makeClassROCplot()   ###
###-###-###-###-###-###-###-##

#DESCRIPTION
#Wrapper function to create a Class ROC plot (a plot with multiple ROC curves
#comparing your positive class (i.e. your cases) against each of your other
#classes (i.e. your controls)) from a pooledDataObject using the genes
#from a filterObject (or multiple filterObjects) to make predictions with the
#geometric mean score model

#PARAMETERS
#pooledDataObject - A MANATEE \code{pooledDataObject}
#filterObject - A MANATEE \code{filterObject}
#upgenes - if desired, a character vector of the upgenes can be provided instead of a filterObject
#downgenes - if desired, a character vector of the downgenes can be provided instead of a filterObject
#caseNames - if left blank, this will be a character vector of the class(es) that correspond to your cases (i.e. samples labeled with a 1 in pooledDataObject$class).
#            If the user wants to exclude any of these case classes, then this can be set as the names of all of the case classes that should be included
#otherNames - if left blank, this will be a character vector of the class(es) that correspond to your controls (i.e. samples labeled with a 0 in pooledDataObject$class).
#             If the user wants to exclude any of these control classes, then this can be set as the names of all of the control classes that should be included
#group.colname - this designates the column name of the group vector within the pooledDataObject$pheno dataframe

#GRAPHICAL PARAMETERS
#title - title of the figure
#textSize - use this to easily increase or decrease the size of all the text in the plot
#full.labels - if TRUE, then the full names of the case categories will be used in the legend, if FALSE then a generic "Cases" will be used instead
#rounding - how many digits to round the AUC and CI to
#colors - vector of colors for the lines
#plot.unweighted - should the unweighted ROC curve be plotted (i.e. the normal ROC curve without class divisions)

#HELPER FUNCTION PARAMETERS
#lwd - vector of line widths. if NULL, the lines will all have a width of 1
#unweighted.lwd - width of the unweighted curve
#unweighted.col - color of the unweighted curve
#background.color - background color of the plot
#grid.marks - increment between grid lines
#grid.color - grid line color
#grid.lty - grid line type
#grid.lwd - grid line width
#diag.lty - diagonal line style
#legend.lty - legend style (0 is no box, 1 is boxed legend)
#legend.color - color of the legend box
#cex.main - title size
#cex.subtitle - subtitle size
#cex.legend - legend text size

#RETURN VALUE
#A ROC plot showing diagnostic performance across multiple designated classes

#REQUIRED PACKAGES: ROCR, MetaIntegrator

makeClassROCplot <- function(pooledDataObject, filterObject=NULL, upgenes=NULL, downgenes=NULL, caseNames=NULL, otherNames=NULL, group.colname="group", title=NULL, textSize=NULL,
                             full.labels=TRUE, rounding=3, colors=c(2,4,6,5,3,7:10), plot.unweighted=TRUE, ...){
  if(!checkManateeObject(pooledDataObject,"pooledDataObject")){stop("Invalid pooledDataObject provided")}
  if(is.null(filterObject)){
    if(is.null(upgenes) && is.null(downgenes)){stop("If filterObject is null, then upgenes and downgenes must be provided")}
    .makeClassROCplotHelper(prediction.type="geometric.smd", class=pooledDataObject$class, lasso.obj=NULL, geneMtx=pooledDataObject$genes, lambda.vals=NULL, pos.genes=upgenes,
                            neg.genes=downgenes, labelVec=pooledDataObject$pheno[,group.colname], caseNames=caseNames, otherNames=otherNames, title=title, textSize=textSize, full.labels = full.labels,
                            rounding=rounding, colors=colors, plot.unweighted=plot.unweighted, ...)
  }else{
    if(!checkManateeObject(filterObject,"filterObject")){stop("Invalid filterObject provided")}
    .makeClassROCplotHelper(prediction.type="geometric.smd", class=pooledDataObject$class, lasso.obj=NULL, geneMtx=pooledDataObject$genes, lambda.vals=NULL, pos.genes=filterObject$upGeneNames,
                            neg.genes=filterObject$downGeneNames, labelVec=pooledDataObject$pheno[,group.colname], caseNames=caseNames, otherNames=otherNames, title=title, textSize=textSize, full.labels = full.labels,
                            rounding=rounding, colors=colors, plot.unweighted=plot.unweighted, ...)
  }
}



###-###-###-###-###-###-###-##
###   makeClassPRCplot()   ###
###-###-###-###-###-###-###-##

#DESCRIPTION
#Wrapper function to create a Class PRC plot (a plot with multiple PRC curves
#comparing your positive class (i.e. your cases) against each of your other
#classes (i.e. your controls)) from a pooledDataObject using the genes
#from a filterObject (or multiple filterObjects) to make predictions with the
#geometric mean score model

#PARAMETERS
#pooledDataObject - A MANATEE \code{pooledDataObject}
#filterObject - A MANATEE \code{filterObject}
#upgenes - if desired, a character vector of the upgenes can be provided instead of a filterObject
#downgenes - if desired, a character vector of the downgenes can be provided instead of a filterObject
#caseNames - if left blank, this will be a character vector of the class(es) that correspond to your cases (i.e. samples labeled with a 1 in pooledDataObject$class).
#            If the user wants to exclude any of these case classes, then this can be set as the names of all of the case classes that should be included
#otherNames - if left blank, this will be a character vector of the class(es) that correspond to your controls (i.e. samples labeled with a 0 in pooledDataObject$class).
#             If the user wants to exclude any of these control classes, then this can be set as the names of all of the control classes that should be included
#group.colname - this designates the column name of the group vector within the pooledDataObject$pheno dataframe

#GRAPHICAL PARAMETERS
#title - title of the figure
#textSize - use this to easily increase or decrease the size of all the text in the plot
#full.labels - if TRUE, then the full names of the case categories will be used in the legend, if FALSE then a generic "Cases" will be used instead
#rounding - how many digits to round the AUC and CI to
#colors - vector of colors for the lines
#plot.unweighted - should the unweighted PRC curve be plotted (i.e. the normal PRC curve without class divisions)

#HELPER FUNCTION PARAMETERS
#lwd - vector of line widths. if NULL, the lines will all have a width of 1
#unweighted.lwd - width of the unweighted curve
#unweighted.col - color of the unweighted curve
#background.color - background color of the plot
#grid.marks - increment between grid lines
#grid.color - grid line color
#grid.lty - grid line type
#grid.lwd - grid line width
#legend.lty - legend style (0 is no box, 1 is boxed legend)
#legend.color - color of the legend box
#cex.main - title size
#cex.subtitle - subtitle size
#cex.legend - legend text size

#RETURN VALUE
#A PRC plot showing diagnostic performance across multiple designated classes

#REQUIRED PACKAGES: ROCR, MetaIntegrator, zoo, pracma

makeClassPRCplot <- function(pooledDataObject, filterObject=NULL, upgenes=NULL, downgenes=NULL, caseNames=NULL, otherNames=NULL, group.colname="group", title=NULL, textSize=NULL,
                             full.labels=TRUE, rounding=3, colors=c(2,4,6,5,3,7:10), plot.unweighted=TRUE, ...){
  if(!checkManateeObject(pooledDataObject,"pooledDataObject")){stop("Invalid pooledDataObject provided")}
  if(is.null(filterObject)){
    if(is.null(upgenes) && is.null(downgenes)){stop("If filterObject is null, then upgenes and downgenes must be provided")}
    .makeClassPRCplotHelper(prediction.type="geometric.smd", class=pooledDataObject$class, lasso.obj=NULL, geneMtx=pooledDataObject$genes, lambda.vals=NULL, pos.genes=upgenes,
                            neg.genes=downgenes, labelVec=pooledDataObject$pheno[,group.colname], caseNames=caseNames, otherNames=otherNames, title=title, textSize=textSize, full.labels = full.labels,
                            rounding=rounding, colors=colors, plot.unweighted=plot.unweighted, ...)
  }else{
    if(!checkManateeObject(filterObject,"filterObject")){stop("Invalid filterObject provided")}
    .makeClassPRCplotHelper(prediction.type="geometric.smd", class=pooledDataObject$class, lasso.obj=NULL, geneMtx=pooledDataObject$genes, lambda.vals=NULL, pos.genes=filterObject$upGeneNames,
                            neg.genes=filterObject$downGeneNames, labelVec=pooledDataObject$pheno[,group.colname], caseNames=caseNames, otherNames=otherNames, title=title, textSize=textSize, full.labels = full.labels,
                            rounding=rounding, colors=colors, plot.unweighted=plot.unweighted, ...)
  }
}



###-###-###-###-###-###-###
###   getNumStudies()   ###
###-###-###-###-###-###-###

#DESCRIPTION
#This function determines how many studies and what percentage of samples each
#of your genes was measured in. It also reports how many of those studies have
#at least one case (numStudies1) or at least one control (numStudies0), to allow
#for filtering of any genes that were only measured in cases or controls.

#PARAMETERS
#pooledDataObject - A MANATEE pooledDataObject Object
#geneNames - Specifies which genes to calculate numStudies statistics for. If left blank, then statistics will be calculated for all genes.

#RETURN VALUE
#returns a data frame that has four statistics for each gene:
#numStudies - The number of studies that the gene was measured in
#numStudies1 - The number of studies with at least one case sample that the gene was measured in
#numStudies0 - The number of studies with at least one control sample that the gene was measured in
#propSamples - The proportion of samples that the gene was measured in

#REQUIRED PACKAGES: NONE
getNumStudies <- function(pooledDataObject, geneNames=NULL){
  #add some error checking for geneNames at some point
  studyNames = pooledDataObject$pheno$name
  studyNames.unique = unique(studyNames)
  total = length(studyNames)
  studyClasses = lapply(studyNames.unique,function(name) pooledDataObject$class[studyNames %in% name])
  names(studyClasses) = studyNames.unique
  if(is.null(geneNames)){
    my.genes = pooledDataObject$genes
  }else{
    my.genes = pooledDataObject$genes[geneNames,]
  }

  numStudiesResults = do.call(rbind,apply(my.genes,1,function(gene){
    count.num = count = count1 = count0 = 0
    for(name in studyNames.unique){
      studyGenes = gene[studyNames %in% name]
      if(!any(is.na(studyGenes))){
        count = count+1
        count.num = count.num+length(studyGenes)
        if(any(studyClasses[[name]]==1)){
          count1 = count1+1
        }
        if(any(studyClasses[[name]]==0)){
          count0 = count0+1
        }
      }
    }
    return(data.frame(numStudies=count,numStudies1=count1,numStudies0=count0,propSamples=count.num/total))
  }))
  rm(my.genes)

  return(numStudiesResults)
}



###-###-###-###-###-###-###-##
###   mostRecentFilter()   ###
###-###-###-###-###-###-###-##

#DESCRIPTION
#Get the most recent filter

#Given a \code{pooledDataObject}, this function will look through
#\code{$filterResults} for the most recent filter used and will either return
#the most recently generated filterObject, or just the name of the most recently
#generated filterObject.

#PARAMETERS
#pooledDataObject - A MANATEE pooledDataObject Object
#returnName - if TRUE, then only the name of the most recent filter will be returned. If FALSE, then the actual filterObject will be returned.

#RETURN VALUE
#Either the most recently generated filterObject or the name of the most recent filter

#REQUIRED PACKAGES: NONE
mostRecentFilter <- function(pooledDataObject, returnName=FALSE){
  positions <- order(unlist(lapply(pooledDataObject$filterResults,
                                   function(i) i$filterDescription$timestamp)), decreasing = TRUE)
  recentName = names(pooledDataObject$filterResults)[positions[1]]
  if(returnName){
    return(recentName)
  }else{
    return(pooledDataObject$filterResults[[recentName]])
  }
}



###-###-###-###-###-###
###   runiLasso()   ###
###-###-###-###-###-###

#DESCRIPTION
#Unpolished explanation (modify later): Since lasso initially automatically
#removes genes, I start off by training lasso over and over again as it steadily
#removes genes (e.g. I start off with 100 genes, see that the first lasso has an
#“optimal” cutoff at a model with 90 genes so I retrain using only those 90
#genes, see that the next lasso has an optimal cutoff at 80 genes, etc) until it
#stops automatically removing genes. I then remove one gene at a time, choosing
#the best n-1 gene lasso model each time (e.g. I start with 50 genes, then
#choose the best 49 gene model, retrain lasso with those 49 genes, then choose
#the best 48 gene model, etc).

#PARAMETERS
#pooledDataObject - A MANATEE pooledDataObject. If weighting is desired, then sample weights must be stored in pooledDataObject$sample.weights
#filterObject - a MANATEE filterObject
#upgenes - character vector of the upregulated genes
#downgenes - character vector of the downregulated genes
#Valid.pooledDataObject - A MANATEE pooledDataObject, to be used for validation of each gene signature
#perf.meas - Which diagnostic performance metric to use. Options are "AUC", "AUPRC", and "AvgAUC" (Default: c("AUC","AUPRC"))
#group.colname - if AvgAUC is to be calculated, this designates the column name of the group vector within the pooledDataObject$pheno dataframe
#caseNames - if left blank, this will be a character vector of the class(es) that correspond to your cases (i.e. samples labeled with a 1 in pooledDataObject$class).
#            If the user wants to exclude any of these case classes, then this can be set as the names of all of the case classes that should be included
#otherNames - if left blank, this will be a character vector of the class(es) that correspond to your controls (i.e. samples labeled with a 0 in pooledDataObject$class).
#             If the user wants to exclude any of these control classes, then this can be set as the names of all of the control classes that should be included
#weighted - if TRUE, weighting is performed using the weights in pooledDataObject$sample.weights. If not specified, then weighting will automatically be performed if \code{pooledDataObject$sample.weights} is populated,
#           but weighting will not be performed if \code{pooledDataObject$sample.weights} is NULL
#nfolds - number of folds; default is 10. Although nfolds can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets. Smallest value allowable is nfolds=3
#alpha - The elasticnet mixing parameter, with 0 ≤ a ≤ 1. alpha=1 is the lasso penalty, and alpha=0 the ridge penalty
#nlambda - The number of lambda values; default is 300
#seed - random seed
#rounding - how many digits to round the performance metrics to
#use.lambda.1se - if TRUE, then lambda.1se will be used as the threshold when selecting new models during the first loop. If FALSE, then lambda.min will be used instead.
#... - other arguments for cv.glmnet()

#RETURN VALUE
#Returns a list containing the results of the analysis:
#   discoResults: the results within the discovery data
#   validResults: the results within the validation data
#   iLassoDescription: information about the parameters of the analysis

#REQUIRED PACKAGES: glmnet, MetaIntegrator, ROCR, zoo, pracma

#maybe add progress bar, or output certain AUCs, or something while it runs to give people a sense of where it's at
runiLasso <- function(pooledDataObject, filterObject=NULL, upgenes=NULL, downgenes=NULL, Valid.pooledDataObject=NULL,
                      perf.meas=c("AUC","AUPRC"), group.colname="group", caseNames=NULL, otherNames=NULL, weighted=NULL,
                      nfolds=10, alpha=1, nlambda=300, seed=1337, rounding=3, use.lambda.1se=TRUE, ...){
  call = sys.call(which = 0)
  #do default MANATEE checking stuff
  if(!checkManateeObject(pooledDataObject,"pooledDataObject")){stop("Invalid pooledDataObject provided")}
  if(!is.null(filterObject) && !checkManateeObject(filterObject,"filterObject")){stop("Invalid filterObject provided")}
  #check if perf.meas is valid (can't contain anything other than AUC/AUPRC/AvgAUC)
  perf.meas = tolower(perf.meas)
  if(any(!perf.meas %in% c("auc","auprc","avgauc"))){
    stop("An invalid perf.meas was provided. The only options are AUC, AUPRC, and AvgAUC (not case-sensitive)")
  }
  
  #check to make sure either filterObject is present, or upgenes/downgenes are present
  if(is.null(filterObject) && is.null(upgenes) && is.null(downgenes)){
    stop("Either a filterObject must be provided or upgenes/downgenes must be provided")
  }
  #if filterObject is NULL, make the upgenes/downgenes into a filterObject
  if(is.null(filterObject)){
    filterObject = list(upGeneNames = upgenes,downGeneNames = downgenes)
  }
  
  #
  #when AvgAUC is included, do the caseNames/otherNames stuff here
  #
  
  valid = TRUE
  if(is.null(Valid.pooledDataObject)){
    valid = FALSE
  }
  
  #First run a loop where I keep on retraining until lasso stops automatically removing genes
  iLassoResults = list()
  if(valid){iLassoResults.valid = list()}
  filter.curr = filterObject
  firstLoop = TRUE
  set.seed(seed)
  
  while(firstLoop){
    loop1.auc = loop1.auprc = loop1.avgauc = NULL
    if(valid){loop1V.auc = loop1V.auprc = loop1V.avgauc = NULL}
    scores = getGeneScores(pooledDataObject$genes,filter.curr$upGeneNames,filter.curr$downGeneNames)
    if(valid){scores.valid = getGeneScores(Valid.pooledDataObject$genes,filter.curr$upGeneNames,filter.curr$downGeneNames)}
    if("auc" %in% perf.meas){
      auc.all = getAUC(pooledDataObject$class,scores,rounding)
      loop1.auc = c(auc.all$auc,auc.all$auc.CI[1],auc.all$auc.CI[2])
      names(loop1.auc) = c("AUC","AUC.ci.lower","AUC.ci.upper")
      if(valid){
        auc.all.v = getAUC(Valid.pooledDataObject$class,scores.valid,rounding)
        loop1V.auc = c(auc.all.v$auc,auc.all.v$auc.CI[1],auc.all.v$auc.CI[2])
        names(loop1V.auc) = c("AUC","AUC.ci.lower","AUC.ci.upper")
      }
    }
    if("auprc" %in% perf.meas){
      auprc.all = getAUPRC(pooledDataObject$class,scores,rounding)
      loop1.auprc = c(auprc.all$auprc,auprc.all$auprc.CI[1],auprc.all$auprc.CI[2])
      names(loop1.auprc) = c("AUPRC","AUPRC.ci.lower","AUPRC.ci.upper")
      if(valid){
        auprc.all.v = getAUPRC(Valid.pooledDataObject$class,scores.valid,rounding)
        loop1V.auprc = c(auprc.all.v$auprc,auprc.all.v$auprc.CI[1],auprc.all.v$auprc.CI[2])
        names(loop1V.auprc) = c("AUPRC","AUPRC.ci.lower","AUPRC.ci.upper")
      }
    }
    if("avgauc" %in% perf.meas){
      #ADD THIS LATER
    }
    
    numGenes = length(filter.curr[[1]])+length(filter.curr[[2]])
    geneList = c(paste(filter.curr[[1]],collapse = " / "),paste(filter.curr[[2]],collapse = " / "))
    names(geneList) = c("upgenes","downgenes")
    loop1Results = c(numGenes,loop1.auc,loop1.auprc,loop1.avgauc,geneList)
    names(loop1Results)[1] = "numGenes"
    iLassoResults[[as.character(numGenes)]] = loop1Results
    if(valid){
      loop1Results.valid = c(numGenes,loop1V.auc,loop1V.auprc,loop1V.avgauc,geneList)
      names(loop1Results.valid)[1] = "numGenes"
      iLassoResults.valid[[as.character(numGenes)]] = loop1Results.valid
    }
    
    lasso.loop1 = runLasso.cv(pooledDataObject,filter.curr,weighted=weighted,nfolds=nfolds,alpha=alpha,nlambda=nlambda,plot=F)
    
    if(use.lambda.1se){
      next.numGenes = min(lasso.loop1$lassoResults$nzero[which(lasso.loop1$lassoResults$lambda == lasso.loop1$lassoResults$lambda.1se)])
    }else{
      next.numGenes = min(lasso.loop1$lassoResults$nzero[which(lasso.loop1$lassoResults$lambda == lasso.loop1$lassoResults$lambda.min)])
    }
    
    topGenes = getTopLassoGenes(lasso.loop1, next.numGenes, "topLassoGenes",plot = F)
    filter.curr = list(upGeneNames = topGenes$upgenes, downGeneNames = topGenes$downgenes)
    if(next.numGenes >= numGenes){firstLoop=FALSE}
  }
  cat(paste0("Lasso stopped automatically removing genes at ",numGenes, " genes; the remaining genes will be manually removed."))
  loop1.finalnumGenes = numGenes
  
  #next run a loop where I remove 1 gene (or however many is the fewest genes I can remove) at a time
  secondLoop = TRUE
  loop2.allNums = lasso.loop1$lassoResults$glmnet.fit$df
  topGenes = getTopLassoGenes(lasso.loop1, max(loop2.allNums[loop2.allNums < loop1.finalnumGenes]), "topLassoGenes",plot = F)
  filter.curr = list(upGeneNames = topGenes$upgenes, downGeneNames = topGenes$downgenes)
  #not worth it to spend time on this now, but eventually should include something
  #to deal with a situation where i've already reached the minimum number of genes
  #before going to the second loop
  
  while(secondLoop){
    loop2.auc = loop2.auprc = loop2.avgauc = NULL
    if(valid){loop2V.auc = loop2V.auprc = loop2V.avgauc = NULL}
    scores = getGeneScores(pooledDataObject$genes,filter.curr$upGeneNames,filter.curr$downGeneNames)
    if(valid){scores.valid = getGeneScores(Valid.pooledDataObject$genes,filter.curr$upGeneNames,filter.curr$downGeneNames)}
    if("auc" %in% perf.meas){
      auc.all = getAUC(pooledDataObject$class,scores,rounding)
      loop2.auc = c(auc.all$auc,auc.all$auc.CI[1],auc.all$auc.CI[2])
      names(loop2.auc) = c("AUC","AUC.ci.lower","AUC.ci.upper")
      if(valid){
        auc.all.v = getAUC(Valid.pooledDataObject$class,scores.valid,rounding)
        loop2V.auc = c(auc.all.v$auc,auc.all.v$auc.CI[1],auc.all.v$auc.CI[2])
        names(loop2V.auc) = c("AUC","AUC.ci.lower","AUC.ci.upper")
      }
    }
    if("auprc" %in% perf.meas){
      auprc.all = getAUPRC(pooledDataObject$class,scores,rounding)
      loop2.auprc = c(auprc.all$auprc,auprc.all$auprc.CI[1],auprc.all$auprc.CI[2])
      names(loop2.auprc) = c("AUPRC","AUPRC.ci.lower","AUPRC.ci.upper")
      if(valid){
        auprc.all.v = getAUPRC(Valid.pooledDataObject$class,scores.valid,rounding)
        loop2V.auprc = c(auprc.all.v$auprc,auprc.all.v$auprc.CI[1],auprc.all.v$auprc.CI[2])
        names(loop2V.auprc) = c("AUPRC","AUPRC.ci.lower","AUPRC.ci.upper")
      }
    }
    if("avgauc" %in% perf.meas){
      #ADD THIS LATER
    }
    
    numGenes = length(filter.curr[[1]])+length(filter.curr[[2]])
    geneList = c(paste(filter.curr[[1]],collapse = " / "),paste(filter.curr[[2]],collapse = " / "))
    names(geneList) = c("upgenes","downgenes")
    loop2Results = c(numGenes,loop2.auc,loop2.auprc,loop2.avgauc,geneList)
    names(loop2Results)[1] = "numGenes"
    iLassoResults[[as.character(numGenes)]] = loop2Results
    if(valid){
      loop2Results.valid = c(numGenes,loop2V.auc,loop2V.auprc,loop2V.avgauc,geneList)
      names(loop2Results.valid)[1] = "numGenes"
      iLassoResults.valid[[as.character(numGenes)]] = loop2Results.valid
    }
    
    lasso.loop2 = runLasso.cv(pooledDataObject,filter.curr,weighted=weighted,nfolds=nfolds,alpha=alpha,nlambda=nlambda,plot=F)
    
    loop2.allNums = lasso.loop2$lassoResults$nzero #replacing df here to see if I can avoid a bug
    loop2.allNums = loop2.allNums[loop2.allNums < numGenes]
    if(max(loop2.allNums) <= 1){
      secondLoop = FALSE
    }else{
      topGenes = getTopLassoGenes(lasso.loop2, max(loop2.allNums), "topLassoGenes",plot = F)
      filter.curr = list(upGeneNames = topGenes$upgenes, downGeneNames = topGenes$downgenes)
    }
  }
  
  #now need to collapse the results list
  iLassoResults = data.table(do.call(rbind,iLassoResults),stringsAsFactors = FALSE)
  iLassoResults.valid = data.table(do.call(rbind,iLassoResults.valid),stringsAsFactors = FALSE)
  
  iLassoDescription = list(perf.meas=perf.meas, loop1_finalnumGenes=loop1.finalnumGenes,weighted=weighted,nfolds=nfolds,
                           alpha=alpha,nlambda=nlambda,lambda_cutoff=ifelse(use.lambda.1se,"lambda.1se","lambda.min"),
                           group.colname=group.colname,caseNames=caseNames,otherNames=otherNames,seed=seed,call=call)
  
  iLassoResults = list(discoResults = iLassoResults, validResults = iLassoResults.valid,
                       iLassoDescription = iLassoDescription)
  return(iLassoResults)
}

#wrapper function with the old function name for compatibility w/ old scripts
runIteratingLasso <- function(pooledDataObject, filterObject=NULL, upgenes=NULL, downgenes=NULL, Valid.pooledDataObject=NULL,
                              perf.meas=c("AUC","AUPRC"), group.colname="group", caseNames=NULL, otherNames=NULL, weighted=NULL,
                              nfolds=10, alpha=1, nlambda=300, seed=1337, rounding=3, use.lambda.1se=TRUE, ...){
  return(runiLasso(pooledDataObject, filterObject, upgenes, downgenes, Valid.pooledDataObject,
                   perf.meas, group.colname, caseNames, otherNames, weighted,
                   nfolds, alpha, nlambda, seed, rounding, use.lambda.1se, ...))
}











































###-###-###-###-###-###-###-##-#
###   checkManateeObject()   ###
###-###-###-###-###-###-###-##-#

#DESCRIPTION
#This function checks whether an object that the user provides is formatted correctly

#PARAMETERS
#object - object to be checked by the function
#objectType - the object type that you want to check (options are: "Dataset", "DatasetList", "pooledDataObject")

#RETURN VALUE
#TRUE if passed error checking, FALSE otherwise

#REQUIRED PACKAGES: None

#maybe add one for lassoObject?

#need to find a good way to build up these errors and release them all at once (maybe even associated w/ each dataset if possible??) rather than just throwing an error every time I see something wrong

checkManateeObject <- function(object, objectType){
  if(objectType=="Dataset"){
    #needs formattedName
    if(is.null(object$formattedName)){
      stop("Datasets must have their formatted name contained in Dataset$formattedName")
    }
    datasetName = object$formattedName
    #check that rownames of $pheno match up with colnames of $genes
    if(!all(rownames(object$pheno) == colnames(object$genes))){
      stop(paste0("Error in Dataset ",datasetName,": The rownames of Dataset$pheno do not match the colnames of Dataset$genes"))
    }
    #make sure to add something that checks for missing data content of the datasets (since some of the downstream functions can't deal with missing data)
    #needs group in $pheno
    #needs control.0.class in $pheno, and $pheno$control.0.class must have both 0s and 1s
    if(!"control.0.class" %in% colnames(object$pheno)){
      stop(paste0("Error in Dataset ",datasetName,": You must have a column indicating which samples are healthy (0) and which are non-healthy (1), stored in Dataset$pheno$control.0.class",
                  " - this column will be automatically added if you run the formatManateeData() function and keep the AddControl0Class option set to TRUE."))
    }
    if(!(any(object$pheno$control.0.class == 0) && any(object$pheno$control.0.class == 1)) || any(!object$pheno$control.0.class %in% c(0,1))){
      stop(paste0("Error in Dataset ",datasetName,": Dataset$pheno$control.0.class must contain both 0s (to indicate healthy samples) and 1s (to indicate non-healthy samples)",
                  " and cannot contain any other digits."))
    }
    #needs $genes
    #class has to be a factor IS THIS STILL TRUE??
    #class has to be 0 and 1 for now (output message saying that multiclass will be added later)
    return(TRUE) #dummy function right now, will update later
  }
  if(objectType=="DatasetList"){
    #first check that the DatasetList is a named list
    if(class(object) != "list"){
      stop("DatasetList must be a list")
    }
    if(is.null(names(object))){
      stop("DatasetList must be named")
    }
    #and then just do an lapply with the check dataset function
    TEMP = lapply(object, function(x){checkManateeObject(x,"Dataset")})
    #TEMPORARY - running a for loop instead to try and preserve the error messages
    # for(i in 1:length(object)){
    #   checkManateeObject(object[[i]],"Dataset")
    # }
    return(TRUE) #dummy function right now, will update later
  }
  if(objectType=="pooledDataObject"){
    return(TRUE) #dummy function right now, will update later
    #should be fairly similar to Dataset
    #and maybe doesn't need some of the things Dataset needs
    #check for weights, if it's missing give a warning that says "option weights vector not present"
  }
  if(objectType=="filterObject"){
    return(TRUE) #dummy function right now, will update later
  }
  if(objectType=="manateeResults"){
    return(TRUE) #dummy function right now, will update later
    #make sure that both upgenes and downgenes aren't empty
    #check colnames (use an %in%) - give warning if the colnames aren't exactly correct but don't return false
      #"Some of the column names in manateeResults don't match the expected names. This is allowed, but has the potential to cause errors."
  }
  if(objectType=="bssResults"){
    return(TRUE) #dummy function right now, will update later
  }
}



