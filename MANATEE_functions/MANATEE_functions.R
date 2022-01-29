# MANATEE functions
# Author: Aditya Rao
# Date: 1/25/18
# Contact: adityamr@stanford.edu

# Functions for the MANATEE framework

#OVERARCHING NOTES
#need to standardize upgenes/posgenes/pos.genes/etc (and just lots of naming things tbh)
#use the getAUC and getAUPRC functions to replace a bunch of unnecessary lines of code!
#add ... notation to a bunch more functions
#make it so that more functions can use either a filterObject or a list of upgenes/downgenes



#............................................................................................................#
##### Functions for Processing Data #####
#............................................................................................................#

###-###-###-###-###-###-###-###-#
###   excelDateCorrection()   ###
###-###-###-###-###-###-###-###-#

#DESCRIPTION
#This function is used to correct an issue that comes up in Excel where the
#names of some genes (e.g. MARCH1, DEC1) are automatically converted to dates
#(e.g. 1-Mar, 1-Dec). If all you have is your list of gene symbols, then it will
#fix all of the genes that could only correspond to one date. However, there are
#some genes that will result in the same date (e.g. MARCH1 and MARC1 will both
#produce 1-Mar) so if you want to fix those, you must also provide a list of
#gene IDs that corresponds to the gene symbols. The IDs that this currently
#supports are listed below.

#NOTE: This has been updated to also change some gene symbols that are commonly
#found under the wrong name (e.g. converting "Septin #" genes to "SEPT#")

#NOTE: This also now uses the HGNChelper package to do some additional gene symbol adjustments

#PARAMETERS
#geneNames - list of gene symbols
#geneIDs - list of gene IDs
#HGNCgeneCorrection - If TRUE, then the HGNChelper package will be used to fix all remaining gene names
#print - if TRUE, messages will be printed to the console

#RETURN VALUE
#A modified version of the input geneNames that has been corrected to remove dates

#REQUIRED PACKAGES: None

#This is the list of identifiers that my code currently checks for:
# 1) HGNC
# 2) Entrez Gene
# 3) Ensembl
# 4) OMIM
# 5) UniProtKB
# 6) UniGene Cluster
# 7) UniGene Representative Sequence

excelDateCorrection <- function(geneNames, geneIDs=NULL, HGNCgeneCorrection = TRUE, print=TRUE){
  
  if(HGNCgeneCorrection){
    mapping = suppressMessages(suppressWarnings(HGNChelper::checkGeneSymbols(geneNames,unmapped.as.na=FALSE)))
    #mapping$Suggested.Symbol[mapping$Suggested.Symbol == "NA"] = NA
    mapping$Suggested.Symbol[is.na(mapping$Suggested.Symbol)] = geneNames[is.na(mapping$Suggested.Symbol)]
    geneNames[1:length(geneNames)] = mapping$Suggested.Symbol
    geneNames = gsub(","," /// ",geneNames)
  }
  
  uniqueMonthNames = c("1-Feb","2-Feb","4-Feb","5-Feb","6-Feb","7-Feb","9-Feb","10-Feb","3-Mar","4-Mar","5-Mar",
                       "6-Mar","7-Mar","8-Mar","9-Mar","10-Mar","11-Mar","1-Sep","2-Sep","3-Sep","4-Sep","5-Sep",
                       "6-Sep","7-Sep","8-Sep","9-Sep","10-Sep","11-Sep","12-Sep","14-Sep","15-Sep","1-Dec",
                       "Septin 1","Septin 2","Septin 3","Septin 4","Septin 5","Septin 6","Septin 7","Septin 8",
                       "Septin 9","Septin 10","Septin 11","Septin 12","Septin 13","Septin 14","Selenoprotein 15",
                       "Gcom1")

  geneSymbols = c("FEB1","FEB2","FEB4","FEB5","FEB6","FEB7","FEB9","FEB10","MARCH3","MARCH4","MARCH5",
                  "MARCH6","MARCH7","MARCH8","MARCH9","MARCH10","MARCH11","SEPT1","SEPT2","SEPT3","SEPT4","SEPT5",
                  "SEPT6","SEPT7","SEPT8","SEPT9","SEPT10","SEPT11","SEPT12","SEPT14","SEP15","DEC1",
                  "SEPT1","SEPT2","SEPT3","SEPT4","SEPT5","SEPT6","SEPT7","SEPT8",
                  "SEPT9","SEPT10","SEPT11","SEPT12","SEPT13","SEPT14","SEP15",
                  "GCOM1")

  badData = c("1-Mar","2-Mar",uniqueMonthNames) %in% geneNames
  if(print && any(badData==TRUE)){
    cat("These dates are currently in your gene names:",c("1-Mar","2-Mar",uniqueMonthNames)[badData],"\n")
  } else{
    if(print){cat("There are no dates in your gene names :)")}
    return(geneNames)
  }

  for(i in 1:length(uniqueMonthNames)){
    geneNames[geneNames==uniqueMonthNames[i]] = geneSymbols[i]
  }

  if(!is.null(geneIDs)){
    if(any(geneNames=="1-Mar")){
      iter = which(geneNames=="1-Mar", na.rm = T)
      for(i in 1:length(iter)){
        if(geneIDs[iter[i]] %in% c("26077","55016","ENSG00000145416","613331","Q8TCQ1","Hs.592804","NM_001166373")){
          geneNames[iter[i]] = "MARCH1"
        } else if(geneIDs[iter[i]] %in% c("26189","64757","ENSG00000186205","614126","Q5VT66","Hs.497816","AK092439")){
          geneNames[iter[i]] = "MARC1"
        }
      }
    }
    if(any(geneNames=="2-Mar")){
      iter = which(geneNames=="2-Mar", na.rm = T)
      for(i in 1:length(iter)){
        if(geneIDs[iter[i]] %in% c("28038","51257","ENSG00000099785","613332","Q9P0N8","Hs.631861","BC032624")){
          geneNames[iter[i]] = "MARCH2"
        } else if(geneIDs[iter[i]] %in% c("26064","54996","ENSG00000117791","614127","Q969Z3","Hs.369042","AK125512")){
          geneNames[iter[i]] = "MARC2"
        }
      }
    }
  }

  if(print){
    badData2 = c("1-Mar","2-Mar",uniqueMonthNames) %in% geneNames
    if(any(badData2==TRUE)){
      cat("These dates are still in your gene names:",c("1-Mar","2-Mar",uniqueMonthNames)[badData2])
    } else{
      cat("All dates were removed from your gene names.")
    }
  }

  return(geneNames)
}



###-###-###-###-###-###-##-#
###   .NAreplaceVals()   ###
###-###-###-###-###-###-##-#

#DESCRIPTION
#This function takes a data frame or matrix and, for each column, replaces all
#of the NAs in that dataframe with either ones or zeroes, based on the pattern
#strings that are provided.

#PARAMETERS
#data - the data frame or matrix that you want to have adjusted
#zero.string - a string with the pattern that identifies the column names for the columns that will have their NAs replaced by zeroes.
#              This will be used with grep, so regular expressions can be used.
#one.string - a string with the pattern that identifies the column names for the columns that will have their NAs replaced by ones.
#             This will be used with grep, so regular expressions can be used.

#RETURN VALUE
#A modified version of the input dataframe with the NAs replaced

#REQUIRED PACKAGES: None

.NAreplaceVals <- function(data, zero.string, one.string){
  zerocols = grep(zero.string,colnames(data))
  onecols = grep(one.string,colnames(data))
  for(col in zerocols){
    data[,col][is.na(data[,col])] = 0
  }
  for(col in onecols){
    data[,col][is.na(data[,col])] = 1
  }
  return(data)
}



###-###-###-###-###-###-###-###-###
###   combineDuplicateGenes()   ###
###-###-###-###-###-###-###-###-###

#DESCRIPTION
#This function takes a gene matrix and combines all of the rows that have
#duplicate genes by taking the average of all the duplicate genes for each
#sample. Any NAs are first replaced with the average of the other duplicates of
#that gene for that sample.

#PARAMETERS
#geneMtx = matrix of gene expression values (genes in rows, samples in columns)

#RETURN VALUE
#A modified version of the input geneMtx that has had its duplicate genes combined

#REQUIRED PACKAGES: None

combineDuplicateGenes <- function(geneMtx){
  if(length(rownames(geneMtx))==length(unique(rownames(geneMtx)))){
    return(geneMtx)
  }
  dupvals = unique(rownames(geneMtx)[duplicated(rownames(geneMtx))])
  duplicateindex = rownames(geneMtx) %in% dupvals
  duplicategenes = geneMtx[duplicateindex,,drop=F]
  uniquegenes = geneMtx[!duplicateindex,,drop=F]

  newvals = do.call(rbind, lapply(dupvals, function(val){
    dupindex = which(rownames(duplicategenes)==val)
    dupmat = rbind(duplicategenes[dupindex,,drop=F])
    #replace NA values with the average of the other values in that column
    if(any(is.na(dupmat))){
      dupmat = apply(dupmat,2,function(x){
        if(any(is.na(x))){
          avg = mean(x[!is.na(x)])
          x[is.na(x)]=avg
        }
        return(x)
      })
    }
    #average all rows
    return(colMeans(dupmat))
  }))

  rownames(newvals)=dupvals
  if(any(is.na(rownames(newvals)))){ #removing the row named "NA"
    NArow = which(is.na(rownames(newvals)))
    newvals = newvals[-NArow,]
  }
  return(rbind(uniquegenes,newvals))
}



###-###-###-###-###-###-###-##
###   checkSparseGenes()   ###
###-###-###-###-###-###-###-##

#DESCRIPTION
#This function takes a gene matrix that has missing values and checks how many
#genes would be removed by the removeSparseGenes() function at a variety of
#cutoffs

#PARAMETERS
#geneMtx - matrix of gene expression values (genes in rows, samples in columns)
#cutoffs - vector of percents of expression values that must be missing for a gene to be removed
#returnGenes - whether or not to display the genes that missed the cutoff

#RETURN VALUE
#Prints a message indicating how many genes would be removed at the provided cutoffs

#REQUIRED PACKAGES: None

checkSparseGenes <- function(geneMtx, cutoffs=c(0.05,0.1,0.15,0.2,0.25,0.3),returnGenes=FALSE){
  total = length(geneMtx[1,])
  generemove = lapply(cutoffs,function(cut){
    NAcutoff = total*cut
    num.remove=0
    geneList=""
    for(i in 1:length(geneMtx[,1])){
      if(sum(is.na(geneMtx[i,]))>NAcutoff){
        num.remove=num.remove+1
        geneList=paste0(geneList,rownames(geneMtx)[i],", ")
      }
    }
    return(list(num.remove,geneList))
  })

  for(i in 1:length(cutoffs)){
    cat(sprintf("At a cutoff of %s, %s genes would be removed\n",cutoffs[i],generemove[[i]][1]))
    if(returnGenes){
      cat(paste0("Gene List: ",generemove[[i]][2],"\n"))
    }
  }
}



###-###-###-###-###-###-###-###
###   removeSparseGenes()   ###
###-###-###-###-###-###-###-###

#DESCRIPTION
#This function takes a gene matrix that has missing values and removes any genes
#that have more missing data than your cutoff

#PARAMETERS
#geneMtx - matrix of gene expression values (genes in rows, samples in columns)
#cutoff - percent of expression values that must be missing for a gene to be removed

#RETURN VALUE
#A modified version of the input geneMtx with any genes below the cutoff removed

#REQUIRED PACKAGES: None

removeSparseGenes <- function(geneMtx, cutoff=0.15){
  total = length(geneMtx[1,])
  NAcutoff = total*cutoff
  newGeneMtx = do.call(rbind,apply(geneMtx,1,function(gene){
    num.NAs = sum(is.na(gene))
    if(num.NAs <= NAcutoff){
      return(gene)
    }
  }))
  cat(paste(length(geneMtx[,1])-length(newGeneMtx[,1]), "genes removed."))
  return(newGeneMtx)
}



###-###-###-###-###-###-###
###   removeSamples()   ###
###-###-###-###-###-###-###

#DESCRIPTION
#This function removes a set of samples from a dataset list by removing it from
#pheno, genes, expr, etc.

#PARAMETERS
#dataset - dataset to be modified (formatted as a list)
#samples - vector of index numbers of the samples you want to remove
#pheno - whether to remove the samples from $pheno
#genes - whether to remove the samples from $genes
#expr - whether to remove the samples from $expr
#class - whether to remove the samples from $class
#other.elements - vector of other elements to remove (currently only works with one-dimensional structures)

#RETURN VALUE
#A modified version of the input dataset with the designated samples removed
#from the indicated fields

#REQUIRED PACKAGES: None

removeSamples <- function(dataset,samples,pheno=TRUE,genes=TRUE,expr=TRUE,class=TRUE,other.elements=NULL){
  if(length(samples) == 0){
    warning("No samples were removed")
    return(dataset)
  }

  checkNames = c("pheno","genes","expr","class")[c(pheno,genes,expr,class)] %in% names(dataset)
  if(any(checkNames==FALSE)){
    missing = c("pheno","genes","expr","class")[c(pheno,genes,expr,class)][!checkNames]
    if(length(missing)==1){err=paste0("Missing ",missing,".")}
    else if(length(missing)==2){err=paste0("Missing ",missing[1]," and ",missing[2],".")}
    else{err=paste0("Missing ",paste(missing[-c(length(missing))],collapse = ", "),", and ",missing[length(missing)],".")}
    stop(err)
  }

  if(pheno){dataset$pheno = dataset$pheno[-samples,,drop=F]}
  if(genes){dataset$genes = dataset$genes[,-samples,drop=F]}
  if(expr){dataset$expr = dataset$expr[,-samples,drop=F]}
  if(class){
    dataset$class = as.numeric(dataset$class)
    dataset$class = dataset$class[-c(samples)]
  }
  if(!is.null(other.elements)){
    for(element in other.elements){
      if(!(element %in% names(dataset))){stop(sprintf("%s is not present in this dataset",element))}
      dataset[[element]] = dataset[[element]][-c(samples)]
    }
  }

  return(dataset)
}



###-###-###-###-###-###-##-#
###   makeMetaObject()   ###
###-###-###-###-###-###-##-#

#DESCRIPTION
#This function takes a list of datasets and converts it into a MetaIntegrator MetaObject

#PARAMETERS
#datasetlist - named list of datasets
#add - if true, only designated samples will be included
#remove - if true, designated samples will be removed
#select.col - if add or remove is true, this indicates the name of the row in $pheno that will be used to identify samples
#select.names - if add or remove is true, this is a vector of labels indicating which samples will be added/removed

#RETURN VALUE
#A modified version of the input datasetlist that is formatted as valid MetaIntegrator MetaObject

#REQUIRED PACKAGES: MetaIntegrator

makeMetaObject <- function(datasetlist,add = FALSE,remove = FALSE, select.col = NULL, select.names = NULL){
  require(MetaIntegrator)
  if(add && remove){
    stop("Cannot add and remove concurrently (if it's necessary, this functionality can be added)")
  }
  if((add || remove) && (is.null(select.col) || is.null(select.names))){
    stop("Either select.col or select.names is missing/null")
  }
  datasetlist = lapply(datasetlist,function(dataset){
    if(!("expr" %in% names(dataset))){
      dataset$expr = data.matrix(dataset$genes)
    }
    if(!("exp_comment" %in% names(dataset))){
      dataset$exp_comment = "placeholder exp_comment"
    }
    if(!("keys" %in% names(dataset))){
      dataset$keys = rownames(dataset$expr)
      names(dataset$keys) = dataset$keys
    }
    if(!("key_comment" %in% names(dataset))){
      dataset$key_comment = "Typical annotation"
    }
    if(!("formattedName" %in% names(dataset))){
      dataset$formattedName = dataset$pheno$name[1]
    }

    dataset$class=as.numeric(as.character(dataset$class))

    if(add){
      addsamples = c()
      select.colnum = which(colnames(dataset$pheno)==select.col)
      for(name in select.names){
        addsamples = append(addsamples,which(dataset$pheno[,select.colnum] == name))
      }
      if(length(addsamples)>0){
        dataset$pheno=dataset$pheno[addsamples,,drop=F]
        dataset$expr=dataset$expr[,addsamples,drop=F]
        if(("genes" %in% names(dataset))){
          dataset$genes=dataset$genes[,addsamples,drop=F]
        }
        dataset$class=dataset$class[addsamples]
      }
    }else if(remove){
      removesamples = c()
      select.colnum = which(colnames(dataset$pheno)==select.col)
      for(name in select.names){
        removesamples = append(removesamples,which(dataset$pheno[,select.colnum] == name))
      }
      if(length(removesamples)>0){
        dataset$pheno=dataset$pheno[-removesamples,,drop=F]
        dataset$expr=dataset$expr[,-removesamples,drop=F]
        if(("genes" %in% names(dataset))){
          dataset$genes=dataset$genes[,-removesamples,drop=F]
        }
        dataset$class=dataset$class[-removesamples,drop=F]
      }
    }

    names(dataset$class)=rownames(dataset$pheno)

    #run make.names on the colnames of expr, the rownames of pheno, and the names of class
    colnames(dataset$expr) = make.names(colnames(dataset$expr))
    rownames(dataset$pheno) = make.names(rownames(dataset$pheno))
    names(dataset$class) = make.names(names(dataset$class))

    if(!checkDataObject(dataset, "Dataset")){
      warning(paste0(dataset$formattedName," is not a valid MetaIntegrator Dataset Object"))
    }
    return(dataset)
  })

  datasetlist.out = list()
  datasetlist.out$originalData = datasetlist
  if(!checkDataObject(datasetlist.out, "Meta", "Pre-Analysis")){
    warning("The output object is not a valid MetaIntegrator Dataset Object")
  }
  return(datasetlist.out)
}



###-###-###-###-###-###-###-#
###   makeClassVector()   ###
###-###-###-###-###-###-###-#

#DESCRIPTION
#This function takes a vector of labels that indicate which group each sample
#belongs to and then creates a vector where every sample corresponding to the
#class that you've designated as your cases is represented by a 1 and every
#other sample is a 0

#PARAMETERS
#labels - either a vector of labels or a "pheno" dataframe with sample names as the rownames and the labels in $group that lists the disease group associated with each sample
#caseNames - name of the class(es) that you're considering to be your case (e.g. "Malaria")
#sampleNames - if labels is a vector, then this is the corresponding vector of sample names

#RETURN VALUE
#A vector of 1s and 0s that map to the provided labels, where the cases are 1 and the controls are 0

#REQUIRED PACKAGES: None

makeClassVector <- function(labels,caseNames,sampleNames = NULL){
  if(class(labels)[1]=="data.frame"){
    if(is.null(sampleNames)){sampleNames=rownames(labels)}
    labels = labels$group
  }
  class = rep(NA,length(labels))
  class[!(as.character(labels) %in% caseNames)] = 0
  class[as.character(labels) %in% caseNames] = 1
  class = as.numeric(class)
  if(!is.null(sampleNames)){
    stopifnot(identical(length(labels),length(sampleNames)))
    names(class) = sampleNames
  }
  return(class)
}



###-###-###-###-###-###-###-###
###   makeControl0Class()   ###
###-###-###-###-###-###-###-###

#DESCRIPTION
#This function takes a vector of labels that indicate which group each sample
#belongs to and then creates a class vector for use with COCONUT, where every
#sample corresponding to the class that you've designated as your control is
#represented by a 0 and every other sample is a 1

#PARAMETERS
#labels - either a vector of labels or a "pheno" dataframe with sample names as the rownames and the labels in $group that lists the disease group associated with each sample
#caseNames - name of the class(es) that you're considering to be your control (e.g. "Healthy")
#sampleNames - if labels is a vector, then this is the corresponding vector of sample names

#RETURN VALUE
#A vector of 1s and 0s that map to the provided labels, where the controls are 1 and the cases are 0

#REQUIRED PACKAGES: None

makeControl0Class <- function(labels,controlNames,sampleNames = NULL){
  if(class(labels)[1]=="data.frame"){
    if(is.null(sampleNames)){sampleNames=rownames(labels)}
    labels = labels$group
  }
  class = as.character(labels)
  class[!(class %in% controlNames)] = 1
  class[class %in% controlNames] = 0
  class = as.numeric(class)
  if(!is.null(sampleNames)){
    stopifnot(identical(length(labels),length(sampleNames)))
    names(class) = sampleNames
  }
  return(class)
}



###-###-###-###-###-###
###   splitData()   ###
###-###-###-###-###-###

#DESCRIPTION
#This function splits data into discovery and validation before COCONUT conormalization.

#PARAMETERS
#GSElist - list of datasets. It needs to have $pheno (matrix with samples in rows, sample info in columns) and $genes (matrix with genes in rows, samples in columns).
#           You need to have $name (vector of dataset names) in $pheno. Also, you need to have $class (numeric vector with 0=healthy,1=cases) in each dataset to be able to run COCONUT afterwards.
#min.healthy - The minimum number of healthy samples needed for a dataset to be split. If there are fewer healthies, then all of the dataset will go in discovery.
#prop.discovery - The proportion of data to put in discovery. Default is 70%
#healthy.names - if your "healthy" samples are not labeled as "Healthy" in the group column of $pheno, you can indicate the alternative "Healthy" name(s) here
#forced.discovery - Any datasets that you want to put only in discovery (using the dataset name stored in the $pheno$name vector)
#forced.validation - Any datasets that you want to put only in validation (using the dataset name stored in the $pheno$name vector)
#seed - the random seed to use

#RETURN VALUE
#A list containing the newly split Discovery data in one field and the Validation data in the other field

#REQUIRED PACKAGES: None

splitData <- function(GSElist, min.healthy=10, prop.discovery=0.7, healthy.names = "Healthy", forced.discovery=NULL, forced.validation=NULL, seed=1337){
  set.seed(seed)
  combo.datasets = lapply(GSElist, function(gse){
    if(gse$pheno$name[1] %in% forced.discovery){
      return(list(disco=gse,valid=NULL))
    } else if(gse$pheno$name[1] %in% forced.validation){
      return(list(disco=NULL,valid=gse))
    } else{
      ishealthy = gse$pheno$group %in% healthy.names
      if(sum(ishealthy)<min.healthy || sum(!ishealthy) <= 1){
        return(list(disco=gse,valid=NULL))
      }
      numDisc.h = round(sum(ishealthy)*prop.discovery, digits = 0)
      numDisc.d = round(sum(!ishealthy)*prop.discovery, digits = 0)
      hindices = sort(sample(c(1:nrow(gse$pheno))[ishealthy], numDisc.h, replace=FALSE))
      dindices = sort(sample(c(1:nrow(gse$pheno))[!ishealthy], numDisc.d, replace=FALSE))
      indices=append(hindices,dindices)
      gse.disco=gse
      gse.valid=gse
      gse.disco$pheno = gse.disco$pheno[indices,,drop=F]
      gse.disco$genes = gse.disco$genes[,indices,drop=F]
      gse.valid$pheno = gse.valid$pheno[-indices,,drop=F]
      gse.valid$genes = gse.valid$genes[,-indices,drop=F]
      #maybe change this to removeSamples at some point?
      return(list(disco=gse.disco,valid=gse.valid))
    }
  })

  disco.datasets = lapply(combo.datasets, function(gse){
    return(gse$disco)
  })

  valid.datasets = lapply(combo.datasets, function(gse){
    return(gse$valid)
  })

  for(i in length(GSElist):1){
    if(is.null(disco.datasets[[i]])){
      disco.datasets[[i]]=NULL
    }
    if(is.null(valid.datasets[[i]])){
      valid.datasets[[i]]=NULL
    }
  }

  return(list(Discovery=disco.datasets,Validation=valid.datasets))
}



###-###-###-###-###-###-###-###
###   bootstrapCocoData()   ###
###-###-###-###-###-###-###-###

#DESCRIPTION
#This function randomly selects a fraction of your COCONUT normalized samples, with replacement.

#PARAMETERS
#data - list of COCONUT normalized data. It needs to have $pheno (matrix with samples in rows, sample info in columns) and $genes (matrix with genes in rows, samples in columns).
#seed - the random seed to use
#case.separate - If true, the cases and controls are bootstrapped separately. Use this when you have a small number of cases to ensure that you select enough cases.
#case.rows - Vector of rows (in $pheno) representing which samples are cases
#prop.choose - The proportion of data to choose each time (default is 1, which is a typical bootstrap)
#forced.inclusion - Vector of sample names that you want to always include at least once (based on rowname in $pheno/colname in $genes)
#forced.exclusion - Vector of sample names that you want to always exclude (based on rowname in $pheno/colname in $genes)

#RETURN VALUE
#A modified version of the input data that now contains only the bootstrapped samples

#REQUIRED PACKAGES: None

bootstrapCocoData <- function(data, seed=1337, case.separate=TRUE, case.rows=NULL, prop.choose=1, forced.inclusion=NULL, forced.exclusion=NULL){
  set.seed(seed)

  if(any(forced.exclusion %in% forced.inclusion)){
    stop("Some samples are listed in both forced.inclusion and forced.exclusion.")
  }
  if(!is.null(forced.exclusion)){
    exclude.samples = !(rownames(data$pheno) %in% forced.exclusion)
    data$pheno = data$pheno[exclude.samples,,drop=F]
    data$genes = data$genes[,exclude.samples,drop=F]
  }
  if(!is.null(forced.inclusion)){
    include.samples = (rownames(data$pheno) %in% forced.inclusion)
    forced.pheno = data$pheno[include.samples,,drop=F]
    forced.genes = data$genes[,include.samples,drop=F]
  }

  if(case.separate){
    if(is.null(case.rows)){
      stop("No case.rows provided.")
    }
    num.choose.case = round(length(case.rows)*prop.choose)
    num.choose.control = round((dim(data$pheno)[1]-length(case.rows))*prop.choose) - length(forced.inclusion)
    choose.case = sample(case.rows,num.choose.case,replace=TRUE)
    choose.control = sample(c(1:dim(data$pheno)[1])[-case.rows],num.choose.control,replace=TRUE)
    pheno.case = data$pheno[choose.case,,drop=F]
    genes.case = data$genes[,choose.case,drop=F]
    pheno.control = data$pheno[choose.control,,drop=F]
    genes.control = data$genes[,choose.control,drop=F]
    pheno=rbind(pheno.case,pheno.control)
    genes=cbind(genes.case,genes.control)
  }else{
    num.choose = round(dim(data$pheno)[1]*prop.choose) - length(forced.inclusion)
    choose = sample(c(1:dim(data$pheno)[1]),num.choose,replace=TRUE)
    pheno = data$pheno[choose,,drop=F]
    genes = data$genes[,choose,drop=F]
  }

  if(!is.null(forced.inclusion)){
    pheno = rbind(forced.pheno,pheno)
    genes = cbind(forced.genes,genes)
  }

  return(list(pheno=pheno,genes=genes))
}




#............................................................................................................#
##### Functions for Analyzing Data #####
#............................................................................................................#

###-###-###-###-###-###-###
###   getGeneScores()   ###
###-###-###-###-###-###-###

#DESCRIPTION
#Given a matrix of genes and a list of the positive and negative genes in your
#signature, this function calculates the signature scores for each sample in
#that gene matrix.

#PARAMETERS
#geneMtx - matrix of gene expression values (genes in rows, samples in columns)
#pos - vector of the positive genes in your signature
#neg - vector of the negative genes in your signature
#makePos - if you are getting an error due to negative values in your matrix, set this as true
#out.missing - if true, the names of the missing genes will be output

#RETURN VALUE
#A vector of signature scores for each sample

#REQUIRED PACKAGES: None

getGeneScores <- function(geneMtx, pos, neg, makePos = TRUE, out.missing=TRUE){
  if(out.missing){
    missingpos = pos[!(pos %in% rownames(geneMtx))]
    missingneg = neg[!(neg %in% rownames(geneMtx))]
    missing = c(missingpos,missingneg)
    if(length(missing)>0){
      cat("Missing these genes:",missing,"\n")
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

  if(any(is.nan(scores))){ #this means all upgenes or all downgenes were missing from some samples
    nanIndex = which(is.nan(scores))
    for(i in nanIndex){
      my.score = c(geomMean(geneMtx[pos,i,drop=F],na.rm=T),geomMean(geneMtx[neg,i,drop=F],na.rm=T))
      my.score = my.score[!is.nan(my.score)]
      if(length(my.score) == 1){
        scores[i] = my.score
      }else if(length(my.score) == 2){
        scores[i] = my.score[1] - my.score[2]
      }
    }
  }

  return(scores)
}



##-###-###-###-###-#
###   getAUC()   ###
##-###-###-###-###-#

#DESCRIPTION
#Given a class vector and a vector of signature scores, this function calculates
#the AUC and 95% confidence interval. NOTE: This is metaintegrator, not ROCR

#PARAMETERS
#class - vector of classifications for your samples, with 0 representing your negative class (controls) and 1 representing your positive class (cases)
#scores - vector of signature score values for your samples
#rounding - how many digits to round the AUC and CI to

#RETURN VALUE
#A list with the AUC in one entry, and the 95% confidence interval in the other entry

#REQUIRED PACKAGES: MetaIntegrator

getAUC <- function(class, scores, rounding=3){
  pred_ROC = prediction(scores, class)
  auc = performance(pred_ROC, "auc")@y.values[[1]]
  auc = round(auc, rounding)
  CI = .getAUC_CI(auc,sum(class==1),sum(class==0))
  lo = max(0,round(CI[1], rounding))
  up = min(1,round(CI[2], rounding))
  
  # MI_ROC = calculateROC(as.numeric(as.character(class)), as.numeric(scores))
  # auc = round(MI_ROC$auc, rounding)
  # lo = max(0,round(MI_ROC$auc.CI[1], rounding))
  # up = min(1,round(MI_ROC$auc.CI[2], rounding))
  
  return(list(auc=auc,auc.CI=c(lo,up)))
}



###-###-###-###-###-##
###   getAUPRC()   ###
###-###-###-###-###-##

#DESCRIPTION
#Given a class vector and a vector of signature scores, this function calculates
#the AUPRC and 95% confidence interval

#PARAMETERS
#class - vector of classifications for your samples, with 0 representing your negative class (controls) and 1 representing your positive class (cases)
#scores - vector of signature score values for your samples
#rounding - how many digits to round the AUPRC and CI to

#RETURN VALUE
#A list with the AUPRC in one entry, and the 95% confidence interval in the other entry

#REQUIRED PACKAGES: ROCR, zoo, pracma

getAUPRC <- function(class, scores, rounding=3){
  pred_PRC = prediction(as.numeric(scores), matrix(class, nrow = length(class), ncol = 1))
  perf = performance(pred_PRC, "prec", "rec")
  if(is.nan(perf@y.values[[1]][1])){
    perf@y.values[[1]][1] = 1
  }
  perf@x.values[[1]] = c(0,perf@x.values[[1]],1)
  perf@y.values[[1]] = c(1,perf@y.values[[1]],0)
  auprc_all = .calcauprc(perf@x.values[[1]],perf@y.values[[1]],class)
  auprc = round(auprc_all$auprc,rounding)
  lo = round(max(0,auprc_all$auprc.CI[1]),rounding)
  up = round(min(1,auprc_all$auprc.CI[2]),rounding)
  return(list(auprc=auprc,auprc.CI=c(lo,up)))
}



###-###-###-###-###-###-###-###
###   weightByCase.half()   ###
###-###-###-###-###-###-###-###

#DESCRIPTION
#This function is used to make a vector of weighted values that corresponds to
#your samples. This weighting scheme makes your case class 50% of the weight,
#and then weights equally amongst the rest of the classes. If you explicitly set
#"otherNames", then any classes that you don't include (i.e. classes not
#included in caseNames OR otherNames) will have a weight of 1 (in theory, this
#should mean that they are "unweighted").

#PARAMETERS
#labelVec - vector of labels that lists the disease group associated with each sample
#caseNames - name of the class(es) that you're considering to be your case (e.g. "Malaria")
#otherNames - vector of all the other classes. If left blank, the function will use the levels of labelVec to figure out what the non-case names are

#RETURN VALUE
#numerical vector of weights for each sample

#REQUIRED PACKAGES: None

weightByCase.half <- function(labelVec, caseNames, otherNames=NULL){
  if(is.null(otherNames)){
    labelVec = factor(labelVec)
    labelVec=droplevels(labelVec)
    otherNames=levels(labelVec)[-c(which(levels(labelVec) %in% caseNames))]
  }

  catNums = catRatio = rep(0,length(otherNames)+1)
  catNums[1] = length(labelVec[labelVec %in% caseNames])
  for(i in 1:length(otherNames)){
    catNums[i+1] = length(labelVec[labelVec == otherNames[i]])
  }

  caseNum=catNums[1]
  otherNum=sum(catNums[2:(length(otherNames)+1)])
  catRatio[1]=otherNum/caseNum
  for(i in 1:length(otherNames)){
    catRatio[i+1] = (1/catNums[i+1])*(otherNum/length(otherNames))
  }

  weightVec=labelVec
  levels(weightVec)=unique(c(levels(weightVec),catRatio,1))
  weightVec[weightVec %in% caseNames] = catRatio[1]
  for(i in 1:length(otherNames)){
    weightVec[weightVec == otherNames[i]] = catRatio[i+1]
  }
  weightVec[!(weightVec %in% catRatio)]=1
  weightVec=as.numeric(as.character(weightVec))

  return(weightVec)
}



##-###-###-###-###-###-###-#
###   weightAllEqual()   ###
##-###-###-###-###-###-###-#

#DESCRIPTION
#This function is used to make a vector of weighted values that corresponds to
#your samples. This weighting scheme weights all classes equally. If you
#explicitly set "classNames", then any classes that you don't include will have
#a weight of 1 (in theory, this should mean that they are "unweighted").

#PARAMETERS
#labelVec - vector of labels that lists the disease group associated with each sample
#classNames - vector of all your classes. If left blank, the function will use the levels of labelVec to figure out what the class names are

#RETURN VALUE
#numerical vector of weights for each sample

#REQUIRED PACKAGES: None

weightAllEqual <- function(labelVec, classNames=NULL){
  if(is.null(classNames)){
    labelVec = factor(labelVec)
    labelVec=droplevels(labelVec)
    classNames=levels(labelVec)
  }

  catNums = catRatio = rep(0,length(classNames))
  for(i in 1:length(classNames)){
    catNums[i] = length(labelVec[labelVec == classNames[i]])
  }

  for(i in 1:length(classNames)){
    catRatio[i] = (1/catNums[i])*(sum(catNums)/length(classNames))
  }

  weightVec=labelVec
  levels(weightVec)=unique(c(levels(weightVec),catRatio,1))
  for(i in 1:length(classNames)){
    weightVec[weightVec == classNames[i]] = catRatio[i]
  }
  weightVec[!(weightVec %in% catRatio)]=1
  weightVec=as.numeric(as.character(weightVec))

  return(weightVec)
}



###-###-###-###-###-###-###-###
###   getGeneNums.lasso()   ###
###-###-###-###-###-###-###-###

#DESCRIPTION
#This function is used to figure out how many genes (i.e. non-zero coefficients)
#are present in a lasso model at various values of lambda

#PARAMETERS
#lasso.obj - the lasso object
#lambda.vals - vector of all your lambda values

#RETURN VALUE
#integer indicating the number of genes that have non-zero coefficients in the
#provided lasso object at the given lambda

#REQUIRED PACKAGES: glmnet

getGeneNums.lasso <- function(lasso.obj,lambda.vals=NULL){
  if(is.null(lambda.vals)){
    lambda_vals = c(lasso.obj$lambda.min,lasso.obj$lambda.1se)
  }
  coef = coef(lasso.obj, s = lambda.vals)
  gene_nums = rep(0,length(lambda.vals))
  for(i in 1:length(lambda.vals)){
    gene_nums[i] = sum(!(coef[-1,i]==0))
  }
  return(gene_nums)
}



###-###-###-###-###-###-###-###-##-#
###   .getNonzeroGenes.lasso()   ###
###-###-###-###-###-###-###-###-##-#

#DESCRIPTION
#This function is used to make a vector of all of the genes with non-zero
#coefficients that are present in a lasso model at a defined value of lambda

#PARAMETERS
#lasso.obj - the lasso object
#lambda - the lambda value for which you want to get the non-zero genes

#RETURN VALUE
#upgenes - character vector of the upregulated genes
#downgenes - character vector of the downregulated genes
#nonzeroGenes - named vector of the coefficients of the non-zero genes

#REQUIRED PACKAGES: glmnet

.getNonzeroGenes.lasso <- function(lasso.obj,lambda=lasso.obj$lambda.1se){
  sigcoef = coef(lasso.obj, lambda)
  nonzeroGenes = sigcoef[-1,1]
  nonzeroGenes = nonzeroGenes[nonzeroGenes!=0]
  upgenes = names(nonzeroGenes)[which(nonzeroGenes>0)]
  downgenes = names(nonzeroGenes)[which(nonzeroGenes<0)]
  return(list(upgenes=upgenes,downgenes=downgenes,nonzeroGenes=nonzeroGenes))
}



###-###-###-###-###-###-###-##
###   metaROCsummaries()   ###
###-###-###-###-###-###-###-##

#DESCRIPTION
#Get relevant statistical read-outs from the output of getMetaROCStats()

#PARAMETERS
#stats - output of getMetaROCStats()
#title - title of the plot
#method - method used to compute summary meta-statistics

#RETURN VALUE
#prints messages and makes plots that show the alpha and beta from the output of getMetaROCStats()

#REQUIRED PACKAGES: rmeta

metaROCsummaries <- function(stats, title, method="random"){
  require(rmeta)
  alpha <- with(stats, meta.summaries(d=tstar_alpha, se=SE_alpha*sqrt(N), method=method))
  print(alpha)
  alpha$names <- rownames(stats)
  plot(alpha, main=paste0("Summary of Alpha, ", title) )

  beta <- with(stats, meta.summaries(d=tstar_beta, se=SE_beta*sqrt(N), method=method))
  print(beta)
  beta$names <- rownames(stats)
  plot(beta, main=paste0("Summary of Beta, ", title), xlim=c(-1,1) )
  # plot(newMetaROC(alpha=alpha$summary, beta=beta$summary), type='l', main=title)
  # lines(newMetaROC(alpha=alpha$summary + 1.96*alpha$se.summary, beta=beta$summary + 1.96*beta$se.summary), type='l', lty=2)
  # lines(newMetaROC(alpha=alpha$summary - 1.96*alpha$se.summary, beta=beta$summary - 1.96*beta$se.summary), type='l', lty=2)
}



#............................................................................................................#
##### Functions for Making Figures #####
#............................................................................................................#

###-###-###-###-###-###-###-##-#
###   .makeROCplotHelper()   ###
###-###-###-###-###-###-###-##-#

#DESCRIPTION
#Creates a Receiver Operating curve for a variety of methods that produce a
#vector of numeric predictions. This is the less user-friendly helper function
#version.

#PARAMETERS
#prediction.type - method used to obtain predictions. Lasso for prediction.type="lasso", cross-validated lasso for prediction.type="cv.lasso",
#                  and standardized difference of geometric means for prediction.type="geometric.smd" or "signature.score"
#class - vector of classifications for your samples, with 0 representing your negative class (controls) and 1 representing your positive class (cases)
#lasso.obj - For prediction.type="lasso" or "cv.lasso", the lasso model used to generate the predictions
#geneMtx - For prediction.type="lasso" or "cv.lasso", matrix of gene expression values (genes in rows, samples in columns)
#lambda.vals - For prediction.type="lasso" or "cv.lasso", the lambda values for each lasso model that you want to use (default is lasso.obj$lambda.min, lasso.obj$lambda.1se)
#scoreList - For prediction.type="geometric.smd", a list of vectors containing each set of scores
#gene.nums - For prediction.type="geometric.smd", the number of genes used to generate each vector of scores

#GRAPHICAL PARAMETERS
#title - title of the figure
#subtitle - subtitle of the figure
#legend.names - if you don't want the default legend labels (i.e. the number of genes in each filterObject) then you can provide alternative names as a character vector
#textSize - use this to easily increase or decrease the size of all the text in the plot
#rounding - how many digits to round the AUC and CI to
#colors - vector of colors for the lines
#legend - if true, a legend will be included
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

#REQUIRED PACKAGES: ROCR, MetaIntegrator (some use cases require glmnet)

.makeROCplotHelper <- function(prediction.type, class, lasso.obj, geneMtx, lambda.vals=NULL, scoreList, gene.nums=NULL,
                               title=NULL, subtitle=NULL, legend.names=NULL, textSize=NULL, rounding=3, colors=c(2,4,6,5,3,7:10), legend=TRUE, ROC.lty=1,
                               ROC.lwd=1, background.color="white", grid.marks=0.1, grid.color="gray93",
                               grid.lty=1, grid.lwd=0.9, diag.lty=2, legend.lty=0, legend.color="transparent", cex.main=NULL, cex.subtitle=NULL, cex.legend=NULL, ...){
  if(prediction.type=="lasso"){
    require(glmnet)
    if(is.null(lambda.vals)){lambda.vals = c(lasso.obj$lambda.min, lasso.obj$lambda.1se)}
    pred = predict(lasso.obj, newx = t(geneMtx), s = lambda.vals, type = "response")
    if(is.null(gene.nums)){gene.nums = getGeneNums.lasso(lasso.obj,lambda.vals)}
  }else if(prediction.type=="cv.lasso"){
    require(glmnet)
    if(is.null(lambda.vals)){lambda.vals = c(lasso.obj$lambda.min, lasso.obj$lambda.1se)}
    pred = predict(lasso.obj, newx = t(geneMtx), s = lambda.vals, type = "response")
    if(is.null(gene.nums)){gene.nums = getGeneNums.lasso(lasso.obj,lambda.vals)}
  }else if(prediction.type=="geometric.smd" || prediction.type=="signature.score"){
    pred = do.call(cbind,(scoreList))
  }else{stop("Prediction type not recognized.")}

  require(ROCR)
  require(MetaIntegrator)
  pred_ROC = prediction(pred, matrix(class, nrow = length(class), ncol = length(pred[1,])))
  perf = performance(pred_ROC, "tpr", "fpr")
  auc.perf = performance(pred_ROC, "auc")
  auc = as.numeric(auc.perf@y.values)
  auc_CI = rep(0,length(pred[1,]))
  for(i in 1:length(pred[1,])){
    CI = .getAUC_CI(auc[i],sum(class==1),sum(class==0))
    lo = max(0,round(CI[1], rounding))
    up = min(1,round(CI[2], rounding))
    auc_CI[i] = paste(" (95% CI ", lo, " - ", up, ")", sep = "")
  }
  if(is.null(legend.names)){
    gene.auc = paste(gene.nums, " genes: AUC=", round(auc, rounding), auc_CI, sep = "")
  }else{
    gene.auc = paste(legend.names, ": AUC=", round(auc, rounding), auc_CI, sep = "")
  }

  if(is.null(cex.main) && is.null(textSize)){
    plot(perf, main = title, ylim=c(0,1),xlim=c(0,1), xlab="False Positive Rate (1-Specificity)", ylab="True Positive Rate (Sensitivity)", ...)
  }else if(!is.null(cex.main) && is.null(textSize)){
    plot(perf, main = title, cex.main = cex.main, ylim=c(0,1),xlim=c(0,1), xlab="False Positive Rate (1-Specificity)", ylab="True Positive Rate (Sensitivity)", ...)
  }else if(!is.null(cex.main) && !is.null(textSize)){
    plot(perf, main = title, ylim = c(0, 1), xlim = c(0,1), cex.lab = textSize, cex.axis = textSize, cex.main = cex.main, cex.sub = textSize, xlab="False Positive Rate (1-Specificity)", ylab="True Positive Rate (Sensitivity)", ...)
  }else{
    plot(perf, main = title, ylim = c(0, 1), xlim = c(0,1), cex.lab = textSize, cex.axis = textSize, cex.main = textSize, cex.sub = textSize, xlab="False Positive Rate (1-Specificity)", ylab="True Positive Rate (Sensitivity)", ...)
  }
  
  if(is.null(cex.subtitle) && is.null(textSize)){
    mtext(subtitle, line = 0.5)
  }else if(!is.null(cex.subtitle)){
    mtext(subtitle, line = 0.5, cex=cex.subtitle)
  }else{
    mtext(subtitle, line = 0.5, cex=textSize*0.85)
  }
  
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = background.color)
  abline(h=seq(0, 1, grid.marks), v=seq(0, 1, grid.marks), col=grid.color, lty=grid.lty, lwd=grid.lwd)
  plot(perf, col = as.list(colors), add = T, lty=ROC.lty, lwd=ROC.lwd)
  lines(c(0,1),c(0,1),col = "black", lty = diag.lty)

  if(legend){
    if(is.null(cex.legend) && is.null(textSize)){
      legend("bottomright", inset = 0.03, legend = gene.auc, fill = colors, box.lty = legend.lty, bg=legend.color)
    }else if(!is.null(cex.legend)){
      legend("bottomright", inset = 0.03, legend = gene.auc, fill = colors, box.lty = legend.lty, cex = cex.legend, bg=legend.color)
    }else{
      legend("bottomright", inset = 0.03, legend = gene.auc, fill = colors, box.lty = legend.lty, cex = textSize * 0.85, bg=legend.color)
    }
  }
}



###-###-###-###-###-###-###-##-#
###   .makePRCplotHelper()   ###
###-###-###-###-###-###-###-##-#

#DESCRIPTION
#Creates a Precision-Recall curve for a variety of methods that produce a vector
#of numeric predictions. This is the less user-friendly helper function version.

#PARAMETERS
#prediction.type - method used to obtain predictions. Lasso for prediction.type="lasso", cross-validated lasso for prediction.type="cv.lasso",
#                  and standardized difference of geometric means for prediction.type="geometric.smd" or "signature.score"
#class - vector of classifications for your samples, with 0 representing your negative class (controls) and 1 representing your positive class (cases)
#lasso.obj - For prediction.type="lasso" or "cv.lasso", the lasso model used to generate the predictions
#geneMtx - For prediction.type="lasso" or "cv.lasso", matrix of gene expression values (genes in rows, samples in columns)
#lambda.vals - For prediction.type="lasso" or "cv.lasso", the lambda values for each lasso model that you want to use (default is lasso.obj$lambda.min, lasso.obj$lambda.1se)
#scoreList - For prediction.type="geometric.smd", a list of vectors containing each set of scores
#gene.nums - For prediction.type="geometric.smd", the number of genes used to generate each vector of scores

#GRAPHICAL PARAMETERS
#title - title of the figure
#subtitle - subtitle of the figure
#legend.names - if you don't want the default legend labels (i.e. the number of genes in each filterObject) then you can provide alternative names as a character vector
#rounding - how many digits to round the AUC and CI to
#colors - vector of colors for the lines
#legend - if true, a legend will be included
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

#REQUIRED PACKAGES: ROCR, zoo, pracma (some use cases require glmnet)

.makePRCplotHelper <- function(prediction.type, class, lasso.obj, geneMtx, lambda.vals, scoreList, gene.nums=NULL,
                               title=NULL, subtitle=NULL, legend.names=NULL, textSize = NULL, rounding=3, colors=c(2,4,6,5,3,7:10), legend=TRUE, PRC.lty=1,
                               PRC.lwd=1, background.color="white", grid.marks=0.1, grid.color="gray93",
                               grid.lty=1, grid.lwd=0.9, legend.lty=0, legend.color="transparent", cex.main=NULL, cex.subtitle=NULL, cex.legend=NULL, ...){
  if(prediction.type=="lasso"){
    require(glmnet)
    if(is.null(lambda.vals)){lambda.vals = c(lasso.obj$lambda.min, lasso.obj$lambda.1se)}
    pred = predict(lasso.obj, newx = t(geneMtx), s = lambda.vals, type = "response")
    if(is.null(gene.nums)){gene.nums = getGeneNums.lasso(lasso.obj,lambda.vals)}
  }else if(prediction.type=="cv.lasso"){
    require(glmnet)
    if(is.null(lambda.vals)){lambda.vals = c(lasso.obj$lambda.min, lasso.obj$lambda.1se)}
    pred = predict(lasso.obj, newx = t(geneMtx), s = lambda.vals, type = "response")
    if(is.null(gene.nums)){gene.nums = getGeneNums.lasso(lasso.obj,lambda.vals)}
  }else if(prediction.type=="geometric.smd" || prediction.type=="signature.score"){
    pred = do.call(cbind,(scoreList))
  }else{stop("Prediction type not recognized.")}

  require(ROCR)
  require(zoo)
  require(pracma)
  pred_PRC = prediction(pred, matrix(class, nrow = length(class), ncol = length(pred[1,])))
  perf = performance(pred_PRC, "prec", "rec")
  auprc_vals = auprc_CI = rep(0,length(pred[1,]))
  for(i in 1:length(pred[1,])){
    #check for NaN bug
    if(is.nan(perf@y.values[[i]][1])){
      perf@y.values[[i]][1] = 1
    }
    #add extra points at (1,0) and (0,1) so that it looks nicer
    perf@x.values[[i]] = c(0,perf@x.values[[i]],1)
    perf@y.values[[i]] = c(1,perf@y.values[[i]],0)
    auprc=.calcauprc(perf@x.values[[i]],perf@y.values[[i]],class)
    auprc_vals[i] = auprc$auprc
    lo = max(0,round(auprc$auprc.CI[1], rounding))
    up = min(1,round(auprc$auprc.CI[2], rounding))
    auprc_CI[i] = paste(" (95% CI ", lo, " - ", up, ")", sep = "")
  }
  if(is.null(legend.names)){
    gene.auprc = paste(gene.nums, " genes: AUPRC=", round(auprc_vals, rounding), auprc_CI, sep = "")
  }else{
    gene.auprc = paste(legend.names, ": AUPRC=", round(auprc_vals, rounding), auprc_CI, sep = "")
  }

  if(is.null(cex.main) && is.null(textSize)){
    plot(perf, main = title, ylim=c(0,1),xlim=c(0,1), ...)
  }else if(!is.null(cex.main) && is.null(textSize)){
    plot(perf, main = title, cex.main = cex.main, ylim=c(0,1),xlim=c(0,1), ...)
  }else if(!is.null(cex.main) && !is.null(textSize)){
    plot(perf, main = title, ylim = c(0, 1), xlim = c(0,1), cex.lab = textSize, cex.axis = textSize, cex.main = cex.main, cex.sub = textSize, ...)
  }else{
    plot(perf, main = title, ylim = c(0, 1), xlim = c(0,1), cex.lab = textSize, cex.axis = textSize, cex.main = textSize, cex.sub = textSize, ...)
  }
  
  if(is.null(cex.subtitle) && is.null(textSize)){
    mtext(subtitle, line = 0.5)
  }else if(!is.null(cex.subtitle)){
    mtext(subtitle, line = 0.5, cex=cex.subtitle)
  }else{
    mtext(subtitle, line = 0.5, cex=textSize*0.85)
  }
  
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = background.color)
  abline(h=seq(0, 1, grid.marks), v=seq(0, 1, grid.marks), col=grid.color, lty=grid.lty, lwd=grid.lwd)
  plot(perf, col = as.list(colors), add = T, lty=PRC.lty, lwd=PRC.lwd)
  #lines(c(1,0),c(0,1),col = "black", lty = diag.lty)
  if(legend){
    if(is.null(cex.legend) && is.null(textSize)){
      legend("bottomleft", inset = 0.03, legend = gene.auprc, fill = colors, box.lty = legend.lty, bg=legend.color)
    }else if(!is.null(cex.legend)){
      legend("bottomleft", inset = 0.03, legend = gene.auprc, fill = colors, box.lty = legend.lty, cex = cex.legend, bg=legend.color)
    }else{
      legend("bottomleft", inset = 0.03, legend = gene.auprc, fill = colors, box.lty = legend.lty, cex = textSize * 0.85, bg=legend.color)
    }
  }
}



###-###-###-###-###-###-###-###-###-#
###   .makeClassROCplotHelper()   ###
###-###-###-###-###-###-###-###-###-#

#DESCRIPTION
#Creates a ROC plot that assesses the performance of your signature and makes a
#ROC curve comparing your positive class (i.e. your cases) against each of your
#other classes (i.e. your controls). This is the less user-friendly helper
#function version.

#PARAMETERS
#prediction.type - method used to obtain predictions. Lasso for prediction.type="lasso", cross-validated lasso for prediction.type="cv.lasso",
#                  and standardized difference of geometric means for prediction.type="geometric.smd" or "signature.score"
#class - vector of classifications for your samples, with 0 representing your negative class (controls) and 1 representing your positive class (cases)
#lasso.obj - For prediction.type="lasso" or "cv.lasso", the lasso model used to generate the predictions
#geneMtx - matrix of gene expression values (genes in rows, samples in columns)
#lambda.vals - For prediction.type="lasso" or "cv.lasso", the lambda value for the lasso model that you want to use (default is lasso.obj$lambda.1se)
#pos.genes - vector of the positive genes in your signature
#neg.genes - vector of the negative genes in your signature
#labelVec - character vector or factor of descriptive class labels (e.g. infecting pathogen)
#caseNames - if left blank, this will be a character vector of the class(es) that correspond to your cases (i.e. samples labeled with a 1 in class).
#            If the user wants to exclude any of these case classes, then this can be set as the names of all of the case classes that should be included
#otherNames - if left blank, this will be a character vector of the class(es) that correspond to your controls (i.e. samples labeled with a 0 in class).
#             If the user wants to exclude any of these control classes, then this can be set as the names of all of the control classes that should be included

#GRAPHICAL PARAMETERS
#title - title of the plot
#textSize - use this to easily increase or decrease the size of all the text in the plot
#full.labels - if TRUE, then the full names of the case categories will be used in the legend, if FALSE then a generic "Cases" will be used instead
#rounding - how many digits to round the AUC and CI to
#colors - colors of the lines
#plot.unweighted - should the unweighted ROC curve be plotted (i.e. the normal ROC curve without class divisions)
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
#An ROC plot showing diagnostic performance across multiple designated classes

#REQUIRED PACKAGES: ROCR, MetaIntegrator (some use cases require glmnet)

.makeClassROCplotHelper <- function(prediction.type, class, lasso.obj, geneMtx, lambda.vals=NULL, pos.genes, neg.genes, labelVec, caseNames=NULL, otherNames=NULL,
                                    title, textSize = NULL, full.labels = TRUE, rounding=3, colors=c(2,4,6,5,3,7:10), plot.unweighted=TRUE, lwd=NULL, unweighted.lwd=1,
                                    unweighted.col="grey", background.color="white", grid.marks=0.1, grid.color="gray93",
                                    grid.lty=1, grid.lwd=0.9, diag.lty=2, legend.lty=0, legend.color="transparent", cex.main = NULL, cex.subtitle=NULL, cex.legend=NULL, ...){
  if(is.factor(labelVec)){labelVec=droplevels(labelVec)}
  if(!is.null(caseNames) && any(!c(caseNames %in% unique(labelVec[class==1])))){
    stop("If caseNames is set by the user, then all of the provided names must be labels that are used for case samples (i.e. samples that are labeled with 1 in the provided class vector)")
  }
  if(!is.null(otherNames) && any(!c(otherNames %in% unique(labelVec[class==0])))){
    stop("If otherNames is set by the user, then all of the provided names must be labels that are used for control samples (i.e. samples that are labeled with 0 in the provided class vector)")
  }
  if(!is.null(lambda.vals) && length(lambda.vals) > 1){
    stop("Only one lambda value can be tested at a time.")
  }
  if((prediction.type=="lasso" || prediction.type=="cv.lasso") && is.null(lambda.vals)){
    lambda.vals = lasso.obj$lambda.1se
  }
  require(ROCR)
  require(MetaIntegrator)

  #figure out your case labels and non-case labels
  labelVec = factor(labelVec)
  labelVec = droplevels(labelVec)
  if(is.null(caseNames)){
    #if(!is.factor(labelVec)){labelVec=factor(labelVec)} #delete this line if there's no bugs
    caseNames=as.character(unique(labelVec[class==1]))
  }
  if(is.null(otherNames)){
    #if(!is.factor(labelVec)){labelVec=factor(labelVec)} #delete this line if there's no bugs
    otherNames=levels(labelVec)[-c(which(levels(labelVec) %in% caseNames))]
  }

  if(length(colors) < length(otherNames)){
    mult = ceiling(length(otherNames)/length(colors))
    colors = rep(colors, mult)
  }
  if(!is.null(lwd) && (length(lwd) != length(otherNames))){stop("If lwd is set by the user, it must be the same length as otherNames")}
  classIndex = lapply(otherNames, function(other){
    return(append(which(labelVec==other),which(labelVec %in% caseNames)))
  })

  #get predictions/scores
  gene.nums=NULL
  if(prediction.type=="lasso"){
    pred = as.numeric(predict(lasso.obj, newx = t(geneMtx), s = lambda.vals, type = "response"))
    gene.nums = getGeneNums.lasso(lasso.obj,lambda.vals)
  }else if(prediction.type=="cv.lasso"){
    pred = as.numeric(predict(lasso.obj, newx = t(geneMtx), s = lambda.vals, type = "response"))
    gene.nums = getGeneNums.lasso(lasso.obj,lambda.vals)
  }else if(prediction.type=="geometric.smd" || prediction.type=="signature.score"){
    pred = as.numeric(getGeneScores(geneMtx,pos.genes,neg.genes))
    gene.nums = length(pos.genes) + length(neg.genes)
  }

  #get AUCs
  aucVec=loVec=upVec=rep(0,length(classIndex))
  for(i in 1:length(classIndex)){
    MI_ROC = calculateROC(as.numeric(as.character(class[classIndex[[i]] ])), pred[classIndex[[i]] ])
    aucVec[i] = round(MI_ROC$auc, rounding)
    loVec[i] = max(0,round(MI_ROC$auc.CI[1], rounding))
    upVec[i] = min(1,round(MI_ROC$auc.CI[2], rounding))
  }
  if(plot.unweighted){
    MI_ROC_unweighted = calculateROC(as.numeric(as.character(class)), pred)
    unweighted.auc = round(MI_ROC_unweighted$auc, rounding)
    unweighted.lo = max(0,round(MI_ROC_unweighted$auc.CI[1], rounding))
    unweighted.up = min(1,round(MI_ROC_unweighted$auc.CI[2], rounding))
  }

  #get plotting values
  class.perf = lapply(classIndex, function(index){
    pred_ROC = prediction(pred[index], class[index])
    perf = performance(pred_ROC, "tpr", "fpr")
    return(perf)
  })
  if(plot.unweighted){
    pred_ROC.unweighted = prediction(pred, class)
    unweighted.perf = performance(pred_ROC.unweighted, "tpr", "fpr")
  }

  #make labels
  blurbVec = rep("",length(classIndex))
  for(i in 1:length(classIndex)){
    if(full.labels){
      blurbVec[i]=sprintf("%s vs. %s: AUC=%s (95%% CI %s-%s)",paste(caseNames,collapse=", "),otherNames[i],aucVec[i],loVec[i],upVec[i])
    }else{
      blurbVec[i]=sprintf("%s vs. %s: AUC=%s (95%% CI %s-%s)","Cases",otherNames[i],aucVec[i],loVec[i],upVec[i])
    }
  }
  if(plot.unweighted){
    blurbVec=append(blurbVec,sprintf("Unweighted AUC=%s (95%% CI %s-%s)",unweighted.auc,unweighted.lo,unweighted.up))
  }

  avgauc=round(sum(aucVec)/length(classIndex),rounding)
  avgauc.CI = .aucCI(avgauc,class)
  avgauc.lo = max(0,round(avgauc.CI[1],rounding))
  avgauc.up = min(1,round(avgauc.CI[2],rounding))

  if(is.null(lwd)){lwd=rep(1,length(classIndex))}
  
  if(is.null(cex.main) && is.null(textSize)){
    plot(class.perf[[1]], col=colors, main = title, xlab="False Positive Rate (1-Specificity)", ylab="True Positive Rate (Sensitivity)", ...)
  }else if(!is.null(cex.main) && is.null(textSize)){
    plot(class.perf[[1]], col=colors, main = title, cex.main=cex.main, xlab="False Positive Rate (1-Specificity)", ylab="True Positive Rate (Sensitivity)", ...)
  }else if(!is.null(cex.main) && !is.null(textSize)){
    plot(class.perf[[1]], col=colors, main = title, cex.main=cex.main, cex.lab = textSize, cex.axis = textSize, cex.sub = textSize, xlab="False Positive Rate (1-Specificity)", ylab="True Positive Rate (Sensitivity)", ...)
  }else{
    plot(class.perf[[1]], col=colors, main = title, cex.main=textSize, cex.lab = textSize, cex.axis = textSize, cex.sub = textSize, xlab="False Positive Rate (1-Specificity)", ylab="True Positive Rate (Sensitivity)", ...)
  }
  
  if(is.null(cex.subtitle) && is.null(textSize)){
    mtext(sprintf("%s Genes; Average AUC=%s (95%% CI %s-%s)",gene.nums,avgauc,avgauc.lo,avgauc.up), line = 0.5)
  }else if(!is.null(cex.subtitle)){
    mtext(sprintf("%s Genes; Average AUC=%s (95%% CI %s-%s)",gene.nums,avgauc,avgauc.lo,avgauc.up), line = 0.5, cex = cex.subtitle)
  }else{
    mtext(sprintf("%s Genes; Average AUC=%s (95%% CI %s-%s)",gene.nums,avgauc,avgauc.lo,avgauc.up), line = 0.5, cex = textSize*0.85)
  }
  
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = background.color)
  abline(h=seq(0, 1, grid.marks), v=seq(0, 1, grid.marks), col=grid.color, lty=grid.lty, lwd=grid.lwd)
  for(i in 1:length(class.perf)){
    plot(class.perf[[i]], col=colors[i], add = T, lwd=lwd[i])
  }
  legend.colors = colors[1:length(otherNames)]
  if(plot.unweighted){
    plot(unweighted.perf, col=unweighted.col, add = T, lwd=unweighted.lwd)
    legend.colors = c(legend.colors,unweighted.col)
  }
  lines(c(0,1),c(0,1),col = "black", lty = diag.lty)
  if(is.null(cex.legend) && is.null(textSize)){
    legend("bottomright", inset = 0.03, legend = blurbVec, fill = legend.colors, box.lty = legend.lty, bg=legend.color)
  }else if(!is.null(cex.legend)){
    legend("bottomright", inset = 0.03, legend = blurbVec, fill = legend.colors, box.lty = legend.lty, cex = cex.legend, bg=legend.color)
  }else{
    legend("bottomright", inset = 0.03, legend = blurbVec, fill = legend.colors, box.lty = legend.lty, cex = textSize * 0.85, bg=legend.color)
  }
}



###-###-###-###-###-###-###-###-###-#
###   .makeClassPRCplotHelper()   ###
###-###-###-###-###-###-###-###-###-#

#DESCRIPTION
#Creates a PRC plot that assesses the performance of your signature and makes a
#PRC curve comparing your positive class (i.e. your cases) against each of your
#other classes (i.e. your controls). This is the less user-friendly helper
#function version.

#NOTE: I should probably make a version of this later that works with multiple caseNames, and compares each case category to all control categories

#PARAMETERS
#prediction.type - method used to obtain predictions. Lasso for prediction.type="lasso", cross-validated lasso for prediction.type="cv.lasso",
#                  and standardized difference of geometric means for prediction.type="geometric.smd" or "signature.score"
#class - vector of classifications for your samples, with 0 representing your negative class (controls) and 1 representing your positive class (cases)
#lasso.obj - For prediction.type="lasso" or "cv.lasso", the lasso model used to generate the predictions
#geneMtx - matrix of gene expression values (genes in rows, samples in columns)
#lambda.vals - For prediction.type="lasso" or "cv.lasso", the lambda value for the lasso model that you want to use (default is lasso.obj$lambda.1se)
#pos.genes - vector of the positive genes in your signature
#neg.genes - vector of the negative genes in your signature
#labelVec - character vector or factor of descriptive class labels (e.g. infecting pathogen)
#caseNames - if left blank, this will be a character vector of the class(es) that correspond to your cases (i.e. samples labeled with a 1 in class).
#            If the user wants to exclude any of these case classes, then this can be set as the names of all of the case classes that should be included
#otherNames - if left blank, this will be a character vector of the class(es) that correspond to your controls (i.e. samples labeled with a 0 in class).
#             If the user wants to exclude any of these control classes, then this can be set as the names of all of the control classes that should be included

#GRAPHICAL PARAMETERS
#title - title of the plot
#textSize - use this to easily increase or decrease the size of all the text in the plot
#full.labels - if TRUE, then the full names of the case categories will be used in the legend, if FALSE then a generic "Cases" will be used instead
#rounding - how many digits to round the AUC and CI to
#colors - colors of the lines
#plot.unweighted - should the unweighted PRC curve be plotted (i.e. the normal PRC curve without class divisions)
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

#REQUIRED PACKAGES: ROCR, MetaIntegrator, zoo, pracma (some use cases require glmnet)

.makeClassPRCplotHelper <- function(prediction.type, class, lasso.obj, geneMtx, lambda.vals=NULL, pos.genes, neg.genes, labelVec, caseNames=NULL, otherNames=NULL,
                                    title, textSize = NULL, full.labels = TRUE, rounding=3, colors=c(2,4,6,5,3,7:10), plot.unweighted=TRUE, lwd=NULL, unweighted.lwd=1,
                                    unweighted.col="grey", background.color="white", grid.marks=0.1, grid.color="gray93",
                                    grid.lty=1, grid.lwd=0.9, legend.lty=0, legend.color="transparent", cex.main = NULL, cex.subtitle=NULL, cex.legend=NULL, ...){
  if(is.factor(labelVec)){labelVec=droplevels(labelVec)}
  if(!is.null(caseNames) && any(!c(caseNames %in% unique(labelVec[class==1])))){
    stop("If caseNames is set by the user, then all of the provided names must be labels that are used for case samples (i.e. samples that are labeled with 1 in the provided class vector)")
  }
  if(!is.null(otherNames) && any(!c(otherNames %in% unique(labelVec[class==0])))){
    stop("If otherNames is set by the user, then all of the provided names must be labels that are used for control samples (i.e. samples that are labeled with 0 in the provided class vector)")
  }
  if(!is.null(lambda.vals) && length(lambda.vals) > 1){
    stop("Only one lambda value can be tested at a time.")
  }
  if((prediction.type=="lasso" || prediction.type=="cv.lasso") && is.null(lambda.vals)){
    lambda.vals = lasso.obj$lambda.1se
  }
  require(ROCR)
  require(MetaIntegrator)
  require(zoo)
  require(pracma)

  #figure out your case labels and non-case labels
  labelVec = factor(labelVec)
  labelVec = droplevels(labelVec)
  if(is.null(caseNames)){
    caseNames=as.character(unique(labelVec[class==1]))
  }
  if(is.null(otherNames)){
    otherNames=levels(labelVec)[-c(which(levels(labelVec) %in% caseNames))]
  }

  if(length(colors) < length(otherNames)){
    mult = ceiling(length(otherNames)/length(colors))
    colors = rep(colors, mult)
  }
  if(!is.null(lwd) && (length(lwd) != length(otherNames))){stop("If lwd is set by the user, it must be the same length as otherNames")}
  classIndex = lapply(otherNames, function(other){
    return(append(which(labelVec==other),which(labelVec %in% caseNames)))
  })

  #get predictions/scores
  gene.nums=NULL
  if(prediction.type=="lasso"){
    pred = as.numeric(predict(lasso.obj, newx = t(geneMtx), s = lambda.vals, type = "response"))
    gene.nums = getGeneNums.lasso(lasso.obj,lambda.vals)
  }else if(prediction.type=="cv.lasso"){
    pred = as.numeric(predict(lasso.obj, newx = t(geneMtx), s = lambda.vals, type = "response"))
    gene.nums = getGeneNums.lasso(lasso.obj,lambda.vals)
  }else if(prediction.type=="geometric.smd" || prediction.type=="signature.score"){
    pred = as.numeric(getGeneScores(geneMtx,pos.genes,neg.genes))
    gene.nums = length(pos.genes) + length(neg.genes)
  }

  #get AUPRCs
  auprcVec=loVec=upVec=rep(0,length(classIndex))
  class.perf = list()
  for(i in 1:length(classIndex)){
    pred_PRC = prediction(pred[classIndex[[i]] ], as.numeric(as.character(class[classIndex[[i]] ])))
    perf = performance(pred_PRC, "prec", "rec")
    #check for NaN bug
    if(is.nan(perf@y.values[[1]][1])){
      perf@y.values[[1]][1] = 1
    }
    #add extra points at (1,0) and (0,1) so that it looks nicer
    perf@x.values[[1]] = c(0,perf@x.values[[1]],1)
    perf@y.values[[1]] = c(1,perf@y.values[[1]],0)
    class.perf[[i]] = perf
    auprc=.calcauprc(perf@x.values[[1]],perf@y.values[[1]],class)
    auprcVec[i] = auprc$auprc
    loVec[i] = max(0,round(auprc$auprc.CI[1], rounding))
    upVec[i] = min(1,round(auprc$auprc.CI[2], rounding))
  }

  if(plot.unweighted){
    pred_PRC.unweighted = prediction(pred, class)
    unweighted.perf = performance(pred_PRC.unweighted, "prec", "rec")
    #check for NaN bug
    if(is.nan(unweighted.perf@y.values[[1]][1])){
      unweighted.perf@y.values[[1]][1] = 1
    }
    #add extra points at (1,0) and (0,1) so that it looks nicer
    unweighted.perf@x.values[[1]] = c(0,unweighted.perf@x.values[[1]],1)
    unweighted.perf@y.values[[1]] = c(1,unweighted.perf@y.values[[1]],0)
    unweighted.info =.calcauprc(unweighted.perf@x.values[[1]],unweighted.perf@y.values[[1]],class)
    unweighted.auprc = round(unweighted.info$auprc,rounding)
    unweighted.lo = max(0,round(unweighted.info$auprc.CI[1], rounding))
    unweighted.up = min(1,round(unweighted.info$auprc.CI[2], rounding))
  }

  avgauprc = sum(auprcVec)/length(classIndex)
  total.n.pos = sum(class==1)
  avgauprc.lo = max(0,round(avgauprc-1.96*sqrt((avgauprc*(1-avgauprc))/total.n.pos),rounding))
  avgauprc.up = min(1,round(avgauprc+1.96*sqrt((avgauprc*(1-avgauprc))/total.n.pos),rounding))
  avgauprc = round(avgauprc,rounding)
  auprcVec = round(auprcVec,rounding)

  #make labels
  blurbVec = rep("",length(classIndex))
  for(i in 1:length(classIndex)){
    if(full.labels){
      blurbVec[i]=sprintf("%s vs. %s: AUPRC=%s (95%% CI %s-%s)",paste(caseNames,collapse=", "),otherNames[i],auprcVec[i],loVec[i],upVec[i])
    }else{
      blurbVec[i]=sprintf("%s vs. %s: AUPRC=%s (95%% CI %s-%s)","Cases",otherNames[i],auprcVec[i],loVec[i],upVec[i])
    }
  }
  if(plot.unweighted){
    blurbVec=append(blurbVec,sprintf("Unweighted AUPRC=%s (95%% CI %s-%s)",unweighted.auprc,unweighted.lo,unweighted.up))
  }

  if(is.null(lwd)){lwd=rep(1,length(classIndex))}
  
  if(is.null(cex.main) && is.null(textSize)){
    plot(class.perf[[1]], col=colors, main = title, ...)
  }else if(!is.null(cex.main) && is.null(textSize)){
    plot(class.perf[[1]], col=colors, main = title, cex.main=cex.main, ...)
  }else if(!is.null(cex.main) && !is.null(textSize)){
    plot(class.perf[[1]], col=colors, main = title, cex.main=cex.main, cex.lab = textSize, cex.axis = textSize, cex.sub = textSize, ...)
  }else{
    plot(class.perf[[1]], col=colors, main = title, cex.main=textSize, cex.lab = textSize, cex.axis = textSize, cex.sub = textSize, ...)
  }
  
  if(is.null(cex.subtitle) && is.null(textSize)){
    mtext(sprintf("%s Genes; Average AUPRC=%s (95%% CI %s-%s)",gene.nums,avgauprc,avgauprc.lo,avgauprc.up), line = 0.5)
  }else if(!is.null(cex.subtitle)){
    mtext(sprintf("%s Genes; Average AUPRC=%s (95%% CI %s-%s)",gene.nums,avgauprc,avgauprc.lo,avgauprc.up), line = 0.5, cex = cex.subtitle)
  }else{
    mtext(sprintf("%s Genes; Average AUPRC=%s (95%% CI %s-%s)",gene.nums,avgauprc,avgauprc.lo,avgauprc.up), line = 0.5, cex = textSize*0.85)
  }
  
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = background.color)
  abline(h=seq(0, 1, grid.marks), v=seq(0, 1, grid.marks), col=grid.color, lty=grid.lty, lwd=grid.lwd)
  for(i in 1:length(class.perf)){
    plot(class.perf[[i]], col=colors[i], add = T, lwd=lwd[i])
  }
  legend.colors = colors[1:length(otherNames)]
  if(plot.unweighted){
    plot(unweighted.perf, col=unweighted.col, add = T, lwd=unweighted.lwd)
    legend.colors = c(legend.colors,unweighted.col)
  }
  #lines(c(0,1),c(0,1),col = "black", lty = diag.lty)
  if(is.null(cex.legend) && is.null(textSize)){
    legend("bottomleft", inset = 0.03, legend = blurbVec, fill = legend.colors, box.lty = legend.lty, bg=legend.color)
  }else if(!is.null(cex.legend)){
    legend("bottomleft", inset = 0.03, legend = blurbVec, fill = legend.colors, box.lty = legend.lty, cex = cex.legend, bg=legend.color)
  }else{
    legend("bottomleft", inset = 0.03, legend = blurbVec, fill = legend.colors, box.lty = legend.lty, cex = textSize * 0.85, bg=legend.color)
  }
}



##-###-###-###-###-###-###-#
###   makeViolinPlot()   ###
##-###-###-###-###-###-###-#

#DESCRIPTION
#Given a gene matrix, a vector of class labels, and the positive and negative
#genes in your signature, this function creates a figure that has violin plots
#(plotted by signature score) for each category you choose and then superimposes
#a dotted line over the violin plots that indicates what the optimal signature
#score cutoff would be (determined using Youden Index).

#PARAMETERS
#geneMtx - matrix of gene expression values (genes in rows, samples in columns)
#class - vector of classifications for your samples, with 0 representing your negative class (controls) and 1 representing your positive class (cases)
#pos.genes - vector of the positive genes in your signature
#neg.genes - vector of the negative genes in your signature
#grouping - factor with labels for each sample indicating what category each sample belongs in
#names - vector containing each label you want to make a violin plot for, in order. if NULL, this is automatically generated from the levels of "grouping"
#labelData - set this as TRUE if you want to see the labels of each of your categories in the colored bars (defined by the "names" parameter)
#title - title of the figure
#display.auc - whether to display the global AUC as part of the figure title
#DZpal - the colors used for your violin plots
#group.min - minimum number of datapoints required for each group to be included in the plot (note that all groups will always be included in the AUC calculation)
#barlayers - number of layers of colored bars
#barwidth - width of the colored bars
#ylim - if you want to explicitly set the ylim
#ylim.lo - use this to increase the upper limit of the y axis
#ylim.hi - use this to decrease the lower limit of the y axis
#labsize - size of the text used for the labels
#xlim - limits of the x axis
#cex.main - size of the title
#ylabadjust - if you want the labels to be moved downwards, because sometimes they are placed too high
#ytitleadjust - if you want the title to be moved downwards
#cutoff.method - what optimal cutoff method to use. Choose from Youden, MinValueSe, MinValueSp, or most other options within OptimalCutpoints (some that don't work can probably be added)

#RETURN VALUE
#A violin plot comparing the signature scores for each class of samples

#REQUIRED PACKAGES: OptimalCutpoints, vioplot

makeViolinPlot <- function(geneMtx, class, pos.genes, neg.genes, grouping, names=NULL, labelData=T, title=NULL, display.auc=T,
                           DZpal = c("firebrick3","#0072B2"), group.min=NULL, barlayers=1, barwidth=0.065, ylim=NULL, ylim.lo=0,
                           ylim.hi=0, labsize=1, xlim=NULL, cex.main=1, ylabadjust=0, ytitleadjust=0, cutoff.method="Youden"){
  test.score <- getGeneScores(geneMtx, pos.genes, neg.genes)
  if(any(is.na(grouping))){grouping[is.na(grouping)] = "NA"} #bugfix
  
  ## get ROC data; this is equivalent to the internal methods
  require(OptimalCutpoints)
  opt <- optimal.cutpoints(X="test.score", status="class", tag.healthy=0,
                           data=data.frame(class, test.score),
                           methods=c("Youden", cutoff.method),
                           control=control.cutpoints(valueSe=0.945))
  AUC <- signif(opt$Youden$Global$measures.acc$AUC, 3)
  cat(sprintf("AUC = %s (90%% CI %s - %s)\n", AUC[1], max(0,AUC[2]), min(1,AUC[3])))
  cat(sprintf("At optimal cutoff, summary sensitivity = %s, summary specificity = %s  @ cutoff %s\n",
              signif(opt$Youden$Global$optimal.cutoff$Se, 3),
              signif(opt$Youden$Global$optimal.cutoff$Sp, 3),
              opt$Youden$Global$optimal.cutoff$cutoff))
  # cat(sprintf("Choosing sensitivity = %s, then specificity = %s @ cutoff %s\n",
  #             signif(opt$MinValueSe[[1]]$optimal.cutoff$Se, 3),
  #             signif(opt$MinValueSe[[1]]$optimal.cutoff$Sp, 3),
  #             opt$MinValueSe[[1]]$optimal.cutoff$cutoff))
  if(cutoff.method=="Youden"){
    cutoff <- opt$Youden$Global$optimal.cutoff$cutoff[1]
  }else{
    cutoff <- opt[[cutoff.method]][[1]]$optimal.cutoff$cutoff[1]
  }
  
  if(!is.null(group.min)){
    remove.groups = names(which(table(grouping)<group.min))
    remove.entries = which(grouping %in% remove.groups)
    geneMtx = geneMtx[,-remove.entries]
    class = class[-remove.entries]
    grouping = grouping[-remove.entries]
  }
  
  #make plot
  extra <- length(class)*0.01
  if(is.null(xlim)){xlim=c(0-extra, length(class)+extra)}
  classTab <- table(class)
  ylim.new = .getYLim(test.score)
  ylim.new[1] = ylim.new[1]-ylim.lo;ylim.new[2] = ylim.new[2]+ylim.hi
  
  plot("",xlim=xlim,xaxs="i",xaxt='n',xlab="",ylim =ylim.new,ylab="Score")
  if(!is.null(title)){
    if(display.auc){title(main = paste(title, sprintf(" Global AUC = %s (95%% CI %s - %s)\n", AUC[1], max(0,AUC[2]), min(1,AUC[3]))), line=0-ytitleadjust, cex.main=cex.main)}
    else{title(main = title, line=1-ytitleadjust, cex.main=cex.main)}
  }
  colBarWidth <- (par("usr")[4]-par("usr")[3])*barwidth
  
  if(is.null(names)){
    grouping = factor(grouping)
    grouping = droplevels(grouping)
    names = levels(grouping)
  }
  
  #trying to fix regular expression bug
  names.regx = names
  #names = gsub("\\\n|\\(|\\)|\\+|\\.|\\^|\\$|\\*|\\?","",names)
  #grouping = gsub("\\\n|\\(|\\)|\\+|\\.|\\^|\\$|\\*|\\?","",grouping)
  names = gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", names) #escape all special characters
  grouping = gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", grouping) #escape all special characters
  
  dataPal <- .getDataPal(names)

  curr.start = 0
  for(i in 1:length(names)){
    ## establish widths
    index <- grep(paste0("^",names[i],"$"), grouping)
    if(length(index) == 0){index <- which(grouping==names[i])}
    x.start <- curr.start
    x.end <- curr.start + length(index)
    curr.start <- x.end
    classTab <- table(class[index])
    classTab <- classTab[as.character(sort(unique(class)))]
    classTab[is.na(classTab)] <- 0
    classTab=classTab[length(classTab):1]
    if(classTab[1]>0 && classTab[2]>0){
      index1=index[class[index]==1]
      index2=index[class[index]==0]
    }else{index1=index2=index}
    
    ## boxes on top showing group name
    col <- dataPal[grep(paste0("^",names[i],"$"), names)]
    odd1 <- (grep(paste0("^",names[i],"$"), names) - 1) %% barlayers  ## how many layers of colorbars
    if(length(odd1) == 0){odd1 <- (which(names==names[i]) - 1) %% barlayers}
    y.upper <- par("usr")[4] - (0.1 + colBarWidth*odd1)
    rect(x.start, y.upper-colBarWidth, x.end, y.upper,
         col=col, lwd=1, lend=2, border=NA)
    ## grey vertical lines
    lines(x=rep(x.start, 2), y=c(par("usr")[3], y.upper-colBarWidth),
          lwd=1, col="grey30", lty=2)
    
    ## violinplots of score
    library(vioplot)
    if(classTab[1]>0){
      capture.output(vioplot(test.score[index1], add=T,
                             at=x.start+(classTab[1])/2, wex=classTab[1]*0.8,
                             range=0, lwd=1.2, col=DZpal[1], pchMed="-", colMed="white"),file="/dev/null")
      if(classTab[1]==1){
        points(x=x.start+(classTab[1])/2,y=test.score[index1],
               pch=21,col="black",bg=DZpal[1])
      }
    }
    if(classTab[2]>0){
      capture.output(vioplot(test.score[index2], add=T,
                             at=x.start+classTab[1]+(classTab[2])/2, wex=classTab[2]*0.8,
                             range=0, lwd=1.2, col=DZpal[2], pchMed="-", colMed="white"),file="/dev/null")
      if(classTab[2]==1){
        points(x=x.start+classTab[1]+(classTab[2])/2,y=test.score[index2],
               pch=21,col="black",bg=DZpal[2])
      }
    }
    if(labelData){
      odd2 <- (grep(paste0("^",names[i],"$"), names) - 1) %% barlayers
      if(length(odd2) == 0){odd2 <- (which(names==names[i]) - 1) %% barlayers}
      text(x=(x.start+x.end)/2, y=par("usr")[4]-colBarWidth*(odd2+0.8)-0.04+ylabadjust,
           labels=names.regx[i],col="black", cex=labsize)#, srt=45) ##rotates labels
    } else {
      my.labels = grep(paste0("^",names[i],"$"), names)
      if(length(my.labels) == 0){my.labels <- which(names==names[i])}
      text(x=mean(c(x.start,x.end)), y=y.upper-colBarWidth/2, labels=my.labels,cex=1.2)
    }
  }
  lines(x=rep(x.end, 2), y=c(par("usr")[3], par("usr")[4]-colBarWidth-0.1),
        lwd=1, col="grey30", lty=2)
  
  ## Global cutoff
  abline(h=cutoff, lwd=2, col="grey10", lty=4)
}



###-###-###-###-###-###-##
###   makeGenePlot()   ###
###-###-###-###-###-###-##

#DESCRIPTION
#Given a gene matrix, a vector of class labels, and the positive and negative
#genes in your signature, this function creates a figure that has the individual
#expression values of each gene in your signature, for each sample in your data

#PARAMETERS
#geneMtx - matrix of gene expression values (genes in rows, samples in columns)
#class - vector of classifications for your samples, with 0 representing your negative class (controls) and 1 representing your positive class (cases)
#pos.genes - vector of the positive genes in your signature
#neg.genes - vector of the negative genes in your signature
#grouping - factor with labels for each sample indicating what category each sample belongs in
#names - vector containing each label you want to have a category for, in order. if NULL, this is automatically generated from the levels of "grouping"
#title - title of the figure
#pointSize - the size of the points
#ylim - use this if you want to manually set the y limits of the plot
#labelData - set this as true if you want to see the labels of each of your categories in the colored bars (defined by the "names" parameter)
#DZpal - the colors used for the bars on the bottom that indicate cases and controls
#barlayers - number of layers of colored bars
#barwidth - width of the colored bars
#lowerbarwidth - width of the bars on the bottom that indicate cases and controls
#greyScaleGenes - use this if you want the genes to be in greyscale
#legend - if true, a small legend will be made showing the color for each gene

#RETURN VALUE
#A plot showing the expression of the genes in your signature across all the samples

#REQUIRED PACKAGES: RColorBrewer

makeGenePlot <- function(geneMtx, class, genes, grouping, names=NULL, title, pointSize=0.6, ylim=NULL, labelData=T,
                         DZpal = c("firebrick3","#0072B2"), barlayers=1, barwidth=0.065, lowerbarwidth=0.75, greyScaleGenes=F, legend=F){
  geneMtx = geneMtx[genes,,drop=F]

  library(RColorBrewer)
  ### Create main plot; this shows individual genes
  #this color stuff isn't completely tested
  if(nrow(geneMtx)>8){
    mypal <- colorRampPalette(brewer.pal(9, "Set1"), bias=0.5)(nrow(geneMtx))
  } else if (greyScaleGenes){
    mypal <- colorRampPalette(c("grey20", "grey70"))(nrow(geneMtx))
  } else {
    mypal <- brewer.pal(nrow(geneMtx), "Set1")
  }

  extra <- length(class)*0.01
  if(is.null(ylim)){ylim <- .getYLim(data.matrix(geneMtx))}
  plot("", xlim=c(1-extra, length(class)+extra),
       xaxs="i", xaxt='n', xlab="", ylab = "Gene Exp.",
       ylim = ylim, main=title)

  if(is.null(names)){
    grouping = factor(grouping)
    grouping = droplevels(grouping)
    names = levels(grouping)
  }
  #reorder geneMtx based on the grouping/names
  geneMtx.new = do.call(cbind,lapply(names,function(name){
    index = which(grouping==name)
    return(geneMtx[,index])
  }))

  dataPal <- .getDataPal(names)
  colBarWidth <- (par("usr")[4]-par("usr")[3])*barwidth

  curr.start=0
  for(i in 1:length(names)){
    ## establish widths
    index <- grep(paste0("^",names[i],"$"), grouping)
    x.start <- curr.start - 0.5
    x.end <- curr.start + length(index)
    curr.start <- x.end

    ## boxes on top showing dataset ID
    col <- dataPal[grep(paste0("^",names[i],"$"), names)]
    odd1 <- (grep(paste0("^",names[i],"$"), names) - 1) %% barlayers
    y.upper <- par("usr")[4] - (0.1 + colBarWidth*odd1)
    rect(x.start, y.upper-colBarWidth, x.end+0.5, y.upper,
         col=col, lwd=1, lend=2, border=NA)

    ## dataset ID for legend
    #untested
    if(labelData){
      text(x=(x.start+x.end)/2, y=y.upper-colBarWidth/2, labels=names[i])
    }

    ## boxes on bottom showing cases and controls
    ## controlled by DZpal
    classTab <- table(class[index])
    classTab <- classTab[as.character(sort(unique(class)))]
    classTab[is.na(classTab)] <- 0
    classTab=classTab[length(classTab):1]
    y.lower <- par("usr")[3]

    ## boxes are automatically expanded to present classes
    for(j in 1:length(unique(class))){
      rect(x.start + sum(classTab[0:(j-1)]), y.lower,
           x.start + sum(classTab[0:j]), y.lower+colBarWidth*lowerbarwidth,
           col=DZpal[j], lwd=1, lend=2, border=NA)
    }

    lines(x=rep(x.start, 2), y=c(par("usr")[3], y.upper-colBarWidth),
          lwd=1, col="grey30", lty=2)
  }

  ## last vertical dividing line
  lines(x=rep(length(class), 2),
        y=c(par("usr")[3], par("usr")[4]-colBarWidth-0.1),
        lwd=1, col="grey50", lty=2)

  ## Actual points
  for(k in 1:nrow(geneMtx.new)){
    points(1:ncol(geneMtx.new), y=geneMtx.new[k,], col=mypal[k], pch=20, cex=pointSize)
  }

  if(legend){
    par(font=3)
    legend("bottomright", legend=rownames(geneMtx.new),
           fill=mypal, bg="white")#, cex=1, bty="n", y.intersp=0.7)
    par(font=1)
  }
}



##-###-###-###-###-###-###-#
###   makeForestPlot()   ###
##-###-###-###-###-###-###-#

#DESCRIPTION
#Given a gene matrix and a vector of class labels, this function creates a
#forest plot for your data. This version of the function is used when you want
#to look at the expression of one gene across a number of groups (e.g. datasets,
#diseases, etc.).

#PARAMETERS
#geneMtx - matrix of gene expression values (genes in rows, samples in columns)
#genes - the gene you want to use
#grouping - factor with labels for each sample indicating what category each sample belongs in
#names - vector containing each label you want to have a category for, in order. if NULL, this is automatically generated from the levels of "grouping"
#labels - vector containing the labels that you want to show up on the plot, in order. if NULL, "names" is used.
#xlim - if you want to manually set the x limits of the plot
#title - if you don't want the title to be the name of your gene, use this to indicate an alternate title
#box - color of the box
#lines - color of the whiskers on the box
#zero - color of the dotted line at 0
#summary - color of the summary diamond
#text - color of the labels

#RETURN VALUE
#A forest plot comparing the effect size of a gene across datasets/diseases/etc.

#REQUIRED PACKAGES: rmeta

makeForestPlot <- function(geneMtx, gene, grouping, names=NULL, labels=NULL, xlim=NULL, title=NULL,
                           box = "blue",lines = "lightblue",zero = "black", summary = "orange", text = "black", ...){
  if(!(gene %in% rownames(geneMtx))){stop("The requested gene does not exist in the given gene matrix")}
  
  if(is.null(names)){
    grouping = factor(grouping)
    grouping = droplevels(grouping)
    names = levels(grouping)
  }
  if(is.null(labels)){labels=names}
  
  # #trying to fix regular expression bug
  # names_regx = gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", names) #escape all special characters
  # grouping_regx = gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", grouping) #escape all special characters
  
  ES.g = rep(0,length(names));ES.se = rep(0,length(names))
  for(i in 1:length(names)){
    #index = grep(paste0("^",names[i],"$"), grouping)
    index = which(grouping == names[i]) #I think this should work - not sure why I was using regex before
    test.class=rep(0,length(grouping))
    test.class[index]=1
    ES = .getES(as.numeric(geneMtx[gene,]),test.class)
    ES.g[i] = ES[9]
    ES.se[i] = ES[10]
  }
  
  if(is.null(title)){
    rmeta::metaplot(ES.g,ES.se,1/ES.se^2,labels=labels,xlab="Effect Size",ylab="",cex=1.2,main=gene,xlim=xlim,
                    colors = rmeta::meta.colors(box = box,lines = lines,zero = zero, summary = summary, text = text), ...)
  } else{
    rmeta::metaplot(ES.g,ES.se,1/ES.se^2,labels=labels,xlab="Effect Size",ylab="",cex=1.2,main=title,xlim=xlim,
                    colors = rmeta::meta.colors(box = box,lines = lines,zero = zero, summary = summary, text = text), ...)
  }
}



###-###-###-###-###-###-###-###-###
###   makeForestPlot.pooled()   ###
###-###-###-###-###-###-###-###-###

#DESCRIPTION
#Given a gene matrix, a vector of class labels, and a list of genes, this
#function creates a forest plot for your data. This version of the function is
#used when all of your data is pooled and you want to look at the expression of
#each of your genes between your cases and controls.

#PARAMETERS
#geneMtx - matrix of gene expression values (genes in rows, samples in columns)
#class - vector of classifications for your samples, with 0 representing your negative class (controls) and 1 representing your positive class (cases)
#chosengenes - list of the genes you want to use
#xlim - if you want to manually set the x limits of the plot
#title - title of the figure
#box - color of the box
#lines - color of the whiskers on the box
#zero - color of the dotted line at 0
#summary - color of the summary diamond
#text - color of the labels
#out.missing - if true, the names of the missing genes will be output

#RETURN VALUE
#A forest plot comparing the effect size of your genes between cases and controls

#REQUIRED PACKAGES: None

makeForestPlot.pooled <- function(geneMtx, class, chosengenes, xlim=NULL, title=NULL, include.summary=TRUE,box = "blue",lines = "lightblue",
                                  zero = "black", summary = "orange", text = "black",out.missing=TRUE, ...){
  if(out.missing){
    missing = chosengenes[!(chosengenes %in% rownames(geneMtx))]
    if(length(missing)>0){
      cat("Missing these genes:",missing,"\n")
    }
  }
  chosengenes = chosengenes[chosengenes %in% rownames(geneMtx)]
  
  all.ES = apply(geneMtx[chosengenes,,drop=F],1,function(x) .getES(as.numeric(x),as.numeric(as.character(class))))
  ES.g = all.ES[9,]
  ES.se = all.ES[10,]
  ES.pooled = pool.inverseVar(t(data.frame(genes=ES.g)),t(data.frame(genes=ES.se)),method="random")
  if(include.summary){
    rmeta::metaplot(ES.g,ES.se,1/ES.se^2,labels=chosengenes,xlab="Effect Size",ylab="",cex=1.2,main=title,xlim=xlim,
                    colors = rmeta::meta.colors(box = box,lines = lines,zero = zero, summary = summary, text = text),
                    summn=ES.pooled[2],sumse=ES.pooled[3],sumnn=1/ES.pooled[3]^2, ...)
  }else{
    rmeta::metaplot(ES.g,ES.se,1/ES.se^2,labels=chosengenes,xlab="Effect Size",ylab="",cex=1.2,main=title,xlim=xlim,
                    colors = rmeta::meta.colors(box = box,lines = lines,zero = zero, summary = summary, text = text), ...)
  }
}



###-###-###-###-###-###-###-###-###-###-#
###   makeForestPlot.multidataset()   ###
###-###-###-###-###-###-###-###-###-###-#
#No description cuz i'm sleepy

makeForestPlot.multidataset <- function(pooledDataObject, gene, title=NULL, include.summary=TRUE, xlim=NULL,box = "#ef4f4f",
                                        lines="#ffcda3", zero = "black", summary = "#74c7b8", text = "black", ...){
  EStab = data.frame(do.call(rbind,lapply(pooledDataObject,function(gse){
    if(gene %in% rownames(gse$genes)){
      stats = .getES(as.numeric(gse$genes[gene,]),gse$class)
      return(c(stats[9:10],name=gse$formattedName))
    }else{
      return(c(g=NA,se.g=NA,name=gse$formattedName))
    }
  })))
  EStab$g = as.numeric(as.character(EStab$g))
  EStab$se.g = as.numeric(as.character(EStab$se.g))
  ES.pooled = pool.inverseVar(t(data.frame(genes=EStab$g)),t(data.frame(genes=EStab$se.g)),method="random")
  
  if(include.summary){
    rmeta::metaplot(EStab$g,EStab$se.g,1/EStab$se.g^2,labels=EStab$name,xlab="Effect Size",ylab="",cex=1.2,main=title,xlim=xlim,
                    colors = rmeta::meta.colors(box = box,lines=lines,zero = zero, summary = summary, text = text),
                    summn=ES.pooled[2],sumse=ES.pooled[3],sumnn=1/ES.pooled[3]^2, ...)
  }else{
    rmeta::metaplot(EStab$g,EStab$se.g,1/EStab$se.g^2,labels=EStab$name,xlab="Effect Size",ylab="",cex=1.2,main=title,xlim=xlim,
                    colors = rmeta::meta.colors(box = box,lines=lines,zero = zero, summary = summary, text = text), ...)
  }
}



##-###-###-###-###-###-###-#
###   discoValidPlot()   ###
##-###-###-###-###-###-###-#

#DESCRIPTION
#After splitting your data into discovery and validation, use this to make sure
#that the distributions of disease samples within each cohort is fairly similar

#PARAMETERS
#discoLabels - labels for all of your samples in discovery
#validLabels - labels for all of your samples in validation
#diseaseNames - list of diseases that you have labeled in discoLabels and validLabels
#colors - the colors for each category

#RETURN VALUE
#A bar plot showing the distribution of samples in discovery and validation

#REQUIRED PACKAGES: None

discoValidPlot <- function(discoLabels, validLabels, diseaseNames, colors=NULL){
  disco_diseases = vector("numeric");valid_diseases = vector("numeric")
  for(i in 1:length(diseaseNames)){
    disco_diseases[i] = length(discoLabels[discoLabels == diseaseNames[i]])
    valid_diseases[i] = length(validLabels[validLabels == diseaseNames[i]])
  }

  if(is.null(colors)){
    colors = sample(rainbow(length(diseaseNames)))
  }
  both_disease = as.matrix(cbind(valid_diseases/sum(valid_diseases),disco_diseases/sum(disco_diseases)))
  barplot(both_disease, beside = F, horiz = T, col = colors, main = "Distribution of Disease Samples in Discovery and Validation",
          names.arg = c("Validation","Discovery"))
  legend("left", legend = diseaseNames, fill = colors, ncol = 7, cex = 0.75, inset = 0.01, bg="white")
}



###-###-###-###-###-###-###-###-###-##
###   makeClassErrorsPlot.pred()   ###
###-###-###-###-###-###-###-###-###-##

#DESCRIPTION
#This function analyzes a set of predictions and determines the class-specific
#error rates at a given cutoff and plots them.

#PARAMETERS
#method - method used to get the optimal cutoff
#pred - vector of predictions
#class - vector of classifications for your samples, with 0 representing your negative class (controls) and 1 representing your positive class (cases)

#predVec - vector of predictions
#labelVec - vector of descriptive class labels (e.g. infecting pathogen)
#class - vector of classifications for your samples, with 0 representing your negative class (controls) and 1 representing your positive class (cases)
#cutoff - cutoff used to convert predictions into classifications
#subsetData - if true, then the data will be subsetted according to keepIndex
#keepIndex - if subsetData is true, then this describes the indices of the samples that you want to include in the plot
#title - title of the plot
#min - minimum number of samples needed for each class to be included in the plot
#vertjust - vertical justifcation of the text that displays the total number of samples in each class
#ylim - limits of the y axis
#ybreaks - breakpoints on the y axis

#RETURN VALUE
#A bar plot showing classification error

#REQUIRED PACKAGES: ggplot2

makeClassErrorsPlot.pred <- function(predVec, labelVec, class, cutoff, subsetData=FALSE, keepIndex,
                                     title=NULL, min=5, vertjust=-1, ylim=c(0,1), ybreaks=seq(0,1,0.1)){
  if(!subsetData){keepIndex=TRUE}
  true.labels = class[keepIndex]
  pred.labels = as.numeric(predVec>cutoff)[keepIndex]
  pred.class = labelVec[keepIndex]

  pred.class = droplevels(pred.class)
  pred.error.mat = lapply(levels(pred.class),function(className){
    pathIndex = which(pred.class==className)
    if(length(pathIndex)>=min){
      errRate=sum(pred.labels[pathIndex]!=true.labels[pathIndex])/length(pathIndex)
      errRate=round(errRate,3)
      return(c(className,paste(length(pathIndex),"samples",sep="\n"),errRate))
    }
  })

  pred.error = data.frame(do.call(rbind,pred.error.mat),stringsAsFactors = FALSE)
  colnames(pred.error) = c("Pathogen","Total","Error")
  pred.error$Error=as.numeric(pred.error$Error)

  ggplot(pred.error,aes(x=Pathogen,y=Error,fill=Pathogen))+geom_bar(stat="identity")+theme(legend.position="none")+
    scale_y_continuous(name="Error Rate",limits=ylim,breaks=ybreaks)+geom_text(aes(label=Total),vjust=vertjust)+ggtitle(title)
}



##-###-###-###-###-###-###-###-#
###   makeSummaryROCplot()   ###
##-###-###-###-###-###-###-###-#

#DESCRIPTION
#This function makes a ROC plot for multiple datasets with a summary ROC curve included

#Statistics are based on the Kester and Buntinx Method, from (Kester and
#Buntnix, Med Decis Making, 2000). Methods have also been added by Tim Sweeney
#(2015) for better estimates in cases with low numbers of tpr/fpr values

#PARAMETERS
#DatasetList - list of datasets. Each dataset must have a matrix of gene expression values with genes in rows
#              and samples in columns (stored in $genes), a character vector of labels for each sample (stored
#              in $pheno$group) and the name of the dataset (stored in $formattedName)
#pos.genes - vector of the positive genes in your signature
#neg.genes - vector of the negative genes in your signature
#caseNames - name of the class(es) that you're considering to be your case (e.g. "Malaria")
#ROC.stats - use this if you want to manually pass in the output of getMetaROCStats()
#title - title of the plot
#size - size of the text/legend/etc.
#rounding - how many digits to round the AUC and CI to
#method - method used to compute summary meta-statistics
#numCores - number of CPUs to use if parallel computing is desired
#minPoints - minimum number of points required for bootstrap to be used
#bootReps - number of bootstrap iterations
#auc1.thresh - if the AUC of a dataset is above this threshold, then it is treated as if the AUC were 1

#RETURN VALUE
#A ROC plot showing diagnostic performance summarized across a number of ROC curves

#REQUIRED PACKAGES: MetaIntegrator, ROCR, boot, ggplot2, rmeta

makeSummaryROCplot <- function(DatasetList,pos.genes,neg.genes,caseNames,ROC.stats=NULL,title=NULL,size=14,rounding=3,method="random",numCores=1,
                               minPoints=5,bootReps=1000,auc1.thresh=0.99){
  require(ggplot2)
  require(rmeta)
  if(is.null(ROC.stats)){
    ROC.stats = getMetaROCStats(DatasetList,pos.genes,neg.genes,caseNames,numCores=numCores,minPoints=minPoints,bootReps=bootReps,auc1.thresh=auc1.thresh)
  }
  legendText = sprintf("%s AUC=%s (95%% CI %s-%s)",rownames(ROC.stats),round(ROC.stats$AUC,rounding),
                       round(ROC.stats$CI.lower,rounding),round(ROC.stats$CI.upper,rounding))
  ROC.stats <- na.omit(ROC.stats)
  if(nrow(ROC.stats)>1){
    alpha <- with(ROC.stats, meta.summaries(d=tstar_alpha, se=SE_alpha*sqrt(N), method=method))
    beta <- with(ROC.stats, meta.summaries(d=tstar_beta, se=SE_beta*sqrt(N), method=method))
  } else {
    alpha <- with(ROC.stats, data.frame(summary=tstar_alpha, se.summary=SE_alpha))
    beta <- with(ROC.stats, data.frame(summary=tstar_beta, se.summary=SE_beta))
  }
  roc.summ <- newMetaROC(alpha=alpha$summary, beta=beta$summary)
  roc.lower.b.hi <- newMetaROC(alpha=alpha$summary - 1.96*alpha$se.summary, beta=beta$summary + 1.96*beta$se.summary, print=F)
  roc.lower.b.lo <- newMetaROC(alpha=alpha$summary - 1.96*alpha$se.summary, beta=beta$summary - 1.96*beta$se.summary, print=F)
  roc.upper.b.hi <- newMetaROC(alpha=alpha$summary + 1.96*alpha$se.summary, beta=beta$summary + 1.96*beta$se.summary, print=F)
  roc.upper.b.lo <- newMetaROC(alpha=alpha$summary + 1.96*alpha$se.summary, beta=beta$summary - 1.96*beta$se.summary, print=F)

  roc.lower <- data.frame(FPR = roc.lower.b.hi$FPR, TPR = pmin(roc.lower.b.hi$TPR, roc.lower.b.lo$TPR, na.rm=T) )
  roc.upper <- data.frame(FPR = roc.upper.b.hi$FPR, TPR = pmax(roc.upper.b.hi$TPR, roc.upper.b.lo$TPR, na.rm=T) )

  auc.summ <- round(aucROCframe(roc.summ),rounding)
  auc.lower <- round(aucROCframe(roc.lower),rounding)
  auc.upper <- round(aucROCframe(roc.upper),rounding)

  plotData = do.call(rbind,lapply(DatasetList, function(dataset){
    scores = getGeneScores(dataset$genes,pos.genes,neg.genes)
    class = makeClassVector(dataset$pheno$group,caseNames)
    pred_ROC = prediction(scores, class)
    perf = performance(pred_ROC, "tpr", "fpr")
    return(data.frame(name=dataset$formattedName,FPR=unlist(perf@x.values),TPR=unlist(perf@y.values)))
  }))
  plotData = rbind(plotData,data.frame(name="Summary",FPR=roc.summ$FPR,TPR=roc.summ$TPR))

  legendText = c(legendText,sprintf("Summary AUC=%s (95%% CI %s-%s)",auc.summ,auc.lower,auc.upper))
  #force color to dark grey on palette
  hues <- seq(15, 375, length=length(DatasetList)+1)
  dataPal = c(hcl(h=hues, l=65, c=100)[1:length(DatasetList)],"grey25")

  ggplot(plotData, aes(x = FPR, y = TPR)) +
    #Have this first to put it on bottom
    geom_ribbon(data=cbind(summ=roc.summ, lower=roc.lower, upper=roc.upper),
                aes(x= lower.FPR, y=upper.TPR, ymin=lower.TPR, ymax=upper.TPR),
                color="gray25", fill="gray75") +
    geom_line(aes(colour = name)) +
    geom_segment(x=0, y=0, xend=1, yend=1, linetype=5, color="grey20") +
    scale_x_continuous("False Positive Rate (1-Specificity)",breaks=seq(0,1,0.2)) +
    scale_y_continuous("True Positive Rate (Sensitivity)",breaks=seq(0,1,0.2)) +
    scale_color_manual(values=dataPal, labels = legendText) +
    ggtitle(title) +
    theme(text = element_text(size=size)) +
    theme(plot.title = element_text(size=size),
          axis.title.x = element_text(size=size),
          axis.title.y = element_text(size=size, angle=90),
          legend.justification=c(1,0),
          legend.position=c(1,0),
          legend.title = element_blank(),
          legend.key = element_blank(),
          legend.text=element_text(size=size),
          axis.text.x = element_text(size=size),
          axis.text.y = element_text(size=size)) +
    guides(colour=guide_legend(override.aes = list(size=2.8))) +
    geom_line(data=roc.summ, aes(x=FPR, y=TPR), size=1.3, color="gray20")
  ## If just want to define upper and lower lines without fill:
  #geom_line(data=roc.upper, aes(x=FPR, y=TPR), size=1.1, linetype=3, color="gray15") +
  #geom_line(data=roc.lower, aes(x=FPR, y=TPR), size=1.1, linetype=3, color="gray15")
}


#............................................................................................................#
##### Background Functions #####
#............................................................................................................#

###-###-###-###-###-##
###   geomMean()   ###
###-###-###-###-###-##

#DESCRIPTION
#Function for calculating the geometric mean of a numeric vector.
#Note: the method used here for dealing with 0s does not perform well
#when the input vector has a large proportion of 0 values.

#PARAMETERS
#x - vector of numeric values
#na.rm - whether or not to remove NA values

#REQUIRED PACKAGES: None

geomMean <- function(x, na.rm = FALSE){
  if (!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
    warning("argument is not numeric or logical: returning NA")
    return(as.numeric(NA))
  }
  if (na.rm){
    x <- x[!is.na(x)]
  }
  if (any(x < 0)){
    stop("'x' contains negative value(s)")
  }
  if(any(x == 0)){ #if I don't this, having any 0s will result in the geometric mean being 0
    x = x[-c(which(x == 0))]
    if(length(x) == 0){ #this means all values were 0
      return(0)
    }
  }

  ### direct method causes overflow errors, use log method instead
  ### return(prod(x)^(1/length(x)))
  return(exp(sum(log(x))/length(x)))
}



###-###-###-###-###-###
###  .calcauprc()   ###
###-###-###-###-###-###

#DESCRIPTION
#function for calculating the auprc of a precision-recall curve, as well as the
#confidence interval

#PARAMETERS
#precision - vector of precision values
#recall - vector of recall values
#class - vector of classifications for your samples, with 0 representing your negative class (controls) and 1 representing your positive class (cases)
#type.CI - method used to calculate the confidence interval. binomial is more balanced, but logit is guaranteed to be between 0 and 1
#conflevel - confidence level of the confidence interval

#REQUIRED PACKAGES: zoo, pracma

.calcauprc <- function(precision, recall, class, type.CI="binomial", conflevel="95%"){
  require(zoo)
  require(pracma)
  NA.vals=union(which(is.na(precision)),which(is.na(recall)))
  if(length(NA.vals)==0){
    auprc = sum(diff(precision)*rollmean(recall,2))
  }else{
    auprc = sum(diff(precision[-NA.vals])*rollmean(recall[-NA.vals],2))
  }

  allconflevels=c("80%","90%","95%","99%")
  allphis=c(1.28,1.64,1.96,2.576)
  if(conflevel %in% allconflevels){
    phi = allphis[which(conflevel==allconflevels)]
  }else{
    return(auprc)
  }

  if(type.CI=="binomial"){
    n.pos = sum(class==1)
    ci.lower = auprc-phi*sqrt((auprc*(1-auprc))/n.pos)
    ci.upper = auprc+phi*sqrt((auprc*(1-auprc))/n.pos)
    if(ci.upper>1){ci.upper=1}
    if(ci.lower<0){ci.lower=0}
  }else if(type.CI=="logit"){
    mu = logit(auprc)
    n.pos = sum(class==1)
    t = 1/sqrt(n.pos*auprc*(1-auprc))
    e=exp(1)
    ci.lower = (e^(mu-phi*t))/(1+e^(mu-phi*t))
    ci.upper = (e^(mu+phi*t))/(1+e^(mu+phi*t))
  }else{
    return(auprc)
  }

  return(list(auprc=auprc,auprc.CI=c(ci.lower,ci.upper)))
}



##-###-###-###-###-###-###-#
###   getCutoff.pred()   ###
##-###-###-###-###-###-###-#

#DESCRIPTION
#This function is used to find the optimal cutoff for a set of predictions.
#It currently supports the Youden cutoff, but will be expanded later

#PARAMETERS
#method - method used to get the optimal cutoff
#pred - vector of predictions
#class - vector of classifications for your samples, with 0 representing your negative class (controls) and 1 representing your positive class (cases)

#REQUIRED PACKAGES: OptimalCutpoints

getCutoff.pred <- function(pred, class, method="Youden"){
  df = data.frame(pred,class)
  colnames(df)=c("scores","class")
  require(OptimalCutpoints)
  opt=optimal.cutpoints(X="scores", status="class", tag.healthy=0, data=df,methods=method)
  cutoff = opt$Youden$Global$optimal.cutoff$cutoff
  cat(sprintf("At optimal cutoff, summary sensitivity = %s, summary specificity = %s  @ cutoff %s\n",
              signif(opt$Youden$Global$optimal.cutoff$Se, 3),signif(opt$Youden$Global$optimal.cutoff$Sp, 3),cutoff))
  return(cutoff)
}



###-###-###-###-###-###-###-#
###   getMetaROCStats()   ###
###-###-###-###-###-###-###-#

#DESCRIPTION
#This function calculates meta-statistics for multiple ROC curves, which can be
#used to make a summary ROC curve for those datasets

#Statistics are based on the Kester and Buntinx Method, from (Kester and
#Buntnix, Med Decis Making, 2000). Methods have also been added by Tim Sweeney
#(2015) for better estimates in cases with low numbers of tpr/fpr values

#PARAMETERS
#DatasetList - list of datasets. Each dataset must have a matrix of gene expression values with genes in rows
#              and samples in columns (stored in $genes), a character vector of labels for each sample (stored
#              in $pheno$group) and the name of the dataset (stored in $formattedName)
#pos.genes - vector of the positive genes in your signature
#neg.genes - vector of the negative genes in your signature
#caseNames - name of the class(es) that you're considering to be your case (e.g. "Malaria")
#auc1.thresh - if the AUC of a dataset is above this threshold, then it is treated as if the AUC were 1
#bootReps - number of bootstrap iterations
#minPoints - minimum number of points required for bootstrap to be used
#numCores - number of CPUs to use if parallel computing is desired

#RETURN VALUE
#data frame listing the sample size, alpha parameter value and standard error,
#beta parameter and standard error, auc value, and the 95% confidence interval
#upper and lower values for each of the ROC curves

#REQUIRED PACKAGES: MetaIntegrator, ROCR, boot

getMetaROCStats <- function(DatasetList, pos.genes, neg.genes, caseNames, auc1.thresh=0.99,
                            bootReps=1000, minPoints=5, numCores=1){
  require(MetaIntegrator)
  require(ROCR)
  nameVec=rep("",length(DatasetList))
  for(i in 1:length(DatasetList)){
    nameVec[i]=DatasetList[[i]]$formattedName
  }

  ROC.stats = do.call(rbind,mclapply(mc.cores=numCores, DatasetList, function(dataset){
    scores = getGeneScores(dataset$genes,pos.genes,neg.genes)
    class = makeClassVector(dataset$pheno$group,caseNames)
    MI_ROC = calculateROC(as.numeric(as.character(class)), as.numeric(scores))
    auc = MI_ROC$auc
    auc.lo = MI_ROC$auc.CI[1]
    auc.up = min(MI_ROC$auc.CI[2],1)
    pred_ROC = prediction(scores, class)
    perf = performance(pred_ROC, "tpr", "fpr")
    P=sum(class==1)
    N=sum(class==0)
    if(auc>auc1.thresh){
      #as defined in the paper. alpha = log((P+0.5)*(N+0.5)/0.25)
      alpha = log((P+0.5)*(N+0.5)/0.25)
      #ASE = sqrt of asymptotic variance over sqrt n
      ASE = (1/0.5 + 1/0.5 + 1/(P+0.5) + 1/(N+0.5))^(-1/2)/sqrt(P+N)
      return(c(N+P, alpha, ASE , 0, 10, auc, auc.lo, auc.up))
    }else if(auc==0){ #not sure why no threshold here
      #added for edge case of auc=0
      alpha = -log((P+0.5)*(N+0.5)/0.25)
      ASE = (1/0.5 + 1/0.5 + 1/(P+0.5) + 1/(N+0.5))^(-1/2)/sqrt(P+N)
      return(c(N+P, alpha, ASE , 0, 10, auc, auc.lo, auc.up))
    }

    sens = unlist(perf@x.values); min1.spec = unlist(perf@y.values)
    sens[min(which(sens==1))] = 0.98 #not convinced this is necessary
    min1.spec[max(which(min1.spec==0))] = 0.02 #not convinced this is necessary
    #weights defined in paper appendix, sends inf/NAN to 0
    weights = (1/unlist(pred_ROC@tp) + 1/unlist(pred_ROC@tn) + 1/unlist(pred_ROC@fp) + 1/unlist(pred_ROC@fn))^(-1)
    weights[is.nan(weights) | is.na(weights)] = 0
    S = log(sens/(1-sens)) + log(min1.spec/(1-min1.spec))
    D = -1*(log(sens/(1-sens)) - log(min1.spec/(1-min1.spec)))

    #if there are enough non-0 weights, then bootstrap
    points = sum(weights != 0)
    #tim made some note about needing points that vary in x and y, not just non-0
    #weights, but then he commented out the corresponding code so it looks like
    #maybe non-0 weights is enough? If there's an issue, can try this code instead:
    #points = sum(!duplicated(sens) & !duplicated(min1.spec))

    if(points > minPoints){
      require(boot)
      bs <- function(formula, data, indices){
        d = data[indices, ]
        fit = lm(formula, data=d, weights=weights)
        return(coef(fit))
      }
      boot.out = boot(data=data.frame(min1.spec,sens,weights,S,D), statistic=bs, formula=D~S, R=bootReps, weights=weights)
      op <- NULL
      for (i in 1:2) op <- rbind(op, imp.moments(boot.out, index = i)$rat)
      std.error <- sqrt(op[, 2]) #this used to be 2L
      return(c(N+P, op[1,1], std.error[1], op[2,1], std.error[2], auc, auc.lo, auc.up))
    } else if (points >= 2){
      #else assign ASE as standard errors for both alpha and beta
      cat("\tCan't bootstrap, points < minSize...\t")
      model <- lm(D~S, data=data.frame(min1.spec,sens,weights,S,D), weights=weights)
      alpha <- model$coefficients[1]
      beta <- 0
      ASE <- (1/0.5 + 1/0.5 + 1/(P+0.5) + 1/(N+0.5))^(-1/2)/sqrt(P+N)
      return(c(N+P, alpha, ASE , beta, ASE, auc, auc.lo, auc.up))
    } else {
      #if only one point in ROC that is not on axis, can't use.
      cat(sprintf(" has <2 points in ROC curve; cannot compute summary stats\n", dataset$formattedName))
      return(c(N+P, rep(NA, 4), auc, auc.lo, auc.up))
    }
  }))

  colnames(ROC.stats)=c("N","tstar_alpha","SE_alpha","tstar_beta","SE_beta","AUC","CI.lower","CI.upper")
  rownames(ROC.stats)=nameVec
  return(data.frame(ROC.stats))
}



##-###-###-###-###-###-#
###   newMetaROC()   ###
##-###-###-###-###-###-#

#DESCRIPTION
#Using the alpha and beta parameters from getMetaROCStats(), this function makes
#a matrix of FPR and TPR value for the summary ROC curve

#PARAMETERS
#alpha - alpha parameter
#beta - beta parameter
#points - how many points to include for the fpr/tpr
#print - whether to print the summary AUC

#RETURN VALUE
#a matrix of FPR and TPR values

#REQUIRED PACKAGES: None

newMetaROC <- function(alpha, beta, points=1000, print=F){
  ## in practice, beta > 0.95 becomes unbounded
  if(abs(beta)>0.95) beta <- 0.95 * sign(beta)
  A <- alpha/(1-beta)
  B <- (1+beta)/(1-beta)

  by <- 1/points
  newspec <- seq(by, 1-by, by)
  newsens <- exp(A+B*log((1-newspec)/newspec)) / (1 + exp(A+B*log((1-newspec)/newspec)) )
  newroc <- data.frame("FPR" = 1-newspec, "TPR" = newsens)

  auc <- aucROCframe(newroc)
  if(print){
    cat("AUC: ", auc, "\n")
  }
  return(newroc)
}



###-###-###-###-###-###-#
###   aucROCframe()   ###
###-###-###-###-###-###-#

#DESCRIPTION
#Extracts the AUC from a matrix of FPR and TPR values

#PARAMETERS
#newroc - matrix where first column is FPR and second column is TPR

#RETURN VALUE
#the AUC value

#REQUIRED PACKAGES: None

aucROCframe <- function(newroc){
  i <- 2:dim(newroc)[1]
  auc <- (newroc[i-1, "FPR"] - newroc[i, "FPR"]) %*% (newroc[i-1, "TPR"] + newroc[i, "TPR"])/2
  return(auc)
}



###-###-###-###-##-#
###   .aucSE()   ###
###-###-###-###-##-#

#DESCRIPTION
#calculates the standard error given an AUC and vector of classification labels

#PARAMETERS
#auc - AUC value
#class - vector of classifications for your samples, with 0 representing your negative class (controls) and 1 representing your positive class (cases)

#RETURN VALUE
#Returns the standard error of the AUC

#REQUIRED PACKAGES: None

.aucSE <- function(auc,class){
  n.pos = sum(class == 1)
  n.neg = sum(class == 0)
  q1 = auc/(2 - auc)
  q2 = (2 * auc^2)/(1 + auc)
  se.auc = sqrt((((auc * (1 - auc))+(n.pos - 1)*(q1 - auc^2))+((n.neg - 1)*(q2 - auc^2)))/(n.pos * n.neg))
  return(se.auc)
}



###-###-###-###-##-#
###   .aucCI()   ###
###-###-###-###-##-#

#DESCRIPTION
#calculates the 95% confidence interval given an AUC and vector of classification labels

#PARAMETERS
#auc - AUC value
#class - vector of classifications for your samples, with 0 representing your negative class (controls) and 1 representing your positive class (cases)

#RETURN VALUE
#Returns the 95% confidence interval of the AUC

#REQUIRED PACKAGES: None

.aucCI <- function(auc,class){
  n.pos = sum(class == 1)
  n.neg = sum(class == 0)
  q1 = auc/(2 - auc)
  q2 = (2 * auc^2)/(1 + auc)
  se.auc = sqrt((((auc * (1 - auc))+(n.pos - 1)*(q1 - auc^2))+((n.neg - 1)*(q2 - auc^2)))/(n.pos * n.neg))
  ci.upper = auc + (se.auc * 1.96)
  ci.lower = auc - (se.auc * 1.96)
  return(c(ci.lower,ci.upper))
}



###-###-###-###-###-###-###-##
###   .getSummROCalpha()   ###
###-###-###-###-###-###-###-##

#DESCRIPTION
#for the Kester and Buntix method, calculates the alpha parameter given the AUC and beta parameter

#PARAMETERS
#auc - AUC value
#beta - beta parameter

#RETURN VALUE
#returns the predicted alpha parameter

#REQUIRED PACKAGES: None

.getSummROCalpha <- function(auc,beta){
  if(auc==0.5){
    alpha = 0
  }
  else if(beta>=0.95){
    alpha = (0.01-log((1-auc)/auc))/(0.51)
  }
  else{
    alpha = (0.01-log((1-auc)/auc))/(0.681-(beta^2)*0.197)
  }

  if(alpha == Inf){
    return(12) #for now, to avoid errors, I return a very high alpha if alpha was Inf
  }else if(alpha == -Inf || alpha < -12)
    return(-12)
  else{
    return(alpha)
  }
}



###-###-###-###-##-#
###   .getES()   ###
###-###-###-###-##-#

#DESCRIPTION
#Gets the Hedges' g effect size given a vector of gene expression values and a vector of class labels

#PARAMETERS
#value - vector of gene expression values (numeric is expected, but if you have a 1D dataframe (e.g. from running apply on a 2D dataframe) then it will be converted to numeric)
#class - vector of classifications for your samples, with 0 representing your negative class (controls) and 1 representing your positive class (cases)

#RETURN VALUE
#n1 - the number of samples in class 1
#m1 - the mean expression value for class 1
#sd1 - the standard deviation of the expression values for class 1
#n2 - the number of samples in class 2
#m2 - the mean expression value for class 2
#sd2 - the standard deviation of the expression values for class 2
#diff - the difference between means
#pooled.sd - the pooled standard deviation
#g - Hedges' g effect size
#se.g - standard error of the effect size

#REQUIRED PACKAGES: None

.getES <- function(value, class){
  ## function to calculate basic statistics for two-class comparison for a gene
  if(class(value)=="data.frame"){value=as.numeric(value)}
  stopifnot(identical(length(value),length(class)))

  na.vals = which(is.na(value))
  if(length(na.vals)>=1){
    value = value[-na.vals]
    class = class[-na.vals]
  }

  x <- value[which(class==1)]
  y <- value[which(class==0)]

  n1 <- length(x); n2 <- length(y)
  if( n1 < 2 | n2 < 2 )
    return( c(n1=NA, m1=NA, sd1=NA,
              n2=NA, m2=NA, sd2=NA,
              diff=NA, pooled.sd=NA,
              g=NA, se.g=NA) )

  m1   <- mean(x); m2 <- mean(y)
  diff <- m1 - m2
  sd1  <- sd(x);  sd2 <- sd(y)
  sp   <- sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2 )/(n1 + n2 - 2))
  J   <- 1 - 3/(4*(n1 + n2) - 9) #bias correction factor
  g    <- J * diff/sp

  #using the unbiased estimate of the variance of g from http://journals.sagepub.com/doi/pdf/10.3102/1076998606298034
  df = n1+n2-2
  se.g = sqrt((n1+n2)/(n1*n2) + (1-(df-2)/(df*J^2)) * g^2)
  #old version (produces a different result, but only very slightly)
  #se.g <- sqrt((n1+n2)/(n1*n2) + 0.5*g^2 /(n1+n2-3.94))

  return(c(n1=n1, m1=m1, sd1=sd1, n2=n2, m2=m2, sd2=sd2,
           diff=diff, pooled.sd=sp, g=g, se.g=se.g))
}



###-###-###-###-##-#
###   .getFC()   ###
###-###-###-###-##-#

#DESCRIPTION
#Gets the fold change given a vector of gene expression values and a vector of class labels

#PARAMETERS
#value - vector of gene expression values (numeric is expected, but if you have a 1D dataframe (e.g. from running apply on a 2D dataframe) then it will be converted to numeric)
#class - vector of classifications for your samples, with 0 representing your negative class (controls) and 1 representing your positive class (cases)
#logged2 - set as TRUE if the gene expression values are log2 normalized
#adjust - if TRUE, then the fold change value will be adjusted to be on a linear scale (i.e fc 1.5 becomes 0.5, fc 0.2 becomes -4)

#RETURN VALUE
#returns the fold change value

#REQUIRED PACKAGES: None

.getFC <- function(value, class, logged2=TRUE, adjust=TRUE){
  if(class(value)=="data.frame"){value=as.numeric(value)}
  stopifnot(identical(length(value),length(class)))

  na.vals = which(is.na(value))
  if(length(na.vals)>=1){
    value = value[-na.vals]
    class = class[-na.vals]
  }

  x <- value[which(class==1)]
  y <- value[which(class==0)]

  n1 <- length(x); n2 <- length(y)
  if( n1 < 1 | n2 < 1 ){return(NA)}

  m1   <- mean(x); m2 <- mean(y)
  if(logged2){
    fc = 2^(m1 - m2)
  }else{
    fc = m1/m2
  }
  if(adjust){
    if(fc>=1){
      fc = fc-1
    }else{
      fc = -1*(1/fc-1)
    }
  }
  return(fc)
}



###-###-###-###-###-###
###   .getSAMs0()   ###
###-###-###-###-###-###

#DESCRIPTION
#Gets the fudge factor (s0) for SAM

#PARAMETERS
#s0.perc - if the user desires, they can specify what percentile of the s_i values they want s0 to be (see SAM technical details)
#sd - vector of standard deviation values for each gene (from output of ttest.func)

#RETURN VALUE
#returns the SAM fudge factor (s0)

#REQUIRED PACKAGES: samr

.getSAMs0 <- function(s0.perc=NULL, sd){
  if(is.null(s0.perc) & nrow(x) < 500){
    s0 = quantile(sd, 0.05)
    s0.perc = 0.05
  }else if(!is.null(s0.perc)) {
    if ((s0.perc != -1 & s0.perc < 0) | s0.perc >
        100) {
      stop("Illegal value for s0.perc: must be between 0 and 100, or equal\nto (-1) (meaning that s0 should be set to zero)")
    }
    if (s0.perc == -1) {
      s0 = 0
    }
    if (s0.perc >= 0) {
      s0 <- quantile(init.fit$sd, s0.perc/100)
    }
  }
  else if (is.null(s0.perc)) {
    s0 = samr:::est.s0(init.fit$tt, init.fit$sd)$s0.hat
    s0.perc = 100 * sum(init.fit$sd < s0)/length(init.fit$sd)
  }
  return(data.frame(s0=s0,s0.perc=s0.perc))
}



###-###-###-###-###-###-###-#
###   .poolFoldChange()   ###
###-###-###-###-###-###-###-#

#DESCRIPTION
#For pooling fold changes that are not adjusted to a linear scale, since you can't just take a simple mean

#PARAMETERS
#fc.vals - vector of fold change values
#weights - vector of weights, if desired

#RETURN VALUE
#returns the pooled fold change value on the original, unadjusted scale

#REQUIRED PACKAGES: None

.poolFoldChange <- function(fc.vals,weights=NULL){
  fc.vec = rep(0,length(fc.vals))
  for(i in 1:length(fc.vals)){ #put fold change values on a linear scale
    if(fc.vals[i]>=1){
      fc.vec[i] = fc.vals[i]-1
    }else{
      fc.vec[i] = -1*(1/fc.vals[i]-1)
    }
  }
  if(is.null(weights)){
    fc.pool = mean(fc.vec)
  }else{
    fc.pool = sum(fc.vec*weights)
  }
  if(fc.pool>=0){ #convert back
    fc.pool = fc.pool+1
  }else{
    fc.pool = 1/((-1*fc.pool)+1)
  }
  return(fc.pool)
}



###-###-###-###-###-##-#
###   .getAUC_CI()   ###
###-###-###-###-###-##-#

#DESCRIPTION
#Get the confidence interval for the AUC

#PARAMETERS
#auc - AUC value
#n.pos - number of positive classifications (i.e. "1")
#n.pos - number of negative classifications (i.e. "0")

#RETURN VALUE
#returns a vector with the lower and then upper confidence intervals

#REQUIRED PACKAGES: None

.getAUC_CI <- function(auc, n.pos, n.neg){
  q1 <- auc/(2 - auc)
  q2 <- (2 * auc^2)/(1 + auc)
  se.auc <- sqrt((((auc * (1 - auc)) + (n.pos - 1) * (q1 - auc^2)) + ((n.neg - 1) * (q2 - auc^2)))/(n.pos * n.neg))
  ci.upper <- auc + (se.auc * 1.96)
  ci.lower <- auc - (se.auc * 1.96)
  return(c(ci.lower,ci.upper))
}



#### Debugged COCONUT functions####

#added drop=F to a bunch of stuff
#this is to fix this particular error: Error in X$genes[common, ] : incorrect number of dimensions
debugCOCONUT <- function(GSEs, control.0.col, disease.col = NULL, byPlatform = FALSE, platformCol,
                         par.prior = TRUE, itConv = 1e-04, parallel = FALSE, mc.cores = 1){
  common <- Reduce(intersect, lapply(GSEs, function(X) rownames(X$genes)))
  common <- common[!(common %in% c(NA, ""))]
  GSEs.control <- lapply(GSEs, function(x){
    x$pheno <- x$pheno[(x$pheno[, control.0.col] %in% 0), , drop=F]
    x$genes <- x$genes[, rownames(x$pheno), drop=F]
    x$class <- rep(0, ncol(x$genes))
    x
  })
  if(byPlatform){
    checkPlatforms <- function(GSElist){
      platforms <- lapply(GSElist, function(x){
        x$pheno[, grep(platformCol, colnames(x$pheno)), drop=F][1]
      })
      sort(unique(unlist(platforms)))
    }
    if(!identical(checkPlatforms(GSEs), checkPlatforms(GSEs.control))){
      stop("byPlatform = T but not all platforms have associated controls")
    }
  }else{
    check <- unlist(lapply(GSEs.control, function(x) length(x$class) > 1))
    if (!all(check)) {
      stop(paste("Datasets with <1 control:", names(GSEs.control[!check])))
    }
  }
  ComBatcontrol <- debug.CombatCustom(GSEs.control, common, params = "get", 
                                      byPlatform = byPlatform, platformCol = platformCol, par.prior = par.prior, 
                                      itConv = itConv, parallel = parallel, mc.cores = mc.cores)
  if(is.null(disease.col)){
    GSEs.disease <- lapply(GSEs, function(x){
      x$pheno <- x$pheno[!(x$pheno[, control.0.col] %in% 0), , drop=F]
      x$genes <- x$genes[, rownames(x$pheno), drop=F]
      x$class <- x$pheno[, control.0.col, drop=F]
      x
    })
  }else{
    GSEs.disease <- lapply(GSEs, function(x){
      x$pheno <- x$pheno[!is.na(x$pheno[, disease.col]), , drop=F]
      x$genes <- x$genes[, rownames(x$pheno), drop=F]
      x$class <- x$pheno[, disease.col, drop=F]
      x
    })
  }
  GSEs.disease.ComBat <- debug.CombatCustom(GSEs.disease, common, 
                                            params = "have", byPlatform = byPlatform, platformCol = platformCol, 
                                            bayesParams = ComBatcontrol$bayesParams, getPlatforms = ComBatcontrol$getPlatforms)
  return(list(COCONUTList = GSEs.disease.ComBat, rawDiseaseList = GSEs.disease, controlList = ComBatcontrol))
}



debug.CombatCustom <- function(GSE.list.genes, common, params, byPlatform, platformCol, bayesParams,  
                               getPlatforms = NULL, par.prior = TRUE, itConv = 1e-04, parallel = FALSE, mc.cores = 1){
  stopifnot(params %in% c("get", "have"))
  commonGenes <- Reduce(cbind, lapply(GSE.list.genes, function(X) X$genes[common, , drop=F]))
  if(byPlatform){
    cat("\nTreating platforms as batches...\nPlatforms found: ")
    platforms <- factor(unlist(lapply(GSE.list.genes, function(X){
      platform <- grep(platformCol, colnames(X$pheno), ignore.case = T)[1]
      GSEplatform <- as.character(X$pheno[, platform])
      cat(GSEplatform[1], ", ")
      GSEplatform
    })))
    index <- order(platforms)
    platforms <- platforms[index]
    commonGenes <- commonGenes[index]
    requireNamespace("limma")
    cat("\nCo-quantile-normalizing datasets from the same platform ")
    invisible(lapply(levels(platforms), function(platform){
      index <- colnames(commonGenes)[platforms == platform]
      samePlatformData <- commonGenes[, index]
      samePlatformData <- limma::normalizeQuantiles(samePlatformData)
      commonGenes[, index] <<- samePlatformData
      cat(" .")
    }))
  }else {
    cat("\nTreating datasets as batches: ")
    platforms <- unlist(lapply(1:length(GSE.list.genes), 
                               function(i){
                                 rep(i, ncol(GSE.list.genes[[i]]$genes))
                               }))
  }
  if(params == "get"){
    ComBatWithParams <- COCONUT:::.ComBatGetParamsNoCov(commonGenes, 
                                                        batch = platforms, par.prior = par.prior, itConv = itConv, 
                                                        parallel = parallel, mc.cores = mc.cores)
    ComBatWithParams$GSEs <- lapply(GSE.list.genes, function(GSE){
      GSE$genes <- data.frame(ComBatWithParams$bayesdata[, colnames(GSE$genes), drop=F])
      return(GSE)
    })
    ComBatWithParams$bayesdata <- NULL
    ComBatWithParams$getPlatforms <- unique(platforms)
    return(ComBatWithParams)
  }else if(params == "have"){
    if(!(identical(unique(platforms), getPlatforms))){
      print(unique(platforms))
      print(getPlatforms)
      stop("Batches not in identical order between controls and cases.")
    }else{
      cat("\nBatches identical between have-params and get-params...\n")
    }
    commonGenesCombat <- COCONUT:::.ComBatApplyParamsNoCov(commonGenes, batch = platforms, bayesParams = bayesParams)
    GSE.list.genes <- lapply(GSE.list.genes, function(GSE) {
      GSE$genes <- data.frame(commonGenesCombat[, colnames(GSE$genes), drop=F])
      return(GSE)
    })
    return(GSE.list.genes)
  }
}



####Functions without descriptions written####
##
##
##
##FIGURE OUT HOW THESE FUNCTIONs WORK AND WRITE A DESCRIPTION FOR THEM
##
##
##

#plotting function that time wrote for getting the y-limit
.getYLim <- function(vec) {
  cushion <- 0.05
  tmp <- quantile(vec, c(0.01, 0.99), na.rm=T)
  tmp <- c(tmp[1]-cushion*diff(tmp), tmp[2]+cushion*diff(tmp)*2)
  tmp
}


#plotting function that time wrote for getting plot colors
.getDataPal <- function(datanames){
  dataPal <- col2rgb(c("red3", "darkorange", "yellow", "green4", "blue1", "cyan",  "darkmagenta"), alpha=T )
  dataPal["alpha", ] <- 150
  dataPal <- apply(dataPal, 2, function(pal) rgb(pal[1], pal[2], pal[3], pal[4], maxColorValue = 255))
  dataPal <- colorRampPalette(dataPal, alpha=T)(length(datanames))
  dataPal
}


#modified imputeSex function to work better with my dataset structure
imputeSex.genes <- function(myDataset, femGenes = NULL, malGenes = NULL){
  # If femGenes not given, use known X-escapees from iSEXS as default
  if(is.null(femGenes)){
    femGenes = c("XIST","RPS4X","CD40LG","ZRSR2",
                 "EFHC2","CA5B","ZFX","EIF1AX",
                 "CA5BP1","UBA1","SYAP1","DDX3X",
                 "FUNDC1","USP9X","SMC1A","NUP62CL","NAA10")
  }
  
  # If malGenes not give, use Y-chromosome genes from iSEXS from iSEXS as default
  if(is.null(malGenes)){
    malGenes = c("KDM5D","RPS4Y1","EIF1AY","USP9Y",
                 "DDX3Y","UTY","PRKY","ZFY","TMSB4Y")
  }
  
  if (!any(c(femGenes, malGenes) %in% rownames(myDataset$genes))) {
    stop("Sex classification genes not present in dataset")
  }
  femGenes = femGenes[which(femGenes %in% rownames(myDataset$genes))]
  malGenes = malGenes[which(malGenes %in% rownames(myDataset$genes))]
  sexGenes = c(femGenes, malGenes)
  sexExpr = myDataset$genes[sexGenes,,drop=F]
  kmeansResults = stats::kmeans(t(sexExpr), centers = 2)
  if (length(femGenes) > 0 & length(malGenes) > 0) {
    femCenter = kmeansResults$centers[, femGenes]
    if (!is.null(ncol(femCenter))) {
      femCenter = rowMeans(femCenter)
    }
    malCenter = kmeansResults$centers[, malGenes]
    if (!is.null(ncol(malCenter))) {
      malCenter = rowMeans(malCenter)
    }
    highFemCentroid = names(femCenter)[which(femCenter ==
                                               max(femCenter))]
    lowMalCentroid = names(malCenter)[which(malCenter ==
                                              min(malCenter))]
    if (highFemCentroid != lowMalCentroid) {
      stop("Strange centroid assignment in kmeans. Try different genes")
    }
    femaleCentroid = highFemCentroid
  }
  if (length(femGenes) == 0) {
    malCenter = kmeansResults$centers[, malGenes]
    if (!is.null(ncol(malCenter))) {
      malCenter = rowMeans(malCenter)
    }
    femaleCentroid = names(malCenter)[which(malCenter ==
                                              min(malCenter))]
  }
  if (length(malGenes) == 0) {
    femCenter = kmeansResults$centers[, femGenes]
    if (!is.null(ncol(femCenter))) {
      femCenter = rowMeans(femCenter)
    }
    femaleCentroid = names(femCenter)[which(femCenter ==
                                              max(femCenter))]
  }
  imputedSex = ifelse(kmeansResults$cluster == as.numeric(femaleCentroid),
                      yes = "F", no = "M")
  return(imputedSex)
}

unfactor_df <- function(df){
  for(i in 1:ncol(df)){
    if(class(df[,i]) == "factor"){
      df[,i] = as.character(df[,i])
    }
  }
  return(df)
}

emptystring_NA <- function(df){
  for(i in 1:ncol(df)){
    if(any(df[,i] == "",na.rm = T)){
      df[,i][df[,i] == ""] = NA
    }
  }
  return(df)
}

capitalize_str <- function(y) {
  if(class(y) == "factor"){y = as.character(y)}
  my.NAs = which(is.na(y))
  y[my.NAs] = "placeholder"
  if(length(y) == 1){
    c = strsplit(y, " ")[[1]]
    str = paste(toupper(substring(c, 1,1)), substring(c, 2), sep="", collapse=" ")
  }else if(length(y) > 1){
    str = sapply(y, function(x){
      c = strsplit(x, " ")[[1]]
      c = paste(toupper(substring(c, 1,1)), substring(c, 2), sep="", collapse=" ")
      c
    })
    names(str) = NULL
  }
  str[my.NAs] = NA
  return(str)
}

sex_MF <- function(sexVec){
  return(ifelse(grepl("fem|Fem|FEM",sexVec),"F","M"))
}

date_dob <- function(dateVec,endDateMDY=NULL,sep="-",dateMDY=TRUE){
  library(eeptools) #note that eeptools needs the date to be in year/month/day
  endDate = as.Date(endDateMDY,format = "%m-%d-%Y")
  if(dateMDY){
    dateList = strsplit(dateVec,split=sep)
    ageVec = unlist(lapply(dateList, function(date){
      date = as.numeric(date)
      #make year into 4 digits
      if(nchar(date[3]) == 2){
        if(date[3]<19){
          date[3] = 2000+date[3]
        }else{
          date[3] = 1900+date[3]
        }
      }
      startDate = as.Date(paste(date,collapse="-"),format = "%m-%d-%Y")
      return(age_calc(startDate,endDate,units = "years"))
    }))
  }else{
    stop("Still need to add support for any format other than day-month-year")
  }
  return(ageVec)
}

#add I^2 to this
pool.inverseVar <- function(g, se.g, method){
  require(rmeta)
  stopifnot( identical( rownames(g), rownames(se.g) ) )
  out <- matrix( nr=nrow(g), nc=8,
                 dimnames=list( rownames(g), c("n.studies", "summary", "se.summary", "tau2", "p.value", "Q", "df", "pval.het") ) )

  cleanNA <- function(x) return( x[!is.na(x) & is.finite(x) ] )
  for(j in 1:nrow(g)){

    e  <- cleanNA(    g[j, ] )
    se <- cleanNA( se.g[j, ] )
    n  <- length(e)

    if(n==1){
      summ <- e;   se.summ <- se;   tau2 <- NA
      Q.het = NA
      df.het = NA
      pval.het = NA
    } else {
      fit <- meta.summaries(e, se, method = method)
      summ <- fit$summary
      se.summ <- fit$se.summary
      tau2 <- ifelse( method=="fixed", NA, fit$tau2 )
      Q.het = fit$het[1]
      df.het = fit$het[2]
      pval.het = fit$het[3]
      rm(fit)
    }

    pval     <- 2*pnorm( abs(summ/se.summ), lower.tail=FALSE )

    out[j, ] <- c(n, summ, se.summ, tau2, pval, Q.het, df.het, pval.het)
    rm(e, se, n, summ, se.summ, tau2, pval, Q.het, df.het, pval.het)
  }
  return(out)
}





combine.effect.sizes <- function (list.of.effects, between.method="random",
                                  within.method="fixed.iv", everything=TRUE,
                                  parallel=F){

  if( is.null(names(list.of.effects)) )
    names(list.of.effects) <- paste("data", 1:length(list.of.effects), sep="")

  studyEffects <- function(effects) {

    effects <- data.frame(effects)
    effects$keys <- as.character(effects$keys)

    ## remove probes that cannot be mapped or have insufficient observations to calculate effect size
    bad <- which( is.na(effects$g) | is.na(effects$keys) | effects$keys=="NA" )
    effects <- effects[ setdiff(1:nrow(effects), bad), ]

    ## expand probes that maps to multiple keys
    effects <- expand.df( effects )

    ## summarize multiple probes within a study
    effects <- summ.eff.within(effects, option = within.method)
    effects
  }

  if(parallel){
    library(parallel)
    study.effects <- mclapply(list.of.effects, mc.cores=length(list.of.effects), studyEffects)
  } else {
    study.effects <- lapply(list.of.effects, studyEffects)
  }

  tmp <- multimerge(study.effects)
  g    <- tmp[, paste(names(study.effects), "_g", sep = ""), drop=FALSE]
  se.g <- tmp[, paste(names(study.effects), "_se.g", sep = ""), drop=FALSE]

  pooled.estimates <- data.frame( pool.inverseVar(g, se.g, method=between.method ) )

  if (everything) {
    return(list(g=g, se.g=se.g, pooled.estimates=pooled.estimates))
  } else {
    return(pooled.estimates)
  }
}
