
##-###-###-###-###-###-###
###   makeBeeswarm()   ###
##-###-###-###-###-###-###

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
#title - title of the figure
#display.auc - whether to display the global AUC as part of the figure title
#group.colors - the colors used for each class (the number of colors must match the number of classes that are plotted. By default this will be the number of classes in pooledDataObject$class, but that will change if the classLabels parameter is utilized.)
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

#RETURN VALUE
#A violin plot comparing the signature scores for each class of samples

#REQUIRED PACKAGES: OptimalCutpoints, ggbeeswarm

#group.colname - this designates the column name of the group vector within the pooledDataObject$pheno dataframe

#classLabels - if you want to use something other pooledDataObject$class (for example, if you want to highlight three classes instead of just cases vs. controls) you can provide an alternate class vector

# pooledDataObject = ucla_train_copd
# filterObject=NULL
# upgenes = up
# downgenes = down
# group.colname = "copd"
# title="UCLA Training"
# group.min = 3
# group.colors = c("firebrick3","#0072B2")
# dotsize=4
# width=0.4
# bandwidth=0.5
# show.legend=T
# alpha=1
# beeswarm.cex=2
# allText.size = NULL
# title.size=18
# xlab.size=12
# xlab.angle=45
# axistitle.size=14
# ylab.size=12
# legend.size=10
# classLabels=NULL
# show.violin = F
# show.meanSD = F
# show.IQR = T
# show.ttest = F
# show.wilcox = F
# IQR.top=F
# type="quasirandom"
# splitGroup=FALSE
# zscore=FALSE
# score.overwrite=NULL

  
makeBeeswarm <- function(pooledDataObject, filterObject=NULL, upgenes=NULL, downgenes=NULL, type="quasirandom", group.colname = "group", names=NULL, title=NULL, display.auc=T,
                         group.colors = c("firebrick3","#0072B2"), group.min=NULL, classLabels=NULL, show.violin = F, show.meanSD = F, show.IQR = T, show.ttest = F, show.wilcox = F, IQR.top=T, dotsize=3, width=0.4,
                         beeswarm.cex=2, bandwidth=0.5, alpha=0.9, IQR.width=0.03, show.legend=TRUE, allText.size = NULL, title.size=18, axistitle.size=14, legend.size=10, xlab.size=12, 
                         xlab.angle=0, ylab.size=10, splitGroup=FALSE, zscore=FALSE, score.overwrite=NULL, ...){
  require(ggpubr)
  
  #internal variables
  rounding=3
  
  #cutoff.method - what optimal cutoff method to use. Choose from Youden, MinValueSe, MinValueSp, or most other options within OptimalCutpoints (some that don't work can probably be added)
  cutoff.method="Youden"
  
  #internal methods
  gg_mean_sd <- function(x) {
    m <- mean(x)
    ymin <- m-sd(x)
    ymax <- m+sd(x)
    return(c(y=m,ymin=ymin,ymax=ymax))
  }
  
  #checking
  if(!checkManateeObject(pooledDataObject,"pooledDataObject")){stop("Invalid pooledDataObject provided")}
  
  if(show.meanSD && show.IQR){stop("mean+SD and IQR cannot both be shown")}
  if(show.ttest && show.wilcox){stop("t-test and wilcoxon cannot both be shown")}
  
  if(is.null(classLabels)){
    class = pooledDataObject$class
  }else{
    class = classLabels
  }
  if(length(group.colors) != length(unique(class))){
    stop("The number of colors in group.colors must match the number of unique classes in pooledDataObject$class")
  }
  
  type = tolower(type)
  if(type %in% c("quasirandom","quasi","quasirandom plot","quasi plot")){
    type = "quasirandom"
  }else if(type %in% c("dot","dot plot","dotplot","stacked dot","stacked dot plot","stacked dotplot")){
    type = "dot"
  }else if(type %in% c("beeswarm","beeswarm plot","bee","bee plot")){
    type = "beeswarm"
  }else if(type %in% c("pseudorandom","pseudorandom plot","pseudo","pseudo plot")){
    type = "pseudorandom"
  }else if(type %in% c("smiley","smiley plot")){
    type = "smiley"
  }else if(type %in% c("frowney","frowney plot","frowny","frowny plot")){
    type = "frowney"
  }else if(type %in% c("tukey","tukey plot")){
    type = "tukey"
  }else{stop("Plot type could not be identified.")}
  #OTHER OPTIONS: stacked dot w/ horizontal jitter
  
  if(is.null(filterObject)){
    if(is.null(upgenes) && is.null(downgenes)){stop("If filterObject is null, then upgenes and downgenes must be provided")}
    if((class(upgenes)=="character" && class(downgenes)=="character") || (is.null(upgenes) && class(downgenes)=="character") ||
       (class(upgenes)=="character" && is.null(downgenes))){
      test.score = getGeneScores(pooledDataObject$genes,upgenes,downgenes)
      num.genes = length(c(upgenes,downgenes))
    }else{
      stop("upgenes and downgenes must both be character vectors")
    }
  }else{
    if(!checkManateeObject(filterObject,"filterObject")){stop("Invalid filterObject provided")}
    test.score = getSigScores(pooledDataObject,filterObject)
    num.genes = length(c(filterObject$upGeneNames,filterObject$downGeneNames))
  }
  
  if(!is.null(score.overwrite)){
    test.score=score.overwrite
  }
  
  if(zscore){
    test.score = scale(test.score)
  }
  
  grouping = pooledDataObject$pheno[,group.colname]
  if(any(is.na(grouping))){grouping[is.na(grouping)] = "NA"} #bugfix
  if(class(grouping) == "numeric"){grouping = as.character(grouping)}
  
  ## get ROC data; this is equivalent to the internal methods
  require(OptimalCutpoints)
  opt <- optimal.cutpoints(X="test.score", status="class", tag.healthy=0,
                           data=data.frame(class, test.score),
                           methods=c("Youden", cutoff.method),
                           control=control.cutpoints(valueSe=0.945))
  AUC <- getAUC(class,test.score,rounding)
  cat(sprintf("AUC = %s (90%% CI %s - %s)\n", AUC$auc, max(0,AUC$auc.CI[1]), min(1,AUC$auc.CI[2])))
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
    if(length(remove.entries) > 0){
      pooledDataObject$genes = pooledDataObject$genes[,-remove.entries]
      class = class[-remove.entries]
      grouping = grouping[-remove.entries]
      test.score = test.score[-remove.entries]
    }
  }
  
  if(display.auc){title = paste0(title, "\n", sprintf(" Global AUC = %s (95%% CI %s - %s)", AUC$auc, max(0,AUC$auc.CI[1]), min(1,AUC$auc.CI[2])))}
  
  #split group by class, if desired
  if(splitGroup){
    if(is.factor(grouping)){
      if(!setequal(levels(grouping),na.omit(unique(grouping)))){
        stop("The levels of grouping must contain all the entries in grouping if you want grouping to be a factor. This error should never occur I think?")
      }
      my.levels = levels(grouping)
      new.levels = c()
      grouping = as.character(grouping)
      for(gname in my.levels){
        if(length(unique(class[grouping == gname])) > 1){
          new.levels = append(new.levels,sort(unique(paste0(grouping[grouping == gname],".",class[grouping == gname])))) #for now i'm just using sort() to make sure the order is the same
          grouping[grouping == gname] = paste0(grouping[grouping == gname],".",class[grouping == gname])
        }else{
          new.levels = append(new.levels,gname)
        }
      }
      grouping = factor(grouping,levels = new.levels)
      if(any(is.na(grouping))){
        grouping[is.na(grouping)] = "NA"
        warning("Some of the entries in grouping are NA - this may be unintended behavior caused by the way this function handles factors, be careful!")
      }
    }else{
      for(gname in unique(grouping)){
        if(length(unique(class[grouping == gname])) > 1){
          grouping[grouping == gname] = paste0(grouping[grouping == gname],".",class[grouping == gname])
        }
      }
    }
  }
  
  #ggbeeswarm version
  plotData = data.frame(scores = test.score,class=class,group=grouping)

  if(type=="quasirandom"){
    bee = ggplot(data=plotData,aes(x=group,y=scores,color=as.character(class)))+
      ggbeeswarm::geom_quasirandom(size=dotsize, width=width, bandwidth=bandwidth, show.legend=show.legend, alpha=alpha, ...)
    #Other options: varwidth
  }else if(type=="dot"){
    bee = ggplot(data=plotData,aes(x=group,y=scores,color=as.character(class),fill=as.character(class)))+
      geom_dotplot(dotsize=(dotsize/5.5), binaxis="y", stackdir="center", show.legend=show.legend, alpha=alpha, ...)+
      scale_fill_manual(values=group.colors)
    #Other options: binpositions="all"
  }else if(type=="beeswarm"){
    bee = ggplot(data=plotData,aes(x=group,y=scores,color=as.character(class)))+
      ggbeeswarm::geom_beeswarm(size=dotsize, cex=beeswarm.cex, show.legend=show.legend, alpha=alpha, ...)
    #Other options: priority
  }else if(type=="pseudorandom"){
    bee = ggplot(data=plotData,aes(x=group,y=scores,color=as.character(class)))+
      ggbeeswarm::geom_quasirandom(size=dotsize, width=width, bandwidth=bandwidth, show.legend=show.legend, alpha=alpha, method = "pseudorandom", ...)
    #Other options: varwidth
  }else if(type=="smiley"){
    bee = ggplot(data=plotData,aes(x=group,y=scores,color=as.character(class)))+
      ggbeeswarm::geom_quasirandom(size=dotsize, width=width, bandwidth=bandwidth, show.legend=show.legend, alpha=alpha, method = "smiley", ...)
    #Other options: varwidth
  }else if(type=="frowney"){
    bee = ggplot(data=plotData,aes(x=group,y=scores,color=as.character(class)))+
      ggbeeswarm::geom_quasirandom(size=dotsize, width=width, bandwidth=bandwidth, show.legend=show.legend, alpha=alpha, method = "frowney", ...)
    #Other options: varwidth
  }else if(type=="tukey"){
    bee = ggplot(data=plotData,aes(x=group,y=scores,color=as.character(class)))+
      ggbeeswarm::geom_quasirandom(size=dotsize, width=width, bandwidth=bandwidth, show.legend=show.legend, alpha=alpha, method = "tukey", ...)
    #Other options: varwidth
  }
  
  if(show.meanSD){
    bee = bee + stat_summary(fun.data=gg_mean_sd,color="black")
  }
  
  if(show.IQR){
    if(IQR.top){
      bee = bee + geom_boxplot(width=IQR.width,color="black",outlier.shape = NA)
    }else{
      #might want to use IQR.width=0.1 if IQR is at the back
      bee$layers = c(geom_boxplot(width=IQR.width,color="black",outlier.shape = NA),bee$layers)
    }
  }
  
  if(show.violin){
    bee$layers = c(geom_violin(trim=F,lwd=1.5),bee$layers)
  }
  
  bee+scale_colour_manual(values=group.colors)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.title = element_blank(),
          axis.line = element_line(colour = "black"),
          plot.title = element_text(size=title.size*ifelse(is.null(allText.size),1,allText.size),hjust = 0.5),
          axis.text.x = element_text(size=xlab.size*ifelse(is.null(allText.size),1,allText.size),angle=xlab.angle,hjust=ifelse(xlab.angle == 0,0.5,1)),
          axis.text.y = element_text(size=ylab.size*ifelse(is.null(allText.size),1,allText.size)),
          axis.title = element_text(size=axistitle.size*ifelse(is.null(allText.size),1,allText.size)),
          legend.text = element_text(size=legend.size*ifelse(is.null(allText.size),1,allText.size)))+
    labs(x=NULL,y=paste0(num.genes," Gene Score"), title=title)+
    {if(display.auc)geom_hline(yintercept=cutoff, linetype="dashed", color = "black")}+
    {if(show.ttest)stat_compare_means(comparisons = combn(unique(grouping),2,simplify=FALSE),method = "t.test")}+
    {if(show.wilcox)stat_compare_means(comparisons = combn(unique(grouping),2,simplify=FALSE),method = "wilcox.test")}

  #parameters that need to be added
  #option to add point outlines
  
  #add option to customize the legend at some point, but it seems annoying so i'll hold off for now
  
  #add violin plot outline and/or gray background
}


