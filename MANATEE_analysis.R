# Bacterial vs Viral proteases for nanoparticle sensors
# Author: Aditya Rao
# Contact: adityamr@stanford.edu
# Script for getting protease targets that distinguish bacterial and viral infections

#NOTE: data are publicly available from NCBI GEO and EMBL-EBI ArrayExpress.
#This script can be used to normalize and analyze the datasets. 

#Source the necessary libraries and functions
library(ROCR)
library(MetaIntegrator)
library(samr)
library(parallel)
library(pbmcapply)
library(data.table)
library(COCONUT)
sapply(list.files(path="MANATEE_functions/", pattern="*.R",full.names=T), source)

#Split data into Discovery/Hold-out Validation and separately COCONUT conormalize each
#BV.GSEs is comprised of publicly available datasets from NCBI GEO and EMBL-EBI ArrayExpress
BV.discovalid = makeDiscoValid(BV.GSEs, use.pheno.class = TRUE, remove.samples=c("Healthy","Remove","C. Albicans"),
                               seed=36362, conorm.type = "COCONUT",forced.discovery = "GSE16129GPL96")

#format the Discovery and Hold-out Validation data
BV.discovalid$Discovery$pheno$group = droplevels(BV.discovalid$Discovery$pheno$group)
BV.discovalid$Validation$pheno$group = droplevels(BV.discovalid$Validation$pheno$group)
BV.discovalid$Discovery$pheno$group2 = BV.discovalid$Discovery$pheno$group
BV.discovalid$Validation$pheno$group2 = BV.discovalid$Validation$pheno$group
icbac.names = c("S. Typhi","Brucella","Salmonella","L. Monocytogenes","B. Pseudomallei","S. Enterica","O. Tsutsugamushi","R. Typhi","Rickettsia",
                "C. Burnetii","Murine typhus","Scrub typhus","C. Trachomatis","T. Whipplei")
ecbac.names = c("E. Coli","MRSA","MSSA","S. Pneumoniae","S. Pyogenes","S. Agalactiae","S. Aureus","N. Meningitidis","P. Aeruginosa",
                "H. Influenzae","K. Pneumoniae","Enterobacter","M. Catarrhalis","Enterococcus","Acinetobacter","C. Difficile","S. Marcescens",
                "Sphingomonas","Micrococcus","Viridans Streptococci","Corynebacterium","Aeromonas","A. Lwoffii","A. Baumannii","S. Suis",
                "S. Epidermidis","CoNS","E. Faecium","E. Faecalis","A. Hydrophila","S. Anginosus","C. Freundii","E. Cloacae","K. Oxytoca",
                "Neisseria","E. Corrodens","N. Gonorrhoeae","M. Genitalium","Staphylococcus","Pseudomonas","Leptospira","L. Interrogans")
idkbac.names = c("Bacterial","Gram Positive Bacteria","Gram Negative Bacteria")
bac.names = c(icbac.names,ecbac.names,idkbac.names)
viral.names = c("DF","DHF","DSS","Adenovirus","HHV6","Enterovirus","Influenza","RSV","Enterovirus/Rhinovirus","Rhinovirus","Influenza/Rhinovirus",
                "Coronavirus","HSV","Parainfluenza","HMPV","BKV","CMV","Dengue","Varicella","Rotavirus","HIV","Measles","HCMV","EBV")

BV.discovalid$Discovery$pheno$group = BV.discovalid$Discovery$pheno$group2
BV.discovalid$Validation$pheno$group = BV.discovalid$Validation$pheno$group2
levels(BV.discovalid$Discovery$pheno$group) = c(levels(BV.discovalid$Discovery$pheno$group),"Bacterial","Viral")
BV.discovalid$Discovery$pheno$group[BV.discovalid$Discovery$pheno$group %in% bac.names]="Bacterial"
BV.discovalid$Discovery$pheno$group[BV.discovalid$Discovery$pheno$group %in% viral.names]="Viral"
BV.discovalid$Discovery$pheno$group = droplevels(BV.discovalid$Discovery$pheno$group)
levels(BV.discovalid$Validation$pheno$group) = c(levels(BV.discovalid$Validation$pheno$group),"Bacterial","Viral")
BV.discovalid$Validation$pheno$group[BV.discovalid$Validation$pheno$group %in% bac.names]="Bacterial"
BV.discovalid$Validation$pheno$group[BV.discovalid$Validation$pheno$group %in% viral.names]="Viral"
BV.discovalid$Validation$pheno$group = droplevels(BV.discovalid$Validation$pheno$group)

BV.discovalid$Discovery$class = makeClassVector(BV.discovalid$Discovery$pheno$group,"Bacterial")
BV.discovalid$Validation$class = makeClassVector(BV.discovalid$Validation$pheno$group,"Bacterial")

#Load the filtered list of proteases
protease_genes = readRDS("data/protease_genes.rds")

#filter the Discovery and Hold-out Validation to have only the filtered protease genes
BV.discovalid$Discovery$genes = BV.discovalid$Discovery$genes[protease_genes[protease_genes %in% rownames(BV.discovalid$Discovery$genes) & protease_genes %in% rownames(WB.COCO$genes)],]
BV.discovalid$Validation$genes = BV.discovalid$Validation$genes[protease_genes[protease_genes %in% rownames(BV.discovalid$Validation$genes) & protease_genes %in% rownames(WB.COCO$genes)],]

#run Basic Manatee with LOSO
BV.Basic = runManatee(BV.discovalid$Discovery, manatee.type = "Basic", seed=36362, runLeaveOneOutAnalysis = TRUE, numCores = 1)

#filter Manatee with effect size threshold of 0.6 and FDR threshold of 0.01
BV.Basic = filterManatee(BV.Basic, EffectSizeThresh = 0.6, ttestFDR = 0.01,isLeaveOneOut = T)

#test signature performance in Discovery data
makeROCplot(BV.Basic,BV.Basic$filterResults$loo_ES0.6_tFDR0.01_tgp_cf,
            title="Bacterial vs. Viral - Discovery")

#Test the signature on the Hold-out Validation data
makeROCplot(BV.discovalid$Validation,BV.Basic$filterResults$loo_ES0.6_tFDR0.01_tgp_cf,
            title="Bacterial vs. Viral - Holdout Validation")

#COCONUT conormalize the Independent Validation datasets
#IV.GSEs is comprised of publicly available datasets from NCBI GEO and EMBL-EBI ArrayExpress
IV.COCO.out = COCONUT(GSEs = IV.GSEs, control.0.col = "control.0.class")
IV.COCO = combineCOCOoutput(IV.COCO.out)

#format the Independent Validation data
remove.sample.names = c("Healthy","Remove","C. Albicans")
IV.COCO = removeSamples(IV.COCO, which(IV.COCO$pheno$group %in% remove.sample.names),expr=FALSE, class=FALSE)
IV.COCO$pheno$group = droplevels(IV.COCO$pheno$group)
IV.COCO$pheno$group2 = IV.COCO$pheno$group
icbac.names = c("S. Typhi","Brucella","Salmonella","L. Monocytogenes","B. Pseudomallei","S. Enterica","O. Tsutsugamushi","R. Typhi","Rickettsia",
                "C. Burnetii","Murine typhus","Scrub typhus","C. Trachomatis","T. Whipplei","Leptospira","L. Interrogans")
ecbac.names = c("E. Coli","MRSA","MSSA","S. Pneumoniae","S. Pyogenes","S. Agalactiae","S. Aureus","N. Meningitidis","P. Aeruginosa",
                "H. Influenzae","K. Pneumoniae","Enterobacter","M. Catarrhalis","Enterococcus","Acinetobacter","C. Difficile","S. Marcescens",
                "Sphingomonas","Micrococcus","Viridans Streptococci","Corynebacterium","Aeromonas","A. Lwoffii","A. Baumannii","S. Suis",
                "S. Epidermidis","CoNS","E. Faecium","E. Faecalis","A. Hydrophila","S. Anginosus","C. Freundii","E. Cloacae","K. Oxytoca",
                "Neisseria","E. Corrodens","N. Gonorrhoeae","M. Genitalium","Staphylococcus","Pseudomonas")
idkbac.names = c("Bacterial","Gram Positive Bacteria","Gram Negative Bacteria")
bac.names = c(icbac.names,ecbac.names,idkbac.names)
viral.names = c("DF","DHF","DSS","Adenovirus","HHV6","Enterovirus","Influenza","RSV","Enterovirus/Rhinovirus","Rhinovirus","Influenza/Rhinovirus",
                "Coronavirus","HSV","Parainfluenza","HMPV","BKV","CMV","Dengue","Varicella","Rotavirus","HIV","Measles","HCMV","EBV")

IV.COCO$pheno$group = IV.COCO$pheno$group2
levels(IV.COCO$pheno$group) = c(levels(IV.COCO$pheno$group),"Bacterial","Viral")
IV.COCO$pheno$group[IV.COCO$pheno$group %in% bac.names]="Bacterial"
IV.COCO$pheno$group[IV.COCO$pheno$group %in% viral.names]="Viral"
IV.COCO$pheno$group = droplevels(IV.COCO$pheno$group)

IV.COCO$class = makeClassVector(IV.COCO$pheno$group,"Bacterial")

#Test the signature on the Independent Validation data
IndependentValidation = readRDS(file="data/IndependentValidation_data.rds")
makeROCplot(IV.COCO,BV.Basic$filterResults$loo_ES0.6_tFDR0.01_tgp_cf,
            title="Bacterial vs. Viral - Independent Validation")
