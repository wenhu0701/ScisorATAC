library(tidyr)
library(dplyr)
library(data.table)
library(parallel)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(RColorBrewer)
library(patchwork)
library(ggplot2)

workingDir <- '../2024_02_29_multiVelo/'
dataDir <- '../2024_02_29_multiVelo/data/'

print("Reading in state files")
stateFile_case <- fread(file.path(dataDir,'AD_vs_Control/multiVeloMatrixFiles/All_9Cases_multivelo_dPSItestedGenes_Cell_States_mtx.csv'))
stateFile_ctrl <- fread(file.path(dataDir,'AD_vs_Control/multiVeloMatrixFiles/All_10Controls_multivelo_dPSItestedGenes_Cell_States_mtx.csv'))
colnames(stateFile_case)[1] <- colnames(stateFile_ctrl)[1] <- "Sample_BC"

print("Reading in all info")
allInfo_case <- fread(file.path(dataDir,'AD_vs_Control/AllInfoFiles/AllInfo_Cases_Corrected_Incomplete_Celltype.gz'))
allInfo_case <- allInfo_case[,c(1:5,9)]

allInfo_ctrl <- fread(file.path(dataDir,'AD_vs_Control/AllInfoFiles/AllInfo_Incomplete_Corrected_Controls_Celltypes.gz'))
allInfo_ctrl <- allInfo_ctrl[,c(1:5,9)]
colnames(allInfo_case) <- colnames(allInfo_ctrl) <- c("Read","Gene","CT","BC","UMI","ExonChain")

print("Reading in extra info")
gene2clear <- read.table('/athena/tilgnerlab/scratch/anj2026/ENSEMBLE_v34_HumanClearGeneNames', header = F)
colnames(gene2clear) <- c("ENS","Clear")

exonList <- read.table(file.path(dataDir,'AD_vs_Control/testedExons_caseVsControls_wClear'), header = F)
colnames(exonList) <- c("Exon","ENS","dPI","pvalue","cellType","Clear")
head(exonList)

#!/bin/R

# By Anoushka Joglekar 
# 03.2024

## Setup ------


## Functions -----

## Subsets out the state vector and known alternative exons for a gene of interest.

getDataPerGene <- function(stateFile,exonList,goi){
    geneID <- gene2clear %>% filter(Clear == goi) %>% .$ENS
    
    goiStates <- stateFile %>% select(Sample_BC,all_of(goi)) %>% 
        separate(Sample_BC, into = c("Sample","Barc"),sep = "_") %>% 
        separate(Barc,into=c("BC","Suffix"),sep = "-")
    
    subsetExons <- exonList %>% filter(Clear == goi) %>% .$Exon
    
    return(list(subsetExons,goiStates))
}

# ## Check that the exon is alternative and the contingency table passes the chi-sq criterion

checkChiSqCrit <- function(mat){
    if (sum(mat) > 0 & max(colSums(mat))/sum(mat) <= 0.9 & 
        min(rowSums(mat))*min(colSums(mat))/sum(mat) > 5){
        return("pass")
    } else {return("fail")}
}

logOddsCalc <- function(mat){## Set 0's to 0.5 and explicitly calculate log odds ratio for 2x2 matrices
    if(nrow(mat)==2){
        mat[mat==0]=0.5
        lor = log2((mat[1,1]*mat[2,2])/(mat[2,1]*mat[1,2]))
        return(lor)
    } else {
        return("NA")}
}


reportStats <- function(mat){
    p = chisq.test(mat)$p.value
    lor = logOddsCalc(mat)
    return(list(p,lor))
}


## Function to calculate psi per state in each matrix

calcPSI = function(mat){
    psiVec = mat[,1]/rowSums(mat)
    return(psiVec)
}

## Given an exon, calculate inclusion and exclusion counts per read. Test for association with accessibility state

perExon <- function(allInfo_mod,subset_states,eoi,goi){
    se <- as.integer(strsplit(eoi,split = "_")[[1]][c(2,3)])

    exonStatus <- allInfo_mod %>% 
        mutate(status = case_when(eoi == ExonChain ~ 0, start <= se[1] & end >= se[2] ~ 1,TRUE ~ 2)) %>% 
        filter(status != 2) %>%
        select(Read,status,BC) %>% group_by(Read) %>% 
        arrange(status,.by_group = T) %>% slice_head(n = 1)

    exonMat <- as.matrix(table(inner_join(exonStatus,subset_states, by = "BC") %>% 
                     ungroup() %>% select(all_of(goi),status)))
    if (nrow(exonMat) >=2 & ncol(exonMat) >=2 & length(which(rowSums(exonMat) >= 5)) >=2){
        if(nrow(exonMat) ==2 & checkChiSqCrit(exonMat) == "pass"){
            s = reportStats(exonMat)
        } else {s = list("x","x")}
        psiVals <- paste(calcPSI(exonMat),collapse = ", ")
        stateNames <-  paste(rownames(exonMat), collapse = ", ")
        totCounts <- sum(exonMat)
        res <- c(eoi,psiVals,stateNames,totCounts,s[[1]],s[[2]])
        return(res)
    }
}
## Go through all alt exons for a gene and cell type and report p-values where possible

performAnalysisForCT <- function(ss,allInfo,geneID,goi,ct){
    allInfo_subset <- allInfo %>% filter(Gene == geneID & CT == ct)
    if (nrow(allInfo_subset) > 10){
        allInfo_mod <- allInfo_subset %>% 
            rowwise() %>%
            mutate(start = unlist(strsplit(ExonChain,"_"))[2],end = rev(unlist(strsplit(ExonChain,"_")))[2]) %>%
            mutate_at(c("start","end"),as.integer) %>%
            separate_rows(ExonChain,sep = ";%;") %>% filter(ExonChain != "")
    
#         subsetExons <- ss[[1]]
        subsetExons <- unique(allInfo_mod %>% .$ExonChain)
        subset_states <- ss[[2]] %>% filter(BC %in% allInfo_mod$BC)

        exonPVals = lapply(subsetExons, function(eoi){
            pe = perExon(allInfo_mod,subset_states,eoi,goi)
            if(length(pe) == 6){
                exonCTList = c(ct,pe)
                return(exonCTList)}})
        
        allExonPVals = do.call('rbind',exonPVals)
        return(allExonPVals)
    }
}

## Function to run analysis for all CTs for a particular gene

analysisPerGene <- function(goi, stateFile, allInfo, allCTs){
    geneID <- gene2clear %>% filter(Clear == goi) %>% .$ENS
    ss = getDataPerGene(stateFile,exonList,goi)
    allRes = do.call('rbind',lapply(allCTs, function(ct) performAnalysisForCT(ss,allInfo,geneID,goi,ct)))
    allRes <- as.data.frame(allRes)
    if(nrow(allRes) > 1){
        colnames(allRes) <- c("CT","Exon","statePSI","stateID","totalCounts","pvalue","LOR")
        allRes$GeneName <- goi
        allRes$GeneID <- geneID
        return(allRes)
    } else{return()}
}

## Just a function to extract counts once you already know the gene, exon, and CT of interest

getMatrix <- function(goi,eoi,ct){
    
    geneID <- gene2clear %>% filter(Clear == goi) %>% .$ENS

    ss = getDataPerGene(stateFile,exonList,goi)
    
    allInfo_mod <- allInfo %>% filter(Gene == geneID & CT == ct) %>% 
        rowwise() %>% 
        mutate(start = unlist(strsplit(ExonChain,"_"))[2],end = rev(unlist(strsplit(ExonChain,"_")))[2]) %>%
        mutate_at(c("start","end"),as.integer) %>%
        separate_rows(ExonChain,sep = ";%;") %>% filter(ExonChain != "")
    
    subsetExons <- ss[[1]]
    subset_states <- ss[[2]] %>% filter(BC %in% allInfo_mod$BC)
    
    se <- as.integer(strsplit(eoi,split = "_")[[1]][c(2,3)])

    exonStatus <- allInfo_mod %>% 
        mutate(status = case_when(eoi == ExonChain ~ 0, start <= se[1] & end >= se[2] ~ 1,TRUE ~ 2)) %>% 
        filter(status != 2) %>%
        select(Read,status,BC) %>% group_by(Read) %>% 
        arrange(status,.by_group = T) %>% slice_head(n = 1)

    exonMat <- as.matrix(table(inner_join(exonStatus,subset_states, by = "BC") %>% 
                     ungroup() %>% select(all_of(goi),status)))
    psiVals <- paste(calcPSI(exonMat),collapse = ", ")
    
    print(exonMat)
}

## Function to identify genes of interest

getGeneList <- function(stateFile){
    allGenes <- colnames(stateFile)[-1]
    print(length(allGenes))                       

#     filteredGenes <- intersect(unique(exonList$Clear),allGenes)
#     print(length(filteredGenes))

    multiStateGenes <- c()

    for(gene in allGenes){
        if(length(table(stateFile[[gene]])) >= 2){
            multiStateGenes <- c(multiStateGenes,gene)
        }
    }

    print(length(multiStateGenes))
    
    print(paste("Percent multistate",length(multiStateGenes)*100/length(allGenes)))

    #print(paste("Total genes from exonList:",length(unique(exonList$Clear))))
    return(multiStateGenes)
}
                                    
                                    
runForSample <- function(allInfo, stateFile){
    
    allCTs <- grep(pattern = "Unknown",unique(allInfo$CT),invert = T,value = T)
    multiStateGenes <- getGeneList(stateFile)
    
    multiStateRes <- as.data.frame(do.call('rbind',mclapply(multiStateGenes, function(gene) {
    print(gene)
    analysisPerGene(gene, stateFile,allInfo, allCTs)},mc.cores = 8)))
    
    multiStateRes <- multiStateRes %>% distinct()
    return(multiStateRes)
}


caseRes <- runForSample(allInfo_case, stateFile_case)
controlRes <- runForSample(allInfo_ctrl, stateFile_ctrl)



caseRes$Sample <- "Case"
controlRes$Sample <- "Control"

caseRes <- caseRes %>% separate_rows(statePSI, stateID, sep = ",\\s*") %>%
mutate_at("statePSI",as.double) %>% mutate_at("stateID", as.factor) %>% distinct()

controlRes <- controlRes %>% separate_rows(statePSI, stateID, sep = ",\\s*") %>%
mutate_at("statePSI",as.double) %>% mutate_at("stateID", as.factor) %>% distinct()

fullMat <- full_join(caseRes,controlRes, by = c("CT","Exon","stateID","GeneName","GeneID")) %>% 
    mutate(statePSI.x = ifelse(is.na(statePSI.x), -1, statePSI.x),
           Sample.x = ifelse(is.na(Sample.x), "Case", Sample.x),
           statePSI.y = ifelse(is.na(statePSI.y), -1, statePSI.y),
          Sample.y = ifelse(is.na(Sample.y), "Control", Sample.y),)
head(fullMat,5)



dir.create('../2024_02_29_multiVelo/output/v2_ADvsControl_newSamples/')
write.table(fullMat,'../2024_02_29_multiVelo/output/v2_ADvsControl_newSamples/AD_vs_control_fullMat_altWorkflow.tab',sep = "\t",
           quote = F, row.name = F, col.names = T)
           
           
           