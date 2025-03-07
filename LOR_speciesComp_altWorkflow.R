library(tidyr)
library(dplyr)
library(data.table)
library(parallel)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(ggplot2)

workingDir <- '../2024_02_29_multiVelo/'
dataDir <- '../2024_02_29_multiVelo/data/'

print("Reading in state files")
stateFile_case <- fread(file.path(dataDir,'Human_vs_Macaque/multiVeloMatrixFiles/Macaque_2PFC_multivelo_SpeciesComparison_dPSItestedGenes_Cell_States_mtx.csv'))
stateFile_ctrl <- fread(file.path(dataDir,'Human_vs_Macaque/multiVeloMatrixFiles/All_6Controls_multivelo_SpeciesComparison_dPSItestedGenes_Cell_States_mtx.csv'))
colnames(stateFile_case)[1] <- colnames(stateFile_ctrl)[1] <- "Sample_BC"

print("Reading in all info")
allInfo_case <- fread(file.path(dataDir,'Human_vs_Macaque/AllInfoFiles/Macaque_AllInfo_speciescomp_edit23_broad.gz'))
allInfo_case <- allInfo_case[,c(1:5,9)]

allInfo_ctrl <- fread(file.path(dataDir,'Human_vs_Macaque/AllInfoFiles/Human_AllInfo_speciescomp_edit23_broad.gz'))
allInfo_ctrl <- allInfo_ctrl[,c(1:5,9)]
colnames(allInfo_case) <- colnames(allInfo_ctrl) <- c("Read","Gene","CT","BC","UMI","ExonChain")

print("Reading in extra info")
g2c <- read.table(file.path(dataDir,'bidirectional_humanMacaque'), header = F) %>% distinct()
g2c <- g2c %>% separate(V1,into=c("ENS","Suffix"),sep = "[.]") %>% select(-Suffix)
colnames(g2c) <- c("ENS_M","Clear_M","Exon_M","ENS_H","Clear_H","Exon_H")

eL <- read.table(file.path(dataDir,'Human_vs_Macaque/testedExons_humanVsMacaque_filtered'), header = F)
colnames(eL) <- c("Exon_H","ENS_H","Clear_H","Exon_M","ENS_M","Clear_M")
eL <- eL %>% separate(ENS_M,into=c("ENS_M","Suffix"),sep = "[.]") %>% select(-Suffix)


#!/bin/R

# By Anoushka Joglekar 
# 03.2024

## Setup ------


## Functions -----

## Subsets out the state vector and known alternative exons for a gene of interest.

getDataPerGene <- function(stateFile,exonList,goi,gene2clear){
    #print("getDataPerGene")
    geneID <- gene2clear %>% filter(Clear == goi) %>% .$ENS
    print(stateFile$Sample_BC[1])
    if(length(unlist(strsplit(stateFile$Sample_BC[1],"_"))) == 4){
        goiStates <- stateFile %>% select(Sample_BC,all_of(goi)) %>% 
            separate(Sample_BC, into = c("Species","Sample","Region","Barc"),sep = "_") %>% 
            separate(Barc,into=c("BC","Suffix"),sep = "-")
    } else {
        goiStates <- stateFile %>% select(Sample_BC,all_of(goi)) %>% 
            separate(Sample_BC, into = c("Sample","Barc"),sep = "_") %>% 
            separate(Barc,into=c("BC","Suffix"),sep = "-")
    }
    
    subsetExons <- exonList %>% filter(Clear == goi) %>% .$Exon
    
    return(list(subsetExons,goiStates))
}

# ## Check that the exon is alternative and the contingency table passes the chi-sq criterion

checkChiSqCrit <- function(mat){
    #print("checkChiSqCrit")
    if (sum(mat) > 0 & max(colSums(mat))/sum(mat) <= 0.9 & 
        min(rowSums(mat))*min(colSums(mat))/sum(mat) > 5){
        return("pass")
    } else {return("fail")}
}

## If it passes the chisq criterion, calculate psi and log odds ration

logOddsCalc <- function(mat){## Set 0's to 0.5 and explicitly calculate log odds ratio for 2x2 matrices
    #print("logOddsCalc")
    if(nrow(mat)==2){
        mat[mat==0]=0.5
        lor = log2((mat[1,1]*mat[2,2])/(mat[2,1]*mat[1,2]))
        return(lor)
    } else {
        return("NA")}
}

reportStats <- function(mat){
    #print("reportStats")
    p = chisq.test(mat)$p.value
    lor = logOddsCalc(mat)
    return(list(p,lor))
}

## Function to calculate psi per state in each matrix

calcPSI = function(mat){
    #print("calcPSI")
    psiVec = mat[,1]/rowSums(mat)
    return(psiVec)
}

## Given an exon, calculate inclusion and exclusion counts per read. Test for association with accessibility state

perExon <- function(allInfo_mod,subset_states,eoi,goi){
    #print("perExon")
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
    #print("performAnalysisForCT")
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

analysisPerGene <- function(goi, stateFile, allInfo, allCTs,exonList,gene2clear){
    #print("analysisPerGene")
    geneID <- unique(gene2clear %>% filter(Clear == goi) %>% .$ENS)
    if(length(geneID) > 0){
        ss = getDataPerGene(stateFile,exonList,goi,gene2clear)
        allRes = do.call('rbind',lapply(allCTs, function(ct) performAnalysisForCT(ss,allInfo,geneID,goi,ct)))
        allRes <- as.data.frame(allRes)
        if(nrow(allRes) > 1){
            colnames(allRes) <- c("CT","Exon","statePSI","stateID","totalCounts","pvalue","LOR")
            allRes$GeneName <- goi
            allRes$GeneID <- geneID
            return(allRes)
        } else{return()}
    } else{return()}
}


## Function to identify genes of interest

getGeneList <- function(stateFile,exonList){
    #print("getGeneList")
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

#     print(paste("Total genes from exonList:",length(unique(exonList$Clear))))
    return(multiStateGenes)
}
                                    
                                    
runForSample <- function(allInfo, stateFile, exonList, gene2clear){
    
    allCTs <- grep(pattern = "Unknown",unique(allInfo$CT),invert = T,value = T)
    multiStateGenes <- getGeneList(stateFile,exonList)
    
    multiStateRes <- as.data.frame(do.call('rbind',mclapply(multiStateGenes, function(gene) {
    print(gene)
    analysisPerGene(gene, stateFile,allInfo, allCTs,exonList, gene2clear)},mc.cores = 8)))
    
#     multiStateRes <- as.data.frame(do.call('rbind',lapply(multiStateGenes[1:5], function(gene) {
#     print(gene)
#     analysisPerGene(gene, stateFile,allInfo, allCTs, exonList, gene2clear)})))
    
    multiStateRes <- multiStateRes %>% distinct()
    return(multiStateRes)
}
                                        
eL_control <- eL %>% select(contains("_H"))
colnames(eL_control) <- gsub("_H","",colnames(eL_control))
g2c_ctrl <- g2c %>% select(contains("_H"))
colnames(g2c_ctrl) <- gsub("_H","",colnames(g2c_ctrl))

eL_case <- eL %>% select(contains("_M"))
colnames(eL_case) <- gsub("_M","",colnames(eL_case))
g2c_case <- g2c %>% select(contains("_M"))
colnames(g2c_case) <- gsub("_M","",colnames(g2c_case))
                                        
                                        
caseRes <- runForSample(allInfo_case, stateFile_case, eL_case, g2c_case)
controlRes <- runForSample(allInfo_ctrl, stateFile_ctrl, eL_control, g2c_ctrl)
                                        
caseRes$Sample <- "MacaquePFC"
controlRes$Sample <- "HumanCortex"

caseRes <- caseRes %>% separate_rows(statePSI, stateID, sep = ",\\s*") %>%
mutate_at("statePSI",as.double) %>% mutate_at("stateID", as.factor) %>% distinct()

controlRes <- controlRes %>% separate_rows(statePSI, stateID, sep = ",\\s*") %>%
mutate_at("statePSI",as.double) %>% mutate_at("stateID", as.factor) %>% distinct()

cr2 = left_join(caseRes %>% rename(Exon_M = Exon), g2c %>% select(Exon_H,Exon_M))
contR2 = left_join(controlRes %>% rename(Exon_H = Exon), g2c %>% select(Exon_H,Exon_M))


fullMat <- full_join(cr2,contR2, by = c("CT","stateID","GeneName","Exon_H","Exon_M")) %>% 
    mutate(statePSI.x = ifelse(is.na(statePSI.x), -1, statePSI.x),
           Sample.x = ifelse(is.na(Sample.x), "MacaquePFC", Sample.x),
           statePSI.y = ifelse(is.na(statePSI.y), -1, statePSI.y),
          Sample.y = ifelse(is.na(Sample.y), "HumanCortex", Sample.y),)
head(fullMat,5)                                        
                                        
dir.create('../2024_02_29_multiVelo/output/v3_human_vs_macaque/')
write.table(fullMat,'../2024_02_29_multiVelo/output/v3_human_vs_macaque/macaque_vs_human_fullMat_altWorkflow.tab',sep = "\t",
           quote = F, row.name = F, col.names = T)                                        