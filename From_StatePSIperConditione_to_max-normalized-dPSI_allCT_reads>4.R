library(RColorBrewer)
library(ggplot2)


celltypes = c("ExN","InN","OLIG","MG","ASC")
dpsi2fdr2state.sig.all = NULL
largedpsi.sig.noConfirm.all = NULL

for (i in 1:length(celltypes))
  
{
ct.path = paste0("/Users/apple/Wen/Neuroscience_papers/Human_AD_Multiome/Multivelo/state_PSI/state_reads_perExon/AD_allCT/ADvsControl_", celltypes[i])
ct.state0.path = paste0(ct.path, "/State0")
ct.state1.path = paste0(ct.path, "/State1")
ct.state2.path = paste0(ct.path, "/State2")
ct.state3.path = paste0(ct.path, "/State3")

############ filter for exons with at least 5 reads for each state ###########

### state 0 
setwd(ct.state0.path)
state0.case  = read.table("state0_CASES.sampleBoth_totalNum_Inc.tab", sep="\t", header = F, col.names= c("exon","case.state0.reads","case.state0.PSI"))
state0.control  = read.table("state0_CONTROLS.sampleBoth_totalNum_Inc.tab", sep="\t", header = F, col.names= c("exon","ctl.state0.reads","ctl.state0.PSI"))

state0 = merge(state0.case, state0.control, by.x = "exon", by.y = "exon", all.x = T, all.y = T)
state0 = subset(state0, subset = case.state0.reads>=5 | ctl.state0.reads>=5)
state0.both = merge(state0.case, state0.control, by.x = "exon", by.y = "exon", all.x = F, all.y = F)
state0.both = subset(state0.both, subset = case.state0.reads>=5 | ctl.state0.reads>=5)


#### state 1
setwd(ct.state1.path)
state1.case  = read.table("state1_CASES.sampleBoth_totalNum_Inc.tab", sep="\t", header = F, col.names= c("exon","case.state1.reads","case.state1.PSI"))
state1.control  = read.table("state1_CONTROLS.sampleBoth_totalNum_Inc.tab", sep="\t", header = F, col.names= c("exon","ctl.state1.reads","ctl.state1.PSI"))

state1 = merge(state1.case, state1.control, by.x = "exon", by.y = "exon", all.x = T, all.y = T)
state1 = subset(state1, subset = case.state1.reads>=5 | ctl.state1.reads>=5)
state1.both = merge(state1.case, state1.control, by.x = "exon", by.y = "exon", all.x = F, all.y = F)
state1.both = subset(state1.both, subset = case.state1.reads>=5 | ctl.state1.reads>=5)


#### state 2
setwd(ct.state2.path)
state2.case  = read.table("state2_CASES.sampleBoth_totalNum_Inc.tab", sep="\t", header = F, col.names= c("exon","case.state2.reads","case.state2.PSI"))
state2.control  = read.table("state2_CONTROLS.sampleBoth_totalNum_Inc.tab", sep="\t", header = F, col.names= c("exon","ctl.state2.reads","ctl.state2.PSI"))

state2 = merge(state2.case, state2.control, by.x = "exon", by.y = "exon", all.x = T, all.y = T)
state2 = subset(state2, subset = case.state2.reads>=5 | ctl.state2.reads>=5)
state2.both = merge(state2.case, state2.control, by.x = "exon", by.y = "exon", all.x = F, all.y = F)
state2.both = subset(state2.both, subset = case.state2.reads>=5 | ctl.state2.reads>=5)


#### state 3
setwd(ct.state3.path)
state3.case  = read.table("state3_CASES.sampleBoth_totalNum_Inc.tab", sep="\t", header = F, col.names= c("exon","case.state3.reads","case.state3.PSI"))
state3.control  = read.table("state3_CONTROLS.sampleBoth_totalNum_Inc.tab", sep="\t", header = F, col.names= c("exon","ctl.state3.reads","ctl.state3.PSI"))

state3 = merge(state3.case, state3.control, by.x = "exon", by.y = "exon", all.x = T, all.y = T)
state3 = subset(state3, subset = case.state3.reads>=5 | ctl.state3.reads>=5)
state3.both = merge(state3.case, state3.control, by.x = "exon", by.y = "exon", all.x = F, all.y = F)
state3.both = subset(state3.both, subset = case.state3.reads>=5 | ctl.state3.reads>=5)

####
state01 = merge(state0, state1, by.x = "exon", by.y = "exon", all.x = T, all.y = T)
state01.paired = merge(state0.both, state1.both, by.x = "exon", by.y = "exon", all.x = T, all.y = T)

state23 = merge(state2, state3, by.x = "exon", by.y = "exon", all.x = T, all.y = T)
state23.paired = merge(state2.both, state3.both, by.x = "exon", by.y = "exon", all.x = T, all.y = T)


state.all = merge(state01, state23, by.x = "exon", by.y = "exon", all.x = T, all.y = T)
state.all.paired = merge(state01.paired, state23.paired, by.x = "exon", by.y = "exon", all.x = T, all.y = T)

state.all.paired$state0.dpsi = state.all.paired$case.state0.PSI- state.all.paired$ctl.state0.PSI
state.all.paired$state1.dpsi = state.all.paired$case.state1.PSI- state.all.paired$ctl.state1.PSI
state.all.paired$state2.dpsi = state.all.paired$case.state2.PSI- state.all.paired$ctl.state2.PSI
state.all.paired$state3.dpsi = state.all.paired$case.state3.PSI- state.all.paired$ctl.state3.PSI

state.all$state0.dpsi = state.all$case.state0.PSI- state.all$ctl.state0.PSI
state.all$state1.dpsi = state.all$case.state1.PSI- state.all$ctl.state1.PSI
state.all$state2.dpsi = state.all$case.state2.PSI- state.all$ctl.state2.PSI
state.all$state3.dpsi = state.all$case.state3.PSI- state.all$ctl.state3.PSI


##### read in overall-dPSI with FDR
dpsi.path = paste0("/Users/apple/Wen/Neuroscience_papers/Human_AD_Multiome/Multivelo/overall_dPSI_FDR/cases_vs_controls.dPSI.FDR_AD.", celltypes[i])
dpsi2fdr = read.table(dpsi.path, header=F, col.names = c("exon","geneID","dPSI","FDR","celltype"))


###### only consider the exons which are have been tested with overall-dPSI (dPSI across all cell states)
dpsi2fdr2state = merge(dpsi2fdr, state.all, by.x = "exon", by.y = "exon", all.x = F, all.y = F)

#### calculate the normalized state-dPSI for each exon per cell type
dpsi2fdr2state$norm.state0.dpsi = dpsi2fdr2state$state0.dpsi/dpsi2fdr2state$dPSI
dpsi2fdr2state$norm.state1.dpsi = dpsi2fdr2state$state1.dpsi/dpsi2fdr2state$dPSI
dpsi2fdr2state$norm.state2.dpsi = dpsi2fdr2state$state2.dpsi/dpsi2fdr2state$dPSI
dpsi2fdr2state$norm.state3.dpsi = dpsi2fdr2state$state3.dpsi/dpsi2fdr2state$dPSI

dpsi2fdr2state$norm.state0.dpsi = as.numeric(dpsi2fdr2state$norm.state0.dpsi)
dpsi2fdr2state$norm.state1.dpsi = as.numeric(dpsi2fdr2state$norm.state1.dpsi)
dpsi2fdr2state$norm.state2.dpsi = as.numeric(dpsi2fdr2state$norm.state2.dpsi)
dpsi2fdr2state$norm.state3.dpsi = as.numeric(dpsi2fdr2state$norm.state3.dpsi)

###### only keep the exons which have  significant 
dpsi2fdr2state.sig = subset(dpsi2fdr2state, subset= FDR<0.05)


###### find the maximum normalized dPSI among all states
for (k in 1:nrow(dpsi2fdr2state.sig))
{
  dpsi2fdr2state.sig$max.norm.state.psi[k] = max(dpsi2fdr2state.sig[k,26:29],na.rm=TRUE)
  
}



###### 
dpsi2fdr2state.sig = subset(dpsi2fdr2state.sig, subset = max.norm.state.psi!="-Inf")
### confirmed normalized state dPSI
dpsi2fdr2state.sig.confirm = subset(dpsi2fdr2state.sig, subset = max.norm.state.psi>=1)
### non-confirmed normalized state dPSI
dpsi2fdr2state.sig.noconfirm = subset(dpsi2fdr2state.sig, subset = max.norm.state.psi<1)

setwd(ct.path)
file.name = paste0("Human.ADvsControl_Normalized.state.PSI_for_significantExons_",celltypes[i],"_reads>4.tsv")
write.table(dpsi2fdr2state.sig, file=file.name, sep = "\t", quote = F, col.names = T, row.names = F)

file.name = paste0("Human.ADvsControl_Normalized.state.PSI_for_significantExons_Confirmeed_",celltypes[i],"_reads>4.tsv")
write.table(dpsi2fdr2state.sig.confirm, file=file.name, sep = "\t", quote = F, col.names = T, row.names = F)

file.name = paste0("Human.ADvsControl_Normalized.state.PSI_for_significantExons_No_Confirmeed_",celltypes[i],"_reads>4.tsv")
write.table(dpsi2fdr2state.sig.noconfirm, file=file.name, sep = "\t", quote = F, col.names = T, row.names = F)


dpsi2fdr2state.sig.all = rbind(dpsi2fdr2state.sig.all,dpsi2fdr2state.sig)
dpsi2fdr2state.sig.confirm.all = rbind(dpsi2fdr2state.sig.confirm.all,dpsi2fdr2state.sig.confirm)
dpsi2fdr2state.sig.noconfirm.all = rbind(dpsi2fdr2state.sig.noconfirm.all,dpsi2fdr2state.sig.noconfirm)


rm(list=setdiff(ls(), c("dpsi2fdr2state.sig.all","celltypes","largedpsi.sig.noConfirm.all",)))

}


setwd("/Users/apple/Wen/Neuroscience_papers/Human_AD_Multiome/Multivelo/state_PSI/state_reads_perExon/AD_allCT/")
file.name1 = paste0("Human.ADvsControl_Normalized.state.PSI_for_significantExons_all_atleast1condition.reads>4_atleast1.tsv")
write.table(dpsi2fdr2state.sig.all, file=file.name1, sep = "\t", quote = F, col.names = T, row.names = F)

setwd("/Users/apple/Wen/Neuroscience_papers/Human_AD_Multiome/Multivelo/state_PSI/state_reads_perExon/AD_allCT/")
file.name2 = paste0("Human.ADvsControl_Normalized.state.PSI_for_significantExons_all_confirmed_atleast1condition.reads>4_atleast1.tsv")
write.table(dpsi2fdr2state.sig.confirm.all, file=file.name2, sep = "\t", quote = F, col.names = T, row.names = F)

setwd("/Users/apple/Wen/Neuroscience_papers/Human_AD_Multiome/Multivelo/state_PSI/state_reads_perExon/AD_allCT/")
file.name3 = paste0("Human.ADvsControl_Normalized.state.PSI_for_significantExons_all_No_confirmed_atleast1condition.reads>4_atleast1.tsv")
write.table(dpsi2fdr2state.sig.noconfirm.all, file=file.name3, sep = "\t", quote = F, col.names = T, row.names = F)



stats = data.frame(matrix(NA,ncol = 7,nrow = length(unique(dpsi2fdr2state.sig.all$celltype))))
colnames(stats) = c("CT","[-Inf,0.9)","[0.9,1)","[1,Inf)","[-Inf,0.9).ratio","[0.9,1).ratio","[1,Inf).ratio")
celltypes = c("ExN","InN","OLIG","MG","ASC")

for (i in 1:nrow(stats))
{
  ct.state = subset(dpsi2fdr2state.sig.all, subset = celltype == celltypes[i])
  stats$CT[i] = celltypes[i]
  ct.num = nrow(ct.state)
  stats$`[-Inf,0.9)`[i] = nrow(subset(ct.state, subset = max.norm.state.psi<0.9))
  stats$`[0.9,1)`[i] = nrow(subset(ct.state, subset = max.norm.state.psi>=0.9 & max.norm.state.psi<1))
  stats$`[1,Inf)`[i] = nrow(subset(ct.state, subset = max.norm.state.psi>=1))
  stats$`[-Inf,0.9).ratio`[i] = nrow(subset(ct.state, subset = max.norm.state.psi<0.9))/ct.num
  stats$`[0.9,1).ratio`[i] = nrow(subset(ct.state, subset = max.norm.state.psi>=0.9 & max.norm.state.psi<1))/ct.num
  stats$`[1,Inf).ratio`[i] = nrow(subset(ct.state, subset = max.norm.state.psi>=1))/ct.num
}



state.plot = data.frame(matrix(NA, nrow = 3*length(celltypes),ncol = 4))
colnames(state.plot) = c("celltype","range.ratio","group","range.num")

for (i in 1:nrow(stats))
{
  start = (i-1)*3+1
  end = (i-1)*3+3
  state.plot$celltype[start:end] = stats[i,1]
  state.plot$range.num[start:end] = c(stats[i,2:4])
  state.plot$range.ratio[start:end] = c(stats[i,5:7])
  state.plot$group[start:end] = colnames(stats)[2:4]
  
}

state.plot$range.ratio =as.numeric(state.plot$range.ratio)
state.plot$range.num = as.numeric(state.plot$range.num)
state.plot$range.ratio = 100*state.plot$range.ratio


setwd("/Users/apple/Wen/Neuroscience_papers/Human_AD_Multiome/Multivelo/state_PSI/state_reads_perExon/AD_allCT/")

write.table(state.plot, file = "For_plot_ADvsControl_sig.exons_max.norm.state_3ranges_Stats_withExonNumber_allCT_atleast1condition.reads>4.tsv",row.names = F, col.names = T, sep = "\t", quote = F)

######### plot

setwd("/Users/apple/Wen/Neuroscience_papers/Human_AD_Multiome/Multivelo/state_PSI/state_reads_perExon/AD_allCT/")

pdf.file2= "Barplot_ADvsControl_sig.exons_max.norm.state_3ranges_withExonNumber_atleast1condition.reads>4_allCT.pdf"
p2.1 <- ggplot(state.plot, aes(x = celltype, y = range.ratio, fill = group ))+
  geom_bar(position="stack", stat="identity") + theme_bw() + geom_text(aes(label = range.num), size = 8, hjust = 0.5, vjust = 1, position="stack") + theme(text=element_text(size=25)) + xlab("") + ylab("% Exon") + scale_fill_brewer(palette="Pastel2")
print(p2.1)
ggsave(filename=pdf.file2,width = 8, height = 8, dpi = 300)

pdf.file3= "Barplot_ADvsControl_sig.exons_max.norm.state_3ranges_atleast1condition.reads>4_allCT.pdf"
p2.2 <- ggplot(state.plot, aes(x = celltype, y = range.ratio, fill = group ))+
  geom_bar(position="stack", stat="identity") + theme_bw() + theme(text=element_text(size=25)) + xlab("") + ylab("% Exon") + scale_fill_brewer(palette="Pastel2")
print(p2.2)
ggsave(filename=pdf.file3,width = 8, height = 8, dpi = 300)



####### density plot of ExN
dpsi2fdr2state.sig.exn = subset(dpsi2fdr2state.sig.all, subset = celltype == "ExN")

pdf("DensityPlot_ADvsControl_significant.exons_max.norm.state.dPSI_ExN_new_atleast1condition.reads>4_v1.pdf",7,3)
p1 <- ggplot(dpsi2fdr2state.sig.exn, aes(x=max.norm.state.psi)) + geom_density(lwd = 1.5) + theme_classic()+ geom_density(fill="#a8c2fb", color="#e9ecef", alpha=0.8)+ xlab("max.norm.state.dPSI") + ylab("Density")  + theme( axis.text=element_text(size=18), axis.title=element_text(size=20)) 
p1
dev.off()


pdf("DensityPlot_ADvsControl_significant.exons_max.norm.state.dPSI_ExN_new_v2.pdf",5,3)
p1 <- ggplot(dpsi2fdr2state.sig.exn, aes(x=max.norm.state.psi)) + geom_density(lwd = 1.5) + theme_classic()+ geom_density(fill="#a8c2fb", color="#e9ecef", alpha=0.8)+ xlab("max.norm.state.dPSI") + ylab("Density")  + theme( axis.text=element_text(size=18), axis.title=element_text(size=20)) 
p1
dev.off()

pdf("DensityPlot_ADvsControl_significant.exons_max.norm.state.dPSI_ExN_new_v3.pdf",8,2)
p1 <- ggplot(dpsi2fdr2state.sig.exn, aes(x=max.norm.state.psi)) + geom_density(lwd = 1.5) + theme_classic()+ geom_density(fill="#a8c2fb", color="#e9ecef", alpha=0.8)+ xlab("max.norm.state.dPSI") + ylab("Density")  + theme( axis.text=element_text(size=18), axis.title=element_text(size=20)) 
p1
dev.off()

