# Multivelo notebooks for 3 comparisons
# region comparison
Macaque_2PFC_dPSI_RegionComparison_genes_ForGithub.ipynb
Macaque_2VIS_dPSI__RegionComparison_genes-ForGithub.ipynb

# species comparison
Macaque_2PFC_dPSI_SpeciesComparison-ForGithub.ipynb
Multivelo_6Controls_SpeciesComparison_dPSI_tested_genes_broad_celltypes-ForGithub.ipynb

# AD vs Control
Multivelo_10Controls_dPSI_tested_genes_broad_celltypes-forGithub.ipynb
Multivelo_9Cases_dPSI_tested_genes_broad_celltypes-ForGithub.ipynb

#############################################################################
# we use the log-odds-ratio (abbreviated as LOR) to quantify the difference in inclusion between both states.
LOR_ADvsControl_withinHuman_altWorkflow.R
LOR_regionComparison_withinMacaque_altWorkflow
LOR_speciesComp_altWorkflow.R

############################################################################
# Alternative splicing analysis per cell state defined by multivelo output

#v1.1a_exonInclusion_CTspecific_case_control.sh
#####################
Compares exon inclusion and exclusion between 2 groups

#From_AllInfo_to_PSIvaluePerConditionPerState.sh
####################
a) Split the allinfo file of each condition by cell state information provided by Multivelo output
b) Perform differential splicing analysis between two conditions per cell state by calling the script 'v1.1a_exonInclusion_CTspecific_case_control.sh'
c) calculate the PSI (percentage spliced in) value per condition per cell state for all the exons

#From_StatePSIperConditione_to_max-normalized-dPSI_allCT_reads>4.R
#######################
a) calculate the state dPSI between conditions  if applicable
b) subset the exons which have been tested significant with an overall dPSI (across all cell states)



c) normalize the state dPSI  by overall dPSI (across all cell states)
d) find the maximum normalized state dPSI 
e) make a density plot of the maximum normalized state dPSI of each tested exon


