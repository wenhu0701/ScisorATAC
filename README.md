# ScisorATAC


#v1.1a_exonInclusion_CTspecific_case_control.sh
Compares exon inclusion and exclusion between 2 groups

#parameter settings:
# numThreads: number of threads to be used.
# min_reads: minimum number of reads for sum of 2 allInfos for a given exon. 
# OL_fraction: the fraction of the reads for a given position must be either inclusion or exclusion.
# PathToAnnotation: gencode annotation
# PathToChromosomeFile: a list of chromosome names in query, which overs chr1~chr22, chrX and chrY
# celltypeList: cell-type in query, should be consistent with the Allinfo file col3.

#From_AllInfo_to_PSIvaluePerConditionPerState.sh
a) Split the allinfo file of each condition by cell state information provided by Multivelo output
b) Perform differential splicing analysis between two conditions per cell state by calling the script 'v1.1a_exonInclusion_CTspecific_case_control.sh'
c) calculate the PSI (percentage spliced in) value per condition per cell state for all the exons

#From_StatePSIperConditione_to_max-normalized-dPSI_allCT_reads>4.R
a) calculate the state dPSI between conditions  if applicable
b) subset the exons which have been tested significant with an overall dPSI (across all cell states)
c) normalize the state dPSI  by overall dPSI (across all cell states)
d) find the maximum normalized state dPSI 
e) make a density plot of the maximum normalized state dPSI of each tested exon


