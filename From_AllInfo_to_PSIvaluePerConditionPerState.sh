
## Allinfo file is a UMI corrected scisorseqr output which include concatenated information for all poly-exonic, spliced, barcoded, mapped, and 5’ and 3’ complete reads
## the format of Allinfo file:
### col1: readID, ensembl_geneID, cell-type, barcode, UMI, intron-chain, TSS, polyA, exon-chain, status, number of introns

#### split Allinfo file of each condition into 4 states #####
##### For region comparison by State ####################

### Input files
# AllInfo_VIS/AllInfo_PFC: UMI corrected Allinfo file derived from scisorseqr pipeline
# cellToStateFile: the output from multivelo, which is a matrix file composed of state values of each gene per cell. The column-names are gene symbols, the row-names are cell barcodes.
# ensembl2GeneName: the annotation file which includes the ensemble-gene-IDs and corresponding gene-symbols.

### output files
## Allinfo files of all cell-states, with a suffix of "state0.0","state1.0","state2.0","state3.0"
###############################
 
###split by state for AD vs Control  For all cell types

cd /athena/tilgnerlab/scratch/weh4002/Multivelo/state_PSI/state_assign/human_AD
cat AllInfo_Cases_Corrected_Incomplete_Celltype| awk -v cellToStateFile=All_9Cases_multivelo_dPSItestedGenes_Cell_States_mtx.csv -v ensembl2GeneName=hs.symbol2ID 'BEGIN {comm2="cat "cellToStateFile; lineCounter=0; while(comm2|getline){lineCounter++; if(lineCounter==1){n=split($1,geneOrder,","); if(geneOrder[1]!=""){print "PARSING ERROR 1" > "/dev/stderr"; exit(0);} } if(lineCounter % 100 == 0){print "lineCounter="lineCounter > "/dev/stderr";}    if(lineCounter>1){ n=split($1,geneState,",");  split(geneState[1],a,"_"); cellID=substr(a[2],1,length(a[2])-2); print cellID > "test.a"; for(i=2;i<=n;i++){geneCellState[cellID"\t"geneOrder[i]]=geneState[i]; print geneState[i] >> "test.a2";}}} comm3="cat "ensembl2GeneName; while(comm3|getline){ensembl2Gname[$2]=$1;}          } {cellState=geneCellState[$4"\t"ensembl2Gname[$2]]; print $4"\t"ensembl2Gname[$2] >> "test.a3"; print "\""cellState"\"" > "test.b";outputFile="human.AD.allCT.allInfoFile.state"cellState;  print $0 >> outputFile;       }'

cd /athena/tilgnerlab/scratch/weh4002/Multivelo/state_PSI/state_assign/human_control
cat AllInfo_Incomplete_Corrected_Controls_Celltypes| awk -v cellToStateFile=All_10Controls_multivelo_dPSItestedGenes_Cell_States_mtx.csv -v ensembl2GeneName=hs.symbol2ID 'BEGIN {comm2="cat "cellToStateFile; lineCounter=0; while(comm2|getline){lineCounter++; if(lineCounter==1){n=split($1,geneOrder,","); if(geneOrder[1]!=""){print "PARSING ERROR 1" > "/dev/stderr"; exit(0);} } if(lineCounter % 100 == 0){print "lineCounter="lineCounter > "/dev/stderr";}    if(lineCounter>1){ n=split($1,geneState,",");  split(geneState[1],a,"_"); cellID=substr(a[2],1,length(a[2])-2); print cellID > "test.a"; for(i=2;i<=n;i++){geneCellState[cellID"\t"geneOrder[i]]=geneState[i]; print geneState[i] >> "test.a2";}}} comm3="cat "ensembl2GeneName; while(comm3|getline){ensembl2Gname[$2]=$1;}          } {cellState=geneCellState[$4"\t"ensembl2Gname[$2]]; print $4"\t"ensembl2Gname[$2] >> "test.a3"; print "\""cellState"\"" > "test.b";outputFile="human.control.allCT.allInfoFile.state"cellState;  print $0 >> outputFile;       }'

####### extract the reads assigned to a specific cell-type in query, we showed the excitaroty neuron (ExN)  

cat human.AD.allCT.allInfoFile.state0.0 | grep 'ExN' | gzip -c  > human.AD.ExN.allInfoFile.state0.0.gz
cat human.AD.allCT.allInfoFile.state1.0 | grep 'ExN' | gzip -c  > human.AD.ExN.allInfoFile.state1.0.gz
cat human.AD.allCT.allInfoFile.state2.0 | grep 'ExN' | gzip -c  > human.AD.ExN.allInfoFile.state2.0.gz
cat human.AD.allCT.allInfoFile.state3.0 | grep 'ExN' | gzip -c  > human.AD.ExN.allInfoFile.state3.0.gz

cat human.control.allCT.allInfoFile.state0.0 | grep 'ExN' | gzip -c  > human.control.ExN.allInfoFile.state0.0.gz
cat human.control.allCT.allInfoFile.state1.0 | grep 'ExN' | gzip -c  > human.control.ExN.allInfoFile.state1.0.gz
cat human.control.allCT.allInfoFile.state2.0 | grep 'ExN' | gzip -c  > human.control.ExN.allInfoFile.state2.0.gz
cat human.control.allCT.allInfoFile.state3.0 | grep 'ExN' | gzip -c  > human.control.ExN.allInfoFile.state3.0.gz


#### For each cell type in query of each cell state, test for exons if they show significant differential usage in two conditions for comparison
#### the parameters are set for a very loose cutoff, which aimed for keeping all the exons if possible

### v1.1a_exonInclusion_CTspecific_case_control.sh is the code for testing if the exons are differentially utilized between two conditions
# numThreads: number of threads to be used.
# min_reads: minimum number of reads for sum of 2 allInfos for a given exon. 
# OL_fraction: the fraction of the reads for a given position must be either inclusion or exclusion.
# PathToAnnotation: gencode annotation
# PathToChromosomeFile: a list of chromosome names in query, which overs chr1~chr22, chrX and chrY
# celltypeList: cell-type in query, should be consistent with the Allinfo file col3.

MinNumReads=1
PathToChromosomeFile="/athena/tilgnerlab/scratch/caf4010/12_3_22_SingleNucleiMacaque/LongRead/all.chroms.1-22.X"
NumThreads=10
PathToAnnotation="/athena/tilgnerlab/scratch/weh4002/Multivelo/state_PSI/state_assign/gencode.v34.annotation.gtf.gz"
PathToCellTypeFile="/athena/tilgnerlab/scratch/weh4002/Multivelo/state_PSI/state_assign/human.AD.celltypeList"
OL_fraction=0


################# state0
caseListPath="/athena/tilgnerlab/scratch/weh4002/Multivelo/state_PSI/state_assign/ADvsControl_ExN/State0_AD/caseList"
controlListPath="/athena/tilgnerlab/scratch/weh4002/Multivelo/state_PSI/state_assign/ADvsControl_ExN/State0_AD/controlList"
echo "/scratchLocal/" > tmpdirs; time bash /athena/tilgnerlab/scratch/nab4004/snisor-new-probes/cassette-exon-scripts/v1.1a_exonInclusion_CTspecific_case_control.sh $caseListPath $controlListPath tmpdirs $PathToChromosomeFile $NumThreads $PathToAnnotation /athena/tilgnerlab/scratch/nab4004/snisor-new-probes/cassette-exon-scripts/other-scripts/ 0.00 1.00 $MinNumReads zcat $OL_fraction $PathToCellTypeFile &> human_ADvsControl_state1_test


################# state1
caseListPath="/athena/tilgnerlab/scratch/weh4002/Multivelo/state_PSI/state_assign/ADvsControl_ExN/State1_AD/caseList"
controlListPath="/athena/tilgnerlab/scratch/weh4002/Multivelo/state_PSI/state_assign/ADvsControl_ExN/State1_AD/controlList"
echo "/scratchLocal/" > tmpdirs; time bash /athena/tilgnerlab/scratch/nab4004/snisor-new-probes/cassette-exon-scripts/v1.1a_exonInclusion_CTspecific_case_control.sh $caseListPath $controlListPath tmpdirs $PathToChromosomeFile $NumThreads $PathToAnnotation /athena/tilgnerlab/scratch/nab4004/snisor-new-probes/cassette-exon-scripts/other-scripts/ 0.00 1.00 $MinNumReads zcat $OL_fraction $PathToCellTypeFile &> human_ADvsControl_state1_test


################# state2
caseListPath="/athena/tilgnerlab/scratch/weh4002/Multivelo/state_PSI/state_assign/ADvsControl_ExN/State2_AD/caseList"
controlListPath="/athena/tilgnerlab/scratch/weh4002/Multivelo/state_PSI/state_assign/ADvsControl_ExN/State2_AD/caseList/controlList"

echo "/scratchLocal/" > tmpdirs; time bash /athena/tilgnerlab/scratch/nab4004/snisor-new-probes/cassette-exon-scripts/v1.1a_exonInclusion_CTspecific_case_control.sh $caseListPath $controlListPath tmpdirs $PathToChromosomeFile $NumThreads $PathToAnnotation /athena/tilgnerlab/scratch/nab4004/snisor-new-probes/cassette-exon-scripts/other-scripts/ 0.00 1.00 $MinNumReads zcat $OL_fraction $PathToCellTypeFile &> human_ADvsControl_state2_test


################# state3
caseListPath="/athena/tilgnerlab/scratch/weh4002/Multivelo/state_PSI/state_assign/ADvsControl_ExN/State3_AD/caseList"
controlListPath="/athena/tilgnerlab/scratch/weh4002/Multivelo/state_PSI/state_assign/ADvsControl_ExN/State3_AD/caseList/controlList"

echo "/scratchLocal/" > tmpdirs; time bash /athena/tilgnerlab/scratch/nab4004/snisor-new-probes/cassette-exon-scripts/v1.1a_exonInclusion_CTspecific_case_control.sh $caseListPath $controlListPath tmpdirs $PathToChromosomeFile $NumThreads $PathToAnnotation /athena/tilgnerlab/scratch/nab4004/snisor-new-probes/cassette-exon-scripts/other-scripts/ 0.00 1.00 $MinNumReads zcat $OL_fraction $PathToCellTypeFile &> human_ADvsControl_state3_test



#####################
# In the output from the last step, we can get  CONTROLS.sampleBoth.inc.exc.tab.gz and CASES.sampleBoth.inc.exc.tab.gz
# col1: exonID, col2-4: perfect included read #, upstream partial included read#, dnstream partial included read#, col5: excluded read#, col6: total raw read#


##### process the output to get #col1: exonID, col2: total reads, col3: PSI (percentage spliced in) 
# state 0  output folder
zcat CONTROLS.sampleBoth.inc.exc.tab.gz | awk -F"\t" '($2+$3+$4+$5)>0 {print $1"\t"($2+$3+$4+$5)"\t"($2+$3+$4)/($2+$3+$4+$5)}' | sort -u  > state0_CONTROLS.sampleBoth_totalNum_Inc.tab
zcat CASES.sampleBoth.inc.exc.tab.gz | awk -F"\t" '($2+$3+$4+$5)>0 {print $1"\t"($2+$3+$4+$5)"\t"($2+$3+$4)/($2+$3+$4+$5)}' | sort -u > state0_CASES.sampleBoth_totalNum_Inc.tab
# state 1  output folder
zcat CONTROLS.sampleBoth.inc.exc.tab.gz | awk -F"\t" '($2+$3+$4+$5)>0 {print $1"\t"($2+$3+$4+$5)"\t"($2+$3+$4)/($2+$3+$4+$5)}' | sort -u  > state1_CONTROLS.sampleBoth_totalNum_Inc.tab
zcat CASES.sampleBoth.inc.exc.tab.gz | awk -F"\t" '($2+$3+$4+$5)>0 {print $1"\t"($2+$3+$4+$5)"\t"($2+$3+$4)/($2+$3+$4+$5)}' | sort -u > state1_CASES.sampleBoth_totalNum_Inc.tab
# state 2  output folder
zcat CONTROLS.sampleBoth.inc.exc.tab.gz | awk -F"\t" '($2+$3+$4+$5)>0 {print $1"\t"($2+$3+$4+$5)"\t"($2+$3+$4)/($2+$3+$4+$5)}' | sort -u  > state2_CONTROLS.sampleBoth_totalNum_Inc.tab
zcat CASES.sampleBoth.inc.exc.tab.gz | awk -F"\t" '($2+$3+$4+$5)>0 {print $1"\t"($2+$3+$4+$5)"\t"($2+$3+$4)/($2+$3+$4+$5)}' | sort -u > state2_CASES.sampleBoth_totalNum_Inc.tab
# state 3  output folder
zcat CONTROLS.sampleBoth.inc.exc.tab.gz | awk -F"\t" '($2+$3+$4+$5)>0 {print $1"\t"($2+$3+$4+$5)"\t"($2+$3+$4)/($2+$3+$4+$5)}' | sort -u  > state3_CONTROLS.sampleBoth_totalNum_Inc.tab
zcat CASES.sampleBoth.inc.exc.tab.gz | awk -F"\t" '($2+$3+$4+$5)>0 {print $1"\t"($2+$3+$4+$5)"\t"($2+$3+$4)/($2+$3+$4+$5)}' | sort -u > state3_CASES.sampleBoth_totalNum_Inc.tab






