#!/bin/bash

PUIN=$1
DATAPATH=/eos/user/g/gauzinge/PUdata/
declare -a list=("0.0" "0.5" "1.0" "1.5" "2.0")

eval $(scramv1 runtime -sh) || echo "The command 'cmsenv' failed!"

function runAnalysis() {
    local DATAPATH=$1
    local PUIN=$2
    echo '####################################################################'
    echo "running BRIL IT Cluster analysis"
    echo 'Picking up data files with PU ' ${PUIN} ' from directory: ' ${DATAPATH}
    echo '####################################################################'
    echo 'Processing files:'
    for filename in ${DATAPATH}*.root; do
        #get the Proc ID
        TAG=$(echo ${filename} | cut -d'.' -f 3)
        tmp=$(echo ${filename} | cut -d'_' -f 4)
        MAIN=$(echo ${tmp} | cut -d'.' -f 1)
        DECIMAL=$(echo ${filename} | cut -d'.' -f 2)
        PU=${MAIN}.${DECIMAL}
        #if [[ "$PUIN" -eq "$PU" ]]; then
        if [ "$PUIN" = "$PU" ]; then
        echo ${filename} 'with PU ' ${PU}
            cmsRun ITclusterAnalyzer/python/ITclusterAnalyzer_cfg.py print \
                inputFiles=file:${filename} \
                outputFile=temp.root \
                tag=${TAG}
        fi
    done
    
    echo "Now merging output histograms"
    command="hadd -f summary_PU_${PUIN}.root"
    for rootfile in ${PWD}/temp_?.root; do
        command+=" ${rootfile}"
    done
    echo $command
    ${command}
}

if test -z "$PUIN" 
then
   echo "No Pileup number specified - please run as source runAnalysis.sh PU"
   for i in "${list[@]}"
   do
       runAnalysis ${DATAPATH} ${i}
   done
else
    runAnalysis ${DATAPATH} ${PUIN}
fi

echo "Cleaning up..."
rm temp*.root
