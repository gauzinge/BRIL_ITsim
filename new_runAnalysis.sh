#!/bin/bash

eval $(scramv1 runtime -sh) || echo "The command 'cmsenv' failed!"

input="full_RelValNuGun_Files_PU10.txt"
PUIN=$(echo ${input} | cut -d'.' -f 1 | cut -d'_' -f 3 | cut --complement -c 1,2)

function runAnalysis() {
  echo '############################################'
  echo "Running BRIL IT simulation analysis code"
  echo '=> Picking up data files with PU' ${PUIN}
  echo '############################################'
  echo 'Processing files...'
  TAG=1
  while IFS= read -r filename
  do
    cmsRun ITclusterAnalyzer/python/ITclusterAnalyzer_cfg.py print \
                inputFiles=file:${filename} \
                outputFile=temp.root \
                tag=${TAG}
    TAG=`expr $TAG + 1`
  done < "${input}"
}

runAnalysis

echo "Now merging output histograms"
command="hadd -f summary_PU_${PUIN}.root"
for rootfile in ${PWD}/temp_?.root; do
  command+=" ${rootfile}"
done
echo $command
${command}

