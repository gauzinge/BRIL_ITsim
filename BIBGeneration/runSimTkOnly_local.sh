#!/bin/bash

# a helper function to parse PU values
function list_include_item {
    local list="$1"
    local item="$2"
    if [[ $list =~ (^|[[:space:]])"$item"($|[[:space:]])  ]] ; then
        # yes, list include item
        result=0
    else
        result=1
    fi
    return $result
}

#assign the command line arguments
PU=$1
NEVENTS=$2
JOBID=$3
#STAGED=$4

################################################################################
##CHANGE ME ACCORDING TO YOUR NEEDS
################################################################################
NTHREADS=100

#FILE=file:
INFILE=/afs/cern.ch/work/p/pkicsiny/private/cmssw/CMSSW_11_2_0_pre6/src/BRIL_ITsim/BIBGeneration/test_generator/BeamHalo_fullcms.root
OUTDIR=/afs/cern.ch/work/p/pkicsiny/private/cmssw/CMSSW_11_2_0_pre6/src/BRIL_ITsim/BIBGeneration/test_simulator

################################################################################
#PRINT THE ARGUMENTS SUMMARY
################################################################################

echo '###########################################################################'
echo 'Configuration: '
echo 'Pileup Average: '${PU}
echo 'Number of Events: '${NEVENTS}
echo 'JobId: '${JOBID}
echo 'Input File (HEPMC generated): ' ${INFILE}
echo 'OutputDirectory: ' ${OUTDIR}
echo 'Number of Threads: ' ${NTHREADS}
echo '###########################################################################'

################################################################################
#SETUP CMSSW FRAMEWORK
################################################################################

#Extract sandbox
tar -xf sandbox.tar.bz2
#Keep track of release sandbox version
basedir=$PWD
rel=$(echo CMSSW_*)
arch=$(ls $rel/.SCRAM/|grep slc) || echo "Failed to determine SL release!"
old_release_top=$(awk -F= '/RELEASETOP/ {print $2}' $rel/.SCRAM/slc*/Environment) || echo "Failed to determine old releasetop!"
 
# Creating new release
# This isdone so e.g CMSSW_BASE and other variables are not hardcoded to the sandbox setting paths
# which will not exist here
  
echo ">>> creating new release $rel"
mkdir tmp
cd tmp
scramv1 project -f CMSSW $rel
new_release_top=$(awk -F= '/RELEASETOP/ {print $2}' $rel/.SCRAM/slc*/Environment)
cd $rel
echo ">>> preparing sandbox release $rel"
  
#for i in bin cfipython config lib module python src; do
for i in bin cfipython config lib python src; do
    rm -rf "$i"
    mv "$basedir/$rel/$i" .
done
     
echo ">>> fixing python paths"
for f in $(find -iname __init__.py); do
   sed -i -e "s@$old_release_top@$new_release_top@" "$f"
done
         
eval $(scramv1 runtime -sh) || echo "The command 'cmsenv' failed!"
cd "$basedir"
echo "[$(date '+%F %T')] wrapper ready"


################################################################################
##RUN THE ACTUAL SIMULATION
################################################################################

echo "Running the full simulation in one step from directory ${PWD}!"
#local
#command="cmsRun python/BH_SimTrigRec.py print \
#cluster
command="cmsRun python/BH_SimTrigRec.py print \
        nEvents=${NEVENTS} \
	inputFile=file:${INFILE} \
        nThreads=${NTHREADS} \
        jobId=${JOBID} \
        outputDirectory=file:${OUTDIR}"

echo 'Command: ' ${command}
${command}

################################################################################
##CLEANING UP BEHIND MYSELF
################################################################################

echo "Done running the generation"
echo "Cleaning up behing me"
rm -rf tmp/
rm -rf CMSSW_11_2_0_pre6/

