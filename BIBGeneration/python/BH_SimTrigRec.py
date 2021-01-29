# Auto generated configuration file
# using: 
# Revision: 1.163 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: Configuration/Generator/python/TTbar_cfi.py -s GEN,SIM,DIGI,L1,DIGI2RAW,HLT --eventcontent RAWSIM --datatier GEN-SIM-RAW --conditions FrontierConditions_GlobalTag,IDEAL_V11::All -n 10 --no_exec
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

from Configuration.StandardSequences.Eras import eras

# In the line below 'analysis' is an instance of VarParsing object 
options = VarParsing ('analysis')

# Here we have defined our own two VarParsing options 
# add a list of strings for events to process
options.register ('nEvents',
                                 10,
                                 #-1,
                                 VarParsing.multiplicity.singleton,
                                 VarParsing.varType.int,
                  "The number of events to generate: 10")
options.register ('inputFile',
		  'file:/afs/cern.ch/work/p/pkicsiny/private/cmssw/CMSSW_11_2_0_pre6/src/BRIL_ITsim/BIBGeneration/BeamHalo.0.root',
                  #'file:/afs/cern.ch/user/g//BIBSim/CMSSW_11_2_0_pre6/src/BRIL_ITsim/BIBGeneration/BeamHalo.0.root',
                                 VarParsing.multiplicity.singleton,
                                 VarParsing.varType.string,
                  "The input file")
options.register ('nThreads',
                                 1,
                                 VarParsing.multiplicity.singleton,
                                 VarParsing.varType.int,
                  "The number of threads to use: 1")
options.register ('jobId',
                                 0,
                                 VarParsing.multiplicity.singleton,
                                 VarParsing.varType.int,
                  "The job Id: 0")
options.register ('outputDirectory',
                  'file:/afs/cern.ch/work/p/pkicsiny/private/cmssw/CMSSW_11_2_0_pre6/src/BRIL_ITsim/BIBGeneration/BeamHalo.0.root',
		  #'file:/afs/cern.ch/user/g/gauzinge/BIBSim/CMSSW_11_2_0_pre6/src/BRIL_ITsim/BIBGeneration/',
                  # 'file:/afs/cern.ch/work/g/gauzinge/public/',
                                 VarParsing.multiplicity.singleton,
                                 VarParsing.varType.string,
                  "The output directory")

options.parseArguments()
options.outputFile=options.outputDirectory+'/BeamHaloReco.'+str(options.jobId)+'.root'
print("Output File: %s" % (options.outputFile))

process = cms.Process('FULLSIM', eras.Phase2)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
# process.load('Configuration.StandardSequences.DigiToRaw_cff')
# process.load('Configuration.StandardSequences.SimL1Emulator_cff')
# process.load('HLTrigger.Configuration.HLT_GRun_cff')
# process.load('Configuration.StandardSequences.RawToDigi_cff')
# process.load('Configuration.StandardSequences.ReconstructionCosmics_cff')
process.load('RecoLocalTracker.Configuration.RecoLocalTracker_cff') # pixelclusterizer process imported from here
process.load('SimTracker.TrackTriggerAssociation.TrackTriggerAssociator_cff')
process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
# process.load("Configuration.StandardSequences.SimulationRandomNumberGeneratorSeeds_cff")

# process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
# process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')

# needed for pixel RecHits (TkPixelCPERecord), from: https://github.com/cms-sw/cmssw/blob/master/EventFilter/SiPixelRawToDigi/test/runRawToClus_cfg.py
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.Geometry.GeometryExtended2026D63_cff') # specific for generation (full cms)
process.load('Configuration.Geometry.GeometryExtended2026D63Reco_cff') # specific for reconstruction (full cms)
print('Running with full CMS geometry & TkOnly Digitisation, Clustering')

#custom BRIL configs like geometry
#process.load('BRIL_ITsim.DataProductionTkOnly.cmsExtendedGeometry2026D999XML_cff')
# process.load('BRIL_ITsim.DataProductionTkOnly.TkOnlyDigiToRaw_cff')
# process.load('BRIL_ITsim.DataProductionTkOnly.TkOnlyRawToDigi_cff')
#print('Running with special BRIL Tk Only Geometry & TkOnly Digitisation, Clustering')

from L1Trigger.TrackTrigger.TTStubAlgorithmRegister_cfi import *
from L1Trigger.TrackTrigger.TTStub_cfi import *

# randomeze the seeds every time cmsRun is invoked
from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper
randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
randSvc.populate()
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.nEvents)
)

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

# define range of events to process in a single job (event numbering starts from 1 and range is limit inclusive)
n_events = 5000
range_min = int(options.jobId)*n_events + 1
range_max = min((int(options.jobId)+1)*n_events, options.nEvents)
print("Processing job {} with event numbers {}:{}".format(options.jobId, range_min, range_max))

# input source
process.source = cms.Source("PoolSource",
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            fileNames = cms.untracked.vstring(options.inputFile),
                            skipEvents = cms.untracked.uint32(0),
			    eventsToProcess = cms.untracked.VEventRange("1:{}-1:{}".format(range_min, range_max))
)

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
     # SelectEvents = cms.untracked.PSet(
         # SelectEvents = cms.vstring('generation_step')
     # ),
     dataset = cms.untracked.PSet(
       dataTier = cms.untracked.string('GEN-SIM-DIGI'),
       filterName = cms.untracked.string('')
     ),
     splitLevel = cms.untracked.int32(0),
     # outputCommands = cms.untracked.vstring( "keep *_*_*_*"),
     outputCommands = process.RECOSIMEventContent.outputCommands,
     fileName = cms.untracked.string(options.outputFile),
)


# process.RECOSIMEventContent.outputCommands.extend(
     # ['keep *_*_*_GEN'])
# process.RECOSIMEventContent.outputCommands.extend(
     # ['keep *_trackingtruthprod_*_*'])
# process.RECOSIMEventContent.outputCommands.extend(
     # ['keep *_mergedtruth_MergedTrackTruth_*'])
# process.RECOSIMEventContent.outputCommands.extend(
     # ['keep *_g4SimHits_*_*'])
# process.RECOSIMEventContent.outputCommands.extend(
     # ['keep *_TriggerResults_*_*'])
process.RECOSIMEventContent.outputCommands.append('keep  *_*_*_*')
process.RECOSIMEventContent.outputCommands.append('drop  *_*UnsuppressedDigis_*_*')
process.RECOSIMEventContent.outputCommands.append('drop  *_simEcal*_*_*')
process.RECOSIMEventContent.outputCommands.append('drop  *_simHcal*_*_*')
process.RECOSIMEventContent.outputCommands.append('drop  *_mix_*_*')
# process.RECOSIMEventContent.outputCommands.append('keep  *_*_*_*')
# process.RECOSIMEventContent.outputCommands.append('drop  *_mix_*_STUBS')
# process.RECOSIMEventContent.outputCommands.append('drop  PCaloHits_*_*_*')
# process.RECOSIMEventContent.outputCommands.append('drop  *_ak*_*_*')
# process.RECOSIMEventContent.outputCommands.append('drop  *_mix_*_*')
# process.RECOSIMEventContent.outputCommands.append('keep  *_simSi*_*_*')
# process.RECOSIMEventContent.outputCommands.append('drop  *_simEcal*_*_*')
# process.RECOSIMEventContent.outputCommands.append('drop  *_simHcal*_*_*')
# process.RECOSIMEventContent.outputCommands.append('drop  *_TkPixelCPERecord*_*_*')
# process.RECOSIMEventContent.outputCommands.append('keep  *_g4SimHits_Tracker*_*')
# process.RECOSIMEventContent.outputCommands.append('drop  *_g4SimHits_*_*')

# Additional output definition
import SimG4Core.Application.g4SimHits_cfi
process.g4SimHits.Generator.ApplyEtaCuts  = cms.bool(False)
process.g4SimHits.Generator.MinPCut  = cms.double(0.0001) #100keV
process.g4SimHits.Generator.BeamBkgdEvent = cms.untracked.bool(True)
process.g4SimHits.StackingAction.SaveFirstLevelSecondary = cms.untracked.bool(True)
process.g4SimHits.StackingAction.SavePrimaryDecayProductsAndConversionsInCalo = cms.untracked.bool(True)
process.g4SimHits.StackingAction.SavePrimaryDecayProductsAndConversionsInMuon = cms.untracked.bool(True)

#not present in CMSSW_11_2_X
#process.g4SimHits.SteppingAction.KillBeamPipe = cms.bool(False)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')
# Other statements

# Path and EndPath definitions
process.simulation_step   = cms.Path(process.psim*process.mix)
process.digitisationTkOnly_step = cms.Path(process.pdigi_valid)
process.PixelClusterizer_step = cms.Path(process.pixeltrackerlocalreco) # pixel clusterization process
process.L1TrackTrigger_step     = cms.Path(process.TrackTriggerClustersStubs)
process.L1TTAssociator_step     = cms.Path(process.TrackTriggerAssociatorClustersStubs)
# process.L1simulation_step = cms.Path(process.SimL1Emulator)
# process.digi2raw_step     = cms.Path(process.DigiToRaw)

# Path and EndPath definitions
# process.raw2digi_step       = cms.Path(process.RawToDigi)
# process.localreco_step      = cms.Path(process.localReconstructionCosmics)
# process.bhalo_step          = cms.Path(process.beamhaloTracksSeq)
# process.reco_track_step     = cms.Path(process.trackerlocalreco)
# process.reco_calo_step      = cms.Path(process.calolocalreco)
# process.reco_muon_step      = cms.Path(process.muonsLocalRecoCosmics)
process.endjob_step         = cms.Path(process.endOfProcess)
process.out_step            = cms.EndPath(process.output)
    


# Schedule definition
#process.schedule = cms.Schedule(process.raw2digi_step,process.reco_track_step,process.reco_calo_step,process.reco_muon_step,process.endjob_step,process.out_step)




# process.endjob_step       = cms.Path(process.endOfProcess)
# process.out_step = cms.EndPath(process.output)

# Schedule definition

# Simulation and L1

process.schedule = cms.Schedule(process.simulation_step,
                                process.digitisationTkOnly_step,
                                process.PixelClusterizer_step,
                                process.L1TrackTrigger_step,
                                process.L1TTAssociator_step)
                                # process.L1simulation_step,
                                # process.digi2raw_step)

# High level trigger

#process.schedule.extend(process.HLTSchedule)


# Reconstruction

# process.schedule.extend([process.raw2digi_step,process.localreco_step,process.bhalo_step])


process.schedule.extend([process.endjob_step,process.out_step])


#Setup FWK for multithreaded
process.options.numberOfThreads=cms.untracked.uint32(options.nThreads)
process.options.numberOfStreams=cms.untracked.uint32(options.nThreads)

# filter all path with the production filter sequence
# for path in process.paths:
	# getattr(process,path)._seq = process.generator * getattr(process,path)._seq

#do not add changes to your config after this point (unless you know what you are doing)
# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)

# Automatic addition of the customisation function
#customisation functions to only run Tracker Digitisation and Pixel Clustering
from BRIL_ITsim.DataProductionTkOnly.TkOnlyDigi_cff import TkOnlyDigi
process = TkOnlyDigi(process)
from BRIL_ITsim.DataProductionTkOnly.PixelClusterizerOnly_cff import PixelClusterizerOnly
process = PixelClusterizerOnly(process)
# End adding early deletion

