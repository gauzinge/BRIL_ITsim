# Auto generated configuration file
# using: 
# Revision: 1.163 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: Configuration/Generator/python/TTbar_cfi.py -s GEN,SIM,DIGI,L1,DIGI2RAW,HLT --eventcontent FEVTDEBUG --datatier GEN-SIM-RAW --conditions FrontierConditions_GlobalTag,IDEAL_V11::All -n 10 --no_exec
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing


# In the line below 'analysis' is an instance of VarParsing object 
options = VarParsing ('analysis')

# Here we have defined our own two VarParsing options 
# add a list of strings for events to process
options.register ('nEvents',
                                 1000,
                                 VarParsing.multiplicity.singleton,
                                 VarParsing.varType.int,
                  "The number of events to generate: 10")
options.register ('inputFile',
                  'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/bib_generation_fullgeo/halo/BeamHalo.0.root',
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
                  'file:/afs/cern.ch/work/p/pkicsiny/private/cmssw/CMSSW_11_2_0_pre6/src/BRIL_ITsim/BIBGeneration/test_simulator',
                  # 'file:/afs/cern.ch/work/g/gauzinge/public/',
                                 VarParsing.multiplicity.singleton,
                                 VarParsing.varType.string,
                  "The output directory")

options.parseArguments()
options.outputFile=options.outputDirectory+'/BeamHaloReco.'+str(options.jobId)+'.root'
print("Output File: %s" % (options.outputFile))

#edit GA: changed to new C11 era
# from Configuration.StandardSequences.Eras import eras
from Configuration.Eras.Era_Phase2C11_cff import Phase2C11
process = cms.Process('FULLSIM', Phase2C11)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
# process.load('SimGeneral.MixingModule.mix_POISSON_average_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D63Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D63_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
# process.load('Configuration.StandardSequences.Generator_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
# process.load('IOMC.EventVertexGenerators.VtxSmearedHLLHC_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
# process.load('RecoLocalTracker.Configuration.RecoLocalTracker_cff')
process.load('Configuration.StandardSequences.RecoSim_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#custom BRIL configs like Geometry
# process.load('BRIL_ITsim.DataProductionTkOnly.cmsExtendedGeometry2026D999XML_cff')
# process.load('Configuration.Geometry.GeometryExtended2026D63Reco_cff')
# process.load('BRIL_ITsim.DataProductionTkOnly.TkOnlyDigiToRaw_cff')
# process.load('BRIL_ITsim.DataProductionTkOnly.TkOnlyRawToDigi_cff')
# print 'Running with special BRIL Tk Only Geometry & TkOnly Digitisation, Clustering'

# from L1Trigger.TrackTrigger.TTStubAlgorithmRegister_cfi import *
# from L1Trigger.TrackTrigger.TTStub_cfi import *

# randomeze the seeds every time cmsRun is invoked
from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper
randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
randSvc.populate()

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.nEvents)
)

process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring(),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(

        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(1)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    makeTriggerResults = cms.obsolete.untracked.bool,
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(1),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Input source
process.source = cms.Source("PoolSource",
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            fileNames = cms.untracked.vstring(options.inputFile),
                            skipEvents = cms.untracked.uint32(0)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step1_step2_step3 nevts:'+str(options.nEvents)),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition
process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
     #SelectEvents = cms.untracked.PSet(
     #    SelectEvents = cms.vstring('generation_step')
     #),
     dataset = cms.untracked.PSet(
       dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW-RECO'),
       filterName = cms.untracked.string('')
     ),
     fileName = cms.untracked.string(options.outputFile),
     outputCommands = process.FEVTDEBUGEventContent.outputCommands,
     splitLevel = cms.untracked.int32(0),
)


process.FEVTDEBUGoutput.outputCommands.append('drop  *')
process.FEVTDEBUGoutput.outputCommands.append('keep  *_g4SimHits__*')
process.FEVTDEBUGoutput.outputCommands.append('keep  *_*_TrackerHitsPixelEndcapHighTof_*')
process.FEVTDEBUGoutput.outputCommands.append('keep  *_*_TrackerHitsPixelEndcapLowTof_*')
process.FEVTDEBUGoutput.outputCommands.append('keep  *_*_TrackerHitsPixelBarrelHighTof_*')
process.FEVTDEBUGoutput.outputCommands.append('keep  *_*_TrackerHitsPixelBarrelLowTof_*')
process.FEVTDEBUGoutput.outputCommands.append('keep  *_simSiPixelDigis_Pixel_*')
process.FEVTDEBUGoutput.outputCommands.append('keep  *_siPixelClusters_*_*')
process.FEVTDEBUGoutput.outputCommands.append('keep  *_siPixelRecHits_*_*')
process.FEVTDEBUGoutput.outputCommands.append('keep  *_SiPixelRecHits_*_*')
process.FEVTDEBUGoutput.outputCommands.append('keep  *_simSiPixelRecHits_*_*')
process.FEVTDEBUGoutput.outputCommands.append('keep  *_*Pixel*_*_*')
process.FEVTDEBUGoutput.outputCommands.append('keep  *_*_Pixel_*')
process.FEVTDEBUGoutput.outputCommands.append('drop  *_mix_Pixel_*')
process.FEVTDEBUGoutput.outputCommands.append('drop  *_simGmtStage2Digis_*_*')
process.FEVTDEBUGoutput.outputCommands.append('drop  *_simGtStage2Digis_*_*')
process.FEVTDEBUGoutput.outputCommands.append('drop  *_gmtStage2Digis_*_*')
process.FEVTDEBUGoutput.outputCommands.append('drop  *_gtStage2Digis_*_*')

# Additional output definition
import SimG4Core.Application.g4SimHits_cfi
process.g4SimHits.Generator.ApplyEtaCuts  = cms.bool(False)
process.g4SimHits.Generator.MinPCut  = cms.double(0.0001) #100keV
process.g4SimHits.Generator.BeamBkgdEvent = cms.untracked.bool(True)
process.g4SimHits.StackingAction.SaveFirstLevelSecondary = cms.untracked.bool(True)
process.g4SimHits.StackingAction.SavePrimaryDecayProductsAndConversionsInCalo = cms.untracked.bool(True)
process.g4SimHits.StackingAction.SavePrimaryDecayProductsAndConversionsInMuon = cms.untracked.bool(True)

process.mix.digitizers = cms.PSet(process.theDigitizersValid)
process.mix.digitizers.pixel.SSDigitizerAlgorithm.HitDetectionMode = cms.int32(2)
process.mix.digitizers.pixel.PixelDigitizerAlgorithm.ApplyTimewalk = cms.bool(False)
process.mix.digitizers.pixel.PixelDigitizerAlgorithm.ElectronPerAdc = cms.double(500.)

#not present in CMSSW_11_2_X
#process.g4SimHits.SteppingAction.KillBeamPipe = cms.bool(False)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T21', '')
# Other statements

# Path and EndPath definitions
process.simulation_step   = cms.Path(process.psim*process.mix)
process.digitisation_step = cms.Path(process.pdigi_valid)
process.L1TrackTrigger_step = cms.Path(process.L1TrackTrigger)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.recosim_step = cms.Path(process.recosim)
# process.PixelClusterizer_step = cms.Path(process.pixeltrackerlocalreco)
# process.digi2raw_step = cms.Path(process.siPixelRawData)#*process.SiStripDigiToRaw)
# process.raw2digi_step = cms.Path(process.RawToDigi_pixelOnly)
# process.PixelClusterizer_step = cms.Path(process.pixeltrackerlocalreco)
# process.L1TrackTrigger_step     = cms.Path(process.TrackTriggerClustersStubs)
# process.L1TTAssociator_step     = cms.Path(process.TrackTriggerAssociatorClustersStubs)
# process.L1simulation_step = cms.Path(process.SimL1Emulator)
# process.digi2raw_step     = cms.Path(process.DigiToRaw)

# Path and EndPath definitions
# process.raw2digi_step       = cms.Path(process.RawToDigi)
# process.localreco_step      = cms.Path(process.localReconstructionCosmics)
# process.bhalo_step          = cms.Path(process.beamhaloTracksSeq)
# process.reco_track_step     = cms.Path(process.trackerlocalreco)
# process.reco_calo_step      = cms.Path(process.calolocalreco)
# process.reco_muon_step      = cms.Path(process.muonsLocalRecoCosmics)
process.endjob_step         = cms.EndPath(process.endOfProcess)
process.out_step            = cms.EndPath(process.FEVTDEBUGoutput)
    


# Schedule definition
#process.schedule = cms.Schedule(process.raw2digi_step,process.reco_track_step,process.reco_calo_step,process.reco_muon_step,process.endjob_step,process.out_step)




# process.endjob_step       = cms.Path(process.endOfProcess)
# process.out_step = cms.EndPath(process.output)

# Schedule definition

# Simulation and L1

process.schedule = cms.Schedule(process.simulation_step,
                                process.digitisation_step,
                                process.L1TrackTrigger_step,
                                process.L1simulation_step,
                                process.digi2raw_step,
                                process.raw2digi_step,
                                process.reconstruction_step,
                                process.recosim_step,
                                )

# High level trigger
#process.schedule.extend(process.HLTSchedule)

process.schedule.extend([process.endjob_step,process.out_step])


from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC 

#call to customisation function customizeHLTforMC imported from HLTrigger.Configuration.customizeHLTforMC
process = customizeHLTforMC(process)

# Automatic addition of the customisation function from SimGeneral.MixingModule.fullMixCustomize_cff
from SimGeneral.MixingModule.fullMixCustomize_cff import setCrossingFrameOn 

#call to customisation function setCrossingFrameOn imported from SimGeneral.MixingModule.fullMixCustomize_cff
process = setCrossingFrameOn(process)

# End of customisation functions

#Setup FWK for multithreaded
process.options.numberOfThreads=cms.untracked.uint32(options.nThreads)
process.options.numberOfStreams=cms.untracked.uint32(options.nThreads)

# filter all path with the production filter sequence
# for path in process.paths:
	# getattr(process,path).insert(0, process.generator)

# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)

# Automatic addition of the customisation function
#customisation functions to only run Tracker Digitisation and Pixel Clustering
# from BRIL_ITsim.DataProductionTkOnly.TkOnlyDigi_cff import TkOnlyDigi
# process = TkOnlyDigi(process)
# from BRIL_ITsim.DataProductionTkOnly.PixelClusterizerOnly_cff import PixelClusterizerOnly
# process = PixelClusterizerOnly(process)
# End adding early deletion

