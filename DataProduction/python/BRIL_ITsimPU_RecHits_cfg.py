# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
##STEP1
# with command line options: SingleNuE10_cfi.py --mc --conditions auto:phase2_realistic -n ${NEVENTS} --era Phase2 --eventcontent ${EVENTCONTENT} --relval 25000,250 -s GEN,SIM --datatier GEN-SIM --beamspot HLLHC --geometry Extended2023D21 --fileout file:step1.root --nThreads ${NTHREADS}
##STEP2
# with command line options: step2 --mc --conditions auto:phase2_realistic --pileup_input file:${PUFILE} -n ${NEVENTS} --era Phase2 --eventcontent ${EVENTCONTENT} -s DIGI:pdigi_valid,L1,DIGI2RAW --datatier GEN-SIM-DIGI-RAW --pileup ${PUSTRING} --geometry Extended2023D21 --filein file:step1.root --fileout file:step2.root --nThreads ${NTHREADS}
##STEP3
# with command line options: step3 --mc --conditions auto:phase2_realistic --pileup_input file:${PUFILE} -n ${NEVENTS} --era Phase2 --eventcontent ${EVENTCONTENT} --runUnscheduled -s RAW2DIGI,L1Reco,RECO,RECOSIM --datatier GEN-SIM-RECO --pileup ${PUSTRING} --geometry Extended2023D21 --filein file:step2.root --fileout file:step3_pixel_PU_${PU}.${JOBID}.root --nThreads ${NTHREADS}

import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

# In the line below 'analysis' is an instance of VarParsing object 
options = VarParsing ('analysis')

# Here we have defined our own two VarParsing options 
# add a list of strings for events to process

nevents = 10
pu = 10
analyzer_plugin = "PuRecHitDigi_cff"

options.register ('nEvents',
                                 10,
                                 VarParsing.multiplicity.singleton,
                                 VarParsing.varType.int,
                  "The number of events to simulate: {}".format(nevents))
options.register ('pileupFile', ('file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_0.root',
				 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_10.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_11.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_12.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_13.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_14.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_15.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_16.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_17.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_18.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_19.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_1.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_20.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_21.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_22.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_23.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_24.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_25.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_26.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_27.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_28.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_29.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_2.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_30.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_31.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_32.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_33.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_34.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_35.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_36.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_37.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_38.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_39.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_3.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_40.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_41.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_42.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_43.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_44.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_45.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_46.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_47.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_48.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_49.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_4.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_50.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_51.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_52.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_53.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_54.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_55.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_56.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_57.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_58.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_59.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_5.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_6.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_7.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_8.root',
                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias5k_9.root'
),
#                                 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/minbias_fullgeo/minBias300k_0.root',
                                 VarParsing.multiplicity.list,
                                 VarParsing.varType.string,
                                 "File with Minimum Bias events to use as PU overlay")
options.register ('pileupAverage',
                                 10,
                                 VarParsing.multiplicity.singleton,
                                 VarParsing.varType.float,
                  "The average PU number: {}".format(pu))
options.register ('bunchSpace',
                                 25,
                                 VarParsing.multiplicity.singleton,
                                 VarParsing.varType.int,
                  "The bunch spacing: 25 ns")
options.register ('minBunch',
                                 -12,
                                 VarParsing.multiplicity.singleton,
                                 VarParsing.varType.int,
                  "The minimum bunch: -12 BX")
options.register ('maxBunch',
                                 3,
                                 VarParsing.multiplicity.singleton,
                                 VarParsing.varType.int,
                  "The maximum bunch: 3 BX")
# options.register ('eventContent',
                                 # 'FEVTDEBUG',
                                 # VarParsing.multiplicity.singleton,
                                 # VarParsing.varType.string,
                  # "The Event content string: FEVTDEBUG")
options.register ('nThreads',
                                 10,
                                 VarParsing.multiplicity.singleton,
                                 VarParsing.varType.int,
                  "The number of threads to use: 1")
options.register ('jobId',
                                 0,
                                 VarParsing.multiplicity.singleton,
                                 VarParsing.varType.int,
                  "The job Id: 0")
options.register ('outputDirectory',
                  'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/pu_simulations_fullgeo/low_stat_examples',
                                 VarParsing.multiplicity.singleton,
                                 VarParsing.varType.string,
                  "The output directory")

options.parseArguments()

# define range of events to process in a single job (event numbering starts from 1 and range is limit inclusive)
n_events = 10
range_min = int(options.jobId)*n_events + 1
range_max = min((int(options.jobId)+1)*n_events, options.nEvents)
print("Processing job {} with event numbers {}:{}".format(options.jobId, range_min, range_max))

options.outputFile=options.outputDirectory+'/step3_pixel_PU_'+str(options.pileupAverage)+'_'+str(n_events)+'.'+str(options.jobId)+'_{}_notimewalk.root'.format(analyzer_plugin)
print("Output File: %s" % (options.outputFile))
print("PU: {}\nNum. events: {}".format(options.pileupAverage, options.nEvents))
# process = cms.Process('FULLSIM',eras.Phase2)
#mod GA:
# from Configuration.StandardSequences.Eras import eras
from Configuration.Eras.Era_Phase2C11_cff import Phase2C11
from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper

process = cms.Process('FULLSIM',Phase2C11)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
# process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('SimGeneral.MixingModule.mix_POISSON_average_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D63Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D63_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedHLLHC_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('HLTrigger.Configuration.HLT_Fake2_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
# process.load('RecoLocalTracker.Configuration.RecoLocalTracker_cff')
process.load('Configuration.StandardSequences.RecoSim_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# randomeze the seeds every time cmsRun is invoked
randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
randSvc.populate()

# set up the event number (change to options.nEvents if no manual event cap required)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(n_events)
)

# Input source (frist event set manually)
process.source = cms.Source("EmptySource",
			    firstEvent = cms.untracked.uint32(range_min),
)

# process options - edit GA
process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring("ProductNotFound"),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(

        ),
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

# process.options = cms.untracked.PSet(
    # SkipEvent = cms.untracked.vstring('ProductNotFound')
# )

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step1_step2_step3 nevts:'+str(options.nEvents)),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition
process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    #edit GA: do we need SelectEvents?
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW-RECO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string(options.outputFile),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition
# include the filtering step
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
# process.FEVTDEBUGoutput.outputCommands.append('keep  *')
process.FEVTDEBUGoutput.outputCommands.append('drop  *_simGmtStage2Digis_*_*')
process.FEVTDEBUGoutput.outputCommands.append('drop  *_simGtStage2Digis_*_*')
process.FEVTDEBUGoutput.outputCommands.append('drop  *_gmtStage2Digis_*_*')
process.FEVTDEBUGoutput.outputCommands.append('drop  *_gtStage2Digis_*_*')

# process.FEVTDEBUGoutput.outputCommands.append('drop  *_*_*')
# print(process.FEVTDEBUGoutput.outputCommands)

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T21', '')

# generator setup
process.generator = cms.EDProducer("FlatRandomEGunProducer",
    AddAntiParticle = cms.bool(False),
    PGunParameters = cms.PSet(
        MaxE = cms.double(10.01),
        MaxEta = cms.double(2.5),
        MaxPhi = cms.double(3.14159265359),
        MinE = cms.double(9.99),
        MinEta = cms.double(-2.5),
        MinPhi = cms.double(-3.14159265359),
        PartID = cms.vint32(12)
    ),
    Verbosity = cms.untracked.int32(0),
    firstRun = cms.untracked.uint32(1),
    psethack = cms.string('single Nu E 10')
)

# mixing module
process.mix.input.nbPileupEvents.averageNumber = cms.double(options.pileupAverage)
process.mix.bunchspace = cms.int32(options.bunchSpace)
process.mix.minBunch = cms.int32(options.minBunch)
process.mix.maxBunch = cms.int32(options.maxBunch)
# process.mix.seed = cms.int32(@SEED@)
# process.mix.input.fileNames = cms.untracked.vstring([options.pileupFile])
process.mix.input.fileNames = cms.untracked.vstring(*options.pileupFile)

#edit GA: needed?
# process.mix.playback = True
# process.mix.digitizers = cms.PSet()
# for a in process.aliases: delattr(process, a)

process.mix.digitizers = cms.PSet(process.theDigitizersValid)
process.mix.digitizers.pixel.PixelDigitizerAlgorithm.ElectronPerAdc = cms.double(500.) # change digi charge granularity
#process.mix.digitizers.pixel.PixelDigitizerAlgorithm.AdcFullScale = cms.int32(5000)
#process.mix.digitizers.pixel.PixelDigitizerAlgorithm.ThresholdInElectrons_Endcap = cms.double(5000.0)
process.mix.digitizers.pixel.SSDigitizerAlgorithm.HitDetectionMode = cms.int32(2)
process.mix.digitizers.pixel.PixelDigitizerAlgorithm.ApplyTimewalk = cms.bool(False) # toggle timewalk curve overlay
# keep simhits
print("keeping psimhits")
process.mix.mixObjects.mixSH.crossingFrames.extend([ 
        'TrackerHitsPixelEndcapHighTof', 
        'TrackerHitsPixelEndcapLowTof'])

process.mix.digitizers.mergedtruth.select.signalOnlyTP = cms.bool(False)

# root file with ntuple
ntuple_filename = '/afs/cern.ch/work/p/pkicsiny/private/cmssw/CMSSW_11_2_0_pre6/src/BRIL_ITsim/DataProduction/psimhits/PSimHits_ntuple_{}.root'.format(options.jobId)
process.TFileService = cms.Service("TFileService", fileName = cms.string(ntuple_filename), closeFileFast = cms.untracked.bool(True))

#    12 occurrences of key PSimHitCrossingFrame + "mix" + "g4SimHitsTrackerHitsPixelEndcapHighTof" "FULLSIM"
#    12 occurrences of key PSimHitCrossingFrame + "mix" + "g4SimHitsTrackerHitsPixelEndcapLowTof" "FULLSIM"

process.dumpPSimHits = cms.EDAnalyzer(analyzer_plugin,
				simhits_low=cms.InputTag("g4SimHits", "TrackerHitsPixelEndcapLowTof", "FULLSIM"),
                                simhits_high=cms.InputTag("g4SimHits", "TrackerHitsPixelEndcapHighTof", "FULLSIM"),
                                crossingframe_simhits_low=cms.InputTag("mix", "g4SimHitsTrackerHitsPixelEndcapLowTof", "FULLSIM"),
                                crossingframe_simhits_high=cms.InputTag("mix", "g4SimHitsTrackerHitsPixelEndcapHighTof", "FULLSIM"),
				digis=cms.InputTag("simSiPixelDigis", "Pixel", "FULLSIM"),
	  			rechits=cms.InputTag("siPixelRecHitsPreSplitting"),
			        associateHitbySimTrack = cms.bool(True),
				associatePixel = cms.bool(True),
                                associateStrip = cms.bool(False),
                                associateRecoTracks = cms.bool(False),
 		  	        ROUList = cms.vstring('g4SimHitsTrackerHitsPixelEndcapLowTof',
						      'g4SimHitsTrackerHitsPixelEndcapHighTof',
                                                               'TrackerHitsPixelEndcapLowTof',
                                                               'TrackerHitsPixelEndcapHighTof'),
    	   			usePhase2Tracker = cms.bool(True),
		  	        pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
			        nQBins=cms.untracked.uint32(100),
                                maxQBin=cms.untracked.uint32(100000),
				jobid=cms.untracked.uint32(options.jobId),
)

process.dump = cms.EDAnalyzer("EventContentAnalyzer")

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.digitisation_step = cms.Path(process.pdigi_valid)
process.L1TrackTrigger_step = cms.Path(process.L1TrackTrigger)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.recosim_step = cms.Path(process.recosim)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

# for saving simhits
process.psimhitdump_step =  cms.Path(process.dumpPSimHits)
process.dump_step = cms.Path(process.dump)

# process.digi2raw_step = cms.Path(process.siPixelRawData)#*process.SiStripDigiToRaw)
# process.raw2digi_step = cms.Path(process.RawToDigi_pixelOnly)
# process.reconstruction_step = cms.Path(process.pixeltrackerlocalreco)
# process.reconstruction_step = cms.Path(process.trackerlocalreco)


# Schedule definition
process.schedule = cms.Schedule(process.generation_step,
                                process.genfiltersummary_step,
       				#process.psimhitdump_step, # before out_step dump psimhits into a root file
                                process.simulation_step,
				#process.dump_step,
				#process.psimhitdump_step, # before out_step dump psimhits into a root file
                                process.digitisation_step, # mixing is done here
                                process.psimhitdump_step, # before out_step dump psimhits into a root file
                                process.L1TrackTrigger_step,
                                process.L1simulation_step,
                                process.digi2raw_step,
                                process.raw2digi_step,
                                process.reconstruction_step,
                                process.recosim_step,
                                #process.psimhitdump_step # before out_step dump psimhits into a root file
				)
                                # process.L1Reco_step,

# process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.endjob_step,process.FEVTDEBUGoutput_step])

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
for path in process.paths:
	getattr(process,path).insert(0, process.generator)

# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
