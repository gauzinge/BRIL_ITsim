# imports
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os

# create a new CMS process
plugin_file = "BIBAnalyzerOccupancy"
process = cms.Process(plugin_file)

# set up the options
options = VarParsing.VarParsing('analysis')

# input file(s)
bib_type = "gas_oxygen"

#BIB
#input_list = "file:/afs/cern.ch/work/p/pkicsiny/private/cmssw/CMSSW_11_2_0_pre6/src/BRIL_ITsim/BIBGeneration/test_simulator/BeamHaloReco.0.root"

#PU central samples
#input_list = "file:root://cmsxrootd.fnal.gov///store/relval/CMSSW_10_6_0_patch2/RelValNuGun/GEN-SIM-RECO/PU25ns_106X_upgrade2023_realistic_v3_2023D42PU100-v2/10000/1FBE9504-4EF5-064C-A3C3-05A3891A6079.root"

#BIB full stats
input_path = '/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/bib_simulations_fullgeo/{}'.format(bib_type)
input_list = [os.path.join("file:", input_path[1:], input_file) for input_file in os.listdir(input_path)  if ".root" in input_file]
options.inputFiles = input_list

print("Input file: ", input_list)
# output file
options.outputFile='results/occupancy/{}_{}_summary.root'.format(bib_type, plugin_file.lower())

# proccess this many events from input (-1 means all events)
options.maxEvents = -1

#get and parse command line arguments
options.parseArguments()

# load standard Geometry
process.load('Configuration.Geometry.GeometryExtended2026D63Reco_cff')

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append(plugin_file)
#process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#    limit=cms.untracked.int32(-1)
#)
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(False),
 	   	 		     FailPath = cms.untracked.vstring('ProductNotFound'),
				     SkipEvent = cms.untracked.vstring('EventCorrupt'))
#				     IgnoreCompletely = cms.untracked.vstring("DataCorrupt"))
#
# set number of events to process
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))

# set input file(s)
process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring(options.inputFiles)
                            )

# set EDAnalyzer type
process.content = cms.EDAnalyzer("EventContentAnalyzer")

# set the python config of my analyzer
process.BRIL_IT_Analysis = cms.EDAnalyzer(plugin_file,
					 clusters=cms.InputTag("siPixelClustersPreSplitting"), # RECO: np preSplitting
                                         digis=cms.InputTag("simSiPixelDigis", "Pixel", "FULLSIM"), # HLT
                                         simhits_low=cms.InputTag("g4SimHits", "TrackerHitsPixelEndcapLowTof", "FULLSIM"), # SIM: empty
                                         simhits_high=cms.InputTag("g4SimHits", "TrackerHitsPixelEndcapHighTof", "FULLSIM"),
 				 	 simlinks=cms.InputTag("simSiPixelDigis", "Pixel", "FULLSIM"), 
	  				 rechits=cms.InputTag("siPixelRecHitsPreSplitting"),
			                 associateHitbySimTrack = cms.bool(True),
					 associatePixel = cms.bool(True),
                                         associateStrip = cms.bool(False),
                                         associateRecoTracks = cms.bool(False),
 					 ROUList = cms.vstring('TrackerHitsPixelBarrelLowTof',
                                                               'TrackerHitsPixelBarrelHighTof',
                                                               'TrackerHitsPixelEndcapLowTof',
                                                               'TrackerHitsPixelEndcapHighTof'),
    	   				 usePhase2Tracker = cms.bool(True),
					 pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
					 nQBins=cms.untracked.uint32(100),
                                         maxQBin=cms.untracked.uint32(100000),
					 bib_type=cms.untracked.string(bib_type),
)

# the TFIleService that produces the output root files
process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string(options.outputFile)
                                   )


# process.p = cms.Path( ... process.content * ...  )
process.p = cms.Path(process.BRIL_IT_Analysis)
