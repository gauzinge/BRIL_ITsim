# imports
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os

# create a new CMS process
plugin_file = "BIBAnalyzerTof"
process = cms.Process(plugin_file)

# set up the options
options = VarParsing.VarParsing('analysis')

# input file(s)
bib_type = "gas_hydrogen"

#input_path = "file:/afs/cern.ch/work/p/pkicsiny/private/cmssw/CMSSW_11_2_0_pre6/src/BRIL_ITsim/DataProductionTkOnly/pkicsiny_pileup/step3_pixel_PU_200.0.0TkOnly.root"
input_path = '/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/bib_simulations/{}'.format(bib_type)
input_list = [os.path.join("file:", input_path[1:], input_file) for input_file in os.listdir(input_path)]
options.inputFiles = input_list[0]

# output file
options.outputFile='results/{}_tof_summary.root'.format(bib_type)

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
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit=cms.untracked.int32(-1)
)
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(False))

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
					 clusters=cms.InputTag("siPixelClustersPreSplitting"),
                                         digis=cms.InputTag("simSiPixelDigis", "Pixel", "FULLSIM"),
                                         simhits_low=cms.InputTag("g4SimHits", "TrackerHitsPixelEndcapLowTof", "FULLSIM"),
                                         simhits_high=cms.InputTag("g4SimHits", "TrackerHitsPixelEndcapHighTof", "FULLSIM"),
					 simlinks=cms.InputTag("simSiPixelDigis", "Pixel", "FULLSIM"),
					 nEBins=cms.untracked.uint32(10000),
                                         maxEBin=cms.untracked.uint32(100),
                                         )

# the TFIleService that produces the output root files
process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string(options.outputFile)
                                   )


# process.p = cms.Path( ... process.content * ...  )
process.p = cms.Path(process.BRIL_IT_Analysis)
