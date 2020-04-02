# imports
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

analyzer = "ITclusterAnalyzer"
#analyzerCode = "ITclusterAnalyzer_hits"
analyzerCode = "ITclusterAnalyzer_ORIGINAL"
#analyzerCode = "ITclusterAnalyzer"


# create a new CMS process
#process = cms.Process("ITclusterAnalyzer_hits")
process = cms.Process(analyzer)

print("--------------------------------1")
# set up the options
options = VarParsing.VarParsing('analysis')
#set up the defaults
options.inputFiles = 'file:/afs/cern.ch/user/c/corderom/public/PUdata/step3_pixel_PU_1.0.0.root'
#options.inputFiles = 'file:/eos/user/g/gauzinge/PUdata/step3_pixel_PU_1.0.0.root'
options.outputFile='/afs/cern.ch/user/c/corderom/private/mySimDir/CMSSW_10_4_0_pre2/src/BRIL_ITsim/AnalysisOut/summary_runSim.root'
#options.inputFiles = 'file:/afs/cern.ch/user/c/corderom/public/PUdata/step3_pixel_PU_10.0.root'
#options.outputFile='/afs/cern.ch/user/c/corderom/private/mySimDir/CMSSW_10_4_0_pre2/src/BRIL_ITsim/AnalysisOut/summary.root'
options.maxEvents = 200 #all events
#options.maxEvents = 5 #all events
#options.maxEvents = -1 #all events
#options.maxEvents = 100000 #all events

#get and parse command line arguments
options.parseArguments()

# load the geomtry that i modified
process.load('Configuration.Geometry.GeometryExtended2023D21Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')

print("--------------------------------9")
# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append(analyzerCode)
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit=cms.untracked.int32(-1)
)
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(False)
                                    # ,SkipEvent = cms.untracked.vstring('ProductNotFound')
                                    )

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))

# the input file
process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring(options.inputFiles)
                            )

process.content = cms.EDAnalyzer("EventContentAnalyzer")
# the config of my analyzer
print("--------------------------------10")
process.BRIL_IT_Analysis = cms.EDAnalyzer(analyzerCode,
                                         clusters=cms.InputTag("siPixelClusters"),
                                         simlinks=cms.InputTag("simSiPixelDigis", "Pixel", "FULLSIM"),
                                         digis=cms.InputTag("simSiPixelDigis", "Pixel", "FULLSIM"),
                                         # simlinks=cms.InputTag("Pixel"),
                                         # simtracks=cms.InputTag("g4SimHits"),
                                         maxBin=cms.untracked.uint32(5000),
                                         docoincidence=cms.untracked.bool(True),
                                         dx_cut=cms.double(.1),
                                         dy_cut=cms.double(.1),
                                         dz_cut=cms.double(0.9)
                                         )

print("--------------------------------20")
# the TFIleService that produces the output root files
process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string(options.outputFile)
                                   )


# process.p = cms.Path( ... process.content * ...  )
process.p = cms.Path(process.BRIL_IT_Analysis)
