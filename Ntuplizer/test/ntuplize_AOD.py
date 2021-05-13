import FWCore.ParameterSet.Config as cms

process = cms.Process('Demo')

process.load('FWCore.MessageService.MessageLogger_cfi')


process.MessageLogger.cerr.FwkReport.reportEvery = 200 

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(400) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.options = cms.untracked.PSet( numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(1),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True),

   )



process.source = cms.Source("PoolSource",
     duplicateCheckMode=cms.untracked.string("noDuplicateCheck"),
    fileNames = cms.untracked.vstring(
   'file:/afs/cern.ch/work/a/athachay/public/BsToMuMuGamma/RunIIAutumn18DRPremix/BsToMuMuGamma_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-evtgen-pythia8/AODSIM/102X_upgrade2018_realistic_v15-v1/606765BD-9BB4-9741-925C-A0C69B933039.root',      
   #'file:DoublePhotonGun/DoublePhoton0To40FlatPtAODSIM_HI_Reco_1.root',      
    )
)
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("muonNtuplizer.root")
)

process.effMeasure = cms.EDAnalyzer("NTuplizer",
	muons=cms.InputTag("muons"),
	beamSpot = cms.InputTag("offlineBeamSpot"),
	vertices = cms.InputTag("offlinePrimaryVertices")
)


process.p = cms.Path(process.effMeasure)
