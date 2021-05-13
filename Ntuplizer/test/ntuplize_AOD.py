import FWCore.ParameterSet.Config as cms

process = cms.Process('Demo')

process.load('FWCore.MessageService.MessageLogger_cfi')


process.MessageLogger.cerr.FwkReport.reportEvery = 200 

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(400) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.options = cms.untracked.PSet( numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(1),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(4),
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True),

   )



process.source = cms.Source("PoolSource",
     duplicateCheckMode=cms.untracked.string("noDuplicateCheck"),
    fileNames = cms.untracked.vstring(
   'file:',      
   #'file:DoublePhotonGun/DoublePhoton0To40FlatPtAODSIM_HI_Reco_1.root',      
    )
)
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("DoublePhotonGun/DoublePhoton0To40FlatPtAODSIM_EffMesure.root")
    #fileName = cms.string("DoublePhotonGun/DoublePhoton0To40FlatPtAODSIM_EffMesure_HI_lowPt.root")
)

process.effMeasure = cms.EDAnalyzer("RecoEfficiencyAODAnalyzer",
    GenParticleSrc  =cms.InputTag("genParticles"),
#    SuperClusterPhotonSrc = cms.InputTag("particleFlowEGamma")
      GedPhotonSrc    =cms.InputTag("gedPhotons"),
    SuperClusterSrc = cms.string("particleFlowSuperClusterECAL"),
    RefinedSuperClusterSrc = cms.string("particleFlowEGamma")
)


process.p = cms.Path(process.effMeasure)
