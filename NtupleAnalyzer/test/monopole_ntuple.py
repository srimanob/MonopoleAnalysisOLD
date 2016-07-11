import FWCore.ParameterSet.Config as cms

process = cms.Process("ntuple")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load('Configuration/Geometry/GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")
process.load("RecoTracker.Configuration.RecoTracker_cff")
#process.load("RecoTracker.TransientTrackingRecHit.TransientTrackingRecHitBuilderWithoutRefit_cfi")
#process.load("TrackingTools.TrackRefitter.TracksToTrajectories_cff")
process.load('Configuration/StandardSequences/EndOfProcess_cff')


process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True) 
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1),
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

# reduce verbosity
#process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10000)  

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'root://eoscms.cern.ch//eos/cms/store/user/srimanob/monopole/13TeV/v3-spinhalf/monopole_1000GeV/step2_AODSIM_1000_100.root'
    )
)

# Global tag
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_v14', '')

process.ntupleMaker = cms.EDAnalyzer(
    "NtupleAnalyzer",
    outputFile   = cms.string('outputfile.root'), 
    genInfo = cms.InputTag('genParticles'),
    generalTracks = cms.InputTag('generalTracks')
)

process.p = cms.Path(process.ntupleMaker)

 
