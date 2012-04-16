import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")

#process.GlobalTag.globaltag = 'IDEAL_21X::All'
process.GlobalTag.globaltag = 'DESIGN_36_V10::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:myfile.root'
    )
)

process.demo = cms.EDAnalyzer('RPCBendingAnalyzer',
    simHitTag = cms.untracked.InputTag("g4SimHits", "MuonRPCHits"),
    simTrackTag = cms.untracked.InputTag("g4SimHits"),
    RPCrecHitTag = cms.untracked.InputTag("rpcRecHits"),
    RPCDigiSimLinkTag = cms.untracked.InputTag("simMuonRPCDigis", "RPCDigiSimLink"),
    RPCLayer = cms.vuint32(0,1,2,3,4,5),
    codeTH = cms.untracked.uint32(63),
    debug = cms.untracked.bool(True),
    theRootFileName = cms.untracked.string("RPCValidationTree.root")
)


process.p = cms.Path(process.demo)
