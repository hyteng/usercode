import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")


#process.load("Geometry.MuonCommonData.muonIdealGeometryXML_cfi");
#process.load("Geometry.RPCGeometry.rpcGeometry_cfi");
#process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi");
#process.load("CalibMuon.Configuration.Muon_FakeAlignment_cff");
#process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi");

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'DESIGN41X_V0::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:muonseed.root'
        #'-input-'
    )
)

process.demo = cms.EDAnalyzer('RPCSeedValidator',
    theRootFileName = cms.untracked.string('RPCSeedValidator.root'),
    SimHitTag = cms.untracked.InputTag("g4SimHits", "MuonRPCHits"),
    SimTrackTag = cms.untracked.InputTag("g4SimHits"),
    RPCRecHitTag = cms.untracked.InputTag("rpcRecHits"),
    RPCDigiSimLinkTag = cms.untracked.InputTag("simMuonRPCDigis", "RPCDigiSimLink"),
    TrajectorySeedCollectionTag = cms.untracked.InputTag("myRPCSeed", "GoodSeeds"),
    RecHitNumberTH = cms.untracked.uint32(5),
    CodeTH = cms.untracked.uint32(15),
    unCodeTH = cms.untracked.uint32(48),
    FilterType = cms.untracked.int32(1),
    SeedPurityTH = cms.untracked.double(0.5),
    isVertexConstraint = cms.untracked.bool(True),
    SegmentBendingPhiTH = cms.untracked.double(4.0),
    MaxBendingPhiTH = cms.untracked.double(0.09),
    debug = cms.untracked.bool(True)
)

process.p = cms.Path(process.demo)
