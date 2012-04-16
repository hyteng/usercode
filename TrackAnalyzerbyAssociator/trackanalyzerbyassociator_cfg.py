import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'DESIGN41X_V0::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:recSTA415_Mix_VertexnoFilternoOverlap_Pt1.0-100.0Gev_key.root'
    )
)

process.STAAna = cms.EDAnalyzer('TrackAnalyzerbyAssociator',
    #recTrackTag = cms.untracked.InputTag("standAloneMuons", "UpdatedAtVtx", "myRPC"),
    #recTrackTag = cms.untracked.InputTag("standAloneMuonsRPCSeed", "UpdatedAtVtx", "myRPC"),
    #recTrackTag = cms.untracked.InputTag("standAloneMuonsRPCSeedOnlyRPC", "UpdatedAtVtx", "myRPC"),
    #recTrackTag = cms.untracked.InputTag("standAloneMuonsMergedAncientandRPCSeed", "UpdatedAtVtx", "myRPC"),
    #recTrackTag = cms.untracked.InputTag("standAloneMuonsOld", "UpdatedAtVtx", "myRPC"),
    #recTrackTag = cms.untracked.InputTag("standAloneMuonsOldRPCSeed", "UpdatedAtVtx", "myRPC"),  
    #recTrackTag = cms.untracked.InputTag("standAloneMuonsOldRPCSeedOnlyRPC", "UpdatedAtVtx", "myRPC"),
    #recTrackTag = cms.untracked.InputTag("standAloneMuonsOldMergedAncientandRPCSeed", "UpdatedAtVtx", "myRPC"),
    #recTrackTag = cms.untracked.InputTag("standAloneSETMuons", "UpdatedAtVtx", "myRPC"),
    #recTrackTag = cms.untracked.InputTag("standAloneSETMuonsMergedSETMuonandRPCSeed", "UpdatedAtVtx", "myRPC"),
    recTrackTag = cms.untracked.InputTag("-branch-", "UpdatedAtVtx", "myRPC"),
    trackingParticleTag = cms.untracked.InputTag("mergedtruth","MergedTrackTruth"),
    PropagatorName = cms.string("SteppingHelixPropagatorAlong"),
    useTracker = cms.untracked.bool(False),
    useMuon = cms.untracked.bool(True),
    recTrackPurityThreshold = cms.untracked.double(0.5),
    debug = cms.untracked.bool(True),
    theRootFileName = cms.untracked.string("anaSTA415_-branch-_key.root"),
    CSClinksTag = cms.InputTag("simMuonCSCDigis","MuonCSCStripDigiSimLinks"),
    CSCwireLinksTag = cms.InputTag("simMuonCSCDigis","MuonCSCWireDigiSimLinks"),
    crossingframe  = cms.bool(False),
    CSCsimHitsTag = cms.InputTag("g4SimHits","MuonCSCHits"),
    CSCsimHitsXFTag = cms.InputTag("mix","g4SimHitsMuonCSCHits"),
    DTdigisimlinkTag = cms.InputTag("simMuonDTDigis"),
    associatorByWire  = cms.bool(False),
    links_exist = cms.bool(True),
    DTsimhitsTag  = cms.InputTag("g4SimHits","MuonDTHits"),
    DTsimhitsXFTag  = cms.InputTag("mix","g4SimHitsMuonDTHits"),
    dumpDT = cms.bool(False),
    DTdigiTag = cms.InputTag("simMuonDTDigis"),
    DTrechitTag = cms.InputTag("dt1DRecHits"),
    RPCdigisimlinkTag = cms.InputTag("simMuonRPCDigis","RPCDigiSimLink"),
    RPCsimhitsTag  = cms.InputTag("g4SimHits","MuonRPCHits"),
    RPCsimhitsXFTag = cms.InputTag("mix","g4SimHitsMuonRPCHits")
)


process.p = cms.Path(process.STAAna)
