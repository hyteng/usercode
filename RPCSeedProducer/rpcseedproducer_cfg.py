import FWCore.ParameterSet.Config as cms

process = cms.Process("RPCSeed")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.GlobalTag.globaltag = 'DESIGN_36_V10::All'


process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    #fileNames = cms.untracked.vstring('file:muongun.root')
    fileNames = cms.untracked.vstring('rfio:/castor/cern.ch/user/h/hyteng/mc/muongun/muongun_362_DESIGN36V10_1pair_3.8T_Pt3.0-1000.0Gev_Eta-1.0-1.0/muongun_362_DESIGN36V10_1pair_3.8T_Pt3.0-1000.0Gev_Eta-1.0-1.0_key.0001.root')
)

process.myRPCSeed = cms.EDProducer('RPCSeedProducer', 
        BarrelLayerRange = cms.vuint32(5),
        EndcapLayerRange = cms.vuint32(0),
        isCosmic = cms.bool(False),
        isSpecialLayers = cms.bool(False), 
        isMixBarrelwithEndcap = cms.bool(False),
        LayersinBarrel = cms.vuint32(1,1,1,1,0,1),
        LayersinEndcap = cms.vuint32(1,1,1,1,1,1),
        ConstraintedBarrelLayer = cms.vuint32(1,1,1,1,0,0),
        ConstraintedNegativeEndcapLayer = cms.vuint32(1,1,1),
        ConstraintedPositiveEndcapLayer = cms.vuint32(1,1,1),
        RPCRecHitsLabel = cms.InputTag("rpcRecHits"),
        BxRange = cms.uint32(10),
        ClusterSet = cms.vint32(),
        MaxDeltaPhi = cms.double(3.14159265359/3),
        ZError = cms.double(130.0),
        AlgorithmType = cms.uint32(0),
        MagnecticFieldThreshold = cms.double(0.5),
        sampleCount = cms.uint32(20),
        ShareRecHitsNumberThreshold = cms.uint32(1),
        isCheckCandidateOverlap = cms.bool(False),
        isCheckGoodOverlap = cms.bool(True),
        BendingPhiLowerTH = cms.vdouble(0.062),
        BendingPhiUpperTH = cms.vdouble(0.18),
        MeanPt_Parameter0 = cms.vdouble(344.529),
        MeanPt_Parameter1 = cms.vdouble(-108.671),
        MeanPt_Parameter2 = cms.vdouble(11.1376),
        SigmaPt_Parameter0 = cms.vdouble(213.289),
        SigmaPt_Parameter1 = cms.vdouble(-77.8504),
        SigmaPt_Parameter2 = cms.vdouble(7.10351),
        SimHitTag = cms.vInputTag("g4SimHits", "MuonRPCHits"),
        RPCDigiSimLinkTag = cms.InputTag("simMuonRPCDigis", "RPCDigiSimLink"),
        SeedPurityTH = cms.double(0.5),
        useSimData = cms.bool(True)
)

process.content = cms.PSet( 
        outputCommands = cms.untracked.vstring(
            'keep *_*_*_RPCSeed',
            'keep *_rpcRecHits_*_*',
            'keep *_source_*_*',
            'keep SimTracks_g4SimHits_*_*',
            'keep *_g4SimHits_MuonRPCHits_*',
            'keep *_simMuonRPCDigis_RPCDigiSimLink_*'
        )
)
process.out = cms.OutputModule("PoolOutputModule",
        process.content, 
        fileName = cms.untracked.string('muonseed.root')
)

  
process.p = cms.Path(process.myRPCSeed)

process.e = cms.EndPath(process.out)
