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
        AlgorithmType = cms.uint32(1),
        isVertexConstraint = cms.bool(True),
        isContinuousFilter = cms.bool(False),
        MagnecticFieldThreshold = cms.double(0.5),
        sampleCount = cms.uint32(20),
        ShareRecHitsNumberThreshold = cms.uint32(1),
        isCheckCandidateOverlap = cms.bool(False),
        isCheckGoodOverlap = cms.bool(False),
        Cut1234 = cms.vdouble(0.06,4.0,0.06,4.0),
        CutMax = cms.vdouble(0.075,0.09,0.08,0.09),
        BendingPhiLowerTH = cms.vdouble(0.06378,0.09,0.06378,0.09),
        BendingPhiUpperTH = cms.vdouble(0.33,0.48,0.19,0.44),
        BendingPhiFitValueUpperLimit = cms.vdouble(0.22,0.3,0.16,0.3),
        BendingPhiFitSigmaUpperLimit = cms.vdouble(0.19,0.22,0.16,0.22),
        MeanPt_Parameter0 = cms.vdouble(182.529,243.02,289.086,249.456),
        MeanPt_Parameter1 = cms.vdouble(-79.2961,-147.261,-91.6061,-151.007),
        MeanPt_Parameter2 = cms.vdouble(10.9864,25.344,9.68018,25.892),
        SigmaPt_Parameter0 = cms.vdouble(141.976,205.345,180.16,233.097),
        SigmaPt_Parameter1 = cms.vdouble(-49.1516,-87.4596,-53.9042,-98.2261),
        SigmaPt_Parameter2 = cms.vdouble(4.49208,9.6771,4.24705,10.7054),
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
