import FWCore.ParameterSet.Config as cms

process = cms.Process("myRPC")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'DESIGN41X_V0::All'

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#'rfio:/castor/cern.ch/user/h/hyteng/mc/muongun/muongun415x2_DESIGN41XV0_3.8T_Pt1.0-100.0Gev_Eta-1.0-1.0_BetaSmear/muongun415x2_DESIGN41XV0_3.8T_Pt1.0-100.0Gev_Eta-1.0-1.0_BetaSmear_key.1998.root',
#'rfio:/castor/cern.ch/user/h/hyteng/mc/muongun/muongun415x2_DESIGN41XV0_3.8T_Pt1.0-100.0Gev_Eta-1.0-1.0_BetaSmear/muongun415x2_DESIGN41XV0_3.8T_Pt1.0-100.0Gev_Eta-1.0-1.0_BetaSmear_key.1999.root'
#'file:/tmp/hyteng/muongun415x2_DESIGN41XV0_3.8T_Pt1.0-100.0Gev_Eta-1.0-1.0_BetaSmear_key.sample.root'
        'rfio:/castor/cern.ch/user/h/hyteng/mc/muongun/muongun415x2_DESIGN41XV0_3.8T_Pt1.0-100.0Gev_Eta-1.0-1.0_BetaSmear_key.sample.root'
    )
)

process.myRPCSeed = cms.EDProducer('RPCSeedProducer', 
        BarrelLayerRange = cms.vuint32(4),
        EndcapLayerRange = cms.vuint32(0),
        isCosmic = cms.bool(False),
        isSpecialLayers = cms.bool(False),
        isMixBarrelwithEndcap = cms.bool(False),
        LayersinBarrel = cms.vuint32(1,1,1,1,0,1),
        LayersinEndcap = cms.vuint32(1,1,1,1,1,1),
        ConstraintedBarrelLayer = cms.vuint32(0,0,0,0,0,0),
        ConstraintedNegativeEndcapLayer = cms.vuint32(1,1,1),
        ConstraintedPositiveEndcapLayer = cms.vuint32(1,1,1),
        RPCRecHitsLabel = cms.InputTag("rpcRecHits"),
        BxRange = cms.uint32(10),
        ClusterSet = cms.vint32(),
        MaxDeltaPhi = cms.double(3.14159265359/6),
        ZError = cms.double(130.0),
        AlgorithmType = cms.uint32(0),
        isVertexConstraint = cms.bool(True),
        isContinuousFilter = cms.bool(True),
        MagnecticFieldThreshold = cms.double(0.5),
        SampleCount = cms.uint32(20),
        ShareRecHitsNumberThreshold = cms.uint32(2),
        isCheckCandidateOverlap = cms.bool(False),
        isCheckGoodOverlap = cms.bool(False),
        Cut0 = cms.vdouble(0.,0.,0.,0.008,0.,0.008,0.,0.001,0.,0.001,0.,0.001,0.,0.,0.,0.001,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.),
        Cut1 = cms.vdouble(0.,0.,0.,0.006,0.,0.001,0.,0.006,0.,0.001,0.,0.001,0.,0.,0.,0.006,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.),
        Cut2 = cms.vdouble(0.,0.0062,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.001,0.,0.,0.,0.001,0.,0.001,0.,0.008,0.,0.002,0.,0.001,0.,0.002,0.,0.008),
        CutMax = cms.vdouble(0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.),
        ProbingPhiLowerTH = cms.vdouble(0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.),
        ProbingPhiUpperTH = cms.vdouble(0.,0.0492,0.,0.046,0.,0.044,0.,0.046,0.,0.046,0.,0.038,0.,0.048,0.,0.044,0.,0.048,0.,0.041,0.,0.044,0.,0.044,0.,0.039,0.,0.044,0.,0.036),
        ExhaustivePhiTH = cms.vdouble(0.,0.02,0.,0.0075,0.,0.0042,0.,0.02,0.,0.0042,0.,0.0042,0.,0.0075,0.,0.0042,0.,0.0042,0.,0.0025,0.,0.025,0.,0.025,0.,0.022,0.,0.022,0.,0.025),
        BendingPhiLowerTH = cms.vdouble(0.,0.1,0.,0.1,0.,0.1,0.,0.1,0.,0.1,0.,0.08,0.,0.08,0.,0.08,0.,0.1,0.,0.08,0.,0.1,0.,0.1,0.,0.08,0.,0.08,0.,0.08),
        BendingPhiUpperTH = cms.vdouble(0.,0.5,0.,0.5,0.,0.5,0.,0.5,0.,0.5,0.,0.4,0.,0.5,0.,0.45,0.,0.45,0.,0.45,0.,0.45,0.,0.45,0.,0.4,0.,0.4,0.,0.35),
        BendingPhiFitValueUpperLimit = cms.vdouble(0.,0.4,0.,0.3,0.,0.24,0.,0.3,0.,0.2,0.,0.19,0.,0.24,0.,0.24,0.,0.22,0.,0.22,0.,0.24,0.,0.21,0.,0.22,0.,0.22,0.,0.2),
        BendingPhiFitSigmaUpperLimit = cms.vdouble(0.,0.2,0.,0.16,0.,0.14,0.,0.2,0.,0.15,0.,0.19,0.,0.18,0.,0.16,0.,0.14,0.,0.2,0.,0.2,0.,0.2,0.,0.2,0.,0.2,0.,0.18),
        MeanPt_Parameter0 = cms.vdouble(0.,103.4,0.,192.5,0.,294.7,0.,199,0.,551,0.,501.7,0.,381.4,0.,415,0.,391.3,0.,437.8,0.,258.8,0.,326.6,0.,390.4,0.,404.2,0.,499.2),
        MeanPt_Parameter1 = cms.vdouble(0.,-79.26,0.,-115.1,0.,-148.2,0.,-117.6,0.,-221.1,0.,-208.2,0.,-184.6,0.,-193.9,0.,-179.9,0.,-195.8,0.,-132.3,0.,-152.3,0.,-172,0.,-180.1,0.,-186),
        MeanPt_Parameter2 = cms.vdouble(0.,16.97,0.,20.14,0.,22.58,0.,20.2,0.,27.73,0.,25.84,0.,26.39,0.,26.96,0.,25.15,0.,26.25,0.,19.77,0.,21.25,0.,21.87,0.,23.69,0.,19.73),
        SigmaPt_Parameter0 = cms.vdouble(0.,77.29,0.,78.8,0.,261.1,0.,49.69,0.,91.77,0.,56.25,0.,115,0.,155.7,0.,270.9,0.,62.58,0.,60.42,0.,63.97,0.,50.6,0.,86.25,0.,51.39),
        SigmaPt_Parameter1 = cms.vdouble(0.,-28.41,0.,-27.75,0.,-72.15,0.,-19.6,0.,-31.34,0.,-19.83,0.,-39.08,0.,-48.14,0.,-74.42,0.,-22.98,0.,-23.27,0.,-24.39,0.,-19.13,0.,-31.49,0.,-20.54),
        SigmaPt_Parameter2 = cms.vdouble(0.,2.96,0.,2.835,0.,5.521,0.,2.251,0.,3.034,0.,2.056,0.,3.694,0.,4.123,0.,5.665,0.,2.408,0.,2.571,0.,2.647,0.,2.135,0.,3.197,0.,2.36),
        SimHitTag = cms.InputTag("g4SimHits", "MuonRPCHits"),
        RPCDigiSimLinkTag = cms.InputTag("simMuonRPCDigis", "RPCDigiSimLink"),
        SeedPurityTH = cms.double(0.5),
        useSimData = cms.bool(False)

)


process.content = cms.PSet( 
    outputCommands = cms.untracked.vstring(
        'keep *_ancientMuonSeed_*_*',
        'keep *_MuonSeed_*_*',
        'keep *_SETMuonSeed_*_*',
        'keep *_myRPCSeed_*_*',
        'keep *_mergedAncientandRPCSeeds_*_*',
        'keep *_mergedMuonandRPCSeeds_*_*',
        'keep *_mergedSETMuonandRPCSeeds_*_*',
        'keep *_standAloneMuons_*_*',
        'keep *_standAloneMuonsOld_*_*',
        'keep *_standAloneSETMuons_*_*',
        'keep *_standAloneMuonsRPCSeed_*_*',
        'keep *_standAloneMuonsOldRPCSeed_*_*',
        'keep *_standAloneMuonsRPCSeedOnlyRPC_*_*',
        'keep *_standAloneMuonsOldRPCSeedOnlyRPC_*_*',
        'keep *_standAloneMuonsMuonSeed_*_*',
        'keep *_standAloneMuonsOldMuonSeed_*_*',
        'keep *_standAloneMuonsMergedAncientandRPCSeed_*_*',
        'keep *_standAloneMuonsOldMergedAncientandRPCSeed_*_*',
        'keep *_standAloneSETMuonsMergedSETMuonandRPCSeed_*_*',
        'keep *_standAloneMuonsMergedMuonandRPCSeed_*_*',
        'keep *_standAloneMuonsOldMergedMuonandRPCSeed_*_*',
        'keep *_mergedtruth_MergedTrackTruth_*',
        'keep *_simMuonCSCDigis_MuonCSCStripDigiSimLinks_*',
        'keep *_simMuonCSCDigis_MuonCSCWireDigiSimLinks_*',
        'keep *_g4SimHits_MuonCSCHits_*',
        'keep *_mix_g4SimHitsMuonCSCHits_*',
        'keep *_simMuonDTDigis_*_*',
        'keep *_g4SimHits_MuonDTHits_*',
        'keep *_mix_g4SimHitsMuonDTHits_*',
        'keep *_dt1DRecHits_*_*',
        'keep *_simMuonRPCDigis_RPCDigiSimLink_*',
        'keep *_g4SimHits_MuonRPCHits_*',
        'keep *_mix_g4SimHitsMuonRPCHits_*',
        'keep *_rpcRecHits_*_*',
        'keep SimTracks_g4SimHits_*_*'
    )
)

process.out = cms.OutputModule("PoolOutputModule",
    process.content, 
    fileName = cms.untracked.string('recSTA415_Mix_VertexnoFilternoOverlap_Pt1.0-100.0Gev_key.root')
)


# for the old mode
process.standAloneMuonsOld = process.standAloneMuons.clone()
process.standAloneMuonsOld.MuonTrajectoryBuilder = cms.string("StandAloneMuonTrajectoryBuilder")

# Stand alone made with RPC seeds as starting states, but it can collect DT/CSC/RPC hits during pattern recognition 
process.standAloneMuonsRPCSeed = process.standAloneMuons.clone()
process.standAloneMuonsRPCSeed.InputObjects = cms.InputTag('myRPCSeed','GoodSeeds')

# Stand alone made with RPC seeds as starting states, but it can collect DT/CSC/RPC hits during pattern recognition 
process.standAloneMuonsOldRPCSeed = process.standAloneMuons.clone()
process.standAloneMuonsOldRPCSeed.InputObjects = cms.InputTag('myRPCSeed','GoodSeeds')
process.standAloneMuonsOldRPCSeed.MuonTrajectoryBuilder = cms.string("StandAloneMuonTrajectoryBuilder")

# Stand alone made with Muon seeds (new DT/CSC seed) as starting states, but it can collect DT/CSC/RPC hits during pattern recognition
process.standAloneMuonsMuonSeed = process.standAloneMuons.clone()
process.standAloneMuonsMuonSeed.InputObjects = cms.InputTag('MuonSeed')

# Stand alone made with RPC&Ancient merged seeds as starting states, but it can collect DT/CSC/RPC hits during pattern recognition
process.standAloneMuonsMergedAncientandRPCSeed = process.standAloneMuons.clone()
process.standAloneMuonsMergedAncientandRPCSeed.InputObjects = cms.InputTag('mergedAncientandRPCSeeds')

# Stand alone made with RPC&Ancient merged seeds as starting states, and use old mode,also it can collect DT/CSC/RPC hits during pattern recognition
process.standAloneMuonsOldMergedAncientandRPCSeed = process.standAloneMuons.clone()
process.standAloneMuonsOldMergedAncientandRPCSeed.InputObjects = cms.InputTag('mergedAncientandRPCSeeds')
process.standAloneMuonsOldMergedAncientandRPCSeed.MuonTrajectoryBuilder = cms.string("StandAloneMuonTrajectoryBuilder")

# Stand alone made with RPC&Muon(new DT/CSC seed) merged seeds as starting states, but it can collect DT/CSC/RPC hits during pattern recognition
process.standAloneMuonsMergedMuonandRPCSeed = process.standAloneMuons.clone()
process.standAloneMuonsMergedMuonandRPCSeed.InputObjects = cms.InputTag('mergedMuonandRPCSeeds')

# Stand alone made with RPC&Muon(new DT/CSC seed) merged seeds as starting states, with old mode, also it can collect DT/CSC/RPC hits during pattern recognition
process.standAloneMuonsOldMergedMuonandRPCSeed = process.standAloneMuons.clone()
process.standAloneMuonsOldMergedMuonandRPCSeed.InputObjects = cms.InputTag('mergedMuonandRPCSeeds')
process.standAloneMuonsOldMergedMuonandRPCSeed.MuonTrajectoryBuilder = cms.string("StandAloneMuonTrajectoryBuilder")

# Stand alone made with RPC&Muon(new DT/CSC seed) merged seeds as starting states, but it can collect DT/CSC/RPC hits during pattern recognition
process.standAloneSETMuonsMergedSETMuonandRPCSeed = process.standAloneSETMuons.clone()
process.standAloneSETMuonsMergedSETMuonandRPCSeed.InputObjects = cms.InputTag('mergedSETMuonandRPCSeeds')

# Stand alone made with RPC hits only
process.standAloneMuonsRPCSeedOnlyRPC = process.standAloneMuonsRPCSeed.clone()
process.standAloneMuonsRPCSeedOnlyRPC.STATrajBuilderParameters.FilterParameters.EnableCSCMeasuremen = cms.bool(False)
process.standAloneMuonsRPCSeedOnlyRPC.STATrajBuilderParameters.BWFilterParameters.EnableCSCMeasuremen = cms.bool(False)
process.standAloneMuonsRPCSeedOnlyRPC.STATrajBuilderParameters.FilterParameters.EnableDTMeasuremen = cms.bool(False)
process.standAloneMuonsRPCSeedOnlyRPC.STATrajBuilderParameters.BWFilterParameters.EnableDTMeasuremen = cms.bool(False)
process.standAloneMuonsRPCSeedOnlyRPC.STATrajBuilderParameters.FilterParameters.MuonTrajectoryUpdatorParameters.ExcludeRPCFromFit = cms.bool(False)
process.standAloneMuonsRPCSeedOnlyRPC.STATrajBuilderParameters.BWFilterParameters.MuonTrajectoryUpdatorParameters.ExcludeRPCFromFit = cms.bool(False)
# Stand alone made with RPC hits only
process.standAloneMuonsOldRPCSeedOnlyRPC = process.standAloneMuonsOldRPCSeed.clone()
process.standAloneMuonsOldRPCSeedOnlyRPC.STATrajBuilderParameters.FilterParameters.EnableCSCMeasuremen = cms.bool(False)
process.standAloneMuonsOldRPCSeedOnlyRPC.STATrajBuilderParameters.BWFilterParameters.EnableCSCMeasuremen = cms.bool(False)
process.standAloneMuonsOldRPCSeedOnlyRPC.STATrajBuilderParameters.FilterParameters.EnableDTMeasuremen = cms.bool(False)
process.standAloneMuonsOldRPCSeedOnlyRPC.STATrajBuilderParameters.BWFilterParameters.EnableDTMeasuremen = cms.bool(False)
process.standAloneMuonsOldRPCSeedOnlyRPC.STATrajBuilderParameters.FilterParameters.MuonTrajectoryUpdatorParameters.ExcludeRPCFromFit = cms.bool(False)
process.standAloneMuonsOldRPCSeedOnlyRPC.STATrajBuilderParameters.BWFilterParameters.MuonTrajectoryUpdatorParameters.ExcludeRPCFromFit = cms.bool(False)

# merged seeds with ancient and RPC
process.mergedAncientandRPCSeeds = process.mergedStandAloneMuonSeeds.clone()
process.mergedAncientandRPCSeeds.SeedCollections = cms.VInputTag(cms.InputTag("ancientMuonSeed"),  cms.InputTag('myRPCSeed','GoodSeeds'))

# merged seeds with Muon and RPC
process.mergedMuonandRPCSeeds = process.mergedStandAloneMuonSeeds.clone()
process.mergedMuonandRPCSeeds.SeedCollections = cms.VInputTag(cms.InputTag("MuonSeed"),  cms.InputTag('myRPCSeed','GoodSeeds'))
    
# merged seeds with SETMuon and RPC
process.mergedSETMuonandRPCSeeds = process.mergedStandAloneMuonSeeds.clone()
process.mergedSETMuonandRPCSeeds.SeedCollections = cms.VInputTag(cms.InputTag("SETMuonSeed"),  cms.InputTag('myRPCSeed','GoodSeeds'))
    
# process.MuonSeed is new seed and will use sta DT/CSC seeds and all DT/CSC/RPC recHits, while StandAloneMuonSeeds is old ancientSeed, SETMuonSeed is with RC/DT/CSC seeds and all RPC/DT/CSC recHits
#process.p = cms.Path(process.myRPCSeed * process.standAloneMuonSeeds * process.MuonSeed* process.SETMuonSeed * process.mergedAncientandRPCSeeds * process.mergedMuonandRPCSeeds * process.mergedSETMuonandRPCSeeds * process.standAloneMuons * process.standAloneSETMuons * process.standAloneMuonsRPCSeed * process.standAloneMuonsMuonSeed * process.standAloneMuonsMergedAncientandRPCSeed * process.standAloneMuonsMergedMuonandRPCSeed * process.standAloneMuonsMergedSETMuonandRPCSeed)

#process.p = cms.Path(process.myRPCSeed * process.standAloneMuonSeeds * process.MuonSeed * process.SETMuonSeed * process.mergedAncientandRPCSeeds * process.mergedMuonandRPCSeeds * process.mergedSETMuonandRPCSeeds * process.standAloneMuons * process.standAloneMuonsMuonSeed * process.standAloneSETMuons * process.standAloneMuonsRPCSeed * process.standAloneMuonsRPCSeedOnlyRPC * process.standAloneMuonsMergedAncientandRPCSeed * process.standAloneMuonsMergedMuonandRPCSeed * process.standAloneSETMuonsMergedSETMuonandRPCSeed)

# MuonSeed branch (process.MuonSeed) is kickout since it has some bugs and will lead to crash in running
process.p = cms.Path(process.myRPCSeed * process.standAloneMuonSeeds * process.SETMuonSeed * process.mergedAncientandRPCSeeds * process.mergedSETMuonandRPCSeeds * process.standAloneMuons * process.standAloneMuonsOld * process.standAloneSETMuons * process.standAloneMuonsRPCSeed * process.standAloneMuonsOldRPCSeed * process.standAloneMuonsRPCSeedOnlyRPC * process.standAloneMuonsOldRPCSeedOnlyRPC * process.standAloneMuonsMergedAncientandRPCSeed * process.standAloneMuonsOldMergedAncientandRPCSeed * process.standAloneSETMuonsMergedSETMuonandRPCSeed)

process.e = cms.EndPath(process.out)
