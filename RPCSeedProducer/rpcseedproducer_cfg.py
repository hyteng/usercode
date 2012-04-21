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
        BarrelLayerRange = cms.vuint32(4,5),
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
        sampleCount = cms.uint32(20),
        ShareRecHitsNumberThreshold = cms.uint32(2),
        isCheckCandidateOverlap = cms.bool(False),
        isCheckGoodOverlap = cms.bool(False),
        Cut1234 = cms.vdouble(0.06,4.0,0.06,4.0,0.065,0.05,-1,0.05,-1,0.05,-1,0.06),
        #CutMax = cms.vdouble(0.075,0.08,0.08,0.09,0.065,0.08,-1,0.09,-1,0.09,-1,0.08),
        CutMax = cms.vdouble(0.05,0.0,0.05,0.0,0.05,4.0,0.05,4.0,0.05,4.0,0.05,4.0),
        ProbingPhiLowerTH = cms.vdouble(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),
        ProbingPhiUpperTH = cms.vdouble(0.0,0.08,0.0,0.09,0.0,0.08,0.0,0.09,0.0,0.09,0.0,0.1),
        ExhaustivePhiTH = cms.vdouble(0.0,0.08,0.0,0.8,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),
        BendingPhiLowerTH = cms.vdouble(0.06378,0.08,0.06378,0.09,0.065,0.08,-1,0.09,-1,0.09,-1,0.103),
        BendingPhiUpperTH = cms.vdouble(0.33,0.48,0.19,0.44,0.3,0.48,-1,0.48,-1,0.48,-1,0.36),
        BendingPhiFitValueUpperLimit = cms.vdouble(0.22,0.3,0.16,0.3,0.2,0.28,-1,0.32,-1,0.3,-1,0.16),
        BendingPhiFitSigmaUpperLimit = cms.vdouble(0.19,0.22,0.16,0.22,0.17,0.21,-1,0.28,-1,0.28,-1,0.16),
        MeanPt_Parameter0 = cms.vdouble(182.529,304.451,289.086,249.456,336.14,280.45,-1,203.385,-1,203.374,-1,648.297),
        MeanPt_Parameter1 = cms.vdouble(-79.2961,-175.507,-91.6061,-151.007,-123.554,-160.707,-1,-126.38,-1,-126.07,-1,-271.243),
        MeanPt_Parameter2 = cms.vdouble(10.9864,28.4037,9.68018,25.892,13.8076,26.2136,-1,22.6179,-1,22.5726,-1,30.0075),
        SigmaPt_Parameter0 = cms.vdouble(141.976,230.736,180.16,233.097,230.982,266.18,-1,108.439,-1,105.957,-1,714.203),
        SigmaPt_Parameter1 = cms.vdouble(-49.1516,-97.1726,-53.9042,-98.2261,-73.1826,-104.452,-1,-51.8928,-1,-51.1105,-1,-224.504),
        SigmaPt_Parameter2 = cms.vdouble(4.49208,10.5827,4.24705,10.7054,6.04733,10.642,-1,6.46826,-1,6.41354,-1,18.3735),
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
