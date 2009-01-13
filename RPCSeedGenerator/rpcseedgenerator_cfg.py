import FWCore.ParameterSet.Config as cms

process = cms.Process("RPCSeed")

#process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring('rfio:/castor/cern.ch/user/h/hyteng/CMSSW_2_1_10/muongun.root')
)

# MuonRecoGeometryRecord could be used when include following es_module
process.load("Geometry.MuonCommonData.muonIdealGeometryXML_cfi");
#process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi");
process.load("Geometry.CSCGeometry.cscGeometry_cfi");
process.load("Geometry.DTGeometry.dtGeometry_cfi");
process.load("Geometry.RPCGeometry.rpcGeometry_cfi");
process.load("CalibMuon.Configuration.Muon_FakeAlignment_cff");

process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi");
process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi");

# MagneticField.Engine.volumeBasedMagneticField_cfi is obsolete
#process.load("MagneticField.Engine.volumeBasedMagneticField_cfi");
process.load("Configuration.StandardSequences.MagneticField_cff");

process.myRPCSeed = cms.EDProducer('RPCSeedGenerator', 
        RangeofLayersinBarrel = cms.vint32(4),
        RangeofLayersinEndcap = cms.vint32(3),
        isSpecialLayers = cms.bool(False), 
        isMixBarrelwithEndcap = cms.bool(False),
        LayersinBarrel = cms.vint32(1,0,1,0,1,1),
        LayersinEndcap = cms.vint32(1,1,1,1,1,1), 
        RPCRecHitsLabel = cms.InputTag("rpcRecHits"),
        BxRange = cms.uint32(0),
        MaxDeltaPhi = cms.double(3.14159265359/6),
        MaxRSD = cms.double(0.1),
        ShareRecHitsNumberThreshold = cms.uint32(1),
        isCheckcandidateOverlap = cms.bool(True),
        isCheckgoodOverlap = cms.bool(True)
)

process.content = cms.PSet( 
        outputCommands = cms.untracked.vstring('keep *_*_*_RPCSeed',
            'keep *_rpcRecHits_*_*',
            'keep *_source_*_*',
            'keep SimTracks_g4SimHits_*_*',
            'keep *_g4SimHits_MuonRPCHits_*',
            'keep *_simMuonRPCDigis_RPCDigiSimLink_*')
        )
process.out = cms.OutputModule("PoolOutputModule",
        process.content, 
        fileName = cms.untracked.string('muonseed.root')
)

  
process.p = cms.Path(process.myRPCSeed)

process.e = cms.EndPath(process.out)
