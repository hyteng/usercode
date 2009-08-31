import FWCore.ParameterSet.Config as cms

process = cms.Process("RPCSegmentEff")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Geometry.MuonCommonData.muonIdealGeometryXML_cfi")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
process.load("Geometry.CSCGeometry.cscGeometry_cfi")
process.load("Geometry.DTGeometry.dtGeometry_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
process.load("DQMServices.Components.MEtoEDMConverter_cfi")
process.load("DQMServices.Core.DQM_cfg")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.connect = "frontier://PromptProd/CMS_COND_21X_GLOBALTAG"
#process.GlobalTag.globaltag = "CRUZET4_V5P::All"
#process.prefer("GlobalTag")
process.GlobalTag.globaltag = "MC_31X_V2::All"

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#'file:muongun.root'
#'rfio:/castor/cern.ch/cms/store/data/Commissioning08/Cosmics/RECO/v1/000/069/912/00AAB58A-1AAD-DD11-8E38-000423D98750.root'
#'/store/data/BeamCommissioning08/BeamHalo/RECO/v1/000/063/440/0058A079-D786-DD11-99CE-000423D952C0.root'
#'//store/mc/Summer08/CosmicMCBOn10GeV/RECO/COSMMC_21X_v4/0005/0043C01D-92AE-DD11-BB89-001EC9AA8E47.root'
'rfio:/castor/cern.ch/user/h/hyteng/muongun_endcap/MuonGun_Endcap_E100Gev_38T_MC31XV2.085.root'
        )
)

process.rpceff = cms.EDAnalyzer('CSCSegment2RPC',
        MaxD = cms.double(80.0),
        DuplicationCorrection = cms.int32(1),
        rangestrips = cms.double(1.),
        EffSaveRootFile = cms.untracked.bool(True),
        EffRootFileName = cms.untracked.string('RPCEfficiency.root'),
        CSCSegmentTag = cms.InputTag("cscSegments"),
        CSCStripDigiTag = cms.InputTag("muonCSCDigis", "MuonCSCStripDigi"),
        CSCWireDigiTag = cms.InputTag("muonCSCDigis", "MuonCSCWireDigi"),
        checkType = cms.uint32(1),
        sampleType = cms.uint32(1),
        TrackType = cms.int32(13),
        isCSCSegmentFilter = cms.bool(True),
        deltaRTH = cms.double(0.2),
        minSegmentEta = cms.double(0.0),
        maxSegmentEta = cms.double(3.0),
        residualDistanceTH = cms.double(5.0),
        RPCsimHitTag = cms.InputTag("g4SimHits", "MuonRPCHits"),
        RPCDigiTag = cms.InputTag("muonRPCDigis"),
        RPCrecHitTag = cms.InputTag("rpcRecHits"),
        SimTrackTag = cms.InputTag("g4SimHits"),
        
)

process.FEVT = cms.OutputModule("PoolOutputModule",
        outputCommands = cms.untracked.vstring('keep *_MEtoEDMConverter_*_*'),
        fileName = cms.untracked.string('output.root')
        )

process.p = cms.Path(process.rpceff*process.MEtoEDMConverter)

process.outpath = cms.EndPath(process.FEVT)

process.DQM.collectorHost = ''
process.DQM.collectorPort = 9090
    
