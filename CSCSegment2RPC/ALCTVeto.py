import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration/StandardSequences/FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'CRAFT_ALL_V11::All'
process.load("Configuration/StandardSequences/RawToDigi_Data_cff")
process.load("Configuration.StandardSequences.ReconstructionCosmics_cff")

process.source = cms.Source ("PoolSource",
fileNames = cms.untracked.vstring (
# CRAFT data, run 69912, RAW, just a few ramdom files...
'/store/data/Commissioning08/Cosmics/RAW/v1/000/069/912/F2015591-1BAD-DD11-B412-001617E30D4A.root',
'/store/data/Commissioning08/Cosmics/RAW/v1/000/069/912/F29C19AC-32AD-DD11-8430-000423D9890C.root',
'/store/data/Commissioning08/Cosmics/RAW/v1/000/069/912/F42B3FD8-13AD-DD11-BE00-000423D6BA18.root',
'/store/data/Commissioning08/Cosmics/RAW/v1/000/069/912/F42CA98E-29AD-DD11-8782-0016177CA7A0.root',
'/store/data/Commissioning08/Cosmics/RAW/v1/000/069/912/F4370715-21AD-DD11-A59E-000423D98868.root',
'/store/data/Commissioning08/Cosmics/RAW/v1/000/069/912/F60B65BF-1FAD-DD11-9262-001617DBD316.root',
'/store/data/Commissioning08/Cosmics/RAW/v1/000/069/912/F61CBAB7-0CAD-DD11-9943-000423D98800.root',
'/store/data/Commissioning08/Cosmics/RAW/v1/000/069/912/F638405B-17AD-DD11-AD07-001617C3B76A.root',
'/store/data/Commissioning08/Cosmics/RAW/v1/000/069/912/F820B3F4-2AAD-DD11-8C18-000423D60FF6.root',
'/store/data/Commissioning08/Cosmics/RAW/v1/000/069/912/F86ABAE7-09AD-DD11-BE3A-001617C3B73A.root',
'/store/data/Commissioning08/Cosmics/RAW/v1/000/069/912/F8B79B06-26AD-DD11-B0DB-001D09F241B9.root',
'/store/data/Commissioning08/Cosmics/RAW/v1/000/069/912/F8BB3211-2BAD-DD11-87A5-001617DBD316.root',
'/store/data/Commissioning08/Cosmics/RAW/v1/000/069/912/F8C08E78-0DAD-DD11-8583-000423D94AA8.root',
'/store/data/Commissioning08/Cosmics/RAW/v1/000/069/912/F8D57320-34AD-DD11-BCCC-001D09F241B9.root',
'/store/data/Commissioning08/Cosmics/RAW/v1/000/069/912/FA1DB50B-13AD-DD11-A582-000423D99A8E.root',
'/store/data/Commissioning08/Cosmics/RAW/v1/000/069/912/FA4EB7BE-13AD-DD11-8F6B-000423D6B2D8.root',
'/store/data/Commissioning08/Cosmics/RAW/v1/000/069/912/FC387206-0CAD-DD11-8E87-000423D98DD4.root',
'/store/data/Commissioning08/Cosmics/RAW/v1/000/069/912/FE3FC648-17AD-DD11-AC3C-000423D94908.root',
'/store/data/Commissioning08/Cosmics/RAW/v1/000/069/912/FE85826D-2CAD-DD11-8C67-000423D98E6C.root',
'/store/data/Commissioning08/Cosmics/RAW/v1/000/069/912/FEA34B38-2FAD-DD11-9CFF-001617E30F4C.root'

),
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
    )


 
#------------------------------------------
# parameters for the ALCTVeto module
#------------------------------------------
process.alctVeto = cms.EDFilter(
    "ALCTVeto",
    alctDigiTag = cms.InputTag("muonCSCDigis","MuonCSCALCTDigi")
)

#### the path
# muonCSCDigis is needed to create the ALCT Digis from the RAW data.
process.study = cms.Path( process.muonCSCDigis * process.alctVeto )


#### output 
process.outputSkim = cms.OutputModule(
        "PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *','drop *_MEtoEDMConverter_*_*'),
    fileName = cms.untracked.string("/tmp/schmittm/test.root"),
    dataset = cms.untracked.PSet(
      dataTier = cms.untracked.string('RAW-RECO'),
      filterName = cms.untracked.string('ALCTVeto')
    ),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('study'))
)

### process.outpath = cms.EndPath(process.outputSkim)
#
