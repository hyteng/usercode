// -*- C++ -*-
//
// Package:    CSCSegment2RPC
// Class:      CSCSegment2RPC
// 
/**\class CSCSegment2RPC CSCSegment2RPC.cc MyMTCCAnalyzer/CSCSegment2RPC/src/CSCSegment2RPC.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
 */
//
// Original Author:  Haiyun Teng
//         Created:  Wed Feb 11 10:40:41 CET 2009
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// user special include files
#include <FWCore/ServiceRegistry/interface/Service.h>
#include <FWCore/ParameterSet/interface/InputTag.h>
#include <DataFormats/RPCRecHit/interface/RPCRecHitCollection.h>
#include <DataFormats/RPCDigi/interface/RPCDigiCollection.h>
#include <DataFormats/MuonDetId/interface/RPCDetId.h>
#include <DataFormats/CSCRecHit/interface/CSCSegmentCollection.h>
#include <DataFormats/CSCDigi/interface/CSCWireDigiCollection.h>
#include <DataFormats/CSCDigi/interface/CSCStripDigiCollection.h>
#include <DataFormats/MuonDetId/interface/CSCDetId.h>
#include <Geometry/RPCGeometry/interface/RPCGeomServ.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include <FWCore/Framework/interface/ESHandle.h>
#include <Geometry/RPCGeometry/interface/RPCGeometry.h>
#include <Geometry/RPCGeometry/interface/RPCRoll.h>
#include <Geometry/CSCGeometry/interface/CSCGeometry.h>
#include <Geometry/CSCGeometry/interface/CSCChamber.h>
#include <DataFormats/GeometrySurface/interface/Bounds.h>
#include <Geometry/CSCGeometry/interface/CSCLayerGeometry.h>
#include <Geometry/CSCGeometry/interface/CSCStripTopology.h>
#include <DQMOffline/Muon/interface/RPCBookFolderStructure.h>
#include <SimDataFormats/TrackingHit/interface/PSimHitContainer.h>
#include <DataFormats/RPCDigi/interface/RPCDigiCollection.h>
#include <SimDataFormats/Track/interface/SimTrackContainer.h>
#include <Geometry/CommonDetUnit/interface/GeomDetUnit.h>

#include <cmath>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TString.h"

#include <TTree.h>
#include <TFile.h>
#include <TMath.h>

#ifndef PI
#define PI 3.1415926535
#endif

using namespace edm;
using namespace std;
//
// class decleration
//
/*
   class CSCStationIndex {
   public:
   CSCStationIndex():_region(0),_station(0),_ring(0),_chamber(0){}
   CSCStationIndex(int region, int station, int ring, int chamber):
   _region(region),
   _station(station),
   _ring(ring),
   _chamber(chamber){}
   ~CSCStationIndex(){}
   int region() const {return _region;}
   int station() const {return _station;}
   int ring() const {return _ring;}
   t
   q
   n
   int chamber() const {return _chamber;}
   bool operator<(const CSCStationIndex& cscind) const
   {
   if(cscind.region()!=this->region())
   return cscind.region()<this->region();
   else if(cscind.station()!=this->station())
   return cscind.station()<this->station();
   else if(cscind.ring()!=this->ring())
   return cscind.ring()<this->ring();
   else if(cscind.chamber()!=this->chamber())
   return cscind.chamber()<this->chamber();
   return false;
   }

   private:
   int _region;
   int _station;
   int _ring;  
   int _chamber;
   };
 */
class CSCSegment2RPC : public edm::EDAnalyzer {

    public:
        explicit CSCSegment2RPC(const edm::ParameterSet& iConfig);
        ~CSCSegment2RPC();

    private:
        virtual void beginRun(const edm::Run& Run, const edm::EventSetup&);
        virtual void beginJob(const edm::EventSetup& iSetup) ;
        virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
        virtual void endRun(const edm::Run& Run, const edm::EventSetup& iSetup);
        virtual void endJob() ;

        std::map<std::string, MonitorElement*> bookDetUnitSeg(RPCDetId& detId,int nstrips, double stripw, double stripl);
        void sampleCSCSegments();
        void sampleRPCsimHits();
        bool filterCSCSegment(const CSCSegment& sampleCSCSegment);
        bool filterSegmentEta(const CSCSegment& sampleCSCSegment);
        double deltaRforSegment(const GlobalPoint& globalPosition1, const GlobalVector& globalVector1, const GlobalPoint& globalPosition2, const GlobalVector& globalVector2);
        void findpeakTime(const CSCSegment& sampleCSCSegment);
        void fillSample();
        void EfficiencybyRPCrecHit();
        void fillEfficiency();

        // ----------member data ---------------------------
        edm::InputTag CSCSegmentsTag_;
        edm::InputTag CSCStripDigisTag_;
        edm::InputTag CSCWireDigisTag_;
        //edm::InputTag CSCStripDigiSimLinkTag_;
        edm::InputTag RPCsimHitTag_;
        edm::InputTag RPCDigiTag_;
        edm::InputTag RPCrecHitTag_;
        //edm::InputTag SimTrackTag_;

        edm::ESHandle<RPCGeometry> pRPCGeom;
        RPCGeometry* rpcGeometry;
        edm::ESHandle<CSCGeometry> pCSCGeom;
        CSCGeometry* cscGeometry;

        edm::Handle<CSCSegmentCollection> pCSCSegments;
        edm::Handle<CSCStripDigiCollection> pCSCStripDigiCollection;
        edm::Handle<CSCWireDigiCollection> pCSCWireDigiCollection;
        edm::Handle<PSimHitContainer> pRPCsimHits;
        edm::Handle<RPCDigiCollection> pRPCDigis;
        edm::Handle<RPCRecHitCollection> pRPCrecHits;
        edm::Handle<SimTrackContainer> pSimTracks;

        bool isCSCSegmentFilter;
        unsigned int sampleType;
        unsigned int checkType;
        double deltaRTH;
        double minSegmentEta;
        double maxSegmentEta;
        double residualDistanceTH;
        int TrackType;

        int dupli;
        double rangestrips;
        double MaxD;
        double RollSize;
        double StripLength;
        double StripWidth;

        double CSCpeakTime2nd;
        double averageCSCpeakTime;
        int CSCWireTimeBin;  
        int RPCBX;
        double deltaBXTime;

        RPCDetId impactRPCRollId;
        LocalVector impactDirection;
        LocalPoint impactPoint;
        double impactStrip;
        double impactResidual;
        double impactClusterSize;

        bool isvalideRPCBX;
        bool anycoincidence;

        bool issampled;
        bool isSamplefilled;
        bool ischecked;
        bool isEfficiencyfilled;

        char meRPCId [128];
        bool EffSaveRootFile;
        std::string EffRootFileName;

        std::map< CSCDetId, std::vector<RPCDetId> > RollswithCSCChamber;
        std::map< RPCDetId, std::map<std::string, MonitorElement*> >  meCollection;

        DQMStore* DBE;
        MonitorElement* Statistics;
        MonitorElement* hGlobalResClu1R3C;
        MonitorElement* hGlobalResClu1R3B;
        MonitorElement* hGlobalResClu1R3A;
        MonitorElement* hGlobalResClu1R2C;
        MonitorElement* hGlobalResClu1R2B; 
        MonitorElement* hGlobalResClu1R2A;

        MonitorElement* hGlobalResClu2R3C;
        MonitorElement* hGlobalResClu2R3B;
        MonitorElement* hGlobalResClu2R3A;
        MonitorElement* hGlobalResClu2R2C;
        MonitorElement* hGlobalResClu2R2B;
        MonitorElement* hGlobalResClu2R2A;

        MonitorElement* hGlobalResClu3R3C;
        MonitorElement* hGlobalResClu3R3B;
        MonitorElement* hGlobalResClu3R3A;
        MonitorElement* hGlobalResClu3R2C;
        MonitorElement* hGlobalResClu3R2B;
        MonitorElement* hGlobalResClu3R2A;

        // Ntuple method
        unsigned int impactRPCId_NTuple;
        int impactCSCBX_NTuple;
        double impactAngle_NTuple;
        int impactStrip_NTuple;
        unsigned int isEfficiency_NTuple;
        double impactResidual_NTuple;
        int impactClusterSize_NTuple;
        int impactRPCBX_NTuple;

        bool SavetheRootFile;
        string theRootFileName;
        TFile *theFile;
        TTree *ExTree;       
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
CSCSegment2RPC::CSCSegment2RPC(const edm::ParameterSet& iConfig) {

    //now do what ever initialization is needed
    CSCSegmentsTag_ = iConfig.getParameter<edm::InputTag>("CSCSegmentTag");
    CSCStripDigisTag_ = iConfig.getParameter<edm::InputTag>("CSCStripDigiTag");
    CSCWireDigisTag_ = iConfig.getParameter<edm::InputTag>("CSCWireDigiTag");
    //CSCStripDigiSimLinkTag_ = iConfig.getParameter<edm::InputTag>("CSCStripDigiSimLinkTag");

    RPCsimHitTag_ = iConfig.getParameter<edm::InputTag>("RPCsimHitTag");
    RPCDigiTag_ = iConfig.getParameter<edm::InputTag>("RPCDigiTag");
    RPCrecHitTag_ = iConfig.getParameter<edm::InputTag>("RPCrecHitTag");

    //SimTrackTag_ = iConfig.getParameter<edm::InputTag>("SimTrackTag");

    isCSCSegmentFilter = iConfig.getParameter<bool>("isCSCSegmentFilter");
    sampleType = iConfig.getParameter<unsigned int>("sampleType");
    checkType = iConfig.getParameter<unsigned int>("checkType");
    TrackType = iConfig.getParameter<int>("TrackType");
    deltaRTH = iConfig.getParameter<double>("deltaRTH");
    minSegmentEta = iConfig.getParameter<double>("minSegmentEta");
    maxSegmentEta = iConfig.getParameter<double>("maxSegmentEta");
    residualDistanceTH = iConfig.getParameter<double>("residualDistanceTH");
    rangestrips = iConfig.getParameter<double>("rangestrips");
    dupli = iConfig.getParameter<int>("DuplicationCorrection");
    MaxD = iConfig.getParameter<double>("MaxD");

    issampled = false;
    isSamplefilled = false;
    ischecked = false;
    isEfficiencyfilled = false;

    EffSaveRootFile  = iConfig.getUntrackedParameter<bool>("EffSaveRootFile", false);
    EffRootFileName  = iConfig.getUntrackedParameter<std::string>("EffRootFileName", "RPCEfficiency.root"); 
    SavetheRootFile = iConfig.getUntrackedParameter<bool>("EffSaveRootFile2", false);
    theRootFileName = iConfig.getUntrackedParameter<std::string>("EffRootFileName2", "RPCEfficiencyTree.root");

    DBE = edm::Service<DQMStore>().operator->();
    std::string Folder = "Muons/MuonSegEff/";
    DBE->setCurrentFolder(Folder);
    Statistics = DBE->book1D("Statistics","All Statistics",33,0.5,33.5);

    std::cout << "booking Global histograms" << std::endl;

    Folder = "Muons/MuonSegEff/Residuals/EndCap";
    DBE->setCurrentFolder(Folder);

    //Endcap  
    hGlobalResClu1R3C = DBE->book1D("GlobalResidualsClu1R3C", "RPC Residuals Ring 3 Roll C Cluster Size 1", 101, -10., 10.);
    hGlobalResClu1R3B = DBE->book1D("GlobalResidualsClu1R3B", "RPC Residuals Ring 3 Roll B Cluster Size 1", 101, -10., 10.);
    hGlobalResClu1R3A = DBE->book1D("GlobalResidualsClu1R3A", "RPC Residuals Ring 3 Roll A Cluster Size 1", 101, -10., 10.);
    hGlobalResClu1R2C = DBE->book1D("GlobalResidualsClu1R2C", "RPC Residuals Ring 2 Roll C Cluster Size 1", 101, -10., 10.);
    hGlobalResClu1R2B = DBE->book1D("GlobalResidualsClu1R2B", "RPC Residuals Ring 2 Roll B Cluster Size 1", 101, -10., 10.);
    hGlobalResClu1R2A = DBE->book1D("GlobalResidualsClu1R2A", "RPC Residuals Ring 2 Roll A Cluster Size 1", 101, -10., 10.);

    hGlobalResClu2R3C = DBE->book1D("GlobalResidualsClu2R3C", "RPC Residuals Ring 3 Roll C Cluster Size 2", 101, -10., 10.);
    hGlobalResClu2R3B = DBE->book1D("GlobalResidualsClu2R3B", "RPC Residuals Ring 3 Roll B Cluster Size 2", 101, -10., 10.);
    hGlobalResClu2R3A = DBE->book1D("GlobalResidualsClu2R3A", "RPC Residuals Ring 3 Roll A Cluster Size 2", 101, -10., 10.);
    hGlobalResClu2R2C = DBE->book1D("GlobalResidualsClu2R2C", "RPC Residuals Ring 2 Roll C Cluster Size 2", 101, -10., 10.);
    hGlobalResClu2R2B = DBE->book1D("GlobalResidualsClu2R2B", "RPC Residuals Ring 2 Roll B Cluster Size 2", 101, -10., 10.);
    hGlobalResClu2R2A = DBE->book1D("GlobalResidualsClu2R2A", "RPC Residuals Ring 2 Roll A Cluster Size 2", 101, -10., 10.);

    hGlobalResClu3R3C = DBE->book1D("GlobalResidualsClu3R3C", "RPC Residuals Ring 3 Roll C Cluster Size 3", 101, -10., 10.);
    hGlobalResClu3R3B = DBE->book1D("GlobalResidualsClu3R3B", "RPC Residuals Ring 3 Roll B Cluster Size 3", 101, -10., 10.);
    hGlobalResClu3R3A = DBE->book1D("GlobalResidualsClu3R3A", "RPC Residuals Ring 3 Roll A Cluster Size 3", 101, -10., 10.);
    hGlobalResClu3R2C = DBE->book1D("GlobalResidualsClu3R2C", "RPC Residuals Ring 2 Roll C Cluster Size 3", 101, -10., 10.);
    hGlobalResClu3R2B = DBE->book1D("GlobalResidualsClu3R2B", "RPC Residuals Ring 2 Roll B Cluster Size 3", 101, -10., 10.);
    hGlobalResClu3R2A = DBE->book1D("GlobalResidualsClu3R2A", "RPC Residuals Ring 2 Roll A Cluster Size 3", 101, -10., 10.);

    theFile = new TFile(theRootFileName.c_str(), "recreate");
    theFile->cd();
    ExTree= new TTree("ExTree", "ExTree");
    ExTree->Branch("impactRPC", &impactRPCId_NTuple, "impactRPCId/i");
    ExTree->Branch("impactCSCBX", &impactCSCBX_NTuple, "impactCSCBX/I");
    ExTree->Branch("impactAngle", &impactAngle_NTuple, "impactAngle/D");
    ExTree->Branch("impactStrip", &impactStrip_NTuple, "impactStrip/I");
    ExTree->Branch("isEfficiency", &isEfficiency_NTuple, "isEfficiency/i");
    ExTree->Branch("impactResidual", &impactResidual_NTuple, "impactResidual/D");
    ExTree->Branch("impactClusterSize", &impactClusterSize_NTuple, "impactClusterSize/I");
    ExTree->Branch("impactRPCBX", &impactRPCBX_NTuple, "impactRPCBX/I");
}


CSCSegment2RPC::~CSCSegment2RPC() {

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//
void CSCSegment2RPC::beginRun(const edm::Run& Run, const edm::EventSetup&) {

}

std::map<std::string, MonitorElement*> CSCSegment2RPC::bookDetUnitSeg(RPCDetId& detId, int nstrips, double stripw, double stripl) {

    std::map<std::string, MonitorElement*> meMap;

    RPCBookFolderStructure*  FolderStr = new RPCBookFolderStructure(); //Anna
    std::string Folder = "Muons/MuonSegEff/" +  FolderStr->folderStructure(detId);

    DBE->setCurrentFolder(Folder);

    RPCGeomServ RPCname(detId);
    std::string nameRoll = RPCname.name();
    char detUnitLabel[128];
    char layerLabel[128];

    sprintf(detUnitLabel ,"%s", nameRoll.c_str());
    sprintf(layerLabel ,"%s", nameRoll.c_str());

    char meId [128];
    char meTitle [128];

    int rawId = detId.rawId();
    cout << "Booking histogram for RPC roll " << rawId << ", name: " << detUnitLabel << endl;

    // old: from 0.5 to nstrips+0.5 new: from 0 to nstrips
    sprintf(meId, "ExpectedOccupancyFromCSC_%s", detUnitLabel);
    sprintf(meTitle, "ExpectedOccupancyFromCSC_for_%s", layerLabel);
    cout << "Booking for " << meId << endl;
    meMap[meId] = DBE->book1D(meId, meTitle, nstrips, 0, nstrips);

    sprintf(meId, "RPCDataOccupancyFromCSC_%s", detUnitLabel);
    sprintf(meTitle, "RPCDataOccupancyFromCSC_for_%s", layerLabel);
    cout << "Booking for " << meId << endl;
    meMap[meId] = DBE->book1D(meId, meTitle, nstrips, 0, nstrips);

    sprintf(meId, "RPCBXDistribution_%s", detUnitLabel);
    sprintf(meTitle, "RPCBXDistribution_for_%s", layerLabel);
    cout << "Booking for " << meId << endl;
    meMap[meId] = DBE->book1D(meId, meTitle, 11, -5.5, 5.5);

    sprintf(meId,"RPCResidualsFromCSC_%s", detUnitLabel);
    sprintf(meTitle,"RPCResidualsFromCSC_for_%s", layerLabel);
    meMap[meId] = DBE->book1D(meId, meTitle, 101, -7.*stripw, 7*stripw);

    sprintf(meId,"RPCStripYOccupancyFromCSC_%s", detUnitLabel);
    sprintf(meTitle,"RPCStripYOccupancyFromCSC_for_%s", layerLabel);
    meMap[meId] = DBE->book1D(meId, meTitle, 101, -0.5*stripl, 0.5*stripl);

    // RPC Timming 
    sprintf(meId, "RPCBX_%s", detUnitLabel);
    sprintf(meTitle, "RPCBX_for_%s", layerLabel);
    cout << "Booking for " << meId << endl;
    meMap[meId] = DBE->book1D(meId, meTitle, 11, -5.5, 5.5);

    // CSC Strip Timming 
    sprintf(meId, "averageCSCpeakTime_%s", detUnitLabel);
    sprintf(meTitle, "averageCSCpeakTime_for_%s", layerLabel);
    cout << "Booking for " << meId << endl;
    meMap[meId] = DBE->book1D(meId, meTitle, 110, -5.5, 5.5);

    // CSC Wire Timming 
    sprintf(meId, "averageCSCWireBX_%s", detUnitLabel);
    sprintf(meTitle, "averageCSCWireBX_for_%s", layerLabel);
    cout << "Booking for " << meId << endl;
    meMap[meId] = DBE->book1D(meId, meTitle, 11, -5.5, 5.5);

    // Delta Timing
    sprintf(meId, "DeltaBXTime_%s", detUnitLabel);
    sprintf(meTitle, "DeltaBXTime_for_%s", layerLabel);
    cout << "Booking for " << meId << endl;
    meMap[meId] = DBE->book1D(meId, meTitle, 11, -5.5, 5.5);

    // New 2D and more
    sprintf(meId,"ExpectedOccupancy2DFromCSC_%s", detUnitLabel);
    sprintf(meTitle,"ExpectedOccupancy2DFromCSC_for_%s", layerLabel);
    meMap[meId] = DBE->book2D(meId, meTitle, nstrips, -0.5*nstrips*stripw, 0.5*nstrips*stripw, nstrips, -0.5*stripl, 0.5*stripl);

    sprintf(meId,"RPCDataOccupancy2DFromCSC_%s", detUnitLabel);
    sprintf(meTitle,"RPCDataOccupancy2DFromCSC_for_%s", layerLabel);
    meMap[meId] = DBE->book2D(meId, meTitle, nstrips, -0.5*nstrips*stripw, 0.5*nstrips*stripw, nstrips, -0.5*stripl, 0.5*stripl);  

    return meMap;
}

// ------------ method called to for each event  ------------
void CSCSegment2RPC::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    using namespace edm;
    using namespace std;
    //char meRPCId [128];

    // Get CSC segment
    iEvent.getByLabel(CSCSegmentsTag_, pCSCSegments);
    // Get CSC strip digi collection
    iEvent.getByLabel(CSCStripDigisTag_, pCSCStripDigiCollection);
    //Get CSC wire digi collection
    iEvent.getByLabel(CSCWireDigisTag_, pCSCWireDigiCollection);
    // Get the RPC geometry
    iSetup.get<MuonGeometryRecord>().get(pRPCGeom);
    rpcGeometry = (RPCGeometry*)&*pRPCGeom;
    // Get the CSC geomtry
    iSetup.get<MuonGeometryRecord>().get(pCSCGeom);
    cscGeometry = (CSCGeometry*)&*pCSCGeom;
    // Get RPC simHits
    iEvent.getByLabel(RPCsimHitTag_, pRPCsimHits);
    // Get RPC Digis
    iEvent.getByLabel(RPCDigiTag_, pRPCDigis);
    // Get RPC recHit
    iEvent.getByLabel(RPCrecHitTag_, pRPCrecHits);
    // Get SimTrack
    //iEvent.getByLabel(SimTrackTag_, pSimTracks);

    issampled = false;
    isSamplefilled = false;
    ischecked = false;
    isEfficiencyfilled = false;

    if(sampleType == 1)
        sampleCSCSegments();

    if(sampleType == 2)
        sampleRPCsimHits();

}

void CSCSegment2RPC::sampleCSCSegments() {
    unsigned int nSegment = pCSCSegments->size();
    cout << "This event has " << nSegment << " CSCSegments" << endl;
    Statistics->Fill(nSegment);

    std::map<CSCDetId, int> CSCSegmentsCounter;
    for(CSCSegmentCollection::const_iterator CSCSegIter = pCSCSegments->begin(); CSCSegIter != pCSCSegments->end(); CSCSegIter++)
        CSCSegmentsCounter[CSCSegIter->cscDetId()]++;
    for(CSCSegmentCollection::const_iterator CSCSegIter = pCSCSegments->begin(); CSCSegIter != pCSCSegments->end(); CSCSegIter++) {
        issampled = false;
        isSamplefilled = false;
        ischecked = false;
        isEfficiencyfilled = false;
                    
        CSCDetId CSCId = CSCSegIter->cscDetId();
        int Endcap = CSCId.endcap();
        int Station = CSCId.station();
        int Ring = CSCId.ring();
        int Chamber = CSCId.chamber();
        int Layer = CSCId.layer();
        cout << "CSC endcap is " << Endcap << ", Station is " << Station << ", Ring is " << Ring << " Chamber is " << Chamber << ", Layer is " << Layer << endl;
        // 1=Forward(+), 2=Backward(-)
        //if(Endcap != 1)
        //continue;
        if(Station == 4)
            continue;
        if(Ring == 1 || Ring == 4)
            continue;
        if(isCSCSegmentFilter == true && !filterCSCSegment(*CSCSegIter))
            continue;
        if(!filterSegmentEta(*CSCSegIter))
            continue;
        // One chamber has only 1 segment to avoid electron shower segments, and require at least 2 segments for this event
        if(CSCSegmentsCounter[CSCId]==1 && pCSCSegments->size()>=2) {    
            LocalPoint SegmentLocalPosition= CSCSegIter->localPosition();
            LocalVector SegmentLocalDirection=CSCSegIter->localDirection();
            cout << "CSC segment local position: " << SegmentLocalPosition << endl;
            cout << "CSC segment local direction: " << SegmentLocalDirection << endl;
            const GeomDet* CSCDet = cscGeometry->idToDet(CSCId);
            if(dynamic_cast<const CSCChamber*>(CSCDet) != 0) {
                const CSCChamber* CSCCh = dynamic_cast<const CSCChamber*>(CSCDet);
                GlobalPoint SegmentGlobalPosition = CSCCh->toGlobal(SegmentLocalPosition);
                GlobalVector SegmentGlobalDirection = CSCCh->toGlobal(SegmentLocalDirection);
                cout << "CSC segment global position: " << SegmentGlobalPosition << endl;
                cout << "CSC segment global direction: " << SegmentGlobalDirection << endl;
                // Find correspond RPC rolls
                if(RollswithCSCChamber.find(CSCId) != RollswithCSCChamber.end()) {
                    std::vector<RPCDetId> RPCRolls = RollswithCSCChamber[CSCId];
                    for(std::vector<RPCDetId>::const_iterator RPCRollIter = RPCRolls.begin(); RPCRollIter != RPCRolls.end(); RPCRollIter++) {
                        issampled = false;
                        const GeomDetUnit* RPCDet = rpcGeometry->idToDetUnit(*RPCRollIter);
                        if(dynamic_cast<const RPCRoll*>(RPCDet) != 0) {
                            const RPCRoll* sampleRPCRoll = dynamic_cast<const RPCRoll*>(RPCDet);
                            const BoundPlane& sampleRPCRollSurface = sampleRPCRoll->surface();
                            const StripTopology& sampleRPCRollTop = sampleRPCRoll->specificTopology();

                            GlobalPoint test(1, 1, 1);
                            LocalPoint test1 = sampleRPCRoll->toLocal(test);
                            LocalPoint test2 = sampleRPCRollSurface.toLocal(test);
                            cout << "For test: RPC roll local position is :" << test1 << endl;
                            cout << "For test: RPC roll surface local position is :" << test2 << endl;

                            LocalPoint SegmentRPCPosition = sampleRPCRollSurface.toLocal(SegmentGlobalPosition);
                            LocalVector SegmentRPCDirection = sampleRPCRollSurface.toLocal(SegmentGlobalDirection);
                            double X0 = SegmentRPCPosition.x();
                            double Y0 = SegmentRPCPosition.y();
                            double Z0 = SegmentRPCPosition.z();
                            double dx = SegmentRPCDirection.x();
                            double dy = SegmentRPCDirection.y();
                            double dz = SegmentRPCDirection.z();
                            //for(unsigned int iStrip= 1; iStrip <= RPCRo->nstrips(); iStrip++)
                            //cout << "Strip " << iStrip << " local position is " << RPCRo->centreOfStrip((int)iStrip) << endl;
                            double X = X0 + (-Z0) * dx / dz;
                            double Y = Y0 + (-Z0) * dy / dz;
                            LocalPoint impactRPCPosition(X, Y, 0);

                            double S = (SegmentRPCPosition - impactRPCPosition).mag();
                            if(fabs(S) > MaxD) {
                                cout << "Extropolate distance is too long!" << endl;
                                continue;
                            }
                            LocalPoint EdgePosition0 = sampleRPCRollTop.localPosition(0.);
                            LocalPoint EdgePosition1 = sampleRPCRollTop.localPosition((double)(sampleRPCRollTop.nstrips()));
                            cout << "For Roll topology EdgePosition 0 local position is: " << EdgePosition0 << endl;
                            cout << "For Roll topology EdgePosition 1 local position is: " << EdgePosition0 << endl;
                            RollSize = fabs(EdgePosition0.x()-EdgePosition1.x());
                            StripLength = sampleRPCRollTop.stripLength();
                            StripWidth = sampleRPCRollTop.pitch();
                            if(fabs(impactRPCPosition.z()) <= 1 && fabs(impactRPCPosition.x()) <= RollSize*0.5 && fabs(impactRPCPosition.y()) <= StripLength*0.5) {
                                cout << "CSC chamber: " << CSCId.rawId() << ", Endcap: " << Endcap << ", Station: " << Station << ", Ring: " << Ring << ", Chamber: " << Chamber << ", layer: " << Layer << endl;

                                impactRPCRollId = *RPCRollIter; 
                                impactDirection = SegmentRPCDirection; 
                                impactPoint = impactRPCPosition; 
                                impactStrip =  sampleRPCRollTop.strip(impactPoint);                                
                                findpeakTime(*CSCSegIter);

                                // For Ntuple
                                impactRPCId_NTuple = impactRPCRollId.rawId();
                                impactCSCBX_NTuple = CSCWireTimeBin;
                                impactAngle_NTuple = ((impactDirection.theta() > PI/2) ? (PI- impactDirection.theta()) : impactDirection.theta());
                                impactStrip_NTuple = (int)impactStrip+1;
                                cout << "impactRPCId_NTuple is " << impactRPCId_NTuple << ", impactCSCBX_NTuple is " << impactCSCBX_NTuple << ", impactAngle_NTuple is " << impactAngle_NTuple << ", impactStrip_NTuple is " << impactStrip_NTuple << endl;
                                // check efficiency
                                issampled = true;
                                fillSample();
                                EfficiencybyRPCrecHit();
                                // avoid in the boarder of 2 chambers efficiency is calculated for 2 times
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
}

void CSCSegment2RPC::findpeakTime(const CSCSegment& sampleCSCSegment) {

    std::vector<CSCRecHit2D> sampleCSCrecHitDGroup = sampleCSCSegment.specificRecHits();
    // Find the 1st and 2nd peak time
    double peakTime1st = 10000.;
    double peakTime2nd = 5000.;
    double peakTimesum = 0.;

    int WireTimeBinSum = 0;
    int averageWireTimeBin = 0;

    for(std::vector<CSCRecHit2D>::iterator CSCrecHitDIter = sampleCSCrecHitDGroup.begin(); CSCrecHitDIter != sampleCSCrecHitDGroup.end(); CSCrecHitDIter++) {
        double temppeakTime = CSCrecHitDIter->tpeak();
        peakTimesum += temppeakTime;
        if(temppeakTime < peakTime1st) {
            peakTime1st = temppeakTime;
            peakTime2nd = peakTime1st;
        }
        if(temppeakTime >= peakTime1st && temppeakTime < peakTime2nd)
            peakTime2nd = temppeakTime;

        const CSCDetId sampleCSCLayerId = CSCrecHitDIter->cscDetId();
        CSCStripDigiCollection::Range sampleCSCStripDigiRange = pCSCStripDigiCollection->get(sampleCSCLayerId);
        const std::vector<int> sampleCSCStripCluster = CSCrecHitDIter->channels();
        for(CSCStripDigiCollection::const_iterator sampleCSCStripDigiIter = sampleCSCStripDigiRange.first; sampleCSCStripDigiIter != sampleCSCStripDigiRange.second; sampleCSCStripDigiIter++) {
            int sampleStrip = sampleCSCStripDigiIter->getStrip();
            const std::vector<int>::const_iterator resultStripIter = find(sampleCSCStripCluster.begin(), sampleCSCStripCluster.end(), sampleStrip);
            if(resultStripIter != sampleCSCStripCluster.end()) {
            }
        }

        // Wire time bin
        int firstWireGroup = 10000;
        int firstWireTimeBin = -1;
        CSCWireDigiCollection::Range sampleCSCWireDigiRange = pCSCWireDigiCollection->get(sampleCSCLayerId);
        const std::vector<int> sampleCSCWireCluster = CSCrecHitDIter->wgroups();
        cout << "Find " << sampleCSCWireCluster.size() << " wire groups in the cluster" << endl;
        for(CSCWireDigiCollection::const_iterator sampleCSCWireDigiIter = sampleCSCWireDigiRange.first; sampleCSCWireDigiIter != sampleCSCWireDigiRange.second; sampleCSCWireDigiIter++) {
            int sampleWireGroup = sampleCSCWireDigiIter->getWireGroup();
            const std::vector<int>::const_iterator resultWireIter = find(sampleCSCWireCluster.begin(), sampleCSCWireCluster.end(), sampleWireGroup);
            if(resultWireIter != sampleCSCWireCluster.end()) {
                std::vector<int> sampleWireTimeBins = sampleCSCWireDigiIter->getTimeBinsOn();
                for(std::vector<int>::iterator Iter = sampleWireTimeBins.begin(); Iter != sampleWireTimeBins.end(); Iter++)
                    cout << "Sample wire BX time bin word is: " << (*Iter) << endl;
                // The 1st time bin in the wire
                int sampleWireTimeBin = (*(sampleWireTimeBins.begin()));
                if(sampleWireGroup < firstWireGroup) {
                    firstWireGroup = sampleWireGroup;
                    firstWireTimeBin = sampleWireTimeBin;
                }
                cout << "Sample wire BX time bin is: " << firstWireTimeBin << endl;
            }
        }
        WireTimeBinSum += firstWireTimeBin;
    }

    CSCpeakTime2nd = peakTime2nd;
    averageCSCpeakTime = peakTimesum / sampleCSCrecHitDGroup.size();
    cout << "CSCpeakTime2nd is " << CSCpeakTime2nd << ", averageCSCpeakTime is " << averageCSCpeakTime << endl;

    averageWireTimeBin = WireTimeBinSum / sampleCSCrecHitDGroup.size();
    CSCWireTimeBin = averageWireTimeBin - 6;
    cout << "CSC average wire time bin is " << averageWireTimeBin << ", CSC wire BX is " << CSCWireTimeBin << endl;

}

void CSCSegment2RPC::sampleRPCsimHits() {

    for(PSimHitContainer::const_iterator simHitIter = pRPCsimHits->begin(); simHitIter != pRPCsimHits->end(); simHitIter++) {
        // The particle type of the hit may differ from the particle type of the SimTrack with id trackId(). 
        // This happends if the hit was created by a secondary track (e.g. a delta ray) originating from the trackId() and not existing as a separate SimTrack.
        int Particletype = simHitIter->particleType();
        if(Particletype != TrackType)
            continue;

        // Check the DetId is a RPCDetId
        unsigned int sampleDetUnitId = simHitIter->detUnitId();
        DetId sampleDetectorId = DetId(sampleDetUnitId);
        GeomDetEnumerators::SubDetector subDet = rpcGeometry->idToDetUnit(sampleDetectorId)->subDetector();
        if(subDet != GeomDetEnumerators::RPCBarrel && subDet != GeomDetEnumerators::RPCEndcap)        
            continue;

        // Check if inside the collection range
        const RPCDetId sampleRPCRollId(sampleDetectorId);
        if(meCollection.find(sampleRPCRollId) == meCollection.end())
            continue;

        issampled = false;
        isSamplefilled = false;
        ischecked = false;
        isEfficiencyfilled = false;

        const GeomDetUnit* RPCDet = rpcGeometry->idToDetUnit(sampleRPCRollId);
        if(dynamic_cast<const RPCRoll*>(RPCDet) != 0) {
            const RPCRoll* sampleRPCRoll = dynamic_cast<const RPCRoll*>(RPCDet);
            const BoundPlane& sampleRPCRollSurface = sampleRPCRoll->surface();
            const StripTopology& sampleRPCRollTop = sampleRPCRoll->specificTopology();
            Local3DPoint sampleLocalPosition = simHitIter->localPosition();
            LocalPoint impactRPCPosition = sampleLocalPosition;
            LocalVector SegmentRPCDirection = simHitIter->momentumAtEntry(); 
            LocalPoint EdgePosition0 = sampleRPCRollTop.localPosition(0.);
            LocalPoint EdgePosition1 = sampleRPCRollTop.localPosition((double)(sampleRPCRollTop.nstrips()));
            cout << "For Roll topology EdgePosition 0 local position is: " << EdgePosition0 << endl;
            cout << "For Roll topology EdgePosition 1 local position is: " << EdgePosition0 << endl;
            RollSize = fabs(EdgePosition0.x()-EdgePosition1.x());
            StripLength = sampleRPCRollTop.stripLength();
            StripWidth = sampleRPCRollTop.pitch();
            if(fabs(impactRPCPosition.z()) <= 1 && fabs(impactRPCPosition.x()) <= RollSize*0.5 && fabs(impactRPCPosition.y()) <= StripLength*0.5) {

                impactRPCRollId = sampleRPCRollId;
                impactDirection = SegmentRPCDirection;
                impactPoint = impactRPCPosition;
                impactStrip =  sampleRPCRollTop.strip(impactPoint);
                CSCWireTimeBin = 0;

                // For Ntuple 
                impactRPCId_NTuple = impactRPCRollId.rawId(); 
                impactCSCBX_NTuple = CSCWireTimeBin; 
                impactAngle_NTuple = ((impactDirection.theta() > PI/2) ? (PI- impactDirection.theta()) : impactDirection.theta());
                impactStrip_NTuple = (int)impactStrip+1;


                // check efficiency
                issampled = true;
                fillSample();
                EfficiencybyRPCrecHit();
            }
        }
    }
}


bool CSCSegment2RPC::filterCSCSegment(const CSCSegment& sampleCSCSegment) {

    bool isTrackSegment = false;
    CSCDetId sampleCSCId = sampleCSCSegment.cscDetId();
    const CSCChamber* sampleCSCChamber = cscGeometry->chamber(sampleCSCId);
    LocalVector localsampleVector = sampleCSCSegment.localDirection();
    LocalPoint localsamplePosition = sampleCSCSegment.localPosition();
    GlobalVector globalsampleVector = sampleCSCChamber->toGlobal(localsampleVector);
    GlobalPoint globalsamplePosition = sampleCSCChamber->toGlobal(localsamplePosition);
    int sampleCSCStation = sampleCSCId.station();
    int sampleCSCzEndcap = sampleCSCId.zendcap();
    // For CSCDetId::endcap(), 1=forward (+Z); 2=backward (-Z) 
    //if(sampleCSCEndcap == 2)
    //sampleCSCEndcap = -1;

    for(CSCSegmentCollection::const_iterator CSCSegIter = pCSCSegments->begin(); CSCSegIter != pCSCSegments->end(); CSCSegIter++) {
        CSCDetId iterCSCId = CSCSegIter->cscDetId();
        int iterCSCStation = iterCSCId.station();
        int iterCSCzEndcap = iterCSCId.zendcap();
        if((abs(sampleCSCStation-iterCSCStation) != 1) && (sampleCSCzEndcap != iterCSCzEndcap))
            continue;

        const CSCChamber* iterCSCChamber = cscGeometry->chamber(iterCSCId);
        LocalVector localIterVector = CSCSegIter->localDirection();
        LocalPoint localIterPosition = CSCSegIter->localPosition();
        GlobalVector globalIterVector = iterCSCChamber->toGlobal(localIterVector);
        GlobalPoint globalIterPosition = iterCSCChamber->toGlobal(localIterPosition);

        if(deltaRforSegment(globalsamplePosition, globalsampleVector, globalIterPosition, globalIterVector) <= deltaRTH)
            isTrackSegment = true;
    }
    if(isTrackSegment == true)
        cout << "Pass Segment filter. successfully." << endl;
    else
        cout << "Pass Segment filter unsuccessfully." << endl;

    return isTrackSegment;
}

double CSCSegment2RPC::deltaRforSegment(const GlobalPoint& globalPosition1, const GlobalVector& globalVector1, const GlobalPoint& globalPosition2, const GlobalVector& globalVector2) {
    GlobalVector VectorPoint2Point = (GlobalVector)(globalPosition2 - globalPosition1);
    double deltaR1 = sqrt((globalVector2.phi() - globalVector1.phi()).value() * (globalVector2.phi() - globalVector1.phi()).value() + (globalVector2.theta() - globalVector1.theta()) * (globalVector2.theta() - globalVector1.theta()));
    double deltaR2 = sqrt((VectorPoint2Point.phi() - globalVector1.phi()).value() * (VectorPoint2Point.phi() - globalVector1.phi()).value() + (VectorPoint2Point.theta() - globalVector1.theta()) * (VectorPoint2Point.theta() - globalVector1.theta()));
    return (deltaR1 > deltaR2) ? deltaR1 : deltaR2;
}

bool CSCSegment2RPC::filterSegmentEta(const CSCSegment& sampleCSCSegment) {

    bool isinsideEtaRange = false;
    CSCDetId sampleCSCId = sampleCSCSegment.cscDetId();
    const CSCChamber* sampleCSCChamber = cscGeometry->chamber(sampleCSCId);
    LocalVector localsampleVector = sampleCSCSegment.localDirection();
    GlobalVector globalsampleVector = sampleCSCChamber->toGlobal(localsampleVector);
    double sampleSegmentEta = globalsampleVector.eta();
    if((sampleSegmentEta >= minSegmentEta) && (sampleSegmentEta <= maxSegmentEta))
        isinsideEtaRange = true;

    if(isinsideEtaRange == true)
        cout << "Pass Eta filter. successfully." << endl;
    else
        cout << "Pass Eta filter unsuccessfully." << endl;

    return isinsideEtaRange;
}

void CSCSegment2RPC::fillSample() {

    if(issampled == false || isSamplefilled == true) {
        cout << "Not yet sample the RPC rolls or sample already filled" << endl;
        return;
    }

    int  sampleEndcap = impactRPCRollId.region();
    int  sampleStation = impactRPCRollId.station();
    int  sampleRing = impactRPCRollId.ring();
    int  sampleRoll = impactRPCRollId.roll();
    cout << "Find a impact point from this CSC segment onto RPC roll: " << impactRPCRollId.rawId() << ", Endcap: " << sampleEndcap << ", Station: " << sampleStation << ", Ring: " << sampleRing << ", Roll: " << sampleRoll << endl;
    cout << "The impact strip is " << impactStrip << ", impact point on RPC frame: " << impactPoint << endl;

    if(meCollection.find(impactRPCRollId) != meCollection.end()) {
        std::map<std::string, MonitorElement*> meMap = meCollection[impactRPCRollId];
        RPCGeomServ RPCname(impactRPCRollId);
        std::string nameRoll = RPCname.name();
        // Fill CSC Strip peak Time
        sprintf(meRPCId,"averageCSCpeakTime_%s", nameRoll.c_str());
        cout <<"Filling to histogram " << meRPCId << endl;
        if(meMap.find(meRPCId) != meMap.end())
            meMap[meRPCId]->Fill(averageCSCpeakTime);
        else
            cout << "Could not find a histogram " << meRPCId << endl;

        // Fill CSC Wire peak Time
        sprintf(meRPCId,"averageCSCWireBX_%s", nameRoll.c_str());
        cout <<"Filling to histogram " << meRPCId << endl;
        if(meMap.find(meRPCId) != meMap.end())
            meMap[meRPCId]->Fill(CSCWireTimeBin);
        else
            cout << "Could not find a histogram " << meRPCId << endl;

        // Fill extrapolation histogram 
        // N.x strip position correspond to N+1 strip in RPC local frame
        sprintf(meRPCId, "ExpectedOccupancyFromCSC_%s", nameRoll.c_str());
        cout <<"Filling to histogram " << meRPCId << endl;
        if(meMap.find(meRPCId) != meMap.end())
            meMap[meRPCId]->Fill((int)impactStrip+1);
        else
            cout << "Could not find a histogram " << meRPCId << endl;

        // Fill extrapolation histogram
        sprintf(meRPCId, "RPCStripYOccupancyFromCSC_%s", nameRoll.c_str());
        cout <<"Filling to histogram " << meRPCId << endl;
        if(meMap.find(meRPCId) != meMap.end())
            meMap[meRPCId]->Fill(impactPoint.y());
        else
            cout << "Could not find a histogram " << meRPCId << endl;

        // Fill extrapolation 2D histogram
        sprintf(meRPCId, "ExpectedOccupancy2DFromCSC_%s", nameRoll.c_str());
        cout <<"Filling to histogram " << meRPCId << endl;
        if(meMap.find(meRPCId) != meMap.end())
            meMap[meRPCId]->Fill(impactPoint.x(), impactPoint.y());
        else
            cout << "Could not find a histogram " << meRPCId << endl;
    }
    else
        cout << "Could not find a RPC roll index from meCollection" << endl;

    isSamplefilled = true;
}


// Try to find the fire strip from impact RPC roll
void CSCSegment2RPC::EfficiencybyRPCrecHit() {

    if(issampled == false || ischecked == true) {
        cout << "Not yet sample the RPC rolls or already checked efficiency" << endl;
        return;        
    }


    double residualDistance = 10000;
    LocalPoint FirePosition;
    int FireClusterSize = 0;
    cout << "Try to find RPC strips closest to impact point" << endl;

    anycoincidence = false;
    isvalideRPCBX = false;
    isEfficiency_NTuple = 0;

    RPCRecHitCollection::range RollrecHitsRange = pRPCrecHits->get(impactRPCRollId);
    RPCBX = 0;
    for(RPCRecHitCollection::const_iterator recHitIter = RollrecHitsRange.first; recHitIter != RollrecHitsRange.second; recHitIter++) {
        //const RPCRecHit& = *recHitIter;
        LocalPoint candidatePosition = recHitIter->localPosition();
        double residualDistance_Temp = candidatePosition.x()-impactPoint.x();
        int candiadteClusterSize = recHitIter->clusterSize();
        int candidateBX = recHitIter->BunchX();
        cout << "ResidualDistance_Temp is " << residualDistance_Temp << endl;
        if(fabs(residualDistance_Temp) < fabs(residualDistance)) {
            residualDistance = residualDistance_Temp;
            FirePosition = candidatePosition;
            FireClusterSize = candiadteClusterSize;
            RPCBX = candidateBX;
            deltaBXTime = CSCWireTimeBin - RPCBX;
            cout << "Step into a candidate strip. residualDistance: " << residualDistance << ", FireClusterSize: " << FireClusterSize << endl;
        }
    }
    if(fabs(residualDistance) <= (rangestrips+FireClusterSize*0.5)*StripWidth) {
        cout << "Find a fire strip with-in the residual range in RPC roll " << impactRPCRollId.rawId() << endl;
        // For hist
        anycoincidence = true;
        isvalideRPCBX = true;
        impactResidual = residualDistance;
        impactClusterSize = FireClusterSize;
        // For tree
        isEfficiency_NTuple = 1;
        impactResidual_NTuple = residualDistance;
        impactClusterSize_NTuple = FireClusterSize;
        impactRPCBX_NTuple = RPCBX;
    }
    else {
        cout <<"Could not find a fire strip with-in the residual range in RPC roll " << impactRPCRollId.rawId() << endl;
        // For hist
        anycoincidence = false;
        isvalideRPCBX = false;
        impactResidual = residualDistance;
        impactClusterSize = FireClusterSize;
        // For tree
        isEfficiency_NTuple = 0;
        impactResidual_NTuple = residualDistance;
        impactClusterSize_NTuple = FireClusterSize;
        impactRPCBX_NTuple = RPCBX;
    }
    ischecked = true;
    fillEfficiency();
}

void CSCSegment2RPC::fillEfficiency() {

    if(ischecked == false || isEfficiencyfilled == true) {
        cout << "Not yet check the RPC roll or efficiency already filled." << endl;
        return;
    }

    if(anycoincidence == true && isvalideRPCBX == true) {
        if(impactRPCRollId.ring()==2 && impactRPCRollId.roll()==1) {
            cout << "Ring 2 Roll A cluster size: " << impactClusterSize << endl;
            if(impactClusterSize == 1*dupli)
                hGlobalResClu1R2A->Fill(impactResidual); 
            if(impactClusterSize == 2*dupli) 
                hGlobalResClu2R2A->Fill(impactResidual); 
            if(impactClusterSize == 3*dupli) 
                hGlobalResClu3R2A->Fill(impactResidual);
        }
        if(impactRPCRollId.ring()==2 && impactRPCRollId.roll()==2) {
            cout << "Ring 2 Roll B cluster size: " << impactClusterSize << endl;
            if(impactClusterSize == 1*dupli) 
                hGlobalResClu1R2B->Fill(impactResidual); 
            if(impactClusterSize == 2*dupli) 
                hGlobalResClu2R2B->Fill(impactResidual); 
            if(impactClusterSize == 3*dupli) 
                hGlobalResClu3R2B->Fill(impactResidual);
        }
        if(impactRPCRollId.ring()==2 && impactRPCRollId.roll()==3) {
            cout << "Ring 2 Roll C cluster size: " << impactClusterSize << endl;
            if(impactClusterSize == 1*dupli) 
                hGlobalResClu1R2C->Fill(impactResidual); 
            if(impactClusterSize == 2*dupli) 
                hGlobalResClu2R2C->Fill(impactResidual); 
            if(impactClusterSize == 3*dupli) 
                hGlobalResClu3R2C->Fill(impactResidual);
        }
        if(impactRPCRollId.ring()==3 && impactRPCRollId.roll()==1) {
            cout << "Ring 3 Roll A cluster size: " << impactClusterSize << endl;
            if(impactClusterSize == 1*dupli) 
                hGlobalResClu1R3A->Fill(impactResidual); 
            if(impactClusterSize == 2*dupli) 
                hGlobalResClu2R3A->Fill(impactResidual); 
            if(impactClusterSize == 3*dupli) 
                hGlobalResClu3R3A->Fill(impactResidual);
        }
        if(impactRPCRollId.ring()==3 && impactRPCRollId.roll()==2) {
            cout << "Ring 3 Roll B cluster size: " << impactClusterSize << endl;
            if(impactClusterSize == 1*dupli) 
                hGlobalResClu1R3B->Fill(impactResidual); 
            if(impactClusterSize == 2*dupli) 
                hGlobalResClu2R3B->Fill(impactResidual); 
            if(impactClusterSize == 3*dupli) 
                hGlobalResClu3R3B->Fill(impactResidual);
        }
        if(impactRPCRollId.ring()==3 && impactRPCRollId.roll()==3) {
            cout << "Ring 3 Roll C cluster size: " << impactClusterSize << endl;
            if(impactClusterSize == 1*dupli) 
                hGlobalResClu1R3C->Fill(impactResidual); 
            if(impactClusterSize == 2*dupli) 
                hGlobalResClu2R3C->Fill(impactResidual); 
            if(impactClusterSize == 3*dupli) 
                hGlobalResClu3R3C->Fill(impactResidual);
        }

        if(meCollection.find(impactRPCRollId) != meCollection.end()) {
            std::map<std::string, MonitorElement*> meMap = meCollection[impactRPCRollId];
            RPCGeomServ RPCname(impactRPCRollId);
            std::string nameRoll = RPCname.name();
            // Fill BX 
            sprintf(meRPCId,"RPCBX_%s", nameRoll.c_str());
            cout <<"Filling to histogram " << meRPCId << endl;
            if(meMap.find(meRPCId) != meMap.end())
                meMap[meRPCId]->Fill(RPCBX);
            else
                cout << "Could not find a histogram " << meRPCId << endl;
            // Fill Delta Timing
            sprintf(meRPCId,"DeltaBXTime_%s", nameRoll.c_str());
            cout <<"Filling to histogram " << meRPCId << endl;
            if(meMap.find(meRPCId) != meMap.end())
                meMap[meRPCId]->Fill(deltaBXTime);
            else
                cout << "Could not find a histogram " << meRPCId << endl;
            // Fill Residuals histogram
            sprintf(meRPCId, "RPCResidualsFromCSC_%s", nameRoll.c_str());
            cout <<"Filling to histogram " << meRPCId << endl;
            if(meMap.find(meRPCId) != meMap.end())
                meMap[meRPCId]->Fill(impactResidual);
            else
                cout << "Could not find a histogram " << meRPCId << endl;
            // Fill Occupancy histogram
            // N.x strip position correspond to N+1 strip in RPC local frame
            sprintf(meRPCId, "RPCDataOccupancyFromCSC_%s", nameRoll.c_str());
            cout <<"Filling to histogram " << meRPCId << endl;
            if(meMap.find(meRPCId) != meMap.end())
                meMap[meRPCId]->Fill((int)impactStrip+1);
            else
                cout << "Could not find a histogram " << meRPCId << endl;
            // Fill Occupancy 2D histogram
            sprintf(meRPCId, "RPCDataOccupancy2DFromCSC_%s", nameRoll.c_str());
            cout <<"Filling to histogram " << meRPCId << endl;
            if(meMap.find(meRPCId) != meMap.end())
                meMap[meRPCId]->Fill(impactPoint.x(), impactPoint.y());
            else
                cout << "Could not find a histogram " << meRPCId << endl;
        }
        else 
            cout << "Could not find a RPC roll index from meCollection" << endl;
    }
    // Fill tree whatever efficiency is 0 or 1;
    ExTree->Fill();

    isEfficiencyfilled = true;
}


// ------------ method called once each job just before starting event loop  ------------
void CSCSegment2RPC::beginJob(const edm::EventSetup& iSetup) {

    edm::ESHandle<RPCGeometry> pRPCGeom;
    iSetup.get<MuonGeometryRecord>().get(pRPCGeom);
    const RPCGeometry* rpcGeometry = (const RPCGeometry*)&*pRPCGeom;
    edm::ESHandle<CSCGeometry> pCSCGeom;
    iSetup.get<MuonGeometryRecord>().get(pCSCGeom);
    const CSCGeometry* cscGeometry = (const CSCGeometry*)&*pCSCGeom;

    const std::vector<RPCRoll*> rpcRolls= rpcGeometry->rolls();
    for(std::vector<RPCRoll*>::const_iterator RPCRollIter = rpcRolls.begin(); RPCRollIter != rpcRolls.end(); RPCRollIter++) {
        RPCDetId rpcRoll = (*RPCRollIter)->id();
        RPCGeomServ rpcsrv(rpcRoll);
        int Region = rpcRoll.region();
        if(Region == 0)
            continue;
        // For RPC Region is +1/-1 for +Z/-Z, while for CSC endcap is 1/2 for +Z/-Z
        if(Region == -1)
            Region = 2;
        int Station = rpcRoll.station();
        if(Station == 4)
            continue;
        int Ring = rpcRoll.ring();
        if(Ring == 1)
            continue;
        if((Station == 2 || Station == 3) && Ring == 3)
            Ring = 2;
        // correspond CSC chamber??
        int Chamber = rpcsrv.segment();
        int Layer = 0;
        CSCDetId CSCIndex(Region, Station, Ring, Chamber, Layer);
        std::vector<RPCDetId> RPCRolls_temp;
        if(RollswithCSCChamber.find(CSCIndex) != RollswithCSCChamber.end())
            RPCRolls_temp = RollswithCSCChamber[CSCIndex];
        RPCRolls_temp.push_back(rpcRoll);
        RollswithCSCChamber[CSCIndex] = RPCRolls_temp;
        meCollection[rpcRoll] = bookDetUnitSeg(rpcRoll, (*RPCRollIter)->nstrips(), (*RPCRollIter)->pitch(), (*RPCRollIter)->specificTopology().stripLength());
    }

    const std::vector<CSCChamber*> cscChambers = cscGeometry->chambers();
    for(std::vector<CSCChamber*>::const_iterator CSCChamberIter = cscChambers.begin(); CSCChamberIter != cscChambers.end(); CSCChamberIter++) {
        CSCDetId cscChamber = (*CSCChamberIter)->id();
        int Endcap = cscChamber.endcap();
        int Station = cscChamber.station();
        int Ring = cscChamber.ring();
        int Chamber = cscChamber.chamber();
        int Layer = cscChamber.layer();
        // 1=Forward(+), 2=Backward(-)
        //if(Endcap != 1)
        //continue;
        if(Station == 4)
            continue;
        if(Ring == 1 || Ring == 4)
            continue;
        cout << "Checking endcap CSC chamber: " << cscChamber.rawId() << ", Endcap: " << Endcap << ", Station: " << Station << ", Ring: " << Ring << ", Chamber: " << Chamber << ", layer: " << Layer << endl;
        const Bounds& cscBounds = (*CSCChamberIter)->surface().bounds();
        const CSCLayer* cscLayer = (*CSCChamberIter)->layer(1);
        const CSCStripTopology* cscTop = cscLayer->geometry()->topology();
        int nCSCStrip0 = cscLayer->geometry()->numberOfStrips();
        int nCSCStrip1 = cscTop->nstrips();
        cout << "CSC layer strips from layer geometry: " << nCSCStrip0 << ", from layer topology: " << nCSCStrip1 << endl;
        GlobalPoint Edge0 = cscLayer->toGlobal(cscTop->localPosition(0));
        GlobalPoint Edge1 = cscLayer->toGlobal(cscTop->localPosition(nCSCStrip1));
        cout << "CSC Phi range from " << Edge0.phi() << " to " << Edge1.phi() << endl;
        cout << "CSC Phi degree range from " << Edge0.phi().degrees() << " to " << Edge1.phi().degrees() << endl;
        std::vector<RPCDetId> RPCRolls_temp;
        if(RollswithCSCChamber.find(cscChamber) != RollswithCSCChamber.end()) {
            std::vector<RPCDetId> RPCRolls;
            RPCRolls = RollswithCSCChamber[cscChamber];
            cout <<"Bind to " << RPCRolls.size() << " RPC Rolls" << endl;
            for(std::vector<RPCDetId>::const_iterator RPCRollIter = RPCRolls.begin(); RPCRollIter != RPCRolls.end(); RPCRollIter++) {
                RPCGeomServ rpcSrv(*RPCRollIter);
                int rEndcap = RPCRollIter->region();
                int rStation = RPCRollIter->station();
                int rRing = RPCRollIter->ring();
                int rSegment = rpcSrv.segment();
                int rRoll = RPCRollIter->roll();
                const RPCRoll* rpcRoll = rpcGeometry->roll(*RPCRollIter);
                int nRPCStrip = rpcRoll->nstrips();
                GlobalPoint insideEdge0 = rpcRoll->toGlobal(rpcRoll->centreOfStrip(0));
                GlobalPoint insideEdge1 = rpcRoll->toGlobal(rpcRoll->centreOfStrip(nRPCStrip));
                cout << "RPC roll: " << RPCRollIter->rawId() <<" Phi range from " << insideEdge0.phi() << " to " << insideEdge1.phi() << endl;
                cout << "RPC roll: " << RPCRollIter->rawId() <<" Phi degree range from " << insideEdge0.phi().degrees() << " to " << insideEdge1.phi().degrees() << endl;
                LocalPoint inside0 = LocalPoint((*CSCChamberIter)->surface().toLocal(insideEdge0).x(), (*CSCChamberIter)->surface().toLocal(insideEdge0).y(), 0);
                LocalPoint inside1 = LocalPoint((*CSCChamberIter)->surface().toLocal(insideEdge1).x(), (*CSCChamberIter)->surface().toLocal(insideEdge1).y(), 0);
                //if(cscBounds.inside(inside0) && cscBounds.inside(inside1))
                if((Edge0.phi()-insideEdge0.phi()).value()*(Edge1.phi()-insideEdge1.phi()).value() < 0)
                    cout << "Well binding for RPC Roll: " << RPCRollIter->rawId() << ", Endcap: " << rEndcap << ", Station: " << rStation << ", Ring: " << rRing << ", Segment: " << rSegment << ", Roll: " << rRoll << endl;
                else
                    cout << "Wrong binding for RPC Roll " << RPCRollIter->rawId() << ", Endcap: " << rEndcap << ", Station: " << rStation << ", Ring: " << rRing << ", Segment: " << rSegment << ", Roll: " << rRoll << endl;
            }
        }
        else
            cout << "Could not find the binding RPC roll" << endl;
    }
}

// ------------ method called once each job just after ending the event loop  ------------
void CSCSegment2RPC::endJob() {
    DBE = 0;
    if(SavetheRootFile) {
        std::cout << "writing information into file: " << theFile->GetName() << std::endl;
        theFile->Write();
        theFile->Close();
    }
}

void CSCSegment2RPC::endRun(const edm::Run& Run, const edm::EventSetup& iSetup)
{
    if(EffSaveRootFile)
        DBE->save(EffRootFileName);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CSCSegment2RPC);
