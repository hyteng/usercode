// -*- C++ -*-
//
// Package:    CSCSegment2RPC
// Class:      CSCSegment2RPC
// 
/**\class CSCSegment2RPC CSCSegment2RPC.cc MyModule/CSCSegment2RPC/src/CSCSegment2RPC.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
 */
//
// Original Author:  Haiyun Teng
//         Created:  Wed Feb 11 10:40:41 CET 2009
// $Id: CSCSegment2RPC.cc,v 1.10 2012/03/02 14:53:57 hyteng Exp $
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
#include <FWCore/Utilities/interface/InputTag.h>
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
#include <DQM/RPCMonitorDigi/interface/RPCBookFolderStructure.h>
#include <SimDataFormats/TrackingHit/interface/PSimHitContainer.h>
#include <DataFormats/RPCDigi/interface/RPCDigiCollection.h>
#include <SimDataFormats/Track/interface/SimTrackContainer.h>
#include <Geometry/CommonDetUnit/interface/GeomDetUnit.h>
#include <DataFormats/Provenance/interface/EventID.h>
#include "MyModule/CSCSegment2RPC/src/CSC2RPCBinder.h"

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

#ifndef debug
#define debug 1
#endif

using namespace edm;
using namespace std;
//
// class decleration
//
class CSCSegment2RPC : public edm::EDAnalyzer {

    public:
        explicit CSCSegment2RPC(const edm::ParameterSet& iConfig);
        ~CSCSegment2RPC();

    private:
        virtual void beginRun(const edm::Run& Run, const edm::EventSetup& iSetup);
        virtual void beginJob() ;
        virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
        virtual void endRun(const edm::Run& Run, const edm::EventSetup& iSetup);
        virtual void endJob() ;

        void sampleCSCSegments();
        void sampleRPCsimHits();
        bool filterCSCSegment(const CSCSegment& sampleCSCSegment);
        bool filterSegmentEta(const CSCSegment& sampleCSCSegment);
        double deltaRforSegment(const GlobalPoint& globalPosition1, const GlobalVector& globalVector1, const GlobalPoint& globalPosition2, const GlobalVector& globalVector2);
        void findPeakTime(const CSCSegment& sampleCSCSegment);
        void fillSample();
        void EfficiencybyRPCrecHit();
        void fillEfficiency();

        // ----------member data ---------------------------
        edm::InputTag CSCSegmentsTag_;
        //edm::InputTag CSCStripDigisTag_;
        //edm::InputTag CSCWireDigisTag_;
        edm::InputTag RPCsimHitTag_;
        //edm::InputTag RPCDigiTag_;
        edm::InputTag RPCrecHitTag_;
        edm::InputTag simTrackTag_;

        edm::ESHandle<RPCGeometry> theRPCGeometry;
        edm::ESHandle<CSCGeometry> theCSCGeometry;

        edm::Handle<CSCSegmentCollection> pCSCSegments;
        edm::Handle<CSCStripDigiCollection> pCSCStripDigiCollection;
        edm::Handle<CSCWireDigiCollection> pCSCWireDigiCollection;
        edm::Handle<PSimHitContainer> pRPCsimHits;
        edm::Handle<RPCDigiCollection> pRPCDigis;
        edm::Handle<RPCRecHitCollection> pRPCrecHits;
        edm::Handle<SimTrackContainer> psimTracks;

        bool isSegmentMatchFilter;
        unsigned int SampleType;
        double deltaRTH;
        double ConeAngleX;
        double ConeAngleY;
        bool isEtaFilter;
        double MinSegmentEta;
        double MaxSegmentEta;
        int TrackType;

        double RangeStrips;
        double MaxD;
        double RollSize;
        double StripLength;
        double StripWidth;
        double ConeRadiusX;
        double ConeRadiusY;

        double CSCPeakTime2nd;
        double CSCPeakTimeAverage;
        int CSCWireTimeBin;  
        int RPCBX;
        double deltaBXTime;

        RPCDetId impactRPCRollId;
        LocalVector impactDirection;
        LocalPoint impactPoint;

        bool isvalideRPCBX;
        bool anycoincidence;

        bool issampled;
        bool isSamplefilled;
        bool ischecked;
        bool isEfficiencyfilled;

        std::map< CSCDetId, std::vector<RPCDetId> > CSC2RPCRolls;
        CSC2RPCBinder CSC2RPCMap;

        // Ntuple method
        int runNumber;
        unsigned int impactEventId;
        unsigned int impactRPCId;
        double impactCSCPeakTime2nd;
        double impactCSCPeakTimeAverage;
        double impactAngle;
        double impactDistanceZ;
        double impactStrip;
        unsigned int isEfficiency;
        double impactResidual;
        int impactClusterSize;
        int impactRPCBX;
        double impactGlobalPosition_X;
        double impactGlobalPosition_Y;
        double impactLocalPosition_X;
        double impactLocalPosition_Y;
        double CSCGlobalPosition_X;
        double CSCGlobalPosition_Y;
        double CSCGlobalPosition_Z;
        double RPCGlobalPosition_X;
        double RPCGlobalPosition_Y;
        double RPCGlobalPosition_Z;
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
    //CSCStripDigisTag_ = iConfig.getParameter<edm::InputTag>("CSCStripDigiTag");
    //CSCWireDigisTag_ = iConfig.getParameter<edm::InputTag>("CSCWireDigiTag");
    //RPCDigiTag_ = iConfig.getParameter<edm::InputTag>("RPCDigiTag");
    RPCrecHitTag_ = iConfig.getParameter<edm::InputTag>("RPCrecHitTag");
    RPCsimHitTag_ = iConfig.getParameter<edm::InputTag>("RPCsimHitTag");
    simTrackTag_ = iConfig.getParameter<edm::InputTag>("simTrackTag");

    isSegmentMatchFilter = iConfig.getParameter<bool>("isSegmentMatchFilter");
    SampleType = iConfig.getParameter<unsigned int>("SampleType");
    TrackType = iConfig.getParameter<int>("TrackType");
    deltaRTH = iConfig.getParameter<double>("deltaRTH");
    ConeAngleX = iConfig.getParameter<double>("ConeAngleX");
    ConeAngleY = iConfig.getParameter<double>("ConeAngleY");
    isEtaFilter = iConfig.getParameter<bool>("isEtaFilter");
    MinSegmentEta = iConfig.getParameter<double>("MinSegmentEta");
    MaxSegmentEta = iConfig.getParameter<double>("MaxSegmentEta");
    RangeStrips = iConfig.getParameter<double>("RangeStrips");
    MaxD = iConfig.getParameter<double>("MaxD");

    theRootFileName = iConfig.getUntrackedParameter<std::string>("EffRootFileName", "RPCEfficiencyTree.root");

    theFile = new TFile(theRootFileName.c_str(), "recreate");
    theFile->cd();
    ExTree= new TTree("ExTree", "ExTree");
    ExTree->Branch("runNumber", &runNumber, "runNumber/I");
    ExTree->Branch("impactEventId", &impactEventId, "impactEventId/i");
    ExTree->Branch("impactRPC", &impactRPCId, "impactRPCId/i");
    ExTree->Branch("impactCSCPeakTime2nd", &impactCSCPeakTime2nd, "impactCSCPeakTime2nd/D");
    ExTree->Branch("impactCSCPeakTimeAverage", &impactCSCPeakTimeAverage, "impactCSCPeakTimeAverage/D");
    ExTree->Branch("impactAngle", &impactAngle, "impactAngle/D");
    ExTree->Branch("impactStrip", &impactStrip, "impactStrip/D");
    ExTree->Branch("isEfficiency", &isEfficiency, "isEfficiency/i");
    ExTree->Branch("impactResidual", &impactResidual, "impactResidual/D");
    ExTree->Branch("impactClusterSize", &impactClusterSize, "impactClusterSize/I");
    ExTree->Branch("impactRPCBX", &impactRPCBX, "impactRPCBX/I");
    ExTree->Branch("impactDistanceZ", &impactDistanceZ, "impactDistanceZ/D");
    ExTree->Branch("impactConeRadiusX", &ConeRadiusX, "impactConeRadiusX/D");
    ExTree->Branch("impactConeRadiusY", &ConeRadiusY, "impactConeRadiusY/D");
    ExTree->Branch("impactGlobalPosition_X", &impactGlobalPosition_X, "impactGlobalPosition_X/D");
    ExTree->Branch("impactGlobalPosition_Y", &impactGlobalPosition_Y, "impactGlobalPosition_Y/D");
    ExTree->Branch("impactLocalPosition_X", &impactLocalPosition_X, "impactLocalPosition_X/D");
    ExTree->Branch("impactLocalPosition_Y", &impactLocalPosition_Y, "impactLocalPosition_Y/D");
    ExTree->Branch("CSCGlobalPosition_X", &CSCGlobalPosition_X, "CSCGlobalPosition_X/D");
    ExTree->Branch("CSCGlobalPosition_Y", &CSCGlobalPosition_Y, "CSCGlobalPosition_Y/D");
    ExTree->Branch("CSCGlobalPosition_Z", &CSCGlobalPosition_Z, "CSCGlobalPosition_Z/D");
    ExTree->Branch("RPCGlobalPosition_X", &RPCGlobalPosition_X, "RPCGlobalPosition_X/D");
    ExTree->Branch("RPCGlobalPosition_Y", &RPCGlobalPosition_Y, "RPCGlobalPosition_Y/D");
    ExTree->Branch("RPCGlobalPosition_Z", &RPCGlobalPosition_Z, "RPCGlobalPosition_Z/D");
}


CSCSegment2RPC::~CSCSegment2RPC() {

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
}


// ------------ method called to for each event  ------------
void CSCSegment2RPC::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    using namespace edm;
    using namespace std;

    // Get the geometry
    iSetup.get<MuonGeometryRecord>().get(theRPCGeometry);
    //theRPCGeometry = (RPCGeometry*)&*pRPCGeom;
    iSetup.get<MuonGeometryRecord>().get(theCSCGeometry);
    //theCSCGeometry = (CSCGeometry*)&*pCSCGeom;

    // Get CSC segment
    iEvent.getByLabel(CSCSegmentsTag_, pCSCSegments);
    // Get CSC strip digi collection
    //iEvent.getByLabel(CSCStripDigisTag_, pCSCStripDigiCollection);
    //Get CSC wire digi collection
    //iEvent.getByLabel(CSCWireDigisTag_, pCSCWireDigiCollection);
    // Get RPC Digis
    //iEvent.getByLabel(RPCDigiTag_, pRPCDigis);
    // Get RPC recHit
    iEvent.getByLabel(RPCrecHitTag_, pRPCrecHits);
    if(SampleType == 2) {
        // Get RPC simHits
        iEvent.getByLabel(RPCsimHitTag_, pRPCsimHits);
        // Get simTrack
        iEvent.getByLabel(simTrackTag_, psimTracks);
    }

    runNumber = iEvent.run();
    impactEventId = iEvent.id().event();
 
    if(SampleType == 1)
        sampleCSCSegments();

    if(SampleType == 2)
        sampleRPCsimHits();
}

void CSCSegment2RPC::sampleCSCSegments() {
    unsigned int nSegment = pCSCSegments->size();
    if(debug) cout << "This event has " << nSegment << " CSCSegments" << endl;
    //Statistics->Fill(nSegment);

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
        if(debug) cout << "CSC endcap is " << Endcap << ", Station is " << Station << ", Ring is " << Ring << " Chamber is " << Chamber << ", Layer is " << Layer << endl;
        if(Station == 4)
            continue;
        if(Ring == 1 || Ring == 4)
            continue;
        
        if(isEtaFilter == true && !filterSegmentEta(*CSCSegIter))
            continue;
        // One chamber has only 1 segment to avoid electron shower segments, and require at least 2 segments for this event
        if(debug) cout << "CSCSegmentsCounter: " << CSCSegmentsCounter[CSCId] << endl;

        if(CSCSegmentsCounter[CSCId] == 1) {
        //if(CSCSegmentsCounter[CSCId] == 1 && pCSCSegments->size() >= 2) {
            LocalPoint SegmentLocalPosition= CSCSegIter->localPosition();
            LocalVector SegmentLocalDirection=CSCSegIter->localDirection();
            if(debug) cout << "CSC segment local position: " << SegmentLocalPosition << endl;
            if(debug) cout << "CSC segment local direction: " << SegmentLocalDirection << endl;
            const GeomDet* CSCDet = theCSCGeometry->idToDet(CSCId);
            if(dynamic_cast<const CSCChamber*>(CSCDet) != 0) {
                const CSCChamber* CSCCh = dynamic_cast<const CSCChamber*>(CSCDet);
                GlobalPoint SegmentGlobalPosition = CSCCh->toGlobal(SegmentLocalPosition);
                GlobalVector SegmentGlobalDirection = CSCCh->toGlobal(SegmentLocalDirection);
                if(debug) cout << "CSC segment global position: " << SegmentGlobalPosition << endl;
                if(debug) cout << "CSC segment global direction: " << SegmentGlobalDirection << endl;
                CSCGlobalPosition_X = SegmentGlobalPosition.x();
                CSCGlobalPosition_Y = SegmentGlobalPosition.y();
                CSCGlobalPosition_Z = SegmentGlobalPosition.z();
                // Find correspond RPC rolls
                if(CSC2RPCMap.nearbyRPCRolls(CSCId).size() != 0) {
                    std::vector<RPCDetId> RPCRolls = CSC2RPCMap.nearbyRPCRolls(CSCId);
                    for(std::vector<RPCDetId>::const_iterator RPCRollIter = RPCRolls.begin(); RPCRollIter != RPCRolls.end(); RPCRollIter++) {
                        if(ischecked == true)
                            continue;
                        
                        ConeRadiusX = 0.;
                        ConeRadiusY = 0.;

                        int RPCStation = RPCRollIter->station();
                        int RPCRing = RPCRollIter->ring();
                        int RPCRollNumber = RPCRollIter->roll();
                        bool isSpecialRoll = false;
                        if((RPCStation == 1 && RPCRing == 3 && (RPCRollNumber == 1 || RPCRollNumber == 2)) || (RPCStation == 2 && RPCRing == 3 && RPCRollNumber == 1)) {
                            isSpecialRoll = true;
                            if(debug) cout << "SpecialRoll: " << RPCStation << ", " << RPCRing << ", " << RPCRollNumber << endl;
                        }
                        //isSpecialRoll = false;
                        bool passCSCFilter = filterCSCSegment(*CSCSegIter);
                        if(isSegmentMatchFilter == true && !isSpecialRoll && !passCSCFilter) {
                            if(debug) cout << "skip the filter for special roll: " << isSegmentMatchFilter << ", " << isSpecialRoll << ", " << passCSCFilter << endl;
                            continue;
                        }

                        const GeomDetUnit* RPCDet = theRPCGeometry->idToDetUnit(*RPCRollIter);
                        if(dynamic_cast<const RPCRoll*>(RPCDet) != 0) {
                            const RPCRoll* sampleRPCRoll = dynamic_cast<const RPCRoll*>(RPCDet);
                            const BoundPlane& sampleRPCRollSurface = sampleRPCRoll->surface();
                            const StripTopology& sampleRPCRollTop = sampleRPCRoll->specificTopology();

                            GlobalPoint test(1, 1, 1);
                            LocalPoint test1 = sampleRPCRoll->toLocal(test);
                            LocalPoint test2 = sampleRPCRollSurface.toLocal(test);
                            if(debug) cout << "For test: RPC roll local position is :" << test1 << endl;
                            if(debug) cout << "For test: RPC roll surface local position is :" << test2 << endl;

                            LocalPoint SegmentRPCPosition = sampleRPCRoll->toLocal(SegmentGlobalPosition);
                            LocalVector SegmentRPCDirection = sampleRPCRoll->toLocal(SegmentGlobalDirection);
                            double X0 = SegmentRPCPosition.x();
                            double Y0 = SegmentRPCPosition.y();
                            double Z0 = SegmentRPCPosition.z();
                            double dx = SegmentRPCDirection.x();
                            double dy = SegmentRPCDirection.y();
                            double dz = SegmentRPCDirection.z();
                            //for(unsigned int iStrip= 1; iStrip <= RPCRo->nstrips(); iStrip++)
                            //if(debug) cout << "Strip " << iStrip << " local position is " << RPCRo->centreOfStrip((int)iStrip) << endl;
                            double X = X0 + (-Z0) * dx / dz;
                            double Y = Y0 + (-Z0) * dy / dz;
                            LocalPoint impactRPCPosition(X, Y, 0);
                            double S = (SegmentRPCPosition - impactRPCPosition).mag();
                            //ConeRadiusX = S * ConeAngleX / cos(((SegmentRPCDirection.theta() > PI/2)?(PI-SegmentRPCDirection.theta()):(SegmentRPCDirection.theta())));
                            //ConeRadiusY = S * ConeAngleY / cos(((SegmentRPCDirection.theta() > PI/2)?(PI-SegmentRPCDirection.theta()):(SegmentRPCDirection.theta())));
                            ConeRadiusX = S * ConeAngleX;
                            ConeRadiusY = S * ConeAngleY;
                            // For RE1 the special layerout doesn't suit the cone algorithm
                            if(Station == 1)
                                ConeRadiusX = 0;

                            if(fabs(S) > MaxD) {
                                if(debug) cout << "Extropolate distance is too long!" << endl;
                                continue;
                            }
                            
                            StripLength = sampleRPCRollTop.stripLength();
                            StripWidth = sampleRPCRollTop.localPitch(impactRPCPosition);
                            int NStrips = sampleRPCRollTop.nstrips();
                            RollSize = (double)NStrips * StripWidth;

                            if(fabs(impactRPCPosition.z()) <= 1 && fabs(impactRPCPosition.x()) <= (RollSize*0.5+ConeRadiusX) && fabs(impactRPCPosition.y()) <= (StripLength*0.5+ConeRadiusY)) {

                                if(debug) cout << "CSC chamber: " << CSCId.rawId() << ", Endcap: " << Endcap << ", Station: " << Station << ", Ring: " << Ring << ", Chamber: " << Chamber << ", layer: " << Layer << ", RPCRing: " << RPCRing << ", " << "RPCRollNumber: " << RPCRollNumber << ", roll size: " << RollSize << ", strip length: " << StripLength << ", strip width: " << StripWidth << ", ConeRadiusX: " << ConeRadiusX <<  ", ConeRadiusY: " << ConeRadiusY << endl;
                                
                                impactRPCRollId = *RPCRollIter; 
                                impactDirection = SegmentRPCDirection; 
                                impactPoint = impactRPCPosition; 
                                findPeakTime(*CSCSegIter);

                                // For Ntuple
                                impactDistanceZ = Z0;
                                impactRPCId = RPCRollIter->rawId();
                                impactCSCPeakTime2nd = CSCPeakTime2nd;
                                impactCSCPeakTimeAverage = CSCPeakTimeAverage;
                                impactAngle = ((SegmentRPCDirection.theta() > PI/2) ? (PI- SegmentRPCDirection.theta()) : SegmentRPCDirection.theta());
                                impactStrip = sampleRPCRollTop.strip(impactRPCPosition); // If outside the strip topology, the StripNumber will be set to 0/N_strip
                                impactLocalPosition_X = impactRPCPosition.x();
                                impactLocalPosition_Y = impactRPCPosition.y();
                                GlobalPoint impactPoint_Global = sampleRPCRoll->toGlobal(impactRPCPosition);
                                impactGlobalPosition_X = impactPoint_Global.x();
                                impactGlobalPosition_Y = impactPoint_Global.y();

                                RPCGlobalPosition_X = impactPoint_Global.x();
                                RPCGlobalPosition_Y = impactPoint_Global.y();
                                RPCGlobalPosition_Z = impactPoint_Global.z();

                                if(debug) cout << "impactRPCId is " << impactRPCId << ", impactCSCPeakTime2nd is " << impactCSCPeakTime2nd << ", impactAngle is " << impactAngle << ", impactStrip is " << impactStrip << ", impactLocalPosition_X is " << impactLocalPosition_X << ", impactLocalPosition_Y is " << impactLocalPosition_Y << ", impactGlobalPosition_X is " << impactGlobalPosition_X << ", impactGlobalPosition_Y is " << impactGlobalPosition_Y << endl;
                                // check efficiency
                                issampled = true;
                                if(ischecked == false)
                                    EfficiencybyRPCrecHit();
                                if(anycoincidence == true && ischecked == false) {
                                    ischecked = true;
                                    fillEfficiency();
                                }
                                // avoid efficiency is calculated for 2 times
                            }
                            else {
                                if(debug) cout << "out of roll size" << endl;
                            }
                        }
                    }
                    if(issampled == true && ischecked == false) {
                        if(debug) cout << "Could not find a Roll has recHit inside window." << endl;
                        ischecked = true;
                        //fillSample();
                        fillEfficiency();
                    }
                }
            }
        }
    }
}

void CSCSegment2RPC::findPeakTime(const CSCSegment& sampleCSCSegment) {

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

        /*
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
        if(debug) cout << "Find " << sampleCSCWireCluster.size() << " wire groups in the cluster" << endl;
        for(CSCWireDigiCollection::const_iterator sampleCSCWireDigiIter = sampleCSCWireDigiRange.first; sampleCSCWireDigiIter != sampleCSCWireDigiRange.second; sampleCSCWireDigiIter++) {
            int sampleWireGroup = sampleCSCWireDigiIter->getWireGroup();
            const std::vector<int>::const_iterator resultWireIter = find(sampleCSCWireCluster.begin(), sampleCSCWireCluster.end(), sampleWireGroup);
            if(resultWireIter != sampleCSCWireCluster.end()) {
                std::vector<int> sampleWireTimeBins = sampleCSCWireDigiIter->getTimeBinsOn();
                for(std::vector<int>::iterator Iter = sampleWireTimeBins.begin(); Iter != sampleWireTimeBins.end(); Iter++)
                    if(debug) cout << "Sample wire BX time bin word is: " << (*Iter) << endl;
                // The 1st time bin in the wire
                int sampleWireTimeBin = (*(sampleWireTimeBins.begin()));
                if(sampleWireGroup < firstWireGroup) {
                    firstWireGroup = sampleWireGroup;
                    firstWireTimeBin = sampleWireTimeBin;
                }
                if(debug) cout << "Sample wire BX time bin is: " << firstWireTimeBin << endl;
            }
        }
        WireTimeBinSum += firstWireTimeBin;
        */
    }

    CSCPeakTime2nd = peakTime2nd;
    CSCPeakTimeAverage = peakTimesum / sampleCSCrecHitDGroup.size();
    if(debug) cout << "CSCPeakTime2nd is " << CSCPeakTime2nd << ", CSCPeakTimeAverage is " << CSCPeakTimeAverage << endl;

    //averageWireTimeBin = WireTimeBinSum / sampleCSCrecHitDGroup.size();
    //CSCWireTimeBin = averageWireTimeBin - 6;
    //if(debug) cout << "CSC average wire time bin is " << averageWireTimeBin << ", CSC wire BX is " << CSCWireTimeBin << endl;

}

void CSCSegment2RPC::sampleRPCsimHits() {

    for(PSimHitContainer::const_iterator simHitIter = pRPCsimHits->begin(); simHitIter != pRPCsimHits->end(); simHitIter++) {
        // The particle type of the hit may differ from the particle type of the SimTrack with id trackId(). 
        // This happends if the hit was created by a secondary track (e.g. a delta ray) originating from the trackId() and not existing as a separate SimTrack.
        int Particletype = simHitIter->particleType();
        if(abs(Particletype) != TrackType)
            continue;

        // Check the DetId is a RPCDetId
        unsigned int sampleDetUnitId = simHitIter->detUnitId();
        DetId sampleDetectorId = DetId(sampleDetUnitId);
        GeomDetEnumerators::SubDetector subDet = theRPCGeometry->idToDetUnit(sampleDetectorId)->subDetector();
        if(subDet != GeomDetEnumerators::RPCEndcap)        
            continue;

        // Check if inside the collection range
        const RPCDetId sampleRPCRollId(sampleDetectorId);
        //if(meCollection.find(sampleRPCRollId) == meCollection.end())
            //continue;

        issampled = false;
        isSamplefilled = false;
        ischecked = false;
        isEfficiencyfilled = false;

        const GeomDetUnit* RPCDet = theRPCGeometry->idToDetUnit(sampleRPCRollId);
        if(dynamic_cast<const RPCRoll*>(RPCDet) != 0) {
            const RPCRoll* sampleRPCRoll = dynamic_cast<const RPCRoll*>(RPCDet);
            //const BoundPlane& sampleRPCRollSurface = sampleRPCRoll->surface();
            const StripTopology& sampleRPCRollTop = sampleRPCRoll->specificTopology();
            Local3DPoint sampleLocalPosition = simHitIter->localPosition();
            LocalPoint impactRPCPosition = sampleLocalPosition;
            LocalVector SegmentRPCDirection = simHitIter->momentumAtEntry(); 
            StripLength = sampleRPCRollTop.stripLength();
            StripWidth = sampleRPCRollTop.localPitch(impactRPCPosition);
            int NStrips = sampleRPCRollTop.nstrips();
            RollSize = (double)NStrips * StripWidth;
            ConeRadiusX = 0.;
            ConeRadiusY = 0.;
            if(fabs(impactRPCPosition.z()) <= 1 && fabs(impactRPCPosition.x()) <= RollSize*0.5 && fabs(impactRPCPosition.y()) <= StripLength*0.5) {
                impactRPCRollId = sampleRPCRollId;
                impactDirection = SegmentRPCDirection;
                impactPoint = impactRPCPosition;
                CSCWireTimeBin = 0;
                // For Ntuple 
                impactRPCId = sampleRPCRollId.rawId(); 
                impactCSCPeakTime2nd = 0.;
                impactCSCPeakTimeAverage = 0.;
                impactAngle = ((SegmentRPCDirection.theta() > PI/2) ? (PI- SegmentRPCDirection.theta()) : SegmentRPCDirection.theta());
                impactStrip = sampleRPCRollTop.strip(impactRPCPosition);
                impactLocalPosition_X = impactRPCPosition.x();
                impactLocalPosition_Y = impactRPCPosition.y();
                GlobalPoint impactPoint_Global = sampleRPCRoll->toGlobal(impactRPCPosition);
                impactGlobalPosition_X = impactPoint_Global.x();
                impactGlobalPosition_Y = impactPoint_Global.y();
                                                                                

                // check efficiency
                issampled = true;
                //fillSample();
                EfficiencybyRPCrecHit();
                ischecked = true;
                fillEfficiency();
            }
        }
    }
}


bool CSCSegment2RPC::filterCSCSegment(const CSCSegment& sampleCSCSegment) {

    bool isTrackSegment = false;
    CSCDetId sampleCSCId = sampleCSCSegment.cscDetId();
    const CSCChamber* sampleCSCChamber = theCSCGeometry->chamber(sampleCSCId);
    LocalVector localsampleVector = sampleCSCSegment.localDirection();
    LocalPoint localsamplePosition = sampleCSCSegment.localPosition();
    GlobalVector globalsampleVector = sampleCSCChamber->toGlobal(localsampleVector);
    GlobalPoint globalsamplePosition = sampleCSCChamber->toGlobal(localsamplePosition);
    int sampleCSCStation = sampleCSCId.station();
    int sampleCSCzEndcap = sampleCSCId.zendcap();

    for(CSCSegmentCollection::const_iterator CSCSegIter = pCSCSegments->begin(); CSCSegIter != pCSCSegments->end(); CSCSegIter++) {
        CSCDetId iterCSCId = CSCSegIter->cscDetId();
        int iterCSCStation = iterCSCId.station();
        int iterCSCzEndcap = iterCSCId.zendcap();
        if((abs(sampleCSCStation-iterCSCStation) == 0) || (sampleCSCzEndcap != iterCSCzEndcap))
            continue;

        if(debug) cout << "sampleCSCSegment:" << sampleCSCStation << ", " << sampleCSCzEndcap << ". CSCSegIter:" << iterCSCStation << ", " << iterCSCzEndcap << endl;

        const CSCChamber* iterCSCChamber = theCSCGeometry->chamber(iterCSCId);
        LocalVector localIterVector = CSCSegIter->localDirection();
        LocalPoint localIterPosition = CSCSegIter->localPosition();
        GlobalVector globalIterVector = iterCSCChamber->toGlobal(localIterVector);
        GlobalPoint globalIterPosition = iterCSCChamber->toGlobal(localIterPosition);

        double deltaR1 = deltaRforSegment(globalsamplePosition, globalsampleVector, globalIterPosition, globalIterVector);
        double deltaR2 = deltaRforSegment(globalIterPosition, globalIterVector, globalsamplePosition, globalsampleVector);
        if(deltaR1 <= deltaRTH || deltaR2 <= deltaRTH)
            isTrackSegment = true;
    }
    if(isTrackSegment == true) {
        if(debug) cout << "Pass Segment filter. successfully." << endl;
    }
    else {
        if(debug) cout << "Pass Segment filter unsuccessfully." << endl;
    }

    return isTrackSegment;
}

double CSCSegment2RPC::deltaRforSegment(const GlobalPoint& globalPosition1, const GlobalVector& globalVector1, const GlobalPoint& globalPosition2, const GlobalVector& globalVector2) {
    
    GlobalVector VectorPoint2Point = (GlobalVector)(globalPosition2 - globalPosition1);
    double deltaR1 = sqrt((globalVector2.phi() - globalVector1.phi()).value() * (globalVector2.phi() - globalVector1.phi()).value() + (globalVector2.theta() - globalVector1.theta()) * (globalVector2.theta() - globalVector1.theta()));
    double deltaR2 = sqrt((VectorPoint2Point.phi() - globalVector1.phi()).value() * (VectorPoint2Point.phi() - globalVector1.phi()).value() + (VectorPoint2Point.theta() - globalVector1.theta()) * (VectorPoint2Point.theta() - globalVector1.theta()));
    if(debug) cout << "deltaR1: " << deltaR1 << ", deltaR2: " << deltaR2 << ". globalVector1: " << globalVector1.phi().value() << ", " << globalVector1.theta() << ". globalVector2: " << globalVector2.phi().value() << ", " << globalVector2.theta() << endl;
    return (deltaR1 > deltaR2) ? deltaR1 : deltaR2;
}

bool CSCSegment2RPC::filterSegmentEta(const CSCSegment& sampleCSCSegment) {

    bool isinsideEtaRange = false;
    CSCDetId sampleCSCId = sampleCSCSegment.cscDetId();
    const CSCChamber* sampleCSCChamber = theCSCGeometry->chamber(sampleCSCId);
    LocalVector localsampleVector = sampleCSCSegment.localDirection();
    GlobalVector globalsampleVector = sampleCSCChamber->toGlobal(localsampleVector);
    double sampleSegmentEta = globalsampleVector.eta();
    if((fabs(sampleSegmentEta) >= MinSegmentEta) && (fabs(sampleSegmentEta) <= MaxSegmentEta))
        isinsideEtaRange = true;

    if(isinsideEtaRange == true) {
        if(debug) cout << "Pass Eta filter. successfully." << endl;
    }
    else {
        if(debug) cout << "Pass Eta filter unsuccessfully." << endl;
    }

    return isinsideEtaRange;
}

void CSCSegment2RPC::fillSample() {

    if(issampled == false || isSamplefilled == true) {
        if(debug) cout << "Not yet sample the RPC rolls or sample already filled" << endl;
        return;
    }

    isSamplefilled = true;
}


// Try to find the fire strip from impact RPC roll
void CSCSegment2RPC::EfficiencybyRPCrecHit() {

    if(ischecked == true) {
        if(debug) cout << "Already checked and find efficiency" << endl;
        return;        
    }


    double residualDistance = 10000;
    LocalPoint FirePosition;
    int FireClusterSize = 0;
    if(debug) cout << "Try to find RPC strips closest to impact point" << endl;

    anycoincidence = false;
    isvalideRPCBX = false;
    isEfficiency = 0;

    RPCRecHitCollection::range RollrecHitsRange = pRPCrecHits->get(impactRPCRollId);
    RPCBX = 0;
    for(RPCRecHitCollection::const_iterator recHitIter = RollrecHitsRange.first; recHitIter != RollrecHitsRange.second; recHitIter++) {
        //const RPCRecHit& = *recHitIter;
        LocalPoint candidatePosition = recHitIter->localPosition();
        double residualDistance_Temp = candidatePosition.x()-impactPoint.x();
        int candiadteClusterSize = recHitIter->clusterSize();
        int candidateBX = recHitIter->BunchX();
        if(debug) cout << "ResidualDistance_Temp is " << residualDistance_Temp << endl;
        if(fabs(residualDistance_Temp) < fabs(residualDistance)) {
            residualDistance = residualDistance_Temp;
            FirePosition = candidatePosition;
            FireClusterSize = candiadteClusterSize;
            RPCBX = candidateBX;
            deltaBXTime = CSCWireTimeBin - RPCBX;
            if(debug) cout << "Step into a candidate strip. residualDistance: " << residualDistance << ", FireClusterSize: " << FireClusterSize << endl;
        }
    }
    if(fabs(residualDistance) <= ((RangeStrips+FireClusterSize*0.5)*StripWidth)) {
        if(debug) cout << "Find a fire strip with-in the residual range in RPC roll " << impactRPCRollId.rawId() << endl;
        anycoincidence = true;
        isvalideRPCBX = true;
        // For tree
        isEfficiency = 1;
        impactResidual = residualDistance;
        impactClusterSize = FireClusterSize;
        impactRPCBX = RPCBX;
    }
    else {
        if(debug) cout <<"Could not find a fire strip with-in the residual range in RPC roll " << impactRPCRollId.rawId() << endl;
        anycoincidence = false;
        isvalideRPCBX = false;
        // For tree
        isEfficiency = 0;
        impactResidual = residualDistance;
        impactClusterSize = FireClusterSize;
        impactRPCBX = RPCBX;
    }

    // prevent multi checkking of efficiency from one CSC segment
    //if(anycoincidence == true)
        //ischecked = true;
}

void CSCSegment2RPC::fillEfficiency() {

    if(ischecked == false || isEfficiencyfilled == true) {
        if(debug) cout << "Not yet check the RPC roll or efficiency already filled." << endl;
        return;
    }
    // Fill tree whatever efficiency is 0 or 1;
    ExTree->Fill();

    isEfficiencyfilled = true;
}


// ------------ method called once each job just before starting event loop  ------------
void CSCSegment2RPC::beginJob() {
}

void CSCSegment2RPC::beginRun(const edm::Run& Run, const edm::EventSetup& iSetup) {

    edm::ESHandle<RPCGeometry> RPCGeometryEH;
    iSetup.get<MuonGeometryRecord>().get(RPCGeometryEH);
    edm::ESHandle<CSCGeometry> CSCGeometryEH;
    iSetup.get<MuonGeometryRecord>().get(CSCGeometryEH);
    CSC2RPCMap.initiate(RPCGeometryEH, CSCGeometryEH);
    CSC2RPCMap.checkBinding(RPCGeometryEH, CSCGeometryEH);

    const std::vector<RPCRoll*> rpcRolls= RPCGeometryEH->rolls();
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
        if(CSC2RPCRolls.find(CSCIndex) != CSC2RPCRolls.end())
            RPCRolls_temp = CSC2RPCRolls[CSCIndex];
        RPCRolls_temp.push_back(rpcRoll);
        CSC2RPCRolls[CSCIndex] = RPCRolls_temp;
        //meCollection[rpcRoll] = bookDetUnitSeg(rpcRoll, (*RPCRollIter)->nstrips(), (*RPCRollIter)->pitch(), (*RPCRollIter)->specificTopology().stripLength());
    }
}

// ------------ method called once each job just after ending the event loop  ------------
void CSCSegment2RPC::endJob() {
    //DBE = 0;
        if(debug) cout << "writing information into file: " << theFile->GetName() << std::endl;
        theFile->Write();
        theFile->Close();
}

void CSCSegment2RPC::endRun(const edm::Run& Run, const edm::EventSetup& iSetup) {
    //if(EffSaveRootFile)
        //DBE->save(EffRootFileName);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CSCSegment2RPC);
