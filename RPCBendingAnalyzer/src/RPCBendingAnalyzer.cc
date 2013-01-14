// -*- C++ -*-
//
// Package:    RPCBendingAnalyzer
// Class:      RPCBendingAnalyzer
// 
/**\class RPCBendingAnalyzer RPCBendingAnalyzer.cc MyModule/RPCBendingAnalyzer/src/RPCBendingAnalyzer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
 */
//
// Original Author:  Haiyun Teng,591 R-005,+41227671371,
//         Created:  Thu Aug 26 02:12:18 CEST 2010
// $Id: RPCBendingAnalyzer.cc,v 1.7 2012/06/06 19:35:02 hyteng Exp $
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

#include <FWCore/Framework/interface/ESHandle.h>
#include <SimDataFormats/Track/interface/SimTrack.h>
#include <SimDataFormats/TrackingHit/interface/PSimHit.h>
#include <SimDataFormats/RPCDigiSimLink/interface/RPCDigiSimLink.h>
#include <DataFormats/RPCRecHit/interface/RPCRecHit.h>
#include <DataFormats/DetId/interface/DetId.h>
#include <DataFormats/MuonDetId/interface/RPCDetId.h>
#include <Geometry/CommonDetUnit/interface/GeomDetUnit.h>
#include <Geometry/RPCGeometry/interface/RPCGeometry.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <Geometry/RPCGeometry/interface/RPCGeomServ.h>
#include <DataFormats/RPCRecHit/interface/RPCRecHitCollection.h>
#include <SimDataFormats/TrackingHit/interface/PSimHitContainer.h>
#include <SimDataFormats/Track/interface/SimTrackContainer.h>
#include <iostream>
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSet.h"
#include <DataFormats/GeometryVector/interface/Phi.h>
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TTree.h"
#include "TFile.h"

#define RPCBarrelLayerNumber 6
#define RPCEndcapLayerNumber 5

// using name space
using namespace std;
using namespace edm;

//
// class declaration
//

class RPCBendingAnalyzer : public EDAnalyzer {
    public:
        explicit RPCBendingAnalyzer(const ParameterSet&);
        ~RPCBendingAnalyzer();


    private:
        virtual void beginJob() ;
        virtual void analyze(const Event&, const EventSetup&);
        virtual void endJob() ;

        void getTracksInfo();
        void findRPCrecHits(const PSimHit& RPCsimHit, vector<RPCRecHit>& theRPCrecHits);
        int getRPCLayer(const RPCDetId& theRPCDetId);
        void analyzeBending();
        void fillSample(unsigned int Index);
        void fillBendingInfo();

        // ----------member data ---------------------------
        InputTag simHitTag_;
        InputTag simTrackTag_;
        InputTag RPCrecHitTag_;
        InputTag RPCDigiSimLinkTag_;
        vector<unsigned int> RPCLayer;

        Handle<PSimHitContainer> simHitCollectionH;
        Handle<SimTrackContainer> simTrackCollectionH;
        Handle<RPCRecHitCollection> RPCrecHitCollectionH;
        Handle< DetSetVector<RPCDigiSimLink> > thelinkDigis;
        ESHandle<RPCGeometry> theRPCGeometry;

        bool debug;
        unsigned int codeTH;
        unsigned int recHitNumberTH;
        vector<PSimHit> simHitsforTrack;
        vector<unsigned int> recHitNumberforsimHits;
        vector< vector<RPCRecHit> > RPCrecHitsfromTrack;
        RPCRecHit recHitSample[10];
        PSimHit simHitSample[10];

        string theRootFileName;
        TFile *theFile;
        TTree *ExTree;

        unsigned int EventNumber;
        int simTrackId;
        int simTrackType;
        double simTrackMomentumPt;
        double simTrackMomentumEta;
        int simTrackCharge;
        int simTrackvalid;
        int SampleLayer[10];
        double simHitBendingPhi[10][10];
        double simPtBendingPhi[10];
        double recHitBendingPhi[10][10];
        unsigned int RPCBendingLayer[10];
        int ClusterSize[10];
        int BX[10];
        double Eta[10];
        double simPtatRef[10]; // Ref recHit is in layer1(RB1out) for Barrel
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
RPCBendingAnalyzer::RPCBendingAnalyzer(const ParameterSet& iConfig) {
    //now do what ever initialization is needed
    simHitTag_ = iConfig.getUntrackedParameter<InputTag>("simHitTag");
    simTrackTag_ = iConfig.getUntrackedParameter<InputTag>("simTrackTag");
    RPCrecHitTag_ = iConfig.getUntrackedParameter<InputTag>("RPCrecHitTag");
    RPCDigiSimLinkTag_ = iConfig.getUntrackedParameter<InputTag>("RPCDigiSimLinkTag");
    RPCLayer = iConfig.getParameter< vector<unsigned int> >("RPCLayer");
    debug = iConfig.getUntrackedParameter<bool>("debug");
    codeTH = iConfig.getUntrackedParameter<unsigned int>("codeTH");
    recHitNumberTH = iConfig.getUntrackedParameter<unsigned int>("recHitNumberTH");
    theRootFileName = iConfig.getUntrackedParameter<string>("theRootFileName", "RPCValidationTree.root");
    // Booking histogrqam
    theFile = new TFile(theRootFileName.c_str(), "recreate");
    theFile->cd();
    ExTree = new TTree("ExTree", "ExTree");
    ExTree->Branch("EventNumber", &EventNumber, "EventNumber/i");
    ExTree->Branch("simTrackId", &simTrackId, "simTrackId/I");
    ExTree->Branch("simTrackType", &simTrackType,  "simTrackType/I");
    ExTree->Branch("simTrackMomentumPt", &simTrackMomentumPt, "simTrackMomentumPt/D");
    ExTree->Branch("simTrackMomentumEta", &simTrackMomentumEta, "simTrackMomentumEta/D");
    ExTree->Branch("simTrackCharge", &simTrackCharge, "simTrackCharge/I");
    ExTree->Branch("simTrackvalid", &simTrackvalid, "simTrackvalid/I");
    ExTree->Branch("SampleLayer", SampleLayer, "SampleLayer[10]/I");
    ExTree->Branch("simHitBendingPhi", simHitBendingPhi, "simHitBendingPhi[10][10]/D");
    ExTree->Branch("simPtBendingPhi", simPtBendingPhi, "simPtBendingPhi[10]/D");
    ExTree->Branch("recHitBendingPhi", recHitBendingPhi, "recHitBendingPhi[10][10]/D");
    ExTree->Branch("RPCBendingLayer", RPCBendingLayer, "RPCBendingLayer[10]/i");
    ExTree->Branch("simPtatRef", simPtatRef, "simPtatRef[10]/D");
    ExTree->Branch("ClusterSize", ClusterSize, "ClusterSize[10]/I");
    ExTree->Branch("BX", BX, "BX[10]/I");
    ExTree->Branch("Eta", Eta, "Eta[10]/D");
    // Set simTrack index
    EventNumber = 0;
}


RPCBendingAnalyzer::~RPCBendingAnalyzer() {

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void RPCBendingAnalyzer::analyze(const Event& iEvent, const EventSetup& iSetup) {
    using namespace edm;

    iEvent.getByLabel(simHitTag_, simHitCollectionH);
    iEvent.getByLabel(simTrackTag_, simTrackCollectionH);
    iEvent.getByLabel(RPCrecHitTag_, RPCrecHitCollectionH);
    iEvent.getByLabel("simMuonRPCDigis", "RPCDigiSimLink", thelinkDigis);

    iSetup.get<MuonGeometryRecord>().get(theRPCGeometry);

    EventNumber++;
    if(debug) cout << "EventNumner: "<< EventNumber << endl;

    getTracksInfo();
}

void RPCBendingAnalyzer::getTracksInfo() {

    unsigned int simTrackIndex = 0;
    for(SimTrackContainer::const_iterator simTrackIter = simTrackCollectionH->begin(); simTrackIter != simTrackCollectionH->end(); simTrackIter++) {
        simTrackIndex++;
        // for test the bending of simHits
        simHitsforTrack.clear();
        recHitNumberforsimHits.clear();
        RPCrecHitsfromTrack.clear();

        int simTrackId = simTrackIter->trackId();
        int simTrackType = simTrackIter->type();
        if(debug) cout << "SimTrack id: " << simTrackId << endl;
        if(debug) cout << "SimTrack type: " << simTrackType << endl;
        //if(debug) cout << "SimTrack momentum: " << simTrackIter->momentum().Pt() << endl;

        if(abs(simTrackType) != 13)
            continue;

        simTrackCharge = simTrackIter->charge();
        simTrackMomentumPt = simTrackIter->momentum().pt();
        simTrackMomentumEta = simTrackIter->momentum().eta();
        if(debug) cout << "SimTrack momentum Pt: " << simTrackMomentumPt << ", Eta: " << simTrackMomentumEta << endl;

        // recHit number for RPCLayer, 0-5 for barrel, 6-10 for endcap both side
        int recHitNumberinRPCLayer[(RPCBarrelLayerNumber+RPCEndcapLayerNumber)];
        for(unsigned int i = 0; i < (RPCBarrelLayerNumber+RPCEndcapLayerNumber); i++)
            recHitNumberinRPCLayer[i] = 0;

        for(PSimHitContainer::const_iterator simHitIter = simHitCollectionH->begin(); simHitIter != simHitCollectionH->end(); simHitIter++) {
            int TrackId1 = simHitIter->trackId();
            if(TrackId1 != simTrackId)
                continue;

            // The particle type of the hit may differ from the particle type of the SimTrack with id trackId(). 
            // This happends if the hit was created by a secondary track (e.g. a delta ray) originating from the trackId() and not existing as a separate SimTrack.
            int ParticleType1 = simHitIter->particleType();
            if(ParticleType1 != simTrackType)
                continue;

            if(debug) cout << "simHit particle type : " << ParticleType1 << endl;
            int DetUnitId1 = simHitIter->detUnitId();
            DetId DetectorId1(DetUnitId1);
            GeomDetEnumerators::SubDetector subDet = theRPCGeometry->idToDetUnit(DetectorId1)->subDetector();
            // For barrel and endcap in case of double layer implenment in endcap
            if(subDet != GeomDetEnumerators::RPCBarrel && subDet != GeomDetEnumerators::RPCEndcap)
                continue;

            if(debug) cout << "simHit DetUnit id: " << DetUnitId1 << endl;
            Local3DPoint HitPoint = simHitIter->entryPoint();
            if(debug) cout << "simHit local hit point: " << HitPoint << endl;
            LocalVector Momentum = simHitIter->momentumAtEntry();
            if(debug) cout << "simHit local momentum: " << Momentum << endl;
            unsigned int processType = simHitIter->processType();
            if(debug) cout << "simHit process type: " << processType << endl;

            const GeomDetUnit* rpcRoll = theRPCGeometry->idToDetUnit(DetectorId1);
            GlobalPoint simHitPosition = rpcRoll->toGlobal(HitPoint);
            GlobalVector globalMomentum = rpcRoll->toGlobal(Momentum);
            if(debug) cout << "At this simHit the global Position is " << simHitPosition << endl;
            if(debug) cout << "At this simHit the global Pt value: " << globalMomentum.perp() << ", Phi: " << globalMomentum.phi() << endl;
            simHitsforTrack.push_back(*simHitIter);
            vector<RPCRecHit> tempRPCrecHits;
            tempRPCrecHits.clear();
            findRPCrecHits(*simHitIter, tempRPCrecHits);
            RPCrecHitsfromTrack.push_back(tempRPCrecHits);
            for(vector<RPCRecHit>::const_iterator recHitIter = tempRPCrecHits.begin(); recHitIter != tempRPCrecHits.end(); recHitIter++) {
                //RPCrecHitsfromTrack.push_back(*recHitIter);
                RPCDetId theRPCDetId = recHitIter->rpcId();
                int theRPCStationLayer = getRPCLayer(theRPCDetId);
                // If StationLayer<0 bending and codeTH will not use this recHit but still save into RPCrecHitsforTrack, here the getRPCLayer() will return 1-6 for barrel and 7-9 for endcap
                if(theRPCStationLayer >= 1)
                    recHitNumberinRPCLayer[theRPCStationLayer-1] += 1;
            }
            recHitNumberforsimHits.push_back(tempRPCrecHits.size());
        }

        // code for checking layer occupancy
        unsigned int code = 0;
        unsigned int recHitNumber = 0;
        for(unsigned int RPCLayerIndex = 0; RPCLayerIndex < (RPCBarrelLayerNumber+RPCEndcapLayerNumber); RPCLayerIndex++) {
            if(recHitNumberinRPCLayer[RPCLayerIndex] != 0) {
                unsigned value = 1;
                recHitNumber++;
                for(unsigned int i = 0; i < RPCLayerIndex; i++)
                    value *= 2;
                code += value;
            }
        }
        if(((int)code & (int)codeTH) == (int)codeTH && recHitNumber >= recHitNumberTH)
            simTrackvalid = 1;
        else
            simTrackvalid = 0;


        if(simTrackvalid == 1)
            analyzeBending();
    }
}

void RPCBendingAnalyzer::findRPCrecHits(const PSimHit& RPCsimHit, vector<RPCRecHit>& theRPCrecHits) {

    theRPCrecHits.clear();
    int TrackId1 = RPCsimHit.trackId();
    // The particle type of the hit may differ from the particle type of the SimTrack with id trackId(). 
    // This happends if the hit was created by a secondary track (e.g. a delta ray) originating from the trackId() and not existing as a separate SimTrack.
    int ParticleType1 = RPCsimHit.particleType();
    int DetUnitId1 = RPCsimHit.detUnitId();
    DetId DetectorId1(DetUnitId1);
    int firstStrip = -1;
    int lastStrip = -1;
    int DigiBX = 0;
    if(debug) cout << "Checking DigiSimLink..." << endl;
    for(DetSetVector<RPCDigiSimLink>::const_iterator itlink = thelinkDigis->begin(); itlink != thelinkDigis->end(); itlink++) {
        for(DetSet<RPCDigiSimLink>::const_iterator digi_iter=itlink->data.begin(); digi_iter != itlink->data.end(); ++digi_iter) {
            int DetUnitId2 = digi_iter->getDetUnitId();
            int TrackId2 = digi_iter->getTrackId();
            Local3DPoint HitPoint2 = digi_iter->getEntryPoint();
            LocalVector Momentum2 = digi_iter->getMomentumAtEntry();
            int ParticleType2 = digi_iter->getParticleType();
            if((DetUnitId1 == DetUnitId2) && (TrackId1 == TrackId2) && (ParticleType1 == ParticleType2)) {
                int Strip = digi_iter->getStrip();
                DigiBX = digi_iter->getBx();
                if(debug) cout << "Find Digi's Strip: " << Strip << endl;
                if(debug) cout << "Find Digi's BX: " << DigiBX << endl;
                if(debug) cout << "Find Digi's hit positin: " << HitPoint2 << endl;
                if(firstStrip == -1) {
                    firstStrip = Strip;
                    lastStrip = Strip;
                }
                if(Strip < firstStrip)
                    firstStrip = Strip;
                if(Strip > lastStrip)
                    lastStrip = Strip;
            }
        }
    }

    RPCRecHitCollection::range RPCrecHitRange = RPCrecHitCollectionH->get(DetectorId1);
    int associatedrecHitforsimHit = 0;
    for(RPCRecHitCollection::const_iterator recHitIter = RPCrecHitRange.first; recHitIter != RPCrecHitRange.second; recHitIter++) {        
        int firstRecStrip = recHitIter->firstClusterStrip();
        int ClusterSize = recHitIter->clusterSize();
        int recBX = recHitIter->BunchX();
        bool isValid = recHitIter->isValid();
        if(debug) cout << "RecHit 1st strip: " << firstRecStrip << endl;
        if(debug) cout << "RecHit ClusterSize: " << ClusterSize << endl;
        if(debug) cout << "RecHit BX: " << recBX << endl;
        if(debug) cout << "RecHit Valid: " << isValid << endl;
        if((firstRecStrip <= lastStrip)  && ((firstRecStrip + ClusterSize - 1) >= firstStrip) && isValid && recBX == DigiBX) {
            if(debug) cout << "Find SimHit-RecHit" << endl;
            RPCDetId theRPCDetId = recHitIter->rpcId();
            RPCGeomServ theRPCService = RPCGeomServ(theRPCDetId);
            string theRPCName = theRPCService.name();
            const GeomDetUnit *rpcRoll = theRPCGeometry->idToDetUnit(theRPCDetId);
            LocalPoint localPosition = recHitIter->localPosition();
            GlobalPoint globalPosition = rpcRoll->toGlobal(localPosition);
            if(debug) cout << "RPC det is: " << theRPCName << " GlobalPosition: " << globalPosition << endl;
            if(debug) cout << "Push this recHit to stack for validation." << endl;
            associatedrecHitforsimHit++;
            theRPCrecHits.push_back(*recHitIter);
        }
    }
    if(debug) cout << "Find " << associatedrecHitforsimHit << " recHits in the same DetUnit" << endl;
}

int RPCBendingAnalyzer::getRPCLayer(const RPCDetId& theRPCDetId) {

    RPCGeomServ theRPCService = RPCGeomServ(theRPCDetId);
    string theRPCName = theRPCService.name();
    int Region = theRPCDetId.region(); // Region=0 for barrel, +1 for RE+, -1 for RE-
    int Station = theRPCDetId.station();
    int Layer = theRPCDetId.layer();
    int RPCStationLayer = -1;
    if(Region == 0) {
        if(Station >= 3) // RB3,RB4
            RPCStationLayer = Station + 2;
        else {
            if(Station == 2) // RB2in/put
                RPCStationLayer = Layer + 2;
            else // RB1in/out
                RPCStationLayer = Layer;
        }
    }
    else { 
        // Before endcap double layer implement, RE1-4
        RPCStationLayer = Station + 6;
        // After endcap double layer implement, RE2in/out may included
    }
    if(debug) cout << "RPCDet: " << theRPCName << " , on RPCStationLayer: " << RPCStationLayer << endl;
    return RPCStationLayer;
}

void RPCBendingAnalyzer::analyzeBending() {

    unsigned int Index = 0;
    fillSample(Index);
}

void RPCBendingAnalyzer::fillSample(unsigned int Index) {

    unsigned int RPCLayerSize = RPCLayer.size();

    if(Index >= RPCLayerSize) {
        if(debug) cout << "RPCBendingLayerIndex is outside the range!" << endl;
        return;
    }

    unsigned int simHitIndex = 0;
    SampleLayer[Index] = -1;
    for(vector<PSimHit>::const_iterator simHitIter = simHitsforTrack.begin(); simHitIter != simHitsforTrack.end(); simHitIter++) {
        int DetUnitId = simHitIter->detUnitId();
        RPCDetId sampleRPCDetId(DetUnitId);
        if(debug) cout << "simHit: ";
        int theRPCStationLayer1 = getRPCLayer(sampleRPCDetId);
        if(theRPCStationLayer1 == (int)RPCLayer[Index]) {
            SampleLayer[Index] = (int)RPCLayer[Index];
            simHitSample[RPCLayer[Index]] = *simHitIter;
            vector<RPCRecHit> theRPCrecHits = RPCrecHitsfromTrack[simHitIndex];
            for(vector<RPCRecHit>::const_iterator recHitIter = theRPCrecHits.begin(); recHitIter != theRPCrecHits.end(); recHitIter++) {
                RPCDetId RPCId = recHitIter->rpcId();
                if(debug) cout << "recHit: ";
                int theRPCStationLayer2 = getRPCLayer(sampleRPCDetId);
                if(theRPCStationLayer2 == (int)RPCLayer[Index]) {
                    recHitSample[RPCLayer[Index]] = *recHitIter;
                    if(Index < (RPCLayerSize - 1))
                        fillSample(Index+1);
                    else
                        fillBendingInfo();
                }
            }
        }
        simHitIndex++;
    }

    if(SampleLayer[Index] == -1) {
        SampleLayer[Index] = -1*(int)RPCLayer[Index];
        if(Index < (RPCLayerSize - 1))
            fillSample(Index+1);
        else
            fillBendingInfo();
    }

}

void RPCBendingAnalyzer::fillBendingInfo() {

    GlobalPoint simHitPosition[10];
    GlobalVector simHitMomentum[10];
    GlobalPoint recHitPosition[10];

    for(unsigned int i = 0; i < RPCLayer.size(); i++) {
        if(SampleLayer[i] >= 0) {
            int theDetUnitId = simHitSample[RPCLayer[i]].detUnitId();
            RPCDetId theRPCDetId(theDetUnitId);
            const GeomDetUnit *theRPCRoll = theRPCGeometry->idToDetUnit(theRPCDetId);
            Local3DPoint LocalPosition = simHitSample[RPCLayer[i]].localPosition();   
            GlobalPoint GlobalPosition = theRPCRoll->toGlobal(LocalPosition);
            simHitPosition[RPCLayer[i]] = GlobalPosition;
            LocalVector LocalMomentum = simHitSample[RPCLayer[i]].momentumAtEntry();
            GlobalVector GlobalMomentum = theRPCRoll->toGlobal(LocalMomentum);
            simHitMomentum[RPCLayer[i]] = GlobalMomentum;
        }
        else {
            simHitPosition[RPCLayer[i]] = GlobalPoint(0,0,0);
            simHitMomentum[RPCLayer[i]] = GlobalVector(0,0,0);
        }
    }

    for(unsigned int i = 0; i < RPCLayer.size(); i++) {
        if(SampleLayer[i] >= 0) {
            ClusterSize[RPCLayer[i]] = recHitSample[RPCLayer[i]].clusterSize();
            BX[RPCLayer[i]] = recHitSample[RPCLayer[i]].BunchX();
            simPtatRef[RPCLayer[i]] = simHitMomentum[RPCLayer[i]].perp();
            RPCDetId theRPCDetId = recHitSample[RPCLayer[i]].rpcId();
            const GeomDetUnit *theRPCRoll = theRPCGeometry->idToDetUnit(theRPCDetId);
            LocalPoint LocalPosition = recHitSample[RPCLayer[i]].localPosition();  
            GlobalPoint GlobalPosition = theRPCRoll->toGlobal(LocalPosition);
            recHitPosition[RPCLayer[i]] = GlobalPosition;
            Eta[RPCLayer[i]] = GlobalPosition.eta();
        }
        else {
            ClusterSize[RPCLayer[i]] = -1;
            BX[RPCLayer[i]] = 0;
            simPtatRef[RPCLayer[i]] = 0.;
            Eta[RPCLayer[i]] = 0.;
            recHitPosition[RPCLayer[i]] = GlobalPoint(0,0,0);
        }
    }

    for(unsigned int i = 0; i < RPCLayer.size(); i++) {
        if(SampleLayer[i] > 0) {
            simPtBendingPhi[RPCLayer[i]] = simHitMomentum[RPCLayer[i]].phi().value();
            simHitBendingPhi[RPCLayer[i]][RPCLayer[i]] = simHitPosition[RPCLayer[i]].phi().value();
            recHitBendingPhi[RPCLayer[i]][RPCLayer[i]] = recHitPosition[RPCLayer[i]].phi().value();
            for(unsigned int j = i+1; j < RPCLayer.size(); j++) {
                if(SampleLayer[j] > 0) {
                    simHitBendingPhi[RPCLayer[i]][RPCLayer[j]] = ((GlobalVector)(simHitPosition[RPCLayer[j]] - simHitPosition[RPCLayer[i]])).phi().value();
                    recHitBendingPhi[RPCLayer[i]][RPCLayer[j]] = ((GlobalVector)(recHitPosition[RPCLayer[j]] - recHitPosition[RPCLayer[i]])).phi().value();
                }
            }
        }
    }
    ExTree->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void RPCBendingAnalyzer::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void RPCBendingAnalyzer::endJob() {
    // save and close tree
    theFile->Write();
    theFile->Close();
}

//define this as a plug-in
DEFINE_FWK_MODULE(RPCBendingAnalyzer);
