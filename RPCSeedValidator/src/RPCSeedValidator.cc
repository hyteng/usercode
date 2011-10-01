// -*- C++ -*-
//
// Package:    RPCSeedValidator
// Class:      RPCSeedValidator
// 
/**\class RPCSeedValidator RPCSeedValidator.cc MyModule/RPCSeedValidator/src/RPCSeedValidator.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
 */
//
// Original Author:  Haiyun Teng
//         Created:  Thu Nov 20 01:40:00 CET 2008
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

// user special files
#include <FWCore/Framework/interface/ESHandle.h>
#include <SimDataFormats/Track/interface/SimTrack.h>
#include <SimDataFormats/TrackingHit/interface/PSimHit.h>
#include <SimDataFormats/RPCDigiSimLink/interface/RPCDigiSimLink.h>
#include <DataFormats/RPCRecHit/interface/RPCRecHit.h>
#include <DataFormats/DetId/interface/DetId.h>
#include <DataFormats/MuonDetId/interface/RPCDetId.h>
#include <Geometry/RPCGeometry/interface/RPCRoll.h>
#include <Geometry/CommonDetUnit/interface/GeomDetUnit.h>
#include <Geometry/RPCGeometry/interface/RPCGeometry.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <DataFormats/RPCRecHit/interface/RPCRecHitCollection.h>
#include <SimDataFormats/TrackingHit/interface/PSimHitContainer.h>
#include <SimDataFormats/Track/interface/SimTrackContainer.h>
#include <iostream>
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include <DataFormats/TrackingRecHit/interface/TrackingRecHit.h>
#include <DataFormats/GeometryVector/interface/Phi.h>
#include "TTree.h"
#include "TFile.h"
//C++
//#include <vector>


// using name space
using namespace std;
using namespace edm;

//
// class decleration
//

class BendingPhiIndexType {
    public:
        BendingPhiIndexType() {m[0] = 0; m[1] = 0; m[2] = 0; m[3] = 0;};
        BendingPhiIndexType(int a0, int a1, int a2, int a3) {m[0] = a0; m[1] = a1; m[2] = a2; m[3] = a3;};
        int m[4];
};

class RPCSeedValidator : public edm::EDAnalyzer {
    public:
        explicit RPCSeedValidator(const edm::ParameterSet&);
        ~RPCSeedValidator();


    private:
        virtual void beginJob();
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob();

        void getTrackInfo(const SimTrack& theSimTrack);
        void findSeedforTrack();
        void findSimHitfromSeedRecHit(const edm::ESHandle<RPCGeometry>& theRPCGeometry, std::vector<GlobalPoint>& SimHitRPCPos, std::vector<GlobalPoint>& RecHitRPCPos, const std::vector<GlobalPoint>& SimHitPositionCollection, const std::vector<unsigned int>& RecHitNumberforSimHits, const std::vector<RPCRecHit>& RPCRecHitsfromTrack, const std::vector<GlobalPoint>& SeedRecHitPositionCollection);
        bool SegmentFilter(int theFilterType);
        void scanOnePointCollection(int RPCLayer);
        void getBendingPhiCandidate();
        int getRPCLayer(const RPCRecHit& theRPCRecHit);
        void compareSeedBending();
        void getSimPtatRef();
        // ----------member data ---------------------------
        edm::InputTag SimHitTag_;
        edm::InputTag SimTrackTag_;
        edm::InputTag RPCRecHitTag_;
        edm::InputTag RPCDigiSimLinkTag_;
        edm::InputTag TrajectorySeedCollectionTag_;

        edm::Handle<PSimHitContainer> pSimHits;
        edm::Handle<SimTrackContainer> pSimTracks;
        edm::Handle<RPCRecHitCollection> pRPCRecHits;
        edm::Handle< edm::DetSetVector<RPCDigiSimLink> > theLinkDigis;
        edm::Handle<TrajectorySeedCollection> pTrajectorySeedCollection;
        edm::ESHandle<RPCGeometry> theRPCGeometry;

        bool debug;
        unsigned int RecHitNumberTH;
        unsigned int CodeTH;
        unsigned int unCodeTH;
        int FilterType;
        bool isVertexConstraint;
        double SegmentBendingPhiTH;
        double MaxBendingPhiTH;
        double SeedPurityTH;

        unsigned int RecHitNumber;
        unsigned int theCode;
        double TrackBendingPhi;
        std::vector<GlobalPoint> SimHitPositionCollection;
        std::vector<GlobalVector> SimHitMomentumCollection;
        std::vector<unsigned int> RecHitNumberforSimHits;
        std::vector<RPCRecHit> RPCRecHitsfromTrack;
        std::vector<GlobalPoint> SeedRecHitPositionCollection;
        int RecHitNumberinRPCLayer[11];
        std::vector< std::vector<GlobalPoint> > RPCRecHitLayer;
        GlobalPoint OnePointCollection[11];
        std::vector<BendingPhiIndexType> BendingFilter;
        std::vector<double> MaxBendingPhiCollection;
        std::vector<double> SegmentBendingPhiCollection;
        double MaxBendingPhi;
        double SegmentBendingPhi;

        std::string theRootFileName;
        TFile *theFile;
        TTree *ExTree;

        unsigned int EventNumber;
        int SimTrackId;
        int SimTrackType;
        double SimTrackMomentum;
        double SimTrackDirectionPhi;
        int SimTrackCharge;
        int SimTrackValid;
        bool PassSegmentFilter;
        DetId RefDet;
        double SimMomentumatRef;
        double SimDirectionPhiatRef;
        double SimBendingEntryPositionX;
        double SimBendingEntryPositionY;
        double SimBendingEntryPositionZ;
        double SimBendingLeavePositionX;
        double SimBendingLeavePositionY;
        double SimBendingLeavePositionZ;
        double SimBendingPhi;
        unsigned int SeedNumber;
        int SeedCharge;
        double SeedPurity;
        double SeedQuality;
        double RecMomentumatRef;
        double RecDirectionPhiatRef;
        double RecBendingEntryPositionX;
        double RecBendingEntryPositionY;
        double RecBendingEntryPositionZ;
        double RecBendingLeavePositionX;
        double RecBendingLeavePositionY;
        double RecBendingLeavePositionZ;
        double RecBendingPhi;
        double RecBendingLastPhi;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//
bool lessR(const GlobalPoint& it1, const GlobalPoint& it2) {
    // Don't need to use value() in Geom::Phi to short
    return (it1.perp() < it2.perp());
}
//
// constructors and destructor
//
RPCSeedValidator::RPCSeedValidator(const edm::ParameterSet& iConfig) {

    //now do what ever initialization is needed
    SimHitTag_ = iConfig.getUntrackedParameter<edm::InputTag>("SimHitTag");
    SimTrackTag_ = iConfig.getUntrackedParameter<edm::InputTag>("SimTrackTag");
    RPCRecHitTag_ = iConfig.getUntrackedParameter<edm::InputTag>("RPCRecHitTag"); 
    RPCDigiSimLinkTag_ = iConfig.getUntrackedParameter<edm::InputTag>("RPCDigiSimLinkTag");
    TrajectorySeedCollectionTag_ = iConfig.getUntrackedParameter<edm::InputTag>("TrajectorySeedCollectionTag");
    debug = iConfig.getUntrackedParameter<bool>("debug");
    RecHitNumberTH = iConfig.getUntrackedParameter<unsigned int>("RecHitNumberTH");
    CodeTH = iConfig.getUntrackedParameter<unsigned int>("CodeTH");
    unCodeTH = iConfig.getUntrackedParameter<unsigned int>("unCodeTH");
    FilterType = iConfig.getUntrackedParameter<int>("FilterType");
    isVertexConstraint = iConfig.getUntrackedParameter<bool>("isVertexConstraint");
    SegmentBendingPhiTH = iConfig.getUntrackedParameter<double>("SegmentBendingPhiTH");
    MaxBendingPhiTH = iConfig.getUntrackedParameter<double>("MaxBendingPhiTH");
    SeedPurityTH = iConfig.getUntrackedParameter<double>("SeedPurityTH");
    theRootFileName = iConfig.getUntrackedParameter<string>("theRootFileName", "RPCValidationTree.root");
    // Booking histogrqam
    theFile = new TFile(theRootFileName.c_str(), "recreate");
    theFile->cd();
    ExTree = new TTree("ExTree", "ExTree");
    ExTree->Branch("EventNumber", &EventNumber, "EventNumber/i");
    ExTree->Branch("SimTrackId", &SimTrackId, "SimTrackId/I");
    ExTree->Branch("SimTrackType", &SimTrackType,  "SimTrackType/I");
    ExTree->Branch("SimTrackMomentum", &SimTrackMomentum, "SimTrackMomentum/D");
    ExTree->Branch("SimTrackDirectionPhi", &SimTrackDirectionPhi, "SimTrackDirectionPhi/D");
    ExTree->Branch("SimTrackCharge", &SimTrackCharge, "SimTrackCharge/I");
    ExTree->Branch("SimTrackValid", &SimTrackValid, "SimTrackValid/I");
    ExTree->Branch("PassSegmentFilter", &PassSegmentFilter, "PassSegmentFilter/O");
    ExTree->Branch("SimMomentumatRef", &SimMomentumatRef, "SimMomentumatRef/D");
    ExTree->Branch("SimDirectionPhiatRef", &SimDirectionPhiatRef, "SimDirectionPhiatRef/D");
    ExTree->Branch("SimBendingEntryPositionX", &SimBendingEntryPositionX, "SimBendingEntryPositionX/D");
    ExTree->Branch("SimBendingEntryPositionY", &SimBendingEntryPositionY, "SimBendingEntryPositionY/D");
    ExTree->Branch("SimBendingEntryPositionZ", &SimBendingEntryPositionZ, "SimBendingEntryPositionZ/D");
    ExTree->Branch("SimBendingLeavePositionX", &SimBendingLeavePositionX, "SimBendingLeavePositionX/D");
    ExTree->Branch("SimBendingLeavePositionY", &SimBendingLeavePositionY, "SimBendingLeavePositionY/D");
    ExTree->Branch("SimBendingLeavePositionZ", &SimBendingLeavePositionZ, "SimBendingLeavePositionZ/D");
    ExTree->Branch("SimBendingPhi", &SimBendingPhi, "SimBendingPhi/D");
    ExTree->Branch("SeedNumber", &SeedNumber, "SeedNumber/i");
    ExTree->Branch("SeedCharge", &SeedCharge, "SeedCharge/I");
    ExTree->Branch("SeedPurity", &SeedPurity, "SeedPurity/D");
    ExTree->Branch("SeedQuality", &SeedQuality, "SeedQuality/D");
    ExTree->Branch("RecMomentumatRef", &RecMomentumatRef ,"RecMomentumatRef/D");
    ExTree->Branch("RecDirectionPhiatRef", &RecDirectionPhiatRef ,"RecDirectionPhiatRef/D");
    ExTree->Branch("RecBendingEntryPositionX", &RecBendingEntryPositionX, "RecBendingEntryPositionX/D");
    ExTree->Branch("RecBendingEntryPositionY", &RecBendingEntryPositionY, "RecBendingEntryPositionY/D");
    ExTree->Branch("RecBendingEntryPositionZ", &RecBendingEntryPositionZ, "RecBendingEntryPositionZ/D");
    ExTree->Branch("RecBendingLeavePositionX", &RecBendingLeavePositionX, "RecBendingLeavePositionX/D");
    ExTree->Branch("RecBendingLeavePositionY", &RecBendingLeavePositionY, "RecBendingLeavePositionY/D");
    ExTree->Branch("RecBendingLeavePositionZ", &RecBendingLeavePositionZ, "RecBendingLeavePositionZ/D");
    ExTree->Branch("RecBendingPhi", &RecBendingPhi ,"RecBendingPhi/D");
    ExTree->Branch("RecBendingLastPhi", &RecBendingLastPhi, "RecBendingLastPhi/D");
    // Set SimTrack index
    EventNumber = 0;
}


RPCSeedValidator::~RPCSeedValidator() {

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void RPCSeedValidator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    using namespace edm;    

    iEvent.getByLabel(SimHitTag_, pSimHits);
    iEvent.getByLabel(SimTrackTag_, pSimTracks);
    iEvent.getByLabel(RPCRecHitTag_, pRPCRecHits);
    iEvent.getByLabel(RPCDigiSimLinkTag_, theLinkDigis);
    iEvent.getByLabel(TrajectorySeedCollectionTag_, pTrajectorySeedCollection);
    iSetup.get<MuonGeometryRecord>().get(theRPCGeometry);

    EventNumber++;
    if(debug) cout << "Event "<< EventNumber << endl;
    unsigned int SimTrackIndex = 0;
    for(SimTrackContainer::const_iterator pSimTrack = pSimTracks->begin(); pSimTrack != pSimTracks->end(); pSimTrack++) {
        SimTrackIndex++;

        int SimTrackId = pSimTrack->trackId();
        int SimTrackType = pSimTrack->type();
        if(debug) cout << "SimTrack id: " << SimTrackId << endl;
        if(debug) cout << "SimTrack type: " << SimTrackType << endl;
        //if(debug) cout << "SimTrack momentum: " << pSimTrack->momentum().Pt() << endl;

        if(abs(SimTrackType) != 13)
            continue;

        getTrackInfo(*pSimTrack);
    }
}

void RPCSeedValidator::getTrackInfo(const SimTrack& theSimTrack) {

    // for test the bending of SimHits
    SimHitPositionCollection.clear();
    SimHitMomentumCollection.clear();
    RecHitNumberforSimHits.clear();
    RPCRecHitsfromTrack.clear();
    SeedRecHitPositionCollection.clear();

    // initiciate layer count array, 0-5 for barrel, 6-10 for endcap
    for(unsigned int i = 0; i < 11; i++)
        RecHitNumberinRPCLayer[i] = 0;

    SimTrackId = theSimTrack.trackId();
    SimTrackType = theSimTrack.type();
    if(debug) cout << "SimTrack id: " << SimTrackId << endl;
    if(debug) cout << "SimTrack type: " << SimTrackType << endl;

    SimTrackCharge = theSimTrack.charge();
    SimTrackMomentum = theSimTrack.momentum().pt();
    SimTrackDirectionPhi = theSimTrack.momentum().phi();
    if(debug) cout << "SimTrack momentum: " << SimTrackMomentum << endl;

    for(PSimHitContainer::const_iterator pSimHit = pSimHits->begin(); pSimHit != pSimHits->end(); pSimHit++) {
        int TrackId1 = pSimHit->trackId();
        if(TrackId1 != SimTrackId)
            continue;

        // The particle type of the hit may differ from the particle type of the SimTrack with id trackId(). 
        // This happends if the hit was created by a secondary track (e.g. a delta ray) originating from the trackId() and not existing as a separate SimTrack.
        int ParticleType1 = pSimHit->particleType();
        if(ParticleType1 != SimTrackType)
            continue;

        if(debug) cout << "SimHit particle type : " << ParticleType1 << endl;
        int DetUnitId1 = pSimHit->detUnitId();
        DetId DetectorId1(DetUnitId1);
        GeomDetEnumerators::SubDetector subDet = theRPCGeometry->idToDetUnit(DetectorId1)->subDetector();
        // Now we only check barrel seed
        if(subDet != GeomDetEnumerators::RPCBarrel)
            continue;

        if(debug) cout << "SimHit DetUnit id: " << DetUnitId1 << endl;
        Local3DPoint HitPoint1 = pSimHit->entryPoint();
        if(debug) cout << "SimHit local hit point: " << HitPoint1 << endl;
        LocalVector Momentum1 = pSimHit->momentumAtEntry();
        if(debug) cout << "SimHit local momentum: " << Momentum1 << endl;
        unsigned int processType = pSimHit->processType();
        if(debug) cout << "SimHit process type: " << processType << endl;

        const RPCRoll * theRPCRoll = theRPCGeometry->roll(DetectorId1);
        GlobalPoint SimHitPosition1 = theRPCRoll->toGlobal(HitPoint1);
        GlobalVector GlobalMomentum1 = theRPCRoll->toGlobal(Momentum1);
        if(debug) cout << "At this SimHit the global Position is " << SimHitPosition1 << endl;
        if(debug) cout << "At this SimHit the global Pt value: " << GlobalMomentum1.perp() << ", Phi: " << GlobalMomentum1.phi() << ", in degree: " << GlobalMomentum1.phi().degrees() << endl;
        SimHitPositionCollection.push_back(SimHitPosition1);
        SimHitMomentumCollection.push_back(GlobalMomentum1);
        int FirstStrip = -1;
        int LastStrip = -1;
        int DigiBX = 0;
        if(debug) cout << "Checking DigiSimLink..." << endl;
        for(edm::DetSetVector<RPCDigiSimLink>::const_iterator LinkIter = theLinkDigis->begin(); LinkIter != theLinkDigis->end(); LinkIter++) {
            for(edm::DetSet<RPCDigiSimLink>::const_iterator DigiIter = LinkIter->data.begin(); DigiIter != LinkIter->data.end(); ++DigiIter) {
                int DetUnitId2 = DigiIter->getDetUnitId();
                int TrackId2 = DigiIter->getTrackId();
                Local3DPoint HitPoint2 = DigiIter->getEntryPoint();
                LocalVector Momentum2 = DigiIter->getMomentumAtEntry();
                int ParticleType2 = DigiIter->getParticleType();
                if((DetUnitId1 == DetUnitId2) && (TrackId1 == TrackId2) && (ParticleType1 == ParticleType2)) {
                    int Strip = DigiIter->getStrip();
                    DigiBX = DigiIter->getBx();
                    if(debug) cout << "Find Digi's Strip: " << Strip << endl;
                    if(debug) cout << "Find Digi's BX: " << DigiBX << endl;
                    if(debug) cout << "Find Digi's hit positin: " << HitPoint2 << endl;
                    if(FirstStrip == -1) {
                        FirstStrip = Strip;
                        LastStrip = Strip;
                    }
                    if(Strip < FirstStrip)
                        FirstStrip = Strip;
                    if(Strip > LastStrip)
                        LastStrip = Strip;
                }
            }
        }

        RPCRecHitCollection::range RPCRecHitRange = pRPCRecHits->get(DetectorId1);
        int AssociatedRecHitforSimHit = 0;
        for(RPCRecHitCollection::const_iterator RPCRecHitIter = RPCRecHitRange.first; RPCRecHitIter != RPCRecHitRange.second; RPCRecHitIter++) {        
            int FirstRecStrip = RPCRecHitIter->firstClusterStrip();
            int ClusterSize = RPCRecHitIter->clusterSize();
            int RecBX = RPCRecHitIter->BunchX();
            bool isValid = RPCRecHitIter->isValid();
            if(debug) cout << "RecHit 1st strip: " << FirstRecStrip << endl;
            if(debug) cout << "RecHit ClusterSize: " << ClusterSize << endl;
            if(debug) cout << "RecHit BX: " << RecBX << endl;
            if(debug) cout << "RecHit Valid: " << isValid << endl;
            if((FirstRecStrip <= LastStrip)  && ((FirstRecStrip + ClusterSize - 1) >= FirstStrip) && isValid && RecBX == DigiBX) {
                if(debug) cout << "Find SimHit-RecHit" << endl;
                LocalPoint RecLocalPosition = RPCRecHitIter->localPosition();
                const GeomDetUnit *theRPCRoll = theRPCGeometry->idToDetUnit(DetectorId1);
                GlobalPoint RecGlobalPosition = theRPCRoll->toGlobal(RecLocalPosition);
                if(debug) cout << " GlobalPosition: " << RecGlobalPosition << endl;
                AssociatedRecHitforSimHit++;
                RPCRecHitsfromTrack.push_back(*RPCRecHitIter);
                int theRPCLayer = getRPCLayer(*RPCRecHitIter);
                if(debug) cout << "Push this RecHit to layer " << theRPCLayer << " for Validation." << endl;
                RecHitNumberinRPCLayer[theRPCLayer] += 1;
            }
        }
        RecHitNumberforSimHits.push_back(AssociatedRecHitforSimHit);

        if(AssociatedRecHitforSimHit > 1)
            if(debug) cout << "Find " << AssociatedRecHitforSimHit << " RecHits in the same DetUnit" << endl;
    }

    // Code for checking layer occupancy
    RecHitNumber = 0;
    theCode = 0;
    for(unsigned int RPCLayerIndex = 0; RPCLayerIndex < 11; RPCLayerIndex++) {
        if(RecHitNumberinRPCLayer[RPCLayerIndex] != 0) {
            RecHitNumber++;
            unsigned value = 1;
            for(unsigned int i = 0; i < RPCLayerIndex; i++)
                value *= 2;
            theCode += value;
        }
    }

    if((RecHitNumber >= RecHitNumberTH) && (((int)theCode & (int)CodeTH) == (int)CodeTH) && ((int)theCode & (int)unCodeTH) != 0) {
        SimTrackValid = 1;
        PassSegmentFilter = SegmentFilter(FilterType);
    }
    else {
        SimTrackValid = 0;
        PassSegmentFilter = false;
    }

    if(debug) cout << "Core for this track is " << theCode << ", SimTrack Valid is " << SimTrackValid << ", PassSegmentFilter: " << PassSegmentFilter << endl;

    findSeedforTrack();
}

void RPCSeedValidator::findSeedforTrack() {

    if(debug) cout << "Finding corresponding RecHits in " << pTrajectorySeedCollection->size() << " seeds..." << endl;

    // Set the default value for no seed case and always fill it for each SimTrack
    SimMomentumatRef = -1;
    SimDirectionPhiatRef = 0;
    SimBendingEntryPositionX = 0;
    SimBendingEntryPositionY = 0;
    SimBendingEntryPositionZ = 0;
    SimBendingLeavePositionX = 0;
    SimBendingLeavePositionY = 0;
    SimBendingLeavePositionZ = 0;
    SimBendingPhi = 0;
    SeedNumber = 0; 
    SeedCharge = 0;
    SeedPurity = 0.;
    SeedQuality = 0;
    RecMomentumatRef = -1;
    RecDirectionPhiatRef = 0;
    RecBendingEntryPositionX = 0;
    RecBendingEntryPositionY = 0;
    RecBendingEntryPositionZ = 0;
    RecBendingLeavePositionX = 0;
    RecBendingLeavePositionY = 0;
    RecBendingLeavePositionZ = 0;
    RecBendingPhi = 0;
    RecBendingLastPhi = 0;
    ExTree->Fill();

    // For Valid track check the good seed and fake seed of this track
    if(SimTrackValid == 1) {
        for(TrajectorySeedCollection::const_iterator RPCSeedIter = pTrajectorySeedCollection->begin(); RPCSeedIter != pTrajectorySeedCollection->end(); RPCSeedIter++) {
            unsigned int NumberofRecHitsinSeed = 0;
            TrajectorySeed::range SeedRecHitsRange = RPCSeedIter->recHits();
            SeedRecHitPositionCollection.clear();
            for(TrajectorySeed::const_iterator RPCSeedRecHitIter = SeedRecHitsRange.first; RPCSeedRecHitIter != SeedRecHitsRange.second; RPCSeedRecHitIter++) {
                LocalPoint SeedRecHitLocalPosition = RPCSeedRecHitIter->localPosition();
                DetId SeedRPCId = RPCSeedRecHitIter->geographicalId();
                const GeomDetUnit *theRPCRoll = theRPCGeometry->idToDetUnit(SeedRPCId);
                GlobalPoint SeedRecHitGlobalPosition = theRPCRoll->toGlobal(SeedRecHitLocalPosition);
                SeedRecHitPositionCollection.push_back(SeedRecHitGlobalPosition);
                if(debug) cout << "GlobalPosition in seed's RecHit: " << SeedRecHitGlobalPosition << ". On DetId: " << SeedRPCId.rawId() << endl;
                for(std::vector<RPCRecHit>::const_iterator it = RPCRecHitsfromTrack.begin(); it != RPCRecHitsfromTrack.end(); it++) {
                    // This is mot implented in RPCRecHit class
                    //TrackingRecHit::SharedInputType type = (TrackingRecHit::SharedInputType)(0);
                    //bool isthesame = RPCSeedRecHitIter->sharesInput(&(*it), type);
                    LocalPoint RecHitLocalPosition = it->localPosition();
                    RPCDetId RecRPCId = it->rpcId();
                    const GeomDetUnit *theRPCRoll = theRPCGeometry->idToDetUnit(RecRPCId);
                    GlobalPoint RecHitGlobalPosition = theRPCRoll->toGlobal(RecHitLocalPosition);
                    if(debug) cout << "GlobalPosition in stack's RecHit: " << RecHitGlobalPosition << ". On RPCRollId: " << RecRPCId.rawId() << endl;
                    if(SeedRecHitGlobalPosition.x() == RecHitGlobalPosition.x() && SeedRecHitGlobalPosition.y() == RecHitGlobalPosition.y() && SeedRecHitGlobalPosition.z() == RecHitGlobalPosition.z()) {
                        if(debug) cout << "Find one RecHit in seed corresponding to track " << SimTrackId << endl;
                        NumberofRecHitsinSeed++;
                    }
                }
            }

            if(debug) cout << "Find " << NumberofRecHitsinSeed << " RecHits in this seed for track " << SimTrackId << endl;
            SeedPurity = (double)NumberofRecHitsinSeed / (double)RPCSeedIter->nHits();
            if(debug) cout << "SeedPurity is " << SeedPurity << endl;

            if(SeedPurity >= SeedPurityTH) {
                SeedNumber++;
                LocalVector SeedMomentum = RPCSeedIter->startingState().parameters().momentum();
                std::vector<float> SeedError = RPCSeedIter->startingState().errorMatrix();
                SeedQuality = SeedError[0];
                if(debug) cout << "Find SeedQuality " << SeedQuality << " in this seed for track " << SimTrackId << endl;
                SeedCharge = RPCSeedIter->startingState().parameters().charge();
                if(debug) cout << "Charge of Rec is: " << SeedCharge << ", charge of Sim is: " << SimTrackCharge << endl;
                RefDet = RPCSeedIter->startingState().detId();
                GlobalVector GP = theRPCGeometry->idToDetUnit(RefDet)->toGlobal(SeedMomentum);
                RecMomentumatRef = GP.perp();
                RecDirectionPhiatRef = GP.phi().value();
                
                if(SeedPurity == 1.)
                    compareSeedBending();
                else {
                    SimBendingEntryPositionX = 0;
                    SimBendingEntryPositionY = 0;
                    SimBendingEntryPositionZ = 0;
                    SimBendingLeavePositionX = 0;
                    SimBendingLeavePositionY = 0;
                    SimBendingLeavePositionZ = 0;
                    SimBendingPhi = 0;
                    RecBendingEntryPositionX = 0;
                    RecBendingEntryPositionY = 0;
                    RecBendingEntryPositionZ = 0;
                    RecBendingLeavePositionX = 0;
                    RecBendingLeavePositionY = 0;
                    RecBendingLeavePositionZ = 0;
                    RecBendingPhi = 0;
                    RecBendingLastPhi = 0;
                }
                
                getSimPtatRef();

                ExTree->Fill();

                if(debug) cout << "SimTrackId: " << SimTrackId << ", SimTrackType: " << SimTrackType << ", SeedNumber: " << SeedNumber << ", SeedPurity: " << SeedPurity << ", SeedCharge: " << SeedCharge << ", SimTrackCharge: " << SimTrackCharge << ", SimBendingPhi: " << SimBendingPhi << ", RecBendingPhi: " << RecBendingPhi << ", RecMomentumatRef: " << RecMomentumatRef << ", SimMomentumatRef: " << SimMomentumatRef << endl;
                if(debug) cout << "Find " << SeedNumber << "th seed in this track " << SimTrackId << endl;
            }
        }
    }
    // invailid seed
    else {
        for(TrajectorySeedCollection::const_iterator RPCSeedIter = pTrajectorySeedCollection->begin(); RPCSeedIter != pTrajectorySeedCollection->end(); RPCSeedIter++) {
            unsigned int NumberofRecHitsinSeed = 0;
            TrajectorySeed::range SeedRecHitsRange = RPCSeedIter->recHits();
            for(TrajectorySeed::const_iterator RPCSeedRecHitIter = SeedRecHitsRange.first; RPCSeedRecHitIter != SeedRecHitsRange.second; RPCSeedRecHitIter++) {
                LocalPoint SeedRecHitLocalPosition = RPCSeedRecHitIter->localPosition();
                DetId SeedRPCId = RPCSeedRecHitIter->geographicalId();
                const GeomDetUnit *theRPCRoll = theRPCGeometry->idToDetUnit(SeedRPCId);
                GlobalPoint SeedRecHitGlobalPosition = theRPCRoll->toGlobal(SeedRecHitLocalPosition);
                if(debug) cout << " GlobalPosition in seed's RecHit: " << SeedRecHitGlobalPosition << endl;
                for(std::vector<RPCRecHit>::const_iterator it = RPCRecHitsfromTrack.begin(); it != RPCRecHitsfromTrack.end(); it++) {
                    // This is mot implented in RPCRecHit class
                    //TrackingRecHit::SharedInputType type = (TrackingRecHit::SharedInputType)(0);
                    //bool isthesame = RPCSeedRecHitIter->sharesInput(&(*it), type);
                    LocalPoint RecHitLocalPosition = it->localPosition();
                    DetId RecRPCId = it->rpcId();
                    const GeomDetUnit *theRPCRoll = theRPCGeometry->idToDetUnit(RecRPCId);
                    GlobalPoint RecHitGlobalPosition = theRPCRoll->toGlobal(RecHitLocalPosition);
                    if(debug) cout << " GlobalPosition in stack's RecHit: " << RecHitGlobalPosition << endl;
                    if(SeedRecHitGlobalPosition.x() == RecHitGlobalPosition.x() && SeedRecHitGlobalPosition.y() == RecHitGlobalPosition.y() && SeedRecHitGlobalPosition.z() == RecHitGlobalPosition.z()) {
                        if(debug) cout << "Find one RecHit in seed corresponding to track " << SimTrackId << endl;
                        NumberofRecHitsinSeed++;
                    }
                }
            }
            SeedPurity = (double)NumberofRecHitsinSeed / (double)RPCSeedIter->nHits();

            if(debug) cout << "Find " << NumberofRecHitsinSeed << " RecHits in this seed for track " << SimTrackId << endl;

            if(SeedPurity > 0) {
                SeedNumber++;
                std::vector<float> SeedError = RPCSeedIter->startingState().errorMatrix();
                SeedQuality = SeedError[0];
                if(debug) cout << "Find SeedQuality " << SeedQuality << " in this seed for invailid track " << SimTrackId << endl;
                SeedCharge = RPCSeedIter->startingState().parameters().charge();
                if(debug) cout << "Charge of Rec is: " << SeedCharge << ", charge of Sim is: " << SimTrackCharge << endl;

                ExTree->Fill();
            }
        }
    }
}

int RPCSeedValidator::getRPCLayer(const RPCRecHit& theRPCRecHit) {
    int theRPCLayer;
    RPCDetId RPCId = theRPCRecHit.rpcId();
    int Region = RPCId.region();
    int Ring = RPCId.ring();
    int Station = RPCId.station();
    int Layer = RPCId.layer();
    if(debug) cout << "RPC det is " << Region << ", " << Ring << ", " << Station << ", " << Layer << endl;
    if(Region == 0) {
        if(Station >= 3)
            theRPCLayer = Station + 1;
        else {
            if(Station == 2)
                theRPCLayer = Layer + 1;
            else
                theRPCLayer = Layer - 1;
        }
    }
    else {
        // Before endcap double layer implement, RE1-4
        theRPCLayer = Station + 5;
        // After endcap double layer implement, RE2in/out may included
    }
    return theRPCLayer;
}

bool RPCSeedValidator::SegmentFilter(int theFilterType) {

    RPCRecHitLayer.clear();
    for(unsigned int i = 0; i < 11; i++) {
        std::vector<GlobalPoint> RPCRecHitTempLayer;
        RPCRecHitTempLayer.resize(0);
        RPCRecHitLayer.push_back(RPCRecHitTempLayer);
    }

    for(std::vector<RPCRecHit>::const_iterator RPCRecHitIter = RPCRecHitsfromTrack.begin(); RPCRecHitIter != RPCRecHitsfromTrack.end(); RPCRecHitIter++) {
        LocalPoint RecLocalPosition = RPCRecHitIter->localPosition();
        RPCDetId RecRPCId = RPCRecHitIter->rpcId();
        const GeomDetUnit *theRPCRoll = theRPCGeometry->idToDetUnit(RecRPCId);
        GlobalPoint RecGlobalPosition = theRPCRoll->toGlobal(RecLocalPosition);
        int RPCLayer = getRPCLayer(*RPCRecHitIter);
        if(debug) cout << "Track's RecHit in RPCLayer: " << RPCLayer << ", GlobalPosition: " << RecGlobalPosition << endl;
        RPCRecHitLayer[RPCLayer].push_back(RecGlobalPosition);
    }

    if((FilterType == 1 || FilterType == 0) && isVertexConstraint == true) {
        BendingFilter.push_back(BendingPhiIndexType(0,1,2,3));
        BendingFilter.push_back(BendingPhiIndexType(0,1,1,4));
        BendingFilter.push_back(BendingPhiIndexType(2,3,3,4));
        BendingFilter.push_back(BendingPhiIndexType(0,1,1,5));
        BendingFilter.push_back(BendingPhiIndexType(2,3,3,5));
        BendingFilter.push_back(BendingPhiIndexType(0,0,0,1));
        BendingFilter.push_back(BendingPhiIndexType(2,2,2,3));
    }
    if((FilterType == 1 || FilterType == 0) && isVertexConstraint == false) {
        BendingFilter.push_back(BendingPhiIndexType(0,1,2,3));
        BendingFilter.push_back(BendingPhiIndexType(0,1,1,4));
        BendingFilter.push_back(BendingPhiIndexType(2,3,3,4));
        BendingFilter.push_back(BendingPhiIndexType(0,1,1,5));
        BendingFilter.push_back(BendingPhiIndexType(2,3,3,5));
    }
    if((FilterType == 2 || FilterType == 0) && isVertexConstraint == true) {
        BendingFilter.push_back(BendingPhiIndexType(0,1,1,4));
        BendingFilter.push_back(BendingPhiIndexType(0,1,1,5));
        BendingFilter.push_back(BendingPhiIndexType(0,0,0,1));
    }
    if((FilterType == 2 || FilterType == 0) && isVertexConstraint == false) {
        BendingFilter.push_back(BendingPhiIndexType(0,1,1,4));
        BendingFilter.push_back(BendingPhiIndexType(0,1,1,5));
    }

    MaxBendingPhiCollection.clear();
    SegmentBendingPhiCollection.clear();
    scanOnePointCollection(0);
    
    MaxBendingPhi = 0.;
    for(unsigned int i = 0; i < MaxBendingPhiCollection.size(); i++) {
        if(fabs(MaxBendingPhi) < MaxBendingPhiCollection[i])
            MaxBendingPhi = MaxBendingPhiCollection[i];        
    }
    SegmentBendingPhi = 0.;
    for(unsigned int i = 0; i < SegmentBendingPhiCollection.size(); i++) {
        if(fabs(SegmentBendingPhi) < SegmentBendingPhiCollection[i])
            SegmentBendingPhi = SegmentBendingPhiCollection[i];
    }

    bool theSegmentFilterResult = false;
    if(MaxBendingPhi >= MaxBendingPhiTH || SegmentBendingPhi >= SegmentBendingPhiTH)
        theSegmentFilterResult = true;
        
    if(FilterType < 0)
        theSegmentFilterResult = true;

    if(debug) cout << "Valid track TrackBendingPhi: " << TrackBendingPhi << ". theSegmentFilterResult: " << theSegmentFilterResult << endl;
    return theSegmentFilterResult;
}

void RPCSeedValidator::scanOnePointCollection(int RPCLayer) {
    if(debug) cout << "Position on layer " << RPCLayer << endl;
    if(RPCRecHitLayer[RPCLayer].size() > 0)
        for(unsigned int i = 0; i < RPCRecHitLayer[RPCLayer].size(); i++) {
            OnePointCollection[RPCLayer] = RPCRecHitLayer[RPCLayer][i];
            if(RPCLayer < 10)
                scanOnePointCollection(RPCLayer+1);
            else
                getBendingPhiCandidate();
        }
    else {
        OnePointCollection[RPCLayer] = GlobalPoint(0,0,0);
        if(RPCLayer < 10)
            scanOnePointCollection(RPCLayer+1);
        else
            getBendingPhiCandidate();
    }
}

void RPCSeedValidator::getBendingPhiCandidate() {
    double TempMaxBendingPhi = 0.;
    double TempBendingPhi;
    double TempSegmentBendingPhi = 0.;
    for(unsigned int Index = 0; Index < BendingFilter.size(); Index++) {
        int i = BendingFilter[Index].m[0];
        int j = BendingFilter[Index].m[1];
        int k = BendingFilter[Index].m[2];
        int l = BendingFilter[Index].m[3];
        if(OnePointCollection[i] == GlobalPoint(0,0,0) || OnePointCollection[j] == GlobalPoint(0,0,0) || OnePointCollection[k] == GlobalPoint(0,0,0) || OnePointCollection[l] == GlobalPoint(0,0,0))
            continue;
        Geom::Phi<float> SegmentPhi1;
        Geom::Phi<float> SegmentPhi2;
        if(i == j)
            SegmentPhi1 = ((GlobalVector)(OnePointCollection[i] - GlobalPoint(0,0,0))).phi();
        else
            SegmentPhi1 = ((GlobalVector)(OnePointCollection[j] - OnePointCollection[i])).phi();
        SegmentPhi2 = ((GlobalVector)(OnePointCollection[l] - OnePointCollection[k])).phi();
        TempBendingPhi = (SegmentPhi2-SegmentPhi1).value();
        if(fabs(TempMaxBendingPhi) < fabs(TempBendingPhi))
            TempMaxBendingPhi = TempBendingPhi;

        if(Index == 0)
            TempSegmentBendingPhi = TempBendingPhi;
    }
    MaxBendingPhiCollection.push_back(TempMaxBendingPhi);
    SegmentBendingPhiCollection.push_back(TempSegmentBendingPhi);
}

void RPCSeedValidator::compareSeedBending() {
    if(debug) cout << "Comparing the Sim and rec bending for " << SeedRecHitPositionCollection.size() << " seed hits." << endl;
    std::vector<GlobalPoint> SimHitRPCPos;
    std::vector<GlobalPoint> RecHitRPCPos;
    SimHitRPCPos.clear();
    RecHitRPCPos.clear();
    sort(SeedRecHitPositionCollection.begin(), SeedRecHitPositionCollection.end(), lessR);
    findSimHitfromSeedRecHit(theRPCGeometry, SimHitRPCPos, RecHitRPCPos, SimHitPositionCollection, RecHitNumberforSimHits, RPCRecHitsfromTrack, SeedRecHitPositionCollection);
    
    if(FilterType == 1) {
        GlobalVector SimHitRPCSeg1 = SimHitRPCPos[1] - SimHitRPCPos[0];
        GlobalVector SimHitRPCSeg2 = SimHitRPCPos[3] - SimHitRPCPos[2];
        SimBendingEntryPositionX = SimHitRPCPos[1].x();
        SimBendingEntryPositionY = SimHitRPCPos[1].y();
        SimBendingEntryPositionZ = SimHitRPCPos[1].z();
        SimBendingLeavePositionX = SimHitRPCPos[2].x();
        SimBendingLeavePositionY = SimHitRPCPos[2].y();
        SimBendingLeavePositionZ = SimHitRPCPos[2].z();
        SimBendingPhi = (SimHitRPCSeg2.phi()-SimHitRPCSeg1.phi()).value();
        GlobalVector RecHitRPCSeg1 = RecHitRPCPos[1] - RecHitRPCPos[0];
        GlobalVector RecHitRPCSeg2 = RecHitRPCPos[3] - RecHitRPCPos[2];
        RecBendingEntryPositionX = RecHitRPCPos[1].x();
        RecBendingEntryPositionY = RecHitRPCPos[1].y();
        RecBendingEntryPositionZ = RecHitRPCPos[1].z();
        RecBendingLeavePositionX = RecHitRPCPos[2].x();
        RecBendingLeavePositionY = RecHitRPCPos[2].y();
        RecBendingLeavePositionZ = RecHitRPCPos[2].z();
        GlobalVector RecHitRPCSegLast = RecHitRPCPos[RecHitRPCPos.size()-1] - RecHitRPCPos[1];
        RecBendingLastPhi = (RecHitRPCSegLast.phi() - RecHitRPCSeg1.phi()).value();
        RecBendingPhi = (RecHitRPCSeg2.phi()-RecHitRPCSeg1.phi()).value();
    }
    if(FilterType == 2) {
        GlobalVector SimHitRPCSeg1 = SimHitRPCPos[1] - SimHitRPCPos[0];
        GlobalVector SimHitRPCSeg2 = SimHitRPCPos[2] - SimHitRPCPos[1];
        SimBendingEntryPositionX = SimHitRPCPos[1].x();
        SimBendingEntryPositionY = SimHitRPCPos[1].y();
        SimBendingEntryPositionZ = SimHitRPCPos[1].z();
        SimBendingLeavePositionX = SimHitRPCPos[2].x();
        SimBendingLeavePositionY = SimHitRPCPos[2].y();
        SimBendingLeavePositionZ = SimHitRPCPos[2].z();
        SimBendingPhi = (SimHitRPCSeg2.phi()-SimHitRPCSeg1.phi()).value();
        GlobalVector RecHitRPCSeg1 = RecHitRPCPos[1] - RecHitRPCPos[0];
        GlobalVector RecHitRPCSeg2 = RecHitRPCPos[2] - RecHitRPCPos[1];
        RecBendingEntryPositionX = RecHitRPCPos[1].x();
        RecBendingEntryPositionY = RecHitRPCPos[1].y();
        RecBendingEntryPositionZ = RecHitRPCPos[1].z();
        RecBendingLeavePositionX = RecHitRPCPos[2].x();
        RecBendingLeavePositionY = RecHitRPCPos[2].y();
        RecBendingLeavePositionZ = RecHitRPCPos[2].z();
        GlobalVector RecHitRPCSegLast = RecHitRPCPos[RecHitRPCPos.size()-1] - RecHitRPCPos[1];
        RecBendingLastPhi = (RecHitRPCSegLast.phi() - RecHitRPCSeg1.phi()).value();
        RecBendingPhi = (RecHitRPCSeg2.phi()-RecHitRPCSeg1.phi()).value();
    }
    
    if(debug) cout << "SimBendingPhi is " << SimBendingPhi << ", RecBendingPhi is " << RecBendingPhi << ", ratio is " << SimBendingPhi/RecBendingPhi  << ", diffbendPhi is " << (RecBendingPhi-SimBendingPhi) << ", RecBendingLastPhi is " << RecBendingLastPhi << endl;
}

void RPCSeedValidator::getSimPtatRef() {
    bool isset = false;
    for(PSimHitContainer::const_iterator pSimHit = pSimHits->begin(); pSimHit != pSimHits->end(); pSimHit++) {
        int TrackId1 = pSimHit->trackId();
        if (TrackId1 != SimTrackId)
            continue;
        int Particletype1 = pSimHit->particleType();
        if(Particletype1 != SimTrackType)
            continue;
        DetId DetectorId1(pSimHit->detUnitId());
        if(RefDet == DetectorId1 && isset == false) {
            LocalVector LocalMomentum = pSimHit->momentumAtEntry();
            const RPCRoll * theRPCRoll = theRPCGeometry->roll(RefDet);
            GlobalVector GlobalMomentum = theRPCRoll->toGlobal(LocalMomentum);
            GlobalPoint SimEntryPosition = theRPCRoll->toGlobal(pSimHit->entryPoint());
            GlobalPoint SimLeavePosition = theRPCRoll->toGlobal(pSimHit->exitPoint());
            SimMomentumatRef = GlobalMomentum.perp();
            SimDirectionPhiatRef = GlobalMomentum.phi();
            if(debug) cout << "@ ref SimHit's entry Position is: " << SimEntryPosition.x() << ", " << SimEntryPosition.y() << ", " << SimEntryPosition.z() << endl;
            if(debug) cout << "@ ref SimHit's leave Position is: " << SimLeavePosition.x() << ", " << SimLeavePosition.y() << ", " << SimLeavePosition.z() << endl;
            if(debug) cout << "@ ref SimHit's Pt and direction is: " << SimMomentumatRef << ", " << SimDirectionPhiatRef << endl;
            if(debug) cout << "@ ref RecHit's Pt and direction is: " << RecMomentumatRef << ", " << RecDirectionPhiatRef << endl;
            isset = true;
        }
    }
}

void RPCSeedValidator::findSimHitfromSeedRecHit(const edm::ESHandle<RPCGeometry>& theRPCGeometry, std::vector<GlobalPoint>& SimHitRPCPos, std::vector<GlobalPoint>& RecHitRPCPos, const std::vector<GlobalPoint>& SimHitPositionCollection, const std::vector<unsigned int>& RecHitNumberforSimHits, const std::vector<RPCRecHit>& RPCRecHitsfromTrack, const std::vector<GlobalPoint>& SeedRecHitPositionCollection) {

    for(std::vector<GlobalPoint>::const_iterator seedRecHitPositionIter = SeedRecHitPositionCollection.begin(); seedRecHitPositionIter != SeedRecHitPositionCollection.end(); seedRecHitPositionIter++) {
        if(debug) cout << "Checking one RecHit position in seed: " << (*seedRecHitPositionIter) << endl;
        //RecHitRPCPos.push_back((*seedRecHitPositionIter));
        unsigned int RecHitIndex = 0;
        unsigned int SimHitIndex = 0;
        unsigned int NumberofRecHitforSimHit = RecHitNumberforSimHits[SimHitIndex];
        for(RecHitIndex = 0; RecHitIndex < RPCRecHitsfromTrack.size(); RecHitIndex++) {
            // Find the SimHit for current RecHit
            while(NumberofRecHitforSimHit == 0) {
                SimHitIndex++;
                NumberofRecHitforSimHit = RecHitNumberforSimHits[SimHitIndex];
            }
            if(debug) cout << "Find SimHit index " << SimHitIndex << " and it has " << NumberofRecHitforSimHit << " RecHits." << endl;
            LocalPoint lpos = RPCRecHitsfromTrack[RecHitIndex].localPosition();
            RPCDetId RPCId = RPCRecHitsfromTrack[RecHitIndex].rpcId();
            const GeomDetUnit *rpcRoll = theRPCGeometry->idToDetUnit(RPCId);
            GlobalPoint gpos = rpcRoll->toGlobal(lpos);
            //if(debug) cout << endl;
            if(debug) cout << "Checking for RecHit from track, position: " << gpos << endl;
            if((gpos.x() == seedRecHitPositionIter->x()) && (gpos.y() == seedRecHitPositionIter->y()) && (gpos.z() == seedRecHitPositionIter->z())) {
                SimHitRPCPos.push_back(SimHitPositionCollection[SimHitIndex]);
                RecHitRPCPos.push_back((*seedRecHitPositionIter));
                if(debug) cout << "Find this SimHit " << SimHitIndex << " and save its position." << endl;
                break;
            }
            NumberofRecHitforSimHit--;
            //if(NumberofRecHitforSimHit == 0) {
            //    SimHitIndex++;
            //    NumberofRecHitforSimHit = RecHitNumberforSimHits[SimHitIndex];
            //}
        }
    }
    if(debug) cout << "All " << RecHitRPCPos.size() << " seed's RecHit and correspond SimHit position is saved in to array" << endl;
}

// ------------ method called once each job just before starting event loop  ------------
void RPCSeedValidator::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void RPCSeedValidator::endJob() {
    // save and close tree
    theFile->Write();
    theFile->Close();
}

//define this as a plug-in
DEFINE_FWK_MODULE(RPCSeedValidator);
