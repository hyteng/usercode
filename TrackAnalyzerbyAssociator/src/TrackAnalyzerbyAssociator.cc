// -*- C++ -*-
//
// Package:    TrackAnalyzerbyAssociator
// Class:      TrackAnalyzerbyAssociator
// 
/**\class TrackAnalyzerbyAssociator TrackAnalyzerbyAssociator.cc MyModule/TrackAnalyzerbyAssociator/src/TrackAnalyzerbyAssociator.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
 */
//
// Original Author:  Haiyun Teng
//         Created:  Tue Nov 10 16:30:59 CET 2009
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

#include "FWCore/Framework/interface/EventSetup.h"
#include <FWCore/Utilities/interface/InputTag.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h>
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include <DataFormats/TrackReco/interface/Track.h>
#include <SimMuon/MCTruth/interface/DTHitAssociator.h>
#include <SimMuon/MCTruth/interface/MuonTruth.h>
#include <SimMuon/MCTruth/interface/RPCHitAssociator.h>
#include <SimDataFormats/EncodedEventId/interface/EncodedEventId.h>
#include <DataFormats/DetId/interface/DetId.h>
#include <Geometry/CommonDetUnit/interface/GeomDetUnit.h>
#include <iostream>
#include <Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h>
#include <Geometry/Records/interface/GlobalTrackingGeometryRecord.h>
#include <TTree.h>
#include <TFile.h>
#include <SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h>
#include <SimDataFormats/TrackingHit/interface/PSimHit.h>
#include <TrackingTools/GeomPropagators/interface/Propagator.h>
#include <TrackingTools/Records/interface/TrackingComponentsRecord.h>
#include <TrackingTools/TransientTrack/interface/TransientTrack.h>
#include <MagneticField/Engine/interface/MagneticField.h>
#include <MagneticField/Records/interface/IdealMagneticFieldRecord.h>
#include <DataFormats/GeometrySurface/interface/BoundPlane.h>
#include <TrackingTools/TrajectoryState/interface/TrajectoryStateClosestToPoint.h>
#include <TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h>
#include <TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h>

using namespace edm;
using namespace std;
//
// class decleration
//

class TrackAnalyzerbyAssociator : public edm::EDAnalyzer {
    public:
        explicit TrackAnalyzerbyAssociator(const edm::ParameterSet&);
        ~TrackAnalyzerbyAssociator();


    private:
        virtual void beginJob();
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob();
        void getParticleSample(TrackingParticleCollection::size_type i);
        void associateTracks();
        void getTrackInformation(const reco::Track& recTrack);
        void validation(const std::vector< std::vector<SimHitIdpr> >& SamplematchIdsCollection, const unsigned int& SampleSelectedrecHitNumber);
        // ----------member data ---------------------------
        edm::InputTag recTrackTag_;
        edm::InputTag trackingParticleTag_;
        std::string thePropagatorName;
        bool useTracker;
        bool useMuon;
        double recTrackPurityThreshold;
        bool debug;

        edm::ESHandle<MagneticField> theMagneticField;
        edm::ESHandle<GlobalTrackingGeometry> theTrackingGeometryEH;
        edm::Handle<reco::TrackCollection> recTracksH;
        edm::Handle<TrackingParticleCollection> trackingParticlesH;
        edm::ParameterSet iConfiguration;

        std::string theRootFileName;
        TFile* theFile;
        TTree* ExTree;
        unsigned int EventNumber;
        unsigned int matchNumber;
        double recTrackrefMomentum;
        double recTrackrefPhi;
        double recTrackrefEta;
        double recTrackinnerMomentum;
        double recTrackinnerPhi;
        double recTrackinnerEta;
        unsigned int recTrackinnerValid;
        double recTrackouterMomentum;
        double recTrackouterPhi;
        double recTrackouterEta;
        unsigned int recTrackouterValid;
        double simTrackrefMomentum;
        double simTrackrefPhi;
        double simTrackrefEta;
        double simTrackinnerMomentum;
        double simTrackinnerPhi;
        double simTrackinnerEta;
        unsigned int simTrackinnerMatch;
        double simTrackouterMomentum;
        double simTrackouterPhi;
        double simTrackouterEta;
        unsigned int simTrackouterMatch;
        double recTrackinnerMomentumofTSOS;
        double recTrackinnerPhiofTSOS;
        double recTrackinnerEtaofTSOS;
        unsigned int recTrackinnerValidofTSOS;
        double recTrackouterMomentumofTSOS;
        double recTrackouterPhiofTSOS;
        double recTrackouterEtaofTSOS;
        unsigned int recTrackouterValidofTSOS;
        double recTrackimpactMomentumofTSOS;
        double recTrackimpactPhiofTSOS;
        double recTrackimpactEtaofTSOS;
        unsigned int recTrackimpactValidofTSOS;
        int recTrackCharge;
        double simTrackMomentumPt;
        double simTrackPhi;
        double simTrackEta;
        int simTrackCharge;

        Propagator* thePropagator;
        TrackerHitAssociator* TrackerTruth;
        DTHitAssociator* DTTruth;
        MuonTruth* CSCTruth;
        RPCHitAssociator* RPCTruth;

        bool interestingTrack;
        bool associatedrecTracksselected;
        TrackingParticleRef trackingParticleSampleR;
        GlobalVector refMomentum;
        GlobalPoint refPosition;
        unsigned int refDetId;
        GlobalVector innerrecMomentum;
        bool innerrecValid;
        GlobalVector outerrecMomentum;
        bool outerrecValid;
        std::vector< std::vector<SimHitIdpr> > TrackermatchIdsCollection;
        std::vector< std::vector<SimHitIdpr> > DTmatchIdsCollection;
        std::vector< std::vector<SimHitIdpr> > CSCmatchIdsCollection;
        std::vector< std::vector<SimHitIdpr> > RPCmatchIdsCollection;
        std::vector< std::vector<SimHitIdpr> > MuonmatchIdsCollection;
        std::vector< std::vector<SimHitIdpr> > TotalmatchIdsCollection;
        unsigned int TrackerValidrecHitNumber;
        unsigned int DTValidrecHitNumber;
        unsigned int CSCValidrecHitNumber;
        unsigned int RPCValidrecHitNumber;
        unsigned int MuonValidrecHitNumber;
        unsigned int TotalValidrecHitNumber;
        unsigned int innerDetId;
        unsigned int outerDetId;
        GlobalVector simMomentum;
        GlobalVector innersimMomentum;
        bool innersimFound;
        GlobalVector outersimMomentum;
        bool outersimFound;
        double Purity;
        reco::TransientTrack transientTrack;
        GlobalVector innerMomentumofTSOS;
        bool innerValidofTSOS;
        GlobalVector outerMomentumofTSOS;
        bool outerValidofTSOS;
        GlobalVector impactMomentumofTSOS;
        bool impactValidofTSOS;
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
TrackAnalyzerbyAssociator::TrackAnalyzerbyAssociator(const edm::ParameterSet& iConfig) {
    //now do what ever initialization is needed
    recTrackTag_ = iConfig.getUntrackedParameter<edm::InputTag>("recTrackTag");
    trackingParticleTag_ = iConfig.getUntrackedParameter<edm::InputTag>("trackingParticleTag");
    thePropagatorName = iConfig.getParameter<std::string>("PropagatorName");
    useTracker = iConfig.getUntrackedParameter<bool>("useTracker");
    useMuon = iConfig.getUntrackedParameter<bool>("useMuon");
    recTrackPurityThreshold = iConfig.getUntrackedParameter<double>("recTrackPurityThreshold");
    debug = iConfig.getUntrackedParameter<bool>("debug");

    iConfiguration = iConfig;

    theRootFileName = iConfig.getUntrackedParameter<std::string>("theRootFileName", "TrackValidation.root");
    theFile = new TFile(theRootFileName.c_str(), "recreate");
    theFile->cd();
    ExTree = new TTree("ExTree", "ExTree");
    //ExTree->Branch("trackingParticleMomentum", &simMomentum, "simMomentum/D");
    ExTree->Branch("EventNumber", &EventNumber, "EventNumber/i");
    ExTree->Branch("trackingParticleMatch", &matchNumber, "matchNumber/i");
    ExTree->Branch("recTrackPurity", &Purity, "Purity/D");
    ExTree->Branch("recTrackrefMomentum", &recTrackrefMomentum, "recTrackrefMomentum/D");
    ExTree->Branch("recTrackrefPhi", &recTrackrefPhi, "recTrackrefPhi/D");
    ExTree->Branch("recTrackrefEta", &recTrackrefEta, "recTrackrefEta/D");
    ExTree->Branch("recTrackinnerMomentum", &recTrackinnerMomentum, "recTrackinnerMomentum/D");
    ExTree->Branch("recTrackinnerPhi", &recTrackinnerPhi, "recTrackinnerPhi/D");
    ExTree->Branch("recTrackinnerEta", &recTrackinnerEta, "recTrackinnerEta/D");
    ExTree->Branch("recTrackinnerValid", &recTrackinnerValid, "recTrackinnerValid/i");
    ExTree->Branch("recTrackouterMomentum", &recTrackouterMomentum, "recTrackouterMomentum/D");
    ExTree->Branch("recTrackouterPhi", &recTrackouterPhi, "recTrackouterPhi/D");
    ExTree->Branch("recTrackouterEta", &recTrackouterEta, "recTrackouterEta/D");
    ExTree->Branch("recTrackouterValid", &recTrackouterValid, "recTrackouterValid/i");
    ExTree->Branch("simTrackrefMomentum", &simTrackrefMomentum, "simTrackrefMomentum/D");
    ExTree->Branch("simTrackrefPhi", &simTrackrefPhi, "simTrackrefPhi/D");
    ExTree->Branch("simTrackrefEta", &simTrackrefEta, "simTrackrefEta/D");
    ExTree->Branch("simTrackinnerMomentum", &simTrackinnerMomentum, "simTrackinnerMomentum/D");
    ExTree->Branch("simTrackinnerPhi", &simTrackinnerPhi, "simTrackinnerPhi/D");
    ExTree->Branch("simTrackinnerEta", &simTrackinnerEta, "simTrackinnerEta/D");
    ExTree->Branch("simTrackinnerMatch", &simTrackinnerMatch, "simTrackinnerMatch/i");
    ExTree->Branch("simTrackouterMomentum", &simTrackouterMomentum, "simTrackouterMomentum/D");
    ExTree->Branch("simTrackouterPhi", &simTrackouterPhi, "simTrackouterPhi/D");
    ExTree->Branch("simTrackouterEta", &simTrackouterEta, "simTrackouterEta/D");
    ExTree->Branch("simTrackouterMatch", &simTrackouterMatch, "simTrackouterMatch/i");
    ExTree->Branch("recTrackinnerMomentumofTSOS", &recTrackinnerMomentumofTSOS, "recTrackinnerMomentumofTSOS/D");
    ExTree->Branch("recTrackinnerPhiofTSOS", &recTrackinnerPhiofTSOS, "recTrackinnerPhiofTSOS/D");
    ExTree->Branch("recTrackinnerEtaofTSOS", &recTrackinnerEtaofTSOS, "recTrackinnerEtaofTSOS/D");
    ExTree->Branch("recTrackinnerValidofTSOS", &recTrackinnerValidofTSOS, "recTrackinnerValidofTSOS/i");
    ExTree->Branch("recTrackouterMomentumofTSOS", &recTrackouterMomentumofTSOS, "recTrackouterMomentumofTSOS/D");
    ExTree->Branch("recTrackouterPhiofTSOS", &recTrackouterPhiofTSOS, "recTrackouterPhiofTSOS/D");
    ExTree->Branch("recTrackouterEtaofTSOS", &recTrackouterEtaofTSOS, "recTrackouterEtaofTSOS/D");
    ExTree->Branch("recTrackouterValidofTSOS", &recTrackouterValidofTSOS, "recTrackouterValidofTSOS/i");
    ExTree->Branch("recTrackimpactMomentumofTSOS", &recTrackimpactMomentumofTSOS, "recTrackimpactMomentumofTSOS/D");
    ExTree->Branch("recTrackimpactPhiofTSOS", &recTrackimpactPhiofTSOS, "recTrackimpactPhiofTSOS/D");
    ExTree->Branch("recTrackimpactEtaofTSOS", &recTrackimpactEtaofTSOS, "recTrackimpactEtaofTSOS/D");
    ExTree->Branch("recTrackimpactValidofTSOS", &recTrackimpactValidofTSOS, "recTrackimpactValidofTSOS/i");
    ExTree->Branch("recTrackCharge", &recTrackCharge, "recTrackCharge/I");
    ExTree->Branch("simTrackMomentumPt", &simTrackMomentumPt, "simTrackMomentumPt/D");
    ExTree->Branch("simTrackPhi", &simTrackPhi, "simTrackPhi/D");
    ExTree->Branch("simTrackEta", &simTrackEta, "simTrackEta/D");
    ExTree->Branch("simTrackCharge", &simTrackCharge, "simTrackCharge/I");
}


TrackAnalyzerbyAssociator::~TrackAnalyzerbyAssociator() {
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void TrackAnalyzerbyAssociator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace std;
    // Get simulation inform
    iEvent.getByLabel(recTrackTag_, recTracksH);
    iEvent.getByLabel(trackingParticleTag_, trackingParticlesH);
    // Tracking Geometry
    iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometryEH);
    // Get the magnetic field for building transient track
    iSetup.get<IdealMagneticFieldRecord>().get(theMagneticField);
    // Get the propagator for building impact point from transient track
    ESHandle<Propagator> prop;
    iSetup.get<TrackingComponentsRecord>().get(thePropagatorName, prop);
    thePropagator = prop->clone();
    thePropagator->setPropagationDirection(alongMomentum);

    //TrackerTruth = new TrackerHitAssociator(iEvent);
    DTTruth = new DTHitAssociator(iEvent, iSetup, iConfiguration, false);
    CSCTruth = new MuonTruth(iEvent, iSetup, iConfiguration);
    RPCTruth = new RPCHitAssociator(iEvent, iSetup, iConfiguration);

    EventNumber = iEvent.id().event();
    // Loop on simTrack and get information
    for(TrackingParticleCollection::size_type i = 0; i < trackingParticlesH->size(); i++) {
        associatedrecTracksselected = false;
        getParticleSample(i);
        if(interestingTrack == true)
            associateTracks();
    }

    //delete TrackerTruth;
    delete thePropagator;
    delete DTTruth;
    delete CSCTruth;
    delete RPCTruth;
}

void TrackAnalyzerbyAssociator::getParticleSample(TrackingParticleCollection::size_type i) {
    trackingParticleSampleR = TrackingParticleRef(trackingParticlesH, i);
    const std::vector<SimTrack> simTracks = trackingParticleSampleR->g4Tracks();
    int thePDGId = trackingParticleSampleR->pdgId();
    simTrackCharge = 0; // we don't care other particles
    if(thePDGId == 13)
        simTrackCharge = -1;
    if(thePDGId == -13)
        simTrackCharge = 1;

    interestingTrack = false;
    if(abs(trackingParticleSampleR->pdgId()) == 13)
        for(std::vector<SimTrack>::const_iterator simTrackIter = simTracks.begin(); simTrackIter != simTracks.end(); simTrackIter++) {
            int TrackType = simTrackIter->type();
            if(debug) cout << "Get simTrack from trackingParticle. Type " << TrackType << ", ID: " << simTrackIter->trackId() << endl;
            if(abs(TrackType) == 13) {
                interestingTrack = true;
                simMomentum = GlobalVector(simTrackIter->momentum().Px(), simTrackIter->momentum().Py(), simTrackIter->momentum().Pz());
                simTrackMomentumPt = simMomentum.perp();
                simTrackPhi = simMomentum.phi().value();
                simTrackEta = simMomentum.eta();
            }
        }
}


void TrackAnalyzerbyAssociator::associateTracks() {

    if(associatedrecTracksselected == true)
        if(debug) cout << "!!!The recTracks have already been selected! Now will re-select them again" << endl;

    matchNumber = 0;
    for(reco::TrackCollection::const_iterator TrackIter = recTracksH->begin(); TrackIter != recTracksH->end(); TrackIter++) {
        //TrackermatchIdsCollection.clear();
        DTmatchIdsCollection.clear();
        CSCmatchIdsCollection.clear();
        RPCmatchIdsCollection.clear();
        MuonmatchIdsCollection.clear();
        TotalmatchIdsCollection.clear();
        TrackerValidrecHitNumber = 0;
        DTValidrecHitNumber = 0;
        CSCValidrecHitNumber = 0;
        RPCValidrecHitNumber = 0;
        MuonValidrecHitNumber = 0;
        TotalValidrecHitNumber = 0;

        getTrackInformation(*TrackIter);

        unsigned int SampleSelectedrecHitNumber = 0;
        std::vector< std::vector<SimHitIdpr> > SamplematchIdsCollection;
        if(useTracker == true && useMuon == true) {
            SamplematchIdsCollection = TotalmatchIdsCollection;
            SampleSelectedrecHitNumber = TotalValidrecHitNumber;
        }
        if(useTracker == true && useMuon == false) {
            SamplematchIdsCollection = TrackermatchIdsCollection;
            SampleSelectedrecHitNumber = TrackerValidrecHitNumber;
        }
        if(useTracker == false && useMuon == true) {
            SamplematchIdsCollection = MuonmatchIdsCollection;
            SampleSelectedrecHitNumber = MuonValidrecHitNumber;
        }

        validation(SamplematchIdsCollection, SampleSelectedrecHitNumber);

        if(debug) cout << "This trackingParticle has TrackerValidrecHitNumber " << TrackerValidrecHitNumber << ", DTValidrecHitNumber " << DTValidrecHitNumber << ", CSCValidrecHitNumber " << CSCValidrecHitNumber << ", RPCValidrecHitNumber " << RPCValidrecHitNumber << endl;
        if(Purity >= recTrackPurityThreshold) {
            if(debug) cout << "The reco Track is fully associated to this trackingParticle." << endl;
            recTrackrefMomentum = refMomentum.mag();
            recTrackrefPhi = refMomentum.phi().value();
            recTrackrefEta = refMomentum.eta();
            recTrackinnerMomentum = innerrecMomentum.mag();
            recTrackinnerPhi = innerrecMomentum.phi().value();
            recTrackinnerEta = innerrecMomentum.eta();
            recTrackinnerValid = innerrecValid;
            recTrackouterMomentum = outerrecMomentum.mag();
            recTrackouterPhi = outerrecMomentum.phi().value();
            recTrackouterEta = outerrecMomentum.eta();
            recTrackouterValid = outerrecValid;
            recTrackinnerMomentumofTSOS = innerMomentumofTSOS.mag();
            recTrackinnerPhiofTSOS = innerMomentumofTSOS.phi().value();
            recTrackinnerEtaofTSOS = innerMomentumofTSOS.eta();
            recTrackinnerValidofTSOS = innerValidofTSOS;
            recTrackouterMomentumofTSOS = outerMomentumofTSOS.mag();
            recTrackouterPhiofTSOS = outerMomentumofTSOS.phi().value();
            recTrackouterEtaofTSOS = outerMomentumofTSOS.eta();
            recTrackouterValidofTSOS = outerValidofTSOS;
            recTrackimpactMomentumofTSOS = impactMomentumofTSOS.mag();
            recTrackimpactPhiofTSOS = impactMomentumofTSOS.phi().value();
            recTrackimpactEtaofTSOS = impactMomentumofTSOS.eta();
            recTrackimpactValidofTSOS = impactValidofTSOS;
            recTrackCharge = TrackIter->charge();
            matchNumber++;
            //simTrackrefMomentum = 
            simTrackMomentumPt = simMomentum.perp();
            simTrackPhi = simMomentum.phi().value();
            simTrackEta = simMomentum.eta();    
            simTrackinnerMomentum = innersimMomentum.mag();
            simTrackinnerPhi = innersimMomentum.phi().value();
            simTrackinnerEta = innersimMomentum.eta();
            simTrackinnerMatch = innersimFound;
            simTrackouterMomentum = outersimMomentum.mag();
            simTrackouterPhi = outersimMomentum.phi().value();
            simTrackouterEta = outersimMomentum.eta();
            simTrackouterMatch = outersimFound;
            ExTree->Fill();
        }
        else {
            if(debug) cout << "The reco Track is not fully associated to this trackingParticle." << endl;
        }
    }
    if(matchNumber == 0) {
        recTrackrefMomentum = 0;
        recTrackrefPhi = 0;
        recTrackrefEta = 0;
        recTrackinnerMomentum = 0;
        recTrackinnerPhi = 0;
        recTrackinnerEta = 0;
        recTrackinnerValid = 0;
        recTrackouterMomentum = 0;
        recTrackouterPhi = 0;
        recTrackouterEta = 0;
        recTrackouterValid = 0;
        recTrackinnerMomentumofTSOS = 0;
        recTrackinnerPhiofTSOS = 0;
        recTrackinnerEtaofTSOS = 0;
        recTrackinnerValidofTSOS = 0;
        recTrackouterMomentumofTSOS = 0;
        recTrackouterPhiofTSOS = 0;
        recTrackouterEtaofTSOS = 0;
        recTrackouterValidofTSOS = 0;
        recTrackimpactMomentumofTSOS = 0;
        recTrackimpactPhiofTSOS = 0;
        recTrackimpactEtaofTSOS = 0;
        recTrackimpactValidofTSOS = 0;
        recTrackCharge = 0;
        matchNumber = 0;
        //simTrackMomentumPt = 0;
        //simTrackPhi = 0;
        //simTrackEta = 0;
        //simTrackCharge = 0;
        simTrackinnerMomentum = 0;
        simTrackinnerPhi = 0;
        simTrackinnerEta = 0;
        simTrackinnerMatch = 0;
        simTrackouterMomentum = 0;
        simTrackouterPhi = 0;
        simTrackouterEta = 0;
        simTrackouterMatch = 0;
        ExTree->Fill();
    }
    associatedrecTracksselected = true;
}

void TrackAnalyzerbyAssociator::getTrackInformation(const reco::Track& recTrack) {

    // matchId and ref points
    TrackermatchIdsCollection.clear();
    DTmatchIdsCollection.clear();
    CSCmatchIdsCollection.clear();
    RPCmatchIdsCollection.clear();
    MuonmatchIdsCollection.clear();
    TotalmatchIdsCollection.clear();
    TrackerValidrecHitNumber = 0;
    DTValidrecHitNumber = 0;
    CSCValidrecHitNumber = 0;
    RPCValidrecHitNumber = 0;
    MuonValidrecHitNumber = 0;
    TotalValidrecHitNumber = 0;

    std::vector<SimHitIdpr> simHitPairs;
    const math::XYZPoint innerPosition = recTrack.innerPosition();
    const math::XYZPoint outerPosition = recTrack.outerPosition();

    innerDetId = recTrack.innerDetId();
    outerDetId = recTrack.outerDetId();

    refMomentum = GlobalVector(recTrack.momentum().x(), recTrack.momentum().y(), recTrack.momentum().z());
    refPosition = GlobalPoint(recTrack.referencePoint().x(), recTrack.referencePoint().y(), recTrack.referencePoint().z());

    innerrecMomentum = GlobalVector(recTrack.extra()->innerMomentum().x(), recTrack.extra()->innerMomentum().y(), recTrack.extra()->innerMomentum().z());
    innerrecValid = recTrack.innerOk();
    outerrecMomentum = GlobalVector(recTrack.extra()->outerMomentum().x(), recTrack.extra()->outerMomentum().y(), recTrack.extra()->outerMomentum().z());
    outerrecValid = recTrack.outerOk();

    if(debug) cout << "innerPosition: " << innerPosition << ", on DetId: " << innerDetId << ", is OK: " << innerrecValid << ", momentum: " << innerrecMomentum << ". outerPosition: " << outerPosition << ", on DetId: " << outerDetId << ", is OK: " << outerrecValid << ", momentum: " << outerrecMomentum << "refPosition: " << refPosition << ", refMomentum: " << refMomentum << ". With recHit size: " << recTrack.recHitsSize() << endl;

    unsigned int i = 0;
    for(trackingRecHit_iterator recHitIter = recTrack.recHitsBegin(); recHitIter != recTrack.recHitsEnd(); recHitIter++) {
        i++;
        if(debug) cout << "Now checking " << i << "th recHit. isValid is " << (*recHitIter)->isValid() << endl;
        if((*recHitIter)->isValid() == false)
            continue;

        simHitPairs.clear();
        DetId recHitDetId = (*recHitIter)->geographicalId();
        DetId::Detector Det = recHitDetId.det();
        GeomDetEnumerators::SubDetector subDet = theTrackingGeometryEH->idToDet(recHitDetId)->subDetector();
        LocalPoint localrecHitPosition = (*recHitIter)->localPosition();
        GlobalPoint globalrecHitPosition = theTrackingGeometryEH->idToDet(recHitDetId)->toGlobal(localrecHitPosition);
        if(debug) cout << "localrecHitPosition: " << localrecHitPosition << ". globalrecHitPosition: " << globalrecHitPosition << endl;
        /*
           if(Det == DetId::Tracker) {
           simHitPairs = TrackerTruth->associateHitId(*(*recHitIter));
           TrackerValidrecHitNumber++;
           TrackermatchIdsCollection.push_back(simHitPairs);
           }
         */
        if(subDet == GeomDetEnumerators::DT) {
            if(debug) cout << "recHit is from DT." << endl;
            simHitPairs = DTTruth->associateHitId(*(*recHitIter));
            DTValidrecHitNumber++;
            DTmatchIdsCollection.push_back(simHitPairs);
            MuonmatchIdsCollection.push_back(simHitPairs);
        }
        if(subDet == GeomDetEnumerators::CSC) {
            if(debug) cout << "recHit is from CSC." << endl;
            simHitPairs = CSCTruth->associateHitId(*(*recHitIter));
            CSCValidrecHitNumber++;
            CSCmatchIdsCollection.push_back(simHitPairs);
            MuonmatchIdsCollection.push_back(simHitPairs);
        }
        if(subDet == GeomDetEnumerators::RPCBarrel || subDet == GeomDetEnumerators::RPCEndcap) {
            simHitPairs = RPCTruth->associateRecHit(*(*recHitIter));
            if(debug) cout << "recHit is from RPC, link to simHit pairs size: " << simHitPairs.size() << endl;
            RPCValidrecHitNumber++;
            RPCmatchIdsCollection.push_back(simHitPairs);
            MuonmatchIdsCollection.push_back(simHitPairs);
        }
        TotalmatchIdsCollection.push_back(simHitPairs);
        // inner and outer point seems not settled on the recHits, may be they are updated hit points and have residual with the original rechits
        /*
           if(globalrecHitPosition.x() == innerPosition.x() && globalrecHitPosition.y() == innerPosition.y() && globalrecHitPosition.z() == innerPosition.z()) {
           innermactchedIds = simHitPairs;
           innerDetId = recHitDetId;
           if(debug) cout << "Find inner mactched recHit. DetId: " << innerDetId << endl;
           }
           if(globalrecHitPosition.x() == outerPosition.x() && globalrecHitPosition.y() == outerPosition.y() && globalrecHitPosition.z() == outerPosition.z()) {
           outermactchedIds = simHitPairs;
           outerDetId = recHitDetId;
           if(debug) cout << "Find outer mactched recHit. DetId: " << outerDetId << endl;
           }
         */
    }
    MuonValidrecHitNumber = DTValidrecHitNumber + CSCValidrecHitNumber + RPCValidrecHitNumber;
    TotalValidrecHitNumber = MuonValidrecHitNumber + TrackerValidrecHitNumber;

    // By transient track
    transientTrack = reco::TransientTrack(recTrack, &*theMagneticField, theTrackingGeometryEH);
    // Get the trajectory state near the track's reference point? 
    TrajectoryStateClosestToPoint refTSCP = transientTrack.impactPointTSCP();
    // Set this trajectory state to free, that is the state only in global frame
    const FreeTrajectoryState& startFTS = refTSCP.theState();
    // Use the propagator to get the state on RPC surface
    const GeomDet* innerDet = theTrackingGeometryEH->idToDet((DetId)innerDetId);
    const BoundPlane& innerBoundPlane = innerDet->surface();
    const GeomDet* outerDet = theTrackingGeometryEH->idToDet((DetId)outerDetId);
    const BoundPlane& outerBoundPlane = outerDet->surface();
    TrajectoryStateOnSurface impactTSOS = transientTrack.impactPointState();
    TrajectoryStateOnSurface innerTSOS = thePropagator->propagate(startFTS, innerBoundPlane);
    TrajectoryStateOnSurface outerTSOS = thePropagator->propagate(startFTS, outerBoundPlane);
    // Must check the valid status of TSOS, or it may get segment fail while running
    if(impactTSOS.isValid()) {
        impactMomentumofTSOS = impactTSOS.globalMomentum();
        impactValidofTSOS = true;    
    }
    else {
        impactMomentumofTSOS = GlobalVector(0, 0, 0);
        impactValidofTSOS = false;
    }
    if(innerTSOS.isValid()) {
        innerMomentumofTSOS = innerTSOS.globalMomentum();
        innerValidofTSOS = true;
    }
    else {
        innerMomentumofTSOS = GlobalVector(0, 0, 0);
        innerValidofTSOS = false;
    }
    if(outerTSOS.isValid()) {
        outerMomentumofTSOS = outerTSOS.globalMomentum();
        outerValidofTSOS = true;
    }
    else {
        outerMomentumofTSOS = GlobalVector(0, 0, 0);
        outerValidofTSOS = false;
    }
    if(debug) cout << "impactMmentumofTSOS: " << impactMomentumofTSOS << ", impactValidofTSOS: " << impactValidofTSOS << ". innerMomentumofTSOS: " << innerMomentumofTSOS << ", innerValidofTSOS: " << innerValidofTSOS << ". outerMomentumofTSOS: " << outerMomentumofTSOS << ", outerValidofTSOS: " << outerValidofTSOS << endl;
}

void TrackAnalyzerbyAssociator::validation(const std::vector< std::vector<SimHitIdpr> >& SamplematchIdsCollection, const unsigned int& SampleSelectedrecHitNumber) {
    unsigned int SampleMatchedrecHitNumber = 0;
    for(std::vector< std::vector<SimHitIdpr> >::const_iterator SamplematchIds = SamplematchIdsCollection.begin(); SamplematchIds != SamplematchIdsCollection.end(); SamplematchIds++) {
        bool matchthisHit = false;
        for(std::vector<SimHitIdpr>::const_iterator matchIdIter = SamplematchIds->begin(); matchIdIter != SamplematchIds->end(); matchIdIter++) {
            for(TrackingParticle::g4t_iterator simTrackIter = trackingParticleSampleR->g4Track_begin(); simTrackIter !=  trackingParticleSampleR->g4Track_end(); ++simTrackIter) {
                int simPID = simTrackIter->type();
                bool GeneratorPart = simTrackIter->noGenpart();
                unsigned int simTrackId = simTrackIter->trackId();
                if(debug) cout << "Sim track PID: " << simPID << ". GeneratorPart: " << GeneratorPart << ". simTrackId: " << simTrackId << endl;
                if(abs(simPID) != 13)
                    continue;

                if(matchIdIter->first == simTrackIter->trackId() && matchIdIter->second == simTrackIter->eventId()) 
                    matchthisHit = true;
            }
        }
        if(matchthisHit == true)
            SampleMatchedrecHitNumber++;
        if(debug) cout << "For this recHit's simHitIdPair matchthisHit is " << matchthisHit << endl;
    }
    Purity = (double)SampleMatchedrecHitNumber / (double)SampleSelectedrecHitNumber;

    // By inner and outer information 
    innersimMomentum = GlobalVector(0, 0, 0);
    outersimMomentum = GlobalVector(0, 0, 0);
    innersimFound = false;
    outersimFound = false;
    for(TrackingParticle::g4t_iterator simTrackIter = trackingParticleSampleR->g4Track_begin(); simTrackIter !=  trackingParticleSampleR->g4Track_end(); ++simTrackIter) {
        int simPID = simTrackIter->type();
        bool GeneratorPart = simTrackIter->noGenpart();
        unsigned int simTrackId = simTrackIter->trackId();
        if(debug) cout << "Sim track PID: " << simPID << ". GeneratorPart: " << GeneratorPart << ". simTrackId: " << simTrackId << endl;
        if(abs(simPID) != 13)
            continue;

        simMomentum = GlobalVector(simTrackIter->momentum().Px(), simTrackIter->momentum().Py(), simTrackIter->momentum().Pz());
        for(std::vector<PSimHit>::const_iterator simHitIter = trackingParticleSampleR->pSimHit_begin(); simHitIter != trackingParticleSampleR->pSimHit_end(); simHitIter++) {
            unsigned int simHitTrackId = simHitIter->trackId();
            EncodedEventId simHitEventId = simHitIter->eventId();
            unsigned int simHitDetUnitId = simHitIter->detUnitId();
            if(debug) cout << "simHit for TrackId " << simHitTrackId << ", on DetId " << simHitDetUnitId << endl;
            if(simHitDetUnitId == innerDetId && simHitTrackId == simTrackId) {
                LocalVector localMomentum = simHitIter->momentumAtEntry();
                GlobalVector globalMomentum = theTrackingGeometryEH->idToDetUnit((DetId)simHitDetUnitId)->toGlobal(localMomentum);
                if(debug) cout << "Inner point sim momentum is " << globalMomentum << endl;
                innersimMomentum = globalMomentum;
                innersimFound = true;
            }
            if(simHitDetUnitId == outerDetId && simHitTrackId == simTrackId) {
                LocalVector localMomentum = simHitIter->momentumAtEntry();
                GlobalVector globalMomentum = theTrackingGeometryEH->idToDetUnit((DetId)simHitDetUnitId)->toGlobal(localMomentum);
                if(debug) cout << "Outer point sim momentum is " << globalMomentum << endl;
                outersimMomentum = globalMomentum;
                outersimFound = true;
            }
        }
    }
}
// ------------ method called once each job just before starting event loop  ------------
void TrackAnalyzerbyAssociator::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void TrackAnalyzerbyAssociator::endJob() {
    theFile->Write();
    theFile->Close();
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackAnalyzerbyAssociator);
