// -*- C++ -*-
//
// Package:    RPCSeedProducer
// Class:      RPCSeedProducer
// 
/**\class RPCSeedProducer RPCSeedProducer.cc MyModule/RPCSeedProducer/src/RPCSeedProducer.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Haiyun Teng
//         Created:  Wed Oct 29 17:24:36 CET 2008
// $Id: RPCSeedProducer.cc,v 1.4 2012/03/25 18:28:14 hyteng Exp $
//
//


// system include files
#include <memory>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
// special include files
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include <vector>
// Using other classes
#include "MyModule/RPCSeedProducer/src/DebugSignal.h"
#include "MyModule/RPCSeedProducer/src/RPCSeedData.h"
#include "MyModule/RPCSeedProducer/src/RPCSeedPattern.h"
#include "MyModule/RPCSeedProducer/src/RPCSeedFinder.h"
#include "MyModule/RPCSeedProducer/src/RPCSeedrecHitFinder.h"
#include "MyModule/RPCSeedProducer/src/RPCCosmicSeedrecHitFinder.h"
#include "MyModule/RPCSeedProducer/src/RPCSeedLayerFinder.h"
#include "MyModule/RPCSeedProducer/src/RPCSeedOverlapper.h"
// Geometry
#include "TrackingTools/DetLayers/interface/DetLayer.h"
#include "RecoMuon/MeasurementDet/interface/MuonDetLayerMeasurements.h"
#include "RecoMuon/DetLayers/interface/MuonDetLayerGeometry.h"
#include "RecoMuon/Records/interface/MuonRecoGeometryRecord.h"
// Framework
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
// Math
#include "math.h"
// C++
#include <vector>

//
// constants, enums and typedefs
//
using namespace std;
using namespace edm;

typedef MuonTransientTrackingRecHit::MuonRecHitPointer MuonRecHitPointer;
typedef MuonTransientTrackingRecHit::ConstMuonRecHitPointer ConstMuonRecHitPointer;
typedef MuonTransientTrackingRecHit::MuonRecHitContainer MuonRecHitContainer;
typedef MuonTransientTrackingRecHit::ConstMuonRecHitContainer ConstMuonRecHitContainer;
typedef RPCSeedPattern::WeightedTrajectorySeed WeightedTrajectorySeed;

//
// class decleration
//

class RPCSeedFinder;

class RPCSeedProducer : public edm::EDProducer {
    public:
        explicit RPCSeedProducer(const edm::ParameterSet& iConfig);
        ~RPCSeedProducer();

    private:
        virtual void beginJob();
        virtual void beginRun(const edm::Run&, const edm::EventSetup& iSetup);
        virtual void produce(edm::Event& iEvent, const edm::EventSetup& iSetup);
        virtual void endJob();

        // ----------member data ---------------------------
        RPCSeedFinder Finder;
        RPCSeedrecHitFinder recHitFinder;
        RPCCosmicSeedrecHitFinder CosmicrecHitFinder;
        RPCSeedLayerFinder LayerFinder;
        RPCSeedOverlapper Overlapper;
        std::vector<WeightedTrajectorySeed> CandidateWeightedSeeds;
        std::vector<WeightedTrajectorySeed> GoodWeightedSeeds;
        edm::InputTag RPCrecHitTag_;
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
RPCSeedProducer::RPCSeedProducer(const edm::ParameterSet& iConfig) {

    // Configure other modules
    Finder.configure(iConfig);
    recHitFinder.configure(iConfig);
    CosmicrecHitFinder.configure(iConfig);
    LayerFinder.configure(iConfig);
    Overlapper.configure(iConfig);
    // Register the production
    produces<TrajectorySeedCollection>("GoodSeeds");
    produces<TrajectorySeedCollection>("CandidateSeeds");
    // Get event data Tag
    RPCrecHitTag_ = iConfig.getParameter<edm::InputTag>("RPCRecHitsLabel");
    
    if(debug) cout << endl << "[RPCSeedProducer] --> Constructor called" << endl;
}


RPCSeedProducer::~RPCSeedProducer() {
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    if(debug) cout << "[RPCSeedProducer] --> Destructor called" << endl;
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void RPCSeedProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    /* This is an event example
    //Read 'ExampleData' from the Event
    Handle<ExampleData> pIn;
    iEvent.getByLabel("example",pIn);

    //Use the ExampleData to create an ExampleData2 which 
    // is put into the Event
    std::auto_ptr<ExampleData2> pOut(new ExampleData2(*pIn));
    iEvent.put(pOut);
    */

    /* this is an EventSetup example
    //Read SetupData from the SetupRecord in the EventSetup
    ESHandle<SetupData> pSetup;
    iSetup.get<SetupRecord>().get(pSetup);
    */

    if(debug) cout << "test0" << endl;

    CosmicrecHitFinder.setEdge(iSetup);
    Overlapper.setEventSetup(iSetup);

    if(debug) cout << "test1" << endl;
    // clear Weighted Seeds from last reconstruction  
    GoodWeightedSeeds.clear();
    CandidateWeightedSeeds.clear();
    if(debug) cout << "test2" << endl;

    // Create the pointer to the Seed container
    auto_ptr<TrajectorySeedCollection> GoodCollection(new TrajectorySeedCollection());
    auto_ptr<TrajectorySeedCollection> CandidateCollection(new TrajectorySeedCollection());

    // Muon Geometry - DT, CSC and RPC 
    edm::ESHandle<MuonDetLayerGeometry> muonLayers;
    iSetup.get<MuonRecoGeometryRecord>().get(muonLayers);

    // Get the RPC layers
    vector<DetLayer*> RPCBarrelLayers = muonLayers->barrelRPCLayers();
    const DetLayer* RB4L  = RPCBarrelLayers[5];
    const DetLayer* RB3L  = RPCBarrelLayers[4];
    const DetLayer* RB22L = RPCBarrelLayers[3];
    const DetLayer* RB21L = RPCBarrelLayers[2];
    const DetLayer* RB12L = RPCBarrelLayers[1];
    const DetLayer* RB11L = RPCBarrelLayers[0];
    vector<DetLayer*> RPCEndcapLayers = muonLayers->endcapRPCLayers();
    const DetLayer* REM3L = RPCEndcapLayers[0];
    const DetLayer* REM2L = RPCEndcapLayers[1];
    const DetLayer* REM1L = RPCEndcapLayers[2];
    const DetLayer* REP1L = RPCEndcapLayers[3];
    const DetLayer* REP2L = RPCEndcapLayers[4];
    const DetLayer* REP3L = RPCEndcapLayers[5];

    // Get RPC recHits by MuonDetLayerMeasurements, while CSC and DT is set to false and with empty InputTag
    MuonDetLayerMeasurements muonMeasurements(edm::InputTag(), edm::InputTag(), RPCrecHitTag_, false, false, true);

    // Dispatch RPC recHits to the corresponding DetLayer, 6 layers for barrel and 3 layers for each endcap
    MuonRecHitContainer recHitsRPC[RPCLayerNumber];
    recHitsRPC[0] = muonMeasurements.recHits(RB11L, iEvent);
    recHitsRPC[1] = muonMeasurements.recHits(RB12L, iEvent);
    recHitsRPC[2] = muonMeasurements.recHits(RB21L, iEvent);
    recHitsRPC[3] = muonMeasurements.recHits(RB22L, iEvent);
    recHitsRPC[4] = muonMeasurements.recHits(RB3L, iEvent);
    recHitsRPC[5] = muonMeasurements.recHits(RB4L, iEvent); 
    recHitsRPC[6] = muonMeasurements.recHits(REM1L, iEvent);
    recHitsRPC[7] = muonMeasurements.recHits(REM2L, iEvent);
    recHitsRPC[8] = muonMeasurements.recHits(REM3L, iEvent);
    recHitsRPC[9] = muonMeasurements.recHits(REP1L, iEvent);
    recHitsRPC[10] = muonMeasurements.recHits(REP2L, iEvent);
    recHitsRPC[11] = muonMeasurements.recHits(REP3L, iEvent);

    // Print the size of recHits in each DetLayer
    if(debug) cout << "RB1in "  << recHitsRPC[0].size()  << " recHits" << endl;
    if(debug) cout << "RB1out " << recHitsRPC[1].size()  << " recHits" << endl;
    if(debug) cout << "RB2in "  << recHitsRPC[2].size()  << " recHits" << endl;
    if(debug) cout << "RB2out " << recHitsRPC[3].size()  << " recHits" << endl;
    if(debug) cout << "RB3 "    << recHitsRPC[4].size()  << " recHits" << endl;
    if(debug) cout << "RB4 "    << recHitsRPC[5].size()  << " recHits" << endl;
    if(debug) cout << "REM1 "   << recHitsRPC[6].size()  << " recHits" << endl;
    if(debug) cout << "REM2 "   << recHitsRPC[7].size()  << " recHits" << endl;
    if(debug) cout << "REM3 "   << recHitsRPC[8].size()  << " recHits" << endl;
    if(debug) cout << "REP1 "   << recHitsRPC[9].size()  << " recHits" << endl;
    if(debug) cout << "REP2 "   << recHitsRPC[10].size() << " recHits" << endl;
    if(debug) cout << "REP3 "   << recHitsRPC[11].size() << " recHits" << endl;

    // Set Input of RPCSeedFinder, PCSeedrecHitFinder, CosmicrecHitFinder, RPCSeedLayerFinder
    recHitFinder.setInput(recHitsRPC);
    CosmicrecHitFinder.setInput(recHitsRPC);
    LayerFinder.setInput(recHitsRPC);
    
    // Set Magnetic Field EventSetup of RPCSeedFinder
    Finder.setEventSetup(iSetup);
    // Set SimData
    Finder.setSimData(iEvent, iSetup);

    // Start from filling layers to filling seeds
    LayerFinder.fill();
    Overlapper.run();

    // Save seeds to event
    for(vector<WeightedTrajectorySeed>::iterator Weightedseed = GoodWeightedSeeds.begin(); Weightedseed != GoodWeightedSeeds.end(); ++Weightedseed)
        GoodCollection->push_back((*Weightedseed).first);
    for(vector<WeightedTrajectorySeed>::iterator Weightedseed = CandidateWeightedSeeds.begin(); Weightedseed != CandidateWeightedSeeds.end(); ++Weightedseed)
        CandidateCollection->push_back((*Weightedseed).first);

    // Put the seed to event
    if(debug) cout << "Putting the seeds in to branch." << endl;
    iEvent.put(GoodCollection, "GoodSeeds");
    iEvent.put(CandidateCollection, "CandidateSeeds");

    // Unset the input of RPCSeedFinder, PCSeedrecHitFinder, RPCSeedLayerFinder
    if(debug) cout << "Unset the input of sub module." << endl;
    Finder.unsetSimData();
    recHitFinder.unsetInput();
    CosmicrecHitFinder.unsetInput();
    LayerFinder.unsetInput();
    if(debug) cout << "Event done." << endl;
}

void RPCSeedProducer::beginJob() {

    // Set link and EventSetup of RPCSeedFinder, PCSeedrecHitFinder, CosmicrecHitFinder, RPCSeedLayerFinder
    if(debug) cout << "Set link of sub modules." << endl;
    
    Finder.setOutput(&GoodWeightedSeeds, &CandidateWeightedSeeds);
    recHitFinder.setOutput(&Finder);
    CosmicrecHitFinder.setOutput(&Finder);
    LayerFinder.setOutput(&recHitFinder, &CosmicrecHitFinder);
    Overlapper.setIO(&GoodWeightedSeeds, &CandidateWeightedSeeds);
}

void RPCSeedProducer::beginRun(const edm::Run&, const edm::EventSetup& iSetup) {
    // only call at the beginning of the whole program running
}

void RPCSeedProducer::endJob() {
    if(debug) cout << "All jobs completed." << endl;
}

//define this as a plug-in
DEFINE_FWK_MODULE(RPCSeedProducer);
