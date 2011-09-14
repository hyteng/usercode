/**
 *  See header file for a description of this class.
 *
 */


#include "MyModule/RPCSeedProducer/src/RPCSeedFinder.h"
#include "MyModule/RPCSeedProducer/src/DebugSignal.h"
#include <iomanip>

using namespace std;
using namespace edm;


RPCSeedFinder::RPCSeedFinder() {

    // Initiate the member
    isRecHitset = false;
    isConfigured = false;
    isOutputset = false;
    isEventSetupset = false;
    OneSeed.clear();
    SimOneSeed.clear();
}

RPCSeedFinder::~RPCSeedFinder() {

}

void RPCSeedFinder::configure(const edm::ParameterSet& iConfig) {

    useSimData = iConfig.getParameter<bool>("useSimData");
    OneSeed.configure(iConfig);
    SimOneSeed.configure(iConfig);
    isConfigured = true;
}

void RPCSeedFinder::setSimData(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    SimOneSeed.setSimData(iEvent, iSetup);
}

void RPCSeedFinder::unsetSimData() {
    SimOneSeed.unsetSimData();
}

void RPCSeedFinder::setOutput(std::vector<WeightedTrajectorySeed> *GoodWeightedRef, std::vector<WeightedTrajectorySeed> *CandidateWeightedRef) {

    GoodWeightedSeedsRef = GoodWeightedRef;
    CandidateWeightedSeedsRef = CandidateWeightedRef;
    isOutputset = true;
}

void RPCSeedFinder::clear() {
    OneSeed.clear();
    SimOneSeed.clear();
}

void RPCSeedFinder::setRecHits(ConstMuonRecHitContainer &RecHits) {
    OneSeed.setRecHits(RecHits);
    SimOneSeed.setRecHits(RecHits);
    isRecHitset = true;
}

/*
void RPCSeedFinder::setLayers(const std::vector<unsigned int>& Layers) {
    OneSeed.setlayer(Layers);
    isLayerset = true;
}
*/

void RPCSeedFinder::setEventSetup(const edm::EventSetup& iSetup) {
    eSetup = &iSetup;
    isEventSetupset = true;
}

void RPCSeedFinder::seed() {

    if(debug) cout << "[RPCSeedFinder] --> seeds called" << endl;

    if(isRecHitset == false || isOutputset == false || isConfigured == false || isEventSetupset == false) {
        if(debug) cout << "SeedFinder has not set Configuration or IO: " << isRecHitset << ", " << isOutputset << ", " << isConfigured << ", " << isEventSetupset << endl;
        return;
    }
    
    WeightedTrajectorySeed theWeightedSeed;
    int isGoodSeed = 0;
    const edm::EventSetup &iSetup = *eSetup;
    if(useSimData == false)
        theWeightedSeed = OneSeed.seed(iSetup, isGoodSeed);
    else
        theWeightedSeed = SimOneSeed.seed(iSetup, isGoodSeed);

    // Push back the Good seed
    if(debug) cout << "the Seed is Good or not: " << isGoodSeed << endl;
    if(isGoodSeed == 1) {
        GoodWeightedSeedsRef->push_back(theWeightedSeed);
    }
    // Push back cadidate seed but not the fake seed
    if(isGoodSeed >= 0) {
        CandidateWeightedSeedsRef->push_back(theWeightedSeed);
    }

    // Unset the signal
    OneSeed.clear();
    SimOneSeed.clear();
    isRecHitset = false;
}
