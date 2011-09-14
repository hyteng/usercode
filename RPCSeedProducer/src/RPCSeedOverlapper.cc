/**
 *  See header file for a description of this class.
 *
 */


#include "MyModule/RPCSeedProducer/src/RPCSeedOverlapper.h"
#include "MyModule/RPCSeedProducer/src/DebugSignal.h"
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <Geometry/CommonDetUnit/interface/GeomDetUnit.h>

using namespace std;
using namespace edm;

RPCSeedOverlapper::RPCSeedOverlapper() {

    isConfigured = false; 
    isIOset = false;
    isEventSetupset = false;
}

RPCSeedOverlapper::~RPCSeedOverlapper() {

}

void RPCSeedOverlapper::configure(const edm::ParameterSet& iConfig) {

    isCheckGoodOverlap = iConfig.getParameter<bool>("isCheckGoodOverlap");
    isCheckCandidateOverlap = iConfig.getParameter<bool>("isCheckCandidateOverlap");
    ShareRecHitsNumberThreshold = iConfig.getParameter<unsigned int>("ShareRecHitsNumberThreshold");
    isConfigured = true;
}

void RPCSeedOverlapper::setIO(std::vector<WeightedTrajectorySeed> *GoodWeightedRef, std::vector<WeightedTrajectorySeed> *CandidateWeightedRef) {
    if(debug) cout << "set RPCSeedOverlapper IO." << endl;
    GoodWeightedSeedsRef = GoodWeightedRef;
    CandidateWeightedSeedsRef = CandidateWeightedRef;
    isIOset = true;
}

void RPCSeedOverlapper::unsetIO() {
    if(debug) cout << "unset RPCSeedOverlapper IO." << endl;
    isIOset = false;
}

void RPCSeedOverlapper::setEventSetup(const edm::EventSetup& iSetup) {

    eSetup = &iSetup;
    isEventSetupset = true;
}

void RPCSeedOverlapper::run() {

    if(isConfigured == false || isIOset == false || isEventSetupset == false) {
        if(debug) cout << "RPCSeedOverlapper has not set Configuration or IO: " << isConfigured << ", " << isIOset << ", " << isEventSetupset << endl;
        return;
    }
    if(isCheckGoodOverlap == true)
        CheckOverlap(*eSetup, GoodWeightedSeedsRef);
    if(isCheckCandidateOverlap == true)
        CheckOverlap(*eSetup, CandidateWeightedSeedsRef);
}

void RPCSeedOverlapper::CheckOverlap(const edm::EventSetup& iSetup, std::vector<WeightedTrajectorySeed> *WeightedSeedsRef) {

    std::vector<WeightedTrajectorySeed> SortWeightedSeeds;
    std::vector<WeightedTrajectorySeed> TempWeightedSeeds;
    edm::OwnVector<TrackingRecHit> TempRecHits;

    edm::ESHandle<RPCGeometry> rpcGeometry;
    iSetup.get<MuonGeometryRecord>().get(rpcGeometry);

    while(WeightedSeedsRef->size() != 0) {
        if(debug) cout << "Finding the Weighted seeds group from " << WeightedSeedsRef->size() << " seeds which share some recHits" << endl; 
        // Take 1st seed in SeedsRef as referrence and find a collection which always share some recHits with some other
        TempRecHits.clear();
        TempWeightedSeeds.clear();
        int N = 0;
        for(vector<WeightedTrajectorySeed>::iterator itWeightedseed = WeightedSeedsRef->begin(); itWeightedseed != WeightedSeedsRef->end(); N++) {
            TrajectorySeed::range RecHitsRange = itWeightedseed->first.recHits();
            if(N == 0) {
                if(debug) cout << "Always take the 1st Weighted seed to be the referrence." << endl;
                for(TrajectorySeed::const_iterator it = RecHitsRange.first; it != RecHitsRange.second; it++) {
                    if(debug) cout << "Put its recHits to TempRecHits" << endl;
                    TempRecHits.push_back(it->clone());
                }
                if(debug) cout << "Put it to TempWeightedSeeds" << endl;
                TempWeightedSeeds.push_back(*itWeightedseed);
                if(debug) cout << "Then erase from WeightedSeedsRef->" << endl;
                itWeightedseed = WeightedSeedsRef->erase(itWeightedseed);
            }
            else {
                if(debug) cout << "Come to other Weighted seed for checking " << itWeightedseed->first.nHits() << " recHits from " << TempRecHits.size() << " Temp recHits" << endl;
                unsigned int ShareRecHitsNumber = 0;
                for(TrajectorySeed::const_iterator it = RecHitsRange.first; it != RecHitsRange.second; it++) {
                    if(isShareHit(TempRecHits, *it, rpcGeometry))
                        ShareRecHitsNumber++;
                }
                if(ShareRecHitsNumber >= ShareRecHitsNumberThreshold) {
                    if(debug) cout <<"This seed is found to belong to current share group" << endl;
                    for(TrajectorySeed::const_iterator it = RecHitsRange.first; it != RecHitsRange.second; it++) {
                        if(!isShareHit(TempRecHits, *it, rpcGeometry)) {
                            if(debug) cout << "Put its extra recHits to TempRecHits" << endl;
                            TempRecHits.push_back(it->clone());
                        }
                    }
                    if(debug) cout << "Put it to TempSeeds" << endl;
                    TempWeightedSeeds.push_back(*itWeightedseed);
                    if(debug) cout << "Then erase from SeedsRef" << endl;
                    itWeightedseed = WeightedSeedsRef->erase(itWeightedseed);
                }
                else
                    itWeightedseed++;
            }
        }
        // Find the Best Weighted seed and kick out those share recHits with it
        // The Best Weighted seed save in SortWeightedSeeds, those don't share recHits with it will be push back to WeightedSeedsRef for next while loop
        WeightedTrajectorySeed BestWeightedSeed;
        vector<WeightedTrajectorySeed>::iterator BestWeightediter;
        // Find the min Spt wrt Pt as the Best Seed
        double Quality = 1000000;
        unsigned NumberofHits = 0;
        if(debug) cout << "Find " << TempWeightedSeeds.size() << " seeds into one trajectory group" << endl;
        for(vector<WeightedTrajectorySeed>::iterator itWeightedseed = TempWeightedSeeds.begin(); itWeightedseed != TempWeightedSeeds.end(); itWeightedseed++) {
            unsigned int nHits = itWeightedseed->first.nHits();
            //std::vector<float> seed_error = itWeightedseed->first.startingState().errorMatrix();
            //double Spt = seed_error[1];
            double WeightedQuality = itWeightedseed->second;
            if(debug) cout << "Find a Weighted seed with quality " << WeightedQuality << endl;
            if((NumberofHits < nHits) || (NumberofHits == nHits && WeightedQuality < Quality)) {
                NumberofHits = nHits;
                Quality = WeightedQuality;
                BestWeightedSeed = *itWeightedseed;
                BestWeightediter = itWeightedseed;
            }
        }
        if(debug) cout << "Best Good Temp seed's quality is " << Quality <<endl;
        SortWeightedSeeds.push_back(BestWeightedSeed);
        TempWeightedSeeds.erase(BestWeightediter);
        TempRecHits.clear();

        for(TrajectorySeed::const_iterator it = BestWeightedSeed.first.recHits().first; it != BestWeightedSeed.first.recHits().second; it++)
            TempRecHits.push_back(it->clone());

        for(vector<WeightedTrajectorySeed>::iterator itWeightedseed = TempWeightedSeeds.begin(); itWeightedseed != TempWeightedSeeds.end(); ) {
            if(debug) cout << "Checking the Temp Weighted seed's " << itWeightedseed->first.nHits() << " hits to " << TempRecHits.size() << " Temp recHits" << endl;
            TrajectorySeed::range RecHitsRange = itWeightedseed->first.recHits();
            bool isShare = false;
            for(TrajectorySeed::const_iterator it = RecHitsRange.first; it != RecHitsRange.second; it++)
                if(isShareHit(TempRecHits, *it, rpcGeometry))
                    isShare = true;

            if(isShare == true) {
                if(debug) cout << "Find one Temp seed share some recHits with Best Weighted seed" << endl;
                itWeightedseed = TempWeightedSeeds.erase(itWeightedseed);
            }
            else {
                if(debug) cout << "This seed has no relation with Best Weighted seed" << endl;
                WeightedSeedsRef->push_back(*itWeightedseed);
                itWeightedseed = TempWeightedSeeds.erase(itWeightedseed);
            }
        }
    }
    // At the end exchange SeedsRef with SortSeeds
    WeightedSeedsRef->clear();
    *WeightedSeedsRef = SortWeightedSeeds;
}

bool RPCSeedOverlapper::isShareHit(const edm::OwnVector<TrackingRecHit> &RecHits, const TrackingRecHit& hit, edm::ESHandle<RPCGeometry> rpcGeometry) {

    bool istheSame = false;
    unsigned int n = 1;
    if(debug) cout << "Checking from " << RecHits.size() << " Temp recHits" << endl;

    LocalPoint lpos1 = hit.localPosition();
    DetId RPCId1 = hit.geographicalId();
    const GeomDetUnit *rpcroll1 = rpcGeometry->idToDetUnit(RPCId1);
    GlobalPoint gpos1 = rpcroll1->toGlobal(lpos1);
    if(debug) cout << "The hit's position: " << gpos1.x() << ", " << gpos1.y() << ", " << gpos1.z() << endl;
    for(edm::OwnVector<TrackingRecHit>::const_iterator it = RecHits.begin(); it !=RecHits.end(); it++, n++) {
        if(debug) cout << "Checking the " << n << " th recHit from TempRecHits" << endl;
        LocalPoint lpos2 = it->localPosition();
        DetId RPCId2 = it->geographicalId();
        const GeomDetUnit *rpcroll2 = rpcGeometry->idToDetUnit(RPCId2);
        GlobalPoint gpos2 = rpcroll2->toGlobal(lpos2);
        if(debug) cout << "The Temp hit's position: " << gpos2.x() << ", " << gpos2.y() << ", " << gpos2.z() << endl;

        if((gpos1.x() == gpos2.x()) && (gpos1.y() == gpos2.y()) && (gpos1.z() == gpos2.z())) {
            if(debug) cout << "This hit is found to be the same" << endl;
            istheSame = true;
        }
    }
    return istheSame;
}
