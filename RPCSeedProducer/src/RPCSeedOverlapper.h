#ifndef MyModule_RPCSeedProducer_RPCSeedOverlapper_H
#define MyModule_RPCSeedProducer_RPCSeedOverlapper_H

/**  \class RPCSeedPattern
 *
 *  \author Haiyun.Teng - Peking University
 *
 *
 */


#include <DataFormats/TrajectorySeed/interface/TrajectorySeed.h>
#include <FWCore/Framework/interface/EventSetup.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <Geometry/RPCGeometry/interface/RPCGeometry.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <DataFormats/Common/interface/OwnVector.h>
#include <DataFormats/TrackingRecHit/interface/TrackingRecHit.h>
#include "MyModule/RPCSeedProducer/src/RPCSeedPattern.h"

class RPCSeedOverlapper {

    typedef RPCSeedPattern::WeightedTrajectorySeed WeightedTrajectorySeed;

    public:
        RPCSeedOverlapper();
        ~RPCSeedOverlapper();
        void setIO(std::vector<WeightedTrajectorySeed> *GoodWeightedRef, std::vector<WeightedTrajectorySeed> *CandidateWeightedRef);
        void unsetIO();
        void run();    
        void configure(const edm::ParameterSet& iConfig);
        void setEventSetup(const edm::EventSetup& iSetup);
    private:
        void CheckOverlap(const edm::EventSetup& iSetup, std::vector<WeightedTrajectorySeed> *SeedsRef);
        bool isShareHit(const edm::OwnVector<TrackingRecHit> &RecHits, const TrackingRecHit& hit, edm::ESHandle<RPCGeometry> rpcGeometry);
        // Signal for call run()
        bool isConfigured;
        bool isIOset;
        bool isEventSetupset;
        // Parameters for configuration
        bool isCheckGoodOverlap;
        bool isCheckCandidateOverlap;
        unsigned int ShareRecHitsNumberThreshold;
        // IO ref
        std::vector<WeightedTrajectorySeed> *GoodWeightedSeedsRef;
        std::vector<WeightedTrajectorySeed> *CandidateWeightedSeedsRef;
        const edm::EventSetup *eSetup;
};

#endif
