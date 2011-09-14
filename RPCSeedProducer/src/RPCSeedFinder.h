#ifndef MyModule_RPCSeedProducer_RPCSeedFinder_H
#define MyModule_RPCSeedProducer_RPCSeedFinder_H


/** \class RPCSeedFinder
 *  
 *   \author Haiyun.Teng - Peking University
 *
 *  
 */


#include "MyModule/RPCSeedProducer/src/RPCSeedPattern.h"
#include "MyModule/RPCSeedProducer/src/SimRPCSeedPattern.h"
#include <DataFormats/TrajectorySeed/interface/TrajectorySeed.h>
#include <RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/EventSetup.h>
#include <vector>
#include <algorithm>

namespace edm {class EventSetup;}

class RPCSeedFinder {
    typedef MuonTransientTrackingRecHit::MuonRecHitPointer MuonRecHitPointer;
    typedef MuonTransientTrackingRecHit::ConstMuonRecHitPointer ConstMuonRecHitPointer;
    typedef MuonTransientTrackingRecHit::MuonRecHitContainer MuonRecHitContainer;
    typedef MuonTransientTrackingRecHit::ConstMuonRecHitContainer ConstMuonRecHitContainer;
    typedef RPCSeedPattern::WeightedTrajectorySeed WeightedTrajectorySeed;

    public:
        RPCSeedFinder();
        ~RPCSeedFinder();
        void configure(const edm::ParameterSet& iConfig);
        void setOutput(std::vector<WeightedTrajectorySeed> *GoodWeightedRef, std::vector<WeightedTrajectorySeed> *CandidateWeightedRef);
        void setRecHits(ConstMuonRecHitContainer &RecHits);
        void setSimData(const edm::Event& iEvent, const edm::EventSetup& iSetup);
        void unsetSimData();
        //void setLayers(const std::vector<unsigned int>& Layers);
        void clear();
        void setEventSetup(const edm::EventSetup& iSetup);
        void seed();

    private:
        bool useSimData;
        // Signal for call fillLayers()
        bool isRecHitset;
        //bool isLayerset;
        bool isConfigured;
        bool isOutputset;
        bool isEventSetupset;
        const edm::EventSetup *eSetup;
        RPCSeedPattern OneSeed;
        SimRPCSeedPattern SimOneSeed;
        //ConstMuonRecHitContainer theRecHits;
        std::vector<WeightedTrajectorySeed> *GoodWeightedSeedsRef;
        std::vector<WeightedTrajectorySeed> *CandidateWeightedSeedsRef;
};
#endif
