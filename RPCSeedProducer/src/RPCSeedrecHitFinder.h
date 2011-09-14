#ifndef MyModule_RPCSeedProducer_RPCSeedrecHitFinder_H
#define MyModule_RPCSeedProducer_RPCSeedrecHitFinder_H

/** \class RPCSeedLayerFinder
 *  
 *   \author Haiyun.Teng - Peking University
 *
 *  
 */


#include "MyModule/RPCSeedProducer/src/RPCSeedFinder.h"
#include "MyModule/RPCSeedProducer/src/RPCSeedData.h"
#include <RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>


class RPCSeedrecHitFinder {

    typedef MuonTransientTrackingRecHit::MuonRecHitPointer MuonRecHitPointer;
    typedef MuonTransientTrackingRecHit::ConstMuonRecHitPointer ConstMuonRecHitPointer;
    typedef MuonTransientTrackingRecHit::MuonRecHitContainer MuonRecHitContainer;
    typedef MuonTransientTrackingRecHit::ConstMuonRecHitContainer ConstMuonRecHitContainer;

    public:
        RPCSeedrecHitFinder();
        ~RPCSeedrecHitFinder();
        void configure(const edm::ParameterSet& iConfig);
        void setInput(MuonRecHitContainer (&recHits)[RPCLayerNumber]);
        void unsetInput();
        void setOutput(RPCSeedFinder *Seed);
        void setLayers(const std::vector<unsigned int>& Layers);
        void fillrecHits();
    private:
        void complete(unsigned int LayerIndex);
        double getdeltaPhifromrecHits();
        void checkandfill();

        // ----------member data ---------------------------

        // parameters for configuration
        unsigned int BxRange;
        double MaxDeltaPhi;
        std::vector<int> ClusterSet;
        // Signal for call fillrecHits()
        bool isLayerset;
        bool isConfigured;
        bool isInputset;
        bool isOutputset;
        // Enable layers in Barrel and Endcap
        std::vector<unsigned int> RPCLayers;
        MuonRecHitContainer* recHitsRPC[RPCLayerNumber];
        GlobalVector theInitialVector;
        ConstMuonRecHitContainer theRecHits;
        RPCSeedFinder *theSeed;
};

#endif
