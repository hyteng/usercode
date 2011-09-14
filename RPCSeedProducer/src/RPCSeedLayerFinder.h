#ifndef MyModule_RPCSeedProducer_RPCSeedLayerFinder_H
#define MyModule_RPCSeedProducer_RPCSeedLayerFinder_H

/** \class RPCSeedLayerFinder
 *  
 *   \author Haiyun.Teng - Peking University
 *
 *  
 */


#include <RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h>
#include "MyModule/RPCSeedProducer/src/RPCSeedData.h"
#include "MyModule/RPCSeedProducer/src/RPCSeedrecHitFinder.h"
#include "MyModule/RPCSeedProducer/src/RPCCosmicSeedrecHitFinder.h"
#include <FWCore/ParameterSet/interface/ParameterSet.h>


class RPCSeedLayerFinder {

    typedef MuonTransientTrackingRecHit::MuonRecHitPointer MuonRecHitPointer;
    typedef MuonTransientTrackingRecHit::ConstMuonRecHitPointer ConstMuonRecHitPointer;
    typedef MuonTransientTrackingRecHit::MuonRecHitContainer MuonRecHitContainer;
    typedef MuonTransientTrackingRecHit::ConstMuonRecHitContainer ConstMuonRecHitContainer;

    public:
        RPCSeedLayerFinder();
        ~RPCSeedLayerFinder();
        void configure(const edm::ParameterSet& iConfig);
        void setInput(MuonRecHitContainer (&RecHitinRPC)[RPCLayerNumber]);
        void unsetInput();
        void setOutput(RPCSeedrecHitFinder* Ref, RPCCosmicSeedrecHitFinder* CosmicRef);
        void fill();

    private:
        void fillLayers();
        void fillCosmicLayers();
        // create special N layers to fill to seeds
        void SpecialLayers(int last, unsigned int NumberofLayers, int type);
        bool checkBarrelConstrain();
        bool checkNegativeEndcapConstrain();
        bool checkPositiveEndcapConstrain();

        // ----------member data ---------------------------

        // The ref of RPCSeedrecHitFinder which will be call after gathering a set of layers 
        RPCSeedrecHitFinder* RPCrecHitFinderRef;
        RPCCosmicSeedrecHitFinder* RPCCosmicrecHitFinderRef;
        // The parameters for configuration
        bool isCosmic;
        bool isMixBarrelwithEndcap;
        std::vector<unsigned int> BarrelLayerRange;
        std::vector<unsigned int> EndcapLayerRange;
        bool isSpecialLayers;
        std::vector<unsigned int> LayersinEndcap;
        std::vector<unsigned int> LayersinBarrel;
        std::vector<unsigned int> ConstraintedBarrelLayer;
        std::vector<unsigned int> ConstraintedNegativeEndcapLayer;
        std::vector<unsigned int> ConstraintedPositiveEndcapLayer;
        // Signal for call fillLayers()
        bool isConfigured;
        bool isInputset;
        bool isOutputset;
        // Enable layers in Barrel and Endcap
        std::vector<unsigned int> RPCLayers;
        // Information of recHits in each layer
        unsigned int RecHitNumberinLayer[RPCLayerNumber];
};

#endif
