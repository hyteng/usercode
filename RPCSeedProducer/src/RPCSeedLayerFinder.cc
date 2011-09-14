/**
 *  See header file for a description of this class.
 *
 */


#include "MyModule/RPCSeedProducer/src/RPCSeedLayerFinder.h"
#include "MyModule/RPCSeedProducer/src/DebugSignal.h"

using namespace std;
using namespace edm;


RPCSeedLayerFinder::RPCSeedLayerFinder() {

    // Initiate the member
    LayersinRPC.clear();  
    isConfigured = false;
    isInputset = false;
    isOutputset = false;
}

RPCSeedLayerFinder::~RPCSeedLayerFinder() {

}

void RPCSeedLayerFinder::configure(const edm::ParameterSet& iConfig) {

    // Set the configuration
    isCosmic = iConfig.getParameter<bool>("isCosmic");
    isMixBarrelwithEndcap = iConfig.getParameter<bool>("isMixBarrelwithEndcap");
    BarrelLayerRange = iConfig.getParameter< std::vector<unsigned int> >("BarrelLayerRange");
    EndcapLayerRange = iConfig.getParameter< std::vector<unsigned int> >("EndcapLayerRange");
    isSpecialLayers = iConfig.getParameter<bool>("isSpecialLayers");
    LayersinBarrel = iConfig.getParameter< std::vector<unsigned int> >("LayersinBarrel");
    LayersinEndcap = iConfig.getParameter< std::vector<unsigned int> >("LayersinEndcap");
    ConstraintedBarrelLayer = iConfig.getParameter< std::vector<unsigned int> >("ConstraintedBarrelLayer");
    ConstraintedNegativeEndcapLayer = iConfig.getParameter< std::vector<unsigned int> >("ConstraintedNegativeEndcapLayer");
    ConstraintedPositiveEndcapLayer = iConfig.getParameter< std::vector<unsigned int> >("ConstraintedPositiveEndcapLayer");
    
    // Set the signal open
    isConfigured = true;
}

void RPCSeedLayerFinder::setInput(MuonRecHitContainer (&RecHitinRPC)[RPCLayerNumber]) {

    for(unsigned int i = 0; i < RPCLayerNumber; i++)
        RecHitNumberinLayer[i] = RecHitinRPC[i].size();
    // Set the signal open
    isInputset = true;
}

void RPCSeedLayerFinder::unsetInput() {

    isInputset = false;
}

void RPCSeedLayerFinder::setOutput(RPCSeedrecHitFinder* Ref = NULL, RPCCosmicSeedrecHitFinder* CosmicRef = NULL) {

    RPCrecHitFinderRef = Ref;
    RPCCosmicrecHitFinderRef = CosmicRef;
    isOutputset = true;
}

void RPCSeedLayerFinder::fill() {

    // Check if already configured
    if(isConfigured == false || isInputset == false || isOutputset == false) {
        if(debug) cout << "RPCSeedLayerFinder has not set Configure or IO: " << isConfigured << ", " << isInputset << ", " << isOutputset << endl;
        return;
    }

    // Clear the vector LayersinRPC
    LayersinRPC.clear();

    // Now fill the Layers
    if(isCosmic == true) {
        if(RPCCosmicrecHitFinderRef != NULL) {
            if(debug) cout << "filling cosmic layers..." << endl;
            fillCosmicLayers();
        }
        else
            if(debug) cout << "RPCCosmicrecHitFinderRef not set" << endl;
    }
    else {
        if(RPCrecHitFinderRef != NULL) {
            if(debug) cout << "filling collision layers..." << endl;
            fillLayers();
        }
        else
            if(debug) cout << "RPCrecHitFinderRef not set" << endl;
    }
    if(debug) cout << "Finish filling layers." << endl;
}

void RPCSeedLayerFinder::fillLayers() {

    if(isSpecialLayers == false && isMixBarrelwithEndcap == false) {
        for(std::vector<unsigned int>::iterator NumberofLayersinBarrel = BarrelLayerRange.begin(); NumberofLayersinBarrel != BarrelLayerRange.end(); NumberofLayersinBarrel++) {
            // find N layers out of 6 Barrel Layers to fill to SeedinRPC
            unsigned int NumberofLayers = *NumberofLayersinBarrel;
            if(NumberofLayers < 1 || NumberofLayers > BarrelLayerNumber)
                continue;
            int type = 0;  // type=0 for barrel
            LayersinRPC.clear();
            SpecialLayers(-1, NumberofLayers, type);
            LayersinRPC.clear();
        }

        for(std::vector<unsigned int>::iterator NumberofLayersinEndcap = EndcapLayerRange.begin(); NumberofLayersinEndcap != EndcapLayerRange.end(); NumberofLayersinEndcap++) {
            unsigned int NumberofLayers = *NumberofLayersinEndcap;
            if(NumberofLayers < 1 || NumberofLayers > EachEndcapLayerNumber)
                continue;
            int type = 1; // type=1 for endcap
            // for -Z layers
            LayersinRPC.clear();
            SpecialLayers(BarrelLayerNumber-1, NumberofLayers, type);
            LayersinRPC.clear();
            //for +Z layers
            LayersinRPC.clear();
            SpecialLayers(BarrelLayerNumber+EachEndcapLayerNumber-1, NumberofLayers, type);
            LayersinRPC.clear();
        }
    }

    if(isSpecialLayers == true && isMixBarrelwithEndcap == false) {
        // Fill barrel layer for seed
        bool EnoughforBarrel = true;
        unsigned int i = 0;
        LayersinRPC.clear();
        for(std::vector<unsigned int>::iterator it = LayersinBarrel.begin(); it != LayersinBarrel.end(); it++, i++) {   
            if((*it) != 0 && i < BarrelLayerNumber) {
                if(RecHitNumberinLayer[i] != 0)
                    LayersinRPC.push_back(i);
                else {
                    if(debug) cout << "Not recHits in special Barrel layer " << i << endl;
                    EnoughforBarrel = false;
                }
            }
        }
        if(EnoughforBarrel && (LayersinRPC.size() != 0)) {
            // Initiate and call recHit Finder
            RPCrecHitFinderRef->setLayers(LayersinRPC);
            RPCrecHitFinderRef->fillrecHits();
        }
        LayersinRPC.clear();

        // Fill -Z and +Z endcap layer
        bool EnoughforEndcap = true;

        // Fill endcap- layer for seed
        i = BarrelLayerNumber;
        EnoughforEndcap = true;
        LayersinRPC.clear();
        for(std::vector<unsigned int>::iterator it = LayersinEndcap.begin(); it != LayersinEndcap.end(); it++, i++) {
            if((*it) != 0 && i < (BarrelLayerNumber+EachEndcapLayerNumber)) {
                if(RecHitNumberinLayer[i] != 0)
                    LayersinRPC.push_back(i);
                else {
                    if(debug) cout << "Not recHits in special Endcap " << (i - BarrelLayerNumber) << endl;
                    EnoughforEndcap = false;
                }
            }
        }
        if(EnoughforEndcap && (LayersinRPC.size() != 0)) {
            // Initiate and call recHit Finder
            RPCrecHitFinderRef->setLayers(LayersinRPC);
            RPCrecHitFinderRef->fillrecHits();
        }
        LayersinRPC.clear();

        //Fill endcap+ layer for seed
        i = BarrelLayerNumber;
        EnoughforEndcap = true;
        LayersinRPC.clear();
        for(std::vector<unsigned int>::iterator it = LayersinEndcap.begin(); it != LayersinEndcap.end(); it++, i++) {
            if((*it) != 0 && i >= (BarrelLayerNumber+EachEndcapLayerNumber) && i < (BarrelLayerNumber+EachEndcapLayerNumber*2)) {
                if(RecHitNumberinLayer[i] != 0)
                    LayersinRPC.push_back(i);
                else {
                    if(debug) cout << "Not recHits in special Endcap " << i << endl;
                    EnoughforEndcap = false;
                }
            }
        }
        if(EnoughforEndcap && (LayersinRPC.size() != 0)) {
            // Initiate and call recHit Finder
            RPCrecHitFinderRef->setLayers(LayersinRPC);
            RPCrecHitFinderRef->fillrecHits();
        }
        LayersinRPC.clear();
    }

    if(isMixBarrelwithEndcap == true) {
        if(debug) cout <<" Mix is not ready for non-cosmic case" << endl;
        LayersinRPC.clear();
    }
}

void RPCSeedLayerFinder::fillCosmicLayers() {

    // For cosmic only handle the SpecialLayers case
    if(isSpecialLayers == true && isMixBarrelwithEndcap == false) {

        // Fill barrel layer for seed
        unsigned int i = 0;
        LayersinRPC.clear();
        for(std::vector<unsigned int>::iterator it = LayersinBarrel.begin(); it != LayersinBarrel.end(); it++, i++) {   
            if((*it) != 0 && i < BarrelLayerNumber)
                if(RecHitNumberinLayer[i] != 0)
                    LayersinRPC.push_back(i);
        }
        if(LayersinRPC.size() != 0) {
            // Initiate and call recHit Finder
            RPCCosmicrecHitFinderRef->setLayers(LayersinRPC);
            RPCCosmicrecHitFinderRef->fillrecHits();
        }
        LayersinRPC.clear();

        // Fill -Z and +Z endcap layer

        // Fill endcap- layer for seed
        i = BarrelLayerNumber;
        LayersinRPC.clear();
        for(std::vector<unsigned int>::iterator it = LayersinEndcap.begin(); it != LayersinEndcap.end(); it++, i++) {
            if((*it) != 0 && i < (BarrelLayerNumber+EachEndcapLayerNumber))
                if(RecHitNumberinLayer[i] != 0)
                    LayersinRPC.push_back(i);
        }
        if(LayersinRPC.size() != 0) {
            // Initiate and call recHit Finder
            RPCCosmicrecHitFinderRef->setLayers(LayersinRPC);
            RPCCosmicrecHitFinderRef->fillrecHits();
        }
        LayersinRPC.clear();

        //Fill endcap+ layer for seed
        i = BarrelLayerNumber;
        LayersinRPC.clear();
        for(std::vector<unsigned int>::iterator it = LayersinEndcap.begin(); it != LayersinEndcap.end(); it++, i++) {
            if((*it) != 0 && i >= (BarrelLayerNumber+EachEndcapLayerNumber) && i < (BarrelLayerNumber+EachEndcapLayerNumber*2))
                if(RecHitNumberinLayer[i] != 0)
                    LayersinRPC.push_back(i);
        }
        if(LayersinRPC.size() != 0) {
            // Initiate and call recHit Finder
            RPCCosmicrecHitFinderRef->setLayers(LayersinRPC);
            RPCCosmicrecHitFinderRef->fillrecHits();
        }
        LayersinRPC.clear();
    }

    if(isSpecialLayers == true && isMixBarrelwithEndcap == true) {

        // Fill all
        unsigned int i = 0;
        LayersinRPC.clear();
        for(std::vector<unsigned int>::iterator it = LayersinBarrel.begin(); it != LayersinBarrel.end(); it++, i++) {   
            if((*it) != 0 && i < BarrelLayerNumber)
                if(RecHitNumberinLayer[i] != 0)
                    LayersinRPC.push_back(i);
        }
        i = BarrelLayerNumber;
        for(std::vector<unsigned int>::iterator it = LayersinEndcap.begin(); it != LayersinEndcap.end(); it++, i++) {
            if((*it) != 0 && i < (BarrelLayerNumber+EachEndcapLayerNumber*2))
                if(RecHitNumberinLayer[i] != 0)
                    LayersinRPC.push_back(i);
        }

        if(LayersinRPC.size() != 0) {
            // Initiate and call recHit Finder
            RPCCosmicrecHitFinderRef->setLayers(LayersinRPC);
            RPCCosmicrecHitFinderRef->fillrecHits();
        }
        LayersinRPC.clear();
    }

    if(isSpecialLayers == false) {
        if(debug) cout << "Not ready for not SpecialLayers for Cosmic case" << endl;
        LayersinRPC.clear();
    }
}

void RPCSeedLayerFinder::SpecialLayers(int last, unsigned int NumberofLayers, int type) {

    // check type, 0=barrel, 1=endcap, 2=mix

    // barrel has 6 layers
    if(type == 0) {
        if(NumberofLayers > BarrelLayerNumber) {
            if(debug) cout << "NumberofLayers larger than max layers in barrel" << endl;
            return;
        }
        for(unsigned int i = (last+1); i <= (BarrelLayerNumber-NumberofLayers+LayersinRPC.size()); i++) {
            if(RecHitNumberinLayer[i] != 0) {
                LayersinRPC.push_back(i);
                last = i;
                if(LayersinRPC.size() < NumberofLayers)
                    SpecialLayers(last, NumberofLayers, type);
                else {
                    if(checkBarrelConstrain()) {
                        if(debug) cout << "Find special barrel layers: ";
                        for(unsigned int k = 0; k < NumberofLayers; k++)
                            if(debug) cout << LayersinRPC[k] <<" ";
                        if(debug) cout << endl;
                        // Initiate and call recHit Finder
                        RPCrecHitFinderRef->setLayers(LayersinRPC);
                        RPCrecHitFinderRef->fillrecHits();
                    }
                    else
                        if(debug) cout << "The layers don't contain all layers in constrain" << endl;
                }
                LayersinRPC.pop_back();
            }
        }
    }

    // endcap has 3 layers for each -Z and +Z
    if(type == 1) {
        if(NumberofLayers > EachEndcapLayerNumber) {
            if(debug) cout << "NumberofLayers larger than max layers in endcap" << endl;
            return;
        }
        if(last < (BarrelLayerNumber+EachEndcapLayerNumber-1) || (last == (BarrelLayerNumber+EachEndcapLayerNumber-1) && LayersinRPC.size() != 0)) {
            // For -Z case
            for(unsigned int i =  (last+1); i <= (BarrelLayerNumber+EachEndcapLayerNumber-NumberofLayers+LayersinRPC.size()); i++) {
                if(RecHitNumberinLayer[i] != 0) {
                    LayersinRPC.push_back(i);
                    last = i;
                    if(LayersinRPC.size() < NumberofLayers)
                        SpecialLayers(last, NumberofLayers, type);
                    else {
                        if(checkNegativeEndcapConstrain()) {
                            if(debug) cout << "Find special -Z endcap layers: ";
                            for(unsigned int k = 0; k < NumberofLayers; k++)
                                if(debug) cout << LayersinRPC[k] <<" ";
                            if(debug) cout << endl;
                            // Initiate and call recHit Finder
                            RPCrecHitFinderRef->setLayers(LayersinRPC);
                            RPCrecHitFinderRef->fillrecHits();
                        }
                        else
                            if(debug) cout << "The layers don't contain all layers in constrain" << endl;
                    }
                    LayersinRPC.pop_back();
                }
            }
        }
        else {
            // For +Z case
            for(unsigned int i = (last+1); i <= (BarrelLayerNumber+EachEndcapLayerNumber*2-NumberofLayers+LayersinRPC.size()); i++) {
                if(RecHitNumberinLayer[i] != 0) {
                    LayersinRPC.push_back(i);
                    last = i;
                    if(LayersinRPC.size() < NumberofLayers)
                        SpecialLayers(last, NumberofLayers, type);
                    else {
                        if(checkPositiveEndcapConstrain()) {
                            if(debug) cout << "Find special +Z endcap layers: ";
                            for(unsigned int k = 0; k < NumberofLayers; k++)
                                if(debug) cout << LayersinRPC[k] <<" ";
                            if(debug) cout << endl;
                            // Initiate and call recHit Finder
                            RPCrecHitFinderRef->setLayers(LayersinRPC);
                            RPCrecHitFinderRef->fillrecHits();
                        }
                        else
                            if(debug) cout << "The layers don't contain all layers in constrain" << endl;
                    }
                    LayersinRPC.pop_back();
                }
            }
        }
    }
}

bool RPCSeedLayerFinder::checkBarrelConstrain() {

    bool pass = true;
    std::vector<unsigned int> fitConstrain = ConstraintedBarrelLayer;
    for(unsigned int i = 0; i < LayersinRPC.size(); i++)
        fitConstrain[LayersinRPC[i]] = 0;
    for(unsigned int i = 0; i < BarrelLayerNumber; i++)
        if(fitConstrain[i] != 0)
            pass = false;
    return pass;
}

bool RPCSeedLayerFinder::checkPositiveEndcapConstrain() {

    bool pass = true;
    std::vector<unsigned int> fitConstrain = ConstraintedPositiveEndcapLayer;
    for(unsigned int i = 0; i < LayersinRPC.size(); i++)
        fitConstrain[LayersinRPC[i]-BarrelLayerNumber-EachEndcapLayerNumber] = 0;
    for(unsigned int i = 0; i < EachEndcapLayerNumber; i++)
        if(fitConstrain[i] != 0)
            pass = false;
    return pass;
}

bool RPCSeedLayerFinder::checkNegativeEndcapConstrain() {

    bool pass = true;
    std::vector<unsigned int> fitConstrain = ConstraintedNegativeEndcapLayer;
    for(unsigned int i = 0; i < LayersinRPC.size(); i++)
        fitConstrain[LayersinRPC[i]-BarrelLayerNumber] = 0;
    for(unsigned int i = 0; i < EachEndcapLayerNumber; i++)
        if(fitConstrain[i] != 0)
            pass = false;
    return pass;
}
