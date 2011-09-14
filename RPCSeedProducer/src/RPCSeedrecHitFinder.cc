/**
 *  See header file for a description of this class.
 *
 */


#include "MyModule/RPCSeedProducer/src/RPCSeedrecHitFinder.h"
#include "MyModule/RPCSeedProducer/src/DebugSignal.h"
#include <DataFormats/TrackingRecHit/interface/TrackingRecHit.h>
#include <DataFormats/RPCRecHit/interface/RPCRecHit.h>

using namespace std;
using namespace edm;


// Comparator function must be global function?
// Could not be included in .h file, or there will be 2 lessPhi() functions in both RPCSeedGenerator and RPCSeedFinder module
bool LessR(const MuonTransientTrackingRecHit::ConstMuonRecHitPointer& it1, const MuonTransientTrackingRecHit::ConstMuonRecHitPointer& it2) {
    // Don't need to use value() in Geom::Phi to short
    return (it1->globalPosition().perp() < it2->globalPosition().perp());
}

RPCSeedrecHitFinder::RPCSeedrecHitFinder() {

    // Initiate the member
    isLayerset = false;
    isConfigured = false;
    isInputset = false;
    isOutputset = false;
    BxRange = 0;
    MaxDeltaPhi = 0;
    ClusterSet.clear();
    LayersinRPC.clear();
    theRecHits.clear();
}

RPCSeedrecHitFinder::~RPCSeedrecHitFinder() {

}

void RPCSeedrecHitFinder::configure(const edm::ParameterSet& iConfig) {

    // Set the configuration
    BxRange = iConfig.getParameter<unsigned int>("BxRange");
    MaxDeltaPhi = iConfig.getParameter<double>("MaxDeltaPhi");
    ClusterSet = iConfig.getParameter< std::vector<int> >("ClusterSet");

    // Set the signal open
    isConfigured = true;
}

void RPCSeedrecHitFinder::setInput(MuonRecHitContainer (&recHits)[RPCLayerNumber]) {

    for(unsigned int i = 0; i < RPCLayerNumber; i++)
        recHitsRPC[i] = &recHits[i];
    isInputset = true;
}

void RPCSeedrecHitFinder::unsetInput() {

    isInputset = false;
}

void RPCSeedrecHitFinder::setOutput(RPCSeedFinder *Seed) {

    theSeed = Seed;
    isOutputset = true;
}

void RPCSeedrecHitFinder::setLayers(const std::vector<unsigned int>& Layers) {

    LayersinRPC = Layers;
    isLayerset = true;
}

void RPCSeedrecHitFinder::fillrecHits() {

    if(isLayerset == false || isConfigured == false || isOutputset == false || isInputset == false) {
        if(debug) cout << "SeedrecHitFinder has not set the IO or not configured yet: " << isLayerset << ", " << isConfigured << ", " << isOutputset << ", " << isInputset << endl;
        return;
    }
   if(debug) cout << "Now fill recHits from Layers: ";
    
    for(unsigned int k = 0; k < LayersinRPC.size(); k++)
        if(debug) cout << LayersinRPC[k] <<" ";
    if(debug) cout << endl;
    unsigned int LayerIndex = 0;
    theRecHits.clear();
    complete(LayerIndex);

    // Unset the signal
    if(debug) cout << "Finish filling recHits. " << endl;
    LayersinRPC.clear();
    isLayerset = false;
    theRecHits.clear();
}

void RPCSeedrecHitFinder::complete(unsigned int LayerIndex) {

    for(MuonRecHitContainer::const_iterator it = recHitsRPC[LayersinRPC[LayerIndex]]->begin(); it != recHitsRPC[LayersinRPC[LayerIndex]]->end(); it++) {
        if(debug) cout << "Completing layer[" << LayersinRPC[LayerIndex] << "]." << endl;

        // Check validation
        if(!(*it)->isValid())
            continue;

        // Check BX range, be sure there is only RPCRecHit in the MuonRecHitContainer when use the dynamic_cast
        TrackingRecHit* thisTrackingRecHit = (*it)->hit()->clone();
        // Should also delete the RPCRecHit object cast by dynamic_cast<> ?
        RPCRecHit* thisRPCRecHit = dynamic_cast<RPCRecHit*>(thisTrackingRecHit);
        int BX = thisRPCRecHit->BunchX();
        int ClusterSize = thisRPCRecHit->clusterSize();
        delete thisTrackingRecHit;
        // Check BX
        if((unsigned int)abs(BX) > BxRange)
            continue;
        // Check cluster size
        bool Clustercheck = false;
        if(ClusterSet.size() == 0)
            Clustercheck = true;
        for(std::vector<int>::const_iterator CluIter = ClusterSet.begin(); CluIter != ClusterSet.end(); CluIter++)
            if(ClusterSize == (*CluIter))
                Clustercheck = true;
        if(Clustercheck != true)
            continue;
        // Check the recHits Phi range
        GlobalPoint pos = (*it)->globalPosition();
        double Phi = pos.phi();
        if(debug) cout << "Phi: " << Phi << endl;
        // The recHits should locate in some phi range
        theRecHits.push_back(*it);
        double deltaPhi = getdeltaPhifromrecHits();
        if(debug) cout << "Delta phi: "<< deltaPhi << endl;
        theRecHits.pop_back();
        
        if(fabs(deltaPhi) > MaxDeltaPhi)
            continue;

        // If pass all, add to the seed
        theRecHits.push_back(*it);
        if(debug) cout << "RecHit's global position: " << pos.x() << ", " << pos.y() << ", " << pos.z() << endl;

        // Check if this recHit is the last one in the seed
        // If it is the last one, calculate the seed
        // If it is not the last one, continue to fill the seed from other layers
        if(LayerIndex == (LayersinRPC.size()-1)) {
            if(debug) cout << "Check and fill one seed." << endl;
            checkandfill();
        }
        else
            complete(LayerIndex+1);

        // Remember to pop the recHit before add another one from the same layer!
        theRecHits.pop_back();
    }
}

double RPCSeedrecHitFinder::getdeltaPhifromrecHits() {

    ConstMuonRecHitContainer SortRecHits = theRecHits;
    sort(SortRecHits.begin(), SortRecHits.end(), LessR);
    if(debug) cout << "Sorted recHit's R: ";
    for(ConstMuonRecHitContainer::const_iterator Iter = SortRecHits.begin(); Iter != SortRecHits.end(); Iter++)
        if(debug) cout << (*Iter)->globalPosition().perp() << ", ";
    if(debug) cout << endl;
    // Calculate the deltaPhi, take care Geom::Phi always in range [-pi,pi)
    // In case of some deltaPhi larger then Pi, use value() in Geom::Phi to get the true value in radians of Phi, then do the calculation
    double DeltaPhi = 0;
    int n = SortRecHits.size();
    if(n == 2) {
        ConstMuonRecHitContainer::const_iterator StartIter = SortRecHits.begin();
        ConstMuonRecHitContainer::const_iterator EndIter = SortRecHits.end();
        EndIter--;
        theInitialVector = (*EndIter)->globalPosition() - (*StartIter)->globalPosition();
        DeltaPhi = ((*EndIter)->globalPosition().phi() - (*StartIter)->globalPosition().phi()).value();
    }
    if(n >= 3) {
        ConstMuonRecHitContainer::const_iterator StartIter = SortRecHits.begin();
        ConstMuonRecHitContainer::const_iterator SampleIter = SortRecHits.end();
        SampleIter--;
        GlobalVector SampleVector = (*SampleIter)->globalPosition() - (*StartIter)->globalPosition();
        DeltaPhi = (SampleVector.phi() - theInitialVector.phi()).value();
        if(debug) cout << "sample DeltaPhi: " << DeltaPhi << endl;
    }
    return DeltaPhi;
    /*
    if(SortRecHits.size() <= 1)
        return deltaPhi;
    if(SortRecHits.size() == 2) {
        ConstMuonRecHitContainer::const_iterator iter1 = SortRecHits.begin();
        ConstMuonRecHitContainer::const_iterator iter2 = SortRecHits.begin();
        iter2++;
        theInitialVector = iter2->globalPosition() - iter1->globalPosition();
        deltaPhi = 0.;
        //deltaPhi = (((*iter2)->globalPosition().phi().value() - (*iter1)->globalPosition().phi().value()) > M_PI) ? (2 * M_PI - ((*iter2)->globalPosition().phi().value() - (*iter1)->globalPosition().phi().value())) : ((*iter2)->globalPosition().phi().value() - (*iter1)->globalPosition().phi().value());
        return deltaPhi;
    }
    else {
        deltaPhi = 2 * M_PI;
        int n = 0;
        for(ConstMuonRecHitContainer::const_iterator iter = SortRecHits.begin(); iter != SortRecHits.end(); iter++) {   
            if(debug) cout << "Before this loop deltaPhi is " << deltaPhi << endl;
            n++;
            double deltaPhi_more = 0;
            double deltaPhi_less = 0;
            if(iter == SortRecHits.begin()) {
                if(debug) cout << "Calculateing frist loop..." << endl;
                ConstMuonRecHitContainer::const_iterator iter_more = ++iter;
                --iter;
                ConstMuonRecHitContainer::const_iterator iter_less = SortRecHits.end();
                --iter_less;
                if(debug) cout << "more_Phi: " << (*iter_more)->globalPosition().phi() << ", less_Phi: " << (*iter_less)->globalPosition().phi() << ", iter_Phi: " << (*iter)->globalPosition().phi() << endl;
                deltaPhi_more = (2 * M_PI) - ((*iter_more)->globalPosition().phi().value() - (*iter)->globalPosition().phi().value());
                deltaPhi_less = (*iter_less)->globalPosition().phi().value() - (*iter)->globalPosition().phi().value();
            }
            else if(iter == (--SortRecHits.end())) {
                if(debug) cout << "Calculateing last loop..." << endl;
                ConstMuonRecHitContainer::const_iterator iter_less = --iter;
                ++iter;
                ConstMuonRecHitContainer::const_iterator iter_more = SortRecHits.begin();
                if(debug) cout << "more_Phi: " << (*iter_more)->globalPosition().phi() << ", less_Phi: " << (*iter_less)->globalPosition().phi() << ", iter_Phi: " << (*iter)->globalPosition().phi() << endl;
                deltaPhi_less = (2 * M_PI) - ((*iter)->globalPosition().phi().value() - (*iter_less)->globalPosition().phi().value());
                deltaPhi_more = (*iter)->globalPosition().phi().value() - (*iter_more)->globalPosition().phi().value();
            }
            else {
                if(debug) cout << "Calculateing " << n << "st loop..." << endl;
                ConstMuonRecHitContainer::const_iterator iter_less = --iter;
                ++iter;
                ConstMuonRecHitContainer::const_iterator iter_more = ++iter;
                --iter;
                if(debug) cout << "more_Phi: " << (*iter_more)->globalPosition().phi() << ", less_Phi: " << (*iter_less)->globalPosition().phi() << ", iter_Phi: " << (*iter)->globalPosition().phi() << endl;
                deltaPhi_less = (2 * M_PI) - ((*iter)->globalPosition().phi().value() - (*iter_less)->globalPosition().phi().value());
                deltaPhi_more = (2 * M_PI) - ((*iter_more)->globalPosition().phi().value() - (*iter)->globalPosition().phi().value());
            }
            if(deltaPhi > deltaPhi_more)
                deltaPhi = deltaPhi_more;
            if(deltaPhi > deltaPhi_less)
                deltaPhi = deltaPhi_less;

            if(debug) cout << "For this loop deltaPhi_more is " << deltaPhi_more << endl;
            if(debug) cout << "For this loop deltaPhi_less is " << deltaPhi_less << endl;
            if(debug) cout << "For this loop deltaPhi is " << deltaPhi << endl;
        }
        return deltaPhi;
    }
    */
}

void RPCSeedrecHitFinder::checkandfill() {

    if(theRecHits.size() >= 3) {
        theSeed->clear();
        theSeed->setRecHits(theRecHits);
        //theSeed->setLayers(LayersinRPC);
        theSeed->seed();
    }
    else
        if(debug) cout << "Layer less than 3, could not fill a RPCSeedFinder" << endl;
}
