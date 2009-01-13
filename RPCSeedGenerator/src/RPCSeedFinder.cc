/**
 *  See header file for a description of this class.
 *
 */


#include "MyMTCCAnalyzer/RPCSeedGenerator/src/RPCSeedFinder.h"
#include "MyMTCCAnalyzer/RPCSeedGenerator/src/RPCSeedHits.h"
#include "RecoMuon/TrackingTools/interface/MuonPatternRecoDumper.h"

#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "DataFormats/GeometryVector/interface/Pi.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "DataFormats/GeometryVector/interface/CoordinateSets.h"
#include "DataFormats/GeometrySurface/interface/BoundPlane.h"
#include "DataFormats/GeometrySurface/interface/RectangularPlaneBounds.h"

#include "DataFormats/Common/interface/OwnVector.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrajectoryState/interface/PTrajectoryStateOnDet.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include <iomanip>

using namespace std;

typedef MuonTransientTrackingRecHit::MuonRecHitPointer MuonRecHitPointer;
typedef MuonTransientTrackingRecHit::ConstMuonRecHitPointer ConstMuonRecHitPointer;
typedef MuonTransientTrackingRecHit::MuonRecHitContainer MuonRecHitContainer;

// Comparator function must be global function?
// Could not in include in .h file, or there will be 2 lessPhi() functions in both RPCSeedGenerator and RPCSeedFinder module
bool lessPhi(const MuonTransientTrackingRecHit::MuonRecHitPointer& it1, const MuonTransientTrackingRecHit::MuonRecHitPointer& it2) {
    
    return (it1->globalPosition().phi() < it2->globalPosition().phi());
}


RPCSeedFinder::RPCSeedFinder(){}

//vector<TrajectorySeed> RPCSeedFinder::seeds(const edm::EventSetup& eSetup, double MaxRSD) const {
void RPCSeedFinder::seeds(const edm::EventSetup& eSetup, double MaxRSD, std::vector<TrajectorySeed>& goodSeeds, std::vector<TrajectorySeed>& candidateSeeds) const {
    
    cout << "[RPCSeedFinder] --> seeds class called" << endl;
    
    //vector<TrajectorySeed> goodSeeds(0);
    //vector<TrajectorySeed> candidateSeeds(0);
    TrajectorySeed theSeed;
    RPCSeedHits oneSeed;
    
    int NumberofHitsinSeed = 0;
    for (MuonRecHitContainer::const_iterator iter = theRecHits.begin(); iter!= theRecHits.end(); iter++)
    {
        oneSeed.add(*iter);
        NumberofHitsinSeed++;
    }
    int isGoodSeed = 0;
    if(NumberofHitsinSeed >= 3)
    {
        theSeed = oneSeed.seed(eSetup, MaxRSD, isGoodSeed);
        // Push back the good seed
        if(isGoodSeed == 1)
        {
            cout << "[RPCSeedFinder] --> Seeds from " << NumberofHitsinSeed << " recHits." << endl;
            goodSeeds.push_back(theSeed);
        }
        // Push back cadidate seed but not the fake seed
        if(isGoodSeed >= 0)
        {
            candidateSeeds.push_back(theSeed);
        }
    }
    //return theSeeds;
}

double RPCSeedFinder::GetdeltaPhifromHits() {
    
    MuonTransientTrackingRecHit::MuonRecHitContainer sortRecHits = theRecHits;
    sort(sortRecHits.begin(), sortRecHits.end(), lessPhi);
    cout << "Sorted recHit's Phi: ";
    for(MuonRecHitContainer::const_iterator iter = sortRecHits.begin(); iter != sortRecHits.end(); iter++)
        cout << (*iter)->globalPosition().phi() << ", ";
    cout << endl;
    // Calculate the deltaPhi, take care Geom::Phi always in range [-pi,pi)
    // In case of some deltaPhi larger then pi, use value() in Geom::Phi to get the true value in radians of Phi, then do the calculation
    double deltaPhi = 0;
    if(sortRecHits.size() <= 1)
        return deltaPhi;
    if(sortRecHits.size() == 2)
    {
        MuonRecHitContainer::const_iterator iter1 = sortRecHits.begin();
        MuonRecHitContainer::const_iterator iter2 = sortRecHits.begin();
        iter2++;
        deltaPhi = (((*iter2)->globalPosition().phi().value() - (*iter1)->globalPosition().phi().value()) > M_PI) ? (2 * M_PI - ((*iter2)->globalPosition().phi().value() - (*iter1)->globalPosition().phi().value())) : ((*iter2)->globalPosition().phi().value() - (*iter1)->globalPosition().phi().value());
        return deltaPhi;
    }
    else
    {
        deltaPhi = 2 * M_PI;
        int n = 0;
        for(MuonRecHitContainer::const_iterator iter = sortRecHits.begin(); iter != sortRecHits.end(); iter++)
        {   
            cout << "Before this loop deltaPhi is " << deltaPhi << endl;
            n++;
            double deltaPhi_more = 0;
            double deltaPhi_less = 0;
            if(iter == sortRecHits.begin())
            {
                cout << "Calculateing frist loop..." << endl;
                MuonRecHitContainer::const_iterator iter_more = ++iter;
                --iter;
                MuonRecHitContainer::const_iterator iter_less = sortRecHits.end();
                --iter_less;
                cout << "more_Phi: " << (*iter_more)->globalPosition().phi() << ", less_Phi: " << (*iter_less)->globalPosition().phi() << ", iter_Phi: " << (*iter)->globalPosition().phi() << endl;
                deltaPhi_more = (2 * M_PI) - ((*iter_more)->globalPosition().phi().value() - (*iter)->globalPosition().phi().value());
                deltaPhi_less = (*iter_less)->globalPosition().phi().value() - (*iter)->globalPosition().phi().value();
            }
            else if(iter == (--sortRecHits.end()))
            {
                cout << "Calculateing last loop..." << endl;
                MuonRecHitContainer::const_iterator iter_less = --iter;
                ++iter;
                MuonRecHitContainer::const_iterator iter_more = sortRecHits.begin();
                cout << "more_Phi: " << (*iter_more)->globalPosition().phi() << ", less_Phi: " << (*iter_less)->globalPosition().phi() << ", iter_Phi: " << (*iter)->globalPosition().phi() << endl;
                deltaPhi_less = (2 * M_PI) - ((*iter)->globalPosition().phi().value() - (*iter_less)->globalPosition().phi().value());
                deltaPhi_more = (*iter)->globalPosition().phi().value() - (*iter_more)->globalPosition().phi().value();
            }
            else
            {
                cout << "Calculateing " << n <<"st loop..." << endl;
                MuonRecHitContainer::const_iterator iter_less = --iter;
                ++iter;
                MuonRecHitContainer::const_iterator iter_more = ++iter;
                --iter;
                cout << "more_Phi: " << (*iter_more)->globalPosition().phi() << ", less_Phi: " << (*iter_less)->globalPosition().phi() << ", iter_Phi: " << (*iter)->globalPosition().phi() << endl;
                deltaPhi_less = (2 * M_PI) - ((*iter)->globalPosition().phi().value() - (*iter_less)->globalPosition().phi().value());
                deltaPhi_more = (2 * M_PI) - ((*iter_more)->globalPosition().phi().value() - (*iter)->globalPosition().phi().value());
            }
            if(deltaPhi > deltaPhi_more)
                deltaPhi = deltaPhi_more;
            if(deltaPhi > deltaPhi_less)
                deltaPhi = deltaPhi_less;
            
            cout << "For this loop deltaPhi_more is " << deltaPhi_more << endl;
            cout << "For this loop deltaPhi_less is " << deltaPhi_less << endl;
            cout << "For this loop deltaPhi is " << deltaPhi << endl;
        }
        return deltaPhi;
    }
}

