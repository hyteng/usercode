/*
 *  See header file for a description of this class.
 *
 *
 *  $Date: 2009/10/31 02:00:48 $
 *  $Revision: 1.2 $
 *  \author Haiyun.Teng - Peking University
 *
 */

#include "MyModule/RPCSeedProducer/src/RPCSeedPattern.h"
#include "MyModule/RPCSeedProducer/src/RPCSeedData.h"
#include "MyModule/RPCSeedProducer/src/DebugSignal.h"
#include <RecoMuon/TrackingTools/interface/MuonPatternRecoDumper.h>
#include <MagneticField/Engine/interface/MagneticField.h>
#include <MagneticField/Records/interface/IdealMagneticFieldRecord.h>
#include <TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h>
#include <TrackingTools/DetLayers/interface/DetLayer.h>
#include <DataFormats/TrajectoryState/interface/PTrajectoryStateOnDet.h>
#include <DataFormats/Common/interface/OwnVector.h>
#include <DataFormats/MuonDetId/interface/RPCDetId.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h>
#include <Geometry/CommonDetUnit/interface/GeomDet.h>
#include <Geometry/RPCGeometry/interface/RPCChamber.h>
#include <Geometry/RPCGeometry/interface/RPCGeometry.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <Geometry/CommonDetUnit/interface/GeomDetUnit.h>

#include "gsl/gsl_statistics.h"
#include "TH1F.h"
#include "math.h"

using namespace std;
using namespace edm;


RPCSeedPattern::RPCSeedPattern() {
    isPatternChecked = false;
    isConfigured = false;
    MagnecticFieldThreshold = 0.5;
}

RPCSeedPattern::~RPCSeedPattern() {
}


void RPCSeedPattern::configure(const edm::ParameterSet& iConfig) {

    ZError = iConfig.getParameter<double>("ZError");
    MagnecticFieldThreshold = iConfig.getParameter<double>("MagnecticFieldThreshold");
    sampleCount = iConfig.getParameter<unsigned int>("sampleCount");
    AlgorithmType = iConfig.getParameter<unsigned int>("AlgorithmType");
    BendingPhiLowerTH = iConfig.getParameter< vector<double> >("BendingPhiLowerTH");
    BendingPhiUpperTH = iConfig.getParameter< vector<double> >("BendingPhiUpperTH");
    MeanPt_Parameter0 = iConfig.getParameter< vector<double> >("MeanPt_Parameter0");
    MeanPt_Parameter1 = iConfig.getParameter< vector<double> >("MeanPt_Parameter1");
    MeanPt_Parameter2 = iConfig.getParameter< vector<double> >("MeanPt_Parameter2");
    SigmaPt_Parameter0 = iConfig.getParameter< vector<double> >("SigmaPt_Parameter0");
    SigmaPt_Parameter1 = iConfig.getParameter< vector<double> >("SigmaPt_Parameter1");
    SigmaPt_Parameter2 = iConfig.getParameter< vector<double> >("SigmaPt_Parameter2");
    isConfigured = true;

    if(BendingPhiLowerTH.size() == 0 || BendingPhiUpperTH.size() == 0 || MeanPt_Parameter0.size() == 0 || MeanPt_Parameter1.size() == 0 || MeanPt_Parameter2.size() == 0 || SigmaPt_Parameter0.size() == 0 || SigmaPt_Parameter1.size() == 0 || SigmaPt_Parameter2.size() == 0 )
        isConfigured = false;
}

void RPCSeedPattern::setRecHits(const ConstMuonRecHitContainer& RecHits) {

    if(debug) cout << "setting RPCPattern recHits." << endl;
    for(ConstMuonRecHitContainer::const_iterator Iter = RecHits.begin(); Iter != RecHits.end(); Iter++) {
        DetId theDetId = (*Iter)->geographicalId();
        const GeomDetUnit* theDetUnit = (*Iter)->detUnit();
        GeomDetEnumerators::SubDetector subDet = theDetUnit->subDetector();
        int Layer = -1;
        if(subDet == GeomDetEnumerators::RPCBarrel) {
            int RPCStation = RPCDetId(theDetId).station(); // For Barrel 1-4
            int RPCLayer = RPCDetId(theDetId).layer(); // Only for RB1/2 the layer is 1/2 for in/out
            Layer = (RPCStation > 2) ? (RPCStation + 1) : (RPCStation * 2 + RPCLayer - 3);
        }
        if(subDet == GeomDetEnumerators::RPCEndcap) {
            int RPCRegion = RPCDetId(theDetId).region();
            int RPCStation = RPCDetId(theDetId).station();
            int LayerShift = (RPCRegion == -1) ? NegativeEndcapLayerShift : PositiveEndcapLayerShift;
            Layer = RPCStation + LayerShift;  // For Endcap 1-3/4/5(while RE4 complete/double layer in RE2)
        }

        if(debug) cout << "RecHitLayer: " << Layer << endl;    

        if(Layer == -1)
            continue;

        theRecHits.push_back(*Iter);
        theRecHitLayers.push_back(Layer);
    }
    isRecHitset = true;
}

void RPCSeedPattern::clear() {
    theRecHits.clear();
    theRecHitLayers.clear();
    isRecHitset = false;
}


RPCSeedPattern::WeightedTrajectorySeed RPCSeedPattern::seed(const edm::EventSetup& eSetup, int& isGoodSeed) {

    if(isConfigured == false || isRecHitset == false) {
        if(debug) cout << "Configuration not set yet: " << isConfigured << ", " << isRecHitset << endl;
        return createFakeSeed(isGoodSeed);
    }

    // Which kind of magnetic field should get?
    eSetup.get<IdealMagneticFieldRecord>().get(theMagneticField);
    measureRecHitandMagneticField();

    isPatternChecked = false;
    isGoodPattern = -1;

    Algorithm = checkAlgorithm();
    if(Algorithm == 1) {
        DoubleSegmentPattern();
        checkDoubleSegmentPattern();
    }
    // For barrel singal segment pattern
    if(Algorithm == 2) { 
        SingleSegmentPattern();
        checkSingleSegmentPattern();
    }
    // For endcap singal segment pattern
    if(Algorithm == 3) {
        //SingleSegmentPattern(2);
        //checkSingleSegmentPattern(2);
    }
    return createSeed(isGoodSeed);
}


void RPCSeedPattern::measureRecHitandMagneticField() {
    // Get distance and magnetice field sampling information, recHit's position is not the border of Chamber and Iron
    sampleMagneticField.clear();
    DistanceXY = 0;
    DistanceZ = 0;
    IntervalMagneticFlux.clear();
    for(ConstMuonRecHitContainer::const_iterator RecHitIter = theRecHits.begin(); RecHitIter != (theRecHits.end()-1); RecHitIter++) {
        GlobalPoint FirstPosition = (*RecHitIter)->globalPosition();
        GlobalPoint LastPosition = (*(RecHitIter+1))->globalPosition();
        DistanceXY += ((GlobalVector)(LastPosition - FirstPosition)).perp();
        DistanceZ += ((GlobalVector)(LastPosition - FirstPosition)).z();
        GlobalPoint* samplePosition = new GlobalPoint[sampleCount];
        double dX = (LastPosition.x() - FirstPosition.x()) / (sampleCount + 1);
        double dY = (LastPosition.y() - FirstPosition.y()) / (sampleCount + 1);
        double dZ = (LastPosition.z() - FirstPosition.z()) / (sampleCount + 1);
        bool isLowMagneticField = true;
        for(unsigned int index = 0; index < sampleCount; index++) {
            samplePosition[index] = GlobalPoint((FirstPosition.x()+dX*(index+1)), (FirstPosition.y()+dY*(index+1)), (FirstPosition.z()+dZ*(index+1)));
            GlobalVector TempMagneticField = theMagneticField->inTesla(samplePosition[index]);
            // for endcap the magnetic field could be complex so use mag() intead of z() while checking
            if(fabs(TempMagneticField.z()) > MagnecticFieldThreshold)
                isLowMagneticField = false;
            if(debug) cout << "Sampling magnetic field : " << TempMagneticField << endl;
            sampleMagneticField.push_back(TempMagneticField);
        }
        if(debug) cout << "IntervalMagneticFlux(0-high,1-low): " << isLowMagneticField << endl;
        IntervalMagneticFlux.push_back(isLowMagneticField);
        delete [] samplePosition;
    }
    int IntervalIndex = findIntervalIndex();
    if(IntervalIndex >= 0)
        MeanMagneticField = getMeanMagneticField(IntervalIndex);
    else 
        MeanMagneticField = GlobalVector(0, 0, 0);

    if(debug) cout << "MeanMagneticField: " << MeanMagneticField << ". DistanceXY: " << DistanceXY << ", DistanceZ: " << DistanceZ << endl;
}


int RPCSeedPattern::findIntervalIndex() {
    int IntervalIndex = -1;
    bool findIntervalIndex = false;
    // find the 1st interval with low magnetic field in full range, which lead to a straight segment. Then check if the next interval are with high magnetic field, where the track start the bending.
    for(int Index = 0; Index < (int)IntervalMagneticFlux.size(); Index++) {
        if(IntervalMagneticFlux[Index] == true && findIntervalIndex == false)
            IntervalIndex = Index;
        if(IntervalMagneticFlux[Index] == false && IntervalIndex != -1 && IntervalIndex == (Index-1) && findIntervalIndex == false) {
            IntervalIndex = Index;
            findIntervalIndex = true;
        }
    }
    if(findIntervalIndex == false)
        IntervalIndex = -1;

    if(debug) cout << "Find IntervalIndex: " << IntervalIndex << endl;

    // The ref point index
    RefIndex = IntervalIndex;

    return IntervalIndex;
}

GlobalVector RPCSeedPattern::getMeanMagneticField(const int IntervalIndex) {
    GlobalVector theMagneticField(0, 0, 0);
    unsigned int sampleNumber = 0;
    if(IntervalIndex < 0 || sampleCount <= 0)
        return theMagneticField;
    for(unsigned int i = IntervalIndex*sampleCount; i < (IntervalIndex+1)*sampleCount; i++) {
        GlobalVector TempMagnaeticField = sampleMagneticField[i];
        if(TempMagnaeticField.mag() > MagnecticFieldThreshold) {
            theMagneticField += TempMagnaeticField;
            sampleNumber++;
        }
    }
    return theMagneticField / sampleNumber;
}

int RPCSeedPattern::checkAlgorithm() {
    std::vector<int> AlgorithmChoice;
    AlgorithmChoice.clear();
    bool isBarrel =false;
    bool isPositiveEndcap = false;
    bool isNegativeEndcap = false;
    int theBarrelOccupancyCode = 0;
    int theEndcapOccupancyCode = 0;
    for(std::vector<unsigned int>::const_iterator LayerIter = theRecHitLayers.begin(); LayerIter != theRecHitLayers.end(); LayerIter++) {
        int RecHitLayer = *LayerIter;
        if(RecHitLayer < BarrelLayerNumber) {
            isBarrel = true;
            theBarrelOccupancyCode |= (int)pow(2, RecHitLayer);
        }
        else {
            if(RecHitLayer < (BarrelLayerNumber + EachEndcapLayerNumber) && RecHitLayer >=BarrelLayerNumber) {
                isNegativeEndcap = true;
                theEndcapOccupancyCode |= (int)pow(2, RecHitLayer-BarrelLayerNumber);
            }
            if(RecHitLayer >=(BarrelLayerNumber + EachEndcapLayerNumber) && RecHitLayer <(BarrelLayerNumber + EachEndcapLayerNumber * 2)) {
                isPositiveEndcap = true;
                theEndcapOccupancyCode |= (int)pow(2, RecHitLayer-BarrelLayerNumber-EachEndcapLayerNumber);
            }
        }
    }
    // Set algorithm choice with priority, 1-barrel_double_segment, 2-barrel_singal_segment, 3-endcap_singal_segment
    if(debug) cout << "theBarrelOccupancyCode: " << theBarrelOccupancyCode << ", theEndcapOccupancyCode: " << theEndcapOccupancyCode << endl;
    if(isBarrel == true && isNegativeEndcap == false && isPositiveEndcap == false && (theBarrelOccupancyCode & BarrelDoubleSegmentCode) == BarrelDoubleSegmentCode && (theBarrelOccupancyCode & BarrelDoubleSegmentOptionalCode) != 0)
        AlgorithmChoice.push_back(1);
    if(isBarrel == true && isNegativeEndcap == false && isPositiveEndcap == false && (theBarrelOccupancyCode & BarrelSingleSegmentCode) == BarrelSingleSegmentCode && (theBarrelOccupancyCode & BarrelSingleSegmentOptionalCode) != 0)
        AlgorithmChoice.push_back(2);
    if(isBarrel == false && isNegativeEndcap == true && isPositiveEndcap == false && (theEndcapOccupancyCode & EndcapSingleSegmentCode) == EndcapSingleSegmentCode && (theEndcapOccupancyCode & EndcapSingleSegmentOptionalCode) != 0)
        AlgorithmChoice.push_back(3);
    if(isBarrel == false && isNegativeEndcap == false && isPositiveEndcap == true && (theEndcapOccupancyCode & EndcapSingleSegmentCode) == EndcapSingleSegmentCode && (theEndcapOccupancyCode & EndcapSingleSegmentOptionalCode) != 0)
        AlgorithmChoice.push_back(3);

    // Auto choice or manual choise
    int FinalAlgorithm = -1;
    if(AlgorithmType == 0 && AlgorithmChoice.size() > 0) {
        if(checkParameters((unsigned int)AlgorithmChoice[0]))
            FinalAlgorithm = AlgorithmChoice[0];
    }
    else {
        for(unsigned int i = 0; i < AlgorithmChoice.size(); i++)
            if(AlgorithmType == AlgorithmChoice[i])
                if(checkParameters((unsigned int)AlgorithmType))
                    FinalAlgorithm = AlgorithmType;
    }
    cout << "Choose FinalAlgorithm: " << FinalAlgorithm << endl;
    return FinalAlgorithm;
}


bool RPCSeedPattern::checkParameters(unsigned int theAlgorithmType) {

    bool isParametersSet = true;
    if(BendingPhiLowerTH.size() < theAlgorithmType || BendingPhiUpperTH.size() < theAlgorithmType || MeanPt_Parameter0.size() < theAlgorithmType || MeanPt_Parameter1.size() < theAlgorithmType || MeanPt_Parameter2.size() < theAlgorithmType || SigmaPt_Parameter0.size() < theAlgorithmType || SigmaPt_Parameter1.size() < theAlgorithmType || SigmaPt_Parameter2.size() < theAlgorithmType )
        isParametersSet = false;

    return isParametersSet;
}


void RPCSeedPattern::DoubleSegmentPattern() {

    if(debug) cout << "estimating RPC Seed by Double Segment Pattern..." << endl;
    // Unset the pattern estimation signal
    isPatternChecked = false;
    isGoodPattern = 0;

    if(Algorithm != 1) {
        if(debug) cout << "Algorithm does not fit this pattern." << endl;
        isPatternChecked = true;
        isGoodPattern = -1;
        return;
    }

    // Form segments
    unsigned int Index = 0;
    unsigned int n  = 0;
    SegmentRB.clear();
    // reserve only prepair the reserve space but not the sie, it still needs the push/insert method to visit the element
    SegmentRB.resize(theRecHits.size()-2);
    for(ConstMuonRecHitContainer::const_iterator RecHitIter = theRecHits.begin(); RecHitIter != theRecHits.end(); RecHitIter++) {
        if(debug) cout << "Index " << Index << ",n "  << n << ", Size: " << SegmentRB.size() <<  endl;
        if((Index%2) == 0 && Index <= 3) {
            SegmentRB[n].first = (*RecHitIter);
        }
        if(debug) cout << "SegmentRB " << n << " recHit: " << (*RecHitIter)->globalPosition() << endl;

        if((Index%2) == 1 || Index > 3) {
            SegmentRB[n].second = (*RecHitIter);
            n++;
        }

        //if(debug) cout << "SegmentRB " << n << " recHit: " << (*RecHitIter)->globalPosition() << endl;

        Index++;
    }
    n--;
    for(; n >= 2; n--)
        SegmentRB[n].first = theRecHits[3];

    computePtwithDoubleSegment();
}


void RPCSeedPattern::computePtwithDoubleSegment() {
    if(SegmentRB.size() < 2) {
        if(debug) cout << "Segment not enough." << endl;
        BendingWise = 0;
        MeanPt = 0.;
        SigmaPt = 0.;
        Charge = 0;
        ZDirection = 0;
        Momentum = GlobalVector(0, 0, 0);
        isPatternChecked = true;
        isGoodPattern = -1;
        return;
    }

    GlobalVector SegmentVector0 = (GlobalVector)((SegmentRB[0].second)->globalPosition() - (SegmentRB[0].first)->globalPosition());
    GlobalVector SegmentVector1 = (GlobalVector)((SegmentRB[1].second)->globalPosition() - (SegmentRB[1].first)->globalPosition());
    double BendingPhi = (SegmentVector1.phi() - SegmentVector0.phi()).value();
    if(fabs(BendingPhi) < BendingPhiLowerTH[Algorithm-1]) {
        BendingWise = 0;
        Charge = 0;
        isGoodPattern = 0;
        isPatternChecked = true;
    }
    else {
        BendingWise = (BendingPhi > 0.) ? 1 : -1;
        Charge = BendingWise * (int)(fabs(MeanMagneticField.z()) / MeanMagneticField.z()) * -1;
    }

    if(fabs(BendingPhi) <= BendingPhiUpperTH[Algorithm-1] && fabs(BendingPhi) >= BendingPhiLowerTH[Algorithm-1]) {
        MeanPt = BendingPhi * BendingPhi * MeanPt_Parameter0[Algorithm-1] + fabs(BendingPhi) * MeanPt_Parameter1[Algorithm-1] + MeanPt_Parameter2[Algorithm-1];
        SigmaPt = BendingPhi * BendingPhi * SigmaPt_Parameter0[Algorithm-1] + fabs(BendingPhi) * SigmaPt_Parameter1[Algorithm-1] + SigmaPt_Parameter2[Algorithm-1];
        SigmaPt *= 2.;
    }
    else {
        isGoodPattern = 0;
        isPatternChecked = true;
        MeanPt = 0.;
        SigmaPt = 0.;
    }

    // Check the Z direction
    if(fabs((*(theRecHits.end()-1))->globalPosition().z() - (*(theRecHits.begin()))->globalPosition().z()) > ZError) {
        if(((*(theRecHits.end()-1))->globalPosition().z() - (*(theRecHits.begin()))->globalPosition().z()) > ZError)
            ZDirection = 1;
        else
            ZDirection = -1;
    }
    else
        ZDirection = 0;

    GlobalVector InitialPt = GlobalVector(SegmentVector0.x(), SegmentVector0.y(), 0);
    GlobalVector InitialMomentum = InitialPt.unit() * DistanceXY + GlobalVector(0, 0, DistanceZ);
    Momentum = InitialMomentum.unit() * InitialMomentum.mag() * MeanPt / InitialMomentum.perp();

    if(debug) cout << "BendingPhi: " << BendingPhi << ", MeanPt: " << MeanPt << "SigmaPt: " << SigmaPt << ", Momentum: " << Momentum << endl;
}


void RPCSeedPattern::checkDoubleSegmentPattern() {

    if(isPatternChecked == true)
        return;

    isGoodPattern = 1;
    GlobalVector SegmentVector1 = (GlobalVector)((SegmentRB[1].second)->globalPosition() - (SegmentRB[1].first)->globalPosition());
    unsigned int n = SegmentRB.size() - 1;
    for(; n >= 2; n--) {
        GlobalVector SegmentVectorFilter2 = (GlobalVector)((SegmentRB[n].second)->globalPosition() - (SegmentRB[n].first)->globalPosition());
        if(debug) cout << "SegmentVector1: " << (SegmentRB[1].first)->globalPosition() << ", " << (SegmentRB[1].second)->globalPosition() << ". SegmentVectorFilter2: " << (SegmentRB[n].first)->globalPosition() << ", " << (SegmentRB[n].second)->globalPosition() << ". SegmentVector1: " << SegmentVector1.phi() << ", SegmentVectorFilter2: " << SegmentVectorFilter2.phi() << endl;
        double BendingPhiFilter2 = (SegmentVectorFilter2.phi() - SegmentVector1.phi()).value();
        if(debug) cout << "BendingPhiFilter2: " << BendingPhiFilter2 << endl;
        if((BendingPhiFilter2 * (double)BendingWise) <= 0.) {
            if(debug) cout << "Bad pattern: bendingPhi of filter2 does not match." << endl;
            isGoodPattern = 0;
        }
    }

    // Check the Z direction
    if(debug) cout << "Check ZDirection is :" << ZDirection << endl;
    for(ConstMuonRecHitContainer::const_iterator iter = theRecHits.begin(); iter != (theRecHits.end()-1); iter++) {
        if(ZDirection == 0) {
            if(fabs((*(iter+1))->globalPosition().z()-(*iter)->globalPosition().z()) > ZError) {
                if(debug) cout << "Pattern find error in Z direction: wrong perpendicular direction" << endl;
                isGoodPattern = 0;
            }
        }
        else {
            if((int)(((*(iter+1))->globalPosition().z()-(*iter)->globalPosition().z())/ZError)*ZDirection < 0) {
                if(debug) cout << "Pattern find error in Z direction: wrong Z direction" << endl;
                isGoodPattern = 0;
            }
        }
    }

    isPatternChecked = true;
}

void RPCSeedPattern::SingleSegmentPattern() {

    if(debug) cout << "estimating RPC Seed by Double Segment Pattern..." << endl;
    // Unset the pattern estimation signal
    isPatternChecked = false;
    isGoodPattern = 0;

    if(Algorithm != 2) {
        if(debug) cout << "Algorithm does not fit this pattern." << endl;
        isPatternChecked = true;
        isGoodPattern = -1;
        return;
    }

    // Form segments
    unsigned int Index = 0;
    unsigned int n  = 0;
    SegmentRB.clear();
    // reserve only prepair the reserve space but not the sie, it still needs the push/insert method to visit the element
    if(Algorithm == 2) {
        SegmentRB.resize(theRecHits.size()-1);
        for(ConstMuonRecHitContainer::const_iterator RecHitIter = theRecHits.begin(); RecHitIter != theRecHits.end(); RecHitIter++) {
            if(debug) cout << "Index " << Index << ",n "  << n << ", Size: " << SegmentRB.size() <<  endl;
            if((Index%2) == 0 && Index <= 0) {
                SegmentRB[n].first = (*RecHitIter);
            }
            if(debug) cout << "SegmentRB " << n << " recHit: " << (*RecHitIter)->globalPosition() << endl;

            if((Index%2) == 1 || Index > 0) {
                SegmentRB[n].second = (*RecHitIter);
                n++;
            }

            //if(debug) cout << "SegmentRB " << n << " recHit: " << (*RecHitIter)->globalPosition() << endl;

            Index++;
        }
        n--;
        for(; n >= 1; n--)
            SegmentRB[n].first = theRecHits[1];
    }

    computePtwithSingleSegment();
}


void RPCSeedPattern::computePtwithSingleSegment() {
    if(SegmentRB.size() < 2) {
        if(debug) cout << "Segment not enough." << endl;
        BendingWise = 0;
        MeanPt = 0.;
        SigmaPt = 0.;
        Charge = 0;
        ZDirection = 0;
        Momentum = GlobalVector(0, 0, 0);
        isPatternChecked = true;
        isGoodPattern = -1;
        return;
    }

    GlobalVector SegmentVector0 = (GlobalVector)((SegmentRB[0].second)->globalPosition() - (SegmentRB[0].first)->globalPosition());
    GlobalVector SegmentVector1 = (GlobalVector)((SegmentRB[1].second)->globalPosition() - (SegmentRB[1].first)->globalPosition());
    double BendingPhi = (SegmentVector1.phi() - SegmentVector0.phi()).value();
    if(fabs(BendingPhi) < BendingPhiLowerTH[Algorithm-1]) {
        BendingWise = 0;
        Charge = 0;
        isGoodPattern = 0;
        isPatternChecked = true;
    }
    else {
        BendingWise = (BendingPhi > 0.) ? 1 : -1;
        Charge = BendingWise * (int)(fabs(MeanMagneticField.z()) / MeanMagneticField.z()) * -1;
    }

    if(fabs(BendingPhi) <= BendingPhiUpperTH[Algorithm-1] && fabs(BendingPhi) >= BendingPhiLowerTH[Algorithm-1]) {
        MeanPt = BendingPhi * BendingPhi * MeanPt_Parameter0[Algorithm-1] + fabs(BendingPhi) * MeanPt_Parameter1[Algorithm-1] + MeanPt_Parameter2[Algorithm-1];
        SigmaPt = BendingPhi * BendingPhi * SigmaPt_Parameter0[Algorithm-1] + fabs(BendingPhi) * SigmaPt_Parameter1[Algorithm-1] + SigmaPt_Parameter2[Algorithm-1];
        SigmaPt *= 2.;
    }
    else {
        isGoodPattern = 0;
        isPatternChecked = true;
        MeanPt = 0.;
        SigmaPt = 0.;
    }

    // Check the Z direction
    if(fabs((*(theRecHits.end()-1))->globalPosition().z() - (*(theRecHits.begin()))->globalPosition().z()) > ZError) {
        if(((*(theRecHits.end()-1))->globalPosition().z() - (*(theRecHits.begin()))->globalPosition().z()) > ZError)
            ZDirection = 1;
        else
            ZDirection = -1;
    }
    else
        ZDirection = 0;

    GlobalVector InitialPt = GlobalVector(SegmentVector0.x(), SegmentVector0.y(), 0);
    GlobalVector InitialMomentum = InitialPt.unit() * DistanceXY + GlobalVector(0, 0, DistanceZ);
    Momentum = InitialMomentum.unit() * InitialMomentum.mag() * MeanPt / InitialMomentum.perp();

    if(debug) cout << "BendingPhi: " << BendingPhi << ", MeanPt: " << MeanPt << "SigmaPt: " << SigmaPt << ", Momentum: " << Momentum << endl;
}


void RPCSeedPattern::checkSingleSegmentPattern() {

    if(isPatternChecked == true)
        return; 

    isGoodPattern = 1;
    GlobalVector SegmentVector1 = (GlobalVector)((SegmentRB[1].second)->globalPosition() - (SegmentRB[1].first)->globalPosition());
    unsigned int n = SegmentRB.size() - 1;
    for(; n >= 2; n--) {
        GlobalVector SegmentVectorFilter2 = (GlobalVector)((SegmentRB[n].second)->globalPosition() - (SegmentRB[n].first)->globalPosition());
        if(debug) cout << "SegmentVector1: " << (SegmentRB[1].first)->globalPosition() << ", " << (SegmentRB[1].second)->globalPosition() << ". SegmentVectorFilter2: " << (SegmentRB[n].first)->globalPosition() << ", " << (SegmentRB[n].second)->globalPosition() << ". SegmentVector1: " << SegmentVector1.phi() << ", SegmentVectorFilter2: " << SegmentVectorFilter2.phi() << endl;
        double BendingPhiFilter2 = (SegmentVectorFilter2.phi() - SegmentVector1.phi()).value();
        if(debug) cout << "BendingPhiFilter2: " << BendingPhiFilter2 << endl;
        if((BendingPhiFilter2 * (double)BendingWise) <= 0.) {
            if(debug) cout << "Bad pattern: bendingPhi of filter2 does not match." << endl;
            isGoodPattern = 0;
        }   
    }   
    
    // Check the Z direction
    if(debug) cout << "Check ZDirection is :" << ZDirection << endl;
    for(ConstMuonRecHitContainer::const_iterator iter = theRecHits.begin(); iter != (theRecHits.end()-1); iter++) {
        if(ZDirection == 0) {
            if(fabs((*(iter+1))->globalPosition().z()-(*iter)->globalPosition().z()) > ZError) {
                if(debug) cout << "Pattern find error in Z direction: wrong perpendicular direction" << endl;
                isGoodPattern = 0;
            }   
        }   
        else {
            if((int)(((*(iter+1))->globalPosition().z()-(*iter)->globalPosition().z())/ZError)*ZDirection < 0) {
                if(debug) cout << "Pattern find error in Z direction: wrong Z direction" << endl;
                isGoodPattern = 0;
            }   
        }   
    }   

    isPatternChecked = true;
}


RPCSeedPattern::WeightedTrajectorySeed RPCSeedPattern::createFakeSeed(int& isGoodSeed) {
    // Create a fake seed and return
    if(debug) cout << "Now create a fake seed" << endl;
    isPatternChecked = true;
    isGoodPattern = -1;
    Charge = 0;
    Momentum = GlobalVector(0., 0., 0.);

    const ConstMuonRecHitPointer theBestRefRecHit = BestRefRecHit();

    LocalPoint RefPosition = theBestRefRecHit->localPosition();
    LocalVector RefMomentum = theBestRefRecHit->det()->toLocal(Momentum);
    LocalTrajectoryParameters theLTP(RefPosition, RefMomentum, Charge);

    AlgebraicSymMatrix theErrorMatrix(5,0);
    theErrorMatrix = theBestRefRecHit->parametersError().similarityT(theBestRefRecHit->projectionMatrix());
    theErrorMatrix[0][0] = 0.;
    LocalTrajectoryError theLTE(theErrorMatrix);

    TrajectoryStateOnSurface theTSOS(theLTP, theLTE, theBestRefRecHit->det()->surface(), &*theMagneticField);

    DetId theDetId = theBestRefRecHit->geographicalId();
    TrajectoryStateTransform theTST;
    PTrajectoryStateOnDet *seedTSOS = theTST.persistentState(theTSOS, theDetId.rawId());

    edm::OwnVector<TrackingRecHit> container;
    for(ConstMuonRecHitContainer::const_iterator iter=theRecHits.begin(); iter!=theRecHits.end(); iter++)
        container.push_back((*iter)->hit()->clone());

    TrajectorySeed theSeed(*seedTSOS, container, alongMomentum);
    WeightedTrajectorySeed theWeightedSeed;
    theWeightedSeed.first = theSeed;
    theWeightedSeed.second = PatternQuality;
    isGoodSeed = isGoodPattern;

    delete seedTSOS;
    return theWeightedSeed;
}

RPCSeedPattern::WeightedTrajectorySeed RPCSeedPattern::createSeed(int& isGoodSeed) {

    if(isPatternChecked == false || isGoodPattern == -1) {
        cout <<"Pattern is not good: " << isPatternChecked << ", " << isGoodPattern << ". Create a fake seed instead!" << endl;
        return createFakeSeed(isGoodSeed);
    }

    // Get the reference recHit, DON'T use the recHit on 1st layer(inner most layer)
    const ConstMuonRecHitPointer theBestRefRecHit = BestRefRecHit();

    LocalPoint RefPosition = theBestRefRecHit->localPosition();
    LocalVector RefMomentum = theBestRefRecHit->det()->toLocal(Momentum);
    LocalTrajectoryParameters theLTP(RefPosition, RefMomentum, Charge);

    LocalTrajectoryError theLTE = getErrorMatrix(theBestRefRecHit);

    TrajectoryStateOnSurface theTSOS(theLTP, theLTE, theBestRefRecHit->det()->surface(), &*theMagneticField);
    DetId theDetId = theBestRefRecHit->geographicalId();
    TrajectoryStateTransform theTST;
    PTrajectoryStateOnDet *seedTSOS = theTST.persistentState(theTSOS, theDetId.rawId());

    edm::OwnVector<TrackingRecHit> container;
    for(ConstMuonRecHitContainer::const_iterator iter=theRecHits.begin(); iter!=theRecHits.end(); iter++) {
        // This casting withou clone will cause memory overflow when used in push_back
        // Since container's deconstructor functiion free the pointer menber!
        //TrackingRecHit* pt = dynamic_cast<TrackingRecHit*>(&*(*iter));
        //cout << "Push recHit type " << pt->getType() << endl;
        container.push_back((*iter)->hit()->clone());
    }

    TrajectorySeed theSeed(*seedTSOS, container, alongMomentum);
    WeightedTrajectorySeed theWeightedSeed;
    theWeightedSeed.first = theSeed;
    theWeightedSeed.second = PatternQuality;
    isGoodSeed = isGoodPattern;

    delete seedTSOS;
    return theWeightedSeed;
}

LocalTrajectoryError RPCSeedPattern::getErrorMatrix(const ConstMuonRecHitPointer& theRefRecHit) {

    // for possable lack of error matrix element[i][j], we use 5 parameter to construct theLTE;
    double dX = 0;
    double dY = 0;
    double dXdZ = 0;
    double dYdZ = 0;
    double dPInv = 0;

    // trackingrecHit from RPC should have projectionMatrix of 5*1 with only P[3]=1 for dx error, but not sure for that, so we extract dX and dY from the error matrix for construct local trajectory error.
    // projectionMatrix is used to project the track's full error vector (1/dP, dXdZ, dYdZ, dX, dY) the sub detector's error dimension(<5), since some detector does not have ability for direction/2D measurement, e.g. CSCSegment only provide direction and position, which has error like (dXdZ, dYdZ, dX, dY), so it's projectionMatrix is a 4*5 matrix which project 5D track error vector to its onw 4D space. Formula should be DetErrorVec=P*TrackErrorVec, while P is the projectionMatrix from track measurement space to det measurement space;
    // while the errorMatrix is the covariance matrix for error parameter, which is under sub detector's measuremant space(<5D, n*n matrix, n<5), so the covariance matrix for track measurement space(5D, 5*5 matrix) should be TrackError=P-1*DetError*P
    AlgebraicSymMatrix theErrorMatrix(5, 0);
    theErrorMatrix = theRefRecHit->parametersError().similarityT(theRefRecHit->projectionMatrix());   
    dX = sqrt(theErrorMatrix[3][3]);
    dY = sqrt(theErrorMatrix[4][4]);


    LocalError theLocalError = (SegmentRB[0].first)->localPositionError();
    double dXatRef = sqrt(theLocalError.xx());
    double dYatRef = sqrt(theLocalError.yy());
    const GeomDetUnit* RefRPCRoll = (SegmentRB[0].second)->detUnit();
    LocalPoint RecHit0 = RefRPCRoll->toLocal((SegmentRB[0].first)->globalPosition());
    LocalPoint RecHit1 = RefRPCRoll->toLocal((SegmentRB[0].second)->globalPosition());
    LocalVector LocalSegment0 = (LocalVector)(RecHit1 - RecHit0);
    LocalVector LocalSegment1 = LocalVector(LocalSegment0.x()+2*dXatRef, LocalSegment0.y(), LocalSegment0.z());
    LocalVector LocalSegment2 = LocalVector(LocalSegment0.x()-2*dXatRef, LocalSegment0.y(), LocalSegment0.z());
    LocalVector LocalSegment3 = LocalVector(LocalSegment0.x(), LocalSegment0.y()+2*dYatRef, LocalSegment0.z());
    LocalVector LocalSegment4 = LocalVector(LocalSegment0.x(), LocalSegment0.y()-2*dYatRef, LocalSegment0.z());
    double dXdZ1 = LocalSegment1.x() / LocalSegment1.z() - LocalSegment0.x() / LocalSegment0.z();
    double dXdZ2 = LocalSegment2.x() / LocalSegment2.z() - LocalSegment0.x() / LocalSegment0.z();
    dXdZ = (fabs(dXdZ1) >= fabs(dXdZ2)) ? dXdZ1 : dXdZ2;
    double dYdZ1 = LocalSegment3.y() / LocalSegment3.z() - LocalSegment0.y() / LocalSegment0.z();
    double dYdZ2 = LocalSegment4.y() / LocalSegment4.z() - LocalSegment0.y() / LocalSegment0.z();
    dYdZ = (fabs(dYdZ1) >= fabs(dYdZ2)) ? dYdZ1 : dYdZ2;

    double SigmaP = SigmaPt * Momentum.mag() / MeanPt;
    dPInv = SigmaP / (Momentum.mag() * Momentum.mag());

    LocalTrajectoryError theLTE(dX, dY, dXdZ, dYdZ, dPInv);
    return theLTE;
}

RPCSeedPattern::ConstMuonRecHitPointer RPCSeedPattern::BestRefRecHit() const {

    if(RefIndex >= 0 && RefIndex < (int)theRecHits.size()-1)
        return theRecHits[RefIndex];
    else
        return theRecHits[0];
}
