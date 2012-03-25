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
    isVertexConstraint = iConfig.getParameter<bool>("isVertexConstraint");
    isContinuousFilter = iConfig.getParameter<bool>("isContinuousFilter");
    Cut1234 = iConfig.getParameter< vector<double> >("Cut1234");
    CutMax = iConfig.getParameter< vector<double> >("CutMax");
    BendingPhiLowerTH = iConfig.getParameter< vector<double> >("BendingPhiLowerTH");
    BendingPhiUpperTH = iConfig.getParameter< vector<double> >("BendingPhiUpperTH");
    BendingPhiFitValueUpperLimit = iConfig.getParameter< vector<double> >("BendingPhiFitValueUpperLimit");
    BendingPhiFitSigmaUpperLimit = iConfig.getParameter< vector<double> >("BendingPhiFitSigmaUpperLimit");
    MeanPt_Parameter0 = iConfig.getParameter< vector<double> >("MeanPt_Parameter0");
    MeanPt_Parameter1 = iConfig.getParameter< vector<double> >("MeanPt_Parameter1");
    MeanPt_Parameter2 = iConfig.getParameter< vector<double> >("MeanPt_Parameter2");
    SigmaPt_Parameter0 = iConfig.getParameter< vector<double> >("SigmaPt_Parameter0");
    SigmaPt_Parameter1 = iConfig.getParameter< vector<double> >("SigmaPt_Parameter1");
    SigmaPt_Parameter2 = iConfig.getParameter< vector<double> >("SigmaPt_Parameter2");
    isConfigured = true;

    if(isVertexConstraint == true && isContinuousFilter == false)
        applyFilter = false;
    else
        applyFilter = true;

    if(Cut1234.size() == 0 || CutMax.size() == 0 || BendingPhiFitValueUpperLimit.size() == 0 || BendingPhiFitSigmaUpperLimit.size() == 0 || BendingPhiLowerTH.size() == 0 || BendingPhiUpperTH.size() == 0 || MeanPt_Parameter0.size() == 0 || MeanPt_Parameter1.size() == 0 || MeanPt_Parameter2.size() == 0 || SigmaPt_Parameter0.size() == 0 || SigmaPt_Parameter1.size() == 0 || SigmaPt_Parameter2.size() == 0 )
        isConfigured = false;
}

void RPCSeedPattern::setRecHits(const ConstMuonRecHitContainer& RecHits) {

    if(debug) cout << "setting RPCPattern recHits." << endl;
    theRecHits.clear();
    theRecHitLayers.clear();
    for(ConstMuonRecHitContainer::const_iterator Iter = RecHits.begin(); Iter != RecHits.end(); Iter++) {
        DetId theDetId = (*Iter)->geographicalId();
        const GeomDetUnit* theDetUnit = (*Iter)->detUnit();
        GeomDetEnumerators::SubDetector subDet = theDetUnit->subDetector();
        int Layer = -1;
        int SeedLayer = -1;
        if(subDet == GeomDetEnumerators::RPCBarrel) {
            int RPCStation = RPCDetId(theDetId).station(); // For Barrel 1-4
            int RPCLayer = RPCDetId(theDetId).layer(); // Only for RB1/2 the layer is 1/2 for in/out
            Layer = (RPCStation > 2) ? (RPCStation + 1) : (RPCStation * 2 + RPCLayer - 3);
            SeedLayer = Layer;
        }
        if(subDet == GeomDetEnumerators::RPCEndcap) {
            int RPCRegion = RPCDetId(theDetId).region();
            int RPCStation = RPCDetId(theDetId).station();
            int LayerShift = (RPCRegion == -1) ? NegativeEndcapLayerShift : PositiveEndcapLayerShift;
            Layer = RPCStation + LayerShift;  // For Endcap 1-3/4/5(while RE4 complete/double layer in RE2)
            SeedLayer = RPCStation - 1;
        }

        if(debug) cout << "RecHitLayer: " << Layer << endl;    

        if(Layer == -1)
            continue;

        theRecHits.push_back(*Iter);
        theRecHitLayers.push_back(Layer);
        theRecHitPosition[SeedLayer] = (*Iter)->globalPosition();
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
    createPattern();
    if(Algorithm <= 4 && Algorithm >= 0)
        checkDoubleSegmentPattern();
    if(Algorithm >= 5 && Algorithm <= 12)
        checkSingleSegmentPattern();

    computeSegmentPattern();

    return createSeed(isGoodSeed);
}


void RPCSeedPattern::measureRecHitandMagneticField() {
    // Get distance and magnetice field sampling information, recHit's position is not the border of Chamber and Iron
    sampleMagneticField.clear();
    DistanceXY = 0;
    IntervalMagneticFlux.clear();
    for(ConstMuonRecHitContainer::const_iterator RecHitIter = theRecHits.begin(); RecHitIter != (theRecHits.end()-1); RecHitIter++) {
        GlobalPoint FirstPosition = (*RecHitIter)->globalPosition();
        GlobalPoint LastPosition = (*(RecHitIter+1))->globalPosition();
        DistanceXY += ((GlobalVector)(LastPosition - FirstPosition)).perp();
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
    if(IntervalIndex > 0) {
        RefIndex = IntervalIndex - 1;
        RefSegment.first = theRecHits[RefIndex];
        RefSegment.second = theRecHits[RefIndex+1];
        theRefRecHit = theRecHits[RefIndex];
        MeanMagneticField = getMeanMagneticField(IntervalIndex);
    }
    else {
        RefIndex = 0;
        RefSegment.first = theRecHits[0];
        RefSegment.second = theRecHits[1];
        theRefRecHit = theRecHits[0];
        MeanMagneticField = getMeanMagneticField(1);
    }

    if(isVertexConstraint == true)
        DistanceZ = theRecHits[theRecHits.size()-1]->globalPosition().z() * (theRecHits[theRecHits.size()-1]->globalPosition().perp()-theRecHits[0]->globalPosition().perp()) / theRecHits[theRecHits.size()-1]->globalPosition().perp();
    else
        DistanceZ = theRecHits[theRecHits.size()-1]->globalPosition().z() - theRecHits[0]->globalPosition().z();

    if(fabs(DistanceZ) > ZError) {
        if(DistanceZ > ZError)
            ZDirection = 1;
        else
            ZDirection = -1;
    }
    else
        ZDirection = 0;

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
    if(sampleNumber == 0)
        sampleNumber = 1.;
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
    for(std::vector<int>::const_iterator LayerIter = theRecHitLayers.begin(); LayerIter != theRecHitLayers.end(); LayerIter++) {
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
    // Set algorithm choice with priority, value/10=index_number_in_algorithm_parameter_vector, value%10=algorithm_type:1-barrel_double_segment, 2-barrel_singal_segment, 3-endcap_singal_segment
    if(debug) cout << "theBarrelOccupancyCode: " << theBarrelOccupancyCode << ", theEndcapOccupancyCode: " << theEndcapOccupancyCode << endl;
    if(isBarrel == true && isNegativeEndcap == false && isPositiveEndcap == false && theBarrelOccupancyCode == BarrelDoubleSegmentCode1)
        AlgorithmChoice.push_back(11);
    if(isBarrel == true && isNegativeEndcap == false && isPositiveEndcap == false && theBarrelOccupancyCode == BarrelDoubleSegmentCode2)
        AlgorithmChoice.push_back(21);
    if(isBarrel == true && isNegativeEndcap == false && isPositiveEndcap == false && theBarrelOccupancyCode == BarrelSingleSegmentCode1)
        AlgorithmChoice.push_back(32);
    if(isBarrel == true && isNegativeEndcap == false && isPositiveEndcap == false && (theBarrelOccupancyCode & BarrelSingleSegmentOptionalCode1) == BarrelSingleSegmentOptionalCode1 && (theBarrelOccupancyCode & BarrelSingleSegmentOptionalCode2) != 0 && (theBarrelOccupancyCode & BarrelSingleSegmentOptionalCode4) == BarrelSingleSegmentOptionalCode4)
        AlgorithmChoice.push_back(42);
    if(isBarrel == true && isNegativeEndcap == false && isPositiveEndcap == false && (theBarrelOccupancyCode & BarrelSingleSegmentOptionalCode1) == BarrelSingleSegmentOptionalCode1 && (theBarrelOccupancyCode & BarrelSingleSegmentOptionalCode2) != 0 && (theBarrelOccupancyCode & BarrelSingleSegmentOptionalCode3) == BarrelSingleSegmentOptionalCode3)
        AlgorithmChoice.push_back(52);
    if(isBarrel == true && isNegativeEndcap == false && isPositiveEndcap == false && theBarrelOccupancyCode == BarrelSingleSegmentCode2)
        AlgorithmChoice.push_back(62);
    if(isBarrel == false && isNegativeEndcap == true && isPositiveEndcap == false && theEndcapOccupancyCode == EndcapSingleSegmentCode1)
        AlgorithmChoice.push_back(73);
    if(isBarrel == false && isNegativeEndcap == false && isPositiveEndcap == true && theEndcapOccupancyCode == EndcapSingleSegmentCode1)
        AlgorithmChoice.push_back(83);

    // Auto choice or manual choise
    int FinalAlgorithm = -1;
    if(AlgorithmType == 0 && AlgorithmChoice.size() > 0) {
        if(checkParameters((unsigned int)(AlgorithmChoice[0]/10)))
            FinalAlgorithm = (int)AlgorithmChoice[0]/10;
    }
    else {
        for(unsigned int i = 0; i < AlgorithmChoice.size(); i++)
            if(AlgorithmType == (int)(AlgorithmChoice[i]%10))
                if(checkParameters((unsigned int)(AlgorithmChoice[i]/10)))
                    FinalAlgorithm = (int)AlgorithmChoice[i]/10;
    }

    PatternQuality = FinalAlgorithm;
    if(isVertexConstraint == false)
        FinalAlgorithm = FinalAlgorithm * 2 - 1;
    else
        FinalAlgorithm = FinalAlgorithm * 2;

    cout << "Choose FinalAlgorithm: " << FinalAlgorithm << endl;
    return FinalAlgorithm;
}

bool RPCSeedPattern::checkParameters(unsigned int theAlgorithmType) {

    unsigned int AlgorithmIndex;
    if(isVertexConstraint == false)
        AlgorithmIndex = theAlgorithmType * 2 - 1;
    else
        AlgorithmIndex = theAlgorithmType * 2;

    bool isParametersSet = true;
    if(Cut1234.size() < AlgorithmIndex || CutMax.size() < AlgorithmIndex || BendingPhiFitValueUpperLimit.size() < AlgorithmIndex || BendingPhiFitSigmaUpperLimit.size() < AlgorithmIndex || BendingPhiLowerTH.size() < AlgorithmIndex || BendingPhiUpperTH.size() < AlgorithmIndex || MeanPt_Parameter0.size() < AlgorithmIndex || MeanPt_Parameter1.size() < AlgorithmIndex || MeanPt_Parameter2.size() < AlgorithmIndex || SigmaPt_Parameter0.size() < AlgorithmIndex || SigmaPt_Parameter1.size() < AlgorithmIndex || SigmaPt_Parameter2.size() < AlgorithmIndex )
        isParametersSet = false;

    return isParametersSet;
}

void RPCSeedPattern::createPattern() {

    BendingFilter.clear();
    if(Algorithm == 1) {
        BendingFilter.push_back(BendingPhiIndexType(0,1,2,3));
        BendingFilter.push_back(BendingPhiIndexType(0,1,1,5));
        BendingFilter.push_back(BendingPhiIndexType(2,3,3,5));
    }
    if(Algorithm == 2) {
        BendingFilter.push_back(BendingPhiIndexType(0,1,2,3));
        BendingFilter.push_back(BendingPhiIndexType(0,1,1,5));
        BendingFilter.push_back(BendingPhiIndexType(2,3,3,5));
        BendingFilter.push_back(BendingPhiIndexType(2,2,2,3));
        BendingFilter.push_back(BendingPhiIndexType(0,0,0,1));
    }
    if(Algorithm == 3) {
        BendingFilter.push_back(BendingPhiIndexType(0,1,2,3));
        BendingFilter.push_back(BendingPhiIndexType(0,1,1,4));
        BendingFilter.push_back(BendingPhiIndexType(2,3,3,4));
    }
    if(Algorithm == 4) {
        BendingFilter.push_back(BendingPhiIndexType(0,1,2,3));
        BendingFilter.push_back(BendingPhiIndexType(0,1,1,4));
        BendingFilter.push_back(BendingPhiIndexType(2,3,3,4));
        BendingFilter.push_back(BendingPhiIndexType(2,2,2,3));
        BendingFilter.push_back(BendingPhiIndexType(0,0,0,1));
    }
    if(Algorithm == 5) {
        BendingFilter.push_back(BendingPhiIndexType(0,1,1,4));
        BendingFilter.push_back(BendingPhiIndexType(0,1,1,5));
    }
    if(Algorithm == 6) {
        BendingFilter.push_back(BendingPhiIndexType(0,1,1,4));
        BendingFilter.push_back(BendingPhiIndexType(0,1,1,5));
        BendingFilter.push_back(BendingPhiIndexType(0,0,0,1));
    }
    if(Algorithm == 7) {
        if(find(theRecHitLayers.begin(), theRecHitLayers.end(), 3) == theRecHitLayers.end())
            BendingFilter.push_back(BendingPhiIndexType(0,1,1,2)); // we take both layer2/3 but layer2 is worse
        else
            BendingFilter.push_back(BendingPhiIndexType(0,1,1,3));
        BendingFilter.push_back(BendingPhiIndexType(0,1,1,5));
    }
    if(Algorithm == 8) {
        if(find(theRecHitLayers.begin(), theRecHitLayers.end(), 3) == theRecHitLayers.end())
            BendingFilter.push_back(BendingPhiIndexType(0,1,1,2)); // we take both layer2/3 but layer2 is worse
        else
            BendingFilter.push_back(BendingPhiIndexType(0,1,1,3));
        BendingFilter.push_back(BendingPhiIndexType(0,1,1,5));
        BendingFilter.push_back(BendingPhiIndexType(0,0,0,1));
    }
    if(Algorithm == 9) {
        if(find(theRecHitLayers.begin(), theRecHitLayers.end(), 3) == theRecHitLayers.end())
            BendingFilter.push_back(BendingPhiIndexType(0,1,1,2)); // we take both layer2/3 but layer2 is worse
        else
            BendingFilter.push_back(BendingPhiIndexType(0,1,1,3));
        BendingFilter.push_back(BendingPhiIndexType(0,1,1,4));
    }
    if(Algorithm == 10) {
        if(find(theRecHitLayers.begin(), theRecHitLayers.end(), 3) == theRecHitLayers.end())
            BendingFilter.push_back(BendingPhiIndexType(0,1,1,2)); // we take both layer2/3 but layer2 is worse
        else
            BendingFilter.push_back(BendingPhiIndexType(0,1,1,3));
        BendingFilter.push_back(BendingPhiIndexType(0,1,1,4));
        BendingFilter.push_back(BendingPhiIndexType(0,0,0,1));
    }
    if(Algorithm == 11) {
        BendingFilter.push_back(BendingPhiIndexType(2,3,3,4));
        BendingFilter.push_back(BendingPhiIndexType(2,3,3,5));
    }
    if(Algorithm == 12) {
        BendingFilter.push_back(BendingPhiIndexType(2,3,3,4));
        BendingFilter.push_back(BendingPhiIndexType(2,3,3,5));
        BendingFilter.push_back(BendingPhiIndexType(2,2,2,3));
    }

    BendingPhiCollection.clear();
    for(unsigned int Index = 0; Index < BendingFilter.size(); Index++) {
        int i = BendingFilter[Index].m[0];
        int j = BendingFilter[Index].m[1];
        int k = BendingFilter[Index].m[2];
        int l = BendingFilter[Index].m[3];

        Geom::Phi<float> PhiI2J;
        Geom::Phi<float> PhiK2L;
        if(i != j)
            PhiI2J = ((GlobalVector)(theRecHitPosition[j] - theRecHitPosition[i])).phi();
        else
            PhiI2J = ((GlobalVector)(theRecHitPosition[j] - GlobalPoint(0,0,0))).phi();
        PhiK2L = ((GlobalVector)(theRecHitPosition[l] - theRecHitPosition[k])).phi();
        double TempBendingPhi = (PhiK2L - PhiI2J).value();
        // vertex bendingPhi is reverse w.r.t normal case
        if(i == j)
            TempBendingPhi *= -1.;
        BendingPhiCollection.push_back(TempBendingPhi);
    }
    BendingPhiMax = findMaxBendingPhi();

    // set the signal
    isGoodPattern = -1;
    isPatternChecked = false;
}

void RPCSeedPattern::checkDoubleSegmentPattern() {
    if(isPatternChecked == true)
        return;

    if(debug) cout << "checkDoubleSegmentPattern." << endl;

    isGoodPattern = 1;

    if(BendingPhiCollection.size() < 3 || Algorithm <= 0 || Algorithm > 4)
        isGoodPattern = -1;
    else {
        if(debug) cout << "BendingPhiCollection[0]: " << BendingPhiCollection[0] << ", BendingPhiMax: " << BendingPhiMax << endl;

        if(fabs(BendingPhiCollection[0]) < Cut1234[Algorithm-1] && fabs(BendingPhiMax) < CutMax[Algorithm-1])
            isGoodPattern = -1;
        if(applyFilter == true) 
            if(BendingPhiCollection[0]*BendingPhiCollection[2] < 0.)
                isGoodPattern = -1;
        
        // Check the Z direction
        if(debug) cout << "Check ZDirection is :" << ZDirection << endl;
        for(ConstMuonRecHitContainer::const_iterator iter = theRecHits.begin(); iter != (theRecHits.end()-1); iter++) {
            if(ZDirection == 0) {
                if(fabs((*(iter+1))->globalPosition().z()-(*iter)->globalPosition().z()) > ZError) {
                    if(debug) cout << "Pattern find error in Z direction: wrong perpendicular direction" << endl;
                    isGoodPattern = -1;
                }
            }
            else {
                if((int)(((*(iter+1))->globalPosition().z()-(*iter)->globalPosition().z())/ZError)*ZDirection < 0) {
                    if(debug) cout << "Pattern find error in Z direction: wrong Z direction" << endl;
                    isGoodPattern = -1;
                }
            }
        }
    }
    isPatternChecked = true;
}

void RPCSeedPattern::checkSingleSegmentPattern() {
    if(isPatternChecked == true)
        return;

    if(debug) cout << "checkSingleSegmentPattern." << endl;
    isGoodPattern = 1;

    if(BendingPhiCollection.size() < 2 || Algorithm <= 4 || Algorithm > 12)
        isGoodPattern = -1;
    else {
        if(debug) cout << "BendingPhiCollection[0]: " << BendingPhiCollection[0] << ", BendingPhiMax: " << BendingPhiMax
            << endl;

        if(fabs(BendingPhiCollection[0]) < Cut1234[Algorithm-1] && fabs(BendingPhiMax) < CutMax[Algorithm-1])
            isGoodPattern = -1;
        if(applyFilter == true)
            if(BendingPhiCollection[0]*(BendingPhiCollection[1]-BendingPhiCollection[0]) < 0.)
                isGoodPattern = -1;

        // Check the Z direction
        if(debug) cout << "Check ZDirection is :" << ZDirection << endl;
        for(ConstMuonRecHitContainer::const_iterator iter = theRecHits.begin(); iter != (theRecHits.end()-1); iter++) {
            if(ZDirection == 0) {
                if(fabs((*(iter+1))->globalPosition().z()-(*iter)->globalPosition().z()) > ZError) {
                    if(debug) cout << "Pattern find error in Z direction: wrong perpendicular direction" << endl;
                    isGoodPattern = -1;
                }
            }
            else {
                if((int)(((*(iter+1))->globalPosition().z()-(*iter)->globalPosition().z())/ZError)*ZDirection < 0) {
                    if(debug) cout << "Pattern find error in Z direction: wrong Z direction" << endl;
                    isGoodPattern = -1;
                }
            }
        }
    }
    isPatternChecked = true;
}

void RPCSeedPattern::computeSegmentPattern() {
    if(isGoodPattern < 0 || isPatternChecked == false) {
        if(debug) cout << "Pattern not pass filter." << endl;
        BendingWise = 0;
        MeanPt = 0.;
        SigmaPt = 0.;
        Charge = 0;
        Momentum = GlobalVector(0, 0, 0);
        isGoodPattern = -1;
        return;
    }

    isGoodPattern = 1;
    //double BendingPhi = findMaxBendingPhi();
    if(fabs(BendingPhiMax) < BendingPhiLowerTH[Algorithm-1]) {
        BendingWise = 0;
        Charge = 0;
        isGoodPattern = 0;
    }
    else {
        BendingWise = (BendingPhiMax > 0.) ? 1 : -1;
        Charge = BendingWise * (int)(fabs(MeanMagneticField.z()) / MeanMagneticField.z()) * -1;
    }

    if(fabs(BendingPhiMax) <= BendingPhiUpperTH[Algorithm-1] && fabs(BendingPhiMax) >= BendingPhiLowerTH[Algorithm-1]) {
        double theBendingPhi = fabs(BendingPhiMax);
        if(fabs(BendingPhiMax) > BendingPhiFitValueUpperLimit[Algorithm-1])
            theBendingPhi = BendingPhiFitValueUpperLimit[Algorithm-1];
        MeanPt = theBendingPhi * theBendingPhi * MeanPt_Parameter0[Algorithm-1] + theBendingPhi * MeanPt_Parameter1[Algorithm-1] + MeanPt_Parameter2[Algorithm-1];
        if(fabs(BendingPhiMax) > BendingPhiFitSigmaUpperLimit[Algorithm-1])
            theBendingPhi = BendingPhiFitSigmaUpperLimit[Algorithm-1];
        SigmaPt = theBendingPhi * theBendingPhi * SigmaPt_Parameter0[Algorithm-1] + theBendingPhi * SigmaPt_Parameter1[Algorithm-1] + SigmaPt_Parameter2[Algorithm-1];
        SigmaPt *= 3.;
    }
    else {
        isGoodPattern = 0;
        MeanPt = 0.;
        SigmaPt = 0.;
    }

    GlobalVector RefSegmentVector = GlobalVector((RefSegment.second)->globalPosition() - (RefSegment.first)->globalPosition());
    GlobalVector InitialPt = GlobalVector(RefSegmentVector.x(), RefSegmentVector.y(), 0);
    GlobalVector InitialMomentum = InitialPt.unit() * DistanceXY + GlobalVector(0, 0, DistanceZ);
    /*
    if(isVertexConstraint == false)
        InitialMomentum = InitialPt.unit() * DistanceXY + GlobalVector(0, 0, DistanceZ);
    else
        InitialMomentum = GlobalVector(InitialPt.x()*theRefRecHit->globalPosition().perp()/(cos(BendingPhiCollection[BendingPhiCollection.size()-1])*InitialPt.perp()), InitialPt.y()*theRefRecHit->globalPosition().perp()/(cos(BendingPhiCollection[BendingPhiCollection.size()-1])*InitialPt.perp()), theRefRecHit->globalPosition().z());
    */
    Momentum = InitialMomentum.unit() * InitialMomentum.mag() * MeanPt / InitialMomentum.perp();

    if(debug) cout << "BendingPhiMax: " << BendingPhiMax << ", MeanPt: " << MeanPt << "SigmaPt: " << SigmaPt << ", Momentum: " << Momentum << endl;
}

double RPCSeedPattern::findMaxBendingPhi() {
    double MaxBendingPhi = 0.;
    double TempBendingPhi;
    for(unsigned int Index = 0; Index < BendingPhiCollection.size(); Index++) {
        TempBendingPhi = BendingPhiCollection[Index];
        if(fabs(TempBendingPhi) > fabs(MaxBendingPhi))
            MaxBendingPhi = TempBendingPhi;
    }
    return MaxBendingPhi;
}

RPCSeedPattern::WeightedTrajectorySeed RPCSeedPattern::createFakeSeed(int& isGoodSeed) {
    // Create a fake seed and return
    if(debug) cout << "Now create a fake seed" << endl;

    LocalPoint RefPosition = theRefRecHit->localPosition();
    LocalVector RefMomentum = theRefRecHit->det()->toLocal(Momentum);
    LocalTrajectoryParameters theLTP(RefPosition, RefMomentum, Charge);

    AlgebraicSymMatrix theErrorMatrix(5,0);
    theErrorMatrix = theRefRecHit->parametersError().similarityT(theRefRecHit->projectionMatrix());
    theErrorMatrix[0][0] = 0.;
    LocalTrajectoryError theLTE(theErrorMatrix);

    TrajectoryStateOnSurface theTSOS(theLTP, theLTE, theRefRecHit->det()->surface(), &*theMagneticField);

    DetId theDetId = theRefRecHit->geographicalId();
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
        if(debug) cout <<"Pattern is not good: " << isPatternChecked << ", " << isGoodPattern << ". Create a fake seed instead!" << endl;
        return createFakeSeed(isGoodSeed);
    }

    //const ConstMuonRecHitPointer theRefRecHit = BestRefRecHit();

    LocalPoint RefPosition = theRefRecHit->localPosition();
    LocalVector RefMomentum = theRefRecHit->det()->toLocal(Momentum);
    LocalTrajectoryParameters theLTP(RefPosition, RefMomentum, Charge);

    LocalTrajectoryError theLTE = getErrorMatrix();

    TrajectoryStateOnSurface theTSOS(theLTP, theLTE, theRefRecHit->det()->surface(), &*theMagneticField);
    DetId theDetId = theRefRecHit->geographicalId();
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

LocalTrajectoryError RPCSeedPattern::getErrorMatrix() {

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


    LocalError theLocalError = theRefRecHit->localPositionError();
    double dXatRef = sqrt(theLocalError.xx());
    double dYatRef = sqrt(theLocalError.yy());
    
    const GeomDetUnit* RefRPCRoll = theRefRecHit->detUnit();
    LocalPoint RecHit0 = RefRPCRoll->toLocal((RefSegment.first)->globalPosition());
    LocalPoint RecHit1 = RefRPCRoll->toLocal((RefSegment.second)->globalPosition());
    LocalPoint theLocalVertex = RefRPCRoll->toLocal(GlobalPoint(0,0,0));
    LocalVector LocalSegment0 = (LocalVector)(RecHit1 - RecHit0);
    LocalVector LocalSegment1 = LocalVector(LocalSegment0.x()+2*dXatRef, LocalSegment0.y(), LocalSegment0.z());
    LocalVector LocalSegment2 = LocalVector(LocalSegment0.x()-2*dXatRef, LocalSegment0.y(), LocalSegment0.z());
    LocalVector LocalSegment3 = LocalVector(LocalSegment0.x(), LocalSegment0.y()+2*dYatRef, LocalSegment0.z());
    LocalVector LocalSegment4 = LocalVector(LocalSegment0.x(), LocalSegment0.y()-2*dYatRef, LocalSegment0.z());
    double dXdZ1 = LocalSegment1.x() / LocalSegment1.z() - LocalSegment0.x() / LocalSegment0.z();
    double dXdZ2 = LocalSegment2.x() / LocalSegment2.z() - LocalSegment0.x() / LocalSegment0.z();
    dXdZ = (fabs(dXdZ1) >= fabs(dXdZ2)) ? dXdZ1 : dXdZ2;
    // the RPCrecHit has large bias along local Y direction(along the Strip), VertexConstraint measurement is no enough for evaluating the dYdZ error, or the kalman filter may miss some accurate DT/CSC recHits.
    
    if(isVertexConstraint == false) {
        double dYdZ1 = LocalSegment3.y() / LocalSegment3.z() - LocalSegment0.y() / LocalSegment0.z();
        double dYdZ2 = LocalSegment4.y() / LocalSegment4.z() - LocalSegment0.y() / LocalSegment0.z();
        dYdZ = (fabs(dYdZ1) >= fabs(dYdZ2)) ? dYdZ1 : dYdZ2;
    }
    else {
        dYdZ = dYatRef / (RecHit0.z()-theLocalVertex.z());
    }
    
    /*
    double dYdZ1 = LocalSegment3.y() / LocalSegment3.z() - LocalSegment0.y() / LocalSegment0.z();
    double dYdZ2 = LocalSegment4.y() / LocalSegment4.z() - LocalSegment0.y() / LocalSegment0.z();
    dYdZ = (fabs(dYdZ1) >= fabs(dYdZ2)) ? dYdZ1 : dYdZ2;
    */
    double SigmaP = SigmaPt * Momentum.mag() / MeanPt;
    dPInv = SigmaP / (Momentum.mag() * Momentum.mag());

    LocalTrajectoryError theLTE(dX, dY, dXdZ, dYdZ, dPInv);
    return theLTE;
}
