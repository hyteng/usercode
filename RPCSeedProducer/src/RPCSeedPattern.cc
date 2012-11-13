/*
 *  See header file for a description of this class.
 *
 *
 *  $Date: 2012/11/03 03:11:21 $
 *  $Revision: 1.21 $
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
    SampleCount = iConfig.getParameter<unsigned int>("SampleCount");
    AlgorithmType = iConfig.getParameter<unsigned int>("AlgorithmType");
    isVertexConstraint = iConfig.getParameter<bool>("isVertexConstraint");
    isContinuousFilter = iConfig.getParameter<bool>("isContinuousFilter");
    Cut0 = iConfig.getParameter< vector<double> >("Cut0");
    Cut1 = iConfig.getParameter< vector<double> >("Cut1");
    Cut2 = iConfig.getParameter< vector<double> >("Cut2");
    CutMax = iConfig.getParameter< vector<double> >("CutMax");
    BendingPhiLowerTH = iConfig.getParameter< vector<double> >("BendingPhiLowerTH");
    BendingPhiUpperTH = iConfig.getParameter< vector<double> >("BendingPhiUpperTH");
    ProbingPhiUpperTH = iConfig.getParameter< vector<double> >("ProbingPhiUpperTH");
    ProbingPhiLowerTH = iConfig.getParameter< vector<double> >("ProbingPhiLowerTH");
    ExhaustivePhiTH = iConfig.getParameter< vector<double> >("ExhaustivePhiTH");
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

    if(Cut0.size() == 0 || Cut1.size() == 0 || Cut2.size() == 0 || CutMax.size() == 0 || BendingPhiFitValueUpperLimit.size() == 0 || BendingPhiFitSigmaUpperLimit.size() == 0 || BendingPhiLowerTH.size() == 0 || BendingPhiUpperTH.size() == 0 || ProbingPhiUpperTH.size() == 0 || ProbingPhiLowerTH.size() == 0 || ExhaustivePhiTH.size() == 0 || MeanPt_Parameter0.size() == 0 || MeanPt_Parameter1.size() == 0 || MeanPt_Parameter2.size() == 0 || SigmaPt_Parameter0.size() == 0 || SigmaPt_Parameter1.size() == 0 || SigmaPt_Parameter2.size() == 0 )
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
        theRecHitPosition[Layer] = (*Iter)->globalPosition();
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
    createRPCPattern();
    if(Algorithm >= 1)
        checkRPCPattern();

    computeRPCPattern();

    return createSeed(isGoodSeed);
}


void RPCSeedPattern::measureRecHitandMagneticField() {
    // Get distance and magnetice field sampling information, recHit's position is not the border of Chamber and Iron
    SampleMagneticField.clear();
    DistanceXY = 0;
    IntervalMagneticFlux.clear();
    for(ConstMuonRecHitContainer::const_iterator RecHitIter = theRecHits.begin(); RecHitIter != (theRecHits.end()-1); RecHitIter++) {
        GlobalPoint FirstPosition = (*RecHitIter)->globalPosition();
        GlobalPoint LastPosition = (*(RecHitIter+1))->globalPosition();
        DistanceXY += ((GlobalVector)(LastPosition - FirstPosition)).perp();
        GlobalPoint* SamplePosition = new GlobalPoint[SampleCount];
        double dX = (LastPosition.x() - FirstPosition.x()) / (SampleCount + 1);
        double dY = (LastPosition.y() - FirstPosition.y()) / (SampleCount + 1);
        double dZ = (LastPosition.z() - FirstPosition.z()) / (SampleCount + 1);
        bool isLowMagneticField = true;
        for(unsigned int index = 0; index < SampleCount; index++) {
            SamplePosition[index] = GlobalPoint((FirstPosition.x()+dX*(index+1)), (FirstPosition.y()+dY*(index+1)), (FirstPosition.z()+dZ*(index+1)));
            GlobalVector TempMagneticField = theMagneticField->inTesla(SamplePosition[index]);
            // for endcap the magnetic field could be complex so use mag() intead of z() while checking
            if(fabs(TempMagneticField.z()) > MagnecticFieldThreshold)
                isLowMagneticField = false;
            if(debug) cout << "Sampling magnetic field : " << TempMagneticField << endl;
            SampleMagneticField.push_back(TempMagneticField);
        }
        if(debug) cout << "IntervalMagneticFlux(0-high,1-low): " << isLowMagneticField << endl;
        IntervalMagneticFlux.push_back(isLowMagneticField);
        delete [] SamplePosition;
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
        MeanMagneticField = getMeanMagneticField(0);
    }

    if(isVertexConstraint == true)
        DistanceZ = theRecHits[theRecHits.size()-1]->globalPosition().z() * (theRecHits[theRecHits.size()-1]->globalPosition().perp()-theRecHits[0]->globalPosition().perp()) / theRecHits[theRecHits.size()-1]->globalPosition().perp();
    else
        DistanceZ = theRecHits[theRecHits.size()-1]->globalPosition().z() - theRecHits[0]->globalPosition().z();

    ZDirection = 0;
    if(fabs(DistanceZ) > ZError || isVertexConstraint == true) {
        if(DistanceZ > ZError)
            ZDirection = 1;
        else
            ZDirection = -1;
    }

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
        IntervalIndex = 0;

    if(debug) cout << "Find IntervalIndex: " << IntervalIndex << endl;

    return IntervalIndex;
}

GlobalVector RPCSeedPattern::getMeanMagneticField(const int IntervalIndex) {
    GlobalVector theMagneticField(0, 0, 0);
    unsigned int SampleNumber = 0;
    if(IntervalIndex < 0 || SampleCount <= 0)
        return theMagneticField;
    for(unsigned int i = IntervalIndex*SampleCount; i < (IntervalIndex+1)*SampleCount; i++) {
        GlobalVector TempMagnaeticField = SampleMagneticField[i];
        if(TempMagnaeticField.mag() > MagnecticFieldThreshold) {
            theMagneticField += TempMagnaeticField;
            SampleNumber++;
        }
    }
    if(SampleNumber == 0)
        SampleNumber = 1.;
    return theMagneticField / SampleNumber;
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
    if(isBarrel == true && isNegativeEndcap == false && isPositiveEndcap == false && theBarrelOccupancyCode == BarrelPatternCode1)
        AlgorithmChoice.push_back(1);
    if(isBarrel == true && isNegativeEndcap == false && isPositiveEndcap == false && theBarrelOccupancyCode == BarrelPatternCode2)
        AlgorithmChoice.push_back(2);
    if(isBarrel == true && isNegativeEndcap == false && isPositiveEndcap == false && theBarrelOccupancyCode == BarrelPatternCode3)
        AlgorithmChoice.push_back(3);
    if(isBarrel == true && isNegativeEndcap == false && isPositiveEndcap == false && theBarrelOccupancyCode == BarrelPatternCode4)
        AlgorithmChoice.push_back(4);
    if(isBarrel == true && isNegativeEndcap == false && isPositiveEndcap == false && theBarrelOccupancyCode == BarrelPatternCode5)
        AlgorithmChoice.push_back(5);
    if(isBarrel == true && isNegativeEndcap == false && isPositiveEndcap == false && theBarrelOccupancyCode == BarrelPatternCode6)
        AlgorithmChoice.push_back(6);
    if(isBarrel == true && isNegativeEndcap == false && isPositiveEndcap == false && theBarrelOccupancyCode == BarrelPatternCode7)
        AlgorithmChoice.push_back(7);
    if(isBarrel == true && isNegativeEndcap == false && isPositiveEndcap == false && theBarrelOccupancyCode == BarrelPatternCode8)
        AlgorithmChoice.push_back(8);
    if(isBarrel == true && isNegativeEndcap == false && isPositiveEndcap == false && theBarrelOccupancyCode == BarrelPatternCode9)
        AlgorithmChoice.push_back(9);
    if(isBarrel == true && isNegativeEndcap == false && isPositiveEndcap == false && theBarrelOccupancyCode == BarrelPatternCode10)
        AlgorithmChoice.push_back(10);
    if(isBarrel == true && isNegativeEndcap == false && isPositiveEndcap == false && theBarrelOccupancyCode == BarrelPatternCode11)
        AlgorithmChoice.push_back(11);
    if(isBarrel == true && isNegativeEndcap == false && isPositiveEndcap == false && theBarrelOccupancyCode == BarrelPatternCode12)
        AlgorithmChoice.push_back(12);
    if(isBarrel == true && isNegativeEndcap == false && isPositiveEndcap == false && theBarrelOccupancyCode == BarrelPatternCode13)
        AlgorithmChoice.push_back(13);
    if(isBarrel == true && isNegativeEndcap == false && isPositiveEndcap == false && theBarrelOccupancyCode == BarrelPatternCode14)
        AlgorithmChoice.push_back(14);
    if(isBarrel == true && isNegativeEndcap == false && isPositiveEndcap == false && theBarrelOccupancyCode == BarrelPatternCode15)
        AlgorithmChoice.push_back(15);
    if(isBarrel == true && isNegativeEndcap == false && isPositiveEndcap == false && theBarrelOccupancyCode == BarrelPatternCode16)
        AlgorithmChoice.push_back(16);
    if(isBarrel == true && isNegativeEndcap == false && isPositiveEndcap == false && theBarrelOccupancyCode == BarrelPatternCode17)
        AlgorithmChoice.push_back(17);
    if(isBarrel == true && isNegativeEndcap == false && isPositiveEndcap == false && theBarrelOccupancyCode == BarrelPatternCode18)
        AlgorithmChoice.push_back(18);
    if(isBarrel == true && isNegativeEndcap == false && isPositiveEndcap == false && theBarrelOccupancyCode == BarrelPatternCode19)
        AlgorithmChoice.push_back(19);

    // Auto choice or manual choise
    int FinalAlgorithm = -1;
    if(AlgorithmType == 0 && AlgorithmChoice.size() > 0) {
        if(checkParameters((unsigned int)(AlgorithmChoice[0])))
            FinalAlgorithm = AlgorithmChoice[0];
    }
    else {
        for(unsigned int i = 0; i < AlgorithmChoice.size(); i++)
            if(AlgorithmType == AlgorithmChoice[i])
                if(checkParameters((unsigned int)(AlgorithmChoice[i])))
                    FinalAlgorithm = AlgorithmChoice[i];
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
    if(Cut0.size() < AlgorithmIndex || Cut1.size() < AlgorithmIndex || Cut2.size() < AlgorithmIndex || CutMax.size() < AlgorithmIndex || BendingPhiFitValueUpperLimit.size() < AlgorithmIndex || BendingPhiFitSigmaUpperLimit.size() < AlgorithmIndex || BendingPhiLowerTH.size() < AlgorithmIndex || BendingPhiUpperTH.size() < AlgorithmIndex || ProbingPhiUpperTH.size() < AlgorithmIndex || ProbingPhiLowerTH.size() < AlgorithmIndex || ExhaustivePhiTH.size() < AlgorithmIndex || MeanPt_Parameter0.size() < AlgorithmIndex || MeanPt_Parameter1.size() < AlgorithmIndex || MeanPt_Parameter2.size() < AlgorithmIndex || SigmaPt_Parameter0.size() < AlgorithmIndex || SigmaPt_Parameter1.size() < AlgorithmIndex || SigmaPt_Parameter2.size() < AlgorithmIndex )
        isParametersSet = false;

    return isParametersSet;
}

void RPCSeedPattern::createRPCPattern() {
    
    if(Algorithm <= 0) {
        isGoodPattern = -1;
        isPatternChecked = true;
        return;
    }

    BendingFilter.clear();

    if(isVertexConstraint == true)
        for(unsigned int i = 0; i < theRecHitLayers.size()-1; i++)
            for(unsigned int j = i; j < theRecHitLayers.size()-1; j++)
                for(unsigned int k = j; k < theRecHitLayers.size()-1; k++)
                    for(unsigned int l = k+1; l < theRecHitLayers.size(); l++) {
                        // we are in vertex constraint mode, so take vertex as much as possible
                        if(i != j)
                            continue;
                        // only RB1/2 could link with vertex, since RB3 are close to RB4 and bring in large widthing
                        if(i == j && ((theRecHitLayers[i] < BarrelLayerNumber && theRecHitLayers[i] >= (BarrelLayerNumber-2)) || (theRecHitLayers[i] >= BarrelLayerNumber && theRecHitLayers[i] < (BarrelLayerNumber+EachEndcapLayerNumber) && theRecHitLayers[i] >= (BarrelLayerNumber+EachEndcapLayerNumber-2)) || (theRecHitLayers[i] >= (BarrelLayerNumber+EachEndcapLayerNumber) && theRecHitLayers[i] < (BarrelLayerNumber+2*EachEndcapLayerNumber) && theRecHitLayers[i] >= (BarrelLayerNumber+2*EachEndcapLayerNumber-2))))
                            continue;
                        // close layers bring in large widthing
                        if((theRecHitLayers[i] < (BarrelLayerNumber-1) && theRecHitLayers[j] == (theRecHitLayers[i]+1)) || (theRecHitLayers[k] < (BarrelLayerNumber-1) && theRecHitLayers[l] == (theRecHitLayers[k]+1)))
                            continue;
                        BendingFilter.push_back(BendingPhiIndexType(theRecHitLayers[i], theRecHitLayers[j], theRecHitLayers[k], theRecHitLayers[l]));
                    }

    if(isVertexConstraint == false)
        for(unsigned int i = 0; i < theRecHitLayers.size()-2; i++)
            for(unsigned int j = i+1; j < theRecHitLayers.size()-1; j++)
                for(unsigned int k = j; k < theRecHitLayers.size()-1; k++)
                    for(unsigned int l = k+1; l < theRecHitLayers.size(); l++) {
                        // close layers bring in large widthing
                        //if((theRecHitLayers[i] < (BarrelLayerNumber-1) && theRecHitLayers[j] == (theRecHitLayers[i]+1)) || (theRecHitLayers[k] < (BarrelLayerNumber-1) && theRecHitLayers[l] == (theRecHitLayers[k]+1)))
                            //continue;
                        BendingFilter.push_back(BendingPhiIndexType(theRecHitLayers[i], theRecHitLayers[j], theRecHitLayers[k], theRecHitLayers[l]));
                    }
                    
    BendingPhiCollection.clear();
    for(unsigned int Index = 0; Index < BendingFilter.size(); Index++) {
        int i = BendingFilter[Index].m[0];
        int j = BendingFilter[Index].m[1];
        int k = BendingFilter[Index].m[2];
        int l = BendingFilter[Index].m[3];
        /*
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
        */
        double TempBendingPhi = getdPhi(i,j,k,l);
        BendingPhiCollection.push_back(TempBendingPhi);
    }
    BendingPhiMax = findMaxBendingPhi();

    // set the signal
    isGoodPattern = -1;
    isPatternChecked = false;
}

double RPCSeedPattern::getdPhi(int i, int j, int k, int l) {

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
    
    return TempBendingPhi;
}

void RPCSeedPattern::checkRPCPattern() {

    if(isPatternChecked == true)
        return;

    if(debug) cout << "checkDoubleSegmentPattern." << endl;

    isGoodPattern = 1;

    if(isVertexConstraint == true) {
        // choose different and far way RPCLayer for small widthing and less reverse bending in validating bendingPhi, Val0 and Val1 are for exhaustive bending correction, wihle Val2 is for patterns not redundent enough to filter a possible wrong bending
        double BendingPhiVal0 = getdPhi(theRecHitLayers[0],theRecHitLayers[0],theRecHitLayers[0],theRecHitLayers[theRecHitLayers.size()-2]);
        double BendingPhiVal1 = getdPhi(theRecHitLayers[1],theRecHitLayers[1],theRecHitLayers[1],theRecHitLayers[theRecHitLayers.size()-1]);
        double BendingPhiVal2 = getdPhi(theRecHitLayers[0],theRecHitLayers[0],theRecHitLayers[0],theRecHitLayers[theRecHitLayers.size()-1]);

        // filter small reverse bending if needed
        if(fabs(BendingPhiMax) < CutMax[Algorithm-1] || fabs(BendingPhiVal0) < Cut0[Algorithm-1] || fabs(BendingPhiVal1) < Cut1[Algorithm-1] || fabs(BendingPhiVal2) < Cut2[Algorithm-1]) {
            if(debug) cout << "block by BendingPhiTH." << endl;
            isGoodPattern = -1;
        }

        if(applyFilter == true)
            if((BendingPhiMax * BendingPhiVal0 < 0. && Cut0[Algorithm-1] > 0.) || (BendingPhiMax * BendingPhiVal1 < 0. && Cut1[Algorithm-1] > 0.) || (BendingPhiMax * BendingPhiVal2 < 0. && Cut2[Algorithm-1] > 0.)) {
                // only below the ExhaustivePhiTH we correct the BendingPhiMax, or it should come from a noise bias
                if(fabs(BendingPhiMax) <= ExhaustivePhiTH[Algorithm-1] && Cut0[Algorithm-1] > 0. && Cut1[Algorithm-1] > 0. && BendingPhiMax * BendingPhiVal0 < 0. && BendingPhiMax * BendingPhiVal1 < 0.)
                    BendingPhiMax *= -1.;
                else {
                    if(debug) cout << "block by Filter." << endl;
                    isGoodPattern = -1;
                }
            }
    }
    else {
        // choose different and far way RPCLayer for small widthing and less reverse bending in validating bendingPhi
        double BendingPhiVal0 = getdPhi(theRecHitLayers[0],theRecHitLayers[1],theRecHitLayers[1],theRecHitLayers[theRecHitLayers.size()-2]);
        double BendingPhiVal1 = getdPhi(theRecHitLayers[1],theRecHitLayers[2],theRecHitLayers[2],theRecHitLayers[theRecHitLayers.size()-1]);
        double BendingPhiVal2 = getdPhi(theRecHitLayers[0],theRecHitLayers[1],theRecHitLayers[1],theRecHitLayers[theRecHitLayers.size()-1]);

        // filter small reverse bending if needed
        if(fabs(BendingPhiMax) < CutMax[Algorithm-1] || fabs(BendingPhiVal0) < Cut0[Algorithm-1] || fabs(BendingPhiVal1) < Cut1[Algorithm-1] || fabs(BendingPhiVal2) < Cut2[Algorithm-1]) {
            if(debug) cout << "block by BendingPhiTH." << endl;
            isGoodPattern = -1;
        }

        if(applyFilter == true)
            //if(fabs(BendingPhiMax) <= ExhaustivePhiTH[Algorithm-1])
                if((BendingPhiMax * BendingPhiVal0 < 0. && Cut0[Algorithm-1] > 0.) || (BendingPhiMax * BendingPhiVal1 < 0. && Cut1[Algorithm-1] > 0.) || (BendingPhiMax * BendingPhiVal2 < 0. && Cut2[Algorithm-1] > 0.)) {
                    if(Cut0[Algorithm-1] > 0. && Cut1[Algorithm-1] > 0. && BendingPhiMax * BendingPhiVal0 < 0. && BendingPhiMax * BendingPhiVal1 < 0.)
                        BendingPhiMax *= -1.;
                    else {
                        if(debug) cout << "block by Filter." << endl;
                        isGoodPattern = -1;
                    }
                }
    }

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

    isPatternChecked = true;
}

void RPCSeedPattern::computeRPCPattern() {
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
    
    BendingWise = (BendingPhiMax > 0.) ? 1 : -1;
    Charge = BendingWise * (int)(fabs(MeanMagneticField.z()) / MeanMagneticField.z()) * -1;

    if(fabs(BendingPhiMax) <= BendingPhiUpperTH[Algorithm-1] && fabs(BendingPhiMax) >= BendingPhiLowerTH[Algorithm-1]) {
        double theBendingPhi = fabs(BendingPhiMax);
        if(fabs(BendingPhiMax) > BendingPhiFitValueUpperLimit[Algorithm-1])
            theBendingPhi = BendingPhiFitValueUpperLimit[Algorithm-1];
        MeanPt = theBendingPhi * theBendingPhi * MeanPt_Parameter0[Algorithm-1] + theBendingPhi * MeanPt_Parameter1[Algorithm-1] + MeanPt_Parameter2[Algorithm-1];
        if(fabs(BendingPhiMax) > BendingPhiFitSigmaUpperLimit[Algorithm-1])
            theBendingPhi = BendingPhiFitSigmaUpperLimit[Algorithm-1];
        SigmaPt = theBendingPhi * theBendingPhi * SigmaPt_Parameter0[Algorithm-1] + theBendingPhi * SigmaPt_Parameter1[Algorithm-1] + SigmaPt_Parameter2[Algorithm-1];
        //SigmaPt *= 3.; //in STA filter it will be applied
    }

    if(fabs(BendingPhiMax) <= ProbingPhiLowerTH[Algorithm-1]) {
        isGoodPattern = -1;
        MeanPt = UpperLimitMeanPt2;
        SigmaPt = UpperLimitSigmaPt2;
        Charge = 0;
    }
    if(fabs(BendingPhiMax) > ProbingPhiLowerTH[Algorithm-1] && fabs(BendingPhiMax) < ProbingPhiUpperTH[Algorithm-1]){
        MeanPt = UpperLimitMeanPt2;
        SigmaPt = UpperLimitSigmaPt2;
    }
    if(fabs(BendingPhiMax) >= ProbingPhiUpperTH[Algorithm-1] && fabs(BendingPhiMax) < BendingPhiLowerTH[Algorithm-1]) {
        double CenterBendingPhi = BendingPhiLowerTH[Algorithm-1];
        double CenterMeanPt = CenterBendingPhi * CenterBendingPhi * MeanPt_Parameter0[Algorithm-1] + CenterBendingPhi * MeanPt_Parameter1[Algorithm-1] + MeanPt_Parameter2[Algorithm-1];
        double CenterSigmaPt = CenterBendingPhi * CenterBendingPhi * SigmaPt_Parameter0[Algorithm-1] + CenterBendingPhi * SigmaPt_Parameter1[Algorithm-1] + SigmaPt_Parameter2[Algorithm-1];
        MeanPt = CenterMeanPt - (CenterBendingPhi - fabs(BendingPhiMax)) * (CenterMeanPt - UpperLimitMeanPt1) / (BendingPhiLowerTH[Algorithm-1] - ProbingPhiUpperTH[Algorithm-1]);
        SigmaPt = CenterSigmaPt - (CenterBendingPhi - fabs(BendingPhiMax)) * (CenterSigmaPt - UpperLimitSigmaPt1) / (BendingPhiLowerTH[Algorithm-1] - ProbingPhiUpperTH[Algorithm-1]);
    }
    if(fabs(BendingPhiMax) > BendingPhiUpperTH[Algorithm-1]) {
        MeanPt = LowerLimitMeanPt;
        SigmaPt = LowerLimitSigmaPt;
    }

    //GlobalVector RefSegmentVector = GlobalVector((RefSegment.second)->globalPosition() - (RefSegment.first)->globalPosition());
    //GlobalVector InitialPt = GlobalVector(RefSegmentVector.x(), RefSegmentVector.y(), 0);
    GlobalVector InitialPt = correctPhiatRef();
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

GlobalVector RPCSeedPattern::correctPhiatRef() {
    GlobalVector RefSegmentVector = GlobalVector((RefSegment.second)->globalPosition() - (RefSegment.first)->globalPosition());
    double LinkedPhi = RefSegmentVector.unit().phi().value();
    double LinkedTheta = RefSegmentVector.unit().theta();
    double BendingPhiCorrection = 0.;
    double SampleHelixS = RefSegmentVector.perp() / SampleCount;
    for(unsigned int i = RefIndex*SampleCount; i < (RefIndex+1)*SampleCount; i++) {
        GlobalVector TempMagnaeticField = SampleMagneticField[i];
        double TempBz = TempMagnaeticField.z();
        if(TempBz == 0.|| MeanPt == 0. || Charge == 0)
            continue;
        double TempBendingRadius = fabs((MeanPt*100.)/(0.3*TempBz));
        double TempBendingPhi = ((double)Charge*-1.*fabs(TempBz)*SampleHelixS)/(TempBendingRadius*TempBz);
        BendingPhiCorrection += TempBendingPhi;
    }
    double CorrectedPhi = LinkedPhi - BendingPhiCorrection;
    double CorrectedTheta = PI/2.;
    if(debug) cout << "LinkedPhi: " << LinkedPhi << ", CorrectedPhi: " << CorrectedPhi << endl;
    GlobalVector RefPtVector(GlobalVector::Polar(CorrectedTheta, CorrectedPhi, 1.) );
    return RefPtVector;
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
    double dX = sqrt(theErrorMatrix[3][3]);
    double dY = sqrt(theErrorMatrix[4][4]);
    double dXdZ = sqrt(theErrorMatrix[1][1]);
    double dYdZ = sqrt(theErrorMatrix[2][2]);
    double dPInv = sqrt(theErrorMatrix[0][0]);
    LocalTrajectoryError theLTE(dX, dY, dXdZ, dYdZ, dPInv);

    TrajectoryStateOnSurface theTSOS(theLTP, theLTE, theRefRecHit->det()->surface(), &*theMagneticField);

    DetId theDetId = theRefRecHit->geographicalId();
    //TrajectoryStateTransform theTST;
    PTrajectoryStateOnDet SeedTSOS = trajectoryStateTransform::persistentState(theTSOS, theDetId.rawId());

    edm::OwnVector<TrackingRecHit> container;
    for(ConstMuonRecHitContainer::const_iterator iter=theRecHits.begin(); iter!=theRecHits.end(); iter++)
        container.push_back((*iter)->hit()->clone());

    TrajectorySeed theSeed(SeedTSOS, container, alongMomentum);
    WeightedTrajectorySeed theWeightedSeed;
    theWeightedSeed.first = theSeed;
    theWeightedSeed.second = PatternQuality;
    isGoodSeed = isGoodPattern;

    //delete SeedTSOS;
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
    //TrajectoryStateTransform theTST;
    PTrajectoryStateOnDet SeedTSOS = trajectoryStateTransform::persistentState(theTSOS, theDetId.rawId());

    edm::OwnVector<TrackingRecHit> container;
    for(ConstMuonRecHitContainer::const_iterator iter=theRecHits.begin(); iter!=theRecHits.end(); iter++) {
        // This casting withou clone will cause memory overflow when used in push_back
        // Since container's deconstructor functiion free the pointer menber!
        //TrackingRecHit* pt = dynamic_cast<TrackingRecHit*>(&*(*iter));
        //cout << "Push recHit type " << pt->getType() << endl;
        container.push_back((*iter)->hit()->clone());
    }

    TrajectorySeed theSeed(SeedTSOS, container, alongMomentum);
    WeightedTrajectorySeed theWeightedSeed;
    theWeightedSeed.first = theSeed;
    theWeightedSeed.second = PatternQuality;
    isGoodSeed = isGoodPattern;

    //delete SeedTSOS;
    return theWeightedSeed;
}

LocalTrajectoryError RPCSeedPattern::getErrorMatrix() {

    // for possable lack of error matrix element[i][j], we use 5 parameter to construct theLTE;
    double dX = 0.;
    double dY = 0.;
    double dXdZ = 0.;
    double dYdZ = 0.;
    double dPInv = 0.;

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
    /*
       LocalVector LocalSegment1 = LocalVector(LocalSegment0.x()+2*dXatRef, LocalSegment0.y(), LocalSegment0.z());
       LocalVector LocalSegment2 = LocalVector(LocalSegment0.x()-2*dXatRef, LocalSegment0.y(), LocalSegment0.z());
       LocalVector LocalSegment3 = LocalVector(LocalSegment0.x(), LocalSegment0.y()+2*dYatRef, LocalSegment0.z());
       LocalVector LocalSegment4 = LocalVector(LocalSegment0.x(), LocalSegment0.y()-2*dYatRef, LocalSegment0.z());
       double dXdZ1 = LocalSegment1.x() / LocalSegment1.z() - LocalSegment0.x() / LocalSegment0.z();
       double dXdZ2 = LocalSegment2.x() / LocalSegment2.z() - LocalSegment0.x() / LocalSegment0.z();
       dXdZ = (fabs(dXdZ1) >= fabs(dXdZ2)) ? dXdZ1 : dXdZ2;
       */
    dXdZ = fabs(2*dXatRef/LocalSegment0.z());
    // the RPCrecHit has large bias along local Y direction(along the Strip), VertexConstraint measurement is no enough for evaluating the dYdZ error, or the kalman filter may miss some accurate DT/CSC recHits.

    if(isVertexConstraint == false) {
        //double dYdZ1 = LocalSegment3.y() / LocalSegment3.z() - LocalSegment0.y() / LocalSegment0.z();
        //double dYdZ2 = LocalSegment4.y() / LocalSegment4.z() - LocalSegment0.y() / LocalSegment0.z();
        //dYdZ = (fabs(dYdZ1) >= fabs(dYdZ2)) ? dYdZ1 : dYdZ2;
        dYdZ =fabs(dYatRef/LocalSegment0.z());
    }
    else {
        dYdZ = fabs(dYatRef/(RecHit0.z()-theLocalVertex.z()));
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
