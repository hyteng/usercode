#ifndef MyModule_RPCSeedProducer_RPCSeedPattern_H
#define MyModule_RPCSeedProducer_RPCSeedPattern_H

/**  \class RPCSeedPattern
 *
 *  \author Haiyun.Teng - Peking University
 *
 *
 */


#include "MyModule/RPCSeedProducer/src/RPCSeedData.h"
#include <FWCore/Framework/interface/EventSetup.h>
#include <DataFormats/TrajectorySeed/interface/TrajectorySeed.h>
#include <RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h>
#include <MagneticField/Engine/interface/MagneticField.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <TrackingTools/TrajectoryParametrization/interface/LocalTrajectoryError.h>
#include <vector>

class RPCSeedPattern {

    public:
    typedef MuonTransientTrackingRecHit::MuonRecHitPointer MuonRecHitPointer;
    typedef MuonTransientTrackingRecHit::ConstMuonRecHitPointer ConstMuonRecHitPointer;
    typedef MuonTransientTrackingRecHit::MuonRecHitContainer MuonRecHitContainer;
    typedef MuonTransientTrackingRecHit::ConstMuonRecHitContainer ConstMuonRecHitContainer;

    class BendingPhiIndexType {
        public:
            BendingPhiIndexType() {m[0] = 0; m[1] = 0; m[2] = 0; m[3] = 0;};
            BendingPhiIndexType(int a0, int a1, int a2, int a3) {m[0] = a0; m[1] = a1; m[2] = a2; m[3] = a3;};
            int m[4];
    };

    typedef std::pair<ConstMuonRecHitPointer, ConstMuonRecHitPointer> RPCSegment;
    typedef std::pair<TrajectorySeed, double> WeightedTrajectorySeed;

    public:
    RPCSeedPattern(); 
    ~RPCSeedPattern();
    void configure(const edm::ParameterSet& iConfig);
    void clear();
    void setRecHits(const ConstMuonRecHitContainer& RecHits);

    private:
    friend class RPCSeedFinder;
    WeightedTrajectorySeed seed(const edm::EventSetup& eSetup, int& isGoodSeed); 
    int checkAlgorithm();
    bool checkParameters(unsigned int theAlgorithmType);
    void createPattern();
    void checkDoubleSegmentPattern();
    void checkSingleSegmentPattern();
    void computeSegmentPattern();
    double findMaxBendingPhi();
    void measureRecHitandMagneticField();
    int findIntervalIndex();
    GlobalVector getMeanMagneticField(const int IntervalIndex);
    LocalTrajectoryError getErrorMatrix();
    WeightedTrajectorySeed createFakeSeed(int& isGoodSeed);
    WeightedTrajectorySeed createSeed(int& isGoodSeed);
    //double getDistance(const ConstMuonRecHitPointer& precHit, const GlobalVector& Center) const;

    // ----------member data ---------------------------

    // parameters for configuration
    double ZError;
    double MagnecticFieldThreshold;
    unsigned int sampleCount;
    int AlgorithmType;
    bool isVertexConstraint;
    bool isContinuousFilter;
    bool applyFilter;
    std::vector<double> Cut1234;
    std::vector<double> CutMax;
    std::vector<double> BendingPhiLowerTH;
    std::vector<double> BendingPhiUpperTH;
    std::vector<double> BendingPhiFitValueUpperLimit;
    std::vector<double> BendingPhiFitSigmaUpperLimit;
    std::vector<double> MeanPt_Parameter0;
    std::vector<double> MeanPt_Parameter1;
    std::vector<double> MeanPt_Parameter2;
    std::vector<double> SigmaPt_Parameter0;
    std::vector<double> SigmaPt_Parameter1;
    std::vector<double> SigmaPt_Parameter2;

    // signals
    bool isConfigured;
    bool isRecHitset;
    bool isPatternChecked;
    // recHits collection
    ConstMuonRecHitContainer theRecHits;
    std::vector<int> theRecHitLayers;
    GlobalPoint theRecHitPosition[BarrelLayerNumber+EachEndcapLayerNumber];
    edm::ESHandle<MagneticField> theMagneticField;
    // magnetic field info
    std::vector<GlobalVector> sampleMagneticField;
    std::vector<bool> IntervalMagneticFlux;
    GlobalVector MeanMagneticField;
    // Pattern parameter
    int Algorithm;
    int RefIndex;
    std::vector<double> BendingPhiCollection;
    std::vector<BendingPhiIndexType> BendingFilter;
    RPCSegment RefSegment;
    ConstMuonRecHitPointer theRefRecHit;
    double BendingPhiMax;
    int BendingWise;
    int ZDirection;
    double MeanPt;
    double SigmaPt;
    int Charge;
    double DistanceXY;
    double DistanceZ;
    GlobalVector Momentum;
    // pattern estimation
    int isGoodPattern;
    double PatternQuality;
};

#endif
