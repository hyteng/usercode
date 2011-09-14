#ifndef MyModule_RPCSeedProducer_RPCSeedPattern_H
#define MyModule_RPCSeedProducer_RPCSeedPattern_H

/**  \class RPCSeedPattern
 *
 *  \author Haiyun.Teng - Peking University
 *
 *
 */


#include "FWCore/Framework/interface/EventSetup.h"
#include <DataFormats/TrajectorySeed/interface/TrajectorySeed.h>
#include <RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h>
#include <MagneticField/Engine/interface/MagneticField.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <TrackingTools/TrajectoryParametrization/interface/LocalTrajectoryError.h>
#include <vector>

#ifndef upper_limit_pt
#define upper_limit_pt 100
#endif

#ifndef lower_limit_pt
#define lower_limit_pt 3.0
#endif


class RPCSeedPattern {

    typedef MuonTransientTrackingRecHit::MuonRecHitPointer MuonRecHitPointer;
    typedef MuonTransientTrackingRecHit::ConstMuonRecHitPointer ConstMuonRecHitPointer;
    typedef MuonTransientTrackingRecHit::MuonRecHitContainer MuonRecHitContainer;
    typedef MuonTransientTrackingRecHit::ConstMuonRecHitContainer ConstMuonRecHitContainer;

    public:
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
    void DoubleSegmentPattern();
    void computePtwithDoubleSegment();
    void checkDoubleSegmentPattern();
    void SingleSegmentPattern();
    void computePtwithSingleSegment();
    void checkSingleSegmentPattern();
    void measureRecHitandMagneticField();
    int findIntervalIndex();
    GlobalVector getMeanMagneticField(const int IntervalIndex);
    ConstMuonRecHitPointer BestRefRecHit() const;
    LocalTrajectoryError getErrorMatrix(const ConstMuonRecHitPointer& theRefRecHit);
    WeightedTrajectorySeed createFakeSeed(int& isGoodSeed);
    WeightedTrajectorySeed createSeed(int& isGoodSeed);
    //double getDistance(const ConstMuonRecHitPointer& precHit, const GlobalVector& Center) const;

    // ----------member data ---------------------------

    // parameters for configuration
    double ZError;
    double MagnecticFieldThreshold;
    unsigned int sampleCount;
    int AlgorithmType;
    std::vector<double> BendingPhiLowerTH;
    std::vector<double> BendingPhiUpperTH;
    std::vector<double> MeanPt_Parameter0;
    std::vector<double> MeanPt_Parameter1;
    std::vector<double> MeanPt_Parameter2;
    std::vector<double> SigmaPt_Parameter0;
    std::vector<double> SigmaPt_Parameter1;
    std::vector<double> SigmaPt_Parameter2;
    
    // signals
    bool isConfigured;
    bool isRecHitset;
    // recHits collection
    ConstMuonRecHitContainer theRecHits;
    std::vector<unsigned int> theRecHitLayers;
    edm::ESHandle<MagneticField> theMagneticField;
    // magnetic field info
    std::vector<GlobalVector> sampleMagneticField;
    std::vector<bool> IntervalMagneticFlux;
    GlobalVector MeanMagneticField;
    // double segment pattern
    int Algorithm;
    int RefIndex;
    std::vector<RPCSegment> SegmentRB;
    int BendingWise;
    int ZDirection;
    double MeanPt;
    double SigmaPt;
    int Charge;
    double DistanceXY;
    double DistanceZ;
    GlobalVector Momentum;
    double PhiResidualRB[4];
    // pattern estimation part
    bool isPatternChecked;
    int isGoodPattern;
    double PatternQuality;
};

#endif
