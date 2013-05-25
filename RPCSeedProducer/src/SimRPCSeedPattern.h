#ifndef MyModule_RPCSeedProducer_SimRPCSeedPattern_H
#define MyModule_RPCSeedProducer_SimRPCSeedPattern_H

/**  \class SimRPCSeedPattern
 *
 *  \author Haiyun.Teng - Peking University
 *
 *
 */

#include <DataFormats/TrajectorySeed/interface/TrajectorySeed.h>
#include <RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h>
#include <MagneticField/Engine/interface/MagneticField.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <TrackingTools/TrajectoryParametrization/interface/LocalTrajectoryError.h>
#include "FWCore/Framework/interface/Event.h"
#include <SimDataFormats/TrackingHit/interface/PSimHitContainer.h>
#include <SimDataFormats/RPCDigiSimLink/interface/RPCDigiSimLink.h>
#include <DataFormats/Common/interface/DetSetVector.h>
#include <MyModule/RPCSeedProducer/src/RPCSeedPattern.h>
#include <random>

class SimRPCSeedPattern {

    typedef MuonTransientTrackingRecHit::MuonRecHitPointer MuonRecHitPointer;
    typedef MuonTransientTrackingRecHit::ConstMuonRecHitPointer ConstMuonRecHitPointer;
    typedef MuonTransientTrackingRecHit::MuonRecHitContainer MuonRecHitContainer;
    typedef MuonTransientTrackingRecHit::ConstMuonRecHitContainer ConstMuonRecHitContainer;

    typedef RPCSeedPattern::BendingPhiIndexType BendingPhiIndexType;
    public:
    typedef std::pair<ConstMuonRecHitPointer, ConstMuonRecHitPointer> RPCSegment;
    typedef std::pair<TrajectorySeed, double> WeightedTrajectorySeed;
    
    public:
    SimRPCSeedPattern(); 
    ~SimRPCSeedPattern();
    void configure(const edm::ParameterSet& iConfig);
    void clear();
    void setRecHits(const ConstMuonRecHitContainer& RecHits);
    void setSimData(const edm::Event& iEvent, const edm::EventSetup& iSetup);
    void unsetSimData();

    private:
    friend class RPCSeedFinder;
    WeightedTrajectorySeed seed(const edm::EventSetup& eSetup, int& isGoodSeed); 
    void mapRecHittoSimHit();
    int checkAlgorithm();
    bool checkParameters(unsigned int theAlgorithmType);
    void createRPCPattern();
    double getdPhi(int i, int j, int k, int l);
    double findMaxBendingPhi();
    void checkRPCPattern();
    void computePatternfromSimData();
    void checkPatternfromSimData();
    void measureRecHitandMagneticField();
    int findIntervalIndex();
    GlobalVector getMeanMagneticField(const int IntervalIndex);
    LocalTrajectoryError getErrorMatrix();
    WeightedTrajectorySeed createFakeSeed(int& isGoodSeed);
    WeightedTrajectorySeed createSeed(int& isGoodSeed);
    //double getDistance(const ConstMuonRecHitPointer& precHit, const GlobalVector& Center) const;

    // ----------member data ---------------------------
    edm::InputTag SimHitTag_;
    edm::InputTag RPCDigiSimLinkTag_;

    edm::Handle<edm::PSimHitContainer> pSimHits;
    edm::Handle< edm::DetSetVector<RPCDigiSimLink> > theLinkDigis;
    
    bool isSimDataset;
    double SeedPurityTH;
    //std::map<PSimHit, ConstMuonRecHitPointer> RecHit2SimHitMap;
    int RefTrackId;
    int RefParticleType;
    std::map<ConstMuonRecHitPointer, PSimHit> RecHit2SimHitMap;
    // parameters for configuration
    double SimSigmaPt;
    double SimBiasPt;
    double SimSigmaPhi;
    double SimBiasPhi;
    double SimSigmaEta;
    double SimBiasEta;
    double ZError;
    double MagnecticFieldThreshold;
    unsigned int SampleCount;
    int AlgorithmType;
    bool isVertexConstraint;
    bool isContinuousFilter;
    bool applyFilter;
    std::vector<double> Cut0;
    std::vector<double> Cut1;
    std::vector<double> Cut2;
    std::vector<double> CutMax;
    std::vector<double> BendingPhiLowerTH;
    std::vector<double> BendingPhiUpperTH;
    std::vector<double> ProbingPhiLowerTH;
    std::vector<double> ProbingPhiUpperTH;
    std::vector<double> ExhaustivePhiTH;
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
    GlobalPoint theRecHitPosition[6];
    edm::ESHandle<MagneticField> theMagneticField;
    // magnetic field info
    std::vector<GlobalVector> SampleMagneticField;
    std::vector<bool> IntervalMagneticFlux;
    GlobalVector MeanMagneticField;
    // double segment pattern
    int Algorithm;
    int RefIndex;
    Geom::Phi<float> BendingPhi[4][6];
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
    // pattern estimation part
    int isGoodPattern;
    double PatternQuality;

    std::default_random_engine theRandomGenerator;
};

#endif
