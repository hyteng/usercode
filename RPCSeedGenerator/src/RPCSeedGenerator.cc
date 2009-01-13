// -*- C++ -*-
//
// Package:    RPCSeedGenerator
// Class:      RPCSeedGenerator
// 
/**\class RPCSeedGenerator RPCSeedGenerator.cc MyMTCCAnalyzer/RPCSeedGenerator/src/RPCSeedGenerator.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Haiyun Teng
//         Created:  Wed Oct 29 17:24:36 CET 2008
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// special include files
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TFile.h"
#include <vector>

// using other classes
#include "MyMTCCAnalyzer/RPCSeedGenerator/src/RPCSeedHits.h"
#include "MyMTCCAnalyzer/RPCSeedGenerator/src/RPCSeedFinder.h"

// Data Formats 
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"

#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"

// Geometry
#include "Geometry/RPCGeometry/interface/RPCRoll.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "TrackingTools/DetLayers/interface/DetLayer.h"

#include "RecoMuon/MeasurementDet/interface/MuonDetLayerMeasurements.h"
#include "RecoMuon/DetLayers/interface/MuonDetLayerGeometry.h"
#include "RecoMuon/Records/interface/MuonRecoGeometryRecord.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

// Framework
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "math.h"

// C++
#include <vector>

//
// constants, enums and typedefs
//
using namespace std;
using namespace edm;

typedef MuonTransientTrackingRecHit::MuonRecHitPointer MuonRecHitPointer;
typedef MuonTransientTrackingRecHit::ConstMuonRecHitPointer ConstMuonRecHitPointer;
typedef MuonTransientTrackingRecHit::MuonRecHitContainer MuonRecHitContainer;

#ifndef EachEndcapLayerNumber
#define EachEndcapLayerNumber 3
#endif

//
// class decleration
//

class RPCSeedFinder;

class RPCSeedGenerator : public edm::EDProducer {
   public:
      explicit RPCSeedGenerator(const edm::ParameterSet& iConfig);
      ~RPCSeedGenerator();

   private:
      virtual void beginJob(const edm::EventSetup& iSetup);
      virtual void produce(edm::Event& iEvent, const edm::EventSetup& iSetup);
      virtual void endJob();

      void complete(RPCSeedFinder& theSeed, int N, int NumberinMuon, int *SeedinMuon, MuonRecHitContainer (&recHitsMuon)[6], const edm::EventSetup& iSetup);
      void CheckandFill(RPCSeedFinder& theSeed, const edm::EventSetup& iSetup);
      void CheckOverlap(const edm::EventSetup& iSetup);
      void CheckcandidateOverlap(const edm::EventSetup& iSetup);
      bool isShareHit(const edm::OwnVector<TrackingRecHit> RecHits, const TrackingRecHit& hit, edm::ESHandle<RPCGeometry> rpcGeometry);

      // fill seeds from special layers
      void FillSeedfromSpecialLayers(const edm::EventSetup& iSetup, int NumberinBarrel, int *SeedinBarrel, MuonRecHitContainer (&recHitsBarrel)[6], int NumberinEndcap, int *SeedinEndcap, MuonRecHitContainer (&recHitsEndcap)[6], bool isMixBarrelwithEndcap);
      // create special N layers to fill to seeds
      void SpecialLayers(const edm::EventSetup& iSetup, int last, int N, int type, int NumberinBarrel, int *SeedinBarrel, MuonRecHitContainer (&recHitsBarrel)[6], int NumberinEndcap, int *SeedinEndcap, MuonRecHitContainer (&recHitsEndcap)[6]);
      
      // ----------member data ---------------------------
      std::vector<TrajectorySeed> candidateSeeds;
      std::vector<TrajectorySeed> goodSeeds;
      edm::InputTag theRPCRecHits;
      std::vector<int> RangeofLayersinBarrel;
      std::vector<int> RangeofLayersinEndcap;
      bool isSpecialLayers;
      bool isMixBarrelwithEndcap;
      unsigned int BxRange;
      double MaxDeltaPhi;
      double MaxRSD;
      bool isCheckgoodOverlap;
      bool isCheckcandidateOverlap;
      unsigned int ShareRecHitsNumberThreshold;
      std::vector<int> LayersinEndcap;
      std::vector<int> LayersinBarrel;
};

/*
//
// constants, enums and typedefs
//
using namespace std;
using namespace edm;

typedef MuonTransientTrackingRecHit::MuonRecHitPointer MuonRecHitPointer;
typedef MuonTransientTrackingRecHit::ConstMuonRecHitPointer ConstMuonRecHitPointer;
typedef MuonTransientTrackingRecHit::MuonRecHitContainer MuonRecHitContainer;
*/

//
// static data member definitions
//

//
// constructors and destructor
//
RPCSeedGenerator::RPCSeedGenerator(const edm::ParameterSet& iConfig)
{
   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
*/
   // now do what ever other initialization is needed
   produces<TrajectorySeedCollection>("goodSeeds");
   produces<TrajectorySeedCollection>("candidateSeeds");
   
   RangeofLayersinBarrel = iConfig.getParameter< std::vector<int> >("RangeofLayersinBarrel");
   RangeofLayersinEndcap = iConfig.getParameter< std::vector<int> >("RangeofLayersinEndcap");
   isCheckgoodOverlap = iConfig.getParameter<bool>("isCheckgoodOverlap");
   isCheckcandidateOverlap = iConfig.getParameter<bool>("isCheckcandidateOverlap");
   isSpecialLayers = iConfig.getParameter<bool>("isSpecialLayers");
   isMixBarrelwithEndcap = iConfig.getParameter<bool>("isMixBarrelwithEndcap");
   BxRange = iConfig.getParameter<unsigned int>("BxRange");
   MaxDeltaPhi = iConfig.getParameter<double>("MaxDeltaPhi");
   MaxRSD = iConfig.getParameter<double>("MaxRSD");
   LayersinBarrel = iConfig.getParameter< std::vector<int> >("LayersinBarrel");
   LayersinEndcap = iConfig.getParameter< std::vector<int> >("LayersinEndcap");
   theRPCRecHits = iConfig.getParameter<edm::InputTag>("RPCRecHitsLabel");
   ShareRecHitsNumberThreshold = iConfig.getParameter<unsigned int>("ShareRecHitsNumberThreshold");
   cout << endl << "[RPCSeedGenerator] --> Constructor called" << endl;
}


RPCSeedGenerator::~RPCSeedGenerator()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   cout << "[RPCSeedGenerator] --> Destructor called" << endl;
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
RPCSeedGenerator::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
/* This is an event example
   //Read 'ExampleData' from the Event
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);

   //Use the ExampleData to create an ExampleData2 which 
   // is put into the Event
   std::auto_ptr<ExampleData2> pOut(new ExampleData2(*pIn));
   iEvent.put(pOut);
*/

/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/
     // clear goodSeeds from last reconstruction  
  goodSeeds.clear();
  candidateSeeds.clear();
  
  // create the pointer to the Seed container
  auto_ptr<TrajectorySeedCollection> goodCollection(new TrajectorySeedCollection());
  auto_ptr<TrajectorySeedCollection> candidateCollection(new TrajectorySeedCollection());
  
  // Muon Geometry - DT, CSC and RPC 
  edm::ESHandle<MuonDetLayerGeometry> muonLayers;
  iSetup.get<MuonRecoGeometryRecord>().get(muonLayers);

  // get the RPC layers
  vector<DetLayer*> RPCBarrelLayers = muonLayers->barrelRPCLayers();
  const DetLayer* RB4L  = RPCBarrelLayers[5];
  const DetLayer* RB3L  = RPCBarrelLayers[4];
  const DetLayer* RB22L = RPCBarrelLayers[3];
  const DetLayer* RB21L = RPCBarrelLayers[2];
  const DetLayer* RB12L = RPCBarrelLayers[1];
  const DetLayer* RB11L = RPCBarrelLayers[0];
  vector<DetLayer*> RPCEndcapLayers = muonLayers->endcapRPCLayers();
  const DetLayer* REM3L = RPCEndcapLayers[0];
  const DetLayer* REM2L = RPCEndcapLayers[1];
  const DetLayer* REM1L = RPCEndcapLayers[2];
  const DetLayer* REP1L = RPCEndcapLayers[3];
  const DetLayer* REP2L = RPCEndcapLayers[4];
  const DetLayer* REP3L = RPCEndcapLayers[5];

  // get RPC recHits by MuonDetLayerMeasurements, while CSC and DT is set to false and with empty InputTag
  MuonDetLayerMeasurements muonMeasurements(edm::InputTag(), edm::InputTag(), theRPCRecHits, false, false, true);

  // dispatch RPC recHits to the corresponding DetLayer, 6 layers for barrel and 3 layers for each endcap
  MuonRecHitContainer recHitsBarrel[6];
  recHitsBarrel[0] = muonMeasurements.recHits(RB11L, iEvent);
  recHitsBarrel[1] = muonMeasurements.recHits(RB12L, iEvent);
  recHitsBarrel[2] = muonMeasurements.recHits(RB21L, iEvent);
  recHitsBarrel[3] = muonMeasurements.recHits(RB22L, iEvent);
  recHitsBarrel[4] = muonMeasurements.recHits(RB3L, iEvent);
  recHitsBarrel[5] = muonMeasurements.recHits(RB4L, iEvent); 
  MuonRecHitContainer recHitsEndcap[6];
  recHitsEndcap[0] = muonMeasurements.recHits(REM3L, iEvent);
  recHitsEndcap[1] = muonMeasurements.recHits(REM2L, iEvent);
  recHitsEndcap[2] = muonMeasurements.recHits(REM1L, iEvent);
  recHitsEndcap[3] = muonMeasurements.recHits(REP1L, iEvent);
  recHitsEndcap[4] = muonMeasurements.recHits(REP2L, iEvent);
  recHitsEndcap[5] = muonMeasurements.recHits(REP3L, iEvent);

  // print the size of recHits in each DetLayer
  cout << "RB1in " << recHitsBarrel[0].size() << " recHits" << endl;
  cout << "RB1out " << recHitsBarrel[1].size() << " recHits" << endl;
  cout << "RB2in " << recHitsBarrel[2].size() << " recHits" << endl;
  cout << "RB2out " << recHitsBarrel[3].size() << " recHits" << endl;
  cout << "RB3 " << recHitsBarrel[4].size() << " recHits" << endl;
  cout << "RB4 " << recHitsBarrel[5].size() << " recHits" << endl;

  cout << "REM3 " << recHitsEndcap[0].size() << " recHits" << endl;
  cout << "REM2 " << recHitsEndcap[1].size() << " recHits" << endl;
  cout << "REM1 " << recHitsEndcap[2].size() << " recHits" << endl;
  cout << "REP1 " << recHitsEndcap[3].size() << " recHits" << endl;
  cout << "REP2 " << recHitsEndcap[4].size() << " recHits" << endl;
  cout << "REP3 " << recHitsEndcap[5].size() << " recHits" << endl;

  // choose enable layers in Barrel and Endcap
  int SeedinBarrel[6];
  int SeedinEndcap[6];
  int NumberinBarrel = 0;
  int NumberinEndcap = 0;
  
  if(isSpecialLayers == false && isMixBarrelwithEndcap == false)
  {
      for(std::vector<int>::iterator NumberofLayersinBarrel = RangeofLayersinBarrel.begin(); NumberofLayersinBarrel != RangeofLayersinBarrel.end(); NumberofLayersinBarrel++)
      {
          // find N layers out of 6 Barrel Layers to fill to SeedinBarrel
          int N = *NumberofLayersinBarrel;
          int type = 0;
          NumberinBarrel = 0;
          NumberinEndcap = 0;
          SpecialLayers(iSetup, -1, N, type, NumberinBarrel, SeedinBarrel, recHitsBarrel, NumberinEndcap, SeedinEndcap, recHitsEndcap);
      }
      for(std::vector<int>::iterator NumberofLayersinEndcap = RangeofLayersinEndcap.begin(); NumberofLayersinEndcap != RangeofLayersinEndcap.end(); NumberofLayersinEndcap++)
      {
          int N = *NumberofLayersinEndcap;
          int type = 1;
          NumberinBarrel = 0;
          NumberinEndcap = 0;
          // for -Z layers
          SpecialLayers(iSetup, -1, N, type, NumberinBarrel, SeedinBarrel, recHitsBarrel, NumberinEndcap, SeedinEndcap, recHitsEndcap);
          //for +Z layers
          SpecialLayers(iSetup, EachEndcapLayerNumber-1, N, type, NumberinBarrel, SeedinBarrel, recHitsBarrel, NumberinEndcap, SeedinEndcap, recHitsEndcap);
      }
  }
  if(isSpecialLayers == true && isMixBarrelwithEndcap == false)
  {
      // Fill barrel layer for seed
      bool EnoughforBarrel = true;
      int i = 0;
      NumberinBarrel = 0;
      NumberinEndcap = 0;     
      for(std::vector<int>::iterator it = LayersinBarrel.begin(); it != LayersinBarrel.end(); it++, i++)
      {   
          if((*it) != false)
          {
              if(recHitsBarrel[i].size() != 0)
              {
                  SeedinBarrel[NumberinBarrel] = i;
                  NumberinBarrel++;
              }
              else
              {
                  cout << "Not recHits in special Barrel " << i << endl;
                  EnoughforBarrel = false;
              }
          }
      }
      if(!EnoughforBarrel)
          return;
      FillSeedfromSpecialLayers(iSetup, NumberinBarrel, SeedinBarrel, recHitsBarrel, NumberinEndcap, SeedinEndcap, recHitsEndcap, isMixBarrelwithEndcap);

      // Fill -Z and +Z endcap layer
      bool EnoughforEndcap = true;
      // Fill endcap- layer for seed
      i = 0;
      NumberinBarrel = 0;
      NumberinEndcap = 0;
      EnoughforEndcap = true;
      for(std::vector<int>::iterator it = LayersinEndcap.begin(); it != LayersinEndcap.end(); it++, i++)
      {
          if((*it) != false && i < EachEndcapLayerNumber)
          {
              if(recHitsEndcap[i].size() != 0)
              {
                  SeedinEndcap[NumberinEndcap] = i;
                  NumberinEndcap++;
              }
              else
              {
                  cout << "Not recHits in special Endcap " << i << endl;
                  EnoughforEndcap = false;
              }
          }
      }
      if(!EnoughforEndcap)
          return;
      FillSeedfromSpecialLayers(iSetup, NumberinBarrel, SeedinBarrel, recHitsBarrel, NumberinEndcap, SeedinEndcap, recHitsEndcap, isMixBarrelwithEndcap);

      //Fill endcap+ layer for seed
      i = 0;
      NumberinBarrel = 0;
      NumberinEndcap = 0;
      EnoughforEndcap = true;
      for(std::vector<int>::iterator it = LayersinEndcap.begin(); it != LayersinEndcap.end(); it++, i++)
      {
          if((*it) != false && i > EachEndcapLayerNumber)
          {
              if(recHitsEndcap[i].size() != 0)
              {
                  SeedinEndcap[NumberinEndcap] = i;
                  NumberinEndcap++;
              }
              else
              {
                  cout << "Not recHits in special Endcap " << i << endl;
                  EnoughforEndcap = false;
              }
          }
      }
      if(!EnoughforEndcap)
          return;
      FillSeedfromSpecialLayers(iSetup, NumberinBarrel, SeedinBarrel, recHitsBarrel, NumberinEndcap, SeedinEndcap, recHitsEndcap, isMixBarrelwithEndcap);
  }

  // check the overlap if required in configuration
  if(isCheckgoodOverlap == true)
      CheckOverlap(iSetup);
  if(isCheckcandidateOverlap == true)
      CheckcandidateOverlap(iSetup);
      
  if(goodSeeds.size() > 1)
  {
      for(vector<TrajectorySeed>::iterator seed = goodSeeds.begin(); seed != goodSeeds.end(); ++seed)
      {
          int counter =0;
          for(vector<TrajectorySeed>::iterator seed2 = seed; seed2 != goodSeeds.end(); ++seed2) 
              if(seed->startingState().parameters().vector() == seed2->startingState().parameters().vector())
                  ++counter;
          if(counter > 1) 
              goodSeeds.erase(seed--);
      }
  }          

  for(vector<TrajectorySeed>::iterator seed = goodSeeds.begin(); seed != goodSeeds.end(); ++seed)
      goodCollection->push_back(*seed);
  for(vector<TrajectorySeed>::iterator seed = candidateSeeds.begin(); seed != candidateSeeds.end(); ++seed)
      candidateCollection->push_back(*seed);

  if(goodSeeds.size() <= 1)
  {
      cout << "Low efficienc for this  event!" << endl;
  }

  // put the seed to event
  iEvent.put(goodCollection, "goodSeeds");
  iEvent.put(candidateCollection, "candidateSeeds");
}

// ------------ method called once each job just before starting event loop  ------------
void 
RPCSeedGenerator::beginJob(const edm::EventSetup& iSetup)
{
    cout << endl << "Begin jobs" << endl;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
RPCSeedGenerator::endJob() {
    cout << endl << "End jobs" << endl;
}

// create special N layers to fill to seeds
void RPCSeedGenerator::SpecialLayers(const edm::EventSetup& iSetup, int last, int N, int type, int NumberinBarrel, int *SeedinBarrel, MuonRecHitContainer (&recHitsBarrel)[6], int NumberinEndcap,int *SeedinEndcap, MuonRecHitContainer (&recHitsEndcap)[6]) {

    // check type, 0=barrel, 1=endcap, 2=mix
    
    // barrel has 6 layers
    if(type == 0)
    {
        if(N > 6)
            return;
        for(int i = (last + 1); i <= (6 - N + NumberinBarrel); i++)
        {
            if(recHitsBarrel[i].size() != 0)
            {
                SeedinBarrel[NumberinBarrel] = i;
                last = i;
                if(NumberinBarrel < (N -1))
                    SpecialLayers(iSetup, last, N, type, (NumberinBarrel + 1), SeedinBarrel, recHitsBarrel, NumberinEndcap, SeedinEndcap, recHitsEndcap);
                else
                {
                    cout << "find special barrel layers: ";
                    for(int k = 0; k < N; k++)
                        cout << SeedinBarrel[k] <<" ";
                    cout << endl;
                    FillSeedfromSpecialLayers(iSetup, (NumberinBarrel + 1), SeedinBarrel, recHitsBarrel, NumberinEndcap, SeedinEndcap, recHitsEndcap, isMixBarrelwithEndcap);
                }
            }
        }
    }

    // endcap has 3 layers for each -Z and +Z
    if(type == 1)
    {
        if(N > EachEndcapLayerNumber)
            return;
        if(last < 2 || (last == (EachEndcapLayerNumber - 1) && NumberinEndcap != 0))
        {
            for(int i = (last + 1); i <= (EachEndcapLayerNumber - N + NumberinEndcap); i++)
            {
                if(recHitsEndcap[i].size() != 0)
                {
                    SeedinEndcap[NumberinEndcap] = i;
                    last = i;
                    if(NumberinEndcap < (N - 1))
                        SpecialLayers(iSetup, last, N, type, NumberinBarrel, SeedinBarrel, recHitsBarrel, (NumberinEndcap + 1), SeedinEndcap, recHitsEndcap);
                    else
                    {
                        cout << "find special -Z endcap layers: ";
                        for(int k = 0; k < N; k++)
                            cout << SeedinEndcap[k] <<" ";
                        cout << endl;
                        FillSeedfromSpecialLayers(iSetup, NumberinBarrel, SeedinBarrel, recHitsBarrel, (NumberinEndcap + 1), SeedinEndcap, recHitsEndcap, isMixBarrelwithEndcap);
                    }
                }
            }
        }
        else
        {
            for(int i = (last + 1); i <= (EachEndcapLayerNumber * 2 - N + NumberinEndcap); i++)
            {
                if(recHitsEndcap[i].size() != 0)
                {
                    SeedinEndcap[NumberinEndcap] = i;
                    last = i;
                    if(NumberinEndcap < (N - 1))
                        SpecialLayers(iSetup, last, N, type, NumberinBarrel, SeedinBarrel, recHitsBarrel, (NumberinEndcap + 1), SeedinEndcap, recHitsEndcap);
                    else
                    {
                        cout << "find special +Z endcap layers: ";
                        for(int k = 0; k < N; k++)
                            cout << SeedinEndcap[k] <<" ";
                        cout << endl;
                        FillSeedfromSpecialLayers(iSetup, NumberinBarrel, SeedinBarrel, recHitsBarrel, (NumberinEndcap + 1), SeedinEndcap, recHitsEndcap, isMixBarrelwithEndcap);
                    }
                }
            }
        }
    }
}

// fill seeds from special layers
void RPCSeedGenerator::FillSeedfromSpecialLayers(const edm::EventSetup& iSetup, int NumberinBarrel, int *SeedinBarrel, MuonRecHitContainer (&recHitsBarrel)[6], int NumberinEndcap, int *SeedinEndcap, MuonRecHitContainer (&recHitsEndcap)[6], bool isMixBarrelwithEndcap) {

  // for non-Mix case, Barrel needs at least 2/3 layers and Endcap needs at least 2/3 layers w/o vertex
  if(isMixBarrelwithEndcap == false && (NumberinBarrel >= 3 || NumberinEndcap >= 3))
  {
      RPCSeedFinder theSeed;
      int N = 0;
      // Fill the Barrel Seed
      complete(theSeed, N, NumberinBarrel, SeedinBarrel, recHitsBarrel, iSetup);
      // Fill -Z and +Z endcap Seed
      complete(theSeed, N, NumberinEndcap, SeedinEndcap, recHitsEndcap, iSetup);
  }
  else
  {
      cout << "Not enough layers for seed reconstruction!" << endl;
  }
}


void RPCSeedGenerator::complete(RPCSeedFinder& theSeed, int N, int NumberinMuon, int *SeedinMuon, MuonRecHitContainer (&recHitsMuon)[6], const edm::EventSetup& iSetup) {

  if(NumberinMuon == 0)
      return;
  for(MuonRecHitContainer::iterator it = recHitsMuon[SeedinMuon[N]].begin(); it != recHitsMuon[SeedinMuon[N]].end(); it++)  
  {
      cout << "Complete layer[" << SeedinMuon[N] << "]." << endl;

      // Check validation
      if(!(*it)->isValid())
          continue;
            
      GlobalPoint pos = (*it)->globalPosition();
      double Phi = pos.phi();
      cout << "Phi: " << Phi << endl;
      // The recHits should locate in some phi range
      theSeed.add(*it);
      double deltaPhi = theSeed.GetdeltaPhifromHits();
      cout << "Delta phi: "<< deltaPhi << endl;
      theSeed.pop();
      if(deltaPhi > MaxDeltaPhi)
          continue;
      
      // If pass all, add to the seed
      theSeed.add(*it);
      cout << "RecHit's global position: " << pos.x() << ", " << pos.y() << ", " << pos.z() << endl;
      
      // Check if this recHit is the last one in the seed
      // If it is the last one, calculate the seed
      if(N == (NumberinMuon - 1))
      {
          cout << "Check and Fill one barrel seed." << endl;
          CheckandFill(theSeed, iSetup);
      }
      // If it is not the last one, continue to fill the seed from other layers
      else
      {
          complete(theSeed, N+1, NumberinMuon, SeedinMuon, recHitsMuon, iSetup);
      }

      // Remember to pop the recHit before add another one from the same layer!
      theSeed.pop();
  }

} 

void RPCSeedGenerator::CheckandFill(RPCSeedFinder& theSeed, const edm::EventSetup& iSetup) {
  
  if(theSeed.nrhit() >= 3)
  {
      /*
      std::vector<TrajectorySeed> the_seeds =  theSeed.seeds(iSetup, MaxRSD);
      for(std::vector<TrajectorySeed>::const_iterator the_seed=the_seeds.begin(); the_seed!=the_seeds.end(); ++the_seed)
      {
          goodSeeds.push_back(*the_seed);
      }
      */
      theSeed.seeds(iSetup, MaxRSD, goodSeeds, candidateSeeds);
  }
}

void RPCSeedGenerator::CheckOverlap(const edm::EventSetup& iSetup) {
    
    std::vector<TrajectorySeed> sortSeeds;
    std::vector<TrajectorySeed> tempSeeds;
    edm::OwnVector<TrackingRecHit> tempRecHits;

    edm::ESHandle<RPCGeometry> rpcGeometry;
    iSetup.get<MuonGeometryRecord>().get(rpcGeometry);
    
    while(goodSeeds.size() != 0)
    {
        cout << "Finding the good seeds group from " << goodSeeds.size() << " seeds which share some recHits" << endl; 
        // Take 1st seed in goodSeeds as referrence and find a collection which always share some recHits with some other
        tempRecHits.clear();
        tempSeeds.clear();
        int N = 0;
        for(vector<TrajectorySeed>::iterator itseed = goodSeeds.begin(); itseed != goodSeeds.end(); N++)
        {
            TrajectorySeed::range RecHitsRange = itseed->recHits();
            if(N == 0)
            {
                cout << "Always take the 1st good seed to be the referrence." << endl;
                for(TrajectorySeed::const_iterator it = RecHitsRange.first; it != RecHitsRange.second; it++)
                {
                    cout << "Put its recHits to tempRecHits" << endl;
                    tempRecHits.push_back(it->clone());
                }
                cout << "Put it to tempSeeds" << endl;
                tempSeeds.push_back(*itseed);
                cout << "Then erase from goodSeeds" << endl;
                itseed = goodSeeds.erase(itseed);
            }
            else
            {
                cout << "Come to other good seed for checking " << itseed->nHits() << " recHits from " << tempRecHits.size() << " temp recHits" << endl;
                unsigned int ShareRecHitsNumber = 0;
                for(TrajectorySeed::const_iterator it = RecHitsRange.first; it != RecHitsRange.second; it++)
                {
                    if(isShareHit(tempRecHits, *it, rpcGeometry))
                        ShareRecHitsNumber++;
                }
                if(ShareRecHitsNumber >= ShareRecHitsNumberThreshold)
                {
                    cout <<"This seed is found to belong to current share group" << endl;
                    for(TrajectorySeed::const_iterator it = RecHitsRange.first; it != RecHitsRange.second; it++)
                    {
                        if(!isShareHit(tempRecHits, *it, rpcGeometry))
                        {
                            cout << "Put its extra recHits to tempRecHits" << endl;
                            tempRecHits.push_back(it->clone());
                        }
                    }
                    cout << "Put it to tempSeeds" << endl;
                    tempSeeds.push_back(*itseed);
                    cout << "Then erase from goodSeeds" << endl;
                    itseed = goodSeeds.erase(itseed);
                }
                else
                    itseed++;
            }
        }
        // Find the best seed and kick out those share recHits with it
        // The best seed save in GoodSeeds, those don't share recHits with it will be push back to goodSeeds for next while loop
        TrajectorySeed bestSeed;
        vector<TrajectorySeed>::iterator bestiter;
        // Find the min Spt wrt Pt as the best Seed
        double minrsd = 1000000;
        for(vector<TrajectorySeed>::iterator itseed = tempSeeds.begin(); itseed != tempSeeds.end(); itseed++)
        {
            std::vector<float> seed_error = itseed->startingState().errorMatrix();
            double rsd = seed_error[0];
            cout << "Find a good seed with Spt_vs_Pt " << rsd << endl;
            if(rsd <= minrsd)
            {
                minrsd = rsd;
                bestSeed = *itseed;
                bestiter = itseed;
            }
        }
        cout << "Best good temp seed's Spt_vs_Pt is " << minrsd <<endl;
        sortSeeds.push_back(bestSeed);
        tempSeeds.erase(bestiter);
        tempRecHits.clear();
        
        for(TrajectorySeed::const_iterator it = bestSeed.recHits().first; it != bestSeed.recHits().second; it++)
        {
            tempRecHits.push_back(it->clone());
        }
        for(vector<TrajectorySeed>::iterator itseed = tempSeeds.begin(); itseed != tempSeeds.end(); )
        {
            cout << "Checking the temp seed's " << itseed->nHits() << " hits to " << tempRecHits.size() << " temp recHits" << endl;
            TrajectorySeed::range RecHitsRange = itseed->recHits();
            bool isShare = false;
            for(TrajectorySeed::const_iterator it = RecHitsRange.first; it != RecHitsRange.second; it++)
            {
                if(isShareHit(tempRecHits, *it, rpcGeometry))
                    isShare = true;
            }
            if(isShare == true)
            {
                cout << "Find one temp seed share some recHits with best seed" << endl;
                itseed = tempSeeds.erase(itseed);
            }
            else
            {
                cout << "This seed has no relation with best seed" << endl;
                goodSeeds.push_back(*itseed);
                itseed = tempSeeds.erase(itseed);
            }
        }
    }
    // At the end exchange goodSeeds with sortSeeds
    goodSeeds = sortSeeds;
}

void RPCSeedGenerator::CheckcandidateOverlap(const edm::EventSetup& iSetup) {

    std::vector<TrajectorySeed> sortSeeds;
    std::vector<TrajectorySeed> tempSeeds;
    edm::OwnVector<TrackingRecHit> tempRecHits;

    edm::ESHandle<RPCGeometry> rpcGeometry;
    iSetup.get<MuonGeometryRecord>().get(rpcGeometry);

    while(candidateSeeds.size() != 0)
    {
        cout << "Finding the candidate seeds group from " << candidateSeeds.size() << " seeds which share some recHits" << endl;
        // Take 1st seed in candidateSeeds as referrence and find a collection which always share some recHits with some other
        tempRecHits.clear();
        tempSeeds.clear();
        int N = 0;
        for(vector<TrajectorySeed>::iterator itseed = candidateSeeds.begin(); itseed != candidateSeeds.end(); N++)
        {
            TrajectorySeed::range RecHitsRange = itseed->recHits();
            if(N == 0)
            {
                cout << "Always take the 1st candidate seed to be the referrence." << endl;
                for(TrajectorySeed::const_iterator it = RecHitsRange.first; it != RecHitsRange.second; it++)
                {
                    cout << "Put its recHits to tempRecHits" << endl;
                    tempRecHits.push_back(it->clone());
                }
                cout << "Put it to tempSeeds" << endl;
                tempSeeds.push_back(*itseed);
                cout << "Then erase from goodSeeds" << endl;
                itseed = candidateSeeds.erase(itseed);
            }
            else
            {
                cout << "Come to other candidate seed for checking " << itseed->nHits() << " recHits from " << tempRecHits.size() << " temp recHits" << endl;
                unsigned int ShareRecHitsNumber = 0;
                for(TrajectorySeed::const_iterator it = RecHitsRange.first; it != RecHitsRange.second; it++)
                {
                    if(isShareHit(tempRecHits, *it, rpcGeometry))
                        ShareRecHitsNumber++;
                }
                if(ShareRecHitsNumber >= ShareRecHitsNumberThreshold)
                {
                    cout <<"This seed is found to belong to current share group" << endl;
                    for(TrajectorySeed::const_iterator it = RecHitsRange.first; it != RecHitsRange.second; it++)
                    {
                        if(!isShareHit(tempRecHits, *it, rpcGeometry))
                        {
                            cout << "Put its extra recHits to tempRecHits" << endl;
                            tempRecHits.push_back(it->clone());
                        }
                    }
                    cout << "Put it to tempSeeds" << endl;
                    tempSeeds.push_back(*itseed);
                    cout << "Then erase from candidateSeeds" << endl;
                    itseed = candidateSeeds.erase(itseed);
                }
                else
                    itseed++;
            }
        }
        // Find the best seed and kick out those share recHits with it
        // The best seed save in candidateSeeds, those don't share recHits with it will be push back to candidateSeeds for next while loop
        TrajectorySeed bestSeed;
        vector<TrajectorySeed>::iterator bestiter;
        // Find the min Spt wrt Pt as the best Seed
        double minrsd = 1000000;
        for(vector<TrajectorySeed>::iterator itseed = tempSeeds.begin(); itseed != tempSeeds.end(); itseed++)
        {
            std::vector<float> seed_error = itseed->startingState().errorMatrix();
            double rsd = seed_error[0];
            cout << "Find a candidate seed with Spt_vs_Pt " << rsd << endl;
            if(rsd <= minrsd)
            {   
                minrsd = rsd;
                bestSeed = *itseed;
                bestiter = itseed;
            }
        }
        cout << "Best candidate temp seed's Spt_vs_Pt is " << minrsd <<endl;
        sortSeeds.push_back(bestSeed);
        tempSeeds.erase(bestiter);
        tempRecHits.clear();

        for(TrajectorySeed::const_iterator it = bestSeed.recHits().first; it != bestSeed.recHits().second; it++)
        {
            tempRecHits.push_back(it->clone());
        }
        for(vector<TrajectorySeed>::iterator itseed = tempSeeds.begin(); itseed != tempSeeds.end(); )
        {
            cout << "Checking the temp seed's " << itseed->nHits() << " hits to " << tempRecHits.size() << " temp recHits" << endl;
            TrajectorySeed::range RecHitsRange = itseed->recHits();
            bool isShare = false;
            for(TrajectorySeed::const_iterator it = RecHitsRange.first; it != RecHitsRange.second; it++)
            {
                if(isShareHit(tempRecHits, *it, rpcGeometry))
                    isShare = true;
            }
            if(isShare == true)
                itseed = tempSeeds.erase(itseed);
            else
            {
                candidateSeeds.push_back(*itseed);
                itseed = tempSeeds.erase(itseed);
            }
        }
    }
    // At the end exchange candidateSeeds with sortSeeds
    candidateSeeds = sortSeeds;
}


bool RPCSeedGenerator::isShareHit(const edm::OwnVector<TrackingRecHit> RecHits, const TrackingRecHit& hit, edm::ESHandle<RPCGeometry> rpcGeometry) {
    
    bool istheSame = false;
    unsigned int n = 1;
    cout << "Checking from " << RecHits.size() << " temp recHits" << endl;
    
    LocalPoint lpos1 = hit.localPosition();
    DetId RPCId1 = hit.geographicalId();
    const GeomDetUnit *rpcroll1 = rpcGeometry->idToDetUnit(RPCId1);
    GlobalPoint gpos1 = rpcroll1->toGlobal(lpos1);
    cout << "The hit's position: " << gpos1.x() << ", " << gpos1.y() << ", " << gpos1.z() << endl;
    for(edm::OwnVector<TrackingRecHit>::const_iterator it = RecHits.begin(); it !=RecHits.end(); it++, n++)
    {
        cout << "Checking the " << n << " th recHit from tempRecHits" << endl;
        LocalPoint lpos2 = it->localPosition();
        DetId RPCId2 = it->geographicalId();
        const GeomDetUnit *rpcroll2 = rpcGeometry->idToDetUnit(RPCId2);
        GlobalPoint gpos2 = rpcroll2->toGlobal(lpos2);
        cout << "The temp hit's position: " << gpos2.x() << ", " << gpos2.y() << ", " << gpos2.z() << endl;
        
        if((gpos1.x() == gpos2.x()) && (gpos1.y() == gpos2.y()) && (gpos1.z() == gpos2.z()))
        {
            cout << "This hit is found to be the same" << endl;
            istheSame = true;
        }
    }
    return istheSame;
}

//define this as a plug-in
DEFINE_FWK_MODULE(RPCSeedGenerator);

