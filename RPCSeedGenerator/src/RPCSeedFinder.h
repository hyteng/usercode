#ifndef RecoMuon_MuonSeedGenerator_RPCSeedFinder_H
#define RecoMuon_MuonSeedGenerator_RPCSeedFinder_H

/** \class RPCSeedFinder
 *  
 *
 *
 *  
 */

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"

#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"

#include <vector>
#include <algorithm>

namespace edm {class EventSetup;}

class RPCSeedFinder {
public:
  
  RPCSeedFinder();

  virtual ~RPCSeedFinder(){};

  void add(MuonTransientTrackingRecHit::MuonRecHitPointer hit) { theRecHits.push_back(hit); }
  void pop() { theRecHits.pop_back(); }
  
  //std::vector<TrajectorySeed> seeds(const edm::EventSetup& eSetup, double MaxRSD) const;
  void seeds(const edm::EventSetup& eSetup, double MaxRSD, std::vector<TrajectorySeed>& goodSeeds, std::vector<TrajectorySeed>& candidateSeeds) const;
  
  MuonTransientTrackingRecHit::ConstMuonRecHitPointer firstRecHit() const { return theRecHits.front(); }
  unsigned int nrhit() const { return theRecHits.size(); }
  double GetdeltaPhifromHits();
private:
    
  MuonTransientTrackingRecHit::MuonRecHitContainer theRecHits;
 
};
#endif
