/*
 *  See header file for a description of this class.
 *
 *
 *  $Date: 2007/06/08 11:59:46 $
 *  $Revision: 1.1 $
 *  \author D. Pagano - University of Pavia & INFN Pavia
 *
 */

#include "MyMTCCAnalyzer/RPCSeedGenerator/src/RPCSeedHits.h"

#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"

#include "RecoMuon/TrackingTools/interface/MuonPatternRecoDumper.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/DetLayers/interface/DetLayer.h"

#include "DataFormats/TrajectoryState/interface/PTrajectoryStateOnDet.h"
#include "DataFormats/Common/interface/OwnVector.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "gsl/gsl_statistics.h"
#include "TH1F.h"
#include "math.h"

using namespace std;

typedef MuonTransientTrackingRecHit::MuonRecHitPointer MuonRecHitPointer;
typedef MuonTransientTrackingRecHit::ConstMuonRecHitPointer ConstMuonRecHitPointer;
typedef MuonTransientTrackingRecHit::MuonRecHitContainer MuonRecHitContainer;

template <class T> T sqr(const T& t) { return t*t; }

TrajectorySeed RPCSeedHits::seed(const edm::EventSetup& eSetup, double MaxRSD, int& isGoodSeed) const
{  

  unsigned int NumberofHitsinSeed = nrhit();
  double *pt = 0;
  double *spt = 0;
  int NumberofPt = 0;
  
  if(NumberofHitsinSeed < 3)
  {
      isGoodSeed = -1;
      // Create a fake seed and return
      ConstMuonRecHitPointer last = (*(theRhits.begin()));
      return createSeed(0, 0, last, eSetup);
  }
  
  cout << "computePtWithoutVtx" << endl;
  NumberofPt = NumberofHitsinSeed * (NumberofHitsinSeed - 1) * (NumberofHitsinSeed - 2) / (3 * 2);
  pt = new double[NumberofPt];
  spt = new double[NumberofPt];
  computePtWithoutVtx(pt, spt);
  
  float ptmean=0.;
  float sptmean=0.;

  computeBestPt(pt, spt, ptmean, sptmean, NumberofPt);
  ConstMuonRecHitPointer last = best_cand();
  
  isGoodSeed = false;
  if(sptmean <= fabs(MaxRSD * ptmean))
  {
      cout << "III--> Seed good Pt : " << ptmean << endl;
      cout << "III--> Seed good Spt : " << ptmean << endl;
      cout << "Spt_vs_Pt : " << sptmean/ptmean << endl;
      isGoodSeed = true;
  }
  else
  {
      cout << "III--> Seed bad Pt : " << ptmean << endl;
      cout << "III--> Seed bad Spt : " << ptmean << endl;
      cout << "Spt_vs_Pt : " << sptmean/ptmean << endl;
      isGoodSeed = false;
  }
  
  delete [] pt;
  delete [] spt;
  return createSeed(ptmean, sptmean, last, eSetup);
}


ConstMuonRecHitPointer RPCSeedHits::best_cand() const {

  MuonRecHitPointer best = 0;

  for (MuonRecHitContainer::const_iterator iter=theRhits.begin();
       iter!=theRhits.end(); iter++) {
       best = (*iter);
  } 

  return best;
}


void RPCSeedHits::computePtWithoutVtx(double* pt, double* spt) const 
{

  int n = 0;
  float *x;
  float *y;
  unsigned int NumberofHitsinSeed = nrhit();
  x = new float[NumberofHitsinSeed];
  y = new float[NumberofHitsinSeed];
  for(MuonRecHitContainer::const_iterator iter=theRhits.begin(); iter!=theRhits.end(); iter++) 
  {
    cout << "X[" << n <<"] = " << (*iter)->globalPosition().x() << ", Y[" << n <<"]= " << (*iter)->globalPosition().y() << endl;
    x[n] = (*iter)->globalPosition().x();
    y[n] = (*iter)->globalPosition().y();
    n++;
  }
  
  n = 0;
  for(int i = 0; i < (NumberofHitsinSeed - 2); i++)
      for(int j = (i + 1); j < (NumberofHitsinSeed - 1); j++)
          for(int k = (j + 1); k < NumberofHitsinSeed; k++)
          {  
              float A = (y[k]-y[j])/(x[k]-x[j]) - (y[j]-y[i])/(x[j]-x[i]);
              float TYO = (x[k]-x[i])/A + (y[k]*y[k]-y[j]*y[j])/((x[k]-x[j])*A) - (y[j]*y[j]-y[i]*y[i])/((x[j]-x[i])*A);
              float TXO = (x[k]+x[j]) + (y[k]*y[k]-y[j]*y[j])/(x[k]-x[j]) - TYO*(y[k]-y[j])/(x[k]-x[j]);
              float XO = 0.5 * TXO;
              float YO = 0.5 * TYO;
              float R2 = (x[i]-XO)*(x[i]-XO) + (y[i]-YO)*(y[i]-YO);
              // How this algorithm get the pt without magnetic field??
              pt[n] = 0.01 * sqrt(R2) * 2 * 0.3;
              n++;
          }
  delete [] x;
  delete [] y;
     
}


void RPCSeedHits::computeBestPt(double* pt, double* spt, float& ptmean, float& sptmean, int NumberofPt) const {
  cout << "[RPCSeedHits] --> computeBestPt class called." << endl;
  double ptall = 0;
  cout << "---< best pt computing >---" << endl;
  for(int i = 0; i < NumberofPt; i++)
  {
      cout << "pt[" << i <<"] = " << pt[i] << endl;
      ptall += pt[i];
  }
  cout << "---------------------------" << endl;

  ptmean = ptall / NumberofPt;
  
  sptmean = 0;
  for(int i = 0; i < NumberofPt; i++)
  {
      sptmean += (pt[i] - ptmean) * (pt[i] - ptmean);
  }
  sptmean /= NumberofPt;
  sptmean = sqrt(sptmean);
  
  cout << "Ptmean : " << ptmean << endl;
  cout << "Sptmean : " << sptmean << endl;
}


TrajectorySeed RPCSeedHits::createSeed(float ptmean, float sptmean, ConstMuonRecHitPointer last, const edm::EventSetup& eSetup) const{
  
  MuonPatternRecoDumper debug;
  
  edm::ESHandle<MagneticField> field;
  eSetup.get<IdealMagneticFieldRecord>().get(field);

  double theMinMomentum = 3.0;
 
  if ( fabs(ptmean) < theMinMomentum ) ptmean = theMinMomentum * ptmean/fabs(ptmean) ;

  AlgebraicVector t(4);
  AlgebraicSymMatrix mat(5,0) ;

  LocalPoint segPos=last->localPosition();
  GlobalVector mom=last->globalPosition()-GlobalPoint();
  GlobalVector polar(GlobalVector::Spherical(mom.theta(), last->globalDirection().phi(), 1.));
  polar *= fabs(ptmean)/polar.perp();
  LocalVector segDirFromPos=last->det()->toLocal(polar);
  int charge=(int)(ptmean/fabs(ptmean));

  LocalTrajectoryParameters param(segPos,segDirFromPos, charge);

  mat = last->parametersError().similarityT( last->projectionMatrix() );
  
  //float p_err = sqr(sptmean/(ptmean*ptmean));
  float p_err = sptmean/abs(ptmean);
  mat[0][0]= p_err;
 
  LocalTrajectoryError error(mat);
  
  TrajectoryStateOnSurface tsos(param, error, last->det()->surface(),&*field);

  cout << "Trajectory State on Surface before the extrapolation"<<endl;
  cout << debug.dumpTSOS(tsos);
  
  DetId id = last->geographicalId();

  cout << "The RecSegment relies on: "<<endl;
  cout << debug.dumpMuonId(id);
  cout << debug.dumpTSOS(tsos);

  TrajectoryStateTransform tsTransform;
  
  PTrajectoryStateOnDet *seedTSOS = tsTransform.persistentState(tsos ,id.rawId());
  
  edm::OwnVector<TrackingRecHit> container;
  for(MuonRecHitContainer::const_iterator iter=theRhits.begin(); iter!=theRhits.end(); iter++)
  {
      // This casting withou clone will cause memory overflow when used in push_back
      // Since container's deconstructor functiion free the pointer menber!
      //TrackingRecHit* pt = dynamic_cast<TrackingRecHit*>(&*(*iter));
      //cout << "Push recHit type " << pt->getType() << endl;
      container.push_back((*iter)->hit()->clone());
  }
  TrajectorySeed theSeed(*seedTSOS, container, oppositeToMomentum);
   
  delete seedTSOS;
  return theSeed;
}
