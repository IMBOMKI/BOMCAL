#include <vector>
#include <iostream>
#include <assert.h>

#include <IOADatabase.hxx>
#include <IHelixTracker.hxx>
#include <IReconTrack.hxx>
#include <IReconTrackCand.hxx>

#include <ICOMETLog.hxx>
#include <ICOMETEvent.hxx>
#include <IHitSelection.hxx>
#include <IMCHit.hxx>
#include <IHit.hxx>
#include <IGeomInfo.hxx>
#include <IFieldManager.hxx>
#include <IGeoField.hxx>

#include <TCanvas.h>
#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>
#include "TEllipse.h"

//TVector3 GetPOCAofTwoWires(TVector3 wireEnd0_lo, TVector3 wireEnd1_lo, TVector3 wireEnd0_up, TVector3 wireEnd1_up);
//TVector3 GetVectorCrossingCenter(TVector3 wireEnd0_lo, TVector3 wireEnd1_lo, TVector3 wireEnd0_up, TVector3 wireEnd1_up, TVector3 POCA);

IHelixTracker::IHelixTracker(const char* name, const char* title)
  :ITracking(name, title)
{;}

IHelixTracker::~IHelixTracker()
{;}

int IHelixTracker::Init(){
  COMETInfo("//----------------------------------------------------------//");
  COMETInfo("// Initialize IHelixTracker");

  return 0;
}

int IHelixTracker::BeginOfEvent(){
  /// Field Map
  return 0;
}

int IHelixTracker::EndOfEvent()
{
  return 0;
}    

void IHelixTracker::LoadHitsAfterHT(COMET::IHandle<COMET::IHitSelection> hitHandle, IHoughTransform* hough){  
  std::vector<int> wireId       = hough->GetRecoWireId();
  std::vector<double> driftDist = hough->GetRecoDriftDist();
  std::vector<int> domain       = hough->GetRecoDomain();
  std::vector<bool> outer       = hough->GetRecoOuterHit();
  std::vector<bool> inner       = hough->GetRecoInnerHit();
  fMaxWireLayerId               = hough->GetMaxWireLayerId();

  fnCALCDCHit=wireId.size();

  // Clear the arrays
  memset (fWireId,-1,sizeof(fWireId));
  memset (fWireLayerId,-1,sizeof(fWireLayerId));
  memset (fWireEnd0X,-1,sizeof(fWireEnd0X));
  memset (fWireEnd0Y,-1,sizeof(fWireEnd0Y));
  memset (fWireEnd0Z,-1,sizeof(fWireEnd0Z));
  memset (fWireEnd1X,-1,sizeof(fWireEnd1X));
  memset (fWireEnd1Y,-1,sizeof(fWireEnd1Y));
  memset (fWireEnd1Z,-1,sizeof(fWireEnd1Z));
  memset (fDriftDist,-1,sizeof(fDriftDist));
  memset (fDomain,-1,sizeof(fDomain));
  memset (fOuter,-1,sizeof(fOuter));
  memset (fInner,-1,sizeof(fInner));

  assert(fnCALCDCHit==hough->GetNumberOfRecognizedHits());
  assert(wireId.size() == driftDist.size());
  assert(wireId.size() == domain.size());

  for (int i=0; i<fnCALCDCHit; i++){
    int wire = wireId.at(i);
    fWireId[i]      = wire;
    fWireLayerId[i] = COMET::IGeomInfo::Get().CDC().GetLayer(wire);
    fWireEnd0X[i]   = COMET::IGeomInfo::Get().CDC().GetWireEnd0(wire).X();
    fWireEnd0Y[i]   = COMET::IGeomInfo::Get().CDC().GetWireEnd0(wire).Y();
    fWireEnd0Z[i]   = COMET::IGeomInfo::Get().CDC().GetWireEnd0(wire).Z();
    fWireEnd1X[i]   = COMET::IGeomInfo::Get().CDC().GetWireEnd1(wire).X();
    fWireEnd1Y[i]   = COMET::IGeomInfo::Get().CDC().GetWireEnd1(wire).Y();
    fWireEnd1Z[i]   = COMET::IGeomInfo::Get().CDC().GetWireEnd1(wire).Z();
    fDriftDist[i]   = driftDist.at(i);
    fDomain[i]      = domain.at(i);
    fOuter[i]       = outer.at(i);
    fInner[i]       = inner.at(i);
    //std::cout << i << "   " << fWireEnd0X[i] << "   " << fWireEnd0Y[i] << "   " << fWireEnd0Z[i] << "   " << fDomain[i] << std::endl;
  }
  std::cout << "Initial Position: " << fCDCEnterX << "  " << fCDCEnterY << "  " << fCDCEnterZ << std::endl;
  std::cout << "Initial Momentum: " << fCDCEnterPx << "  " << fCDCEnterPy << "  " << fCDCEnterPz << std::endl;
}  

void IHelixTracker::AddSideHitPairs(int n, int domain){
  assert(n<=fMaxWireLayerId);
  for (int i=0; i<n; i++){ 
   
    struct SingleHit hit_lo, hit_up;
    struct HitPair hitPair;

    int i_lo=-1;
    int i_up=-1;
    
    for (int i_hit=0; i_hit<fnCALCDCHit; i_hit++){
      if      (fWireLayerId[i_hit]==i   && fDomain[i_hit]==domain && fInner[i_hit]==1) i_lo=i_hit;
      else if (fWireLayerId[i_hit]==i+1 && fDomain[i_hit]==domain && fOuter[i_hit]==1) i_up=i_hit;
    }
    
    if (i_lo==-1 || i_up==-1) continue;

    TVector3 wireEnd0_lo = TVector3(fWireEnd0X[i_lo], fWireEnd0Y[i_lo], fWireEnd0Z[i_lo]);
    TVector3 wireEnd1_lo = TVector3(fWireEnd1X[i_lo], fWireEnd1Y[i_lo], fWireEnd1Z[i_lo]);

    TVector3 wireEnd0_up = TVector3(fWireEnd0X[i_up], fWireEnd0Y[i_up], fWireEnd0Z[i_up]);
    TVector3 wireEnd1_up = TVector3(fWireEnd1X[i_up], fWireEnd1Y[i_up], fWireEnd1Z[i_up]);

    hit_lo.wireEnd0  = wireEnd0_lo;
    hit_lo.wireEnd1  = wireEnd1_lo;
    hit_lo.driftDist = fDriftDist[i_lo];
    hit_lo.wireId    = fWireId[i_lo];
    hit_lo.layerId   = fWireLayerId[i_lo];
    hit_lo.domain    = fDomain[i_lo];

    hit_up.wireEnd0  = wireEnd0_up;
    hit_up.wireEnd1  = wireEnd1_up;
    hit_up.driftDist = fDriftDist[i_up];
    hit_up.wireId    = fWireId[i_up];
    hit_up.layerId   = fWireLayerId[i_up];
    hit_up.domain    = fDomain[i_up];

    TVector3 POCA = GetPOCAofTwoWires(hit_lo.wireEnd0, hit_lo.wireEnd1, hit_up.wireEnd0, hit_up.wireEnd1);
    TVector3 cVec = GetVectorCrossingCenter(hit_lo.wireEnd0, hit_lo.wireEnd1, hit_up.wireEnd0, hit_up.wireEnd1,POCA);

    hitPair.h1 = hit_lo;
    hitPair.h2 = hit_up;
    hitPair.cV = cVec;

    fHitPairs.push_back(hitPair);
    fPOCAs.push_back(POCA);

    std::cout << "Layer: <" << i << " " << i+1 << ">   Domain: " << domain << std::endl;
    std::cout << "POCA:  " << POCA(0) << "  " << POCA(1) << "  " << POCA(2) << std::endl;
  }
}

