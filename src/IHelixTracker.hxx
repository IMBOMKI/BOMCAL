#ifndef THelixTracker_hxx_seen
#define THelixTracker_hxx_seen

#include <ITracking.hxx>
#include <IHoughTransform.hxx>
#include <IFictitiousPlane.hxx>

#include <vector>
#include <map>

#include <TVirtualFitter.h>
#include <TObject.h>
#include <TVector3.h>
#include <TRandom.h>
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TFile.h"

#include <ICOMETLog.hxx>

#include <IG4Trajectory.hxx>
#include <IG4VHit.hxx>
#include <IHandle.hxx>
#include <ICOMETEvent.hxx>
#include <IReconTrack.hxx>
#include <IReconTrackCand.hxx>
#include <IFieldManager.hxx>

class IHelixTracker: public ITracking {
public:
  IHelixTracker(const char*name, const char* title);
  virtual ~IHelixTracker();

  /// called at the begin of run or else (should not be in event-by-event)
  int  Init();

  int BeginOfEvent();

  int EndOfEvent();

  void LoadMCHits(COMET::IHandle<COMET::IHitSelection> hitHandle, COMET::IHandle<COMET::IG4TrajectoryContainer> trajectories,COMET::IHandle<COMET::IG4HitContainer> cdcHits){
    ITracking::LoadMCHits(hitHandle, trajectories, cdcHits);
  }

  void LoadHitsAfterHT(COMET::IHandle<COMET::IHitSelection> hitHandle, IHoughTransform* hough);
  void AddSideHitPairs(int n, int domain);
  std::vector <TVector3> GetPOCAs() { return fPOCAs; }
  TVector3 GetEnterPos() {return TVector3(fCDCEnterX,fCDCEnterY,fCDCEnterZ);}
  TVector3 GetEnterMom() {return TVector3(fCDCEnterPx,fCDCEnterPy,fCDCEnterPz);}  
  

private:
  std::vector <HitPair> fHitPairs;
  std::vector <TVector3> fPOCAs;

  //std::vector <HitPair> fHitPairsDomain2;
};
#endif
