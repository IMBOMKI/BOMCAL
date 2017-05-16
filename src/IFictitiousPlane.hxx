#ifndef TFictitiousPlane_hxx_seen
#define TFictitiousPlane_hxx_seen

#include <ITracking.hxx>
#include <IHoughTransform.hxx>

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


struct SingleHit{
  TVector3 wireEnd0;
  TVector3 wireEnd1;
  double driftDist;
  int wireId;
  int layerId;
  int domain;
};

struct HitPair{
  SingleHit h1;
  SingleHit h2;
  TVector3  cV;
};


class IFictitiousPlane: public ITracking {
public:
  IFictitiousPlane(const char*name, const char* title);
  virtual ~IFictitiousPlane();

  /// called at the begin of run or else (should not be in event-by-event)
  int  Init();

  int BeginOfEvent();

  int EndOfEvent();

  void LoadMCHits(COMET::IHandle<COMET::IHitSelection> hitHandle, COMET::IHandle<COMET::IG4TrajectoryContainer> trajectories,COMET::IHandle<COMET::IG4HitContainer> cdcHits){
    ITracking::LoadMCHits(hitHandle, trajectories, cdcHits);
  }

  void LoadHitsAfterHT(COMET::IHandle<COMET::IHitSelection> hitHandle, IHoughTransform* hough);
  void AddRandomHitPairs(int n, int domain);
  void AddSideHitPairs(int n, int domain);
  void DrawHitsOnFictitiousPlane(TCanvas* canvas);

private:
  std::vector <HitPair> fHitPairs;
  
  //std::vector <HitPair> fHitPairsDomain2;
};
#endif
