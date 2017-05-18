#ifndef TTracking_hxx_seen
#define TTracking_hxx_seen

#include <vector>
#include <map>

#include <TObject.h>
#include <TVector3.h>
#include <TRandom.h>

#include <IAlgorithm.hxx>
#include <ICOMETLog.hxx>

#include <IG4Trajectory.hxx>
#include <IG4VHit.hxx>
#include <IHandle.hxx>
#include <ICOMETEvent.hxx>

#include "TGeoManager.h"
#include "TGeoNode.h"

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

TVector3 GetPOCAofTwoWires(TVector3 wireEnd0_lo, TVector3 wireEnd1_lo, TVector3 wireEnd0_up, TVector3 wireEnd1_up);
TVector3 GetVectorCrossingCenter(TVector3 wireEnd0_lo, TVector3 wireEnd1_lo, TVector3 wireEnd0_up, TVector3 wireEnd1_up, TVector3 POCA);

class ITracking{
private:
  
protected:
  // Input Members
  double fGenTrX;
  double fGenTrY;
  double fGenTrZ;
  double fGenTrT;
  double fGenTrPx;
  double fGenTrPy;
  double fGenTrPz;
  double fGenTrE;

  double fCDCEnterX;
  double fCDCEnterY;
  double fCDCEnterZ;
  double fCDCEnterT;
  double fCDCEnterPx;
  double fCDCEnterPy;
  double fCDCEnterPz;

  int fnCDCHit;
  double fCDCHitX[10000];
  double fCDCHitY[10000];
  double fCDCHitZ[10000];
  double fCDCHitT[10000];
  double fCDCEDep[10000];

  int    fnCALCDCHit;
  double fDriftDist[10000];
  int    fCDCCharge[10000];
  double fWireEnd0X[10000];
  double fWireEnd0Y[10000];
  double fWireEnd0Z[10000];
  double fWireEnd1X[10000];
  double fWireEnd1Y[10000];
  double fWireEnd1Z[10000];
  int    fWireLayerId[10000];
  int    fWireId[10000];
  int    fMaxWireLayerId;
  int    fTurnId[10000];
  int    fTurnNumber;

  // Temporal member
  int    fDomain[10000];
  bool   fSide[10000];
  bool   fOuter[10000];
  bool   fInner[10000];

public:
  ITracking(const char*name, const char* title);
  virtual ~ITracking();
  
  /// called at the begin of run or else (should not be in event-by-event)
  int  Init();

  void LoadMCHits(COMET::IHandle<COMET::IHitSelection> hitHandle, COMET::IHandle<COMET::IG4TrajectoryContainer> trajectories, COMET::IHandle<COMET::IG4HitContainer> cdcHits);
  
  void ShuffleMCHits();

  void PrintMCStatus();

  int GetTurnNumber() { return fTurnNumber; }

  /// called at the end of run or else (should not be in event-by-event)
  int  Finish();

};
#endif
