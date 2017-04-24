#ifndef TMCTrigger_hxx_seen
#define TMCTrigger_hxx_seen

#include <vector>

#include <TObject.h>
#include <TVector3.h>
#include <TRandom.h>

#include <IAlgorithm.hxx>
#include <ICOMETLog.hxx>

#include <IG4VHit.hxx>
#include <IHandle.hxx>
#include <ICOMETEvent.hxx>

#include "TGeoManager.h"
#include "TGeoNode.h"


class IMCTrigger{
private:
  

public:
  IMCTrigger(const char*name, const char* title);
  ~IMCTrigger();
  
  /// called at the begin of run or else (should not be in event-by-event)
  int  Init();

  TGeoNode* GetNode(TVector3 position);

  /// Calculate Number of Photons
  int GetPhotonNumber(int pdgNum, double ene, double trLen);

  /// called at the end of run or else (should not be in event-by-event)
  int  Finish();

  COMET::IG4HitContainer* MakeCTHHitSelection(COMET::IHandle<COMET::IG4HitContainer> & cthhits);
};
#endif
