#ifndef THoughTransform_hxx_seen
#define THoughTransform_hxx_seen

#include <ITracking.hxx>

#include <vector>
#include <map>

#include <TObject.h>
#include <TVector3.h>
#include <TRandom.h>

#include <ICOMETLog.hxx>

#include <IG4Trajectory.hxx>
#include <IG4VHit.hxx>
#include <IHandle.hxx>
#include <ICOMETEvent.hxx>

#include "TGeoManager.h"
#include "TGeoNode.h"


class IHoughTransform: public ITracking {
private:
  
public:
  IHoughTransform(const char*name, const char* title);
  virtual ~IHoughTransform();
  
  /// called at the begin of run or else (should not be in event-by-event)
  int  Init();

  void Clear();

  /// called at the end of run or else (should not be in event-by-event)
  int  Finish();

};
#endif
