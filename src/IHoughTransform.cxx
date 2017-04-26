#include <IHoughTransform.hxx>

#include "HEPUnits.hxx"
#include <vector>
#include <iostream>

#include <IOADatabase.hxx>
#include <IReconBase.hxx>
#include <IReconTrackCand.hxx>
#include <ITrackState.hxx>

#include <ICOMETEvent.hxx>
#include <ICOMETLog.hxx>

#include <IHitSelection.hxx>
#include <IMCHit.hxx>
#include <IG4HitSegment.hxx>
#include <IG4HitCalo.hxx>
#include <IHit.hxx>
#include <IGeomInfo.hxx>
#include <IGeometryId.hxx>
#include <COMETGeomId.hxx>
#include <IGeomIdManager.hxx>
#include <ICTHGeomId.hxx>
#include <ICTHChannelId.hxx>
#include <IChannelId.hxx>
#include <IGeometryDatabase.hxx>
#include <ICDCGeom.hxx>

IHoughTransform::IHoughTransform(const char* name = "IHoughTransform", const char* title="houghtransform"):ITracking(name,title)
{;}
IHoughTransform::~IHoughTransform(){;}

int IHoughTransform::Init()
{
  COMETNamedInfo("IHoughTransform", "//----------------------------------------------------------//");
  COMETNamedInfo("IHoughTransform", "//   Initialize IHoughTransform                               ");
  COMETNamedInfo("IHoughTransform", "//----------------------------------------------------------//");
  return 1;
}

void IHoughTransform::Clear(){
  /*
  fTime.clear();  
  fIndex.clear();
  fScint.clear();
  fModule.clear();
  fDSTimeCluster.clear();
  fUSTimeCluster.clear();
  fPairCandidates.clear();
  */
}

int IHoughTransform::Finish(){
  return 1;
}
