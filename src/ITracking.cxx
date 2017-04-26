#include <ITracking.hxx>

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

ITracking::ITracking(const char* name = "ITracking", const char* title="tracking")
{;}
ITracking::~ITracking(){;}

int ITracking::Init()
{
  COMETNamedInfo("ITracking", "//----------------------------------------------------------//");
  COMETNamedInfo("ITracking", "//   Initialize ITracking                               ");
  COMETNamedInfo("ITracking", "//----------------------------------------------------------//");
  return 1;
}

void ITracking::LoadMCHits(COMET::IHandle<COMET::IHitSelection> hitHandle){

  COMET::IChannelId tmpchanId;

  if (hitHandle){
    COMET::IHitSelection *hits = GetPointer(hitHandle);
    for (COMET::IHitSelection::const_iterator hitSeg = hits->begin(); hitSeg != hits->end(); hitSeg++){
      
      COMET::IChannelId chanId = (*hitSeg)->GetChannelId();
      COMET::IGeometryId geomId = (*hitSeg)->GetGeomId();
      int wire = COMET::IGeomInfo::Get().CDC().GeomIdToWire(geomId);
      
      TVectorD wireMes(7);
      wireMes[0] = COMET::IGeomInfo::Get().CDC().GetWireEnd0(wire).X();
      wireMes[1] = COMET::IGeomInfo::Get().CDC().GetWireEnd0(wire).Y();
      wireMes[2] = COMET::IGeomInfo::Get().CDC().GetWireEnd0(wire).Z();
      wireMes[3] = COMET::IGeomInfo::Get().CDC().GetWireEnd1(wire).X();
      wireMes[4] = COMET::IGeomInfo::Get().CDC().GetWireEnd1(wire).Y();
      wireMes[5] = COMET::IGeomInfo::Get().CDC().GetWireEnd1(wire).Z();
      wireMes[6] = (*hitSeg)->GetDriftDistance();
      
      TVector3 wireend0(wireMes[0],wireMes[1],wireMes[2]);	
      TVector3 wireend1(wireMes[3],wireMes[4],wireMes[5]);
      Double_t Drfit = (*hitSeg)->GetDriftDistance();
      if(!COMET::IGeomInfo::DetectorSolenoid().GetDetPositionInDSCoordinate(wireend0, wireend0)){
	continue;
	std::cout << "MisIdentifided wire is detected" << std::endl;}
      if(!COMET::IGeomInfo::DetectorSolenoid().GetDetPositionInDSCoordinate(wireend1, wireend1)){
	continue;
	std::cout << "MisIdentifided wire is detected" << std::endl;}
      
      int layer = COMET::IGeomInfo::Get().CDC().GetLayer(wire);
      int fWireMaxLayerId;
      if (fWireMaxLayerId<layer){
	fWireMaxLayerId = layer;}
      
      int wireid;
      if(COMET::GeomId::CDC::IsActiveSenseWire(geomId)==1){
	wireid = COMET::IGeomInfo::Get().CDC().GeomIdToWire(geomId);
      }
      else {
	std::cout << "Not SenseWire" << std::endl;
      }
      
      ///////////////////////////CALIBRATED HIT//////////////////////////////////////                
      
      if (hitSeg == hits->begin()){
	tmpchanId = chanId;
	fCDCCharge[fnCALCDCHit]++;
	fWireEnd0X[fnCALCDCHit]=wireend0(0);
	fWireEnd0Y[fnCALCDCHit]=wireend0(1);
	fWireEnd0Z[fnCALCDCHit]=wireend0(2);
	fWireEnd1X[fnCALCDCHit]=wireend1(0);
	fWireEnd1Y[fnCALCDCHit]=wireend1(1);
	fWireEnd1Z[fnCALCDCHit]=wireend1(2);
	fDriftDist[fnCALCDCHit]=wireMes[6];
	fWireLayerId[fnCALCDCHit]=layer;
	fWireId[fnCALCDCHit]=wireid;
      }
      
      // When the next hit is at same Id                                                             
      if (chanId == tmpchanId){
	fCDCCharge[fnCALCDCHit]++;
      }
      
      // When the next hit generates at another Channel Id                                           
      else if (chanId != tmpchanId){
	tmpchanId = chanId;
	fnCALCDCHit++;
	fCDCCharge[fnCALCDCHit]++;
	fWireEnd0X[fnCALCDCHit]=wireend0(0);
	fWireEnd0Y[fnCALCDCHit]=wireend0(1);
	fWireEnd0Z[fnCALCDCHit]=wireend0(2);
	fWireEnd1X[fnCALCDCHit]=wireend1(0);
	fWireEnd1Y[fnCALCDCHit]=wireend1(1);
	fWireEnd1Z[fnCALCDCHit]=wireend1(2);
	fDriftDist[fnCALCDCHit]=wireMes[6];
	fWireLayerId[fnCALCDCHit]=layer;
	fWireId[fnCALCDCHit]=wireid;
      }
    }
  }
}

void ITracking::PrintMCStatus(){
  std::cout << "Print Something..." << std::endl;
}

void ITracking::Clear(){
  fnCALCDCHit=0;
  memset(fDriftDist,0,sizeof(fDriftDist));
  memset(fCDCCharge,0,sizeof(fCDCCharge));
  memset(fWireEnd0X,0,sizeof(fWireEnd0X));  
  memset(fWireEnd0Y,0,sizeof(fWireEnd0Y));  
  memset(fWireEnd0Z,0,sizeof(fWireEnd0Z));  
  memset(fWireEnd1X,0,sizeof(fWireEnd1X));  
  memset(fWireEnd1Y,0,sizeof(fWireEnd1Y));  
  memset(fWireEnd1Z,0,sizeof(fWireEnd1Z));  
  memset(fWireLayerId,0,sizeof(fWireLayerId));  
  memset(fWireId,0,sizeof(fWireId));  
  fWireMaxLayerId=0;
}

int ITracking::Finish(){
  return 1;
}
