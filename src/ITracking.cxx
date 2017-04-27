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

#include <TGeoManager.h>
#include <TGeoNode.h>

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

void ITracking::LoadMCHits(COMET::IHandle<COMET::IHitSelection> hitHandle, COMET::IHandle<COMET::IG4TrajectoryContainer> trajectories){

  TGeoManager* geom = COMET::IOADatabase::Get().Geometry();
  COMET::IG4TrajectoryContainer *TrajCont = GetPointer(trajectories);
  
  if(TrajCont->empty()){
    std::cout<< "Result is not found" << std::endl;
  }
  
  if(!TrajCont->empty()){
    
    for(COMET::IG4TrajectoryContainer::const_iterator seg = TrajCont->begin(); seg != TrajCont->end(); seg++){      
      COMET::IG4Trajectory traj = (*seg).second;      
      /*------  Primary Particle -----*/                 
      if (traj.GetTrackId() == 1){
	TVector3 iniPos = traj.GetInitialPosition().Vect()*(1/unit::cm);
	TVector3 iniMom = traj.GetInitialMomentum().Vect()*(1/unit::MeV);	
	fGenTrX=iniPos(0);
	fGenTrY=iniPos(1);
	fGenTrZ=iniPos(2);
	fGenTrT=traj.GetInitialPosition()(3);	
	fGenTrPx=iniMom(0);
	fGenTrPy=iniMom(1);
	fGenTrPz=iniMom(2);
	fGenTrE=traj.GetInitialMomentum()(3);	      
      }
    }
  }
  
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

int ITracking::Finish(){
  return 1;
}
