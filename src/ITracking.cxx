#include <ITracking.hxx>

#include "HEPUnits.hxx"
#include <vector>
#include <iostream>
#include <string>
#include <algorithm>

#include <IOADatabase.hxx>
#include <IReconBase.hxx>
#include <IReconTrackCand.hxx>
#include <ITrackState.hxx>

#include <ICOMETEvent.hxx>
#include <ICOMETLog.hxx>

#include <IHitSelection.hxx>
#include <IMCHit.hxx>
#include <IG4VHit.hxx>
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
#include <ICDCChannelMap.hxx>
#include <ICDCWireManager.hxx>

#include <TGeoManager.h>
#include <TGeoNode.h>

ITracking::ITracking(const char* name = "ITracking", const char* title="tracking")
:fnCALCDCHit(0), 
  fnCDCHit(0),
 fTurnNumber(0)

{;}
ITracking::~ITracking(){;}

int ITracking::Init()
{
  COMETNamedInfo("ITracking", "//----------------------------------------------------------//");
  COMETNamedInfo("ITracking", "//   Initialize ITracking                               ");
  COMETNamedInfo("ITracking", "//----------------------------------------------------------//");
  return 1;
}

TGeoNode* GetNode(TVector3 position) {
  gGeoManager->PushPath();	
  gGeoManager->GetTopNode()->cd();	 	  
  TGeoNode* volume = gGeoManager->FindNode(position(0),position(1),position(2));
  gGeoManager->PopPath();
  return volume;
}

void ITracking::LoadMCHits(COMET::IHandle<COMET::IHitSelection> hitHandle, COMET::IHandle<COMET::IG4TrajectoryContainer> trajectories, COMET::IHandle<COMET::IG4HitContainer> cdcHits){

  TGeoManager* geom = COMET::IOADatabase::Get().Geometry();
  COMET::ICDCWireManager* WireManager = new COMET::ICDCWireManager();

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

	COMET::IG4Trajectory::Points trajPointSet = traj.GetTrajectoryPoints();
	for(COMET::IG4Trajectory::Points::iterator trajIter = trajPointSet.begin(); trajIter!=trajPointSet.end(); trajIter++ ){
	  COMET::IG4TrajectoryPoint trajPoint = *trajIter;
	  TVector3 tmpPosition = trajPoint.GetPosition().Vect()*(1/unit::cm);
	  //std::cout << GetNode(tmpPosition)->GetName() << std::endl;
	  if (TString(GetNode(tmpPosition)->GetName())=="CDCSenseLayer_0_0") {
	    TVector3 enterPos = trajPoint.GetPosition().Vect()*(1/unit::cm);
	    TVector3 enterMom = trajPoint.GetMomentum()*(1/unit::MeV);
	    fCDCEnterX = enterPos(0);
	    fCDCEnterY = enterPos(1);
	    fCDCEnterZ = enterPos(2);
	    fCDCEnterT = trajPoint.GetPosition()(3)*(1/unit::ns);
	    fCDCEnterPx= enterMom(0);
	    fCDCEnterPy= enterMom(1);
	    fCDCEnterPz= enterMom(2);	    
	    //std::cout << fCDCEnterX << "   " << fCDCEnterY << "   " << fCDCEnterZ << "   " << fCDCEnterT << std::endl;
	    //std::cout << fCDCEnterPx << "   " << fCDCEnterPy << "   " << fCDCEnterPz << "   "<< std::endl;
	      break;
	  }
	}	  		
      }
    }
  }

  std::vector<TString> CDCHitGeometry;
  CDCHitGeometry.push_back("Default");
  int NumOfCDC_0=0;
  int NumOfCDC_1=0;

  if (cdcHits){
    for(COMET::IG4HitContainer::const_iterator hitSeg = cdcHits->begin(); hitSeg != cdcHits->end(); ++hitSeg) {
      COMET::IG4HitSegment* tmpSeg = dynamic_cast<COMET::IG4HitSegment*>(*hitSeg);
            
      if (tmpSeg){

	COMET::IHandle<COMET::IG4Trajectory>  trajectory = TrajCont->GetTrajectory(tmpSeg->GetPrimaryId());
	std::vector <Int_t> trajContributors = tmpSeg->GetContributors();
		
	TVector3 hitPos;
	Double_t hitT;
	
	hitPos(0)=0.5 * (tmpSeg->GetStopX()*(1/unit::cm)+tmpSeg->GetStartX()*(1/unit::cm));
	hitPos(1)=0.5 * (tmpSeg->GetStopY()*(1/unit::cm)+tmpSeg->GetStartY()*(1/unit::cm));
	hitPos(2)=0.5 * (tmpSeg->GetStopZ()*(1/unit::cm)+tmpSeg->GetStartZ()*(1/unit::cm));
	hitT = 0.5 * (tmpSeg->GetStopT()*(1/unit::ns)+tmpSeg->GetStartT()*(1/unit::ns));
	
	TGeoNode* volume = GetNode(hitPos);
	TString geoName = TString(volume->GetName());
	
	if (geoName.Contains("CDCSenseLayer")){
	  fCDCHitX[fnCDCHit]=hitPos(0);
	  fCDCHitY[fnCDCHit]=hitPos(1);
	  fCDCHitZ[fnCDCHit]=hitPos(2);
	  fCDCHitT[fnCDCHit]=hitT;
	  fCDCEDep[fnCDCHit]=tmpSeg->GetEnergyDeposit();
	  	  
	  if (*trajContributors.begin()==1 && trajContributors.size()==1){ // If it is Primary (Signal)
	    if (CDCHitGeometry.back() != geoName){
	      CDCHitGeometry.push_back(geoName);
	      if (geoName=="CDCSenseLayer_0_0"){
		NumOfCDC_0++;
	      }
	      if (geoName=="CDCSenseLayer_1_0"){
		NumOfCDC_1++;
	      }
	    }
	  }	  
	  fTurnId[fnCDCHit]=NumOfCDC_0-1;
	  fnCDCHit++;
	}	  	       	
      }	
    }
  }
  //std::cout << "NumOfCDC_1: " << NumOfCDC_1 << std::endl;
  fTurnNumber=NumOfCDC_1/2;
  

  
  std::map <COMET::IG4HitSegment*, std::vector<COMET::IMCHit*> > HitMap;
  std::vector <COMET::IG4HitSegment*> g4HitOrder;
  std::vector <int> wireIdList;
  
  if (hitHandle){
    COMET::IHitSelection *hits = GetPointer(hitHandle);
    for (COMET::IHitSelection::const_iterator hitSeg = hits->begin(); hitSeg != hits->end(); hitSeg++){
      /// Find SimG4 Hit contributors
      COMET::IMCHit *mcHit = dynamic_cast<COMET::IMCHit*>(GetPointer(*hitSeg));
      std::vector<COMET::IG4VHit*> hitContributors = mcHit->GetContributors();	
      COMET::IG4HitSegment* g4Contributor = dynamic_cast<COMET::IG4HitSegment*>(hitContributors.at(0)); // Currently use only one contributor...
      
      if ( HitMap.find(g4Contributor) == HitMap.end()){      // if NOT exist in map
	std::vector<COMET::IMCHit*> HitVector; HitVector.push_back(mcHit);
	HitMap.insert(std::make_pair(g4Contributor,HitVector));
	g4HitOrder.push_back(g4Contributor);
      }
      
      else if ( HitMap.find(g4Contributor) != HitMap.end()){ // if exist in map
	HitMap[g4Contributor].push_back(mcHit);
	//std::cout << HitMap[g4Contributor].size() << std::endl;
      }
    }
  }
  
  int idOverlap=0;

  for (int i_hit=0 ; i_hit<g4HitOrder.size(); i_hit++){
    
    COMET::IG4HitSegment* g4HitSeg = g4HitOrder.at(i_hit);
    std::vector<COMET::IMCHit*> MCHitVector = HitMap[g4HitSeg];

    TVectorD wireMes(7);	
    
    COMET::IGeometryId geomId = MCHitVector[0]->GetGeomId();
    int wire = COMET::IGeomInfo::Get().CDC().GeomIdToWire(geomId);
    
    if (std::find(wireIdList.begin(), wireIdList.end(), wire) != wireIdList.end() ) {
      idOverlap++;
      continue;      
    }
    wireIdList.push_back(wire);

    fCDCCharge[fnCALCDCHit]    = MCHitVector.size();    
    fWireEnd0X[fnCALCDCHit]    = COMET::IGeomInfo::Get().CDC().GetWireEnd0(wire).X();
    fWireEnd0Y[fnCALCDCHit]    = COMET::IGeomInfo::Get().CDC().GetWireEnd0(wire).Y();
    fWireEnd0Z[fnCALCDCHit]    = COMET::IGeomInfo::Get().CDC().GetWireEnd0(wire).Z();
    fWireEnd1X[fnCALCDCHit]    = COMET::IGeomInfo::Get().CDC().GetWireEnd1(wire).X();
    fWireEnd1Y[fnCALCDCHit]    = COMET::IGeomInfo::Get().CDC().GetWireEnd1(wire).Y();
    fWireEnd1Z[fnCALCDCHit]    = COMET::IGeomInfo::Get().CDC().GetWireEnd1(wire).Z();
    fWireId[fnCALCDCHit]       = wire;	  
    int layer = COMET::IGeomInfo::Get().CDC().GetLayer(wire);

    fWireLayerId[fnCALCDCHit]  = layer;
    if (fMaxWireLayerId<layer) fMaxWireLayerId = layer;
    
    TVector3 g4HitPos;
    g4HitPos(0)=0.5 * (g4HitSeg->GetStopX()*(1/unit::cm)+g4HitSeg->GetStartX()*(1/unit::cm));
    g4HitPos(1)=0.5 * (g4HitSeg->GetStopY()*(1/unit::cm)+g4HitSeg->GetStartY()*(1/unit::cm));
    g4HitPos(2)=0.5 * (g4HitSeg->GetStopZ()*(1/unit::cm)+g4HitSeg->GetStartZ()*(1/unit::cm));
    
    TVector3 local;
    if (!COMET::IGeomInfo::Get().CDC().GetDistanceFromWire(g4HitPos, wire, local)) continue;
    fDriftDist[fnCALCDCHit]=hypot(local.x(),local.y());
        
    //std::cout << fDriftDist[fnCALCDCHit] << std::endl;
    fnCALCDCHit++;
  }

  std::cout << "id Overlap: " << idOverlap << std::endl;
  
  /*
  
  COMET::IChannelId tmpchanId;

  if (hitHandle){
    COMET::IHitSelection *hits = GetPointer(hitHandle);
  
    for (COMET::IHitSelection::const_iterator hitSeg = hits->begin(); hitSeg != hits->end(); hitSeg++){
      
      COMET::IChannelId chanId = (*hitSeg)->GetChannelId();      
      COMET::IGeometryId geomId = (*hitSeg)->GetGeomId();      

      //std::cout << geomId.GetName()  << "   " << geomId.GetSubsystemName() << "   " << geomId.AsInt() << std::endl;
      //std::cout << COMET::GeomId::CDC::IsActiveSenseWire(geomId) << std::endl;
      //std::cout << COMET::GeomId::CDC::GetWireId(geomId,sense,active) << std::endl;

      COMET::IMCHit *mcHit = dynamic_cast<COMET::IMCHit*>(GetPointer(*hitSeg));
      std::vector<COMET::IG4VHit*> hitContributors = mcHit->GetContributors();
      
      int sense=-1; int active=-1;
      int wire = COMET::GeomId::CDC::GetWireId(geomId,sense,active);
      //std::cout << wire << std::endl;
            
      TVectorD wireMes(7);
      wireMes[0] = COMET::IGeomInfo::Get().CDC().GetWireEnd0(wire).X();
      wireMes[1] = COMET::IGeomInfo::Get().CDC().GetWireEnd0(wire).Y();
      wireMes[2] = COMET::IGeomInfo::Get().CDC().GetWireEnd0(wire).Z();
      wireMes[3] = COMET::IGeomInfo::Get().CDC().GetWireEnd1(wire).X();
      wireMes[4] = COMET::IGeomInfo::Get().CDC().GetWireEnd1(wire).Y();
      wireMes[5] = COMET::IGeomInfo::Get().CDC().GetWireEnd1(wire).Z();
      wireMes[6] = (*hitSeg)->GetDriftDistance()*(1/unit::cm);
      
      TVector3 wireend0(wireMes[0],wireMes[1],wireMes[2]);	
      TVector3 wireend1(wireMes[3],wireMes[4],wireMes[5]);

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
      
      //std::cout <<  "X Pos: " << wireMes[0] << std::endl;

      // Merge ionized electrons from same hits
      
      if (hitSeg == hits->begin()){
	tmpchanId = chanId;
	fCDCCharge[fnCALCDCHit]=0; fCDCCharge[fnCALCDCHit]++;
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
	fDriftDist[fnCALCDCHit]+=wireMes[6];
	fCDCCharge[fnCALCDCHit]++;
      }
      
      // When the next hit generates at another Channel Id                                           
      else if (chanId != tmpchanId){

	//std::cout << fWireEnd0X[fnCALCDCHit] << "  " << fWireEnd0Y[fnCALCDCHit] << "   " << fWireEnd0X[fnCALCDCHit] << "   " << fWireEnd1X[fnCALCDCHit] << "   " << fWireEnd1Y[fnCALCDCHit] << "   " << fWireEnd1Z[fnCALCDCHit] << "   " << fCDCCharge[fnCALCDCHit] << "   " << fDriftDist[fnCALCDCHit] << "   ";

	// Average drift distance
	fDriftDist[fnCALCDCHit]=fDriftDist[fnCALCDCHit]/fCDCCharge[fnCALCDCHit];

	// Print to test values
	//	std::cout << fDriftDist[fnCALCDCHit] << std::endl;

	tmpchanId = chanId;
	fnCALCDCHit++;
	fCDCCharge[fnCALCDCHit]=0; fCDCCharge[fnCALCDCHit]++;
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

  */
  
}

void ITracking::ShuffleMCHits(){
  for (int i_swap=0; i_swap<1000; i_swap++){
    int i = rand() % fnCDCHit;
    int j = rand() % fnCDCHit;
    if (i!=j){
      std::swap(fCDCHitX[i], fCDCHitX[j]);
      std::swap(fCDCHitY[i], fCDCHitY[j]);
      std::swap(fCDCHitZ[i], fCDCHitZ[j]);
      std::swap(fCDCHitT[i], fCDCHitY[j]);
      std::swap(fCDCEDep[i], fCDCEDep[j]);
    }
  }
}

void ITracking::PrintMCStatus(){
  std::cout << "----- MC Status -----" << std::endl;
  std::cout << "Number Of Hits: " << fnCALCDCHit << std::endl;
  std::cout << "Initial Energy: " << fGenTrE << std::endl;
  std::cout << std::endl;
}

int ITracking::Finish(){
  return 1;
}


TVector3 GetPOCAofTwoWires(TVector3 wireEnd0_lo, TVector3 wireEnd1_lo, TVector3 wireEnd0_up, TVector3 wireEnd1_up){
  TVector3 u = wireEnd1_lo-wireEnd0_lo;
  TVector3 v = wireEnd1_up-wireEnd0_up;
  TVector3 w = wireEnd0_lo-wireEnd0_up;
  double a,b,c,d,e,s,t;
  a = u * u; 
  b = u * v;
  c = v * v;
  d = u * w;
  e = v * w;
  
  s = (b*e-c*d)/(a*c-b*b);
  t = (a*e-b*d)/(a*c-b*b);

  TVector3 pC=wireEnd0_lo + s * u;
  TVector3 qC=wireEnd0_up + t * v;

  // DOCA 
  //std::cout << "DOCA: " << (w+s*u-t*v).Mag() << std::endl;
  //std::cout << "pC: " << pC(0) << "  " << pC(1) << "  " << pC(2) << std::endl;
  //std::cout << "qC: " << qC(0) << "  " << qC(1) << "  " << qC(2) << std::endl;
  
  return (pC+qC)*0.5;
}

TVector3 GetVectorCrossingCenter(TVector3 wireEnd0_lo, TVector3 wireEnd1_lo, TVector3 wireEnd0_up, TVector3 wireEnd1_up, TVector3 POCA){
  TVector3 u = wireEnd1_lo-wireEnd0_lo;
  TVector3 v = wireEnd1_up-wireEnd0_up;
  
  double u_t   = (POCA(0)-wireEnd0_lo(0))/(wireEnd1_lo(0)-wireEnd0_lo(0));
  double v_t   = (POCA(0)-wireEnd0_up(0))/(wireEnd1_up(0)-wireEnd0_up(0));

  TVector3 c1 = wireEnd0_lo + u_t*u;
  TVector3 c2 = wireEnd0_up + v_t*v;

  //std::cout << "CVector: " << c1(0)-c2(0) << "  " << c1(1)-c2(1) << "  " << c1(2)-c2(2) << std::endl;
  return c2-c1;
}
