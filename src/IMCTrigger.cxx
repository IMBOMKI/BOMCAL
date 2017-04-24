#include "HEPUnits.hxx"
#include <vector>
#include <iostream>

#include <IOADatabase.hxx>
#include <IMCTrigger.hxx>
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

IMCTrigger::IMCTrigger(const char* name = "IMCTrigger", const char* title="mctrigger"){;}
IMCTrigger::~IMCTrigger(){;}

int IMCTrigger::Init()
{
  COMETNamedInfo("IMCTrigger", "//----------------------------------------------------------//");
  COMETNamedInfo("IMCTrigger", "//   Initialize IMCTrigger                               ");
  COMETNamedInfo("IMCTrigger", "//----------------------------------------------------------//");
  return 1;
}

TGeoNode* IMCTrigger::GetNode(TVector3 position){
    gGeoManager->PushPath();	
    gGeoManager->GetTopNode()->cd();	 	  
    TGeoNode* volume = gGeoManager->FindNode(position(0),position(1),position(2));
    gGeoManager->PopPath();
    return volume;
}

int IMCTrigger::GetPhotonNumber(int pdgNum, double ene, double trLen){
    Double_t mass;
    if (pdgNum==11 || pdgNum==-11) mass=0.511;
    else if (pdgNum==13 || pdgNum ==-13) mass=105.658;
    else if (pdgNum==2212 || pdgNum == -2212) mass=938.27;
    else if (pdgNum==211 || pdgNum==-211) mass=139.57;
    else return 0;
   
    Double_t beta=sqrt(ene*ene-mass*mass)/ene;
    Double_t NumOfPhoton = 400 * (1-1/pow(beta*1.5,2)) * trLen;
    return NumOfPhoton;
}

int IMCTrigger::Finish(){
  return 1;
}

COMET::IG4HitContainer* IMCTrigger::MakeCTHHitSelection(COMET::IHandle<COMET::IG4HitContainer> & cthhits){
  COMETNamedDebug("IMCTrigger", "Start MakeCTHHitSelection");
  COMET::IG4HitContainer* hits = new COMET::IG4HitContainer("mchits", "Hits found by IMCTrigger");
  
  if (cthhits){
    for (COMET::IG4HitContainer::const_iterator hitSeg = cthhits->begin(); hitSeg != cthhits->end(); ++ hitSeg){

      COMET::IG4HitSegment* tmpSeg = dynamic_cast<COMET::IG4HitSegment*>(*hitSeg);
      if (tmpSeg){

	TVector3 hitPos;
	Double_t hitT;

	hitPos(0)=0.5 * (tmpSeg->GetStopX()*(1/unit::cm)+tmpSeg->GetStartX()*(1/unit::cm));
	hitPos(1)=0.5 * (tmpSeg->GetStopY()*(1/unit::cm)+tmpSeg->GetStartY()*(1/unit::cm));
	hitPos(2)=0.5 * (tmpSeg->GetStopZ()*(1/unit::cm)+tmpSeg->GetStartZ()*(1/unit::cm));
	hitT = 0.5 * (tmpSeg->GetStopT()*(1/unit::ns)+tmpSeg->GetStartT()*(1/unit::ns));
	TGeoNode* volume = GetNode(hitPos);
	TString geoName = TString(volume->GetName());

	if (geoName.Contains("Cherenkov_pv"))
	  {
	    /*
	    CRKHitX[nCRKHit]=hitPos(0);
	    CRKHitY[nCRKHit]=hitPos(1);
	    CRKHitZ[nCRKHit]=hitPos(2);
	    CRKHitT[nCRKHit]=hitT;
	    CRKEDep[nCRKHit]=tmpSeg->GetEnergyDeposit();
	    
	    COMET::IHandle<COMET::IG4Trajectory> tmpHandle = TrajCont->GetTrajectory(trajContributors.back());
	    COMET::IG4Trajectory *tmpTraj = GetPointer(tmpHandle);
	    Double_t TrLen=tmpSeg->GetTrackLength()*(1/unit::cm);
	    Int_t pdg = tmpTraj->GetPDGEncoding();
	    Double_t ene = tmpTraj->GetInitialMomentum()(3)*(1/unit::MeV);
	    
	    if (GetPhotonNum(pdg,ene,TrLen)>=1 ){
	      
	      Int_t delim = (geoName).Last('_');
	      Int_t CRKIndex= ((geoName).Remove(0,delim+1)).Atoi();
	      
	      TrigCRK[nTrigCRK]=CRKIndex;
	      TrigCRKTime[nTrigCRK]=hitT;
	      nTrigCRK++;
	      
	    }
	    
	    nCRKHit++;	    
	    */
	  }
	
	
	else if (geoName.Contains("Scintillator_pv"))
	  {
	    /*
	    STLHitX[nSTLHit]=hitPos(0);
	    STLHitY[nSTLHit]=hitPos(1);
	    STLHitZ[nSTLHit]=hitPos(2);
	    STLHitT[nSTLHit]=hitT;
	    STLEDep[nSTLHit]=tmpSeg->GetEnergyDeposit();	      
	    
	    if (STLEDep[nSTLHit]*(1/unit::MeV)>0.063){   // Energy deposit is higher than 63keV                        
	      
	      Int_t delim = (geoName).Last('_');
	      Int_t STLIndex= ((geoName).Remove(0,delim+1)).Atoi();
	      
	      TrigSTL[nTrigSTL]=STLIndex;
	      TrigSTLTime[nTrigSTL]=hitT;
	      nTrigSTL++;
	      
	    }
	    
	    nSTLHit++;
	    */
	  }
	
      }  
    }
  }
  
}
