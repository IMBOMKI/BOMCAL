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
#include <IGeomInfo.hxx>
#include <IGeometryId.hxx>
#include <COMETGeomId.hxx>
#include <IGeomIdManager.hxx>
#include <ICTHGeomId.hxx>
#include <ICTHChannelId.hxx>
#include <IChannelId.hxx>
#include <IGeometryDatabase.hxx>
#include <ICDCGeom.hxx>

IMCTrigger::IMCTrigger(const char* name = "IMCTrigger", const char* title="mctrigger")
:fStartT(700), fEndT(1170)
{;}
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

void IMCTrigger::MakeCTHMap(COMET::IHandle<COMET::IG4HitContainer> & cthhits, COMET::IHandle<COMET::IG4TrajectoryContainer> trajectories)
{
  COMET::IG4TrajectoryContainer *TrajCont = GetPointer(trajectories);
  COMETNamedDebug("IMCTrigger", "Start MakeCTHHitSelection");

  if (cthhits){
    int count=0;
    for (COMET::IG4HitContainer::const_iterator hitSeg = cthhits->begin(); hitSeg != cthhits->end(); ++ hitSeg){
      COMET::IG4HitSegment* tmpSeg = dynamic_cast<COMET::IG4HitSegment*>(*hitSeg);
      if (tmpSeg){
	double x =0.5 * (tmpSeg->GetStopX()*(1/unit::cm)+tmpSeg->GetStartX()*(1/unit::cm));
	double y =0.5 * (tmpSeg->GetStopY()*(1/unit::cm)+tmpSeg->GetStartY()*(1/unit::cm));
	double z =0.5 * (tmpSeg->GetStopZ()*(1/unit::cm)+tmpSeg->GetStartZ()*(1/unit::cm));
	double t = 0.5 * (tmpSeg->GetStopT()*(1/unit::ns)+tmpSeg->GetStartT()*(1/unit::ns));

	COMET::IGeometryId CTHId;
	COMET::IOADatabase::Get().GeomId().GetGeometryId(x, y, z, CTHId);
	COMET::ICTHChannelId CTHChannelId;
	COMET::IGeometryDatabase::Get().GetChannelId(CTHChannelId, CTHId);

        UInt_t scint   = CTHChannelId.GetScint();        // Scintillator(1) or Cherenkov(0)     
	UInt_t lightG  = CTHChannelId.GetLightGuide();   // Light Guide             
        UInt_t index   = CTHChannelId.GetCounter();      // segment ID (1-64)         
        UInt_t module  = CTHChannelId.GetModule();       // Upstream (0), Downstream(1)  

	if (scint==0 && lightG==0) //Cherenkov
	  {
	    COMET::IHandle<COMET::IG4Trajectory> tmpHandle = TrajCont->GetTrajectory(tmpSeg->GetContributors().back());
	    COMET::IG4Trajectory *tmpTraj = GetPointer(tmpHandle);
	    double TrLen=tmpSeg->GetTrackLength()*(1/unit::cm);
	    int pdg = tmpTraj->GetPDGEncoding();
            double ene = tmpTraj->GetInitialMomentum()(3)*(1/unit::MeV);

	    if (GetPhotonNumber(pdg,ene,TrLen)>=1 ){	      
	      fTime.insert (std::make_pair(count, t));
	      fIndex.insert(std::make_pair(count, index));
	      fScint.insert(std::make_pair(count, scint));
	      fModule.insert(std::make_pair(count, module));
	      count++;
	    }
	  }
		
	else if (scint==1 && lightG==0)  //Scintillator
	  {	    
	    fTime.insert(std::make_pair(count, t));
	    fIndex.insert(std::make_pair(count, index));
	    fScint.insert(std::make_pair(count, scint));
	    fModule.insert(std::make_pair(count, module));
	    count++;	   
	  }	
      }  
    }
  }  
}

bool TimeComparison(const std::pair<int,double> &A,const std::pair<int,double> &B)
{
  return A.second<B.second;
}

void IMCTrigger::MakeTimeCluster(int Module){
  std::vector <std::pair <int, double> > TimeOrdering;  // (int->key Value, double->time)
  std::map< int, int >::iterator iter;
  int cluCount=0;

  for( iter = fModule.begin(); iter != fModule.end(); ++iter)
    {
      double hitTime = fTime.at(iter->first);
      if (iter->second == Module){      // Down stream
	
	TimeOrdering.push_back(std::make_pair(iter->first,hitTime));
	sort(TimeOrdering.begin(), TimeOrdering.end(), TimeComparison);       
      }      
    }
  
  for (int i=0 ; i<TimeOrdering.size(); i++){
    if ((TimeOrdering.at(i+1).second - TimeOrdering.at(i).second < 10 ) && (i+1<TimeOrdering.size())){
      double refT = TimeOrdering.at(i).second;
      std::vector < int > ClusterKeyVal;
      std::vector < int > ScintClusterKeyVal;
      std::vector < int > CherenClusterKeyVal;

      ClusterKeyVal.push_back(TimeOrdering.at(i).first);
      int j=i+1;
      while ( TimeOrdering.at(j).second-refT < 10 && j<TimeOrdering.size()){	
	bool Push=1;
	for (int i_clu=0; i_clu<ClusterKeyVal.size(); i_clu++){
	  int keyVal  = TimeOrdering.at(j).first;
	  int keyVal2 = ClusterKeyVal.at(i_clu);	  
	  if (fScint.at(keyVal)==fScint.at(keyVal2) && fIndex.at(keyVal)==fIndex.at(keyVal2)){
	    Push=0;
	  }
	}
	if (Push==1) ClusterKeyVal.push_back(TimeOrdering.at(j).first);
	j++;	
      }

      // Sort with Scintillator and Cherenkov
      for (int ii=0; ii < ClusterKeyVal.size(); ii++){      
	int keyVal = ClusterKeyVal.at(ii);
	if (fScint.at(keyVal)==1){
	  ScintClusterKeyVal.push_back(keyVal);
	}
	else if (fScint.at(keyVal)==0){
	  CherenClusterKeyVal.push_back(keyVal);
	}
      }

      if (ScintClusterKeyVal.size()>0 && CherenClusterKeyVal.size()>0){
	if (Module==1) fDSTimeCluster.push_back(std::make_pair(ScintClusterKeyVal,CherenClusterKeyVal));
        else if (Module==1) fUSTimeCluster.push_back(std::make_pair(ScintClusterKeyVal,CherenClusterKeyVal));
      }

      i=j;
      cluCount++;
    }
  }
}

void IMCTrigger::PrintTimeCluster(){
  for (int i=0 ; i<fDSTimeCluster.size(); i++){
    
    

  }

  for (int i=0 ; i<fUSTimeCluster.size(); i++){

  }

}

void IMCTrigger::Clear(){
  fTime.clear();  
  fIndex.clear();
  fScint.clear();
  fModule.clear();
  fDSTimeCluster.clear();
  fUSTimeCluster.clear();
}

int IMCTrigger::Finish(){
  return 1;
}
