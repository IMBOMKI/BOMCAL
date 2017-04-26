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
:fStartT(700), fEndT(1170), fCTHSegNum(64)
{;}
IMCTrigger::~IMCTrigger(){;}

int IMCTrigger::Init()
{
  COMETNamedInfo("IMCTrigger", "//----------------------------------------------------------//");
  COMETNamedInfo("IMCTrigger", "//   Initialize IMCTrigger                               ");
  COMETNamedInfo("IMCTrigger", "//----------------------------------------------------------//");
  return 1;
}

bool DoubleComparison(const std::pair<int,double> &A,const std::pair<int,double> &B)
{
  return A.second<B.second;
}

bool IntComparison(const std::pair<int,int> &A,const std::pair<int,int> &B)
{
  return A.second<B.second;
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

void IMCTrigger::MakeTimeCluster(int Module){
  std::vector <std::pair <int, double> > TimeOrdering;  // (int->key Value, double->time)
  std::map< int, int >::iterator iter;

  if (fModule.size()>0){
    for( iter = fModule.begin(); iter != fModule.end(); ++iter)
      {
	double hitTime = fTime.at(iter->first);
	if (iter->second == Module){      // Down stream
	  
	  TimeOrdering.push_back(std::make_pair(iter->first,hitTime));
	  sort(TimeOrdering.begin(), TimeOrdering.end(), DoubleComparison);       
	}      
      }
  }  

  

  if (TimeOrdering.size()>1){
    for (int i=0 ; i<TimeOrdering.size(); i++){
      if ((i+1<TimeOrdering.size()) && (TimeOrdering.at(i+1).second - TimeOrdering.at(i).second < 10 )){	
	double refT = TimeOrdering.at(i).second;
	std::vector < int > ClusterKeyVal;
	std::vector < int > ScintClusterKeyVal;
	std::vector < int > CherenClusterKeyVal;
	
	ClusterKeyVal.push_back(TimeOrdering.at(i).first);
	int j=i+1;
	while ( j<TimeOrdering.size() && TimeOrdering.at(j).second-refT < 10){	
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
	  if (j == TimeOrdering.size()) break;
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
	  else if (Module==0) fUSTimeCluster.push_back(std::make_pair(ScintClusterKeyVal,CherenClusterKeyVal));
	}
	
	i=j;	
      }
    }
  }
}

void IMCTrigger::PrintTimeCluster(){
  
  if (fDSTimeCluster.size() > 0){

    std::cout << "----------------------Down Stream--------------------" << std::endl;

    for (int i=0 ; i<fDSTimeCluster.size(); i++){
      std::vector< int >  ScintCluster = fDSTimeCluster.at(i).first;
      std::vector< int >  CherenCluster = fDSTimeCluster.at(i).second;
      
      if (ScintCluster.size()>0){
	for (int i_sci = 0; i_sci<ScintCluster.size(); i_sci++){
	  std::cout << fIndex.at(ScintCluster.at(i_sci)) << "   ";
	}      
      }
      
      std::cout << "|   ";
      
      if (CherenCluster.size()>0){
	for (int i_che = 0; i_che<CherenCluster.size(); i_che++){
	  std::cout << fIndex.at(CherenCluster.at(i_che)) << "   ";
	}      
      }

    }

    std::cout << std::endl;
  }
  
  if (fUSTimeCluster.size() > 0){

    std::cout << "----------------------Up Stream----------------------" << std::endl;

    for (int i=0 ; i<fUSTimeCluster.size(); i++){
      std::vector< int >  ScintCluster = fUSTimeCluster.at(i).first;
      std::vector< int >  CherenCluster = fUSTimeCluster.at(i).second;
      
      if (ScintCluster.size()>0){
	for (int i_sci = 0; i_sci<ScintCluster.size(); i_sci++){
	  std::cout << fIndex.at(ScintCluster.at(i_sci)) << "   ";
	}      
      }
      
      std::cout << "|   ";
      
      if (CherenCluster.size()>0){
	for (int i_che = 0; i_che<CherenCluster.size(); i_che++){
	  std::cout << fIndex.at(CherenCluster.at(i_che)) << "   ";
	}      
      }	
      
    }   
    std::cout << std::endl;
  }
}

void IMCTrigger::ApplyShiftCondition(int Module, int shift){
  
  std::vector < std::pair< std::vector< int >, std::vector< int > > > TimeCluster;
  if (Module==1) TimeCluster=fDSTimeCluster;
  if (Module==0) TimeCluster=fUSTimeCluster;

  if (TimeCluster.size()>0){
    
    for (int i=0 ; i<TimeCluster.size() ; i++){
      
      std::vector< std::vector <int> > ScintDistCluster;  // key Value
      std::vector< std::vector <int> > CherenDistCluster; // key Value

      std::pair< std::vector< int >, std::vector< int> > Cluster = TimeCluster.at(i);
      std::vector< int > ScintKey  = Cluster.first;
      std::vector< int > CherenKey = Cluster.second;

      // Make DistanceCluster for Scintillator

      if (ScintKey.size()>1){
	std::vector< std::pair<int, int> > ScintIndex;

	for (int i_sci=0; i_sci<ScintKey.size(); i_sci++){
	  int keyVal = ScintKey.at(i_sci);
	  int indexVal = fIndex.at(keyVal);
	  ScintIndex.push_back(std::make_pair(keyVal,indexVal));
	}
	sort(ScintIndex.begin(), ScintIndex.end(), IntComparison);             
		
	for (int i_sci=0; i_sci<ScintIndex.size(); i_sci++){
	  
	  if (i_sci+1<ScintIndex.size() && (ScintIndex.at(i_sci+1).second - ScintIndex.at(i_sci).second)%fCTHSegNum <= 1){
	    
	    std::vector< int > SingleCluster;
	    SingleCluster.push_back(ScintIndex.at(i_sci).first);
	    
	    while (i_sci+1<ScintIndex.size() && (ScintIndex.at(i_sci+1).second - ScintIndex.at(i_sci).second)%fCTHSegNum <= 1){
	      SingleCluster.push_back(ScintIndex.at(i_sci+1).first);
	      i_sci++;	      	      
	      if (i_sci+1==ScintIndex.size()) break;
	    }	    
	    ScintDistCluster.push_back(SingleCluster);	    
	  }       	  
	}       	
      }
      
      // Make DistanceCluster for Cherenkov
      
      if (CherenKey.size()>1){
	std::vector< std::pair<int, int> > CherenIndex;

	for (int i_che=0; i_che<CherenKey.size(); i_che++){
	  int keyVal = CherenKey.at(i_che);
	  int indexVal = fIndex.at(keyVal);
	  CherenIndex.push_back(std::make_pair(keyVal,indexVal));
	}
	sort(CherenIndex.begin(), CherenIndex.end(), IntComparison);             

	for (int i_che=0; i_che<CherenIndex.size(); i_che++){
	  if (i_che+1<CherenIndex.size() && (CherenIndex.at(i_che+1).second - CherenIndex.at(i_che).second)%fCTHSegNum <= 1){
	    std::vector< int > SingleCluster;
	    SingleCluster.push_back(CherenIndex.at(i_che).first);
	    
	    while (i_che+1<CherenIndex.size() && (CherenIndex.at(i_che+1).second - CherenIndex.at(i_che).second)%fCTHSegNum <= 1){
	      SingleCluster.push_back(CherenIndex.at(i_che+1).first);
	      i_che++;	      	      
	      if (i_che+1==CherenIndex.size()) break;
	    }	    
	    CherenDistCluster.push_back(SingleCluster);
	  }       
	}
      }
      
      ///////////////////////////
      // Apply Shift Condition //
      ///////////////////////////
      
      if (ScintDistCluster.size()>0 && CherenDistCluster.size()>0){
	for (int i_che=0; i_che<CherenDistCluster.size(); i_che++){
	  for (int i_sci=0; i_sci<ScintDistCluster.size(); i_sci++){

	    std::vector < int > CherenCluster_key = CherenDistCluster.at(i_che); // key value
	    std::vector < int > ScintCluster_key = ScintDistCluster.at(i_sci);   // key value
	    
	    std::vector < int >  ScintCluster;   // index
	    std::vector < int > CherenCluster;   // index
	    
	    for (int i_ele=0; i_ele<ScintCluster_key.size(); i_ele++){
	      ScintCluster.push_back(fIndex.at(ScintCluster_key.at(i_ele)));
	    }
	    for (int i_ele=0; i_ele<CherenCluster_key.size(); i_ele++){
	      CherenCluster.push_back(fIndex.at(CherenCluster_key.at(i_ele)));
	    }
	    
	    
	    for (int i_shift=-shift ; i_shift <= shift; i_shift++){
	      int NumOfOverlap=0;
	      for (int i_clu=0; i_clu<ScintCluster.size(); i_clu++){
		
		if ( std::find(CherenCluster.begin(),
			       CherenCluster.end(),
      			       (ScintCluster.at(i_clu)+i_shift)%fCTHSegNum) 
		     != CherenCluster.end() ){
		  NumOfOverlap++;
		}		
	      }
	      if (NumOfOverlap>=2) {
		fPairCandidates.push_back(make_pair(ScintCluster_key,CherenCluster_key));
		break;
	      }
	    }
	  }
	}
      }           
    }
  }
}

void IMCTrigger::PrintPairCandidates(){
  if (fPairCandidates.size()>0){    
    std::cout << "----------------- Pair Candidates -----------------" << std::endl;
    for (int i_cand=0; i_cand<fPairCandidates.size(); i_cand++){
      std::vector< int > ScintCandidate = fPairCandidates.at(i_cand).first;
      std::vector< int > CherenCandidate = fPairCandidates.at(i_cand).second;
      
      for (int i_sci=0; i_sci<ScintCandidate.size(); i_sci++){
	int index;      
	int keyVal=ScintCandidate.at(i_sci);        
	if (fModule.at(keyVal)==0) 
	  {index=fIndex.at(keyVal)+fCTHSegNum;}       
	else if (fModule.at(keyVal)==1) 
	  {index=fIndex.at(keyVal);}	
	std::cout << index << "  ";// << fScint.at(keyVal) << " " << fModule.at(keyVal) << " " << fTime.at(keyVal) << "   ";
      }
      std::cout << "|  ";
      
      for (int i_che=0; i_che<CherenCandidate.size(); i_che++){
	int index;      
	int keyVal=CherenCandidate.at(i_che);        
	if (fModule.at(keyVal)==0) {index=fIndex.at(keyVal)+fCTHSegNum;}
	else {index=fIndex.at(keyVal);}
	std::cout << index << "  ";//<< fScint.at(keyVal) << " " << fModule.at(keyVal) << " " << fTime.at(keyVal) << "   ";
      }      
      std::cout << std::endl;      
    }
  }
}

void IMCTrigger::Process(int shift){
  MakeTimeCluster(1);
  MakeTimeCluster(0);
  ApplyShiftCondition(1,shift);
  ApplyShiftCondition(0,shift);  
}

bool IMCTrigger::GetFourFoldCoincidence(){
  if (fPairCandidates.size()>0) return 1;
  return 0;
}

void IMCTrigger::Clear(){
  fTime.clear();  
  fIndex.clear();
  fScint.clear();
  fModule.clear();
  fDSTimeCluster.clear();
  fUSTimeCluster.clear();
  fPairCandidates.clear();
}

int IMCTrigger::Finish(){
  return 1;
}
