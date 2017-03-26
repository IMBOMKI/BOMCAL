#include "HEPUnits.hxx"
#include <IOADatabase.hxx>
#include <cometEventLoop.hxx>
#include <IG4VHit.hxx>
#include <IG4HitSegment.hxx>
#include <IG4Trajectory.hxx>
#include <ICTHChannelId.hxx>
#include <IGeometryId.hxx>
#include <IGeomInfo.hxx>
#include <IDetectorSolenoidGeom.hxx>
#include <memory>
#include <vector>
#include <algorithm>
#include <TTree.h>
#include <TH1F.h>
#include <TPad.h>
#include <TLorentzVector.h>
#include <TGraph.h>
#include "TROOT.h"
#include "TObject.h"
#include "TFile.h"
#include "TApplication.h"
#include "TSystem.h"
#include "TBrowser.h"
#include "TGeoManager.h"
#include "TGeoNode.h"
#include <TGeoNavigator.h>
#include "TMath.h"
#include "TDirectory.h"
#include <TString.h>
#include "TRandom.h"
#include "TRandom3.h"
#include <math.h>
#include <sstream>


////////////////// How to Use this Code //////////////////////////////////
/*

1. After you run the make, 
          you will get a executable file at Linux~ folder

2. By typing following command,
    $./cthHits3.exe -O filename=<name> -O directory=<name> inputfile1.root inputfile2.root ....
    You will get an output files (filename) at the specified path (directory)

*/
//////////////////////////////////////////////////////////////////////////



class SignalTracking: public COMET::ICOMETEventLoopFunction {
public:
  SignalTracking(): fileName("default.root"), fileMode("recreate"), outputDir("../anal") {}           
  virtual ~SignalTracking() {}  
 
  void Usage(void){
    std::cout << "-O filename=<name> Specify the output file name [default=" << fileName << "]" << std::endl;
    std::cout << "-O directory=<name> Specify the output directory [default=" << outputDir << "]" << std::endl;
    std::cout << "-O filemode=<name> Specify the filemode [default=" << fileMode << "]" << std::endl;
    //std::cout << "-O eventNum=<Num> Specify the eventNum [default=" << eventNum << "]" << std::endl;
    
  }

  virtual bool SetOption(std::string option, std::string value=""){
    if(option== "filename") fileName=value;
    else if(option== "directory") outputDir=value;
    else if(option== "filemode") fileMode=value;
    else return false;
    return true;
  }

  TGeoNode* GetNode(TVector3 position) {
    gGeoManager->PushPath();	
    gGeoManager->GetTopNode()->cd();	 	  
    TGeoNode* volume = gGeoManager->FindNode(position(0),position(1),position(2));
    gGeoManager->PopPath();
    return volume;
  }

  Float_t GetEnergy(TVector3 mom){
    Float_t energy = sqrt(pow(mom(0),2)+pow(mom(1),2)+pow(mom(2),2)+pow(0.511,2));
    return energy;
  }

  Float_t GetPhotonNum(Int_t pdg, Float_t Ene, Float_t TrLen){
    Float_t M;
    if (pdg==11 || pdg==-11){
      M=0.511; 
    }
    else if (pdg==13 || pdg ==-13){
      M=105.658;
    }
    else if (pdg==2212 || pdg == -2212){ //proto
      M=938.27;
    }

    else if (pdg==211 || pdg==-211){
      M=139.57;
    }
    else {
      return 0;
    }

    Float_t beta=sqrt(Ene*Ene-M*M)/Ene;
    //std::cout << beta << std::endl;
    Float_t NumOfPhoton = 400 * (1-1/pow(beta*1.5,2)) * TrLen;
    //std::cout << NumOfPhoton << std::endl;
    
    return NumOfPhoton;

  }

  void Initialize() {
        
    strcpy(fullName, outputDir.c_str());
    strcat(fullName, "/");
    strcat(fullName, fileName.c_str());

    outputFile = TFile::Open(fullName, fileMode.c_str());
    outputFile->cd();
    trdata = new TTree("trdata", "Tree Data");
    
    trdata->Branch("nCDCHit", &nCDCHit, "nCDCHit/I");
    trdata->Branch("CDCHitX", CDCHitX, "CDCHitX[nCDCHit]/F");
    trdata->Branch("CDCHitY", CDCHitY, "CDCHitY[nCDCHit]/F");
    trdata->Branch("CDCHitZ", CDCHitZ, "CDCHitZ[nCDCHit]/F");
    trdata->Branch("CDCHitT", CDCHitT, "CDCHitT[nCDCHit]/F");
    trdata->Branch("CDCLayerId", CDCLayerId, "CDCLayerId[nCDCHit]/I");
    trdata->Branch("CDCEDep", CDCEDep, "CDCEDep[nCDCHit]/F");

    trdata->Branch("nCRKHit", &nCRKHit, "nCRKHit/I");
    trdata->Branch("CRKHitX", CRKHitX, "CRKHitX[nCRKHit]/F");
    trdata->Branch("CRKHitY", CRKHitY, "CRKHitY[nCRKHit]/F");
    trdata->Branch("CRKHitZ", CRKHitZ, "CRKHitZ[nCRKHit]/F");
    trdata->Branch("CRKHitT", CRKHitT, "CRKHitT[nCRKHit]/F");
    trdata->Branch("CRKEDep", CRKEDep, "CRKEDep[nCRKHit]/F");

    trdata->Branch("nSTLHit", &nSTLHit, "nSTLHit/I");
    trdata->Branch("STLHitX", STLHitX, "STLHitX[nSTLHit]/F");
    trdata->Branch("STLHitY", STLHitY, "STLHitY[nSTLHit]/F");
    trdata->Branch("STLHitZ", STLHitZ, "STLHitZ[nSTLHit]/F");
    trdata->Branch("STLHitT", STLHitT, "STLHitT[nSTLHit]/F");
    trdata->Branch("STLEDep", STLEDep, "STLEDep[nSTLHit]/F");

    trdata->Branch("nTrigSTL", &nTrigSTL, "nTrigSTL/I");
    trdata->Branch("nTrigCRK", &nTrigCRK, "nTrigCRK/I");    
    trdata->Branch("TrigSTL", TrigSTL, "TrigSTL[nTrigSTL]/I");
    trdata->Branch("TrigCRK", TrigCRK, "TrigCRK[nTrigCRK]/I");
    trdata->Branch("TrigSTLTime", TrigSTLTime, "TrigSTLTime[nTrigSTL]/F");
    trdata->Branch("TrigCRKTime", TrigCRKTime, "TrigCRKTime[nTrigCRK]/F");
    
    //////////////////////////////////////////////////////////////////////////

    trdata->Branch("NumOfTrack", &NumOfTrack, "NumOfTrack/I");
    trdata->Branch("genTrX", genTrX, "genTrX[NumOfTrack]/F");
    trdata->Branch("genTrY", genTrY, "genTrY[NumOfTrack]/F");
    trdata->Branch("genTrZ", genTrZ, "genTrZ[NumOfTrack]/F");
    trdata->Branch("finTrX", finTrX, "finTrX[NumOfTrack]/F");
    trdata->Branch("finTrY", finTrY, "finTrY[NumOfTrack]/F");
    trdata->Branch("finTrZ", finTrZ, "finTrZ[NumOfTrack]/F");
    trdata->Branch("genTrT", genTrT, "genTrT[NumOfTrack]/F");
    trdata->Branch("genTrPx", genTrPx, "genTrPx[NumOfTrack]/F");
    trdata->Branch("genTrPy", genTrPy, "genTrPy[NumOfTrack]/F");
    trdata->Branch("genTrPz", genTrPz, "genTrPz[NumOfTrack]/F");
    trdata->Branch("genTrE", genTrE, "genTrE[NumOfTrack]/F");
    trdata->Branch("PDGNum", PDGNum, "PDGNum[NumOfTrack]/I");
    trdata->Branch("Tr_ifHitCDC", Tr_ifHitCDC, "Tr_ifHitCDC[NumOfTrack]/O");
    trdata->Branch("Tr_ifHitCRK", Tr_ifHitCRK, "Tr_ifHitCRK[NumOfTrack]/O");
    trdata->Branch("Tr_ifHitSTL", Tr_ifHitSTL, "Tr_ifHitSTL[NumOfTrack]/O");
    trdata->Branch("Tr_ifPrimary", Tr_ifPrimary, "Tr_ifPrimary[NumOfTrack]/O");
    trdata->Branch("Tr_ifStartedAtTarget", Tr_ifStartedAtTarget, "Tr_ifStartedAtTarget[NumOfTrack]/O");
    trdata->Branch("Tr_ifStoppedAtTarget", Tr_ifStoppedAtTarget, "Tr_ifStoppedAtTarget[NumOfTrack]/O");
    trdata->Branch("TrackId", TrackId, "TrackId[NumOfTrack]/I");
    trdata->Branch("ParentId", ParentId, "ParentId[NumOfTrack]/I");
    trdata->Branch("eventId", &eventId, "eventId/I");    
  }

  bool operator () (COMET::ICOMETEvent& event) 
  {

    if (event.GetEventId()!=0) {
      return false;}
    // Hit member Nullifying
    
    nCDCHit=0;
    nCRKHit=0;
    nSTLHit=0;

    /*
    ifHitCDC=false;
    ifSingleTurn=false;
    ifMultiTurn=false;
    ifTwoFoldCoin=false;
    ifFourFoldCoin=false;
    if1011FoldCoin=false;
    if1101FoldCoin=false;
    ifHitCDC5thLayer=false;
    */

    memset (CDCHitX, 0, sizeof (CDCHitX));
    memset (CDCHitY, 0, sizeof (CDCHitY));
    memset (CDCHitZ, 0, sizeof (CDCHitZ));
    memset (CDCHitT, 0, sizeof (CDCHitT));
    memset (CDCEDep, 0, sizeof (CDCEDep));
    memset (CDCLayerId, 0, sizeof (CDCLayerId));
    nTrigCRK=0;
    nTrigSTL=0;                                                                            
    memset (TrigCRK, 0, sizeof (TrigCRK));                               
    memset (TrigSTL, 0, sizeof (TrigSTL));         
    memset (TrigCRKTime, 0, sizeof (TrigCRKTime));                               
    memset (TrigSTLTime, 0, sizeof (TrigSTLTime));
    memset (CRKHitX, 0, sizeof (CRKHitX));                                  
    memset (CRKHitY, 0, sizeof (CRKHitY));               
    memset (CRKHitZ, 0, sizeof (CRKHitZ));                    
    memset (CRKHitT, 0, sizeof (CRKHitT));                               
                                              
    memset (STLHitX, 0, sizeof (STLHitX));                
    memset (STLHitY, 0, sizeof (STLHitY));                
    memset (STLHitZ, 0, sizeof (STLHitZ));                
    memset (STLHitT, 0, sizeof (STLHitT));                 
        
    // Track member Nullifying 

    memset (genTrX, 0, sizeof (genTrX));
    memset (genTrY, 0, sizeof (genTrY));
    memset (genTrZ, 0, sizeof (genTrZ));
    memset (genTrT, 0, sizeof (genTrT));
    memset (finTrX, 0, sizeof (finTrX));
    memset (finTrY, 0, sizeof (finTrY));
    memset (finTrZ, 0, sizeof (finTrZ));
    memset (genTrPx, 0, sizeof (genTrPx));
    memset (genTrPy, 0, sizeof (genTrPy));
    memset (genTrPz, 0, sizeof (genTrPz));
    memset (genTrE, 0, sizeof (genTrE));
    memset (Tr_ifHitCDC, 0, sizeof(Tr_ifHitCDC));
    memset (Tr_ifHitCRK, 0, sizeof(Tr_ifHitCRK));
    memset (Tr_ifHitSTL, 0, sizeof(Tr_ifHitSTL));
    memset (Tr_ifPrimary, 0, sizeof(Tr_ifPrimary));
    memset (Tr_ifStartedAtTarget, 0, sizeof(Tr_ifStartedAtTarget));
    memset (Tr_ifStoppedAtTarget, 0, sizeof(Tr_ifStoppedAtTarget));
    memset (TrackId,0,sizeof(TrackId));
    memset (ParentId,0,sizeof(ParentId));
    NumOfTrack=0;
    eventId=NULL;

    /*------------------------------------------------------------------*/    

    eventId = event.GetEventId();  
   
    COMET::IHandle<COMET::IG4TrajectoryContainer> Trajectories = event.Get<COMET::IG4TrajectoryContainer>("truth/G4Trajectories");
    COMET::IG4TrajectoryContainer *TrajCont = GetPointer(Trajectories);
  
    std::vector<TString> CDCHitGeometry;
    CDCHitGeometry.push_back("Default");


    /*----------------- CDC Data ---------------- */


    COMET::IHandle<COMET::IG4HitContainer> CDCHitCont = event.Get<COMET::IG4HitContainer>("truth/g4Hits/CDC");
    if (CDCHitCont){
      for(COMET::IG4HitContainer::const_iterator hitSeg = CDCHitCont->begin(); hitSeg != CDCHitCont->end(); ++hitSeg) {
	COMET::IG4HitSegment* tmpSeg = dynamic_cast<COMET::IG4HitSegment*>(*hitSeg);


        if (tmpSeg){

	  COMET::IHandle<COMET::IG4Trajectory>  trajectory = TrajCont->GetTrajectory(tmpSeg->GetPrimaryId());
	  std::vector <Int_t> trajContributors = tmpSeg->GetContributors();

          TVector3 hitPos;
          Float_t hitT;

          hitPos(0)=0.5 * (tmpSeg->GetStopX()*(1/unit::cm)+tmpSeg->GetStartX()*(1/unit::cm));
          hitPos(1)=0.5 * (tmpSeg->GetStopY()*(1/unit::cm)+tmpSeg->GetStartY()*(1/unit::cm));
          hitPos(2)=0.5 * (tmpSeg->GetStopZ()*(1/unit::cm)+tmpSeg->GetStartZ()*(1/unit::cm));
          hitT = 0.5 * (tmpSeg->GetStopT()*(1/unit::ns)+tmpSeg->GetStartT()*(1/unit::ns));

          gGeoManager=COMET::IOADatabase::Get().Geometry();
          TGeoNode* volume = GetNode(hitPos);
          TString geoName = TString(volume->GetName());

          if (*trajContributors.begin()==1 && trajContributors.size()==1){ // If it is Primary (Signal)                                                                                                                                                                                                                                                                                      
            if (CDCHitGeometry.back() != geoName){
              CDCHitGeometry.push_back(geoName);
            }
          }

          if (geoName.Contains("CDCSenseLayer")){
		CDCHitX[nCDCHit]=hitPos(0);
		CDCHitY[nCDCHit]=hitPos(1);
		CDCHitZ[nCDCHit]=hitPos(2);
		CDCHitT[nCDCHit]=hitT;
		CDCEDep[nCDCHit]=tmpSeg->GetEnergyDeposit();

		Int_t delimFirst= (geoName).First('_');
		Int_t delimLast = (geoName).Last('_');
		Int_t CDCIndex= (((geoName).Remove(delimLast,100)).Remove(0,delimFirst+1)).Atoi();
		CDCLayerId[nCDCHit]=CDCIndex;
		nCDCHit++;
	  }
        }
      }
    }

    /*----------------- CTH Data ---------------- */
    
    COMET::IHandle<COMET::IG4HitContainer> CTHHitCont = event.Get<COMET::IG4HitContainer>("truth/g4Hits/CTH");
    if (CTHHitCont){
      for(COMET::IG4HitContainer::const_iterator hitSeg = CTHHitCont->begin(); hitSeg != CTHHitCont->end(); ++hitSeg) {
	COMET::IG4HitSegment* tmpSeg = dynamic_cast<COMET::IG4HitSegment*>(*hitSeg);

        if (tmpSeg){

	  COMET::IHandle<COMET::IG4Trajectory>  trajectory = TrajCont->GetTrajectory(tmpSeg->GetPrimaryId());
	  std::vector <Int_t> trajContributors = tmpSeg->GetContributors();
	 
          TVector3 hitPos;
          Float_t hitT;

          hitPos(0)=0.5 * (tmpSeg->GetStopX()*(1/unit::cm)+tmpSeg->GetStartX()*(1/unit::cm));
          hitPos(1)=0.5 * (tmpSeg->GetStopY()*(1/unit::cm)+tmpSeg->GetStartY()*(1/unit::cm));
          hitPos(2)=0.5 * (tmpSeg->GetStopZ()*(1/unit::cm)+tmpSeg->GetStartZ()*(1/unit::cm));
          hitT = 0.5 * (tmpSeg->GetStopT()*(1/unit::ns)+tmpSeg->GetStartT()*(1/unit::ns));

          gGeoManager=COMET::IOADatabase::Get().Geometry();
          TGeoNode* volume = GetNode(hitPos);
          TString geoName = TString(volume->GetName());

          if (geoName.Contains("Cherenkov_pv")){
		CRKHitX[nCRKHit]=hitPos(0);
		CRKHitY[nCRKHit]=hitPos(1);
		CRKHitZ[nCRKHit]=hitPos(2);
		CRKHitT[nCRKHit]=hitT;
		CRKEDep[nCRKHit]=tmpSeg->GetEnergyDeposit();

		COMET::IHandle<COMET::IG4Trajectory> tmpHandle = TrajCont->GetTrajectory(trajContributors.back());
		COMET::IG4Trajectory *tmpTraj = GetPointer(tmpHandle);
		Float_t TrLen=tmpSeg->GetTrackLength()*(1/unit::cm);
		Int_t pdg = tmpTraj->GetPDGEncoding();
		Float_t ene = tmpTraj->GetInitialMomentum()(3)*(1/unit::MeV);

		if (GetPhotonNum(pdg,ene,TrLen)>=1 ){ // beta of (electronor positron)is higher than 2/3                                                                                                                                     
                Int_t delim = (geoName).Last('_');
                Int_t CRKIndex= ((geoName).Remove(0,delim+1)).Atoi();
		TrigCRK[nTrigCRK]=CRKIndex;
		TrigCRKTime[nTrigCRK]=hitT;
		nTrigCRK++;
		
		}

		nCRKHit++;
	  }
	  
	  else if (geoName.Contains("Scintillator_pv")){
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
	  }
	}
      }   
    }

    /*--------------------------- Trajectory ------------------------


      This Trajectory analysis is not perfect. (and It should not be)
      It is just used for conservative study to investigate Bunch Beam's acceptance


     ----------------------------------------------------------------*/

    if(TrajCont->empty()){                                                        
      std::cout<< "Result is not found" << std::endl;                               
    } 
                                                                                       
    if(!TrajCont->empty()){
      

      for(COMET::IG4TrajectoryContainer::const_iterator seg = TrajCont->begin(); seg != TrajCont->end(); seg++){
	
	COMET::IG4Trajectory traj = (*seg).second;     
	
	TVector3 iniPos = traj.GetInitialPosition().Vect()*(1/unit::cm);
	TVector3 finPos = traj.GetFinalPosition().Vect()*(1/unit::cm);
	TVector3 iniMom = traj.GetInitialMomentum().Vect()*(1/unit::MeV);

	TGeoNode* iniVolume=GetNode(iniPos);
	TGeoNode* finVolume=GetNode(finPos);
	
	if (TString(iniVolume->GetName()).Contains("Target")){
	  Tr_ifStartedAtTarget[NumOfTrack]=true;
	}

	if (TString(finVolume->GetName()).Contains("Target")){
	  Tr_ifStoppedAtTarget[NumOfTrack]=true;
	}

        if(!COMET::IGeomInfo::DetectorSolenoid().GetDetPositionInDSCoordinate(iniPos, iniPos)){
          continue;
	  std::cout << "Coordinate change does not work well" << std::endl;}

        if(!COMET::IGeomInfo::DetectorSolenoid().GetDetPositionInDSCoordinate(finPos, finPos)){
          continue;
	  std::cout << "Coordinate change does not work well" << std::endl;}
		
	genTrX[NumOfTrack] = iniPos(0);
	genTrY[NumOfTrack] = iniPos(1);
	genTrZ[NumOfTrack] = iniPos(2);
	genTrT[NumOfTrack] = traj.GetInitialPosition()(3);
	finTrX[NumOfTrack] = finPos(0);
	finTrY[NumOfTrack] = finPos(1);
	finTrZ[NumOfTrack] = finPos(2);
	genTrPx[NumOfTrack] = iniMom(0);
	genTrPy[NumOfTrack] = iniMom(1);
	genTrPz[NumOfTrack] = iniMom(2);	
	genTrE[NumOfTrack] = traj.GetInitialMomentum()(3);
	PDGNum[NumOfTrack] = traj.GetPDGEncoding();
	TrackId[NumOfTrack] = traj.GetTrackId();
	ParentId[NumOfTrack] = traj.GetParentId();
	if (traj.GetParentId()==0) Tr_ifPrimary[NumOfTrack]=1;

	/*
	if (PDGNum[NumOfTrack]==22){
	  COMET::IHandle<COMET::IG4Trajectory> parentTraj_Handle = TrajCont->GetTrajectory(traj.GetParentId());  
	  COMET::IG4Trajectory* parentTraj = GetPointer(parentTraj_Handle);
	  TVector3 parent_iniPos = parentTraj->GetInitialPosition().Vect()*(1/unit::cm);
	  TVector3 parent_finPos = parentTraj->GetFinalPosition().Vect()*(1/unit::cm);
	  TGeoNode* parent_iniVolume=GetNode(parent_iniPos);
	  TGeoNode* parent_finVolume=GetNode(parent_finPos);
	  
	  
	  if (Tr_ifPrimary[NumOfTrack]==false && parentTraj->GetParticleName()=="pi-"){
	  std::cout << parentTraj->GetPDGEncoding() << "   " 
		    << parentTraj->GetParticleName() << "   "
		    << TString(parent_iniVolume->GetName()) << "   "
		    << TString(parent_finVolume->GetName()) << "   " 
		    << parentTraj->GetTrackId() << "   "
		    << std::endl;
	  }	  
	}
	*/	

	COMET::IG4Trajectory::Points trajPointSet = traj.GetTrajectoryPoints();
	
	for(COMET::IG4Trajectory::Points::iterator trajIter = trajPointSet.begin(); trajIter!=trajPointSet.end(); trajIter++ ){
	  COMET::IG4TrajectoryPoint trajPoint = *trajIter;

          TVector3 TrPos = trajPoint.GetPosition().Vect()*(1/unit::cm);
          Float_t TrTime = trajPoint.GetPosition()(3)*(1/unit::ns);
	 
          gGeoManager=COMET::IOADatabase::Get().Geometry();
          TGeoNode* volume = GetNode(TrPos);
	  TString geoName = volume->GetName();

	  if (geoName.Contains("CDCSenseLayer")){
	    Tr_ifHitCDC[NumOfTrack]=true;
	  }

	  else if (geoName.Contains("Cherenkov_pv")){
	    Tr_ifHitCRK[NumOfTrack]=true;
	  }

	  else if (geoName.Contains("Scintillator_pv")){
	    Tr_ifHitSTL[NumOfTrack]=true;
	  }
	  
	}
	
	NumOfTrack++;
      
      } // Trajectory for loop close      
      
    } // if(!Trajectories.empty()) close
    
    trdata->Fill();
    return true;
  }
  
  void Finalize(COMET::ICOMETOutput * const output) {    
    outputFile->cd();    
    trdata->Write();
    outputFile->Close();
  }

private:

  /*-------------- File Variables --------------*/

  TFile* outputFile;
  TTree* trdata;
  TBranch* CDCHit;
  TBranch* evtId;

  std::string outputDir;
  std::string fileName;
  std::string fileMode;
  char fullName[100];

  /*------------- Hit Class Root Member -----------------*/
  
  Int_t nCDCHit;
  Float_t CDCHitX[200000];
  Float_t CDCHitY[200000];
  Float_t CDCHitZ[200000];
  Float_t CDCHitT[200000];
  Float_t CDCEDep[200000];
  Int_t CDCLayerId[200000];

  Int_t nCRKHit;
  Float_t CRKHitX[30000];
  Float_t CRKHitY[30000];
  Float_t CRKHitZ[30000];
  Float_t CRKHitT[30000];
  Float_t CRKEDep[30000];

  Int_t nSTLHit;
  Float_t STLHitX[30000];
  Float_t STLHitY[30000];
  Float_t STLHitZ[30000];
  Float_t STLHitT[30000];
  Float_t STLEDep[30000];

  Int_t nTrigSTL;
  Int_t nTrigCRK;
  Int_t TrigSTL[10000];
  Int_t TrigCRK[10000];
  Float_t TrigSTLTime[10000];
  Float_t TrigCRKTime[10000];
  
  /*------------- Trajectory Class Root Member -----------------*/

  Float_t genTrX[40000];
  Float_t genTrY[40000];
  Float_t genTrZ[40000];
  Float_t finTrX[40000];
  Float_t finTrY[40000];
  Float_t finTrZ[40000];
  Float_t genTrT[40000];
  Float_t genTrPx[40000];
  Float_t genTrPy[40000];
  Float_t genTrPz[40000];
  Float_t genTrE[40000];
  Int_t PDGNum[40000];
  Int_t NumOfTrack;
  Bool_t Tr_ifHitCDC[40000];
  Bool_t Tr_ifHitCRK[40000];
  Bool_t Tr_ifHitSTL[40000];
  Bool_t Tr_ifPrimary[40000];
  Bool_t Tr_ifStartedAtTarget[40000];
  Bool_t Tr_ifStoppedAtTarget[40000];
  Int_t TrackId[40000];
  Int_t ParentId[40000];
  
  TGeoManager* gGeoManager;

  /*-------------------------------------------*/

  Int_t eventId;

};

int main(int argc, char **argv) {
  SignalTracking userCode;
  cometEventLoop(argc,argv,userCode);
}
