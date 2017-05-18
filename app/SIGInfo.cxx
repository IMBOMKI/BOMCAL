#include "HEPUnits.hxx"
#include <IOADatabase.hxx>
#include <cometEventLoop.hxx>
#include <IG4VHit.hxx>
#include <IG4HitSegment.hxx>
#include <IMCHit.hxx>
#include <IHit.hxx>
#include <IG4Trajectory.hxx>

#include <IGeometryId.hxx>
#include <IGeometryDatabase.hxx>
#include <IGeomInfo.hxx>
#include <IGeomIdManager.hxx>
#include <IChannelId.hxx>
#include <ICTHGeomId.hxx>
#include <ICTHChannelId.hxx>
#include <IDetectorSolenoidGeom.hxx>
#include <ICDCmcDigit.hxx>

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
#include <map>

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
    std::cout << "-O directory=<name> Specify the output directory [default=" << fileName << "]" << std::endl;
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

  Double_t GetEnergy(TVector3 mom){
    Double_t energy = sqrt(pow(mom(0),2)+pow(mom(1),2)+pow(mom(2),2)+pow(0.511,2));
    return energy;
  }

  Double_t GetPhotonNum(Int_t pdg, Double_t Ene, Double_t TrLen){
    Double_t M;
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

    Double_t beta=sqrt(Ene*Ene-M*M)/Ene;
    Double_t NumOfPhoton = 400 * (1-1/pow(beta*1.5,2)) * TrLen;

    return NumOfPhoton;

  }

  void Initialize() {
    
    strcpy(fullName, outputDir.c_str());
    strcat(fullName, "/");
    strcat(fullName, fileName.c_str());

    outputFile = TFile::Open(fullName, fileMode.c_str());
    outputFile->cd();
    trdata = new TTree("trdata", "Tree Data");
        
    trdata->Branch("eventId", &eventId, "eventId/I");

    trdata->Branch("nCDCHit", &nCDCHit, "nCDCHit/I");
    trdata->Branch("nCDCHit_Pairep", &nCDCHit_Pairep, "nCDCHit_Pairep/I");
    trdata->Branch("nCDCHit_Pairem", &nCDCHit_Pairem, "nCDCHit_Pairem/I");
    trdata->Branch("nCDCHit_Comptonem", &nCDCHit_Comptonem, "nCDCHit_Comptonem/I");
    trdata->Branch("ifCDC_Pairep", ifCDC_Pairep, "ifCDC_Pairep[nCDCHit]/O");
    trdata->Branch("ifCDC_Pairem", ifCDC_Pairem, "ifCDC_Pairem[nCDCHit]/O");
    trdata->Branch("ifCDC_Comptonem", ifCDC_Comptonem, "ifCDC_Comptonem[nCDCHit]/O");
    trdata->Branch("CDCLayerId_Pairep", CDCLayerId_Pairep, "CDCLayerId_Pairep[nCDCHit_Pairep]/I");
    trdata->Branch("MaxCDCLayerId_Pairep", &MaxCDCLayerId_Pairep, "MaxCDCLayerId_Pairep/I");
    trdata->Branch("CDCLayerId_Pairem", CDCLayerId_Pairem, "CDCLayerId_Pairem[nCDCHit_Pairem]/I");
    trdata->Branch("MaxCDCLayerId_Pairem", &MaxCDCLayerId_Pairem, "MaxCDCLayerId_Pairem/I");

    trdata->Branch("ifSingleTurn_Pairep", &ifSingleTurn_Pairep, "ifSingleTurn_Pairep/O");
    trdata->Branch("ifMultiTurn_Pairep", &ifMultiTurn_Pairep, "ifMultiTurn_Pairep/O");

    trdata->Branch("CDCHitX", CDCHitX, "CDCHitX[nCDCHit]/D");
    trdata->Branch("CDCHitY", CDCHitY, "CDCHitY[nCDCHit]/D");
    trdata->Branch("CDCHitZ", CDCHitZ, "CDCHitZ[nCDCHit]/D");
    trdata->Branch("CDCHitT", CDCHitT, "CDCHitT[nCDCHit]/D");
    trdata->Branch("CDCEDep", CDCEDep, "CDCEDep[nCDCHit]/D");
    trdata->Branch("TurnId",  TurnId, "TurnId[nCDCHit]/I");
    trdata->Branch("TurnNumber",  &TurnNumber, "TurnNumber/I");

    trdata->Branch("nCRKHit", &nCRKHit, "nCRKHit/I");
    trdata->Branch("CRKHitX", CRKHitX, "CRKHitX[nCRKHit]/D");
    trdata->Branch("CRKHitY", CRKHitY, "CRKHitY[nCRKHit]/D");
    trdata->Branch("CRKHitZ", CRKHitZ, "CRKHitZ[nCRKHit]/D");
    trdata->Branch("CRKHitT", CRKHitT, "CRKHitT[nCRKHit]/D");
    trdata->Branch("CRKEDep", CRKEDep, "CRKEDep[nCRKHit]/D");

    trdata->Branch("nSTLHit", &nSTLHit, "nSTLHit/I");
    trdata->Branch("STLHitX", STLHitX, "STLHitX[nSTLHit]/D");
    trdata->Branch("STLHitY", STLHitY, "STLHitY[nSTLHit]/D");
    trdata->Branch("STLHitZ", STLHitZ, "STLHitZ[nSTLHit]/D");
    trdata->Branch("STLHitT", STLHitT, "STLHitT[nSTLHit]/D");
    trdata->Branch("STLEDep", STLEDep, "STLEDep[nSTLHit]/D");

    trdata->Branch("nTrigCRK", &nTrigCRK, "nTrigCRK/I");
    trdata->Branch("TrigCRK", TrigCRK, "TrigCRK[nTrigCRK]/I");
    trdata->Branch("TrigCRKTime", TrigCRKTime, "TrigCRK[nTrigCRK]/D");
    trdata->Branch("nTrigSTL", &nTrigSTL, "nTrigSTL/I");
    trdata->Branch("TrigSTL", TrigSTL, "TrigSTL[nTrigSTL]/I");
    trdata->Branch("TrigSTLTime", TrigSTLTime, "TrigSTL[nTrigSTL]/D");

    trdata->Branch("nCALCDCHit", &nCALCDCHit, "nCALCDCHit/I");
    trdata->Branch("CDCDriftDist", CDCDriftDist, "CDCDriftDist[nCALCDCHit]/D");
    trdata->Branch("CDCCharge", CDCCharge, "CDCCharge[nCALCDCHit]/I");
    trdata->Branch("WireEnd0X", WireEnd0X, "WireEnd0X[nCALCDCHit]/D");
    trdata->Branch("WireEnd0Y", WireEnd0Y, "WireEnd0Y[nCALCDCHit]/D");
    trdata->Branch("WireEnd0Z", WireEnd0Z, "WireEnd0Z[nCALCDCHit]/D");
    trdata->Branch("WireEnd1X", WireEnd1X, "WireEnd1X[nCALCDCHit]/D");
    trdata->Branch("WireEnd1Y", WireEnd1Y, "WireEnd1Y[nCALCDCHit]/D");
    trdata->Branch("WireEnd1Z", WireEnd1Z, "WireEnd1Z[nCALCDCHit]/D");
    trdata->Branch("WireLayerId", WireLayerId, "WireLayerId[nCALCDCHit]/I");
    trdata->Branch("WireMaxLayerId", &WireMaxLayerId, "WireMaxLayerId/I");
    trdata->Branch("WireId", WireId, "WireId[nCALCDCHit]/I");
    
    trdata->Branch("genTrX", &genTrX, "genTrX/D");
    trdata->Branch("genTrY", &genTrY, "genTrY/D");
    trdata->Branch("genTrZ", &genTrZ, "genTrZ/D");
    trdata->Branch("genTrT", &genTrT, "genTrT/D");
    trdata->Branch("genTrPx", &genTrPx, "genTrPx/D");
    trdata->Branch("genTrPy", &genTrPy, "genTrPy/D");
    trdata->Branch("genTrPz", &genTrPz, "genTrPz/D");
    trdata->Branch("genTrE", &genTrE, "genTrE/D");

    trdata->Branch("CDCEnterX", &CDCEnterX, "CDCEnterX/D");
    trdata->Branch("CDCEnterY", &CDCEnterY, "CDCEnterY/D");
    trdata->Branch("CDCEnterZ", &CDCEnterZ, "CDCEnterZ/D");
    trdata->Branch("CDCEnterT", &CDCEnterT, "CDCEnterT/D");
    trdata->Branch("CDCEnterPx", &CDCEnterPx, "CDCEnterPx/D");
    trdata->Branch("CDCEnterPy", &CDCEnterPy, "CDCEnterPy/D");
    trdata->Branch("CDCEnterPz", &CDCEnterPz, "CDCEnterPz/D");
    trdata->Branch("CDCEnterE",  &CDCEnterE,  "CDCEnterE/D");

    trdata->Branch("PairVertexX", &PairVertexX, "PairVertexX/D");
    trdata->Branch("PairVertexY", &PairVertexY, "PairVertexY/D");
    trdata->Branch("PairVertexZ", &PairVertexZ, "PairVertexZ/D");
    trdata->Branch("PairVertexT", &PairVertexT, "PairVertexT/D");

    trdata->Branch("Pairep_TrackId", &Pairep_TrackId, "Pairep_TrackId/I");
    trdata->Branch("Pairep_genTrPx", &Pairep_genTrPx, "Pairep_genTrPx/D");
    trdata->Branch("Pairep_genTrPy", &Pairep_genTrPy, "Pairep_genTrPy/D");
    trdata->Branch("Pairep_genTrPz", &Pairep_genTrPz, "Pairep_genTrPz/D");
    trdata->Branch("Pairep_genTrE", &Pairep_genTrE, "Pairep_genTrE/D");
    trdata->Branch("Pairep_TurnNumber", &Pairep_TurnNumber, "Pairep_TurnNumber/I");

    trdata->Branch("Pairem_TrackId", &Pairem_TrackId, "Pairem_TrackId/I");
    trdata->Branch("Pairem_genTrPx", &Pairem_genTrPx, "Pairem_genTrPx/D");
    trdata->Branch("Pairem_genTrPy", &Pairem_genTrPy, "Pairem_genTrPy/D");
    trdata->Branch("Pairem_genTrPz", &Pairem_genTrPz, "Pairem_genTrPz/D");
    trdata->Branch("Pairem_genTrE", &Pairem_genTrE, "Pairem_genTrE/D");
    trdata->Branch("ifPairVertexAtTarget", &ifPairVertexAtTarget, "ifPairVertexAtTarget/O");
    trdata->Branch("ifPairProdOccurs", &ifPairProdOccurs, "ifPairProdOccurs/O");

    trdata->Branch("ComptonVertexX", &ComptonVertexX, "ComptonVertexX/D");
    trdata->Branch("ComptonVertexY", &ComptonVertexY, "ComptonVertexY/D");
    trdata->Branch("ComptonVertexZ", &ComptonVertexZ, "ComptonVertexZ/D");
    trdata->Branch("ComptonVertexT", &ComptonVertexT, "ComptonVertexT/D");
    trdata->Branch("Comptonem_genTrPx", &Comptonem_genTrPx, "Comptonem_genTrPx/D");
    trdata->Branch("Comptonem_genTrPy", &Comptonem_genTrPy, "Comptonem_genTrPy/D");
    trdata->Branch("Comptonem_genTrPz", &Comptonem_genTrPz, "Comptonem_genTrPz/D");
    trdata->Branch("Comptonem_genTrE", &Comptonem_genTrE, "Comptonem_genTrE/D");
    trdata->Branch("ifComptonVertexAtTarget", &ifComptonVertexAtTarget, "ifComptonVertexAtTarget/O");
    trdata->Branch("ifComptonId2Occurs", &ifComptonId2Occurs, "ifComptonId2Occurs/O");

  }

  bool operator () (COMET::ICOMETEvent& event) 
  {
    /*------------------------- Nullifying --------------------------------*/

    nCDCHit_Pairep=0;
    nCDCHit_Pairem=0;
    nCDCHit_Comptonem=0;
    memset (ifCDC_Pairep, 0, sizeof (ifCDC_Pairep));
    memset (ifCDC_Pairem, 0, sizeof (ifCDC_Pairem));
    memset (ifCDC_Comptonem, 0, sizeof (ifCDC_Comptonem));
    memset (CDCLayerId_Pairep, 0, sizeof (CDCLayerId_Pairep));
    MaxCDCLayerId_Pairep=0;
    memset (CDCLayerId_Pairem, 0, sizeof (CDCLayerId_Pairem));
    MaxCDCLayerId_Pairem=0;
    ifSingleTurn_Pairep=0;
    ifMultiTurn_Pairep=0;

    nCDCHit=0;
    memset (CDCHitX, 0, sizeof (CDCHitX));
    memset (CDCHitY, 0, sizeof (CDCHitY));
    memset (CDCHitZ, 0, sizeof (CDCHitZ));
    memset (CDCHitT, 0, sizeof (CDCHitT));
    memset (CDCEDep, 0, sizeof (CDCEDep));
    memset (TurnId, 0, sizeof (TurnId));
    TurnNumber=0;

    nCRKHit=0;
    memset (CRKHitX, 0, sizeof (CRKHitX));
    memset (CRKHitY, 0, sizeof (CRKHitY));
    memset (CRKHitZ, 0, sizeof (CRKHitZ));
    memset (CRKHitT, 0, sizeof (CRKHitT));
    memset (CRKEDep, 0, sizeof (CRKEDep));

    nSTLHit=0;
    memset (STLHitX, 0, sizeof (STLHitX));
    memset (STLHitY, 0, sizeof (STLHitY));
    memset (STLHitZ, 0, sizeof (STLHitZ));
    memset (STLHitT, 0, sizeof (STLHitT));
    memset (STLEDep, 0, sizeof (STLEDep));

    nTrigCRK=0;
    memset (TrigCRK, 0, sizeof (TrigCRK));
    memset (TrigCRKTime, 0, sizeof (TrigCRKTime));    
    nTrigSTL=0;
    memset (TrigSTL, 0, sizeof (TrigSTL));
    memset (TrigSTLTime, 0, sizeof (TrigSTLTime));

    nCALCDCHit=0;
    memset (CDCDriftDist, 0, sizeof(CDCDriftDist));
    memset (CDCCharge, 0, sizeof(CDCCharge));
    memset (WireEnd0X, 0, sizeof(WireEnd0X));
    memset (WireEnd0Y, 0, sizeof(WireEnd0Y));
    memset (WireEnd0Z, 0, sizeof(WireEnd0Z));
    memset (WireEnd1X, 0, sizeof(WireEnd1X));
    memset (WireEnd1Y, 0, sizeof(WireEnd1Y));
    memset (WireEnd1Z, 0, sizeof(WireEnd1Z));
    memset (WireLayerId, 0, sizeof(WireLayerId));
    WireMaxLayerId=0;
    memset (WireId, 0, sizeof(WireId));

    genTrX=0;
    genTrY=0;
    genTrZ=0;
    genTrT=0;
    genTrPx=0;
    genTrPy=0;
    genTrPz=0;
    genTrE=0;

    CDCEnterX=0;
    CDCEnterY=0;
    CDCEnterZ=0;
    CDCEnterY=0;
    CDCEnterPx=0;
    CDCEnterPy=0;
    CDCEnterPz=0;
    CDCEnterE=0;

    PairVertexX=0;
    PairVertexY=0;
    PairVertexZ=0;
    PairVertexT=0;
    
    Pairep_TrackId=0;
    Pairep_genTrPx=0;
    Pairep_genTrPy=0;
    Pairep_genTrPz=0;
    Pairep_genTrE=0;
    Pairep_TurnNumber=0;

    Pairem_TrackId=0;
    Pairem_genTrPx=0;
    Pairem_genTrPx=0;
    Pairem_genTrPx=0;
    Pairem_genTrE=0;
    ifPairVertexAtTarget=0;
    ifPairProdOccurs=0;

    ComptonVertexX=0;
    ComptonVertexY=0;
    ComptonVertexZ=0;
    ComptonVertexT=0;
    Comptonem_genTrPx=0;
    Comptonem_genTrPy=0;
    Comptonem_genTrPz=0;
    Comptonem_genTrE=0;
    ifComptonVertexAtTarget=0;
    ifComptonId2Occurs=0;

    /*------------------------------------------------------------------*/

    eventId = event.GetEventId();  

    gGeoManager=COMET::IOADatabase::Get().Geometry();
    COMET::IHandle<COMET::IG4TrajectoryContainer> Trajectories = event.Get<COMET::IG4TrajectoryContainer>("truth/G4Trajectories");
    COMET::IG4TrajectoryContainer *TrajCont = GetPointer(Trajectories);

    if(TrajCont->empty()){
      std::cout<< "Result is not found" << std::endl;
    }
    
    if(!TrajCont->empty()){
      
      for(COMET::IG4TrajectoryContainer::const_iterator seg = TrajCont->begin(); seg != TrajCont->end(); seg++){
	
	COMET::IG4Trajectory traj = (*seg).second;


	/*------  Primary Gamma Information -----*/

	if (traj.GetTrackId() == 1){
	  TVector3 iniPos = traj.GetInitialPosition().Vect()*(1/unit::cm);
	  TVector3 iniMom = traj.GetInitialMomentum().Vect()*(1/unit::MeV);
	  
	  genTrX=iniPos(0);
	  genTrY=iniPos(1);
	  genTrZ=iniPos(2);
	  genTrT=traj.GetInitialPosition()(3);

	  genTrPx=iniMom(0);
	  genTrPy=iniMom(1);
	  genTrPz=iniMom(2);
	  genTrE=traj.GetInitialMomentum()(3);	  

	  COMET::IG4Trajectory::Points trajPointSet = traj.GetTrajectoryPoints();
	  for(COMET::IG4Trajectory::Points::iterator trajIter = trajPointSet.begin(); trajIter!=trajPointSet.end(); trajIter++ ){
	    COMET::IG4TrajectoryPoint trajPoint = *trajIter;
	    TVector3 tmpPosition = trajPoint.GetPosition().Vect()*(1/unit::cm);
	    //std::cout << GetNode(tmpPosition)->GetName() << std::endl;
	    if (TString(GetNode(tmpPosition)->GetName())=="CDCSenseLayer_0_0") {
	      TVector3 enterPos = trajPoint.GetPosition().Vect()*(1/unit::cm);
	      TVector3 enterMom = trajPoint.GetMomentum()*(1/unit::MeV);
	      CDCEnterX = enterPos(0);
	      CDCEnterY = enterPos(1);
	      CDCEnterZ = enterPos(2);
	      CDCEnterT = trajPoint.GetPosition()(3)*(1/unit::ns);
	      CDCEnterPx= enterMom(0);
	      CDCEnterPy= enterMom(1);
	      CDCEnterPz= enterMom(2);

	      //std::cout << CDCEnterX << "   " << CDCEnterY << "   " << CDCEnterZ << "   " << CDCEnterT << std::endl;
	      //std::cout << CDCEnterPx << "   " << CDCEnterPy << "   " << CDCEnterPz << "   "<< std::endl;
	      break;
	    }
	  }	  
	}

	/*------  Pair Production Information -----*/

	/* Pair Producted e+ */

	if (traj.GetParentId()==1 && traj.GetPDGEncoding()==-11){

	  Pairep_TrackId = traj.GetTrackId();
	  Pairem_TrackId = Pairep_TrackId-1;

	  //	  std::cout << Pairem_TrackId << "  " << Pairep_TrackId << std::endl;

	  ifPairProdOccurs=true;

	  TVector3 iniPos = traj.GetInitialPosition().Vect()*(1/unit::cm);
	  TVector3 iniMom = traj.GetInitialMomentum().Vect()*(1/unit::MeV);
	  TGeoNode* iniVolume = GetNode(iniPos);

	  if (TString(iniVolume->GetName()).Contains("Target")){
            ifPairVertexAtTarget=true;
	  }

	  PairVertexX=iniPos(0);
	  PairVertexY=iniPos(1);
	  PairVertexZ=iniPos(2);
	  PairVertexT=traj.GetInitialPosition()(3);

	  Pairep_genTrPx=iniMom(0);
	  Pairep_genTrPy=iniMom(1);
	  Pairep_genTrPz=iniMom(2);
	  Pairep_genTrE=traj.GetInitialMomentum()(3);	  	  
	 
	}

      }	
    }

    /* Pair Producted e- */

    if (ifPairProdOccurs==true){
      COMET::IHandle<COMET::IG4Trajectory> Pairem_traj_Handle = TrajCont->GetTrajectory(Pairem_TrackId);
      COMET::IG4Trajectory* Pairem_traj = GetPointer(Pairem_traj_Handle);

      TVector3 iniPos = Pairem_traj->GetInitialPosition().Vect()*(1/unit::cm);
      TVector3 iniMom = Pairem_traj->GetInitialMomentum().Vect()*(1/unit::MeV);

      Pairem_genTrPx=iniMom(0);
      Pairem_genTrPy=iniMom(1);
      Pairem_genTrPz=iniMom(2);
      Pairem_genTrE=Pairem_traj->GetInitialMomentum()(3);	  	        
    }

    /* Compton Scattering e- (Only for TrackId==2) */
    
    if (!TrajCont->empty()){
      if (ifPairProdOccurs==false || (ifPairProdOccurs==true && Pairem_TrackId>2)){
	COMET::IHandle<COMET::IG4Trajectory> Compton_traj_Handle = TrajCont->GetTrajectory(2);
	COMET::IG4Trajectory* Compton_traj = GetPointer(Compton_traj_Handle);
	
	if (Compton_traj){
	  
	  if (Compton_traj->GetParentId()==1 && Compton_traj->GetPDGEncoding()==11){
	    ifComptonId2Occurs=true;
	    
	    TVector3 iniPos = Compton_traj->GetInitialPosition().Vect()*(1/unit::cm);
	    TVector3 iniMom = Compton_traj->GetInitialMomentum().Vect()*(1/unit::MeV);
	    TGeoNode* iniVolume = GetNode(iniPos);
	    
	    if (TString(iniVolume->GetName()).Contains("Target")){
	      ifComptonVertexAtTarget=true;
	    }
	    
	    ComptonVertexX=iniPos(0);
	    ComptonVertexY=iniPos(1);
	    ComptonVertexZ=iniPos(2);
	    ComptonVertexT=Compton_traj->GetInitialPosition()(3);
	    
	    Comptonem_genTrPx=iniMom(0);
	    Comptonem_genTrPy=iniMom(1);
	    Comptonem_genTrPz=iniMom(2);
	    Comptonem_genTrE=Compton_traj->GetInitialMomentum()(3);	  	  
	    
	  }
	}
      }
    }    
    /* ----------------------------
     |                              |
     |      CDC Hit Information     |
     |                              |
       ---------------------------- */
    
    std::vector<TString> CDCHit_PairepGeometry;
    CDCHit_PairepGeometry.push_back("Default");

    std::vector<TString> CDCHitGeometry_Prim;
    CDCHitGeometry_Prim.push_back("Default");
    int NumOfCDC_0_Prim=0;
    int NumOfCDC_1_Prim=0;
    
    COMET::IHandle<COMET::IG4HitContainer> CDCHitCont = event.Get<COMET::IG4HitContainer>("truth/g4Hits/CDC");
    if (CDCHitCont){
      for(COMET::IG4HitContainer::const_iterator hitSeg = CDCHitCont->begin(); hitSeg != CDCHitCont->end(); ++hitSeg) {
	COMET::IG4HitSegment* tmpSeg = dynamic_cast<COMET::IG4HitSegment*>(*hitSeg);
	
        if (tmpSeg){
	  //Double_t TrLen=tmpSeg->GetTrackLength()*(1/unit::cm);	      
	  //std::cout << TrLen << std::endl;
	  
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

	  /*
	  if(!COMET::IGeomInfo::DetectorSolenoid().GetDetPositionInDSCoordinate(hitPos, hitPos)){
	    continue;
	    std::cout << "Mislocated hit is detected" << std::endl;
	  }
	  */
	  
	  if (ifPairProdOccurs==1 && trajContributors[0]==Pairep_TrackId && trajContributors.size()==1){
	    if (CDCHit_PairepGeometry.back() != geoName){
	      CDCHit_PairepGeometry.push_back(geoName);
	      
	      //if (eventId==8503) std::cout << geoName << std::endl;
	    } 
	  }
	  
	    
          if (geoName.Contains("CDCSenseLayer")){
	    CDCHitX[nCDCHit]=hitPos(0);
	    CDCHitY[nCDCHit]=hitPos(1);
	    CDCHitZ[nCDCHit]=hitPos(2);
	    CDCHitT[nCDCHit]=hitT;
	    CDCEDep[nCDCHit]=tmpSeg->GetEnergyDeposit();

	    if (*trajContributors.begin()==1 && trajContributors.size()==1){ // If it is Primary (Signal)
	      if (CDCHitGeometry_Prim.back() != geoName){
		CDCHitGeometry_Prim.push_back(geoName);
		if (geoName=="CDCSenseLayer_0_0"){
		  NumOfCDC_0_Prim++;
		}
		if (geoName=="CDCSenseLayer_1_0"){
		  NumOfCDC_1_Prim++;
		}
	      }
	    }
	    	    
	    TurnId[nCDCHit]=NumOfCDC_0_Prim-1;

	    if (ifPairProdOccurs==1){
	      
	      if (trajContributors[0]==Pairep_TrackId && trajContributors.size()==1){  
		TString ep_geoName = geoName;
		Int_t delimFirst= (ep_geoName).First('_');
		Int_t delimLast = (ep_geoName).Last('_');
		Int_t CDCIndex= (((ep_geoName).Remove(delimLast,100)).Remove(0,delimFirst+1)).Atoi();
		CDCLayerId_Pairep[nCDCHit_Pairep]=CDCIndex+1;
		nCDCHit_Pairep++;
		ifCDC_Pairep[nCDCHit]=true;
	      }
	      
	      else if (trajContributors[0]==Pairem_TrackId && trajContributors.size()==1){  
		TString em_geoName = geoName;
		Int_t delimFirst= (em_geoName).First('_');
		Int_t delimLast = (em_geoName).Last('_');
		Int_t CDCIndex= (((em_geoName).Remove(delimLast,100)).Remove(0,delimFirst+1)).Atoi();
		CDCLayerId_Pairem[nCDCHit_Pairem]=CDCIndex+1;
		nCDCHit_Pairem++;
		ifCDC_Pairem[nCDCHit]=true;
		}
	      
	    }
	    
	    if (ifComptonId2Occurs==1){
	      if (trajContributors[0]==2 && trajContributors.size()==1){
		nCDCHit_Comptonem++;
		ifCDC_Comptonem[nCDCHit]=true;
	      }
	    }


	    nCDCHit++;
	  }	  
        }	
      }
    }

    TurnNumber = NumOfCDC_1_Prim/2;
        
    if (ifPairProdOccurs==1 && nCDCHit_Pairem>0){
      for (Int_t i=0; i<nCDCHit_Pairem; i++){
	if (CDCLayerId_Pairem[i]>MaxCDCLayerId_Pairem){
	  MaxCDCLayerId_Pairem=CDCLayerId_Pairem[i];
	}
      }    
    }

    if (ifPairProdOccurs==1 && nCDCHit_Pairep>0){
      for (Int_t i=0; i<nCDCHit_Pairep; i++){
	if (CDCLayerId_Pairep[i]>MaxCDCLayerId_Pairep){
	  MaxCDCLayerId_Pairep=CDCLayerId_Pairep[i];
	}
      }
       
      int NumOfCDC_0=0;
      int NumOfCDC_1=0;
      int NumOfCDC_2=0;
      
      
      for (std::vector<TString>::iterator geo=CDCHit_PairepGeometry.begin(); geo!=CDCHit_PairepGeometry.end(); geo++)
	{
	  if (*geo=="CDCSenseLayer_0_0")
	    {NumOfCDC_0++;}
	  else if (*geo=="CDCSenseLayer_1_0")
	    {NumOfCDC_1++;}
	  else if (*geo=="CDCSenseLayer_2_0")
	    {NumOfCDC_2++;}
	  
	}
      
      if (NumOfCDC_0==2 && NumOfCDC_1==2 && NumOfCDC_2==2)     {ifSingleTurn_Pairep=true; }
      else if (NumOfCDC_0>=3 && NumOfCDC_1>=4 && NumOfCDC_2>=4) {ifMultiTurn_Pairep=true; }

      Pairep_TurnNumber = (NumOfCDC_1+1)/2;
    }
        
    /* -----------------------------
     |                               |
     |  CTH Trigger Hit Information  |
     |                               |
       ----------------------------- */
    
    COMET::IHandle<COMET::IG4HitContainer> CTHHitCont = event.Get<COMET::IG4HitContainer>("truth/g4Hits/CTH");
    if (CTHHitCont){
      for(COMET::IG4HitContainer::const_iterator hitSeg = CTHHitCont->begin(); hitSeg != CTHHitCont->end(); ++hitSeg) {
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


	  COMET::IGeometryId CTHId;
	  COMET::IOADatabase::Get().GeomId().GetGeometryId(hitPos(0), hitPos(1), hitPos(2), CTHId);
	  COMET::ICTHChannelId CTHChannelId;
	  COMET::IGeometryDatabase::Get().GetChannelId(CTHChannelId, CTHId);

	  /*
	  if(!COMET::IGeomInfo::DetectorSolenoid().GetDetPositionInDSCoordinate(hitPos, hitPos)){
	    continue;
	    std::cout << "Mislocated hit is detected" << std::endl;
	  }
	  */

	  //std::cout << CTHChannelId.GetScint() << "   " << CTHChannelId.GetLightGuide() << "   " << geoName << "   " << CTHChannelId.GetCounter() << "   " << CTHChannelId.GetModule() << std::endl;
       
          if (geoName.Contains("Cherenkov_pv"))
	    {

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
	    }


	  else if (geoName.Contains("Scintillator_pv"))
	    {
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


     /* ------------------------------
     |                               |
     |  CDC Detector Response        |
     |                               |
     --------------------------------*/

    // IMCHit Container

    std::map <COMET::IG4HitSegment*, std::vector<COMET::IMCHit*> > HitMap;

    COMET::IHandle<COMET::IHitSelection> hitHandle = event.Get<COMET::IHitSelection>("./hits/mcCDC");
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
	}

	else if ( HitMap.find(g4Contributor) != HitMap.end()){ // if exist in map
	  HitMap[g4Contributor].push_back(mcHit);
	  //std::cout << HitMap[g4Contributor].size() << std::endl;
	}
      }
    }

    for (std::map< COMET::IG4HitSegment* ,std::vector<COMET::IMCHit* > >::iterator it=HitMap.begin(); it!=HitMap.end(); ++it){
      COMET::IG4HitSegment* g4HitSeg = it->first;
      std::vector<COMET::IMCHit*> MCHitVector = it->second; 
      CDCCharge[nCALCDCHit]=MCHitVector.size();

      for (int i_hit=0; i_hit<MCHitVector.size(); i_hit++){
	TVectorD wireMes(7);	
	if (i_hit==0){
	  COMET::IGeometryId geomId = MCHitVector[i_hit]->GetGeomId();
	  int wire = COMET::IGeomInfo::Get().CDC().GeomIdToWire(geomId);
	  WireEnd0X[nCALCDCHit]    = COMET::IGeomInfo::Get().CDC().GetWireEnd0(wire).X();
	  WireEnd0Y[nCALCDCHit]    = COMET::IGeomInfo::Get().CDC().GetWireEnd0(wire).Y();
	  WireEnd0Z[nCALCDCHit]    = COMET::IGeomInfo::Get().CDC().GetWireEnd0(wire).Z();
	  WireEnd1X[nCALCDCHit]    = COMET::IGeomInfo::Get().CDC().GetWireEnd1(wire).X();
	  WireEnd1Y[nCALCDCHit]    = COMET::IGeomInfo::Get().CDC().GetWireEnd1(wire).Y();
	  WireEnd1Z[nCALCDCHit]    = COMET::IGeomInfo::Get().CDC().GetWireEnd1(wire).Z();
	  WireId[nCALCDCHit]       = wire;	  
	  int layer = COMET::IGeomInfo::Get().CDC().GetLayer(wire);
	  WireLayerId[nCALCDCHit]  = layer;
	  if (WireMaxLayerId<layer) WireMaxLayerId = layer;

	  TVector3 g4HitPos;
          g4HitPos(0)=0.5 * (g4HitSeg->GetStopX()*(1/unit::cm)+g4HitSeg->GetStartX()*(1/unit::cm));
          g4HitPos(1)=0.5 * (g4HitSeg->GetStopY()*(1/unit::cm)+g4HitSeg->GetStartY()*(1/unit::cm));
          g4HitPos(2)=0.5 * (g4HitSeg->GetStopZ()*(1/unit::cm)+g4HitSeg->GetStartZ()*(1/unit::cm));

	  TVector3 local;
	  if (!COMET::IGeomInfo::Get().CDC().GetDistanceFromWire(g4HitPos, wire, local)) continue;
	  CDCDriftDist[nCALCDCHit]=hypot(local.x(),local.y());
	}	
	else if (i_hit != 0 ){	  
	  //CDCDriftDist[nCALCDCHit]+= MCHitVector[i_hit]->GetDriftDistance()*(1/unit::cm);
	}
      }
      //std::cout << CDCDriftDist[nCALCDCHit] << std::endl;
      nCALCDCHit++;
    }


    /*
    COMET::IChannelId tmpchanId;
    COMET::IHandle<COMET::IHitSelection> hitHandle = event.Get<COMET::IHitSelection>("./hits/mcCDC");
    if (hitHandle){
      COMET::IHitSelection *hits = GetPointer(hitHandle);
      for (COMET::IHitSelection::const_iterator hitSeg = hits->begin(); hitSeg != hits->end(); hitSeg++){
	COMET::IChannelId chanId = (*hitSeg)->GetChannelId();
	COMET::IGeometryId geomId = (*hitSeg)->GetGeomId();

	/// Find SimG4 Hit contributors
	COMET::IMCHit *mcHit = dynamic_cast<COMET::IMCHit*>(GetPointer(*hitSeg));
	std::vector<COMET::IG4VHit*> hitContributors = mcHit->GetContributors();
	COMET::IG4HitSegment* g4Contributor = dynamic_cast<COMET::IG4HitSegment*>(hitContributors.at(0));

	int wire = COMET::IGeomInfo::Get().CDC().GeomIdToWire(geomId);

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
	Double_t Drfit = (*hitSeg)->GetDriftDistance();
	
	
	//if(!COMET::IGeomInfo::DetectorSolenoid().GetDetPositionInDSCoordinate(wireend0, wireend0)){
	//  continue;
	//  std::cout << "MisIdentifided wire is detected" << std::endl;}
	//if(!COMET::IGeomInfo::DetectorSolenoid().GetDetPositionInDSCoordinate(wireend1, wireend1)){
	//  continue;
	//  std::cout << "MisIdentifided wire is detected" << std::endl;}
	

	int layer = COMET::IGeomInfo::Get().CDC().GetLayer(wire);
	if (WireMaxLayerId<layer){
	  WireMaxLayerId = layer;}

        int wireid;
        if(COMET::GeomId::CDC::IsActiveSenseWire(geomId)==1){
          wireid = COMET::IGeomInfo::Get().CDC().GeomIdToWire(geomId);
        }
        else {
	  std::cout << "Not SenseWire" << std::endl;
        }
	
	///Merge ionized electrons from same hits

	if (hitSeg == hits->begin()){
	  tmpchanId = chanId;
	  CDCCharge[nCALCDCHit]++;
	  WireEnd0X[nCALCDCHit]=wireend0(0);
	  WireEnd0Y[nCALCDCHit]=wireend0(1);
	  WireEnd0Z[nCALCDCHit]=wireend0(2);
	  WireEnd1X[nCALCDCHit]=wireend1(0);
	  WireEnd1Y[nCALCDCHit]=wireend1(1);
	  WireEnd1Z[nCALCDCHit]=wireend1(2);
	  CDCDriftDist[nCALCDCHit]=wireMes[6];
	  WireLayerId[nCALCDCHit]=layer;
	  WireId[nCALCDCHit]=wireid;
	}

	// When the next hit is at same Id                                                             
	if (chanId == tmpchanId){
	  CDCDriftDist[nCALCDCHit]+=wireMes[6];
	  CDCCharge[nCALCDCHit]++;
	}

	// When the next hit generates at another Channel Id                                           
	else if (chanId != tmpchanId){
	  // Average drift distance
	  CDCDriftDist[nCALCDCHit]=CDCDriftDist[nCALCDCHit]/CDCCharge[nCALCDCHit];

	  // Print to test values
	  //std::cout << WireEnd0X[nCALCDCHit] << "  " << WireEnd0Y[nCALCDCHit] << "   " <<WireEnd0X[nCALCDCHit] << "   " << WireEnd1X[nCALCDCHit] << "   " << WireEnd1Y[nCALCDCHit] << "   " << WireEnd1Z[nCALCDCHit] << "   " << CDCDriftDist[nCALCDCHit] <<  std::endl;

	  
	  // Change channel Id
	  tmpchanId = chanId;
	  nCALCDCHit++;
	  CDCCharge[nCALCDCHit]++;
	  WireEnd0X[nCALCDCHit]=wireend0(0);
	  WireEnd0Y[nCALCDCHit]=wireend0(1);
	  WireEnd0Z[nCALCDCHit]=wireend0(2);
	  WireEnd1X[nCALCDCHit]=wireend1(0);
	  WireEnd1Y[nCALCDCHit]=wireend1(1);
	  WireEnd1Z[nCALCDCHit]=wireend1(2);
	  CDCDriftDist[nCALCDCHit]=wireMes[6];
	  WireLayerId[nCALCDCHit]=layer;
	  WireId[nCALCDCHit]=wireid;
	}
      }
    }
*/

    // IDigit Container

    COMET::IHandle<COMET::IDigitContainer> digitContainer = event.Get<COMET::IDigitContainer>("./digits/CDC");
    if (digitContainer){
      for (COMET::IDigitContainer::const_iterator digitIt = digitContainer->begin(); digitIt != digitContainer->end(); digitIt++){
	COMET::ICDCmcDigit* CDCmcDigit = dynamic_cast<COMET::ICDCmcDigit*>(*digitIt);
	COMET::IChannelId chanId = CDCmcDigit->GetChannelId();
	std::vector<short> adcs  = CDCmcDigit->GetADCs();
	/*
	for (int i=0; i<adcs.size(); i++){
	  std::cout << adcs.at(i) << "   ";
	}
	std::cout << std::endl;
	*/
      }
    }

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
 
  Int_t eventId;    

  Int_t nCDCHit_Pairep;
  Int_t nCDCHit_Pairem;
  Int_t nCDCHit_Comptonem;
  Bool_t ifCDC_Pairep[10000];
  Bool_t ifCDC_Pairem[10000];
  Bool_t ifCDC_Comptonem[10000];
  Int_t CDCLayerId_Pairep[10000];
  Int_t MaxCDCLayerId_Pairep;
  Int_t CDCLayerId_Pairem[10000];
  Int_t MaxCDCLayerId_Pairem;
  Bool_t ifSingleTurn_Pairep;
  Bool_t ifMultiTurn_Pairep;

  Int_t nCDCHit;
  Double_t CDCHitX[100000];  
  Double_t CDCHitY[100000]; 
  Double_t CDCHitZ[100000];  
  Double_t CDCHitT[100000]; 
  Double_t CDCEDep[100000];
  Int_t TurnId[100000];
  Int_t TurnNumber;

  Int_t nCRKHit;
  Double_t CRKHitX[10000];  
  Double_t CRKHitY[10000];  
  Double_t CRKHitZ[10000];  
  Double_t CRKHitT[10000];
  Double_t CRKEDep[10000];

  Int_t nSTLHit;
  Double_t STLHitX[10000];
  Double_t STLHitY[10000];
  Double_t STLHitZ[10000];
  Double_t STLHitT[10000];
  Double_t STLEDep[10000];

  Int_t nTrigSTL;
  Int_t nTrigCRK;
  Int_t TrigSTL[1000];
  Double_t TrigSTLTime[1000];
  Int_t TrigCRK[1000];
  Double_t TrigCRKTime[1000];

  /*------------- Detector Response Part ------------*/

  Int_t nCALCDCHit;
  Double_t CDCDriftDist[10000];
  Int_t CDCCharge[10000];
  Double_t WireEnd0X[10000];
  Double_t WireEnd0Y[10000];
  Double_t WireEnd0Z[10000];
  Double_t WireEnd1X[10000];
  Double_t WireEnd1Y[10000];
  Double_t WireEnd1Z[10000];
  Int_t WireLayerId[10000];
  Int_t WireMaxLayerId;
  Int_t WireId[10000];

  /*------------- Trajectory Class Root Member -----------------*/

  Double_t genTrX;
  Double_t genTrY;
  Double_t genTrZ;
  Double_t genTrT;
  Double_t genTrPx;
  Double_t genTrPy;
  Double_t genTrPz;
  Double_t genTrE;

  Double_t CDCEnterX;
  Double_t CDCEnterY;
  Double_t CDCEnterZ;
  Double_t CDCEnterT;
  Double_t CDCEnterPx;
  Double_t CDCEnterPy;
  Double_t CDCEnterPz;
  Double_t CDCEnterE;
    
  Double_t PairVertexX;
  Double_t PairVertexY;
  Double_t PairVertexZ;
  Double_t PairVertexT;

  Bool_t ifPairProdOccurs;
  Bool_t ifPairVertexAtTarget;

  Int_t Pairep_TrackId;
  Double_t Pairep_genTrPx;
  Double_t Pairep_genTrPy;
  Double_t Pairep_genTrPz;
  Double_t Pairep_genTrE;
  Int_t Pairep_TurnNumber;

  Int_t Pairem_TrackId;
  Double_t Pairem_genTrPx;
  Double_t Pairem_genTrPy;
  Double_t Pairem_genTrPz;
  Double_t Pairem_genTrE;

  Bool_t ifComptonVertexAtTarget;
  Bool_t ifComptonId2Occurs;
  Double_t ComptonVertexX;
  Double_t ComptonVertexY;
  Double_t ComptonVertexZ;
  Double_t ComptonVertexT;
  Double_t Comptonem_genTrPx;
  Double_t Comptonem_genTrPy;
  Double_t Comptonem_genTrPz;
  Double_t Comptonem_genTrE;
   
  /*-------------------------------------------*/

  TGeoManager* gGeoManager;

};

int main(int argc, char **argv) {
  SignalTracking userCode;
  cometEventLoop(argc,argv,userCode);
}
