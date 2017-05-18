#include <IMCTrigger.hxx>
#include <ITracking.hxx>
#include <IHoughTransform.hxx>
#include <IGenFitting.hxx>
#include <IHelixTracker.hxx>

#include <cometEventLoop.hxx>
#include <IMCHit.hxx>
#include <COMETGeomId.hxx>
#include <COMETGeomIdDef.hxx>
#include <IGeomInfo.hxx>
#include <IGeoField.hxx>
#include <TFile.h>
#include <TTree.h>
#include <ICDCWireManager.hxx>
#include <IFieldManager.hxx>
#include <IOADatabase.hxx>

#include <vector>
#include <iostream>
#include <string>
#include <sstream>

#include <EventDisplay.h>
#include <Track.h>
#include <FieldManager.h>

#include <TGeoGlobalMagField.h>
#include <TGeoManager.h>

class TMyEventLoop: public COMET::ICOMETEventLoopFunction {
public:
  TMyEventLoop(): 
    fFileName("default.root"), 
    fFileMode("recreate"), 
    fOutputDir("../anal"),

    fEventsNumForAnaly(150),
    fCoincidenceCount(0),
    fSingleTurnCount(0),
    fMultiTurnCount(0),
    fSaveHoughTransform(0),
    entryId(0)
  {}
  virtual ~TMyEventLoop() {}
 
  void Usage(void){
    std::cout << "-O filename=<name> Specify the output file name [default=" << fFileName << "]" << std::endl;
    std::cout << "-O directory=<name> Specify the output directory [default=" << fFileName << "]" << std::endl;
  }

  virtual bool SetOption(std::string option, std::string value=""){
    if(option== "filename") fFileName=value;
    else if(option== "directory") fOutputDir=value;
    else if(option== "filemode") fFileMode=value;
    else return false;
    return true;
  }

  virtual void Initialize(void) {
    std::cout << "Initialize" << std::endl;

    strcpy(fullName, fOutputDir.c_str());
    strcat(fullName, "/");
    strcat(fullName, fFileName.c_str());
    fOutputFile = TFile::Open(fullName, fFileMode.c_str());
    fOutputFile->cd();
    fTrdata = new TTree("trdata", "Tree Data");
        
    fTrdata->Branch("eventId", &eventId, "eventId/I");
    fTrdata->Branch("entryId", &entryId, "entryId/I");
    fTrdata->Branch("nPOCA", &nPOCA, "nPOCA/I");
    fTrdata->Branch("POCAx", POCAx, "POCAx[nPOCA]/D");
    fTrdata->Branch("POCAy", POCAy, "POCAy[nPOCA]/D");
    fTrdata->Branch("POCAz", POCAz, "POCAz[nPOCA]/D");
    fTrdata->Branch("CDCEnterX", &CDCEnterX, "CDCEnterX/D");
    fTrdata->Branch("CDCEnterY", &CDCEnterY, "CDCEnterY/D");
    fTrdata->Branch("CDCEnterZ", &CDCEnterZ, "CDCEnterZ/D");
    fTrdata->Branch("CDCEnterPx", &CDCEnterPx, "CDCEnterPx/D");
    fTrdata->Branch("CDCEnterPy", &CDCEnterPy, "CDCEnterPy/D");
    fTrdata->Branch("CDCEnterPz", &CDCEnterPz, "CDCEnterPz/D");
    fTrdata->Branch("CDCEnterPz", &CDCEnterPz, "CDCEnterPz/D");
    fTrdata->Branch("RecoPt", &RecoPt, "RecoPt/D");	
  }
  
  bool operator () (COMET::ICOMETEvent& event) {

    eventId = event.GetEventId();
    if (eventId>fEventsNumForAnaly) return true;
        
    std::cout << std::endl;
    std::cout << "*Event Id: " << eventId << std::endl;
    std::cout << std::endl;   

    COMET::IOADatabase::Get().Geometry();

    COMET::IHandle<COMET::IG4HitContainer> CTHHits = event.Get<COMET::IG4HitContainer>("truth/g4Hits/CTH");
    COMET::IHandle<COMET::IG4HitContainer> CDCHits = event.Get<COMET::IG4HitContainer>("truth/g4Hits/CDC");
    COMET::IHandle<COMET::IG4TrajectoryContainer> Trajectories = event.Get<COMET::IG4TrajectoryContainer>("truth/G4Trajectories");
    
    COMET::IHandle<COMET::IHitSelection> CDCHits_DetResp = event.Get<COMET::IHitSelection>("./hits/mcCDC");

    IMCTrigger* MCTrigger = new IMCTrigger("mctrigger", "MC trigger");
    IHoughTransform* HoughTransform = new IHoughTransform("houghtransform","Hough Transform");
    IHelixTracker* HelixTracker = new IHelixTracker("helixtracker","Helix Tracker");

    ////////////////////////////////////////////////////
    //////              Trigger                  ///////               
    ////////////////////////////////////////////////////

    bool FourFoldCoincidence;
    std::vector < std::pair < std::vector< int >, std::vector< int> > > PairCandidates;   
    
    MCTrigger->MakeCTHMap(CTHHits, Trajectories);
    MCTrigger->Process();
    FourFoldCoincidence = MCTrigger->GetFourFoldCoincidence();
    PairCandidates = MCTrigger->GetPairCandidates();
    MCTrigger->PrintResults();
    
    ////////////////////////////////////////////////////
    //////              Tracking                 ///////               
    ////////////////////////////////////////////////////

    TCanvas *c_hits;
    if (fSaveHoughTransform) c_hits = new TCanvas("c_hits", "c_hits", 1000,1000);  

    HoughTransform->LoadMCHits(CDCHits_DetResp, Trajectories, CDCHits);
    HoughTransform->PrintMCStatus();    

    if (FourFoldCoincidence==1 && HoughTransform->PreTrackCut()==1) {     

      HoughTransform->ImportTriggerInfo(PairCandidates);
      HoughTransform->Process();      
      HoughTransform->RecognizeHits();

      int nRecoHit          =HoughTransform->GetNumberOfRecognizedHits();
      bool RecoCL3          =HoughTransform->GetCL3();     
      bool RecoCL5          =HoughTransform->GetCL5();
      int Reco2DCharge      =HoughTransform->Get2DCharge();
      int RecoMaxWireLayerId=HoughTransform->GetMaxWireLayerId();    
      int MCTurnNumber      =HoughTransform->GetTurnNumber();
  
      HoughTransform->PrintResults();

      if (fSaveHoughTransform){
	HoughTransform->DrawEvent(c_hits);
	c_hits->SaveAs("./HoughTransform"+TString(Form ("%d", eventId))+".png");
      }

      std::cout << "Turn Number is: " << MCTurnNumber << std::endl;
      fCoincidenceCount++;

      //////////////////////////////
      ///      Single Turn       ///
      //////////////////////////////

      if (RecoCL5==1 && Reco2DCharge==-1 && nRecoHit>30 && MCTurnNumber==1){
	fSingleTurnCount++;		

	// Helix Tracker
	HelixTracker->LoadMCHits(CDCHits_DetResp, Trajectories, CDCHits);
	HelixTracker->LoadHitsAfterHT(CDCHits_DetResp, HoughTransform);
	HelixTracker->AddSideHitPairs(RecoMaxWireLayerId,2);	 // domain 2
	HelixTracker->AddSideHitPairs(RecoMaxWireLayerId,1);	 // domain 1

	// input POCAs to ROOT file

	std::vector <TVector3> POCAs = HelixTracker->GetPOCAs();
	nPOCA=0;
	for (Int_t i=0; i<POCAs.size(); i++){
	  POCAx[nPOCA]=(POCAs.at(i))(0);
	  POCAy[nPOCA]=(POCAs.at(i))(1);
	  POCAz[nPOCA]=(POCAs.at(i))(2);
	  nPOCA++;
	}

	// MC Truth Value
	TVector3 CDCEnterPos = HelixTracker->GetEnterPos();
	TVector3 CDCEnterMom = HelixTracker->GetEnterMom();
	CDCEnterX = CDCEnterPos(0);
	CDCEnterY = CDCEnterPos(1);
	CDCEnterZ = CDCEnterPos(2);
	CDCEnterPx = CDCEnterMom(0);
	CDCEnterPy = CDCEnterMom(1);
	CDCEnterPz = CDCEnterMom(2);
      	fTrdata->Fill();	

	entryId++;
      }  

      //////////////////////////////
      ///       Multi Turn       ///
      //////////////////////////////
 
      else if (RecoMaxWireLayerId>=4 && Reco2DCharge==-1 && RecoCL3==1 && nRecoHit>30 && MCTurnNumber>1){
	fMultiTurnCount++;
      }
    }
    
    delete MCTrigger;
    delete HoughTransform;
    if (fSaveHoughTransform) {c_hits->Close(); delete c_hits;}

    return true;
  }
  
  void Finalize(COMET::ICOMETOutput* output) {
    std::cout << "Finalize" << std::endl;
    std::cout << "Coincidence Count             : "  << fCoincidenceCount << std::endl;
    std::cout << "SingleTurn Count after Trigger & Track Cut: "  << fSingleTurnCount << std::endl;
    std::cout << "Multi-Turn Count after Trigger & Track Cut: "  << fMultiTurnCount << std::endl;
    fOutputFile->cd();    
    fTrdata->Write();
    fOutputFile->Close();
   
    return;
  }
  
private:

  TFile* fOutputFile;
  TTree* fTrdata;

  std::string fOutputDir;
  std::string fFileName;
  std::string fFileMode;  
  char fullName[100];

  int fEventsNumForAnaly;

  /*** Trigger ***/
  int fCoincidenceCount;
  int fSingleTurnCount;
  int fMultiTurnCount;

  /*** HoughTransform ***/
  bool fSaveHoughTransform;


  /*** Branch Variables ***/
  Int_t    eventId;
  Int_t    entryId;
  Int_t nPOCA;
  Double_t POCAx[1000];
  Double_t POCAy[1000];
  Double_t POCAz[1000];
  Double_t CDCEnterX;
  Double_t CDCEnterY;
  Double_t CDCEnterZ;
  Double_t CDCEnterPx;
  Double_t CDCEnterPy;
  Double_t CDCEnterPz;
  Double_t RecoPt;

};

int main(int argc, char **argv) {
  TMyEventLoop userCode;
  cometEventLoop(argc, argv, userCode);
}
