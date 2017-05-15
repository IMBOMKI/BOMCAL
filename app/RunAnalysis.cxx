#include <IMCTrigger.hxx>
#include <ITracking.hxx>
#include <IHoughTransform.hxx>
#include <IGenFitting.hxx>

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

    fEventsNumForAnaly(10000),
    fCoincidenceCount(0),
    fSingleTurnCount(0),
    fMultiTurnCount(0),
    fSaveHoughTransform(0)
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
    display = genfit::EventDisplay::getInstance();


    strcpy(fullName, fOutputDir.c_str());
    strcat(fullName, "/");
    strcat(fullName, fFileName.c_str());
    fOutputFile = TFile::Open(fullName, fFileMode.c_str());
    fOutputFile->cd();
    fTrdata = new TTree("trdata", "Tree Data");
        
    fTrdata->Branch("eventId", &eventId, "eventId/I");
    fTrdata->Branch("fpFit", &fpFit, "fpFit/D");
    fTrdata->Branch("fChi2", &fChi2, "fChi2/D");
    fTrdata->Branch("fNdf", &fNdf, "fNdf/I");
    fTrdata->Branch("fChi2Ndf", &fChi2Ndf, "fChi2Ndf/D");
    fTrdata->Branch("fpEnter", &fpEnter, "fpEnter/D");
    fTrdata->Branch("fpIni", &fpIni, "fpIni/D");
  }
  
  bool operator () (COMET::ICOMETEvent& event) {
    
    eventId = event.GetEventId();
    if (eventId>fEventsNumForAnaly) return true;
        
    std::cout << std::endl;
    std::cout << "*Event Id: " << eventId << std::endl;
    std::cout << std::endl;
    
//if (!fGeoManager)   fGeoManager   = COMET::IOADatabase::Get().Geometry();
//if (!fFieldManager) fFieldManager = COMET::IFieldManager::Import(); // Open current field maps

    COMET::IOADatabase::Get().Geometry();

    COMET::IHandle<COMET::IG4HitContainer> CTHHits = event.Get<COMET::IG4HitContainer>("truth/g4Hits/CTH");
    COMET::IHandle<COMET::IG4HitContainer> CDCHits = event.Get<COMET::IG4HitContainer>("truth/g4Hits/CDC");
    COMET::IHandle<COMET::IG4TrajectoryContainer> Trajectories = event.Get<COMET::IG4TrajectoryContainer>("truth/G4Trajectories");
    
    COMET::IHandle<COMET::IHitSelection> CDCHits_DetResp = event.Get<COMET::IHitSelection>("./hits/mcCDC");

    IMCTrigger* MCTrigger = new IMCTrigger("mctrigger", "MC trigger");
    IHoughTransform* HoughTransform = new IHoughTransform("houghtransform","Hough Transform");
    IGenFitting* GenFitting = new IGenFitting("genfitting", "GenFitting");

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
      int Reco2DCharge      =HoughTransform->Get2DCharge();
      int RecoMaxWireLayerId=HoughTransform->GetMaxWireLayerId();    
      int MCTurnNumber      =HoughTransform->GetTurnNumber();
  
      HoughTransform->PrintResults();

      if (fSaveHoughTransform){
	HoughTransform->DrawEvent(c_hits);
	c_hits->SaveAs("./EventId"+TString(Form ("%d", eventId))+".png");
      }

      std::cout << "Turn Number is: " << MCTurnNumber << std::endl;
      fCoincidenceCount++;
      
      if (RecoMaxWireLayerId>=4 && Reco2DCharge==-1 && RecoCL3==1 && nRecoHit>30 && MCTurnNumber==1){
	fSingleTurnCount++;		

	GenFitting->Init();
	GenFitting->LoadMCHits(CDCHits_DetResp, Trajectories, CDCHits);       
	//GenFitting->ShuffleMCHits();
	GenFitting->LoadHitsAfterHT(CDCHits_DetResp, HoughTransform);
	
	if (!GenFitting->DoFit()) {
	  return true;
	}
	GenFitting->AddEvent(display);
	fpFit=GenFitting->GetFittedMom(); fpFit *= 1e3;
	fChi2=GenFitting->GetChi2();
	fNdf =GenFitting->GetNdf();
	fChi2Ndf=GenFitting->GetChi2Ndf();
	fpIni   =GenFitting->GetInitialMom();
	fpEnter =GenFitting->GetCDCEntranceMom();

	fTrdata->Fill();	
      }   
      else if (RecoMaxWireLayerId>=4 && Reco2DCharge==-1 && RecoCL3==1 && nRecoHit>30 && MCTurnNumber>1){
	fMultiTurnCount++;
      }
    }
    
    delete MCTrigger;
    delete HoughTransform;
    if (fSaveHoughTransform) {c_hits->Close(); delete c_hits;}

    delete GenFitting;

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
    display->open();
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

  /*** GenFitting ***/
  genfit::EventDisplay* display;
  //TVector3 fFittedMom;

  /*** Branch Variables ***/
  Int_t    eventId;
  Double_t fpFit;
  Double_t fChi2;
  Int_t    fNdf;
  Double_t fChi2Ndf;

  Double_t fpEnter;
  Double_t fpIni;

};

int main(int argc, char **argv) {
  TMyEventLoop userCode;
  cometEventLoop(argc, argv, userCode);
}
