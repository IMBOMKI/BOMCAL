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

#include "time.h"
#include "sys/time.h"

class TMyEventLoop: public COMET::ICOMETEventLoopFunction {
public:
  TMyEventLoop(): 
    fFileName("default.root"), 
    fFileMode("recreate"), 
    fOutputDir("../anal"),

    fEventsNumForAnaly(1000),
    fCoincidenceCount(0),
    fSingleTurnCount(0),
    fMultiTurnCount(0),
    fSaveHoughTransform(0),
    fMaxGenFitTry(5)
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

    fTrdata->Branch("fTriggerActivated", &fTriggerActivated, "fTriggerActivated/O");
    fTrdata->Branch("fnRecoHit", &fnRecoHit, "fnRecoHit/I");
    fTrdata->Branch("fRecoCL3", &fRecoCL3, "fRecoCL3/O");
    fTrdata->Branch("fReco2DCharge", &fReco2DCharge, "fReco2DCharge/I");
    fTrdata->Branch("fRecoMaxWireLayerId", &fRecoMaxWireLayerId, "fRecoMaxWireLayerId/I");
    fTrdata->Branch("fMCTurnNumber", &fMCTurnNumber, "fMCTurnNumber/I");

    fTrdata->Branch("fFitSuccess", &fFitSuccess, "fFitSuccess/O");
    fTrdata->Branch("fpFit", &fpFit, "fpFit/D");
    fTrdata->Branch("fpxFit", &fpxFit, "fpxFit/D");
    fTrdata->Branch("fpyFit", &fpyFit, "fpyFit/D");
    fTrdata->Branch("fpzFit", &fpzFit, "fpzFit/D");
    fTrdata->Branch("fxFit", &fxFit, "fxFit/D");
    fTrdata->Branch("fyFit", &fyFit, "fyFit/D");
    fTrdata->Branch("fzFit", &fzFit, "fzFit/D");

    fTrdata->Branch("fChi2", &fChi2, "fChi2/D");
    fTrdata->Branch("fNdf", &fNdf, "fNdf/I");
    fTrdata->Branch("fChi2Ndf", &fChi2Ndf, "fChi2Ndf/D");
    fTrdata->Branch("fpEnter", &fpEnter, "fpEnter/D");
    fTrdata->Branch("fpIni", &fpIni, "fpIni/D");
    fTrdata->Branch("fpT_HT", &fpT_HT, "fpT_HT/D");
    fTrdata->Branch("fpT_Reseeded", &fpT_Reseeded, "fpT_Reseeded/D");
    fTrdata->Branch("fpT_Truth",&fpT_Truth, "fpT_Truth/D");

    fTrdata->Branch("pVal",&pVal, "pVal/D");
    fTrdata->Branch("hqopPu",&hqopPu, "hqopPu/D");
    fTrdata->Branch("hupPu",&hupPu, "hupPu/D");
    fTrdata->Branch("hvpPu",&hvpPu, "hvpPu/D");
    fTrdata->Branch("huPu",&huPu, "huPu/D");
    fTrdata->Branch("hvPu",&hvPu, "hvPu/D");

    // for Testing purpose

    // MC Initial Value
    fTrdata->Branch("fCDCEnterX", &fCDCEnterX, "fCDCEnterX/D");
    fTrdata->Branch("fCDCEnterY", &fCDCEnterY, "fCDCEnterY/D");
    fTrdata->Branch("fCDCEnterZ", &fCDCEnterZ, "fCDCEnterZ/D");
    fTrdata->Branch("fCDCEnterPx", &fCDCEnterPx, "fCDCEnterPx/D");
    fTrdata->Branch("fCDCEnterPy", &fCDCEnterPy, "fCDCEnterPy/D");
    fTrdata->Branch("fCDCEnterPz", &fCDCEnterPz, "fCDCEnterPz/D");

    // Hough Transform Transverse Mom & Pos Initial Seed
    fTrdata->Branch("fFitEnterY",  &fFitEnterY, "fFitEnterY/D");
    fTrdata->Branch("fFitEnterZ",  &fFitEnterZ, "fFitEnterZ/D");
    fTrdata->Branch("fFitEnterPy", &fFitEnterPy, "fFitEnterPy/D");
    fTrdata->Branch("fFitEnterPz", &fFitEnterPz, "fFitEnterPz/D");

    fTrdata->Branch("fAnalysisTime_HT", &fAnalysisTime_HT, "fAnalysisTime_HT/D");
    fTrdata->Branch("fAnalysisTime_GENFIT", &fAnalysisTime_GENFIT, "fAnalysisTime_GENFIT/D");
  }
  
  bool operator () (COMET::ICOMETEvent& event) {
    
    // Nullifying ////////////////////////
    fTriggerActivated=0;
    fnRecoHit=0;       
    fRecoCL3=0;        
    fReco2DCharge=0;
    fRecoMaxWireLayerId=0;
    fMCTurnNumber=0;    
    fFitSuccess=0;
    fpFit=0;
    fChi2=0;
    fNdf=0;
    fChi2Ndf=0;
    fpEnter=0;
    fpIni=0;
    fpT_HT=0;
    fpT_Reseeded=0;
    fpT_Truth=0;

    fAnalysisTime_HT=0;
    fAnalysisTime_GENFIT=0;
    struct timespec HTstart, HTend;
    struct timespec GENFITstart, GENFITend;
    /////////////////////////////////////

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
    IGenFitting* GenFitting = new IGenFitting("genfitting", "GenFitting");

    ////////////////////////////////////////////////////
    //////              Trigger                  ///////               
    ////////////////////////////////////////////////////

    //bool FourFoldCoincidence;
    std::vector < std::pair < std::vector< int >, std::vector< int> > > PairCandidates;   
    
    MCTrigger->MakeCTHMap(CTHHits, Trajectories);
    MCTrigger->Process();
    fTriggerActivated = MCTrigger->GetFourFoldCoincidence();
    PairCandidates = MCTrigger->GetPairCandidates();
    MCTrigger->PrintResults();
    
    ////////////////////////////////////////////////////
    //////              Tracking                 ///////               
    ////////////////////////////////////////////////////

    /////////////////////
    //                 //
    // Hough Transform // 
    //                 //
    /////////////////////

    HoughTransform->LoadMCHits(CDCHits_DetResp, Trajectories, CDCHits);
    HoughTransform->PrintMCStatus();    

    if (fTriggerActivated==1 && HoughTransform->PreTrackCut()==1) {     

      // Start Time measureing && Event Drawing //////////////////////////////////////
      clock_gettime(CLOCK_MONOTONIC_RAW, &HTstart); 
      TCanvas *c_hits;
      if (fSaveHoughTransform) c_hits = new TCanvas("c_hits", "c_hits", 1000,1000);  
      ////////////////////////////////////////////////////////////////////////////////

      HoughTransform->ImportTriggerInfo(PairCandidates);
      HoughTransform->Process();      
      HoughTransform->RecognizeHits();
      HoughTransform->TuneRadiusWithPOCAs();

      // Extract Transverse Position & Momentum ////////////////////////////
      fpT_HT       = HoughTransform->GetpT_HT();
      fpT_Truth    = HoughTransform->GetpT_Truth();
      fpT_Reseeded = HoughTransform->GetpT_Reseeded();
      TVector3 fCDCEnterPos = HoughTransform->GetEnterPos_Truth();
      TVector3 fCDCEnterMom = HoughTransform->GetEnterMom_Truth();
      fCDCEnterX  = fCDCEnterPos(0);
      fCDCEnterY  = fCDCEnterPos(1);
      fCDCEnterZ  = fCDCEnterPos(2);
      fCDCEnterPx = fCDCEnterMom(0);
      fCDCEnterPy = fCDCEnterMom(1);
      fCDCEnterPz = fCDCEnterMom(2);

      TVector3 fFitEnterPos = HoughTransform->GetEnterXYPair_Reseeded().second;
      TVector3 fFitEnterMom = HoughTransform->GetEnterPxPyPair_Reseeded().second;
      if(!COMET::IGeomInfo::DetectorSolenoid().GetDetPositionInGlobalCoordinate(fFitEnterPos, fFitEnterPos)){std::cout << "Coordinate change fails (Local to Master)" << std::endl; return true;}
      //fFitEnterX  = fFitEnterPos(0);
      fFitEnterY  = fFitEnterPos(1);
      fFitEnterZ  = fFitEnterPos(2);
      fFitEnterPz = -fFitEnterMom(0);
      fFitEnterPy = fFitEnterMom(1);
      //fFitEnterPx = fFitEnterMom(2);
      //////////////////////////////////////////////////////////////////////

      // Extract Track Cut /////////////////////////////////////////////////
      fnRecoHit          =HoughTransform->GetNumberOfRecognizedHits();
      fRecoCL3          =HoughTransform->GetCL3();     
      fReco2DCharge      =HoughTransform->Get2DCharge();
      fRecoMaxWireLayerId=HoughTransform->GetMaxWireLayerId();    
      fMCTurnNumber      =HoughTransform->GetTurnNumber();
      //////////////////////////////////////////////////////////////////////


      HoughTransform->PrintResults();

      if (fSaveHoughTransform){
	HoughTransform->DrawEvent(c_hits);
	c_hits->SaveAs("./EventId"+TString(Form ("%d", eventId))+".png");
	c_hits->Close(); 
	delete c_hits;
      }
      fCoincidenceCount++;
      clock_gettime(CLOCK_MONOTONIC_RAW, &HTend); 
      fAnalysisTime_HT = (HTend.tv_sec - HTstart.tv_sec)*1000000 + (HTend.tv_nsec - HTstart.tv_nsec) / 1000; 

      //////////////// 
      //            //
      // GENFITTING // 
      //            //
      ////////////////

      // Single Turn
      if (fRecoMaxWireLayerId>=4 && fReco2DCharge==-1 && fRecoCL3==1 && fnRecoHit>30 && fMCTurnNumber==1){
	fSingleTurnCount++;		
	GenFitting->Init();
	GenFitting->LoadMCHits(CDCHits_DetResp, Trajectories, CDCHits);       
	//GenFitting->ShuffleMCHits(); // This makes fitting suck...
	GenFitting->LoadHitsAfterHT(CDCHits_DetResp, HoughTransform);
	Int_t nTry=0;	
	
	clock_gettime(CLOCK_MONOTONIC_RAW, &GENFITstart); 

	if (!GenFitting->DoFit(HoughTransform)) {
	  fFitSuccess=0;

	  clock_gettime(CLOCK_MONOTONIC_RAW, &GENFITend); 
	  fAnalysisTime_GENFIT = (GENFITend.tv_sec - GENFITstart.tv_sec)*1000000 + (GENFITend.tv_nsec - GENFITstart.tv_nsec) / 1000;

	  fTrdata->Fill();
	  return true;
	}
	fFitSuccess=1;	
	GenFitting->AddEvent(display);
	fpFit=GenFitting->GetFittedMom();
	TVector3 pFit_vec = GenFitting->GetFittedpVec();
	TVector3 xFit_vec = GenFitting->GetFittedxVec();	
	fpxFit = pFit_vec(0);  fpyFit = pFit_vec(1);   fpzFit = pFit_vec(2);
	fxFit  = xFit_vec(0);  fyFit  = xFit_vec(1);   fzFit  = xFit_vec(2);

	fChi2=GenFitting->GetChi2();
	fNdf =GenFitting->GetNdf();
	fChi2Ndf=GenFitting->GetChi2Ndf();
	pVal    =GenFitting->GetpVal();
	Double_t *pullVal=GenFitting->GetPullValue();
	hqopPu = pullVal[0];
	hupPu  = pullVal[1];
        hvpPu  = *(pullVal+2);
	huPu   = *(pullVal+3);
	hvPu   = *(pullVal+4);

	std::cout << "hqopPu: " << hqopPu << std::endl;
	std::cout << "hupPu: " << hupPu << std::endl;
	std::cout << "hvpPu: " << hvpPu << std::endl;
	std::cout << "huPu: " << huPu << std::endl;
	std::cout << "hvPu: " << hvPu << std::endl;

	fpIni   =GenFitting->GetInitialMom();
	fpEnter =GenFitting->GetCDCEntranceMom();	

	clock_gettime(CLOCK_MONOTONIC_RAW, &GENFITend); 
	fAnalysisTime_GENFIT = (GENFITend.tv_sec - GENFITstart.tv_sec)*1000000 + (GENFITend.tv_nsec - GENFITstart.tv_nsec) / 1000;

	fTrdata->Fill();
	return true;
      }   

      // Multi Turn
      else if (fRecoMaxWireLayerId>=4 && fReco2DCharge==-1 && fRecoCL3==1 && fnRecoHit>30 && fMCTurnNumber>1){
	fMultiTurnCount++;
	fTrdata->Fill();
	return true;
      }
    }
    
    fTrdata->Fill();
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
  Int_t fMaxGenFitTry;
  //TVector3 fFittedMom;

  /*** Branch Variables ***/
  Int_t    eventId;

  Int_t    fTriggerActivated;
  Int_t    fnRecoHit;       
  Bool_t   fRecoCL3;        
  Int_t    fReco2DCharge;
  Int_t    fRecoMaxWireLayerId;
  Int_t    fMCTurnNumber;

  Bool_t   fFitSuccess;
  Double_t fpFit;
  Double_t fpxFit;
  Double_t fpyFit;
  Double_t fpzFit;
  Double_t fxFit;
  Double_t fyFit;
  Double_t fzFit;

  Double_t fChi2;
  Int_t    fNdf;
  Double_t fChi2Ndf;
  Double_t fpEnter;
  Double_t fpIni;
  Double_t fpT_HT;
  Double_t fpT_Reseeded;
  Double_t fpT_Truth;

  Double_t hqopPu;
  Double_t pVal;
  Double_t hupPu;
  Double_t hvpPu;
  Double_t huPu;
  Double_t hvPu;

  Double_t fAnalysisTime_HT;     // Time in mircosecond unit
  Double_t fAnalysisTime_GENFIT;

  // For testing purpose
  Double_t fCDCEnterX;
  Double_t fCDCEnterY;
  Double_t fCDCEnterZ;
  Double_t fCDCEnterPx;
  Double_t fCDCEnterPy;
  Double_t fCDCEnterPz;
  Double_t fFitEnterX;
  Double_t fFitEnterY;
  Double_t fFitEnterZ;
  Double_t fFitEnterPx;
  Double_t fFitEnterPy;
  Double_t fFitEnterPz;


};

int main(int argc, char **argv) {
  TMyEventLoop userCode;
  cometEventLoop(argc, argv, userCode);
}
