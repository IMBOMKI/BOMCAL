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
    fOutputDir("../anal")
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
  }
  
  bool operator () (COMET::ICOMETEvent& event) {
    
    int EventId = event.GetEventId();
    if (EventId>20) return true;
    
    /******** Trigger *******/
    bool FourFoldCoincidence;
    std::vector < std::pair < std::vector< int >, std::vector< int> > > PairCandidates;
    
    /******** Tracking ********/
    bool RecoCL3;
    int  Reco2DCharge;
    int  RecoMaxWireLayerId;

    ////////////////////////////////////////////////////
    //////     SimG4/SimDetResp Information      ///////               
    ////////////////////////////////////////////////////
    
    
    std::cout << std::endl;
    std::cout << "*Event Id: " << EventId << std::endl;
    std::cout << std::endl;
    
//if (!fGeoManager)   fGeoManager   = COMET::IOADatabase::Get().Geometry();
//if (!fFieldManager) fFieldManager = COMET::IFieldManager::Import(); // Open current field maps

    COMET::IOADatabase::Get().Geometry();

    COMET::IHandle<COMET::IG4HitContainer> CTHHits = event.Get<COMET::IG4HitContainer>("truth/g4Hits/CTH");
    COMET::IHandle<COMET::IG4HitContainer> CDCHits = event.Get<COMET::IG4HitContainer>("truth/g4Hits/CDC");
    COMET::IHandle<COMET::IG4TrajectoryContainer> Trajectories = event.Get<COMET::IG4TrajectoryContainer>("truth/G4Trajectories");
    
    COMET::IHandle<COMET::IHitSelection> CDCHits_DetResp = event.Get<COMET::IHitSelection>("./hits/mcCDC");

    ////////////////////////////////////////////////////
    //////         Analysis Class Object         ///////               
    ////////////////////////////////////////////////////

    IMCTrigger* MCTrigger = new IMCTrigger("mctrigger", "MC trigger");
    IHoughTransform* HoughTransform = new IHoughTransform("houghtransform","Hough Transform");
    IGenFitting* GenFitting = new IGenFitting("genfitting", "GenFitting");
    //GenFitting->ImportEnvironments(fGeoManager,fFieldManager);
    GenFitting->Init();

    ////////////////////////////////////////////////////
    //////              Trigger                  ///////               
    ////////////////////////////////////////////////////
    
    MCTrigger->MakeCTHMap(CTHHits, Trajectories);
    MCTrigger->Process();
    FourFoldCoincidence = MCTrigger->GetFourFoldCoincidence();
    PairCandidates = MCTrigger->GetPairCandidates();
    MCTrigger->PrintResults();
    
    ////////////////////////////////////////////////////
    //////              Tracking                 ///////               
    ////////////////////////////////////////////////////

    TCanvas *c_hits = new TCanvas("c_hits", "c_hits", 1000,1000);  
    HoughTransform->LoadMCHits(CDCHits_DetResp, Trajectories, CDCHits);
    HoughTransform->PrintMCStatus();    

    if (FourFoldCoincidence==1 && HoughTransform->PreTrackCut()==1) {     

      HoughTransform->ImportTriggerInfo(PairCandidates);
      HoughTransform->Process();      
      HoughTransform->RecognizeHits();

      RecoCL3           =HoughTransform->GetCL3();
      Reco2DCharge      =HoughTransform->Get2DCharge();
      RecoMaxWireLayerId=HoughTransform->GetMaxWireLayerId();      
      HoughTransform->PrintResults();
      //HoughTransform->DrawEvent(c_hits);
      //c_hits->SaveAs("./EventId"+TString(Form ("%d", EventId))+".png");
      fCoincidenceCount++;

      
      if (RecoMaxWireLayerId>=4 && Reco2DCharge==-1 && RecoCL3==1){
	//GenFitting->LoadHitsAfterHT(CDCHits_DetResp, HoughTransform);
	GenFitting->LoadMCHits(CDCHits_DetResp, Trajectories, CDCHits);
	if (!GenFitting->doFit()) {
	  std::cout << "Fitting Failed" << std::endl; 
	  return true;
	}
	GenFitting->AddEvent(display);
      }      
    }
    
    delete MCTrigger;
    delete HoughTransform;
    delete GenFitting;

    delete c_hits;
    return true;
  }
  
  void Finalize(COMET::ICOMETOutput* output) {
    std::cout << "Finalize" << std::endl;
    std::cout << "Coincidence Count: " << fCoincidenceCount << std::endl;
    display->open();
    return;
  }
  
private:

  std::string fOutputDir;
  std::string fFileName;
  std::string fFileMode;  

  int fCoincidenceCount=0;

  COMET::IFieldManager* fFieldManager;
  TGeoManager *fGeoManager ;
  genfit::EventDisplay* display;
};

int main(int argc, char **argv) {
  TMyEventLoop userCode;
  cometEventLoop(argc, argv, userCode);
}
