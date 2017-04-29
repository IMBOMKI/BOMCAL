#include <IMCTrigger.hxx>
#include <ITracking.hxx>
#include <IHoughTransform.hxx>

#include <cometEventLoop.hxx>
#include <IMCHit.hxx>
#include <COMETGeomId.hxx>
#include <COMETGeomIdDef.hxx>
#include <IGeomInfo.hxx>
#include <TFile.h>
#include <TTree.h>
#include <ICDCWireManager.hxx>

#include <vector>
#include <iostream>
#include <sstream>

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
  }
  
  bool operator () (COMET::ICOMETEvent& event) {

    

    ////////////////////////////////////////////////////
    //////     SimG4/SimDetResp Information      ///////               
    ////////////////////////////////////////////////////
    
    fEventId = event.GetEventId();

    std::cout << "*Event Id: " << fEventId << std::endl;
    std::cout <<  std::endl;

    COMET::IOADatabase::Get().Geometry();
    COMET::IHandle<COMET::IG4HitContainer> CTHHits = event.Get<COMET::IG4HitContainer>("truth/g4Hits/CTH");
    COMET::IHandle<COMET::IG4HitContainer> CDCHits = event.Get<COMET::IG4HitContainer>("truth/g4Hits/CDC");
    COMET::IHandle<COMET::IG4TrajectoryContainer> Trajectories = event.Get<COMET::IG4TrajectoryContainer>("truth/G4Trajectories");
    
    COMET::IHandle<COMET::IHitSelection> CDCHits_DetResp = event.Get<COMET::IHitSelection>("./hits/mcCDC");

    ////////////////////////////////////////////////////
    //////         Analysis Class Object         ///////               
    ////////////////////////////////////////////////////

    fMCTrigger = new IMCTrigger("mctrigger", "MC trigger");
    fHoughTransform = new IHoughTransform("houghtransform","Hough Transform");

    ////////////////////////////////////////////////////
    //////              Trigger                  ///////               
    ////////////////////////////////////////////////////
    
    fMCTrigger->MakeCTHMap(CTHHits, Trajectories);
    //fMCTrigger->SetMCTriggerVariable(1); // Set shift tolerance
    fMCTrigger->Process();
    fFourFoldCoincidence = fMCTrigger->GetFourFoldCoincidence();
    fPairCandidates = fMCTrigger->GetPairCandidates();
    fMCTrigger->PrintResults();
    
    ////////////////////////////////////////////////////
    //////              Tracking                 ///////               
    ////////////////////////////////////////////////////

    fHoughTransform->LoadMCHits(CDCHits_DetResp, Trajectories);
    fHoughTransform->PrintMCStatus();    

    if (fFourFoldCoincidence==1 && fHoughTransform->PreTrackCut()==1) {     
      //fHoughTransform->SetHoughTransformVariables(3,100,300.0,0.02,-0.02,8,5,0,0);
      fHoughTransform->Process();      
      fCoincidenceCount++;
    }
    

    delete fMCTrigger;
    delete fHoughTransform;
    return true;
  }
  
  void Finalize(COMET::ICOMETOutput* output) {
    std::cout << "Finalize" << std::endl;
    std::cout << "Coincidence Count: " << fCoincidenceCount << std::endl;
    return;
  }
  
private:

  std::string fOutputDir;
  std::string fFileName;
  std::string fFileMode;


  int fEventId;

  /******** Trigger *******/
  IMCTrigger* fMCTrigger;
  bool fFourFoldCoincidence;
  std::vector < std::pair < std::vector< int >, std::vector< int> > > fPairCandidates;
  int fCoincidenceCount;
 
  /******** Tracking ********/
  IHoughTransform* fHoughTransform;
  
};

int main(int argc, char **argv) {
  TMyEventLoop userCode;
  cometEventLoop(argc, argv, userCode);
}
