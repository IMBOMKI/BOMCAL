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
    fileName("default.root"), 
    fileMode("recreate"), 
    outputDir("../anal")
  {}
  virtual ~TMyEventLoop() {}
 
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

  virtual void Initialize(void) {
    std::cout << "Initialize" << std::endl;
    //fMCTrigger->Init();

    //fHoughTransform->Init();
  }
  
  bool operator () (COMET::ICOMETEvent& event) {

    ////////////////////////////////////////////////////
    //////     SimG4/SimDetResp Information      ///////               
    ////////////////////////////////////////////////////
    
    fEventId = event.GetEventId();
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
    fMCTrigger->SetMCTriggerVariable(1); // Set shift tolerance
    fMCTrigger->Process();
    fFourFoldCoincidence = fMCTrigger->GetFourFoldCoincidence();
    fPairCandidates = fMCTrigger->GetPairCandidates();
    fMCTrigger->PrintResults();
    
    ////////////////////////////////////////////////////
    //////              Tracking                 ///////               
    ////////////////////////////////////////////////////
    
    if (fFourFoldCoincidence==1) {
      //fHoughTransform->LoadMCHits(CDCHits_DetResp, Trajectories);
      fHoughTransform->SetHoughTransformVariables(3,100,300.0,0.02,-0.02,8);
      FourFoldCount++;
    }
    

    delete fMCTrigger;
    delete fHoughTransform;

    return true;

  }
  
  void Finalize(COMET::ICOMETOutput* output) {
    std::cout << "Finalize" << std::endl;
    std::cout << "FourFold Count: " << FourFoldCount << std::endl;
    return;
  }
  
private:

  std::string outputDir;
  std::string fileName;
  std::string fileMode;
  char fullName[100];

  int fEventId;

  /******** Trigger *******/
  IMCTrigger* fMCTrigger;
  bool fFourFoldCoincidence;
  std::vector < std::pair < std::vector< int >, std::vector< int> > > fPairCandidates;
  int FourFoldCount;
 
  /******** Tracking ********/
  IHoughTransform* fHoughTransform;

};

int main(int argc, char **argv) {
  TMyEventLoop userCode;
  cometEventLoop(argc, argv, userCode);
}
