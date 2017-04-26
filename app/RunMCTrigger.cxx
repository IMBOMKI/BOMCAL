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

#include <vector>
#include <iostream>
#include <sstream>

class TMyEventLoop: public COMET::ICOMETEventLoopFunction {
public:
  TMyEventLoop(): fileName("default.root"), fileMode("recreate"), outputDir("../anal") {}
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
    fMCTrigger = new IMCTrigger("mctrigger", "MC trigger");
    fMCTrigger->Init();

    fHoughTransform = new IHoughTransform("houghtransform","Hough Transform");
    fHoughTransform->Init();

  }
  
  bool operator () (COMET::ICOMETEvent& event) {
    eventNumber = event.GetEventId();
    COMET::IOADatabase::Get().Geometry();
    COMET::IHandle<COMET::IG4HitContainer> CTHHits = event.Get<COMET::IG4HitContainer>("truth/g4Hits/CTH");
    COMET::IHandle<COMET::IG4TrajectoryContainer> CTHTrajectories = event.Get<COMET::IG4TrajectoryContainer>("truth/G4Trajectories");
    
    COMET::IHandle<COMET::IHitSelection> CDCHits_DetResp = event.Get<COMET::IHitSelection>("./hits/mcCDC");


    ////////////////////////////////////////////////////
    //////              Trigger                  ///////               
    ////////////////////////////////////////////////////

    fMCTrigger->MakeCTHMap(CTHHits, CTHTrajectories);
    fMCTrigger->Process(1);
    fFourFoldCoincidence = fMCTrigger->GetFourFoldCoincidence();
    fPairCandidates = fMCTrigger->GetPairCandidates();
    fMCTrigger->PrintPairCandidates();
    fMCTrigger->Clear();

    ////////////////////////////////////////////////////
    //////              Tracking                 ///////               
    ////////////////////////////////////////////////////
    
    if (fFourFoldCoincidence==1) {
      fHoughTransform->LoadMCHits(CDCHits_DetResp);
      FourFoldCount++;
    }
        
    return true;
  }
  
  void Finalize(COMET::ICOMETOutput* output) {
    std::cout << "Finalize" << std::endl;
    std::cout << "FourFold Count: " << FourFoldCount << std::endl;
    delete fMCTrigger;
    return;
  }
  
private:

  std::string outputDir;
  std::string fileName;
  std::string fileMode;
  char fullName[100];

  int eventNumber;

  /********  Trigger  *******/
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
