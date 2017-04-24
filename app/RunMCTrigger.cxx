#include <cometEventLoop.hxx>
#include <IMCHit.hxx>
#include <COMETGeomId.hxx>
#include <COMETGeomIdDef.hxx>
#include <IGeomInfo.hxx>
#include <IMCTrigger.hxx>
#include <TFile.h>
#include <TTree.h>

#include <vector>
#include <iostream>
#include <sstream>

class TMyEventLoop: public COMET::ICOMETEventLoopFunction {
public:
    TMyEventLoop() {}
    virtual ~TMyEventLoop() {}

    virtual void Initialize(void) {
      std::cout << "Initialize" << std::endl;
      fMCTrigger = new IMCTrigger("mctrigger", "MC trigger");
      fMCTrigger->Init();
    }

    bool operator () (COMET::ICOMETEvent& event) {
      eventNumber = event.GetEventId();
      COMET::IOADatabase::Get().Geometry();
      COMET::IHandle<COMET::IG4HitContainer> hits = event.Get<COMET::IG4HitContainer>("truth/g4Hits/CTH");
      return true;
    }

    void Finalize(COMET::ICOMETOutput* output) {
      std::cout << "Finalize" << std::endl;
      delete fMCTrigger;
      return;
    }

private:
    IMCTrigger* fMCTrigger;
    int eventNumber;
};

int main(int argc, char **argv) {
    TMyEventLoop userCode;
    cometEventLoop(argc, argv, userCode);
}
