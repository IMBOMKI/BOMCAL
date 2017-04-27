#ifndef TMCTrigger_hxx_seen
#define TMCTrigger_hxx_seen

#include <vector>
#include <map>

#include <TObject.h>
#include <TVector3.h>
#include <TRandom.h>

#include <IAlgorithm.hxx>
#include <ICOMETLog.hxx>

#include <IG4Trajectory.hxx>
#include <IG4VHit.hxx>
#include <IHandle.hxx>
#include <ICOMETEvent.hxx>

#include "TGeoManager.h"
#include "TGeoNode.h"



class IMCTrigger{
private:

  // Low level Information
  std::map <int, double> fTime; 
  std::map <int, int>    fIndex;      // Segment Index: 1-64
  std::map <int, int>    fScint;      // 1 for Scint, 0 for Cherenkov
  std::map <int, int>    fModule;     // 1 for Down, 0 for Up

  // High level Information (Left/Right vector for Scint/Cherenkov)
  std::vector < std::pair< std::vector< int >, std::vector< int > > > fDSTimeCluster; // int->key Value
  std::vector < std::pair< std::vector< int >, std::vector< int > > > fUSTimeCluster;

  // Pair Candidates
  std::vector < std::pair < std::vector< int >, std::vector< int> > > fPairCandidates_key; // int->key Value
  std::vector < std::pair < std::vector< int >, std::vector< int> > > fPairCandidates;     // int->Index Value

  int fCTHSegNum;
  double fStartT;
  double fEndT;
  int fShift;
  
public:
  IMCTrigger(const char*name, const char* title);
  ~IMCTrigger();
  
  /// called at the begin of run or else (should not be in event-by-event)
  int  Init();

  /// Get Volume node
  TGeoNode* GetNode(TVector3 position);

  /// Calculate number of photons in Cherenkov detector
  int GetPhotonNumber(int pdgNum, double ene, double trLen);

  // Make map for CTH hits' index & time
  void MakeCTHMap(COMET::IHandle<COMET::IG4HitContainer> & cthhits, COMET::IHandle<COMET::IG4TrajectoryContainer> trajectories);

  void SetMCTriggerVariable(int shift){ fShift = shift; }
  void MakeTimeCluster(int Module);
  void PrintTimeCluster();
  void ApplyShiftCondition(int Module, int shift);
  void PrintResults();
  bool GetFourFoldCoincidence();
  std::vector < std::pair < std::vector< int >, std::vector< int> > > GetPairCandidates();
  void Process();

  void Clear();

  /// called at the end of run or else (should not be in event-by-event)
  int  Finish();

};
#endif
