#ifndef TGenFitting_hxx_seen
#define TGenFitting_hxx_seen

#include <ITracking.hxx>
#include <IHoughTransform.hxx>

#include <vector>
#include <map>

#include <TVirtualFitter.h>
#include <TObject.h>
#include <TVector3.h>
#include <TRandom.h>
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TFile.h"

#include <ICOMETLog.hxx>

#include <IG4Trajectory.hxx>
#include <IG4VHit.hxx>
#include <IHandle.hxx>
#include <ICOMETEvent.hxx>
#include <IReconTrack.hxx>
#include <IReconTrackCand.hxx>
#include <IFieldManager.hxx>


////// GenFit
#include <Track.h>
#include <EventDisplay.h>

class IGenFitting: public ITracking {
public:
  IGenFitting(const char*name, const char* title);
  virtual ~IGenFitting();

  /// called at the begin of run or else (should not be in event-by-event)
  int  Init();

  int BeginOfEvent();

  int EndOfEvent();

  void ImportEnvironments(TGeoManager* geom, COMET::IFieldManager* field){
    fGeoManager   = geom;
    fFieldManager = field;
  }

  //TVector3 SmearSeed(TVector3 seed, Double_t smear);

  void LoadMCHits(COMET::IHandle<COMET::IHitSelection> hitHandle, COMET::IHandle<COMET::IG4TrajectoryContainer> trajectories,COMET::IHandle<COMET::IG4HitContainer> cdcHits){
    ITracking::LoadMCHits(hitHandle, trajectories, cdcHits);
  }  
  void ShuffleMCHits(){ ITracking::ShuffleMCHits(); }

  void LoadHitsAfterHT(COMET::IHandle<COMET::IHitSelection> hitHandle, IHoughTransform* hough);
  void AddEvent(genfit::EventDisplay* display) {  display->addEvent(fitTrack); } 
  
  int DoFit(IHoughTransform *hough);

  double GetFittedMom() { return fpFit; }
  double GetChi2() { return fChi2; }
  int GetNdf()  { return fNdf;  }
  double GetChi2Ndf() { return fChi2Ndf; }
  double GetInitialMom() { return sqrt(pow(fGenTrPx,2)+pow(fGenTrPy,2)+pow(fGenTrPz,2)); }
  double GetCDCEntranceMom() { return sqrt(pow(fCDCEnterPx,2)+pow(fCDCEnterPy,2)+pow(fCDCEnterPz,2)); }



private:
  TTree   *fTree;
  Bool_t fUseBetheBloch; /// flag to turn on/off the BetheBloch  
  Bool_t fUseBrems;      /// flag to turn on/off the Brems  
  std::string fMethod; /// fitting method  
  Int_t  fPID;           /// Input Particle ID in PDGEncoding to fit the track  
  UInt_t fMinIterations; /// minimum number of iterations  
  UInt_t fMaxIterations; /// maximum number of iterations  
  UInt_t fMinHitsInTrack;   /// minimum number of hits in track  
  Int_t  fMinNDF;           /// minimum value for NDF  
  Double_t fMaxMomentum;   /// maximum momentum [MeV]  
  Double_t fMinMomentum;   /// minimum momentum [MeV]  
  Double_t fMaxMomDiff;    /// maximum allowed difference between fitted and initial momenta [MeV]  
  Bool_t   fUseExtGeomFile;/// flag to use geometry in external file  
  std::string fGeometry; /// name of the geometry  
  Bool_t   fUseExtFieldFile; /// flag to use field maps in external file  
  std::string fFieldMap; /// name of the fieldmap  
  Bool_t   fUseMCTruth;    /// Flag to use MC true hit position  
  Bool_t   fSmearing;      /// enable Gaussian smearing for the drift distance using given fSigmaD  
  Double_t fSigmaD;        /// position resolution  
  Double_t fSigmaWP;        /// Wire Position resolution
  Bool_t   fSaveHistogram; /// Flage to save Hitogram
  Bool_t   fUseTransverseSeed; /// flag to use field maps in external file  

  TVector3 posInit;
  TVector3 momInit;
  genfit::Track* fitTrack;
  TGeoManager*   fGeoManager;
  COMET::IFieldManager* fFieldManager;
  double fpFit;
  double fChi2;
  int    fNdf;
  double fChi2Ndf;

};
#endif
