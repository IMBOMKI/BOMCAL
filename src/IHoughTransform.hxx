#ifndef THoughTransform_hxx_seen
#define THoughTransform_hxx_seen

#include <ITracking.hxx>

#include <vector>
#include <map>

#include <TObject.h>
#include <TVector3.h>
#include <TRandom.h>

#include <ICOMETLog.hxx>

#include <IG4Trajectory.hxx>
#include <IG4VHit.hxx>
#include <IHandle.hxx>
#include <ICOMETEvent.hxx>

#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TCanvas.h"


class IHoughTransform: public ITracking {
private:
  /*** Trigger Information ***/
  std::vector < std::pair < std::vector< int >, std::vector< int> > > fPairCandidates;
  std::vector < Int_t > fScintUpIndex;
  std::vector < Int_t > fScintDownIndex;
  std::vector < Int_t > fCherenUpIndex;
  std::vector < Int_t > fCherenDownIndex;  // Up & Down Indices Normalized into fCTHSegNum Scale (<64)

  /*** Geometry Variables ***/
  double fDiskRad=10.0;
  int    fCTHSegNum=64;
  int    fNumOfLayers=18;
  double fScintRad=48.28;
  double fScintWidth=9.;
  double fScintHeight=0.5;
  double fScintTiltAngle=13.;
  double fCherenRad=44.78;
  double fCherenWidth=9.;
  double fCherenHeight=1.;
  double fCherenTiltAngle=20.;
  int    fNumOfWiresPerLayer[18]
  ={198,204,210,216,222,228,234,240,246,252,258,264,270,276,282,288,294,300}; // Count only "Actual" Sense Wires
  double fLayerRadius[18]
  ={53.0, 54.6, 56.2, 57.8, 59.4, 61.0, 62.6, 64.2, 65.8, 67.4, 69.0, 70.6, 72.2, 73.8, 75.4, 77.0, 78.6, 80.2};
  double fLayerStereo_TDR[18]
  ={-67.899, 67.640, -67.384, 67.129, -66.876, 66.625, -66.376, 66.129, -65.884, 65.640, -65.398, 65.158, -64.920, 64.683, -75.132, 74.862, -74.593, 74.326}; // After ICEDUST modified we should use this one!
  double fLayerStereo_CAL[18]
  ={-67.004, 66.778, -66.553, 66.327, -66.102, 65.877, -65.652, 65.428, -65.205, 64.982, -64.761, 64.540, -64.319, 64.100, -74.472, 74.220, -73.969, 73.719}; // Up to now, we use this one...
  //double fZOffset[18]
  //={-73.6843, -73.9348,-74.2353,-74.5358, -74.7863, -75.0868, 0,0,0,0,0,0,0,0,0,0,0,0};
  

  /*** HoughTransform Variables ***/
  int fnIter;
  int fnBins;
  double fnPt;
  double fRhoMax;
  double fRhoMin;
  int fBandWidth;
  double fRadUncertainty;
  std::pair <double, double> fRef;

  /*** Hit Variables ***/
  int fnEvenhits=0;
  int fnOddhits=0;    
  double fWireEnd0X_even[10000];
  double fWireEnd0Y_even[10000];  
  double fWireEnd0X_odd[10000];
  double fWireEnd0Y_odd[10000];
  
  double fConfX[10000];
  double fConfY[10000];
  
  int fnConfEven=0;
  int fnConfOdd=0;
  double fConfX_even[3000];
  double fConfX_odd[3000];
  double fConfY_even[3000];
  double fConfY_odd[3000];
  
  std::vector<Int_t> fWireIdsPerLayer[18];
  std::vector< std::vector <Int_t> > fTmpClusterSet;
  std::vector< std::vector <Int_t> > fClusterSet;

  /*** Process Varaibles ***/
  std::pair <Double_t,Double_t> fRef_even;
  std::pair <Double_t,Double_t> fRef_odd;
  double fCx_even;
  double fCy_even;
  double fCx_odd;
  double fCy_odd;
  double fRad_even;
  double fRad_odd;
  double fFitpT;
  double fTruthpT;
  double fOuterR_even;
  double fInnerR_even;
  double fOuterR_odd;
  double fInnerR_odd;
  
  double fAbsCx_even;
  double fAbsCy_even;
  double fAbsCx_odd;
  double fAbsCy_odd;

  int    fnRecoHit_even;
  double fRecoWireEnd0X_even[10000];
  double fRecoWireEnd0Y_even[10000];
  int    fnRecoHit_odd;
  double fRecoWireEnd0X_odd[10000];
  double fRecoWireEnd0Y_odd[10000];


  int    fnRecoHit;
  double fRecoWireEnd0X[10000];
  double fRecoWireEnd0Y[10000];
  double fRecoWireEnd0Z[10000];
  double fRecoWireEnd1X[10000];
  double fRecoWireEnd1Y[10000];
  double fRecoWireEnd1Z[10000];
  double fRecoDriftDist[10000];
  int    fRecoWireLayerId[10000];
  int    fRecoWireId[10000];
  int    fReco2DCharge;
  bool   fRecoCL3;
  int    fRecoMaxWireLayerId;



public:
  IHoughTransform(const char*name, const char* title);
  virtual ~IHoughTransform();
  
  /// called at the begin of run or else (should not be in event-by-event)
  int  Init();

  void LoadMCHits(COMET::IHandle<COMET::IHitSelection> hitHandle, COMET::IHandle<COMET::IG4TrajectoryContainer> trajectories, COMET::IHandle<COMET::IG4HitContainer> cdcHits){
    ITracking::LoadMCHits(hitHandle, trajectories, cdcHits);
  }  
  void ImportTriggerInfo(std::vector < std::pair < std::vector< int >, std::vector< int > > > PairCandidates){
    fPairCandidates = PairCandidates;
  }
  void SetHoughTransformVariables(int nIter, int nBins, double nPt, double rhoMax, double rhoMin, int bandWidth, double rad_uncertainty, double refX, double refY);
  bool PreTrackCut();
  void GetLocalCoordinate();
  void SortHits(); 
  void ConfigureHitClusters();
  double ConfTransX(double x, double y) { return x/(pow(x,2)+pow(y,2));}
  double ConfTransY(double x, double y) { return y/(pow(x,2)+pow(y,2));}
  double HoughTrans(double x, double y, double theta){ return x*TMath::Cos(theta*TMath::Pi()/180)+y*TMath::Sin(theta*TMath::Pi()/180);}
  std::vector<std::pair<double,double> > MakeOrigins(double rad_uncertainty, std::pair<double,double> ref);
  void Process();
  void RecognizeHits();
  bool GetCL3()            {return fRecoCL3;}
  int  GetMaxWireLayerId() {return fRecoMaxWireLayerId;}
  int  Get2DCharge()       {return fReco2DCharge;}
  void DrawEvent(TCanvas* canvas);
  void GetMasterCoordinate(); // Not used
  std::vector<int> GetRecoWireId();
  void PrintMCStatus();
  void PrintResults();

  /// called at the end of run or else (should not be in event-by-event)
  int  Finish();

};
#endif
