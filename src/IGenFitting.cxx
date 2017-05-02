/////////////////////////////////////////////////////////////////
//
// This class used to perform the GenFit for all track candidates.
//
/////////////////////////////////////////////////////////////////
#include <vector>
#include <iostream>

#include <IOADatabase.hxx>
#include <IGenFitting.hxx>
#include <IReconTrack.hxx>
#include <IReconTrackCand.hxx>

#include <ICOMETLog.hxx>
#include <ICOMETEvent.hxx>
#include <IHitSelection.hxx>
#include <IMCHit.hxx>
#include <IHit.hxx>
#include <IGeomInfo.hxx>
#include <IFieldManager.hxx>
#include <IGeoField.hxx>

#include <TCanvas.h>
#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>

////// GenFit
#include <FieldManager.h>
#include <MaterialEffects.h>
#include <TGeoMaterialInterface.h>
#include <AbsMeasurement.h>
#include <WirePointMeasurement.h>
#include <PlanarMeasurement.h>
#include <SpacepointMeasurement.h>
#include <AbsTrackRep.h>
#include <RKTrackRep.h>
#include <StateOnPlane.h>
#include <MeasuredStateOnPlane.h>
#include <Track.h>
#include <KalmanFitter.h>
#include <KalmanFitterRefTrack.h>
#include <KalmanFitterInfo.h>
#include <KalmanFittedStateOnPlane.h>
#include <KalmanFitStatus.h>
#include <DAF.h>
#include <ConstField.h>
#include <EventDisplay.h>

#include <TGeoManager.h>
#include <TGeoGlobalMagField.h>


/// For translating COMET Field into genfit::AbsBField


class GFGeoField : public genfit::AbsBField{
public:
  TVector3 get(const TVector3& pos) const{
    double B[3];
    get(pos.X(),pos.Y(),pos.Z(),B[0],B[1],B[2]);
    return TVector3(B[0],B[1],B[2]);
  }
  void get(const double& posX, const double& posY, const double& posZ, double& Bx, double& By, double& Bz) const{
    double xyz[3]={posX,posY,posZ},B[3];
    TGeoGlobalMagField::Instance()->Field(xyz,B);
    Bx=B[0];
    By=B[1];
    Bz=B[2];
    //    std::cout << posX << "  " << posY << "   " << posZ << "   " << Bx << "   " << By <<  "   " << Bz << std::endl;
  }
};


IGenFitting::IGenFitting(const char* name, const char* title)
  :ITracking(name, title)
  ,fUseBetheBloch(true)
  ,fUseBrems(false)
  ,fMethod("KalmanFitter")
  ,fPID(11)
  ,fMinIterations(10)
  ,fMaxIterations(50)
  ,fMinHitsInTrack(10)
  ,fMinNDF(4)
  ,fMaxMomentum(250)
  ,fMinMomentum(50)
  ,fMaxMomDiff(20)
  ,fUseExtGeomFile(false)
  ,fGeometry("/home/bomki/ICEDUST/BOMKI_analysis/v999/utilities/Geometry/COMETGeometry_Full.root")
  ,fUseExtFieldFile(false)
  ,fFieldMap("/home/bomki/ICEDUST/local_storage/fieldmaps/150630_defaultFieldmap/load_fieldmaps.mac")
  ,fUseMCTruth(true)
  ,fSmearing(true)
  ,fSigmaD(0.02)
{
  ; //if (fSaveTree) fTree = new TTree("gftree", "GenFit tree");
}

IGenFitting::~IGenFitting()
{
  //if (fSaveTree) delete fTree;
}

int IGenFitting::Init(){
  COMETInfo("//----------------------------------------------------------//");
  COMETInfo("// Initialize IGenFitter");

  /// Geometry & materials
  
  /// Set material effects
  COMET::IOADatabase::Get().Geometry();
  genfit::MaterialEffects* mateff = genfit::MaterialEffects::getInstance();
  mateff->setEnergyLossBetheBloch(fUseBetheBloch);
  mateff->setEnergyLossBrems(fUseBrems);  
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());

  // Uniform Field Setup
  TGeoGlobalMagField::Instance()->SetField(new TGeoUniformMagField(10,0,0));
  genfit::FieldManager::getInstance()->init(new GFGeoField());

  return 0;
}

int IGenFitting::BeginOfEvent(){
  /// Field Map
  return 0;
}

int IGenFitting::EndOfEvent()
{
  return 0;
}    

void IGenFitting::LoadHitsAfterHT(COMET::IHandle<COMET::IHitSelection> hitHandle, IHoughTransform* hough){  
  std::vector<int> wireId;
  wireId=hough->GetRecoWireId();
  fnCALCDCHit=wireId.size();
  for (int i=0; i<wireId.size(); i++){
    int wire = wireId.at(i);
    fWireId[i]=wire;
    fWireEnd0X[i] = COMET::IGeomInfo::Get().CDC().GetWireEnd0(wire).X();
    fWireEnd0Y[i] = COMET::IGeomInfo::Get().CDC().GetWireEnd0(wire).Y();
    fWireEnd0Z[i] = COMET::IGeomInfo::Get().CDC().GetWireEnd0(wire).Z();
    fWireEnd1X[i] = COMET::IGeomInfo::Get().CDC().GetWireEnd1(wire).X();
    fWireEnd1Y[i] = COMET::IGeomInfo::Get().CDC().GetWireEnd1(wire).Y();
    fWireEnd1Z[i] = COMET::IGeomInfo::Get().CDC().GetWireEnd1(wire).Z();
  }
}  

int IGenFitting::doFit(){

  //genfit::Track * tmpTrack(NULL);
  //fFitTrack = tmpTrack;
  //genfit::Track *fitTrack(NULL);

  genfit::AbsTrackRep *rep = new genfit::RKTrackRep(fPID);

  TRandom rand;
  Double_t momMag = 105;
  Double_t px, py, pz;
  rand.Sphere(px,py,pz,105);
  /// Set initial position and state
  TVector3 posInit = TVector3(640.,0.,765.);// cm unit
  TVector3 momInit = TVector3(px,py,pz); //    MeV unit
  
  /// Temporal solution to convert the unit
  //posInit *= 0.1;   /// Convert from mm to cm
  momInit *= 0.001; /// Convert from MeV to GeV

  COMETNamedInfo("IGenFitter", "Initial position = (" <<  posInit.x() << ", " << posInit.y() << ", " << posInit.z() << ") cm" <<
                               ", momentum = (" << 1e3*momInit.x() << ", " << 1e3*momInit.y() << ". " << 1e3*momInit.z() << ") MeV" );

  genfit::MeasuredStateOnPlane stateInit(rep);
  TMatrixDSym covMInit(6);
  /// Position resolutions
  double resPos[3] = {0,0,0};
  resPos[0] = resPos[1] = resPos[2] = 0.1;  /// 1 mm
  for (int ii=0; ii<3; ii++) { covMInit(ii, ii) = resPos[ii]*resPos[ii]; }

  bool localCoord=false;

  //set covariant matrix from Momentum and angle resolution
  double sigma_p0     = 0.001; ///   1 MeV
  double sigma_angle0 = 0.1;   /// 100 mrad

  TVector3 dirP    = momInit.Unit();
  TVector3 zaxis   = localCoord ? TVector3(0,0,1) : TVector3(1,0,0);
  TVector3 v_axis  = dirP.Cross(zaxis).Unit();
  TVector3 u_axis  = dirP.Cross(v_axis).Unit();
  TVector3 axis[3] = {dirP, v_axis, u_axis};
  TMatrixD hrot(3, 3);
  for(int ii=0; ii<3; ii++) for(int jj=0; jj<3; jj++) hrot(ii,jj) = axis[jj](ii);
  TMatrixD cov0(3, 3);
  cov0(0,0) = pow(sigma_p0, 2.);
  cov0(1,1) = pow(sigma_angle0*momInit.Mag(), 2.);
  cov0(2,2) = pow(sigma_angle0*momInit.Mag(), 2.);
  TMatrixD covrot(3, 3);
  covrot    = hrot*cov0;
  covrot   *= hrot.T();

  for(int ii=0; ii<3; ii++) for(int jj=0; jj<3; jj++) covMInit[ii+3][jj+3] = covrot(ii, jj);

  rep->setPosMomCov(stateInit, posInit, momInit, covMInit);

  /// Set Track
  TVectorD seedState(6);
  TMatrixDSym seedCov(6);
  rep->get6DStateCov(stateInit, seedState, seedCov);
  fitTrack = new genfit::Track(rep, seedState, seedCov);

  std::vector<genfit::TrackPoint*>     trackPoints; // each track point
  std::vector<genfit::AbsMeasurement*> mesHits;     // each measured hit

  int nTotalHits = 0;
  int nVirtualPlanes =  0;
  
  if(!fUseMCTruth){

    for (int i_hit=0; i_hit<fnCALCDCHit ; i_hit++){
      
      TVectorD wireMes(7);       /// Wire end0(x,y,z), end1(x,y,z), drift distance
      TMatrixDSym wireMatrix(7); /// Uncertainties for wireMes values
      
      //std::cout << fWireEnd0X[i_hit] << "  " << fWireEnd0Y[i_hit] << "  " << fWireEnd0Z[i_hit] << std::endl;
      wireMes[0] = fWireEnd0X[i_hit];
      wireMes[1] = fWireEnd0Y[i_hit];
      wireMes[2] = fWireEnd0Z[i_hit];
      wireMes[3] = fWireEnd1X[i_hit];
      wireMes[4] = fWireEnd1Y[i_hit];
      wireMes[5] = fWireEnd1Z[i_hit];
      wireMes[6] = 0.;
      
      trackPoints.push_back(new genfit::TrackPoint());
      if (fSmearing) wireMes[6] += gRandom->Gaus(0, fSigmaD); /// drift distance is smeared by fSigmaD
     
      for (int row = 0; row < 7; row++) {
	for (int col = 0; col < 7; col++) {
	  if (row != col) wireMatrix(row,col) = 0;
	  if (row < 6)    wireMatrix(row,col) = std::pow(0.001,2);   /// Wire position uncertainties (temporary 10um)
	  else            wireMatrix(row,col) = std::pow(fSigmaD,2); /// Resolution of drift distance
	}
      }
      mesHits.push_back( new genfit::WireMeasurement(wireMes, wireMatrix, fWireId[i_hit], i_hit, trackPoints.at(i_hit)) );  
      
      trackPoints.at(i_hit)->addRawMeasurement(mesHits.at(i_hit));
      fitTrack->insertPoint(trackPoints.at(i_hit), nTotalHits++);  
    }
    
  }

  else if (fUseMCTruth){

    double Sigma[3]={0.02,0.02,0.3};

    for (Int_t i_hit=0; i_hit<fnCDCHit; i_hit++){
      
      TVectorD pos(3);
      TMatrixDSym posMatrix(3);
      pos[0]=fCDCHitX[i_hit];
      pos[1]=fCDCHitY[i_hit];
      pos[2]=fCDCHitZ[i_hit];
      trackPoints.push_back(new genfit::TrackPoint());

      //smearing                                                  
      
      for (Int_t i_xyz=0; i_xyz<3; i_xyz++){
        pos[i_xyz] += gRandom->Gaus(0, Sigma[i_xyz]/TMath::Sqrt(3.));
      }

      for (int row = 0; row < 3; row++) {
        for (int col = 0; col < 3; col++) {
          if (row==col) posMatrix(row,col) = std::pow(Sigma[row]/TMath::Sqrt(3.),2);
          else posMatrix(row,col) = 0;
        }
      }
      
      mesHits.push_back( new genfit::SpacepointMeasurement(pos, posMatrix, i_hit, i_hit, trackPoints.at(i_hit)));
      trackPoints.at(i_hit)->addRawMeasurement(mesHits.at(i_hit));
      fitTrack->insertPoint(trackPoints.at(i_hit), nTotalHits++);                                
    }
  }

  /// Set a fitter
  genfit::AbsKalmanFitter *kalman = NULL;
  if (fMethod=="DAF") {
    kalman = new genfit::DAF();
  } else if (fMethod=="KalmanFitter") {
    kalman = new genfit::KalmanFitter(fMaxIterations);
    kalman->setMinIterations(fMinIterations);
  } else if (fMethod=="KalmanFitterRefTrack") {
    kalman = new genfit::KalmanFitterRefTrack(fMaxIterations);
    kalman->setMinIterations(fMinIterations);
  } else {
    COMETNamedInfo("IGenFitter", "This method is invalid: " << fMethod);
    delete fitTrack;
    return NULL;
  }

  /// Do the fitting
  try{
    kalman->processTrack(fitTrack);
  }catch(genfit::Exception& e){e.what();}

  if (!fitTrack->getFitStatus(rep)->isFitConverged()) {
    kalman->processTrack(fitTrack);
  }
  /*
  if (!fitTrack->getFitStatus(rep)->isFitted()||!fitTrack->getFitStatus(rep)->isFitConverged()) {

    COMETNamedInfo("IGenFitter", "Fitting is failed...");
    std::cout << "Fitting is failed..." << std::endl;
    delete kalman;
    delete fitTrack;
    return NULL;
  }
  */
  /*
  if (fitTrack->getFitStatus(rep)->getChi2()<=0 ||
      (fitTrack->getFitStatus(rep)->getNdf()-2*nVirtualPlanes)<fMinNDF) {
    COMETNamedInfo("IGenFitter", "Fit result might be wrong... (chi2,ndf) = (" << 
                   fitTrack->getFitStatus(rep)->getChi2() << "," << 
                   (fitTrack->getFitStatus(rep)->getNdf()-2*nVirtualPlanes) << ")");
    std::cout << "Fit result might be wrong..." << std::endl;
    delete kalman;
    delete fitTrack;
    return NULL;
  }
  */

  genfit::TrackPoint* tp = fitTrack->getPointWithMeasurementAndFitterInfo(0, rep);
  genfit::KalmanFittedStateOnPlane kfsop(*(static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(rep))->getBackwardUpdate()));
  const TVectorD& state = kfsop.getState();
  const TMatrixDSym& cov = kfsop.getCov();

  double pFit = TMath::Abs(1./state[0]);
  std::cout << "Fitted Momentum is: " << pFit << std::endl;
  std::cout << std::endl;

  //fFitTrack = fitTrack;

  delete kalman;
  //delete fitTrack;
  
  return 1;
}

