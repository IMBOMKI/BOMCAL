/////////////////////////////////////////////////////////////////
//
// This class used to perform the GenFit for all track candidates.
//
/////////////////////////////////////////////////////////////////
#include <vector>
#include <iostream>
#include <assert.h>

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
#include <ProlateSpacepointMeasurement.h>
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
  ,fUseNoiseBetheBloch(true)
  ,fUseNoiseCoulomb(false)
  ,fUseBrems(false)
  ,fUseNoiseBrems(false)
  ,fUseNoEffects(false)
  ,fUseTransverseSeed(false)
  ,fUseMCTruth(false)
  ,fDriftSmearing(true)
  ,fUseDiscreteDriftDistance(false)
  ,fDriftDistanceResolution(0.2)
  ,fSetMaxDriftDistance(1.2)
  ,fSigmaD(0.02)
  ,fSigmaWP(0.001)
   //,fMethod("KalmanFitterRefTrack")
   //,fMethod("KalmanFitter") // "KalmanFitter -> Good for Spacepoint Measurement
  ,fMethod("DAF")       // DAF -> Work for WireMeasurement
  ,fPID(11)
  ,fMinIterations(10)
  ,fMaxIterations(60)
  ,fMinNDF(4)
  ,fSaveHistogram(true)
{
}

IGenFitting::~IGenFitting()
{
}

int IGenFitting::Init(){
  COMETInfo("//----------------------------------------------------------//");
  COMETInfo("// Initialize IGenFitter");

  /// Geometry & materials
  
  /// Set material effects
  COMET::IOADatabase::Get().Geometry();
  genfit::MaterialEffects* mateff = genfit::MaterialEffects::getInstance();
  mateff->setEnergyLossBetheBloch(fUseBetheBloch);
  mateff->setNoiseBetheBloch(fUseNoiseBetheBloch);
  mateff->setNoiseCoulomb(fUseNoiseCoulomb); 
  mateff->setEnergyLossBrems(fUseBrems);  
  mateff->setNoiseBrems(fUseNoiseBrems);  
  mateff->setNoEffects(fUseNoEffects);  

  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());

  // Field Setup

    
  //TGeoGlobalMagField::Instance()->SetField(new TGeoUniformMagField(10,0,0));

  // Check Material 
  genfit::TGeoMaterialInterface* geoMat = new genfit::TGeoMaterialInterface();
  double density;
  double Z;
  double A;
  double radLength;
  double mEE;

  if (geoMat->initTrack(640,0,760,0,1,0)){
    geoMat->getMaterialParameters(density,Z,A,radLength,mEE);
    std::cout << "densitiy: " << density << std::endl;
    std::cout << "Z: " << Z << std::endl;
    std::cout << "A: " << A << std::endl;
    std::cout << "radLength: " << radLength << std::endl;
    std::cout << "mEE: " << mEE << std::endl;
  }



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

void IGenFitting::LoadFieldMap(bool hasCalledFieldMap){
  if (hasCalledFieldMap==0) COMET::IFieldManager::Import();  
  if (COMET::IFieldManager::getObject()) {
    TGeoGlobalMagField::Instance()->SetField(new COMET::IGeoField());
  }
  //fFieldManager = genfit::FieldManager::getInstance();
  genfit::FieldManager::getInstance()->init(new GFGeoField());
}

void IGenFitting::LoadHitsAfterHT(COMET::IHandle<COMET::IHitSelection> hitHandle, IHoughTransform* hough){  
  std::vector<int> wireId      = hough->GetRecoWireId();
  std::vector<double> driftDist= hough->GetRecoDriftDist();
  std::vector<int> domain       = hough->GetRecoDomain();

  fnCALCDCHit=wireId.size();
  assert(fnCALCDCHit==hough->GetNumberOfRecognizedHits());
  assert(wireId.size() == driftDist.size());
  assert(wireId.size() == domain.size());

  for (int i=0; i<fnCALCDCHit; i++){
    int wire = wireId.at(i);
    fWireId[i]=wire;
    fWireEnd0X[i] = COMET::IGeomInfo::Get().CDC().GetWireEnd0(wire).X();
    fWireEnd0Y[i] = COMET::IGeomInfo::Get().CDC().GetWireEnd0(wire).Y();
    fWireEnd0Z[i] = COMET::IGeomInfo::Get().CDC().GetWireEnd0(wire).Z();
    fWireEnd1X[i] = COMET::IGeomInfo::Get().CDC().GetWireEnd1(wire).X();
    fWireEnd1Y[i] = COMET::IGeomInfo::Get().CDC().GetWireEnd1(wire).Y();
    fWireEnd1Z[i] = COMET::IGeomInfo::Get().CDC().GetWireEnd1(wire).Z();
    fDriftDist[i] = driftDist.at(i);
    fDomain[i]      = domain.at(i);
    fWireLayerId[i]=COMET::IGeomInfo::Get().CDC().GetLayer(wire);
    //std::cout << i << "   " << fWireEnd0X[i] << "   " << fWireEnd0Y[i] << "   " << fWireEnd0Z[i] << std::endl;
  }
}  

int IGenFitting::DoFit(IHoughTransform* hough){

  gRandom->SetSeed(1);
  
  genfit::AbsTrackRep *rep = new genfit::RKTrackRep(fPID);

  /// Set initial position and state

  // 1. Random Momentum
  //TVector3 posInit = TVector3(640.,0.,765.);// cm unit
  //TVector3 momInit = TVector3(px,py,pz); //    MeV unit

  // 2. True Birth value  
  //TVector3 posInit = TVector3(fGenTrX,fGenTrY,fGenTrZ);     // cm unit
  //TVector3 momInit = TVector3(fGenTrPx,fGenTrPy,fGenTrPz); // MeV unit

  // 3. True Entering value  
  posInit = TVector3(fCDCEnterX,fCDCEnterY,fCDCEnterZ);     // cm unit
  //posInit = TVector3(fCDCHitX[0],fCDCHitY[0],fCDCHitZ[0]);     // cm unit
  momInit = TVector3(fCDCEnterPx,fCDCEnterPy,fCDCEnterPz); // MeV unit

  // (Optional) Smearing TVector3
  //TVector3 posSmear(gRandom->Gaus(0.3),gRandom->Gaus(0.02),gRandom->Gaus(0.02));
  //TVector3 momSmear(gRandom->Gaus(0.1),gRandom->Gaus(0.1),gRandom->Gaus(0.1));  
  //posInit += posSmear;
  //momInit += momSmear;
  

  if (fUseTransverseSeed==1){
    TVector3 fFitEnterPos = hough->GetEnterXYPair_Reseeded().second;
    TVector3 fFitEnterMom = hough->GetEnterPxPyPair_Reseeded().second;
    if(!COMET::IGeomInfo::DetectorSolenoid().GetDetPositionInGlobalCoordinate(fFitEnterPos, fFitEnterPos)){std::cout << "Coordinate change fails (Local to Master)" << std::endl; return true;}
    Double_t fFitEnterY  = fFitEnterPos(1);
    Double_t fFitEnterZ  = fFitEnterPos(2);
    Double_t fFitEnterPz = -fFitEnterMom(0);
    Double_t fFitEnterPy = fFitEnterMom(1);

    posInit = TVector3(640,fFitEnterY,fFitEnterZ);
    momInit = TVector3(0  ,fFitEnterPy,fFitEnterPz);
    //posInit = TVector3(fCDCEnterX,fFitEnterY,fFitEnterZ);
    //momInit = TVector3(fCDCEnterPx  ,fFitEnterPy,fFitEnterPz);
  }
  /// Temporal solution to convert the unit
  momInit *= 0.001; /// Convert from MeV to GeV

  COMETNamedInfo("IGenFitter", "Initial position = (" <<  posInit.x() << ", " << posInit.y() << ", " << posInit.z() << ") cm" <<
                               ", momentum = (" << 1e3*momInit.x() << ", " << 1e3*momInit.y() << ". " << 1e3*momInit.z() << ") MeV" );

  genfit::MeasuredStateOnPlane stateInit(rep);

  TMatrixDSym covMInit(6);
  /// Position resolutions
  double resPos[3] = {2.,2.,2.};
  for (int ii=0; ii<3; ii++) { covMInit(ii, ii) = resPos[ii]*resPos[ii]; }
  for (int ii=0; ii<3; ii++) { covMInit(ii+3, ii+3) = pow(momInit[ii]*0.2,2);} // 0.1 %


  //set covariant matrix from Momentum and angle resolution
  //double sigma_p0     = 0.001; ///   1 MeV
  /*
  double sigma_angle0 = 0.1;   /// 100 mrad
  TVector3 dirP    = momInit.Unit();
  TVector3 zaxis   = TVector3(1,0,0);
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
  */

  rep->setPosMomCov(stateInit, posInit, momInit, covMInit);

  /// Set Track
  TVectorD seedState(6);
  TMatrixDSym seedCov(6);
  rep->get6DStateCov(stateInit, seedState, seedCov);

  // remember original initial state
  const genfit::StateOnPlane stateRefOrig(stateInit);

  fitTrack = new genfit::Track(rep, seedState, seedCov);
  std::vector<genfit::TrackPoint*>     trackPoints; // each track point
  std::vector<genfit::AbsMeasurement*> mesHits;     // each measured hit

  int nTotalHits = 0;
  int nVirtualPlanes =  0;
  
  if(!fUseMCTruth){

    for (int i_hit=0; i_hit<fnCALCDCHit ; i_hit++){

      //if (fWireLayerId[i_hit]<=1) continue;
      
      TVectorD wireMes(7);       /// Wire end0(x,y,z), end1(x,y,z), drift distance
      TMatrixDSym wireMatrix(7); /// Uncertainties for wireMes values
      
      //std::cout << "Hit Info " << i_hit << "  " << fWireEnd0X[i_hit] << "  " << fWireEnd0Y[i_hit] << "  " << fWireEnd0Z[i_hit] << "  " << fDriftDist[i_hit] <<  std::endl;

      wireMes[0] = fWireEnd0X[i_hit];
      wireMes[1] = fWireEnd0Y[i_hit];
      wireMes[2] = fWireEnd0Z[i_hit];
      wireMes[3] = fWireEnd1X[i_hit];
      wireMes[4] = fWireEnd1Y[i_hit];
      wireMes[5] = fWireEnd1Z[i_hit];
      wireMes[6] = fDriftDist[i_hit];
      if (wireMes[6]>fSetMaxDriftDistance) wireMes[6]=fSetMaxDriftDistance-0.001;
      if (fUseDiscreteDriftDistance){
	Int_t div = wireMes[6]/fDriftDistanceResolution;
	wireMes[6] = fDriftDistanceResolution * div + fDriftDistanceResolution/2.;	
      }
      
      trackPoints.push_back(new genfit::TrackPoint());


      if (fDriftSmearing) wireMes[6] += gRandom->Gaus(0, fSigmaD); /// drift distance is smeared by fSigmaD
     
      for (int row = 0; row < 7; row++) {
	for (int col = 0; col < 7; col++) {
	  if (row != col) wireMatrix(row,col) = 0.;
	  if (row < 6)    wireMatrix(row,col) = std::pow(fSigmaWP,2);   /// Wire position uncertainties (temporary 10um)
	  //if (row < 6)    wireMatrix(row,col) = 0.;
	  else            wireMatrix(row,col) = std::pow(fSigmaD,2); /// Resolution of drift distance
	}
      }

      mesHits.push_back( new genfit::WireMeasurement(wireMes, wireMatrix, fWireId[i_hit], nTotalHits, trackPoints.at(nTotalHits)) );        
      trackPoints.at(nTotalHits)->addRawMeasurement(mesHits.at(nTotalHits));
      fitTrack->insertPoint(trackPoints.at(nTotalHits-1), nTotalHits++);  

      //mesHits.push_back( new genfit::WireMeasurement(wireMes, wireMatrix, fWireId[i_hit], i_hit, trackPoints.at(i_hit)) );        
      //trackPoints.at(i_hit)->addRawMeasurement(mesHits.at(i_hit));
      //fitTrack->insertPoint(trackPoints.at(i_hit), nTotalHits++);  

      //genfit::WireMeasurement *measure = new genfit::WireMeasurement(wireMes, wireMatrix, fWireId[i_hit], i_hit, NULL);
      //measure->setMaxDistance(1.2);
      //fitTrack->insertMeasurement(measure);
    }    
  }

  else if (fUseMCTruth){

    double Sigma[3]={0.01,0.01,0.01};

    for (Int_t i_hit=0; i_hit<fnCDCHit; i_hit++){
      if (fPrimary[i_hit]!=1) continue;

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
      
      ////mesHits.push_back( new genfit::ProlateSpacepointMeasurement(pos, posMatrix, i_hit, nTotalHits, trackPoints.at(nTotalHits)));

      genfit::ProlateSpacepointMeasurement *measure = new genfit::ProlateSpacepointMeasurement(pos, posMatrix, i_hit,nTotalHits, NULL);
      measure->setLargestErrorDirection(TVector3(1,0,0));
      fitTrack->insertMeasurement(measure);

      //mesHits.push_back( new genfit::SpacepointMeasurement(pos, posMatrix, i_hit, nTotalHits, trackPoints.at(nTotalHits)));
      //trackPoints.at(nTotalHits)->addRawMeasurement(mesHits.at(nTotalHits));
      //fitTrack->insertPoint(trackPoints.at(nTotalHits-1), nTotalHits++);                                
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
    return 0;
  }
  kalman->setDeltaPval(3e-3);

  /// Do the fitting
  try{
    fitTrack->checkConsistency();
    gGeoManager->ResetState();
    kalman->processTrack(fitTrack);
    fitTrack->checkConsistency();
  }catch(genfit::Exception& e){e.what();}

  if (!fitTrack->getFitStatus(rep)->isFitConverged()) {
    kalman->processTrack(fitTrack);
  }
  
  if (!fitTrack->getFitStatus(rep)->isFitted()||!fitTrack->getFitStatus(rep)->isFitConverged()) {

    COMETNamedInfo("IGenFitter", "Fitting is failed...");
    std::cout << "Fitting is failed..." << std::endl;
    delete kalman;
    return 0;
  }
    
  
  if (fitTrack->getFitStatus(rep)->getChi2()<=0 ||
      (fitTrack->getFitStatus(rep)->getNdf()-2*nVirtualPlanes)<fMinNDF) {
    COMETNamedInfo("IGenFitter", "Fit result might be wrong... (chi2,ndf) = (" << 
                   fitTrack->getFitStatus(rep)->getChi2() << "," << 
                   (fitTrack->getFitStatus(rep)->getNdf()-2*nVirtualPlanes) << ")");
    std::cout << "Fit result might be wrong..." << std::endl;
    delete kalman;
    //delete fitTrack;
    return 0;
  }
  

  genfit::TrackPoint* tp = fitTrack->getPointWithMeasurementAndFitterInfo(0, rep);
  genfit::KalmanFittedStateOnPlane kfsop(*(static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(rep))->getBackwardUpdate()));
  //genfit::KalmanFittedStateOnPlane kfsop(*(static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(rep))->getUpdate(1)));  

  try{
    rep->extrapolateToPlane(kfsop, stateRefOrig.getPlane());
  }
  catch(genfit::Exception& e){ 
    std::cerr<<"Exception, next track"<<std::endl;
    std::cerr << e.what();
    delete kalman;
    return 0;
  }

  const TVectorD& referenceState = stateRefOrig.getState();
  const TVectorD& state = kfsop.getState();
  const TMatrixDSym& cov = kfsop.getCov();
  
  fpFit = TMath::Abs(1./state[0]);
  pVal = kalman->getPVal(fitTrack, rep);
  hqopPu = (state[0]-referenceState[0]) / sqrt(cov[0][0]);
  hupPu  = (state[1]-referenceState[1]) / sqrt(cov[1][1]);
  hvpPu  = (state[2]-referenceState[2]) / sqrt(cov[2][2]);
  huPu   = (state[3]-referenceState[3]) / sqrt(cov[3][3]);
  hvPu   = (state[4]-referenceState[4]) / sqrt(cov[4][4]);

  fChi2 = fitTrack->getFitStatus(rep)->getChi2();
  fNdf = fitTrack->getFitStatus(rep)->getNdf();  
  fChi2Ndf = fChi2/fNdf;

  delete kalman;

  fxFit_vec = rep->getPos(kfsop);  
  fpFit_vec = rep->getMom(kfsop)*1000;
  fpFit = fpFit_vec.Mag();

  std::cout << "Fitted Momentum is: " << fpFit << std::endl;
  std::cout << std::endl;

  return 1;
}

Double_t* IGenFitting::GetPullValue() { 
  //Double_t* pointer;
  static Double_t PullValue[5];
  //pointer = PullValue;
  PullValue[0]=hqopPu;
  PullValue[1]=hupPu;
  PullValue[2]=hvpPu;
  PullValue[3]=huPu;
  PullValue[4]=hvPu; 
  return PullValue; 
}
