#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "THStack.h"
#include "TMath.h"
#include "TChain.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <tuple>

string IntToString (int a);

Double_t GetAcceptance(Double_t NumberOfTrig,Int_t NumberOfEvent);

Bool_t IfInTimeWindow (Double_t time, Double_t DetTimeWindow);
 

template <class T1, class T2, class Pred = std::less<T2> >
struct sort_pair_second {
  bool operator()(const std::pair<T1,T2>&left, const std::pair<T1,T2>&right) {
    Pred p;
    return p(left.second, right.second);
  }
};

template <class T1, class T2, class Pred = std::less<T1> >
struct sort_pair_first {
  bool operator()(const std::pair<T1,T2>&left, const std::pair<T1,T2>&right) {
    Pred p;
    return p(left.first, right.first);
  }
};

template<int M, template<typename> class F = std::less>
struct TupleCompare
{
  template<typename T>
  bool operator()(T const &t1, T const &t2)
  {
    return F<typename tuple_element<M, T>::type>()(std::get<M>(t1), std::get<M>(t2));
  }
};

std::vector < std::tuple <Int_t, Bool_t, Double_t > > MergeSTLandCRK(std::vector < std::pair<Int_t, Double_t> > STLIndexAndTime, std::vector < std::pair<Int_t, Double_t> > CRKIndexAndTime);

std::vector < std::vector< std::tuple<Int_t, Bool_t, Double_t> > > MakeTimeCluster(std::vector < std::tuple < Int_t,Bool_t, Double_t > > IndexAndTime);

std::vector < std::tuple<std::vector <Int_t>, std::vector <Int_t>, Double_t > > ApplyIntervalCond(std::vector < std::vector< std::tuple<Int_t, Bool_t, Double_t> > > TimeCluster, Int_t Interval=1, Int_t tolerance=1);

std::pair < std::vector <Int_t>, std::vector<Int_t> > ApplyShiftCond(std::vector < std::vector<Int_t> > STLIntervalAppliedCluster, std::vector < std::vector<Int_t> > CRKIntervalAppliedCluster, Int_t tolerance);

void TrigCoincidence(std::string filedir, std::string filename, std::string outputdir){

  /*
  TCanvas *c1 = new TCanvas("c1", "CTH Algorithm", 600, 600);
  c1->SetTitle("CTH Algorithm");  
  c1->cd();
  */

  /*---- Plotting Information ----*/

  Int_t FourFold=0;
  Int_t FourFold_WithECut=0;
  Int_t plusTrackCut=0;
  Double_t ECut=88.32;

  Int_t plusnCDCHit=0;
  Int_t plusCL3=0;
  Int_t plusMaxLayer=0;


  /*---- Trigerring Information ----*/

  Int_t NumOfCTH=64;
  Int_t s_tol=1;

  Double_t DetTimeWindow=700; // <- After 700 // Not used here
  
  /*---- Input File Information ----*/

  std::string InputDir = filedir;
  std::string InputFileName = filename;

  TFile *f = TFile::Open(TString(filedir+filename));
  TTree *t = (TTree*)f->Get("trdata");

  /*---- Output File Information ----*/

  std::string OutputDir = outputdir;
  std::string OutputFileName = "trig_"+filename;
  
  TFile *f_new = new TFile(TString(OutputDir+OutputFileName), "recreate");
  TTree *t_new = t->CloneTree(0);

  /*---- Input File Set Branch Address ----*/

  Int_t nTrigSTL;
  Int_t nTrigCRK;
  Int_t TrigSTL[1000];
  Int_t TrigCRK[1000];
  Double_t TrigSTLTime[1000];
  Double_t TrigCRKTime[1000];


  Int_t nCDCHit_Pairep;
  Int_t nCDCHit_Pairem;
  Int_t nCDCHit_Comptonem;
  Bool_t ifCDC_Pairep[10000];
  Bool_t ifCDC_Pairem[10000];
  Bool_t ifCDC_Comptonem[10000];
  Int_t CDCLayerId_Pairep[10000];
  Int_t CDCLayerId_Pairem[10000];
  Int_t MaxCDCLayerId_Pairep;
  Int_t MaxCDCLayerId_Pairem;
  Bool_t ifSingleTurn_Pairep;
  Bool_t ifMultiTurn_Pairep;

  Int_t nCDCHit;
  Double_t CDCHitX[100000];
  Double_t CDCHitY[100000];
  Double_t CDCHitZ[100000];
  Double_t CDCHitT[100000];
  Double_t CDCEDep[100000];

  Bool_t ifPairProdOccurs;
  Bool_t ifPairVertexAtTarget;
  Double_t PairVertexX;
  Double_t PairVertexY;
  Double_t PairVertexZ;
  Double_t PairVertexT;
  Double_t Pairep_genTrPx;
  Double_t Pairep_genTrPy;
  Double_t Pairep_genTrPz;
  Double_t Pairep_genTrE;
  Double_t Pairem_genTrPx;
  Double_t Pairem_genTrPy;
  Double_t Pairem_genTrPz;
  Double_t Pairem_genTrE;

  Bool_t ifComptonId2Occurs;
  Bool_t ifComptonVertexAtTarget;
  Double_t ComptonVertexX;
  Double_t ComptonVertexY;
  Double_t ComptonVertexZ;
  Double_t ComptonVertexT;
  Double_t Comptonem_genTrPx;
  Double_t Comptonem_genTrPy;
  Double_t Comptonem_genTrPz;
  Double_t Comptonem_genTrE;

  t->SetBranchAddress("nTrigSTL", &nTrigSTL);
  t->SetBranchAddress("nTrigCRK", &nTrigCRK);    
  t->SetBranchAddress("TrigSTL", TrigSTL);
  t->SetBranchAddress("TrigCRK", TrigCRK);
  t->SetBranchAddress("TrigSTLTime", TrigSTLTime);
  t->SetBranchAddress("TrigCRKTime", TrigCRKTime);

  t->SetBranchAddress("nCDCHit_Pairep", &nCDCHit_Pairep);
  t->SetBranchAddress("nCDCHit_Pairem", &nCDCHit_Pairem);
  t->SetBranchAddress("nCDCHit_Comptonem", &nCDCHit_Comptonem);  
  t->SetBranchAddress("ifCDC_Pairep", ifCDC_Pairep);
  t->SetBranchAddress("ifCDC_Pairem", ifCDC_Pairem);
  t->SetBranchAddress("ifCDC_Comptonem", ifCDC_Comptonem);
  t->SetBranchAddress("CDCLayerId_Pairep", CDCLayerId_Pairep);    
  t->SetBranchAddress("CDCLayerId_Pairem", CDCLayerId_Pairem);
  t->SetBranchAddress("MaxCDCLayerId_Pairep", &MaxCDCLayerId_Pairep);  
  t->SetBranchAddress("MaxCDCLayerId_Pairem", &MaxCDCLayerId_Pairem);  
  t->SetBranchAddress("ifSingleTurn_Pairep", &ifSingleTurn_Pairep);
  t->SetBranchAddress("ifMultiTurn_Pairep", &ifMultiTurn_Pairep);


  t->SetBranchAddress("nCDCHit", &nCDCHit);
  t->SetBranchAddress("CDCHitX", CDCHitX);
  t->SetBranchAddress("CDCHitY", CDCHitY);
  t->SetBranchAddress("CDCHitZ", CDCHitZ);
  t->SetBranchAddress("CDCHitT", CDCHitT);
  t->SetBranchAddress("CDCEDep", CDCEDep);
  
            /* Pair e+ and e- */

  t->SetBranchAddress("ifPairProdOccurs", &ifPairProdOccurs); 
  t->SetBranchAddress("ifPairVertexAtTarget", &ifPairVertexAtTarget);
  t->SetBranchAddress("PairVertexX", &PairVertexX);
  t->SetBranchAddress("PairVertexY", &PairVertexY);
  t->SetBranchAddress("PairVertexZ", &PairVertexZ);
  t->SetBranchAddress("PairVertexT", &PairVertexT);
  t->SetBranchAddress("Pairep_genTrPx", &Pairep_genTrPx);
  t->SetBranchAddress("Pairep_genTrPy", &Pairep_genTrPy);
  t->SetBranchAddress("Pairep_genTrPz", &Pairep_genTrPz);
  t->SetBranchAddress("Pairep_genTrE", &Pairep_genTrE);
  t->SetBranchAddress("Pairem_genTrPx", &Pairem_genTrPx);
  t->SetBranchAddress("Pairem_genTrPy", &Pairem_genTrPy);
  t->SetBranchAddress("Pairem_genTrPz", &Pairem_genTrPz);
  t->SetBranchAddress("Pairem_genTrE", &Pairem_genTrE);

          /* Compton Scattering e- */

  t->SetBranchAddress("ifComptonId2Occurs", &ifComptonId2Occurs); 
  t->SetBranchAddress("ifComptonVertexAtTarget", &ifComptonVertexAtTarget);
  t->SetBranchAddress("ComptonVertexX", &ComptonVertexX);
  t->SetBranchAddress("ComptonVertexY", &ComptonVertexY);
  t->SetBranchAddress("ComptonVertexZ", &ComptonVertexZ);
  t->SetBranchAddress("ComptonVertexT", &ComptonVertexT);
  t->SetBranchAddress("Comptonem_genTrPx", &Comptonem_genTrPx);
  t->SetBranchAddress("Comptonem_genTrPy", &Comptonem_genTrPy);
  t->SetBranchAddress("Comptonem_genTrPz", &Comptonem_genTrPz);
  t->SetBranchAddress("Comptonem_genTrE", &Comptonem_genTrE);


  /*-------------------------------------------------------------

    Trigger Algorithm

   -------------------------------------------------------------*/


  for (Int_t i=0; i < t->GetEntries(); i++){

    t->GetEntry(i);

    // Stage 0 - Make std::pair with Up Down STL & CRK

    std::vector <std::pair<Int_t, Double_t> > Up_STLIndexAndTime;
    std::vector <std::pair<Int_t, Double_t> > Down_STLIndexAndTime;
   
    for (Int_t i_STL=0; i_STL<nTrigSTL; i_STL++){
      //      if ( IfInTimeWindow(TrigSTLTime[i_STL],DetTimeWindow)==1){
	if (TrigSTL[i_STL]<=(NumOfCTH-1)){
	  Up_STLIndexAndTime.push_back(std::make_pair(TrigSTL[i_STL],TrigSTLTime[i_STL]));}    
	else if (TrigSTL[i_STL]>(NumOfCTH-1)){
	  Down_STLIndexAndTime.push_back(std::make_pair(TrigSTL[i_STL],TrigSTLTime[i_STL]));}    
	//  }
    }
    
    std::vector <std::pair<Int_t, Double_t> > Up_CRKIndexAndTime;
    std::vector <std::pair<Int_t, Double_t> > Down_CRKIndexAndTime;

    for (Int_t i_CRK=0; i_CRK<nTrigCRK; i_CRK++){
      //  if ( IfInTimeWindow(TrigCRKTime[i_CRK],DetTimeWindow)==1){
      if (TrigCRK[i_CRK]<=(NumOfCTH-1)){
	Up_CRKIndexAndTime.push_back(std::make_pair(TrigCRK[i_CRK],TrigCRKTime[i_CRK])); }
      else if (TrigCRK[i_CRK]>(NumOfCTH-1)){
	Down_CRKIndexAndTime.push_back(std::make_pair(TrigCRK[i_CRK],TrigCRKTime[i_CRK])); }
      // }            
    }
    
    // Stage 0-a - Merge STL & CRK information and sort with Timing Order (0 -> STL , 1 -> CRK)

    std::vector < std::tuple <Int_t, Bool_t, Double_t > > Up_CTHIndexAndTime = MergeSTLandCRK(Up_STLIndexAndTime, Up_CRKIndexAndTime);
    std::vector < std::tuple <Int_t, Bool_t, Double_t > > Down_CTHIndexAndTime = MergeSTLandCRK(Down_STLIndexAndTime, Down_CRKIndexAndTime);
      
    std::sort(Up_CTHIndexAndTime.begin(), Up_CTHIndexAndTime.end(), TupleCompare<2>());
    std::sort(Down_CTHIndexAndTime.begin(), Down_CTHIndexAndTime.end(), TupleCompare<2>());
   
    // Stage 1 - Cluster elements within 10 ns and sort with Index Number Order and consider STL CRK index separately
        
    std::vector < std::vector< std::tuple<Int_t, Bool_t, Double_t> > > Up_CTHTimeCluster = MakeTimeCluster(Up_CTHIndexAndTime);
    std::vector < std::vector< std::tuple<Int_t, Bool_t, Double_t> > > Down_CTHTimeCluster= MakeTimeCluster(Down_CTHIndexAndTime);

    // Stage 2 - Cluster elements' interval // Stage 3 - STL & CRK Cluster shift
    

    std::vector < std::tuple<std::vector <Int_t>, std::vector <Int_t>, Double_t > > Up_CTHIntervalApplied = ApplyIntervalCond(Up_CTHTimeCluster,1,s_tol);
    std::vector < std::tuple<std::vector <Int_t>, std::vector <Int_t>, Double_t > > Down_CTHIntervalApplied = ApplyIntervalCond(Down_CTHTimeCluster,1,s_tol);
    
    if (!Up_CTHIntervalApplied.empty() || !Down_CTHIntervalApplied.empty()){
      FourFold++;      
      if (Pairep_genTrE>ECut){

	FourFold_WithECut++;	    
	t_new->Fill();

	if (nCDCHit_Pairep>30 && (ifSingleTurn_Pairep==1 || ifMultiTurn_Pairep==1) && MaxCDCLayerId_Pairep>=5){
	  plusTrackCut++;
	}
	/*
	if (nCDCHit_Pairep>=30){
	  plusnCDCHit++;
	  if (ifSingleTurn_Pairep==1 || ifMultiTurn_Pairep==1){
	    plusCL3++;
	    if (MaxCDCLayerId_Pairep>=5){
	      plusMaxLayer++;
	    }
	  }
	}
	*/
      }   
    }
    
    if (i%100000==0)
      std::cout << i <<"-th Event is Analyzed"<< std::endl;
    
  }

  std::cout << "FourFold                   : " << FourFold << std::endl;
  std::cout << "FourFold With E>"<< ECut <<" Cut  : " << FourFold_WithECut << std::endl;
  std::cout << "plusTrackCut               : " << plusTrackCut << std::endl;
  /*
  std::cout << "plusnCDCHit                : " << plusnCDCHit << std::endl;
  std::cout << "plusCL3                    : " << plusCL3 << std::endl;
  std::cout << "plusMaxLayer               : " << plusMaxLayer << std::endl;
  */

  f_new->cd();
  t_new->Write();
  f_new->Close();
  //t_new->Print();

  std::cout << "Finish!" << std::endl;
  f->Close();

  
}



















string IntToString (int a)
{
  ostringstream temp;
  temp<<a;
  return temp.str();
}

Double_t GetAcceptance(Double_t NumberOfTrig,Int_t NumberOfEvent){
  return NumberOfTrig/NumberOfEvent*100;
}

Bool_t IfInTimeWindow (Double_t time, Double_t DetTimeWindow){
  if ((time<1170*1 && time>1170*0+DetTimeWindow)){ // || (time<1170*6 && time>1170*5+DetTimeWindow) ){ // Consider 5th (Saturated Bunch)
    return true;
  }
  else {
    return false;
  }
}

std::vector < std::tuple <Int_t, Bool_t, Double_t > > MergeSTLandCRK(std::vector < std::pair<Int_t, Double_t> > STLIndexAndTime, std::vector < std::pair<Int_t, Double_t> > CRKIndexAndTime){

  std::vector < std::tuple < Int_t, Bool_t, Double_t> > Merge;
  for (int i=0; i<STLIndexAndTime.size(); i++){
    std::tuple < Int_t, Bool_t, Double_t> element = std::make_tuple (STLIndexAndTime[i].first, 0, STLIndexAndTime[i].second);
    Merge.push_back(element);
  }
  for (int i=0; i<CRKIndexAndTime.size(); i++){
    std::tuple < Int_t, Bool_t, Double_t> element = std::make_tuple (CRKIndexAndTime[i].first, 1, CRKIndexAndTime[i].second);
    Merge.push_back(element);
  }
  return Merge;
}



std::vector < std::vector< std::tuple<Int_t, Bool_t, Double_t> > > MakeTimeCluster(std::vector < std::tuple < Int_t,Bool_t, Double_t > > IndexAndTime){
  std::vector < std::vector< std::tuple<Int_t, Bool_t, Double_t> > > TimeCluster;
  for (Int_t i=0; i<IndexAndTime.size(); i++){
    if ((std::get<2>(IndexAndTime[i+1])-std::get<2>(IndexAndTime[i])<10) && (i+1<IndexAndTime.size()) ){
      Double_t refTime=std::get<2>(IndexAndTime[i]);

      std::vector< std::tuple<Int_t, Bool_t, Double_t> > ClusterIndex;
      Int_t j=i;
      
      while (std::get<2>(IndexAndTime[j])-refTime<10 && j<IndexAndTime.size()){

	///////////////// Find If there is an element with same Index and STL/CRK info /////////////////////////////////
	///////////////// If there is not (PushOrNot==1), Push it                      /////////////////////////////////
	Bool_t PushOrNot=1;
	if (ClusterIndex.size()==0) ClusterIndex.push_back(IndexAndTime[j]);

	for (Int_t i_Clu=0; i_Clu<ClusterIndex.size(); i_Clu++){
	  if (std::get<0>(IndexAndTime[j])==std::get<0>(ClusterIndex[i_Clu]) && std::get<1>(IndexAndTime[j])==std::get<1>(ClusterIndex[i_Clu]) ){
	    PushOrNot=0;
	  }
	}
	
	if (PushOrNot==1){
	   ClusterIndex.push_back(IndexAndTime[j]);
	}
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	j++;
      }
      i=j;
     
      std::sort(ClusterIndex.begin(), ClusterIndex.end(), TupleCompare<1>());     

      /*
      for (Int_t i=0; i<ClusterIndex.size(); i++){
	std::cout << std::get<0>(ClusterIndex[i]) << " " << std::get<0>(ClusterIndex[i]) << "   ";
      }
      std::cout << std::endl;
      */

      int Boundary=0;
      for (Int_t i=0; i<ClusterIndex.size(); i++){
	if (std::get<1>(ClusterIndex[i])==0){
	  Boundary++;
	}
      }
      std::sort(ClusterIndex.begin(), ClusterIndex.begin()+Boundary, TupleCompare<0>());      
      std::sort(ClusterIndex.begin()+Boundary, ClusterIndex.end(), TupleCompare<0>());      

      /*
      for (Int_t i=0; i<ClusterIndex.size(); i++){
	std::cout << std::get<0>(ClusterIndex[i]) << " " << std::get<1>(ClusterIndex[i]) << "   ";
      }
      std::cout << std::endl;
      */

      TimeCluster.push_back(ClusterIndex);						   
    }      
  }
  //std::cout << TimeCluster.size() << std::endl;
  return TimeCluster;
}


std::vector < std::tuple<std::vector <Int_t>, std::vector <Int_t>, Double_t > > ApplyIntervalCond(std::vector < std::vector< std::tuple<Int_t, Bool_t, Double_t> > > TimeCluster, Int_t Interval=1, Int_t tolerance=1){
  std::vector < std::tuple<std::vector <Int_t>, std::vector <Int_t>, Double_t > > ShiftAppliedCluster;


  for (Int_t i_vec=0 ; i_vec<TimeCluster.size() ; i_vec++){

    std::vector < std::vector<Int_t> > STLIntervalAppliedCluster;
    std::vector < std::vector<Int_t> > CRKIntervalAppliedCluster;
    
    std::vector< std::tuple <Int_t,Bool_t,Double_t> > ClusterIndex = TimeCluster[i_vec];


    std::vector<Double_t> TimeOfClusterElement;
    for (Int_t i=0; i<ClusterIndex.size(); i++){
      TimeOfClusterElement.push_back(std::get<2>(ClusterIndex[i]));
    }
    Double_t TrigTime=*std::min_element(std::begin(TimeOfClusterElement),std::end(TimeOfClusterElement));

    /*
    for (int p=0; p<ClusterIndex.size(); p++){
      std::cout << ClusterIndex[p].first << " " << ClusterIndex[p].second << "   "  ;
    }
    std::cout << " || ";
    */

    int Boundary=0;
    for (Int_t i=0; i<ClusterIndex.size(); i++){
      if (std::get<1>(ClusterIndex[i])==0){
	  Boundary++;
	}
    }

    ///////////////////////////////////// For STL /////////////////////////////////////////
    for (Int_t i=0; i<Boundary; i++){
      if (std::get<0>(ClusterIndex[i+1])-std::get<0>(ClusterIndex[i])<=Interval && i+1<Boundary){
	std::vector< Int_t > AppliedClusterIndex;
	AppliedClusterIndex.push_back(std::get<0>(ClusterIndex[i]));

	/*
	std::cout << std::get<0>(ClusterIndex[i]) << " " << std::get<1>(ClusterIndex[i]) << "  ";
	*/

	while (std::get<0>(ClusterIndex[i+1])-std::get<0>(ClusterIndex[i])<=Interval && i+1<Boundary){
	  AppliedClusterIndex.push_back(std::get<0>(ClusterIndex[i+1]));
	  i++;
	  /*
	  std::cout << std::get<0>(ClusterIndex[i]) << " " << std::get<1>(ClusterIndex[i]) << "  ";
	  */
	}
	/*
	std::cout << " | "; 
	*/
        STLIntervalAppliedCluster.push_back(AppliedClusterIndex);
      }
    }

    ///////////////////////////////////// For CRK /////////////////////////////////////////
    for (Int_t i=Boundary; i<ClusterIndex.size(); i++){
      if (std::get<0>(ClusterIndex[i+1])-std::get<0>(ClusterIndex[i])<=Interval && i+1<ClusterIndex.size()){
	std::vector< Int_t > AppliedClusterIndex;
	AppliedClusterIndex.push_back(std::get<0>(ClusterIndex[i]));
	/*
	std::cout << std::get<0>(ClusterIndex[i]) << " " << std::get<1>(ClusterIndex[i]) << "  ";
	*/
	while (std::get<0>(ClusterIndex[i+1])-std::get<0>(ClusterIndex[i])<=Interval && i+1<ClusterIndex.size()){
	  AppliedClusterIndex.push_back(std::get<0>(ClusterIndex[i+1]));
	  i++;
	  /*
	  std::cout << std::get<0>(ClusterIndex[i]) << " " << std::get<1>(ClusterIndex[i]) << "  ";
	  */
	}		
	/*
	std::cout << " | "; 
	*/
	CRKIntervalAppliedCluster.push_back(AppliedClusterIndex);
      }
    }
    /*
    std::cout << std::endl;
    */

    std::pair < std::vector <Int_t>, std::vector<Int_t> > ShiftAppliedPair;
    ShiftAppliedPair = ApplyShiftCond(STLIntervalAppliedCluster,CRKIntervalAppliedCluster, tolerance);

    /*
    std::cout << "Pair: ";
    for (int i=0; i<(ShiftAppliedPair.first).size(); i++){
      std::cout << (ShiftAppliedPair.first)[i] << "  ";
    }
    std::cout << " | ";
    for (int i=0; i<(ShiftAppliedPair.second).size(); i++){
      std::cout << (ShiftAppliedPair.second)[i] << "  ";
    }
    std::cout << " | " << std::endl;      
    */
    if (!ShiftAppliedPair.first.empty() && !ShiftAppliedPair.second.empty()){
      std::tuple <std::vector <Int_t>, std::vector <Int_t>, Double_t > tmp_tuple = std::make_tuple(ShiftAppliedPair.first, ShiftAppliedPair.second, TrigTime);
      ShiftAppliedCluster.push_back(tmp_tuple);}
  }  
  return ShiftAppliedCluster;
}

std::pair < std::vector <Int_t>, std::vector<Int_t> > ApplyShiftCond(std::vector < std::vector<Int_t> > STLIntervalAppliedCluster, std::vector < std::vector<Int_t> > CRKIntervalAppliedCluster, Int_t tolerance){
  std::pair<std::vector <Int_t>, std::vector <Int_t> >  ShiftAppliedPair;
  std::vector< std::pair<std::vector <Int_t>, std::vector <Int_t> > > PairCandidates;
  if (!CRKIntervalAppliedCluster.empty() && !STLIntervalAppliedCluster.empty()){

    for (Int_t i=0; i<CRKIntervalAppliedCluster.size(); i++){
      for (Int_t j=0; j<STLIntervalAppliedCluster.size(); j++){

        for (Int_t shift=-tolerance; shift<=tolerance; shift++){

          Int_t NumOfOverlap=0;
	  std::vector<Int_t> TempCluster;

          for (Int_t i_clu=0; i_clu < STLIntervalAppliedCluster[j].size(); i_clu++){
            if(std::find(CRKIntervalAppliedCluster[i].begin(),CRKIntervalAppliedCluster[i].end(),STLIntervalAppliedCluster[j][i_clu]+shift) != CRKIntervalAppliedCluster[i]\
	       .end()){
              NumOfOverlap++;
            }
          }

          if (NumOfOverlap >= 2){
            PairCandidates.push_back(std::make_pair(STLIntervalAppliedCluster[j],CRKIntervalAppliedCluster[i]));
          }

        }


      }
    }
  }
  ///////// No Candidate selection Algorithm yet.../////////
  if (!PairCandidates.empty()){
    //if (PairCandidates.size()>1) std::cout << "Candidates are more than 1!!" << std::endl;
    ShiftAppliedPair = PairCandidates.front();
  }
  ///////////////////////////////////////////////
  return ShiftAppliedPair;
}




