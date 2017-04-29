#include "TCanvas.h"
#include "TGraph.h"
#include "TPad.h"
#include "TAxis.h"
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
#include <iostream>
#include <sstream>
#include <TMath.h>

void SignificanceLevel(){
  
  Double_t OnShellPlusOffShell[8]={0,0.876,1.911,2.957,3.965,5.131,6.091,7.003};
  Double_t OnShellOnly[8]={0.319,1.791,3.373,4.877,6.375,7.610,9.083,10.631};
  Double_t SigBR[8]={1,2,3,4,5,6,7,8};

  TGraph *gr = new TGraph(8,SigBR,OnShellOnly);
  TGraph *gr2 = new TGraph(8,SigBR,OnShellPlusOffShell);
  gr->SetMarkerStyle(20);
  gr2->SetMarkerStyle(20);
  gr->Draw("APL");

  gr2->Draw("same");

}
