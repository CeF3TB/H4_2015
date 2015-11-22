#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <sstream> 

#include "DrawTools.h"
#include "RecoTree.h"
#include "channelInfo.h"

#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooFFTConvPdf.h"
#include "RooPlot.h"
#include "RooCBShape.h"
#include "RooCruijff.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TVectorD.h"
#include "TF1.h"

int main( int argc, char* argv[] ) {

  DrawTools::setStyle();



   if( argc<2 ) {
     std::cout << "ERROR. You need to specify the name of the run you want to process." << std::endl;  
     exit(1);
   }

   std::string runName = "";
   std::string tag = "V00";

   if( argc>1 ) {
     std::string runName_str(argv[1]);
     runName = runName_str;
     if( argc>2 ) {
       std::string tag_str(argv[2]);
       tag = tag_str;
     }
   } else {
     std::cout << "Usage:" << std::endl;
     std::cout << "./drawCherenkov [runName] ([tag])" << std::endl;
     exit(12345);
   }

  TString runNumberString(runName);

  std::string fileName;
  fileName = "analysisTrees_"+tag+"/Reco_"+runName+".root";

  TFile* file = TFile::Open(fileName.c_str());
  if( file==0 ) {
    std::cout << "ERROR! Din't find file " << fileName << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }

  TH1F* reso_histo = new TH1F("reso_histo","reso_histo",200,-2,12);

  
  TTree* recoTree=(TTree*)file->Get("recoTree");
  RecoTree t(recoTree);
  Long64_t nentries = t.fChain->GetEntries();
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = t.LoadTree(jentry);
      if (ientry < 0) break;
      nb = t.fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if( t.nino_LEtime<30 || t.nino_LEtime > 50 )continue;
      reso_histo->Fill(t.mcp_time_at_150-t.nino_LEtime);
   }

   float reso_mean = reso_histo->GetMean();
   float reso_sigma = reso_histo->GetRMS();



   TH1F* reso_histo_corr = new TH1F("reso_histo_corr","reso_histo",200,-5*reso_sigma,5*reso_sigma);

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = t.LoadTree(jentry);
      if (ientry < 0) break;
      nb = t.fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if( t.nino_LEtime<30 || t.nino_LEtime > 50)continue;
      reso_histo_corr->Fill(t.mcp_time_at_150-t.nino_LEtime-reso_mean);
   }



   // this is the dir in which the plots will be saved:
  std::string constDirName = "plots_timing";
  system(Form("mkdir -p %s", constDirName.c_str()));
  TString dir(constDirName);

  TCanvas c1;
  reso_histo->Draw();
  c1.SaveAs(dir+"/reso_histo_"+runNumberString+".png");
  c1.SaveAs(dir+"/reso_histo_"+runNumberString+".pdf");

  c1.Clear();
  TF1 *fb = new TF1("fb","[0]*exp(-0.5*((x-[1])/[2])**2)",-1,1.5);
  fb->SetParameter(0,reso_histo_corr->Integral());
  fb->SetParameter(1,0.1);
  fb->SetParameter(2,2);
  reso_histo_corr->Fit(fb,"RN");
  fb->SetLineWidth(2.);
  fb->SetLineColor(kBlue);
  reso_histo_corr->GetXaxis()->SetTitle("time_{nino}-time_{mcp} [ns]");
  float binWidth =reso_histo_corr->GetXaxis()->GetBinWidth(1);
  //  std::string ytitle = Form("Events / %.0f ps",binWidth*1.e3); 
  std::string ytitle = Form("Events");
  reso_histo_corr->SetYTitle(ytitle.c_str());
  reso_histo_corr->Draw();
  fb->Draw("L same");

  TLegend* leg_gauss = new TLegend(0.57, 0.7, 0.8, 0.92);
  leg_gauss->SetTextSize(0.038);
  leg_gauss->AddEntry(  (TObject*)0 ,Form("#sigma = %.0f #pm %.0f ps", fb->GetParameter(2)*1.e3, fb->GetParError(2)*1.e3), "");
  leg_gauss->SetFillColor(0);
  leg_gauss->Draw("same");



  c1.SaveAs(dir+"/reso_histo_corr_gaus_"+runNumberString+".png");
  c1.SaveAs(dir+"/reso_histo_corr_gaus_"+runNumberString+".pdf");


  //fit with cruijff and plot reso
  TH1F* histo;
  histo=reso_histo_corr;
  double peakpos = histo->GetBinCenter(histo->GetMaximumBin());
  double sigma = histo->GetRMS();

  double fitmin;
  double fitmax;
  
  
  fitmin = peakpos-4*sigma;
  fitmax = peakpos+4*sigma;
  
  RooRealVar x("x","deltaT", fitmin, fitmax);
  RooDataHist data("data","dataset with x",x,RooFit::Import(*histo) );

  RooRealVar meanr("meanr","Mean",peakpos-sigma,peakpos-3*sigma, peakpos+3*sigma);
  RooRealVar widthL("widthL","#sigmaL",sigma , 0, 5*sigma);
  RooRealVar widthR("widthR","#sigmaR",sigma , 0, 5*sigma);
  RooRealVar alphaL("alphaL","#alpha",5.08615e-02 , 0., 1.);
  RooRealVar alphaR("alphaR","#alpha",5.08615e-02, 0., 1.);
  int ndf;

  RooPlot* frame;

  RooCruijff fit_fct("fit_fct","fit_fct",x,meanr,widthL,widthR,alphaL,alphaR); ndf = 5;
  fit_fct.fitTo(data);
  
  frame = x.frame("Title");
  frame->SetXTitle("time_{nino}-time_{mcp} [ns]");
  frame->SetYTitle(ytitle.c_str());
  
  data.plotOn(frame);  //this will show histogram data points on canvas 
  fit_fct.plotOn(frame);//this will show fit overlay on canvas  

  double rms,rmsErr;
  rms = (widthL.getVal()+widthR.getVal())/2;
  rmsErr = 0.5*sqrt(widthL.getError()*widthL.getError()+widthR.getError()*widthR.getError());

  TCanvas* cans = new TCanvas();
  cans->cd();
  frame->Draw();
  TLegend* lego = new TLegend(0.57, 0.7, 0.8, 0.92);
  lego->SetTextSize(0.038);
  lego->AddEntry(  (TObject*)0 ,Form("#sigma = %.1f #pm %.1f ps", rms*1.e3, rmsErr*1.e3), "");

  lego->SetFillColor(0);
  lego->Draw("same");
  
  std::string energy(Form("%.0f", t.beamEnergy));
  TPaveText* pave = DrawTools::getLabelTop_expOnXaxis(energy+" GeV Electron Beam");
  pave->Draw("same");

  cans->SaveAs(dir+"/reso_histo_corr"+runNumberString+".png");
  cans->SaveAs(dir+"/reso_histo_corr"+runNumberString+".pdf");


}
