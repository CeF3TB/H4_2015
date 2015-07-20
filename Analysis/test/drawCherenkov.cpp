#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>


#include "DrawTools.h"
#include "RecoTree.h"
#include "channelInfo.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TGaxis.h"

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

  std::string fileName = "analysisTrees_"+tag+"/Reco_"+runName+".root";
  TFile* file = TFile::Open(fileName.c_str());
  if( file==0 ) {
    std::cout << "ERROR! Din't find file " << fileName << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }
  


  TH1F* totalHistos[4];
  TH1F* totalHistosGainCorr[4];
  TH1F* cherHistos[4];
  TH1F* wlsHistos[4];

  for(int i=0;i<CEF3_CHANNELS;++i){
    TString fibre;
    fibre.Form("%d",i); 
    totalHistosGainCorr[i] = new TH1F ("totalHistoGainCorr_"+fibre,"",200,0,0.5);
    totalHistos[i] = new TH1F ("totalHisto_"+fibre,"",200,0,300000);
    cherHistos[i] = new TH1F ("cherHisto_"+fibre,"",200,0,55e3);
    wlsHistos[i] = new TH1F ("wlsHisto_"+fibre,"",200,0,500000);
  }

  TTree* recoTree=(TTree*)file->Get("recoTree");
  RecoTree t(recoTree);

   if (t.fChain == 0) return 0;
   
   Long64_t nentries = t.fChain->GetEntries();

   float gainR1450=1.5e+6;
   float gainR5380=7e+4;
   
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = t.LoadTree(jentry);
      if (ientry < 0) break;
      nb = t.fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      for(int i=0;i<CEF3_CHANNELS;++i){
	totalHistos[i]->Fill(t.cef3_chaInt->at(i));
	if(i==2)	totalHistosGainCorr[i]->Fill(t.cef3_chaInt->at(i)/gainR5380);
	else totalHistosGainCorr[i]->Fill(t.cef3_chaInt->at(i)/gainR1450);
	cherHistos[i]->Fill(t.cef3_chaInt_cher->at(i));
	wlsHistos[i]->Fill(t.cef3_chaInt_wls->at(i));
      }
   }




  TGaxis::SetMaxDigits(3);



  TFile* outFile = TFile::Open("CherenkovPlots_"+runNumberString+".root","recreate");

  //plots of ch int
  TPaveText * pave= new TPaveText(0.75,0.65,0.85,0.85,"NDC");
  pave->SetFillColor(kWhite);
  pave->SetTextSize(0.030);
  pave->SetTextAlign(kHAlignLeft);
  pave->SetTextFont(62);

  TCanvas c1;
  c1.SetLogy();
  totalHistos[1]->SetLineColor(kBlue);
  totalHistos[2]->SetLineColor(kRed);
  totalHistos[3]->SetLineColor(kViolet);

  float max=-1;
  for(int i=0;i<4;++i){
      if(max<totalHistos[i]->GetMaximum()){
      max=totalHistos[i]->GetMaximum();
    }
  }


  for(int i=0;i<4;++i){
    totalHistos[i]->GetYaxis()->SetRangeUser(1,max*1.1);
    totalHistos[i]->GetXaxis()->SetTitle("Charge Integrated");
    totalHistos[i]->GetYaxis()->SetTitle("Events");
    if (i==0)    totalHistos[i]->Draw();
    else totalHistos[i]->Draw("same");
    TString fibre;
    fibre.Form("%d",i); 
    pave->AddText("Fibre "+fibre);  
    pave->SetAllWith("Fibre "+fibre,"color",totalHistos[i]->GetLineColor());
  }
  pave->Draw("same");
  c1.SaveAs("plots_drawCherenkov/chargeIntTotal_"+runNumberString+".png");
  c1.SaveAs("plots_drawCherenkov/chargeIntTotal_"+runNumberString+".pdf");
  c1.Write("chargeIntTotal");


  c1.Clear();
  totalHistosGainCorr[1]->SetLineColor(kBlue);
  totalHistosGainCorr[2]->SetLineColor(kRed);
  totalHistosGainCorr[3]->SetLineColor(kViolet);
  max=-1;
  for(int i=0;i<4;++i){
    if(max<totalHistosGainCorr[i]->GetMaximum()){
      max=totalHistosGainCorr[i]->GetMaximum();
    }
  }
  
  for(int i=0;i<4;++i){
    totalHistosGainCorr[i]->GetYaxis()->SetRangeUser(1,max*1.1);
    totalHistosGainCorr[i]->GetXaxis()->SetTitle("Charge Integrated");
    totalHistosGainCorr[i]->GetYaxis()->SetTitle("Events");
    if (i==0)    totalHistosGainCorr[i]->Draw();
    else totalHistosGainCorr[i]->Draw("same");
  }
  pave->Draw("same");
  c1.SaveAs("plots_drawCherenkov/chargeIntTotalGainCorr_"+runNumberString+".png");
  c1.SaveAs("plots_drawCherenkov/chargeIntTotalGainCorr_"+runNumberString+".pdf");
  c1.Write("chargeIntTotal");


  c1.Clear();
  cherHistos[1]->SetLineColor(kBlue);
  cherHistos[2]->SetLineColor(kRed);
  cherHistos[3]->SetLineColor(kViolet);

  max=-1;
  for(int i=0;i<4;++i){
    if(max<cherHistos[i]->GetMaximum()){
      max=cherHistos[i]->GetMaximum();
    }
  }


  for(int i=0;i<4;++i){
    cherHistos[i]->GetYaxis()->SetRangeUser(1,max*1.1);
    cherHistos[i]->GetXaxis()->SetTitle("Charge Integrated");
    cherHistos[i]->GetYaxis()->SetTitle("Events");
    if (i==0)    cherHistos[i]->Draw();
    else cherHistos[i]->Draw("same");
  }
  pave->Draw("same");
  c1.SaveAs("plots_drawCherenkov/chargeIntCher_"+runNumberString+".png");
  c1.SaveAs("plots_drawCherenkov/chargeIntCher_"+runNumberString+".pdf");
  c1.Write("chargeIntCher");

  c1.Clear();
  wlsHistos[1]->SetLineColor(kBlue);
  wlsHistos[2]->SetLineColor(kRed);
  wlsHistos[3]->SetLineColor(kViolet);
  max=-1;
  for(int i=0;i<4;++i){
    if(max<wlsHistos[i]->GetMaximum()){
      max=wlsHistos[i]->GetMaximum();
    }
  }

  for(int i=0;i<4;++i){
    wlsHistos[i]->GetYaxis()->SetRangeUser(1,max*1.1);
    wlsHistos[i]->GetXaxis()->SetTitle("Charge Integrated");
    wlsHistos[i]->GetYaxis()->SetTitle("Events");
    if (i==0)    wlsHistos[i]->Draw();
    else wlsHistos[i]->Draw("same");
  }
  pave->Draw("same");
  c1.SaveAs("plots_drawCherenkov/chargeIntWls_"+runNumberString+".png");
  c1.SaveAs("plots_drawCherenkov/chargeIntWls_"+runNumberString+".pdf");
  c1.Write("chargeIntWls");


  for(int i=0;i<4;++i){
    totalHistos[i]->Write();
    totalHistosGainCorr[i]->Write();
    cherHistos[i] ->Write();
    wlsHistos[i] ->Write();
  }

  outFile->Write();
  outFile->Close();



  return 0;

}
