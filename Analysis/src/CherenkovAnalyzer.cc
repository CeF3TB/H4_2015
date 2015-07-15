#define CherenkovAnalyzer_cxx
#include "CherenkovAnalyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TRandom3.h>
#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooFFTConvPdf.h"
#include "RooPlot.h"
#include "TVectorD.h"
#include "TPaveText.h"
#include "TGaxis.h"

void CherenkovAnalyzer::Loop(){

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();
  std::cout<<"nentries"<<nentries<<std::endl;
  //    nentries=7000;

  TString runNumberString;
  int digiFreq;


  gROOT->cd();
  TTree* outputTree=new TTree("recoTree","recoTree");
  int nChan=4;
  float digi_charge_integrated_cher[nChan];
  float digi_charge_integrated_wls[nChan];
  float digi_charge_integrated_total[nChan];
  outputTree->Branch( "nChan", &nChan, "nChan/I" );
  outputTree->Branch( "digi_charge_integrated_total", digi_charge_integrated_total, "digi_charge_integrated_total[nChan]/F" );
  outputTree->Branch( "digi_charge_integrated_cher", digi_charge_integrated_cher, "digi_charge_integrated_cher[nChan]/F" );
  outputTree->Branch( "digi_charge_integrated_wls", digi_charge_integrated_wls, "digi_charge_integrated_wls[nChan]/F" );

  TH1F* totalHistos[4];
  TH1F* totalHistosGainCorr[4];
  TH1F* cherHistos[4];
  TH1F* wlsHistos[4];

  for(int i=0;i<4;++i){
    TString fibre;
    fibre.Form("%d",i); 
    totalHistosGainCorr[i] = new TH1F ("totalHistoGainCorr_"+fibre,"",200,0,0.5);
    totalHistos[i] = new TH1F ("totalHisto_"+fibre,"",200,0,300000);
    cherHistos[i] = new TH1F ("cherHisto_"+fibre,"",200,0,55e3);
    wlsHistos[i] = new TH1F ("wlsHisto_"+fibre,"",200,0,500000);
  }
  TGaxis::SetMaxDigits(3);

  TVectorD integrals(4);
  integrals[0]=0.0707221;
  integrals[1]=0.0578497;
  integrals[2]=0.346977;
  integrals[3]=0.0622282;

  float gainR1450=1.5e+6;
  float gainR5380=7e+4;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;      
    if(jentry==0){
      runNumberString.Form("%d",runNumber);
      digiFreq=digi_frequency;
    }
    if(jentry%1000 == 0)std::cout<<"Processing entry:"<<jentry<<std::endl;
      for (int i=0;i<NFIBERS;++i){
	if(digi_max_amplitude->at(i)>10000 || digi_max_amplitude->at(i)<0)continue;
	float chargeWlsUnderCher=(digi_charge_integrated_bare_slow->at(i)*integrals[i])/(1-integrals[i]);//amount of wls charge under cherenkov peak. estimated by fit
	digi_charge_integrated_cher[i]=digi_charge_integrated_bare_fast->at(i) - chargeWlsUnderCher;
	digi_charge_integrated_wls[i]=digi_charge_integrated_bare_slow->at(i) + chargeWlsUnderCher;
	digi_charge_integrated_total[i]=digi_charge_integrated_bare->at(i);
	totalHistos[i]->Fill(digi_charge_integrated_total[i]);
	if(i==2)	totalHistosGainCorr[i]->Fill(digi_charge_integrated_total[i]/gainR5380);
	else totalHistosGainCorr[i]->Fill(digi_charge_integrated_total[i]/gainR1450);
	cherHistos[i]->Fill(digi_charge_integrated_cher[i]);
	wlsHistos[i]->Fill(digi_charge_integrated_wls[i]);
	outputTree->Fill();
      }

  }

  TFile* outFile = TFile::Open("CherenkovAnalyzer_"+runNumberString+".root","recreate");
  outputTree->Write();

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
  c1.SaveAs("plots_analyzeCherenkov/chargeIntTotal_"+runNumberString+".png");
  c1.SaveAs("plots_analyzeCherenkov/chargeIntTotal_"+runNumberString+".pdf");
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
  c1.SaveAs("plots_analyzeCherenkov/chargeIntTotalGainCorr_"+runNumberString+".png");
  c1.SaveAs("plots_analyzeCherenkov/chargeIntTotalGainCorr_"+runNumberString+".pdf");
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
  c1.SaveAs("plots_analyzeCherenkov/chargeIntCher_"+runNumberString+".png");
  c1.SaveAs("plots_analyzeCherenkov/chargeIntCher_"+runNumberString+".pdf");
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
  c1.SaveAs("plots_analyzeCherenkov/chargeIntWls_"+runNumberString+".png");
  c1.SaveAs("plots_analyzeCherenkov/chargeIntWls_"+runNumberString+".pdf");
  c1.Write("chargeIntWls");


  for(int i=0;i<4;++i){
    totalHistos[i]->Write();
    totalHistosGainCorr[i]->Write();
    cherHistos[i] ->Write();
    wlsHistos[i] ->Write();
  }

  outFile->Write();
  outFile->Close();

}
