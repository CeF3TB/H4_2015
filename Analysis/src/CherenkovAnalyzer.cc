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
  outputTree->Branch( "nChan", &nChan, "nChan/I" );
  outputTree->Branch( "digi_charge_integrated_cher", digi_charge_integrated_cher, "digi_charge_integrated_cher[nChan]/F" );
  outputTree->Branch( "digi_charge_integrated_wls", digi_charge_integrated_wls, "digi_charge_integrated_wls[nChan]/F" );


  TVectorD integrals(4);
  integrals[0]=0.0707221;
  integrals[1]=0.0578497;
  integrals[2]=0.346977;
  integrals[3]=0.0622282;


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
	outputTree->Fill();
      }

  }
  TFile* outFile = TFile::Open("CherenkovAnalyzer_"+runNumberString+".root","recreate");
  outputTree->Write();
  outFile->Write("outCherenkovAnalyzer_"+runNumberString+".root");
  outFile->Close();

}
