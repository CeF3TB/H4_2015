#define WaveformUtil_cxx
#include "WaveformUtil.h"
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


void WaveformUtil::Loop(){

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();

  float  mean[NFIBERS][NDIGISAMPLES];
  float  time[NDIGISAMPLES];

  for (int i=0;i<NFIBERS;++i){
    for (int j=0;j<NDIGISAMPLES;++j){
      mean[i][j]=0;
    }
  }

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;      
    for (int i=0;i<1024*4;++i){
      //	std::cout<<"channel:"<<digi_value_ch->at(i)<<digi_value->at(i)<<" "<<digi_value_time->at(i)<<std::endl;
      if(digi_value_ch->at(i) > 3)continue;
      mean[digi_value_ch->at(i)][i-1024*digi_value_ch->at(i)]+=(float)(digi_value->at(i)/nentries);
      if(i<1024)time[i]=digi_value_time->at(i);
    }
  }

//  for (int i=0;i<NFIBERS;++i){
//    for (int j=0;j<NDIGISAMPLES;++j){
//      //std::cout<<"(i,j):"<<i<<","<<j<<" "<<mean[i][j]/nentries<<std::endl;
//      mean[i][j]/=nentries;
//    }
//  }



  TFile* outFile = TFile::Open("outWaveform.root","recreate");

  TRandom3 *r = new TRandom3(0);

  for (int i=0;i<NFIBERS;++i){
    meanWaveGraphs[i]=new TGraph(1024, time, mean[i]);
    TString fiber;
    fiber.Form("%d",i);
    meanWaveGraphs[i]->SetName("waveform_"+fiber);

    meanWaveHistos[i]=new TH1F("waveform_histo_"+fiber,"",1024,0,time[1023]);
    for (int j=0;j<NDIGISAMPLES;++j){  
      for(int k=0;k<mean[i][j]*nentries;k++){//this nentries here is just to get correct errors, slow method
	meanWaveHistos[i]->Fill(time[j]);
      }
    }
    meanWaveHistos[i]->Sumw2();
    meanWaveHistos[i]->Scale(1./nentries);
    meanWaveHistos[i]->Write();
    meanWaveGraphs[i]->Write();

  }


    RooRealVar t("t","time", time[0], time[1023]);
    RooDataHist data("data","dataset with t",t,RooFit::Import(*meanWaveHistos[1]) );
    //    RooDataSet data("data","dataset with t",t);



    // Construct landau(t,ml,sl) ;
    RooRealVar ml("ml","mean landau",0.05e-6,0,0.12e-6) ;
    RooRealVar sl("sl","sigma landau",0.05e-6,0.,0.12e-6) ;
    RooLandau landau("lx","lx",t,ml,sl) ;
  
    // Construct gauss(t,mg,sg)
    RooRealVar mg("mg","mg",1.06706e-07,0.,0.12e-6) ;
    RooRealVar sg("sg","sg",3.91767e-08,0.,0.12e-6) ;
    RooGaussian gauss("gauss","gauss",t,mg,sg) ;

    // Construct landau (x) gauss
    RooFFTConvPdf lxg("lxg","landau (X) gauss",t,landau,gauss) ;

    // Fit gxlx to data
    lxg.fitTo(data,RooFit::Range(time[100],time[1000])) ;
    //    landau.fitTo(data) ;

    // Plot data, landau pdf, landau (X) gauss pdf
    RooPlot* frame = t.frame(RooFit::Title("")) ;
    t.setBins(10);
    data.plotOn(frame,RooFit::DataError(RooAbsData::SumW2)) ;

    lxg.plotOn(frame) ;
    //landau.plotOn(frame) ;
    landau.plotOn(frame,RooFit::LineStyle(kDashed)) ;

    TCanvas c1;
    frame->Draw();
    c1.Write();



  outFile->Write();
  outFile->Close();

}




