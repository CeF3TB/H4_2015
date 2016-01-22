#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <sstream> 

#include "DrawTools.h"
#include "dummyRecoTree.h"
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

#include "interface/Configurator.h"
#include "plotOverallPerformancesConfigurator.h"

plotOverallPerformances_Config_t readConfiguration(std::string configName);
bool passesFibreTopologicalSelection(dummyRecoTree &t);
bool passesChannelTopologicalSelection(dummyRecoTree &t);


int main( int argc, char* argv[] ) {

  DrawTools::setStyle();
  
  
  std::string tag = "V00";
  std::string config = "SiPM2015Config";
  
  if( argc>1 ) {
    std::string tag_str(argv[1]);
    tag = tag_str;
    if( argc>2 ) {
      std::string config_str(argv[2]);
      config=config_str;
    }
  } else {
    std::cout << "Usage:" << std::endl;
    std::cout << "./plotOverallPerformances ([tag]) ([config])" << std::endl;
    exit(12345);
  }

  std::cout<<config<<std::endl;
  
  theConfiguration_=readConfiguration(config);
  
  TH1F* reso_histo_fibre_corr_Amplitude[theConfiguration_.nMaxAmplCuts];//to see behaviour of resolution as a function of ampl
   TH1F* reso_histo_channel_corr_Amplitude[theConfiguration_.nMaxAmplCuts];
   TVectorD resValueAmplitude_fibre(theConfiguration_.nMaxAmplCuts);
   TVectorD resErrValueAmplitude_fibre(theConfiguration_.nMaxAmplCuts);

   TVectorD resValueAmplitude_channel(theConfiguration_.nMaxAmplCuts);
   TVectorD resErrValueAmplitude_channel(theConfiguration_.nMaxAmplCuts);


   for (int i=0;i<theConfiguration_.nMaxAmplCuts;++i){
    TString icut;
    icut.Form("%d",i); 
    reso_histo_fibre_corr_Amplitude[i] = new TH1F("reso_histo_fibre_corr_Amplitude_"+icut,"reso_histo_fibre_corr_Amplitude_"+icut,200,theConfiguration_.rangeXLow,theConfiguration_.rangeXUp);

    reso_histo_channel_corr_Amplitude[i] = new TH1F("reso_histo_channel_corr_Amplitude_"+icut,"reso_histo_channel_corr_Amplitude_"+icut,200,theConfiguration_.rangeXLow,theConfiguration_.rangeXUp);

   }


  std::vector<int> runs;

  TString haddFile = "analysisTrees_"+tag+"/"+theConfiguration_.setup+"_runs.root";

  TString haddString = "hadd analysisTrees_"+tag+"/"+theConfiguration_.setup+"_runs.root ";

  for (int i=0;i<theConfiguration_.runs.size();++i){
    TString fileNameString = "analysisTrees_"+tag+"/Reco_";
    runs.push_back(theConfiguration_.runs[i]);
    fileNameString+=runs[i];
    haddString+=fileNameString;
    haddString+= ".root ";
  }

  system(haddString.Data());

  TFile* file = TFile::Open(haddFile.Data());

  TTree* recoTree=(TTree*)file->Get("recoTree");
  dummyRecoTree t(recoTree);
  Long64_t nentries = t.fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = t.LoadTree(jentry);
      if (ientry < 0) break;
      nb = t.fChain->GetEntry(jentry);   nbytes += nb;
      if( t.mcp_time_frac50<0 ||t.mcp_max_amplitude<200)continue;
      for (int i=0;i<theConfiguration_.nMaxAmplCuts;++i){
      if(passesFibreTopologicalSelection(t) &&  t.cef3_maxAmpl->at(2) < theConfiguration_.channel2CutFibre){ 
	  if(t.cef3_maxAmpl->at(1)>i*theConfiguration_.stepAmplFibre && t.cef3_maxAmpl->at(1)<(i+1)*theConfiguration_.stepAmplFibre){
	    reso_histo_fibre_corr_Amplitude[i]->Fill(t.cef3_time_at_frac50->at(1)-t.mcp_time_frac50);
	    break;
	  }
      }
	  if(passesChannelTopologicalSelection(t)  && t.cef3_maxAmpl->at(2) > theConfiguration_.channel2CutChannel && t.cef3_maxAmpl->at(1)<theConfiguration_.channel1CutChannel){   //cuts on fibre position with hodos and wc. cut on cef3_maxAmpl[2] to reduce hadron contamination, cef3_maxAmpl[1] cut to reduce remaining events hitting the fibre
	    if(t.cef3_maxAmpl->at(1)>i*theConfiguration_.stepAmplChannel && t.cef3_maxAmpl->at(1)<(i+1)*theConfiguration_.stepAmplChannel){
	    reso_histo_channel_corr_Amplitude[i]->Fill(t.cef3_time_at_frac50->at(1)-t.mcp_time_frac50);
	    break;
	    }
	  }
      }
      
   }


   // this is the dir in which the plots will be saved:
   std::string constDirName = "plots_timingPerformance_";
   constDirName+=theConfiguration_.setup;
   system(Form("mkdir -p %s", constDirName.c_str()));
   TString dir(constDirName);
   
   //fit histos
   for (int i=0;i<theConfiguration_.nMaxAmplCuts;++i){
     //fibre
     if(reso_histo_fibre_corr_Amplitude[i]->GetEntries()>=20){
     TH1F* histo;
     histo=reso_histo_fibre_corr_Amplitude[i];
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
    
     std::string ytitle = Form("Events");
     frame = x.frame("Title");
     frame->SetXTitle("time_{fibre}-time_{mcp} [ns]");
     frame->SetYTitle(ytitle.c_str());
     
     data.plotOn(frame);  //this will show histogram data points on canvas 
     fit_fct.plotOn(frame);//this will show fit overlay on canvas  

     double rms,rmsErr;
     rms = (widthL.getVal()+widthR.getVal())/2;
     rmsErr = 0.5*sqrt(widthL.getError()*widthL.getError()+widthR.getError()*widthR.getError());

     resValueAmplitude_fibre[i]=rms*1.e3; 
     resErrValueAmplitude_fibre[i]=rmsErr*1.e3; 


     TCanvas* cans = new TCanvas();
     cans->cd();
     frame->Draw();
     TLegend* lego = new TLegend(0.6, 0.8, 0.8, 0.9);
     lego->SetTextSize(0.036);
     lego->AddEntry(  (TObject*)0 ,Form("#sigma = %.1f #pm %.1f ps", rms*1.e3, rmsErr*1.e3), "");
     
     lego->SetFillColor(0);
     lego->Draw("same");
     
     TPaveText* pave = DrawTools::getLabelTop_expOnXaxis("Electron Beam");
     pave->Draw("same");
     
     TString icut;
     icut.Form("%d",i); 

     cans->SaveAs(dir+"/reso_histo_fibre_cruijff_"+icut+".png");
     cans->SaveAs(dir+"/reso_histo_fibre_cruijff_"+icut+".pdf");
     }

     //channel
     if(reso_histo_channel_corr_Amplitude[i]->GetEntries()>15){
       TH1F* histo;
       histo=reso_histo_channel_corr_Amplitude[i];
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
    
       std::string ytitle = Form("Events");
       frame = x.frame("Title");
       frame->SetXTitle("time_{fibre}-time_{mcp} [ns]");
       frame->SetYTitle(ytitle.c_str());
     
       data.plotOn(frame);  //this will show histogram data points on canvas 
       fit_fct.plotOn(frame);//this will show fit overlay on canvas  

       double rms,rmsErr;
       rms = (widthL.getVal()+widthR.getVal())/2;
       rmsErr = 0.5*sqrt(widthL.getError()*widthL.getError()+widthR.getError()*widthR.getError());

       resValueAmplitude_channel[i]=rms*1.e3; 
       resErrValueAmplitude_channel[i]=rmsErr*1.e3; 


       TCanvas* cans = new TCanvas();
       cans->cd();
       frame->Draw();
       TLegend* lego = new TLegend(0.6, 0.8, 0.8, 0.9);
       lego->SetTextSize(0.036);
       lego->AddEntry(  (TObject*)0 ,Form("#sigma = %.1f #pm %.1f ps", rms*1.e3, rmsErr*1.e3), "");
     
       lego->SetFillColor(0);
       lego->Draw("same");
     
       TPaveText* pave = DrawTools::getLabelTop_expOnXaxis("Electron Beam");
       pave->Draw("same");
     
       TString icut;
       icut.Form("%d",i); 

       cans->SaveAs(dir+"/reso_histo_channel_cruijff_"+icut+".png");
       cans->SaveAs(dir+"/reso_histo_channel_cruijff_"+icut+".pdf");
     }

   }
	  


  return 0;

}


plotOverallPerformances_Config_t readConfiguration(std::string configName){

  Configurator* configurator_ = new Configurator();
  std::string fileName = Form ("./config_timing/plots_performances_%s.xml",configName.c_str());
  configurator_->xmlFileName=fileName.c_str();
  configurator_->Init();

  plotOverallPerformances_Config_t conf;
 
  conf.setup=Configurable::getElementContent (*configurator_, "setup",configurator_->root_element) ;

  string  content = Configurable::getElementContent (*configurator_, "runs",configurator_->root_element) ;
  Configurator::GetVecInt(content,conf.runs);

   conf.nMaxAmplCuts= Configurator::GetInt(Configurable::getElementContent(*configurator_,"nMaxAmplCuts",configurator_->root_element));

   conf.stepAmplFibre= Configurator::GetInt(Configurable::getElementContent(*configurator_,"stepAmplFibre",configurator_->root_element));

   conf.stepAmplChannel= Configurator::GetInt(Configurable::getElementContent(*configurator_,"stepAmplChannel",configurator_->root_element));

   conf.rangeXLow= Configurator::GetDouble(Configurable::getElementContent(*configurator_,"rangeXLow",configurator_->root_element));

   conf.rangeXUp= Configurator::GetDouble(Configurable::getElementContent(*configurator_,"rangeXUp",configurator_->root_element));

   conf.channel2CutFibre= Configurator::GetDouble(Configurable::getElementContent(*configurator_,"channel2CutFibre",configurator_->root_element));
   conf.channel2CutChannel= Configurator::GetDouble(Configurable::getElementContent(*configurator_,"channel2CutChannel",configurator_->root_element));
   conf.channel1CutChannel= Configurator::GetDouble(Configurable::getElementContent(*configurator_,"channel1CutChannel",configurator_->root_element));

   conf.nBinsFibre= Configurator::GetInt(Configurable::getElementContent(*configurator_,"nBinsFibre",configurator_->root_element));
   conf.nBinsChannel= Configurator::GetInt(Configurable::getElementContent(*configurator_,"nBinsChannel",configurator_->root_element));


   return conf;

}

bool passesFibreTopologicalSelection(dummyRecoTree &t){
  float  X=(t.cluster_pos_corr_hodoY1+t.cluster_pos_corr_hodoY2+t.wc_x_corr)/3.;
  float  Y=(t.cluster_pos_corr_hodoY1+t.cluster_pos_corr_hodoY2+t.wc_y_corr)/3.;

  if(Y>-3.5 && Y<-0.5 && X<-2.5 && X>-6 && t.nClusters_hodoX1>0 && t.cluster_pos_corr_hodoX1>-100 && t.cluster_pos_corr_hodoY1>-100&& t.cluster_pos_corr_hodoX2>-100 && t.cluster_pos_corr_hodoY2>-100 && t.wc_x_corr>-20){
    if(0.66*X+4.16+Y<0){//line passing through (-5.5,-0.5) and (-2.5,2.5)
      return true; 
    }
}

  return false;
}

bool passesChannelTopologicalSelection(dummyRecoTree &t){
  float  X=(t.cluster_pos_corr_hodoY1+t.cluster_pos_corr_hodoY2+t.wc_x_corr)/3.;
  float  Y=(t.cluster_pos_corr_hodoY1+t.cluster_pos_corr_hodoY2+t.wc_y_corr)/3.;

  if(Y>-3 && Y<5.5 && X<5.5 && X>-6 && t.nClusters_hodoX1>0 && t.cluster_pos_corr_hodoX1>-100 && t.cluster_pos_corr_hodoY1>-100&& t.cluster_pos_corr_hodoX2>-100 && t.cluster_pos_corr_hodoY2>-100 && t.wc_x_corr>-20){
    if(0.66*X+4.16+Y>0){//line passing through (-5.5,-0.5) and (-2.5,2.5)
    return true; 
  }
}
}
