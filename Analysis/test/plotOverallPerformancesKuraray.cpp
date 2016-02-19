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
#include "RooGenericPdf.h"
#include "RooChi2Var.h"
#include "RooMinuit.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TVectorD.h"
#include "TF1.h"
#include "TGraphErrors.h"

#include "interface/Configurator.h"
#include "plotOverallPerformancesConfigurator.h"
#include "interface/TopologicalSelectionHelper.h"

plotOverallPerformances_Config_t readConfiguration(std::string configName);


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
    std::cout << "./plotOverallPerformancesKuraray ([tag]) ([config])" << std::endl;
    exit(12345);
  }

  std::cout<<config<<std::endl;
  
  theConfiguration_=readConfiguration(config);
  
  TH1F* reso_histo_channel_corr_Amplitude[theConfiguration_.nMaxAmplCuts];

  int isMAPDDiff=0;

  if(theConfiguration_.setup=="MAPDDiff"){
    MAPDDiff=1;
  }

  TH1F* maxAmplHisto=new TH1F("maxAmplHisto","maxAmplHisto",1000,0,1000+3000*isMAPDDiff);

  std::vector<TH1F*> reso_histo_channel_total;

  TVectorD resValueAmplitude_channel(theConfiguration_.nMaxAmplCuts);
  TVectorD resErrValueAmplitude_channel(theConfiguration_.nMaxAmplCuts);

  TVectorD resValueEnergy_channel(theConfiguration_.nMaxAmplCuts);
  TVectorD resErrValueEnergy_channel(theConfiguration_.nMaxAmplCuts);

  for (int i=0;i<theConfiguration_.nMaxAmplCuts;++i){
    TString icut;
    icut.Form("%d",i); 
    reso_histo_channel_corr_Amplitude[i] = new TH1F("reso_histo_channel_corr_Amplitude_"+icut,"reso_histo_channel_corr_Amplitude_"+icut,200,theConfiguration_.rangeXLow,theConfiguration_.rangeXUp);
  }
   




   
  std::vector<int> runs;
  std::vector<float> energies;

  TString haddFile = "analysisTrees_"+tag+"/"+theConfiguration_.setup+"_centralRuns.root";

  if(theConfiguration_.addTagFileName){
    haddFile = "analysisTrees_"+tag+"/"+theConfiguration_.setup+"_"+theConfiguration_.tagFileName+"_centralRuns.root";
  }

  TString haddString = "hadd analysisTrees_"+tag+"/"+theConfiguration_.setup+"_centralRuns.root ";

  if(theConfiguration_.addTagFileName){
    haddString = "hadd analysisTrees_"+tag+"/"+theConfiguration_.setup+"_"+theConfiguration_.tagFileName+"_centralRuns.root ";
  }

  for (int i=0;i<theConfiguration_.runs.size();++i){
    TString fileNameString = "analysisTrees_"+tag+"/Reco_";
    runs.push_back(theConfiguration_.runs[i]);
    energies.push_back(theConfiguration_.energies[i]);
    fileNameString+=runs[i];
    haddString+=fileNameString;
    haddString+= ".root ";
    TString icut;
    icut.Form("%d",(int)energies[i]); 
    reso_histo_channel_total.push_back(new TH1F("reso_histo_channel_total_"+icut,"reso_histo_channel_total_"+icut,200,theConfiguration_.rangeXLow,theConfiguration_.rangeXUp));
  }

  system(haddString.Data());

  TFile* file = TFile::Open(haddFile.Data());

  TTree* recoTree=(TTree*)file->Get("recoTree");
  RecoTree t(recoTree);
  Long64_t nentries = t.fChain->GetEntries();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = t.LoadTree(jentry);
    if (ientry < 0) break;
    nb = t.fChain->GetEntry(jentry);   nbytes += nb;
    if( t.mcp_time_frac50<0 ||t.mcp_max_amplitude<200)continue;
      
    maxAmplHisto->Fill(t.cef3_maxAmpl->at(2));

    for (int i=0;i<theConfiguration_.nMaxAmplCuts;++i){
      if(t.cef3_maxAmpl->at(2)<theConfiguration_.startCutChannel)continue;
	
      float deltaTNoCorr=t.cef3_time_at_frac50->at(2)-t.mcp_time_frac50;

      if(t.cef3_maxAmpl->at(2)>theConfiguration_.startCutChannel+i*theConfiguration_.stepAmplChannel && t.cef3_maxAmpl->at(2)<theConfiguration_.startCutChannel+(i+1)*theConfiguration_.stepAmplChannel){
	if((t.nClusters_hodoX1==1||t.pos_2FibClust_hodoX1>-999) && (t.nClusters_hodoX2==1||t.pos_2FibClust_hodoX2>-999) && (t.nClusters_hodoY1==1||t.pos_2FibClust_hodoY1>-999) && (t.nClusters_hodoY1==1||t.pos_2FibClust_hodoY1>-999)){//exactly one cluster, or, if there are multiple clusters, exactly one 2-fiber cluster
	  if(TMath::Abs(0.5* (t.cluster_pos_corr_hodoX1+t.cluster_pos_corr_hodoX2))< 3 && TMath::Abs( 0.5* (t.cluster_pos_corr_hodoY1+t.cluster_pos_corr_hodoY2))< 3 && (t.wc_x_corr-t.cluster_pos_corr_hodoX2)<4 && (t.wc_y_corr-t.cluster_pos_corr_hodoY2)< 4  && TMath::Abs( (t.cluster_pos_corr_hodoX1-t.cluster_pos_corr_hodoX2))<1.5 &&TMath::Abs( (t.cluster_pos_corr_hodoY1-t.cluster_pos_corr_hodoY2))<1.5){
	    
	    reso_histo_channel_corr_Amplitude[i]->Fill(deltaTNoCorr);
	    if(t.cef3_maxAmpl->at(2)>theConfiguration_.amplCut)	   {
	      for(int j=0;j<energies.size();j++){
		  if(t.beamEnergy==150 && t.cef3_maxAmpl->at(2)<theConfiguration_.channel2CutChannel)continue; //FIXME we use channel for 150 and fibre for 200 to not add variables
		  if(t.beamEnergy==200 && t.cef3_maxAmpl->at(2)<theConfiguration_.channel2CutFibre)continue;
		  if(t.beamEnergy==energies[j]) 
		    reso_histo_channel_total[j]->Fill(deltaTNoCorr);
		  continue;
	      }
	    }
	  }
	}
      }
      
      
    }
    
    
  }

   std::string constDirName = "plots_timingPerformanceKur_";
   constDirName+=theConfiguration_.setup;
   if(theConfiguration_.addTagFileName){
     constDirName+="_";
     constDirName+=theConfiguration_.tagFileName;
   }
   system(Form("mkdir -p %s", constDirName.c_str()));
   TString dir(constDirName);

   for (int i=0;i<theConfiguration_.nMaxAmplCuts;++i){
     std::cout<<"###### entries:"<<i<<" "<<reso_histo_channel_corr_Amplitude[i]->GetEntries()<<std::endl;
     //channel
     if(reso_histo_channel_corr_Amplitude[i]->GetEntries()>25 ){
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

       RooRealVar meanr("meanr","Mean",peakpos,peakpos-3*sigma, peakpos+3*sigma);
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
       //       rms = (widthL.getVal()+widthR.getVal())/2; FIXME
       //       rmsErr = 0.5*sqrt(widthL.getError()*widthL.getError()+widthR.getError()*widthR.getError());

       rms = (widthL.getVal()+widthR.getVal())/2; 
       rmsErr = 0.5*sqrt(widthL.getError()*widthL.getError()+widthR.getError()*widthR.getError());


       resValueAmplitude_channel[i]=rms*1.e3; 
       resErrValueAmplitude_channel[i]=rmsErr*1.e3; 


       TCanvas* cans = new TCanvas();
       cans->cd();
       frame->Draw();
       TLegend* lego = new TLegend(0.57, 0.8, 0.89, 0.9);;
       lego->SetTextSize(0.036);
       lego->SetTextAlign(32); // align right
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


   //total with reasonable cut on ampl

     //channel
   for(int j=0;j<energies.size();j++){
   if(reso_histo_channel_total[j]->GetEntries()>25 ){
     TH1F* histo;
     histo=reso_histo_channel_total[j];
     double peakpos = histo->GetBinCenter(histo->GetMaximumBin());
     double sigma = histo->GetRMS();
     
     double fitmin;
     double fitmax;
     
  
     fitmin = peakpos-4*sigma;
     fitmax = peakpos+4*sigma;
     
     RooRealVar x("x","deltaT", fitmin, fitmax);
     RooDataHist data("data","dataset with x",x,RooFit::Import(*histo) );
     
     RooRealVar meanr("meanr","Mean",peakpos,peakpos-3*sigma, peakpos+3*sigma);
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
     //       rms = (widthL.getVal()+widthR.getVal())/2; FIXME
     //       rmsErr = 0.5*sqrt(widthL.getError()*widthL.getError()+widthR.getError()*widthR.getError());
     
     rms = (widthL.getVal()+widthR.getVal())/2; 
     rmsErr = 0.5*sqrt(widthL.getError()*widthL.getError()+widthR.getError()*widthR.getError());
     
     
     resValueEnergy_channel[j]=rms*1.e3; 
     resErrValueEnergy_channel[j]=rmsErr*1.e3; 
     
     
     TCanvas* cans = new TCanvas();
     cans->cd();
     frame->Draw();
     TLegend* lego = new TLegend(0.57, 0.8, 0.89, 0.9);;
     lego->SetTextSize(0.036);
     lego->SetTextAlign(32); // align right
     lego->AddEntry(  (TObject*)0 ,Form("#sigma = %.1f #pm %.1f ps", rms*1.e3, rmsErr*1.e3), "");
     
     lego->SetFillColor(0);
     lego->Draw("same");
     
     TPaveText* pave = DrawTools::getLabelTop_expOnXaxis("Electron Beam");
     pave->Draw("same");
     
    TString icut;
    icut.Form("%d",(int)energies[j]); 

     cans->SaveAs(dir+"/reso_histo_channel_cruijff_total"+icut+".png");
     cans->SaveAs(dir+"/reso_histo_channel_cruijff_total"+icut+".pdf");
   }
   }

   TCanvas* cans = new TCanvas();
   cans->cd();
   //   reso_histo_channel_total->GetXaxis()->SetTitle("time_{fibre}-time_{mcp} [ns]");
   //   reso_histo_channel_total->Draw();
   //   cans->SaveAs(dir+"/reso_histo_channel_total.png");
   //   cans->SaveAs(dir+"/reso_histo_channel_total.pdf");
  
   cans->Clear();
   cans->cd();
   maxAmplHisto->GetXaxis()->SetTitle("Amplitude [ADC]");
   maxAmplHisto->Draw();
   cans->SaveAs(dir+"/maxAmplHisto.png");
   cans->SaveAs(dir+"/maxAmplHisto.pdf");


   //res vs ampl
   TGraphErrors*   resVsAmplitude_channel = new TGraphErrors(0);
   TGraphErrors*   resVsEnergy_channel = new TGraphErrors(0);

   for(int i=0;i<resValueAmplitude_channel.GetNoElements();i++){
     if(resValueAmplitude_channel[i]!=0){ 
       resVsAmplitude_channel->SetPoint(i,theConfiguration_.startCutChannel+(i+1)*theConfiguration_.stepAmplChannel-theConfiguration_.stepAmplChannel/2.,resValueAmplitude_channel[i]); 
       resVsAmplitude_channel->SetPointError(i,0,resErrValueAmplitude_channel[i]);
     }
   }


   for(int i=0;i<resValueEnergy_channel.GetNoElements();i++){
     if(resValueEnergy_channel[i]!=0){ 
       resVsEnergy_channel->SetPoint(i,energies[i],resValueEnergy_channel[i]); 
       resVsEnergy_channel->SetPointError(i,0,resErrValueEnergy_channel[i]);
     }
   }


   cans->Clear();
   float yup=1.1*(resVsAmplitude_channel->GetY()[0]+resVsAmplitude_channel->GetEY()[0]);
   float xup=(resVsAmplitude_channel->GetX())[resVsAmplitude_channel->GetN()-1]+10; 
   TH2D* h2_axes = new TH2D( "axes", "", 100,resVsAmplitude_channel->GetX()[0]-10 ,xup , 110, 50, 500);
   
   h2_axes->SetYTitle("#sigma_{t} [ps]");
   h2_axes->GetXaxis()->SetTitle("Amplitude [ADC]");
   h2_axes->Draw(""); 
   resVsAmplitude_channel->SetName("resVsAmplitude_channel");
   resVsAmplitude_channel->SetMarkerStyle(20);
   resVsAmplitude_channel->SetMarkerSize(1.6);
   resVsAmplitude_channel->SetMarkerColor(kBlue);
   

   TF1* f= new TF1("fun","sqrt([1]*[1]/(x*x)+[0]*[0])",(resVsAmplitude_channel->GetX())[0]-2,(resVsAmplitude_channel->GetX())[resVsAmplitude_channel->GetN()-1] +10);
   
   
   resVsAmplitude_channel->Fit("fun","R");
   
   resVsAmplitude_channel->Draw("p same");
   
   
   TPaveText* pave;
   if(runs.size()>1) pave= DrawTools::getLabelTop_expOnXaxis("Electron Beam");
   else {
     std::string energy(Form("%.0f", energies[0]));
     pave= DrawTools::getLabelTop_expOnXaxis(energy+" GeV Electron Beam");
   }

   pave->Draw("same");


   cans->Clear();
   float yup_2=1.1*(resVsEnergy_channel->GetY()[0]+resVsEnergy_channel->GetEY()[0]);
   float xup_2=(resVsEnergy_channel->GetX())[resVsEnergy_channel->GetN()-1]+10; 
   TH2D* h2_axes_2 = new TH2D( "axes_2", "", 100,resVsEnergy_channel->GetX()[0]-10 ,xup_2 , 110, 50, yup_2);
   
   h2_axes_2->SetYTitle("#sigma_{t} [ps]");
   h2_axes_2->GetXaxis()->SetTitle("Energy [GeV]");
   h2_axes_2->Draw(""); 
   resVsEnergy_channel->SetName("resVsEnergy_channel");
   resVsEnergy_channel->SetMarkerStyle(20);
   resVsEnergy_channel->SetMarkerSize(1.6);
   resVsEnergy_channel->SetMarkerColor(kBlue);
   

   TF1* f2= new TF1("fun2","sqrt([1]*[1]/(x*x)+[0]*[0])",(resVsEnergy_channel->GetX())[0]-2,(resVsEnergy_channel->GetX())[resVsEnergy_channel->GetN()-1] +10);
   
   
   resVsEnergy_channel->Fit("fun2","R");
   
   resVsEnergy_channel->Draw("p same");
   
   
   if(runs.size()>1) pave= DrawTools::getLabelTop_expOnXaxis("Electron Beam");
   else {
     std::string energy(Form("%.0f", energies[0]));
     pave= DrawTools::getLabelTop_expOnXaxis(energy+" GeV Electron Beam");
   }

   pave->Draw("same");

   
   
   
   TLegend* lego = new TLegend(0.47, 0.8, 0.8, 0.92);
   lego->SetTextSize(0.038);
   //    lego->AddEntry(  (TObject*)0 ,"f(x) = p0 + p1/x", "");
   lego->AddEntry(  (TObject*)0 ,Form("C = %.0f #pm %.0f ps", f2->GetParameter(0), f2->GetParError(0) ), "");
   //    lego->AddEntry(  (TObject*)0 ,Form("N = %.0f #pm %.0f", f->GetParameter(1), f->GetParError(1) ), "");
   lego->SetFillColor(0);
   lego->Draw("same");
   
   cans->SaveAs(dir+"/timingResolutionVsEnergyMAPD_channel.png");
   cans->SaveAs(dir+"/timingResolutionVsEnergyMAPD_channel.pdf");  
   

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

   conf.startCutFibre= Configurator::GetDouble(Configurable::getElementContent(*configurator_,"startCutFibre",configurator_->root_element));

   conf.startCutChannel= Configurator::GetDouble(Configurable::getElementContent(*configurator_,"startCutChannel",configurator_->root_element));

   conf.amplCut= Configurator::GetDouble(Configurable::getElementContent(*configurator_,"amplCut",configurator_->root_element));

   conf.addTagFileName= Configurator::GetInt(Configurable::getElementContent(*configurator_,"addTagFileName",configurator_->root_element));

   conf.tagFileName=Configurable::getElementContent (*configurator_, "tagFileName",configurator_->root_element) ;

  content = Configurable::getElementContent (*configurator_, "energies",configurator_->root_element) ;
  Configurator::GetVecFloat(content,conf.energies);


   return conf;

}


