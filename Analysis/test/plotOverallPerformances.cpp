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
    std::cout << "./plotOverallPerformances ([tag]) ([config])" << std::endl;
    exit(12345);
  }

  std::cout<<config<<std::endl;
  
  theConfiguration_=readConfiguration(config);
  
  TH1F* reso_histo_fibre_corr_Amplitude[theConfiguration_.nMaxAmplCuts];//to see behaviour of resolution as a function of ampl
  TH1F* reso_histo_channel_corr_Amplitude[theConfiguration_.nMaxAmplCuts];

  TH1F* reso_histo_channelPlusFibre_corr_Amplitude[theConfiguration_.nMaxAmplCuts];

  TH1F* reso_histo_fibre_total;//all events at different energies all together with reasonable cut on amplitude
  TH1F* reso_histo_channel_total;

  TVectorD resValueAmplitude_fibre(theConfiguration_.nMaxAmplCuts);
  TVectorD resErrValueAmplitude_fibre(theConfiguration_.nMaxAmplCuts);
  
  TVectorD resValueAmplitude_channel(theConfiguration_.nMaxAmplCuts);
  TVectorD resErrValueAmplitude_channel(theConfiguration_.nMaxAmplCuts);
  TVectorD resValueAmplitude_channelPlusFibre(theConfiguration_.nMaxAmplCuts);
  TVectorD resErrValueAmplitude_channelPlusFibre(theConfiguration_.nMaxAmplCuts);
  //  TVectorD resValueAmplitudeSingleEnergy_channel(theConfiguration_.runs.size());
  //  TVectorD resValueAmplitudeSingleEnergy_fibre(theConfiguration_.runs.size());

  float shift=10;
  bool isNino=false;

  if(theConfiguration_.setup=="Nino"){
    shift=11.5;
    isNino=true;
  }


   for (int i=0;i<theConfiguration_.nMaxAmplCuts;++i){
    TString icut;
    icut.Form("%d",i); 
    reso_histo_fibre_corr_Amplitude[i] = new TH1F("reso_histo_fibre_corr_Amplitude_"+icut,"reso_histo_fibre_corr_Amplitude_"+icut,200,theConfiguration_.rangeXLow,theConfiguration_.rangeXUp);
    reso_histo_channel_corr_Amplitude[i] = new TH1F("reso_histo_channel_corr_Amplitude_"+icut,"reso_histo_channel_corr_Amplitude_"+icut,200,theConfiguration_.rangeXLow,theConfiguration_.rangeXUp);
    reso_histo_channelPlusFibre_corr_Amplitude[i] = new TH1F("reso_histo_channelPlusFibre_corr_Amplitude_"+icut,"reso_histo_channelPlusFibre_corr_Amplitude_"+icut,200,theConfiguration_.rangeXLow,theConfiguration_.rangeXUp);
   }
   
   reso_histo_fibre_total = new TH1F("reso_histo_fibre_total","reso_histo_fibre_total",200,theConfiguration_.rangeXLow+shift,theConfiguration_.rangeXUp+shift);
   reso_histo_channel_total = new TH1F("reso_histo_channel_total","reso_histo_channel_total",200,theConfiguration_.rangeXLow+shift,theConfiguration_.rangeXUp+shift);
   
   
   std::vector<int> runs;
   std::vector<int> energies;

  TString haddFile = "analysisTrees_"+tag+"/"+theConfiguration_.setup+"_runs.root";

  if(theConfiguration_.addTagFileName){
    haddFile = "analysisTrees_"+tag+"/"+theConfiguration_.setup+"_"+theConfiguration_.tagFileName+"_runs.root";
  }

  TString haddString = "hadd analysisTrees_"+tag+"/"+theConfiguration_.setup+"_runs.root ";

  if(theConfiguration_.addTagFileName){
    haddString = "hadd analysisTrees_"+tag+"/"+theConfiguration_.setup+"_"+theConfiguration_.tagFileName+"_runs.root ";
  }

  for (int i=0;i<theConfiguration_.runs.size();++i){
    TString fileNameString = "analysisTrees_"+tag+"/Reco_";
    runs.push_back(theConfiguration_.runs[i]);
    energies.push_back(theConfiguration_.energies[i]);
    fileNameString+=runs[i];
    haddString+=fileNameString;
    haddString+= ".root ";
  }

  system(haddString.Data());

  TFile* file = TFile::Open(haddFile.Data());

  TTree* recoTree=(TTree*)file->Get("recoTree");
  RecoTree t(recoTree);
  Long64_t nentries = t.fChain->GetEntries();



  int dummyCounter=0;
  float lowerBeamEnergy=0.;//cuts are sliding in energy
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = t.LoadTree(jentry);
      if (ientry < 0) break;
      nb = t.fChain->GetEntry(jentry);   nbytes += nb;
      if(jentry==0)lowerBeamEnergy=t.beamEnergy;
      if( t.mcp_time_frac50<0 ||t.mcp_max_amplitude<200)continue;
      for (int i=0;i<theConfiguration_.nMaxAmplCuts;++i){
	bool filled_channel=false;
	bool filled_channelPlusFibre=false;
	if(TopologicalSelectionHelper::passesFibreTopologicalSelection(t,isNino) &&  t.cef3_maxAmpl->at(2) < theConfiguration_.channel2CutFibre*t.beamEnergy/lowerBeamEnergy){ 
	  if(t.cef3_maxAmpl->at(1)>theConfiguration_.startCutFibre+i*theConfiguration_.stepAmplFibre && t.cef3_maxAmpl->at(1)<theConfiguration_.startCutFibre+(i+1)*theConfiguration_.stepAmplFibre){
	    reso_histo_fibre_corr_Amplitude[i]->Fill(t.cef3_time_at_frac50->at(1)-t.mcp_time_frac50);
	    if(t.cef3_maxAmpl->at(1)>theConfiguration_.amplCut)reso_histo_fibre_total->Fill(t.cef3_time_at_frac50->at(1)-t.mcp_time_frac50+shift);
	    reso_histo_channelPlusFibre_corr_Amplitude[i]->Fill(t.cef3_time_at_frac50->at(1)-t.mcp_time_frac50);
	    break;
	  }
      }
	if(TopologicalSelectionHelper::passesChannelTopologicalSelection(t,isNino)  && t.cef3_maxAmpl->at(2) > theConfiguration_.channel2CutChannel*t.beamEnergy/lowerBeamEnergy && t.cef3_maxAmpl->at(1)<theConfiguration_.channel1CutChannel){   //cuts on fibre position with hodos and wc. cut on cef3_maxAmpl[2] to reduce hadron contamination, cef3_maxAmpl[1] cut to reduce remaining events hitting the fibre
	    if(t.cef3_maxAmpl->at(1)>i*theConfiguration_.stepAmplFibre && t.cef3_maxAmpl->at(1)<(i+1)*theConfiguration_.stepAmplFibre){
	      reso_histo_channelPlusFibre_corr_Amplitude[i]->Fill(t.cef3_time_at_frac50->at(1)-t.mcp_time_frac50);
	      filled_channelPlusFibre=true;
	      if(filled_channel)break;
	    }

	    if(t.cef3_maxAmpl->at(1)>theConfiguration_.startCutChannel+i*theConfiguration_.stepAmplChannel && t.cef3_maxAmpl->at(1)<theConfiguration_.startCutChannel+(i+1)*theConfiguration_.stepAmplChannel){
	    reso_histo_channel_corr_Amplitude[i]->Fill(t.cef3_time_at_frac50->at(1)-t.mcp_time_frac50);
	    if(t.cef3_maxAmpl->at(1)>theConfiguration_.amplCut)reso_histo_channel_total->Fill(t.cef3_time_at_frac50->at(1)-t.mcp_time_frac50+shift);
	    filled_channel=true;
	    if(filled_channelPlusFibre) break;
	    }
	  }
      }
      
   }


   // this is the dir in which the plots will be saved:
   std::string constDirName = "plots_timingPerformance_";
   constDirName+=theConfiguration_.setup;
   if(theConfiguration_.addTagFileName){
     constDirName+="_";
     constDirName+=theConfiguration_.tagFileName;
   }
   system(Form("mkdir -p %s", constDirName.c_str()));
   TString dir(constDirName);
   
   //fit histos
   for (int i=0;i<theConfiguration_.nMaxAmplCuts;++i){
     //fibre
     std::cout<<"#####################################"<<i<<" "<<reso_histo_fibre_corr_Amplitude[i]->GetEntries()<<std::endl;
     if(reso_histo_fibre_corr_Amplitude[i]->GetEntries()>=25 && !isNino){
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
     rms = (widthL.getVal()+widthR.getVal())/2;
     rmsErr = 0.5*sqrt(widthL.getError()*widthL.getError()+widthR.getError()*widthR.getError());

     resValueAmplitude_fibre[i]=rms*1.e3; 
     resErrValueAmplitude_fibre[i]=rmsErr*1.e3; 


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

     cans->SaveAs(dir+"/reso_histo_fibre_cruijff_"+icut+".png");
     cans->SaveAs(dir+"/reso_histo_fibre_cruijff_"+icut+".pdf");

     }

     //channel
     if(reso_histo_channel_corr_Amplitude[i]->GetEntries()>25  && !isNino){
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


     //channelPlusFibre
     if(reso_histo_channelPlusFibre_corr_Amplitude[i]->GetEntries()>25  && !isNino){
       TH1F* histo;
       histo=reso_histo_channelPlusFibre_corr_Amplitude[i];
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
       rms = (widthL.getVal()+widthR.getVal())/2;
       rmsErr = 0.5*sqrt(widthL.getError()*widthL.getError()+widthR.getError()*widthR.getError());

       resValueAmplitude_channelPlusFibre[i]=rms*1.e3; 
       resErrValueAmplitude_channelPlusFibre[i]=rmsErr*1.e3; 


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

       cans->SaveAs(dir+"/reso_histo_channelPlusFibre_cruijff_"+icut+".png");
       cans->SaveAs(dir+"/reso_histo_channelPlusFibre_cruijff_"+icut+".pdf");
     }



   }
	  

   TGraphErrors*   resVsAmplitude_fibre = new TGraphErrors(0);
   TGraphErrors*   resVsAmplitude_channel = new TGraphErrors(0);
   TGraphErrors*   resVsAmplitude_channelPlusFibre = new TGraphErrors(0);


   if(!isNino){
     //plot res vs amplitude

     resVsAmplitude_channel->SetName("resVsAmplitude_channel");
     resVsAmplitude_fibre->SetName("resVsAmplitude_fibre");
     resVsAmplitude_channelPlusFibre->SetName("resVsAmplitude_channelPlusFibre");
   
     for(int i=0;i<resValueAmplitude_channel.GetNoElements();i++){
       if(resValueAmplitude_channel[i]!=0){ 
	 //     std::cout<<"setting value"<<i<<" "<<(i+1)*theConfiguration_.stepAmplChannel-theConfiguration_.stepAmplChannel/2.<<" "<<resValueAmplitude_channel[i]<<std::endl;
	 resVsAmplitude_channel->SetPoint(i,(i+1)*theConfiguration_.stepAmplChannel-theConfiguration_.stepAmplChannel/2.,resValueAmplitude_channel[i]); 
	 resVsAmplitude_channel->SetPointError(i,0,resErrValueAmplitude_channel[i]);
       }
       if(resValueAmplitude_fibre[i]!=0){
	 //     std::cout<<"setting value"<<i<<" "<<(i+1)*theConfiguration_.stepAmplFibre/2.<<" "<<resValueAmplitude_fibre[i]<<std::endl;
	 resVsAmplitude_fibre->SetPoint(i,theConfiguration_.startCutFibre+(i+1)*theConfiguration_.stepAmplFibre-theConfiguration_.stepAmplFibre/2.,resValueAmplitude_fibre[i]);
	 resVsAmplitude_fibre->SetPointError(i,0,resErrValueAmplitude_fibre[i]);
       }
       if(resValueAmplitude_channelPlusFibre[i]!=0){ 
	 //     std::cout<<"setting value"<<i<<" "<<(i+1)*theConfiguration_.stepAmplChannel-theConfiguration_.stepAmplChannel/2.<<" "<<resValueAmplitude_channel[i]<<std::endl;
	 resVsAmplitude_channelPlusFibre->SetPoint(i,(i+1)*theConfiguration_.stepAmplFibre-theConfiguration_.stepAmplFibre/2.,resValueAmplitude_channelPlusFibre[i]); 
	 resVsAmplitude_channelPlusFibre->SetPointError(i,0,resErrValueAmplitude_channelPlusFibre[i]);
       }

     }
   }
	 
	 
   std::string outFileName;
   outFileName = dir+"/plotOverallPerformances_"+tag+".root";
   TFile* outFile = TFile::Open(outFileName.c_str(),"recreate");

   if(!isNino){
     TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
     float yup=1.1*(resVsAmplitude_channel->GetY()[2]+resVsAmplitude_channel->GetEY()[0]);
     float xup=(resVsAmplitude_channel->GetX())[resVsAmplitude_channel->GetN()-5]+10; 
     TH2D* h2_axes = new TH2D( "axes", "", 100,resVsAmplitude_channel->GetX()[2]-10 ,xup , 110, 0., yup);


     //channel   
     resVsAmplitude_channel->SetName("resVsAmplitude_channel");
     resVsAmplitude_channel->SetMarkerStyle(20);
     resVsAmplitude_channel->SetMarkerSize(1.6);
     resVsAmplitude_channel->SetMarkerColor(kBlue);


     resVsAmplitude_channel->Write();
     h2_axes->SetYTitle("#sigma_{t} [ps]");
     h2_axes->GetXaxis()->SetTitle("Amplitude [ADC]");
     h2_axes->Draw(""); 

     TF1* f= new TF1("fun","sqrt([1]*[1]/(x*x)+[0]*[0])",(resVsAmplitude_channel->GetX())[2]-2,(resVsAmplitude_channel->GetX())[resVsAmplitude_channel->GetN()-5] +10);
     //      TF1* f= new TF1("fun","sqrt([1]*[1]/(x*x)+[0]*[0])",(resVsAmplitude_channel->GetX())[1]-2,(resVsAmplitude_channel->GetX())[resVsAmplitude_channel->GetN()-7] +10);
     f->SetParLimits(0,50,180);
     f->SetParameter(0,100);
     f->SetParameter(1,1.12135e+04);
     //   f->SetParameter(1,300);
     //     f->SetParameter(2,1.41901e+03);
     //   f->SetParameter(2,100);
     //   f->SetParLimits(2,100,1000);
     std::cout<<" ROOT FIT::::"<<f->Eval(60.)<<std::endl;

     //     f->SetParLimits(1,10,1.19257e+02);
     resVsAmplitude_channel->Fit("fun","R");
   
     resVsAmplitude_channel->Draw("p same");


     TPaveText* pave = DrawTools::getLabelTop_expOnXaxis("Electron Beam");
     pave->Draw("same");



     TLegend* lego = new TLegend(0.47, 0.7, 0.8, 0.92);
     lego->SetTextSize(0.038);
     //    lego->AddEntry(  (TObject*)0 ,"f(x) = p0 + p1/x", "");
     lego->AddEntry(  (TObject*)0 ,Form("C = %.0f #pm %.0f ps", f->GetParameter(0), f->GetParError(0) ), "");
     //    lego->AddEntry(  (TObject*)0 ,Form("N = %.0f #pm %.0f", f->GetParameter(1), f->GetParError(1) ), "");
     lego->SetFillColor(0);
     lego->Draw("same");

     c1->SaveAs(dir+"/timingResolutionVsAmplitudeSiPM_channel.png");
     c1->SaveAs(dir+"/timingResolutionVsAmplitudeSiPM_channel.pdf");  

     std::cout<<resVsAmplitude_channel->GetN()<<" "<<(resVsAmplitude_channel->GetX())[resVsAmplitude_channel->GetN()-1]<< " "<< (resVsAmplitude_channel->GetX())[resVsAmplitude_channel->GetN()]<<std::endl;

     //fibre
     c1->Clear();
     yup=1.1*(resVsAmplitude_fibre->GetY()[0]+resVsAmplitude_fibre->GetEY()[1]);
     xup=(resVsAmplitude_fibre->GetX())[resVsAmplitude_fibre->GetN()-1]+10; 
     TH2D* h2_axes_2 = new TH2D( "axes_2", "", 100, resVsAmplitude_channel->GetX()[1],xup , 110, 0., yup);


     //fibre   
     resVsAmplitude_fibre->SetName("resVsAmplitude_fibre");
     resVsAmplitude_fibre->SetMarkerStyle(20);
     resVsAmplitude_fibre->SetMarkerSize(1.6);
     resVsAmplitude_fibre->SetMarkerColor(kBlue);
     resVsAmplitude_fibre->Write();
     h2_axes_2->SetYTitle("#sigma_{t} [ps]");
     h2_axes_2->GetXaxis()->SetTitle("Amplitude [ADC]");
     h2_axes_2->Draw(""); 
   
     //   TF1* f2= new TF1("fun2","[1]/x+[0]",(resVsAmplitude_fibre->GetX())[1]-10,(resVsAmplitude_fibre->GetX())[resVsAmplitude_fibre->GetN()-1] +10);
     TF1* f2= new TF1("fun2","sqrt([1]*[1]/(x*x)+[0]*[0])",(resVsAmplitude_fibre->GetX())[1]-10,(resVsAmplitude_fibre->GetX())[resVsAmplitude_fibre->GetN()-1] +10);
     //    f2->SetParLimits(0,50,150);
     //  f2->FixParameter(1,1.66061e+04);
     //  f2->FixParameter(2,8.02308e+02);
     resVsAmplitude_fibre->Fit("fun2","R");
  
     resVsAmplitude_fibre->Draw("p same");


     pave->Draw("same");

   
   
     TLegend* lego_2 = new TLegend(0.47, 0.7, 0.8, 0.92);
     lego_2->SetTextSize(0.038);
     //   lego_2->AddEntry(  (TObject*)0 ,"f(x) = p0 + p1/x", "");
     lego_2->AddEntry(  (TObject*)0 ,Form("C = %.0f #pm %.0f ps", f2->GetParameter(0), f2->GetParError(0) ), "");
     //   lego_2->AddEntry(  (TObject*)0 ,Form("N = %.0f #pm %.0f", f2->GetParameter(1), f2->GetParError(1) ), "");
     lego_2->SetFillColor(0);
     lego_2->Draw("same");
   
     c1->SaveAs(dir+"/timingResolutionVsAmplitudeSiPM_fibre.png");
     c1->SaveAs(dir+"/timingResolutionVsAmplitudeSiPM_fibre.pdf");  



     c1->Clear();
     h2_axes_2->Draw("");
     resVsAmplitude_channelPlusFibre->SetName("resVsAmplitude_channelPlusFibre");
     resVsAmplitude_channelPlusFibre->SetMarkerStyle(20);
     resVsAmplitude_channelPlusFibre->SetMarkerSize(1.6);
     resVsAmplitude_channelPlusFibre->SetMarkerColor(kBlue);
     resVsAmplitude_channelPlusFibre->Write();

     resVsAmplitude_channelPlusFibre->Fit("fun2","R");

     resVsAmplitude_channelPlusFibre->Draw("p same");
   

     pave->Draw("same");
   
     TLegend* lego_3 = new TLegend(0.47, 0.7, 0.8, 0.92);
     lego_3->SetTextSize(0.038);
     //   lego_3->AddEntry(  (TObject*)0 ,"f(x) = p0 + p1/x", "");
     lego_3->AddEntry(  (TObject*)0 ,Form("C = %.0f #pm %.0f ps", f2->GetParameter(0), f2->GetParError(0) ), "");
     //   lego_3->AddEntry(  (TObject*)0 ,Form("N = %.0f #pm %.0f", f2->GetParameter(1), f2->GetParError(1) ), "");
     lego_3->SetFillColor(0);
     lego_3->Draw("same");
   
     c1->SaveAs(dir+"/timingResolutionVsAmplitudeSiPM_channelPlusFibre.png");
     c1->SaveAs(dir+"/timingResolutionVsAmplitudeSiPM_channelPlusFibre.pdf");  



     c1->Clear();
     h2_axes_2->Draw(""); 
     resVsAmplitude_fibre->GetFunction("fun2")->SetBit(TF1::kNotDraw);
     resVsAmplitude_fibre->Draw("p same");
     resVsAmplitude_channel->SetLineColor(kRed);
     resVsAmplitude_channel->GetFunction("fun")->SetBit(TF1::kNotDraw);
     resVsAmplitude_channel->SetMarkerColor(kRed);
     resVsAmplitude_channel->Draw("p same");
     c1->SaveAs(dir+"/timingResolutionVsAmplitudeSiPM_fibreAndChannel.png");
     c1->SaveAs(dir+"/timingResolutionVsAmplitudeSiPM_fibreAndChannel.pdf");  
   }

   //fit histos for all energies with reasonable cut on ampl
   //fibre
    std::cout<<reso_histo_fibre_total->GetEntries();
   if(reso_histo_fibre_total->GetEntries()>=25){
     TH1F* histo;
     histo=reso_histo_fibre_total;
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
     rms = (widthL.getVal()+widthR.getVal())/2;
     rmsErr = 0.5*sqrt(widthL.getError()*widthL.getError()+widthR.getError()*widthR.getError());

     TCanvas* cans = new TCanvas();
     cans->cd();
     frame->Draw();
     TLegend* lego = new TLegend(0.57, 0.8, 0.89, 0.9);
     lego->SetTextAlign(32); // align right
     lego->SetTextSize(0.036);
     lego->AddEntry(  (TObject*)0 ,Form("#sigma = %.1f #pm %.1f ps", rms*1.e3, rmsErr*1.e3), "");
     
     lego->SetFillColor(0);
     lego->Draw("same");
     
     TPaveText* pave = DrawTools::getLabelTop_expOnXaxis("Electron Beam");
     pave->Draw("same");
     

     cans->SaveAs(dir+"/reso_histo_fibre_total.png");
     cans->SaveAs(dir+"/reso_histo_fibre_total.pdf");

     }

     //channel
     if(reso_histo_channel_total->GetEntries()>25){
       TH1F* histo;
       histo=reso_histo_channel_total;
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
       rms = (widthL.getVal()+widthR.getVal())/2;
       rmsErr = 0.5*sqrt(widthL.getError()*widthL.getError()+widthR.getError()*widthR.getError());


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
     
       cans->SaveAs(dir+"/reso_histo_channel_total.png");
       cans->SaveAs(dir+"/reso_histo_channel_total.pdf");
     }



   
   outFile->Write();
   outFile->Close();

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


