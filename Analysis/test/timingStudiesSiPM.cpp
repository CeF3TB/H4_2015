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

#include "interface/Configurator.h"
#include "timingPlotsConfigurator.h"
#include "interface/TopologicalSelectionHelper.h"

timingPlots_Config_t readConfiguration(std::string configName, float beamEnergy);


int main( int argc, char* argv[] ) {

  DrawTools::setStyle();



   if( argc<2 ) {
     std::cout << "ERROR. You need to specify the name of the run you want to process." << std::endl;  
     exit(1);
   }

   std::string runName = "";
   std::string tag = "V00";
   std::string config = "SiPM2015Config";

   if( argc>1 ) {
     std::string runName_str(argv[1]);
     runName = runName_str;
     if( argc>2 ) {
       std::string tag_str(argv[2]);
       tag = tag_str;
       if(argc>3){
	 std::string config_str(argv[3]);
	 config=config_str;
       }
     }
   } else {
     std::cout << "Usage:" << std::endl;
     std::cout << "./timingStudies [runName] ([tag]) ([config])" << std::endl;
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

  TTree* recoTree=(TTree*)file->Get("recoTree");
  RecoTree t(recoTree);
  Long64_t nentries = t.fChain->GetEntries();

  t.fChain->GetEntry(0);
  theConfiguration_=readConfiguration(config,t.beamEnergy);

  float xLowChannel=-5.5, xUpChannel=4, yLowChannel=-2, yUpChannel=4.5;
  float xLowFibre=-5.5,xUpFibre=-2.5, yLowFibre=-3, yUpFibre=0;

  if(theConfiguration_.setup=="Nino"){
    xLowChannel=-4; xUpChannel=4; yLowChannel=0.5, yUpChannel=4.5;
    xLowFibre=-5; xUpFibre=0.3; yLowFibre=0.5, yUpFibre=3.9;
  }


  TH1F* reso_histo_fibre = new TH1F("reso_histo_fibre","reso_histo_fibre",2000,theConfiguration_.rangeXLow,theConfiguration_.rangeXUp);  
  TH1F* reso_histo_channel = new TH1F("reso_histo_channel","reso_histo_channel",2000,theConfiguration_.rangeXLow,theConfiguration_.rangeXUp);  
  //maps
  TH2F* timing_map_fibre = new TH2F("timing_map_fibre","timing_map_fibre",16,xLowFibre,xUpFibre,16,yLowFibre,yUpFibre);
  TH2F* amplitude_map_fibre = new TH2F("amplitude_map_fibre","amplitude_map_fibre",16,xLowFibre,xUpFibre,16,yLowFibre,yUpFibre);
  TH2F* timing_map_fibre_norm = new TH2F("timing_map_fibre_norm","timing_map_fibre_norm",16,xLowFibre,xUpFibre,16,yLowFibre,yUpFibre);

  TH2F* timing_map_channel = new TH2F("timing_map_channel","timing_map_channel",22,xLowChannel,xUpChannel,22,-5.5,4.5);
  TH2F* amplitude_map_channel = new TH2F("amplitude_map_channel","amplitude_map_channel",22,xLowChannel,xUpChannel,22,-5.5,4.5);
  TH2F* timing_map_channel_norm = new TH2F("timing_map_channel_norm","timing_map_channel_norm",22,xLowChannel,xUpChannel,22,-5.5,4.5);

  TH2F* amplitude_map_fibre2 = new TH2F("amplitude_map_fibre2","amplitude_map_fibre2",22,xLowChannel,xUpChannel,22,-5.5,4.5);

  TH2F* amplitude_map_fibre0 = new TH2F("amplitude_map_fibre0","amplitude_map_fibre0",22,xLowChannel,xUpChannel,22,-5.5,4.5);
  TH2F* amplitude_map_fibre0and2 = new TH2F("amplitude_map_fibre0and2","amplitude_map_fibre0and2",22,xLowChannel,xUpChannel,22,-5.5,4.5);
  
  //plots with all the selection applied
  TH1F* maxAmpl_sel_fibre = new TH1F("maxAmpl_sel_fibre","maxAmpl_sel_fibre",500,0,4000);
  TH1F* maxAmpl_sel_channel = new TH1F("maxAmpl_sel_channel","maxAmpl_sel_channel",500,0,4000);

  int nBinsChannel=30, nBinsFibre=12;

  TH2F* amplitude_map_sel_channel = new TH2F("amplitude_map_sel_channel","amplitude_map_sel_channel",nBinsChannel,xLowChannel,xUpChannel,nBinsChannel,yLowChannel,yUpChannel);
  TH2F* amplitude_map_sel_channel_ampl_cut = new TH2F("amplitude_map_sel_channel_ampl_cut","amplitude_map_sel_channel_ampl_cut",nBinsChannel,xLowChannel,xUpChannel,nBinsChannel,yLowChannel,yUpChannel);
  TH2F* amplitude_map_fibre2_sel_channel_ampl_cut = new TH2F("amplitude_map_fibre2_sel_channel_ampl_cut","amplitude_map_fibre2_sel_channel_ampl_cut",nBinsChannel,xLowChannel,xUpChannel,nBinsChannel,yLowChannel,yUpChannel);
  TH2F* amplitude_map_sel_fibre = new TH2F("amplitude_map_sel_fibre","amplitude_map_sel_fibre",nBinsFibre,xLowFibre,xUpFibre,nBinsFibre,yLowFibre,yUpFibre);

  TH2F* timing_map_sel_channel = new TH2F("timing_map_sel_channel","timing_map_sel_channel",nBinsChannel,xLowChannel,xUpChannel,nBinsChannel,yLowChannel,yUpChannel);
  TH2F* timing_map_sel_channel_ampl_cut = new TH2F("timing_map_sel_channel_ampl_cut","timing_map_sel_channel_ampl_cut",nBinsChannel,xLowChannel,xUpChannel,nBinsChannel,yLowChannel,yUpChannel);
  TH2F* timing_map_sel_fibre = new TH2F("timing_map_sel_fibre","timing_map_sel_fibre",nBinsFibre,xLowFibre,xUpFibre,nBinsFibre,yLowFibre,yUpFibre);

  TH2F* sel_channel_norm = new TH2F("sel_channel_norm","sel_channel_norm",nBinsChannel,xLowChannel,xUpChannel,nBinsChannel,yLowChannel,yUpChannel);
  TH2F* sel_channel_ampl_cut_norm = new TH2F("sel_channel_ampl_cut_norm","sel_channel_ampl_cut_norm",nBinsChannel,xLowChannel,xUpChannel,nBinsChannel,yLowChannel,yUpChannel);
  TH2F* sel_channel_ampl_and_time_cut_norm = new TH2F("sel_channel_ampl_and_time_cut_norm","sel_channel_ampl_and_time_cut_norm",nBinsChannel,xLowChannel,xUpChannel,nBinsChannel,yLowChannel,yUpChannel);
  TH2F* sel_fibre_norm = new TH2F("sel_fibre_norm","sel_fibre_norm",nBinsFibre,xLowFibre,xUpFibre,nBinsFibre,yLowFibre,yUpFibre);

  if(theConfiguration_.setup=="Nino"){
    std::cout<<"using NINO"<<std::endl;
  }


   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = t.LoadTree(jentry);
      if (ientry < 0) break;
      nb = t.fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if( t.mcp_time_frac50<0 ||t.mcp_max_amplitude<200)continue;
      float deltaT=t.cef3_time_at_frac50->at(1)-t.mcp_time_frac50;
      if(theConfiguration_.setup=="Nino")deltaT=t.nino_LEtime-t.mcp_time_frac50;
      if(TopologicalSelectionHelper::passesFibreTopologicalSelection(t,(theConfiguration_.setup=="Nino")) &&  t.cef3_maxAmpl->at(2) < theConfiguration_.channel2CutFibre)   //cuts on fibre position with hodos and wc. cut on cef3_maxAmpl[2] to reduce remaining events hitting the channel
	reso_histo_fibre->Fill(deltaT);
      if(TopologicalSelectionHelper::passesChannelTopologicalSelection(t,(theConfiguration_.setup=="Nino"))  && t.cef3_maxAmpl->at(2) > theConfiguration_.channel2CutChannel && t.cef3_maxAmpl->at(1)<theConfiguration_.channel1CutChannel)   //cuts on fibre position with hodos and wc. cut on cef3_maxAmpl[2] to reduce hadron contamination, cef3_maxAmpl[1] cut to reduce remaining events hitting the fibre
	reso_histo_channel->Fill(deltaT);
   }


   float reso_mean_fibre = reso_histo_fibre->GetMean();
   float reso_sigma_fibre = reso_histo_fibre->GetRMS();

   float reso_mean_channel = reso_histo_channel->GetMean();
   float reso_sigma_channel = reso_histo_channel->GetRMS();

   TH1F* reso_histo_fibre_corr[theConfiguration_.nMaxAmplCuts];
   TH1F* reso_histo_fibre_corr_Amplitude[theConfiguration_.nMaxAmplCuts];//to see behaviour of resolution as a function of ampl
   TH1F* reso_histo_channel_corr[theConfiguration_.nMaxAmplCuts];
   TH1F* reso_histo_channel_corr_Amplitude[theConfiguration_.nMaxAmplCuts];//to see behaviour of resolution as a function of ampl

   TVectorD resValueTime_fibre(theConfiguration_.nMaxAmplCuts);
   TVectorD resErrValueTime_fibre(theConfiguration_.nMaxAmplCuts);
   TVectorD resErrRelativeValueTime_fibre(theConfiguration_.nMaxAmplCuts);
   TVectorD nEntriesTime_fibre(theConfiguration_.nMaxAmplCuts);  
   TVectorD resValueAmplitude_fibre(theConfiguration_.nMaxAmplCuts);//to see behaviour of resolution as a function of ampl
   TVectorD resErrValueAmplitude_fibre(theConfiguration_.nMaxAmplCuts);//to see behaviour of resolution as a function of ampl

   TVectorD resValueTime_channel(theConfiguration_.nMaxAmplCuts);
   TVectorD resErrValueTime_channel(theConfiguration_.nMaxAmplCuts);
   TVectorD resErrRelativeValueTime_channel(theConfiguration_.nMaxAmplCuts);
   TVectorD nEntriesTime_channel(theConfiguration_.nMaxAmplCuts);  
   TVectorD resValueAmplitude_channel(theConfiguration_.nMaxAmplCuts);//to see behaviour of resolution as a function of ampl
   TVectorD resErrValueAmplitude_channel(theConfiguration_.nMaxAmplCuts);//to see behaviour of resolution as a function of ampl


   for (int i=0;i<theConfiguration_.nMaxAmplCuts;++i){
    TString icut;
    icut.Form("%d",i); 
    reso_histo_fibre_corr[i] = new TH1F("reso_histo_fibre_corr_"+icut,"reso_histo_fibre_corr_"+icut,theConfiguration_.nBinsFibre,-3*reso_sigma_fibre,3*reso_sigma_fibre);
    reso_histo_fibre_corr_Amplitude[i] = new TH1F("reso_histo_fibre_corr_Amplitude_"+icut,"reso_histo_fibre_corr_Amplitude_"+icut,200,-3*reso_sigma_fibre,3*reso_sigma_fibre);


    reso_histo_channel_corr[i] = new TH1F("reso_histo_channel_corr_"+icut,"reso_histo_channel_corr_"+icut,theConfiguration_.nBinsChannel,-3*reso_sigma_channel,3*reso_sigma_channel);
    reso_histo_channel_corr_Amplitude[i] = new TH1F("reso_histo_channel_corr_Amplitude_"+icut,"reso_histo_channel_corr_Amplitude_"+icut,200,-2*reso_sigma_channel,2*reso_sigma_channel);

   }

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = t.LoadTree(jentry);
      if (ientry < 0) break;
      nb = t.fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if( t.mcp_time_frac50<0 || t.mcp_time_frac50>400 || t.cef3_time_at_frac50->at(1)<0 || t.cef3_time_at_frac50->at(1) > 400 || t.mcp_max_amplitude<200)continue;
      //map_fibre of timing vs x and y
      float deltaTNoCorr=t.cef3_time_at_frac50->at(1)-t.mcp_time_frac50;
      if(theConfiguration_.setup=="Nino")deltaTNoCorr=t.nino_LEtime-t.mcp_time_frac50;
      float deltaT=deltaTNoCorr-(reso_mean_channel+reso_mean_fibre)/2.+reso_sigma_channel/2.;//for total maps we do an average of mean time of channel and fibre


      float X=(t.cluster_pos_corr_hodoX1+t.cluster_pos_corr_hodoX2)/2.;
      float Y=(t.cluster_pos_corr_hodoY1+t.cluster_pos_corr_hodoY2)/2.;
      float ampl=t.cef3_maxAmpl->at(1);
      if(theConfiguration_.setup=="Nino")ampl=t.nino_maxAmpl;

	if(deltaT>-2 && deltaT<2 && ampl>0 && ampl<4000){
	   timing_map_fibre_norm->Fill(X,Y);
	   timing_map_fibre->Fill(X,Y,deltaTNoCorr+10);
	   amplitude_map_fibre->Fill(X,Y,ampl);
	   timing_map_channel_norm->Fill(X,Y);
	   timing_map_channel->Fill(X,Y,deltaTNoCorr+10);
	   amplitude_map_channel->Fill(X,Y,ampl);
	   //	  std::cout<<deltaT<<std::endl;
	   amplitude_map_fibre2->Fill(X,Y,t.cef3_maxAmpl->at(2));
	   amplitude_map_fibre0->Fill(X,Y,t.cef3_maxAmpl->at(0));
	   amplitude_map_fibre0and2->Fill(X,Y,t.cef3_maxAmpl->at(0)+t.cef3_maxAmpl->at(2));
	}


	if(TopologicalSelectionHelper::passesFibreTopologicalSelection(t,(theConfiguration_.setup=="Nino"))){
	float deltaTCorr=deltaTNoCorr-reso_mean_fibre+reso_sigma_fibre/2;
	if(t.cef3_maxAmpl->at(2) < theConfiguration_.channel2CutFibre){   //cuts on fibre position with hodos and wc. cut on cef3_maxAmpl[2] to reduce remaining events hitting the channel

	if(deltaTCorr>-2 && deltaTCorr<2 && ampl>0 && ampl<4000){
	  amplitude_map_sel_fibre->Fill(X,Y,ampl);
	  timing_map_sel_fibre->Fill(X,Y,deltaTNoCorr+10);
	  sel_fibre_norm->Fill(X,Y);
	  maxAmpl_sel_fibre->Fill(ampl);
	}
	
	for (int i=0;i<theConfiguration_.nMaxAmplCuts;++i){
	  if(t.cef3_maxAmpl->at(1)>(i>0)*(theConfiguration_.startCutFibre+(i>1)*i*theConfiguration_.stepFibre))	    reso_histo_fibre_corr[i]->Fill(deltaTCorr);
	  if(t.cef3_maxAmpl->at(1)>theConfiguration_.startCutFibre+i*theConfiguration_.stepAmplFibre && t.cef3_maxAmpl->at(1)<(theConfiguration_.startCutFibre+(i+1)*theConfiguration_.stepAmplFibre)) reso_histo_fibre_corr_Amplitude[i]->Fill(deltaTCorr);
	}
	}
	}else if(TopologicalSelectionHelper::passesChannelTopologicalSelection(t,(theConfiguration_.setup=="Nino"))){
	float deltaTCorr=deltaTNoCorr-reso_mean_channel+reso_sigma_channel/2;
	if(t.cef3_maxAmpl->at(2) > theConfiguration_.channel2CutChannel && t.cef3_maxAmpl->at(1)<theConfiguration_.channel1CutChannel){ //cuts on channel position with hodos and wc. cut on cef3_maxAmpl[2] to reduce hadron contamination
	if(deltaTCorr>-2 && deltaTCorr<2 && ampl>0 && ampl<4000){
	  amplitude_map_sel_channel->Fill(X,Y,ampl);  
	  timing_map_sel_channel->Fill(X,Y,deltaTNoCorr+10);
	  sel_channel_norm->Fill(X,Y);
	  maxAmpl_sel_channel->Fill(ampl);
	  if(ampl>100.){
	    amplitude_map_sel_channel_ampl_cut->Fill(X,Y,ampl);
	    timing_map_sel_channel_ampl_cut->Fill(X,Y,deltaTNoCorr+10);
	    amplitude_map_fibre2_sel_channel_ampl_cut->Fill(X,Y,t.cef3_maxAmpl->at(2));
	    sel_channel_ampl_cut_norm->Fill(X,Y);
	    if(deltaTNoCorr>-5.5)sel_channel_ampl_and_time_cut_norm->Fill(X,Y);
	  }

	}

	for (int i=0;i<theConfiguration_.nMaxAmplCuts;++i){
	  if(t.cef3_maxAmpl->at(1)>(i>0)*(theConfiguration_.startCutChannel+(i>1)*i*theConfiguration_.stepChannel))reso_histo_channel_corr[i]->Fill(deltaTNoCorr-reso_mean_channel+reso_sigma_channel/2);
	  if(t.cef3_maxAmpl->at(1)>theConfiguration_.startCutChannel+i*theConfiguration_.stepAmplChannel && t.cef3_maxAmpl->at(1)<(theConfiguration_.startCutChannel+(i+1)*theConfiguration_.stepAmplChannel)) reso_histo_channel_corr_Amplitude[i]->Fill(deltaTNoCorr-reso_mean_channel+reso_sigma_channel/2);
	}

      }
      }
      
   }//nentries

   timing_map_fibre->Divide(timing_map_fibre_norm);
   //   timing_map_fibre->SetAxisRange(theConfiguration_.rangeXLow+1.5,theConfiguration_.rangeXUp-1,"Z")
   //   timing_map_fibre->SetAxisRange(-2,2,"Z");
   timing_map_fibre->GetYaxis()->SetTitle("Y [mm]");
   timing_map_fibre->GetXaxis()->SetTitle("X [mm]");
   amplitude_map_fibre->Divide(timing_map_fibre_norm);
   amplitude_map_fibre->GetYaxis()->SetTitle("Y [mm]");
   amplitude_map_fibre->GetXaxis()->SetTitle("X [mm]");

   amplitude_map_sel_fibre->Divide(sel_fibre_norm);
   amplitude_map_sel_fibre->GetYaxis()->SetTitle("Y [mm]");
   amplitude_map_sel_fibre->GetXaxis()->SetTitle("X [mm]");

   timing_map_sel_fibre->Divide(sel_fibre_norm);
   timing_map_sel_fibre->GetYaxis()->SetTitle("Y [mm]");
   timing_map_sel_fibre->GetXaxis()->SetTitle("X [mm]");


   amplitude_map_sel_channel->Divide(sel_channel_norm);
   amplitude_map_sel_channel->GetYaxis()->SetTitle("Y [mm]");
   amplitude_map_sel_channel->GetXaxis()->SetTitle("X [mm]");

   amplitude_map_sel_channel_ampl_cut->Divide(sel_channel_ampl_cut_norm);
   amplitude_map_sel_channel_ampl_cut->GetYaxis()->SetTitle("Y [mm]");
   amplitude_map_sel_channel_ampl_cut->GetXaxis()->SetTitle("X [mm]");


   timing_map_sel_channel->Divide(sel_channel_norm);
   timing_map_sel_channel->GetYaxis()->SetTitle("Y [mm]");
   timing_map_sel_channel->GetXaxis()->SetTitle("X [mm]");

   timing_map_sel_channel_ampl_cut->Divide(sel_channel_ampl_cut_norm);
   timing_map_sel_channel_ampl_cut->GetYaxis()->SetTitle("Y [mm]");
   timing_map_sel_channel_ampl_cut->GetXaxis()->SetTitle("X [mm]");


   amplitude_map_fibre2->Divide(timing_map_channel_norm);
   amplitude_map_fibre2->GetYaxis()->SetTitle("Y [mm]");
   amplitude_map_fibre2->GetXaxis()->SetTitle("X [mm]");

   amplitude_map_fibre2_sel_channel_ampl_cut->Divide(sel_channel_ampl_cut_norm);
   amplitude_map_fibre2_sel_channel_ampl_cut->GetYaxis()->SetTitle("Y [mm]");
   amplitude_map_fibre2_sel_channel_ampl_cut->GetXaxis()->SetTitle("X [mm]");


   amplitude_map_fibre0->Divide(timing_map_channel_norm);
   amplitude_map_fibre0->GetYaxis()->SetTitle("Y [mm]");
   amplitude_map_fibre0->GetXaxis()->SetTitle("X [mm]");

   amplitude_map_fibre0and2->Divide(timing_map_channel_norm);
   amplitude_map_fibre0and2->GetYaxis()->SetTitle("Y [mm]");
   amplitude_map_fibre0and2->GetXaxis()->SetTitle("X [mm]");


   timing_map_channel->Divide(timing_map_channel_norm);
   //   timing_map_channel->SetAxisRange(theConfiguration_.rangeXLow+1.5,theConfiguration_.rangeXUp-1,"Z");
   //   timing_map_channel->SetAxisRange(-2,2,"Z");
   timing_map_channel->GetYaxis()->SetTitle("Y [mm]");
   timing_map_channel->GetXaxis()->SetTitle("X [mm]");
   amplitude_map_channel->Divide(timing_map_channel_norm);
   amplitude_map_channel->GetYaxis()->SetTitle("Y [mm]");
   amplitude_map_channel->GetXaxis()->SetTitle("X [mm]");

   maxAmpl_sel_channel->SetLineWidth(2);
   maxAmpl_sel_channel->GetXaxis()->SetTitle("Signal Amplitude [ADC channel]");

   maxAmpl_sel_channel->SetLineWidth(2);
   maxAmpl_sel_channel->SetLineColor(kRed);
   maxAmpl_sel_channel->GetXaxis()->SetTitle("Signal Amplitude [ADC channel]");


   // this is the dir in which the plots will be saved:
  std::string constDirName = "plots_timing_";
  constDirName+=theConfiguration_.setup;
  system(Form("mkdir -p %s", constDirName.c_str()));
  TString dir(constDirName);


  ////////////////////////////////////////////fibre plots//////////////////////////////////////////////////////

  TCanvas c1;
  reso_histo_fibre->Draw();
  c1.SaveAs(dir+"/reso_histo_fibre_"+runNumberString+".png");
  c1.SaveAs(dir+"/reso_histo_fibre_"+runNumberString+".pdf");

  for (int i=0;i<theConfiguration_.nMaxAmplCuts;++i){
    std::cout<<i<<" "<<reso_histo_fibre_corr[i]->GetEntries()<<" "<<reso_histo_fibre_corr[i]->Integral(reso_histo_fibre_corr[i]->FindBin(-0.5),reso_histo_fibre_corr[i]->FindBin(0.3))<<std::endl;
    if(reso_histo_fibre_corr[i]->Integral(reso_histo_fibre_corr[i]->FindBin(-0.5),reso_histo_fibre_corr[i]->FindBin(0.3))<5) continue;


    TString icut;
    icut.Form("%d",i); 
    

    c1.Clear();
    bool fitWithGaussian = false;

    std::string ytitle = Form("Events");
      std::string energy(Form("%.0f", t.beamEnergy));
    TPaveText* pave = DrawTools::getLabelTop_expOnXaxis(energy+" GeV Electron Beam");

    if(fitWithGaussian){
      reso_histo_fibre_corr[i]->Fit("gaus","","",-0.5,0.3);
      std::cout<<"fit done"<<std::endl;
      //  double *par=reso_histo_fibre_corr[i]->GetFunction("gaus")->GetParameters();
      TF1* fgaus=reso_histo_fibre_corr[i]->GetFunction("gaus");
      
      fgaus->SetLineWidth(2.);
      fgaus->SetLineColor(kBlue);
      reso_histo_fibre_corr[i]->GetXaxis()->SetTitle("time_{Fibre}-time_{mcp} [ns]");
      float binWidth =reso_histo_fibre_corr[i]->GetXaxis()->GetBinWidth(1);
      //  std::string ytitle = Form("Events / %.0f ps",binWidth*1.e3); 

      reso_histo_fibre_corr[i]->SetYTitle(ytitle.c_str());
      reso_histo_fibre_corr[i]->Draw();


      pave->Draw("same");
      
      fgaus->Draw("same");
      TLegend* leg_gauss = new TLegend(0.6, 0.8, 0.8, 0.9);
      leg_gauss->SetTextSize(0.036);
      leg_gauss->AddEntry(  (TObject*)0 ,Form("#sigma = %.0f #pm %.0f ps", fgaus->GetParameter(2)*1.e3, fgaus->GetParError(2)*1.e3), "");
      leg_gauss->SetFillColor(0);
      leg_gauss->Draw("same");
      

      
      c1.SaveAs(dir+"/reso_histo_fibre_corr_gaus_"+icut+"_"+runNumberString+".png");
      c1.SaveAs(dir+"/reso_histo_fibre_corr_gaus_"+icut+"_"+runNumberString+".pdf");
    }

  //-----------------fit with cruijff ------------------------
  TH1F* histo;

  histo=reso_histo_fibre_corr[i];
  nEntriesTime_fibre[i]=histo->GetEntries();
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
  TLegend* lego = new TLegend(0.57, 0.8, 0.89, 0.9); lego->SetTextAlign(32);
  lego->SetTextSize(0.036);
  lego->AddEntry(  (TObject*)0 ,Form("#sigma = %.1f #pm %.1f ps", rms*1.e3, rmsErr*1.e3), "");

  lego->SetFillColor(0);
  lego->Draw("same");
  
  pave->Draw("same");

  cans->SaveAs(dir+"/reso_histo_fibre_cruijff_"+icut+"_"+runNumberString+".png");
  cans->SaveAs(dir+"/reso_histo_fibre_cruijff_"+icut+"_"+runNumberString+".pdf");


  resValueTime_fibre[i]=rms*1.e3;
  resErrValueTime_fibre[i]=rmsErr*1.e3;

  if(rms!=0)  resErrRelativeValueTime_fibre[i]=rmsErr/rms;

  resValueAmplitude_fibre[i]=reso_histo_fibre_corr_Amplitude[i]->GetRMS()*1.e3; 
  resErrValueAmplitude_fibre[i]=reso_histo_fibre_corr_Amplitude[i]->GetRMSError()*1.e3; 

  //////------------ fit with crystal ball-------
  bool fitWithCB=false;
  if(fitWithCB){
    RooRealVar x("x","deltaT", fitmin, fitmax);
    RooDataHist data("data","dataset with x",x,RooFit::Import(*histo) );
    
    RooRealVar meanr("meanr","Mean",peakpos,peakpos-3*sigma, peakpos+3*sigma);
    RooRealVar width("width","#sigma",sigma , 0, 5.*sigma);
    RooRealVar A("A","Dist",2., 0.0, 7.0);
    RooRealVar N("N","Deg",5, 0.0, 10);
    int ndf;

    RooCBShape fit_fct("fit_fct","fit_fct",x,meanr,width,A,N); ndf = 4;
    fit_fct.fitTo(data);
    fit_fct.plotOn(frame,RooFit::LineColor(4));//this will show fit overlay 

    RooPlot* frame;
    std::string ytitle = Form("Events");
    frame = x.frame("Title");
    frame->SetXTitle("time_{fibre}-time_{mcp} [ns]");
    frame->SetYTitle(ytitle.c_str());
    
    data.plotOn(frame);  //this will show histogram data points on canvas 
    fit_fct.plotOn(frame);//this will show fit overlay on canvas  

    rms = width.getVal();
    rmsErr = width.getError();
    TCanvas* cans = new TCanvas();
    cans->cd();
    frame->Draw();

    TLegend* lego = new TLegend(0.57, 0.7, 0.8, 0.9);
    lego->SetTextSize(0.036);
    lego->AddEntry(  (TObject*)0 ,Form("#sigma = %.1f #pm %.1f ps", rms*1.e3, rmsErr*1.e3), "");
    lego->SetFillColor(0);
    lego->Draw("same");

    pave->Draw("same");

  cans->SaveAs(dir+"/reso_histo_fibre_cb_"+icut+"_"+runNumberString+".png");
  cans->SaveAs(dir+"/reso_histo_fibre_cb_"+icut+"_"+runNumberString+".pdf");


  }


  }

  ////////////////////////////////////////////channel plots//////////////////////////////////////////////////////

  c1.Clear();
  reso_histo_channel->Draw();
  c1.SaveAs(dir+"/reso_histo_channel_"+runNumberString+".png");
  c1.SaveAs(dir+"/reso_histo_channel_"+runNumberString+".pdf");

  for (int i=0;i<theConfiguration_.nMaxAmplCuts;++i){
    if(reso_histo_channel_corr[i]->Integral(reso_histo_channel_corr[i]->FindBin(-0.5),reso_histo_channel_corr[i]->FindBin(0.3))<5) continue;
    
    TString icut;
    icut.Form("%d",i); 
    

    c1.Clear();
    //-----------------fit with cruijff ------------------------
    TH1F* histo;
    
    histo=reso_histo_channel_corr[i];
    nEntriesTime_channel[i]=histo->GetEntries();
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
  
    frame = x.frame("Title");
    frame->SetXTitle("time_{fibre}-time_{mcp} [ns]");

    std::string ytitle = Form("Events");
    frame->SetYTitle(ytitle.c_str());
  
    data.plotOn(frame);  //this will show histogram data points on canvas 
    fit_fct.plotOn(frame);//this will show fit overlay on canvas  

    double rms,rmsErr;
    rms = (widthL.getVal()+widthR.getVal())/2;
    rmsErr = 0.5*sqrt(widthL.getError()*widthL.getError()+widthR.getError()*widthR.getError());

    TCanvas* cans = new TCanvas();
    cans->cd();
    frame->Draw();
    TLegend* lego = new TLegend(0.57, 0.8, 0.89, 0.9); lego->SetTextAlign(32);
    lego->SetTextSize(0.036);
    lego->AddEntry(  (TObject*)0 ,Form("#sigma = %.1f #pm %.1f ps", rms*1.e3, rmsErr*1.e3), "");

    lego->SetFillColor(0);
    lego->Draw("same");

    std::string energy(Form("%.0f", t.beamEnergy));  
    TPaveText* pave = DrawTools::getLabelTop_expOnXaxis(energy+" GeV Electron Beam");
    pave->Draw("same");

    cans->SaveAs(dir+"/reso_histo_channel_cruijff_"+icut+"_"+runNumberString+".png");
    cans->SaveAs(dir+"/reso_histo_channel_cruijff_"+icut+"_"+runNumberString+".pdf");


    resValueTime_channel[i]=rms*1.e3;
    resErrValueTime_channel[i]=rmsErr*1.e3;

    if(rms!=0)  resErrRelativeValueTime_channel[i]=rmsErr/rms;

    resValueAmplitude_channel[i]=reso_histo_channel_corr_Amplitude[i]->GetRMS()*1.e3; 
    resErrValueAmplitude_channel[i]=reso_histo_channel_corr_Amplitude[i]->GetRMSError()*1.e3; 


  }



  std::string outFileName;
  outFileName = dir+"/timingStudiesSiPM_"+tag+"_"+runName+".root";
  TFile* outFile = TFile::Open(outFileName.c_str(),"recreate");

  timing_map_fibre->Write();
  amplitude_map_fibre->Write();
  amplitude_map_fibre2->Write();
  amplitude_map_fibre0->Write();
  amplitude_map_fibre0and2->Write();
  timing_map_fibre_norm->Write();

  timing_map_channel->Write();
  amplitude_map_channel->Write();
  timing_map_channel_norm->Write();

  amplitude_map_sel_fibre->Write();
  amplitude_map_sel_channel->Write();

  sel_fibre_norm->Write();
  sel_channel_norm->Write();
  sel_channel_ampl_cut_norm->Write();
  sel_channel_ampl_and_time_cut_norm->Write();
  timing_map_sel_fibre->Write();
  timing_map_sel_channel->Write();

  maxAmpl_sel_fibre->Write();
  maxAmpl_sel_channel->Write();

  amplitude_map_sel_channel_ampl_cut->Write();
  timing_map_sel_channel_ampl_cut->Write();
  amplitude_map_fibre2_sel_channel_ampl_cut->Write();

  resValueTime_fibre.Write("resValueTime_fibre");
  resValueAmplitude_fibre.Write("resValueAmplitude_fibre");
  resErrValueAmplitude_fibre.Write("resErrValueAmplitude_fibre");
  resErrValueTime_fibre.Write("resErrValueTime_fibre");
  resErrRelativeValueTime_fibre.Write("resErrRelativeValueTime_fibre");
  nEntriesTime_fibre.Write("nEntriesTime_fibre");

  resValueTime_channel.Write("resValueTime_channel");
  resValueAmplitude_channel.Write("resValueAmplitude_channel");
  resErrValueAmplitude_channel.Write("resErrValueAmplitude_channel");
  resErrValueTime_channel.Write("resErrValueTime_channel");
  resErrRelativeValueTime_channel.Write("resErrRelativeValueTime_channel");
  nEntriesTime_channel.Write("nEntriesTime_channel");


  outFile->Write();
  outFile->Close();

  return 0;

}

timingPlots_Config_t readConfiguration(std::string configName, float beamEnergy){

  std::string energy(Form("%.0f", beamEnergy));  

  Configurator* configurator_ = new Configurator();
  std::string fileName = Form ("./config_timing/%s_%sGeV.xml",configName.c_str(),energy.c_str());
  configurator_->xmlFileName=fileName.c_str();
  configurator_->Init();

  timingPlots_Config_t conf;
  
   conf.nMaxAmplCuts= Configurator::GetInt(Configurable::getElementContent(*configurator_,"nMaxAmplCuts",configurator_->root_element));
   conf.setup=Configurable::getElementContent (*configurator_, "setup",configurator_->root_element) ;

   conf.startCutFibre= Configurator::GetInt(Configurable::getElementContent(*configurator_,"startCutFibre",configurator_->root_element));
   conf.stepFibre= Configurator::GetInt(Configurable::getElementContent(*configurator_,"stepFibre",configurator_->root_element));
   conf.stepAmplFibre= Configurator::GetInt(Configurable::getElementContent(*configurator_,"stepAmplFibre",configurator_->root_element));

   conf.startCutChannel= Configurator::GetInt(Configurable::getElementContent(*configurator_,"startCutChannel",configurator_->root_element));
   conf.stepChannel= Configurator::GetInt(Configurable::getElementContent(*configurator_,"stepChannel",configurator_->root_element));
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

//  LocalWords:  maxAmpl
