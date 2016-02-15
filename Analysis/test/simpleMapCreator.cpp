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
    std::cout << "./simpleMapCreator ([tag]) ([config])" << std::endl;
    exit(12345);
  }

  std::cout<<config<<std::endl;
  
  theConfiguration_=readConfiguration(config);
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

  float xLowChannel=-5.5, xUpChannel=4, yLowChannel=-2, yUpChannel=4.5;
  float xLowFibre=-5.5,xUpFibre=-2.5, yLowFibre=-3, yUpFibre=0;
  float amplChannelLow=0, amplChannelUp=300.;
  float amplFibreLow=0, amplFibreUp=700.;
  float timeLow=0, timeUp=100;
  float timeMcpLow=0, timeMcpUp=100;

  float shift=10;
  bool isNino=false;

  if(theConfiguration_.setup=="Nino"){
    xLowChannel=-4.5; xUpChannel=4; yLowChannel=1; yUpChannel=4.5;
    xLowFibre=-5; xUpFibre=0.3; yLowFibre=0.5; yUpFibre=4.5;
    shift=11.5;
    isNino=true;
    amplChannelLow=0, amplChannelUp=65.;
    amplFibreLow=0, amplFibreUp=65.;
    timeLow=30;
    timeMcpLow=30;
  }



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

  int nBinsChannel=30, nBinsFibre=20;

  TH2F* amplitude_map_sel_channel = new TH2F("amplitude_map_sel_channel","amplitude_map_sel_channel",nBinsChannel,xLowChannel,xUpChannel,nBinsChannel,yLowChannel,yUpChannel);
  TH2F* amplitude_map_sel_channel_ampl_cut = new TH2F("amplitude_map_sel_channel_ampl_cut","amplitude_map_sel_channel_ampl_cut",nBinsChannel,xLowChannel,xUpChannel,nBinsChannel,yLowChannel,yUpChannel);
  TH2F* amplitude_map_fibre2_sel_channel_ampl_cut = new TH2F("amplitude_map_fibre2_sel_channel_ampl_cut","amplitude_map_fibre2_sel_channel_ampl_cut",nBinsChannel,xLowChannel,xUpChannel,nBinsChannel,yLowChannel,yUpChannel);
  TH2F* amplitude_map_sel_fibre = new TH2F("amplitude_map_sel_fibre","amplitude_map_sel_fibre",nBinsFibre,xLowFibre,xUpFibre,nBinsFibre,yLowFibre,yUpFibre);

  TH2F* timing_map_sel_channel = new TH2F("timing_map_sel_channel","timing_map_sel_channel",nBinsChannel,xLowChannel,xUpChannel,nBinsChannel,yLowChannel,yUpChannel);
  TH2F* timing_map_sel_channel_ampl_cut = new TH2F("timing_map_sel_channel_ampl_cut","timing_map_sel_channel_ampl_cut",nBinsChannel,xLowChannel,xUpChannel,nBinsChannel,yLowChannel,yUpChannel);
  TH2F* timing_map_sel_channel_ampl_and_time_cut = new TH2F("timing_map_sel_channel_ampl_and_time_cut","timing_map_sel_channel_ampl_and_time_cut",nBinsChannel,xLowChannel,xUpChannel,nBinsChannel,yLowChannel,yUpChannel);
  TH2F* timing_map_sel_fibre = new TH2F("timing_map_sel_fibre","timing_map_sel_fibre",nBinsFibre,xLowFibre,xUpFibre,nBinsFibre,yLowFibre,yUpFibre);
  TH2F* timing_map_sel_fibre_ampl_cut = new TH2F("timing_map_sel_fibre_ampl_cut","timing_map_sel_fibre_ampl_cut",nBinsFibre,xLowFibre,xUpFibre,nBinsFibre,yLowFibre,yUpFibre);
        
  TH2F* timing_map_sel_fibre_ampl_and_time_cut = new TH2F("timing_map_sel_fibre_ampl_and_time_cut","timing_map_sel_fibre_ampl_and_time_cut",nBinsFibre,xLowFibre,xUpFibre,nBinsFibre,yLowFibre,yUpFibre);
  TH2F* sel_fibre_ampl_cut_norm = new TH2F("sel_fibre_ampl_cut_norm","sel_fibre_ampl_cut_norm",nBinsFibre,xLowFibre,xUpFibre,nBinsFibre,yLowFibre,yUpFibre);
  TH2F* sel_fibre_ampl_and_time_cut_norm = new TH2F("sel_fibre_ampl_and_time_cut_norm","sel_fibre_ampl_and_time_cut_norm",nBinsFibre,xLowFibre,xUpFibre,nBinsFibre,yLowFibre,yUpFibre);

  TH2F* sel_channel_norm = new TH2F("sel_channel_norm","sel_channel_norm",nBinsChannel,xLowChannel,xUpChannel,nBinsChannel,yLowChannel,yUpChannel);
  TH2F* sel_channel_ampl_cut_norm = new TH2F("sel_channel_ampl_cut_norm","sel_channel_ampl_cut_norm",nBinsChannel,xLowChannel,xUpChannel,nBinsChannel,yLowChannel,yUpChannel);
  TH2F* sel_channel_ampl_and_time_cut_norm = new TH2F("sel_channel_ampl_and_time_cut_norm","sel_channel_ampl_and_time_cut_norm",nBinsChannel,xLowChannel,xUpChannel,nBinsChannel,yLowChannel,yUpChannel);
  TH2F* sel_fibre_norm = new TH2F("sel_fibre_norm","sel_fibre_norm",nBinsFibre,xLowFibre,xUpFibre,nBinsFibre,yLowFibre,yUpFibre);

  //correlation
  TH2F* timeDiffVsAmpl = new TH2F("timeDiffVsAmpl","timeDiffVsAmpl",300,amplChannelLow,amplChannelUp,100,theConfiguration_.rangeXLow+10, theConfiguration_.rangeXUp+10);
  TH2F* timeDiffVsAmpl_fibre = new TH2F("timeDiffVsAmpl_fibre","timeDiffVsAmpl_fibre",500,amplFibreLow,amplFibreUp,100,theConfiguration_.rangeXLow+10, theConfiguration_.rangeXUp+10);
  TH2F* timeDiffVsTime = new TH2F("timeDiffVsTime","timeDiffVsTime",100,timeLow,timeUp,100, theConfiguration_.rangeXLow+10,theConfiguration_.rangeXUp+10);
  TH2F* timeDiffVsTimeMcp = new TH2F("timeDiffVsTimeMcp","timeDiffVsTimeMcp",100,timeMcpLow,timeMcpUp,100, theConfiguration_.rangeXLow+10,theConfiguration_.rangeXUp+10);
  TH2F* timeDiffVsTimeMcp_amplCut = new TH2F("timeDiffVsTimeMcp_amplCut","timeDiffVsTimeMcp_amplCut",100,timeMcpLow,timeMcpUp,100, theConfiguration_.rangeXLow+10,theConfiguration_.rangeXUp+10);
  TH2F* timeDiffVsTime_amplCut = new TH2F("timeDiffVsTime_amplCut","timeDiffVsTime_amplCut",100,timeLow,timeUp,100, theConfiguration_.rangeXLow+10,theConfiguration_.rangeXUp+10);
  TH2F* timeDiffVsTimeMaxMinusFifty = new TH2F("timeDiffVsTimeMaxMinusFifty","timeDiffVsTimeMaxMinusFifty",100,-20,20,100, theConfiguration_.rangeXLow+10,theConfiguration_.rangeXUp+10);
  TH2F* timeDiffVsTimeMaxMinusFifty_amplCut = new TH2F("timeDiffVsTimeMaxMinusFifty_amplCut","timeDiffVsTimeMaxMinusFifty_amplCut",100,-20,20,100, theConfiguration_.rangeXLow+10,theConfiguration_.rangeXUp+10);
  TH2F* timeDiffVsTimeMaxMinusFifty_amplCut_zoom = new TH2F("timeDiffVsTimeMaxMinusFifty_amplCut_zoom","timeDiffVsTimeMaxMinusFifty_amplCut_zoom",100,2,6,100, theConfiguration_.rangeXLow+10,theConfiguration_.rangeXUp+10);
  TH1F* timeMaxMinusFifty = new TH1F("timeMaxMinusFifty","timeMaxMinusFifty",100,-20,20);
  TH1F* timeMaxMinusFifty_amplCut = new TH1F("timeMaxMinusFifty_amplCut","timeMaxMinusFifty_amplCut",100,-20,20);
  TH1F* timeMaxMinusFifty_timeCut = new TH1F("timeMaxMinusFifty_timeCut","timeMaxMinusFifty_timeCut",100,-20,20);





  int dummyCounter=0;
  float lowerBeamEnergy=0.;//cuts are sliding in energy
   Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = t.LoadTree(jentry);
      if (ientry < 0) break;
      nb = t.fChain->GetEntry(jentry);   nbytes += nb;
      if(jentry==0)lowerBeamEnergy=t.beamEnergy;
      if( t.mcp_time_frac50<0 || t.mcp_time_frac50>400 || t.cef3_time_at_frac50->at(1)<0 || t.cef3_time_at_frac50->at(1) > 400 || t.mcp_max_amplitude<200)continue;
      float deltaTNoCorr=t.cef3_time_at_frac50->at(1)-t.mcp_time_frac50;
      if(theConfiguration_.setup=="Nino")deltaTNoCorr=t.nino_LEtime-t.mcp_time_frac50;

      float X=(t.cluster_pos_corr_hodoX1+t.cluster_pos_corr_hodoX2)/2.;
      float Y=(t.cluster_pos_corr_hodoY1+t.cluster_pos_corr_hodoY2)/2.;
      float ampl=t.cef3_maxAmpl->at(1);
      float amplForCut=t.cef3_maxAmpl->at(1);
      float amplCut=theConfiguration_.amplCut;
      if(theConfiguration_.setup=="Nino"){
	ampl=t.nino_maxAmpl;
      }
      
      if(deltaTNoCorr>theConfiguration_.rangeXLow && deltaTNoCorr<theConfiguration_.rangeXUp && ampl>0 && ampl<4000){
	timing_map_fibre_norm->Fill(X,Y);
	timing_map_fibre->Fill(X,Y,deltaTNoCorr+shift);
	amplitude_map_fibre->Fill(X,Y,ampl);
	timing_map_channel_norm->Fill(X,Y);
	timing_map_channel->Fill(X,Y,deltaTNoCorr+shift);
	amplitude_map_channel->Fill(X,Y,ampl);
	//	  std::cout<<deltaT<<std::endl;
	amplitude_map_fibre2->Fill(X,Y,t.cef3_maxAmpl->at(2));
	amplitude_map_fibre0->Fill(X,Y,t.cef3_maxAmpl->at(0));
	amplitude_map_fibre0and2->Fill(X,Y,t.cef3_maxAmpl->at(0)+t.cef3_maxAmpl->at(2));
      }


	
      if(TopologicalSelectionHelper::passesFibreTopologicalSelection(t,isNino) &&  t.cef3_maxAmpl->at(2) < theConfiguration_.channel2CutFibre*t.beamEnergy/lowerBeamEnergy){ 
	if(deltaTNoCorr>theConfiguration_.rangeXLow && deltaTNoCorr<theConfiguration_.rangeXUp && ampl>0 && ampl<4000){
	  timeDiffVsAmpl_fibre->Fill(ampl,deltaTNoCorr+shift);
	  if(isNino && t.nino_maxAmpl<32)continue;
	  amplitude_map_sel_fibre->Fill(X,Y,ampl);
	  timing_map_sel_fibre->Fill(X,Y,deltaTNoCorr+shift);
	  sel_fibre_norm->Fill(X,Y);
	  maxAmpl_sel_fibre->Fill(ampl);
	  if(amplForCut>amplCut){
	    timing_map_sel_fibre_ampl_cut->Fill(X,Y,deltaTNoCorr+shift);
	    sel_fibre_ampl_cut_norm->Fill(X,Y);
	    if(deltaTNoCorr>-5.4){
	    timing_map_sel_fibre_ampl_and_time_cut->Fill(X,Y,deltaTNoCorr+shift);
	    sel_fibre_ampl_and_time_cut_norm->Fill(X,Y);
	    }
	  }

	}

	  
      }
      if(TopologicalSelectionHelper::passesChannelTopologicalSelection(t,isNino)  && t.cef3_maxAmpl->at(2) > theConfiguration_.channel2CutChannel*t.beamEnergy/lowerBeamEnergy && t.cef3_maxAmpl->at(1)<theConfiguration_.channel1CutChannel){   //cuts on fibre position with hodos and wc. cut on cef3_maxAmpl[2] to reduce hadron contamination, cef3_maxAmpl[1] cut to reduce remaining events hitting the fibre
	    if(deltaTNoCorr>theConfiguration_.rangeXLow && deltaTNoCorr<theConfiguration_.rangeXUp && ampl>0 && ampl<4000){
	      timeDiffVsAmpl->Fill(ampl,deltaTNoCorr+shift);
	      timeDiffVsTime->Fill(t.cef3_time_at_frac50->at(1),deltaTNoCorr+shift);
	      timeDiffVsTimeMcp->Fill(t.mcp_time_frac50,deltaTNoCorr+shift);
	      timeMaxMinusFifty->Fill(t.cef3_maxAmpl_time->at(1)-t.cef3_time_at_frac50->at(1));
	      timeDiffVsTimeMaxMinusFifty->Fill(t.cef3_maxAmpl_time->at(1)-t.cef3_time_at_frac50->at(1),deltaTNoCorr+shift);
	      if(isNino && t.nino_maxAmpl<32)continue;
	      amplitude_map_sel_channel->Fill(X,Y,ampl);  
	      timing_map_sel_channel->Fill(X,Y,deltaTNoCorr+shift);
	      sel_channel_norm->Fill(X,Y);
	      maxAmpl_sel_channel->Fill(ampl);

	      if(amplForCut>amplCut){
		timeDiffVsTimeMcp_amplCut->Fill(t.mcp_time_frac50,deltaTNoCorr+shift);
		timeDiffVsTime_amplCut->Fill(t.cef3_time_at_frac50->at(1),deltaTNoCorr+shift);
		timeMaxMinusFifty_amplCut->Fill(t.cef3_maxAmpl_time->at(1)-t.cef3_time_at_frac50->at(1));
		amplitude_map_sel_channel_ampl_cut->Fill(X,Y,ampl);
		timing_map_sel_channel_ampl_cut->Fill(X,Y,deltaTNoCorr+shift);
		amplitude_map_fibre2_sel_channel_ampl_cut->Fill(X,Y,t.cef3_maxAmpl->at(2));
		sel_channel_ampl_cut_norm->Fill(X,Y);
		timeDiffVsTimeMaxMinusFifty_amplCut->Fill(t.cef3_maxAmpl_time->at(1)-t.cef3_time_at_frac50->at(1),deltaTNoCorr+shift);
		timeDiffVsTimeMaxMinusFifty_amplCut_zoom->Fill(t.cef3_maxAmpl_time->at(1)-t.cef3_time_at_frac50->at(1),deltaTNoCorr+shift);
		if(deltaTNoCorr>-5.4){
		  timeMaxMinusFifty_timeCut->Fill(t.cef3_maxAmpl_time->at(1)-t.cef3_time_at_frac50->at(1));
		  timing_map_sel_channel_ampl_and_time_cut->Fill(X,Y,deltaTNoCorr+shift);
		  sel_channel_ampl_and_time_cut_norm->Fill(X,Y);
		}
	      }
	    }

	  
	  }

      
   }

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
   if(isNino) amplitude_map_sel_fibre->SetAxisRange(30,60,"Z");

   timing_map_sel_fibre->Divide(sel_fibre_norm);
   timing_map_sel_fibre->GetYaxis()->SetTitle("Y [mm]");
   timing_map_sel_fibre->GetXaxis()->SetTitle("X [mm]");
   timing_map_sel_fibre->SetAxisRange(3.5,5.5,"Z");

   amplitude_map_sel_channel->Divide(sel_channel_norm);
   amplitude_map_sel_channel->GetYaxis()->SetTitle("Y [mm]");
   amplitude_map_sel_channel->GetXaxis()->SetTitle("X [mm]");
   if(isNino) amplitude_map_sel_channel->SetAxisRange(30,60,"Z");

   amplitude_map_sel_channel_ampl_cut->Divide(sel_channel_ampl_cut_norm);
   amplitude_map_sel_channel_ampl_cut->GetYaxis()->SetTitle("Y [mm]");
   amplitude_map_sel_channel_ampl_cut->GetXaxis()->SetTitle("X [mm]");
   if(!isNino)   amplitude_map_sel_channel_ampl_cut->SetAxisRange(100,250,"Z");
   else amplitude_map_sel_channel_ampl_cut->SetAxisRange(30,60,"Z");

   timing_map_sel_channel->Divide(sel_channel_norm);
   timing_map_sel_channel->GetYaxis()->SetTitle("Y [mm]");
   timing_map_sel_channel->GetXaxis()->SetTitle("X [mm]");
    timing_map_sel_channel->SetAxisRange(3.5,5.5,"Z");

   timing_map_sel_channel_ampl_cut->Divide(sel_channel_ampl_cut_norm);
   timing_map_sel_channel_ampl_cut->GetYaxis()->SetTitle("Y [mm]");
   timing_map_sel_channel_ampl_cut->GetXaxis()->SetTitle("X [mm]");
   timing_map_sel_channel_ampl_cut->SetAxisRange(3.5,5.5,"Z");

   timing_map_sel_fibre_ampl_cut->Divide(sel_fibre_ampl_cut_norm);
   timing_map_sel_fibre_ampl_cut->GetYaxis()->SetTitle("Y [mm]");
   timing_map_sel_fibre_ampl_cut->GetXaxis()->SetTitle("X [mm]");
   timing_map_sel_fibre_ampl_cut->SetAxisRange(3.5,5,"Z");

   timing_map_sel_fibre_ampl_and_time_cut->Divide(sel_fibre_ampl_and_time_cut_norm);
   timing_map_sel_fibre_ampl_and_time_cut->GetYaxis()->SetTitle("Y [mm]");
   timing_map_sel_fibre_ampl_and_time_cut->GetXaxis()->SetTitle("X [mm]");
   timing_map_sel_fibre_ampl_and_time_cut->SetAxisRange(3.5,5,"Z");

   sel_channel_ampl_cut_norm->GetYaxis()->SetTitle("Y [mm]");
   sel_channel_ampl_cut_norm->GetXaxis()->SetTitle("X [mm]");

   sel_channel_ampl_and_time_cut_norm->GetYaxis()->SetTitle("Y [mm]");
   sel_channel_ampl_and_time_cut_norm->GetXaxis()->SetTitle("X [mm]");


   timing_map_sel_channel_ampl_and_time_cut->Divide(sel_channel_ampl_and_time_cut_norm);
   timing_map_sel_channel_ampl_and_time_cut->GetYaxis()->SetTitle("Y [mm]");
   timing_map_sel_channel_ampl_and_time_cut->GetXaxis()->SetTitle("X [mm]");
   timing_map_sel_channel_ampl_and_time_cut->SetAxisRange(3.5,5.5,"Z");


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

   timeDiffVsAmpl->GetXaxis()->SetTitle("Signal Amplitude [ADC channel]");
   timeDiffVsAmpl->GetYaxis()->SetTitle("time_{fibre}-time_{mcp}");

   timeDiffVsAmpl_fibre->GetXaxis()->SetTitle("Signal Amplitude [ADC channel]");
   timeDiffVsAmpl_fibre->GetYaxis()->SetTitle("time_{fibre}-time_{mcp}");


   timeDiffVsTime->GetXaxis()->SetTitle("time_{fibre}");
   timeDiffVsTime->GetYaxis()->SetTitle("time_{fibre}-time_{mcp}");
   if(!isNino)   timeDiffVsTime->SetAxisRange(10,40,"X");

   timeDiffVsTimeMcp->GetXaxis()->SetTitle("time_{mcp}");
   timeDiffVsTimeMcp->GetYaxis()->SetTitle("time_{fibre}-time_{mcp}");
   if(!isNino)    timeDiffVsTimeMcp->SetAxisRange(10,40,"X");

   timeDiffVsTimeMcp_amplCut->GetXaxis()->SetTitle("time_{mcp}");
   timeDiffVsTimeMcp_amplCut->GetYaxis()->SetTitle("time_{fibre}-time_{mcp}");
   if(!isNino)    timeDiffVsTimeMcp_amplCut->SetAxisRange(10,40,"X");

   timeDiffVsTime_amplCut->GetXaxis()->SetTitle("time_{fibre}");
   timeDiffVsTime_amplCut->GetYaxis()->SetTitle("time_{fibre}-time_{mcp}");
   if(!isNino)     timeDiffVsTime_amplCut->SetAxisRange(10,40,"X");

   timeDiffVsTimeMaxMinusFifty->GetXaxis()->SetTitle("time_{max}-time_{frac50}");
   timeDiffVsTimeMaxMinusFifty->GetYaxis()->SetTitle("time_{fibre}-time_{mcp}");
   timeDiffVsTimeMaxMinusFifty->SetAxisRange(-2,15,"X");

   timeDiffVsTimeMaxMinusFifty_amplCut->GetXaxis()->SetTitle("time_{max}-time_{frac50}");
   timeDiffVsTimeMaxMinusFifty_amplCut->GetYaxis()->SetTitle("time_{fibre}-time_{mcp}");
   timeDiffVsTimeMaxMinusFifty_amplCut->SetAxisRange(-2,15,"X");

   timeDiffVsTimeMaxMinusFifty_amplCut_zoom->GetXaxis()->SetTitle("time_{max}-time_{frac50}");
   timeDiffVsTimeMaxMinusFifty_amplCut_zoom->GetYaxis()->SetTitle("time_{fibre}-time_{mcp}");


   std::string constDirName = "plots_timingPerformance_";
   constDirName+=theConfiguration_.setup;
   if(theConfiguration_.addTagFileName){
     constDirName+="_";
     constDirName+=theConfiguration_.tagFileName;
   }
   system(Form("mkdir -p %s", constDirName.c_str()));
   TString dir(constDirName);


  std::string outFileName;
  outFileName = dir+"/mapCreatorOutput_"+tag+".root";
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

  sel_channel_ampl_and_time_cut_norm->Divide(sel_channel_ampl_cut_norm);
  sel_channel_ampl_cut_norm->Divide(sel_channel_norm);


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
  timing_map_sel_fibre_ampl_cut->Write();
  timing_map_sel_fibre_ampl_and_time_cut->Write();
  timing_map_sel_channel_ampl_and_time_cut->Write();
  amplitude_map_fibre2_sel_channel_ampl_cut->Write();

  timeDiffVsAmpl->Write();
  timeDiffVsAmpl_fibre->Write();
  timeDiffVsTime->Write();
  timeDiffVsTimeMcp->Write();
  timeDiffVsTimeMcp_amplCut->Write();
  timeDiffVsTime_amplCut->Write();
  timeMaxMinusFifty->Write();
  timeMaxMinusFifty_timeCut->Write();
  timeMaxMinusFifty_amplCut->Write();
  timeDiffVsTimeMaxMinusFifty->Write();
  timeDiffVsTimeMaxMinusFifty_amplCut->Write();
  timeDiffVsTimeMaxMinusFifty_amplCut_zoom->Write();
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
