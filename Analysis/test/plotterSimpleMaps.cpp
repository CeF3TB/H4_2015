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
#include "TGraphErrors.h"

#include "interface/Configurator.h"
#include "plotOverallPerformancesConfigurator.h"

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
    std::cout << "./plotterMaps ([tag]) ([config])" << std::endl;
    exit(12345);
  }

  std::cout<<config<<std::endl;
  
  theConfiguration_=readConfiguration(config);

  std::string constDirName = "plots_Simplemaps_";
  constDirName+=theConfiguration_.setup;
  if(theConfiguration_.addTagFileName){
    constDirName+="_";
    constDirName+=theConfiguration_.tagFileName;
  }
  system(Form("mkdir -p %s", constDirName.c_str()));

  
  TString dir(constDirName);



  std::vector<int> runs;
  std::vector<float> energies;


    TString filename= "plots_timingPerformance_"+theConfiguration_.setup+"/mapCreatorOutput_"+tag+".root";

    TFile* file = TFile::Open(filename.Data());
    
    TH2F* timing_map_total=(TH2F*)file->Get("timing_map_channel");
    TH2F* amplitude_map_total=(TH2F*)file->Get("amplitude_map_channel");
    TH2F* amplitude_map_fibre2=(TH2F*)file->Get("amplitude_map_fibre2");

    TH2F* timing_map_sel_channel=(TH2F*)file->Get("timing_map_sel_channel");
    TH2F* timing_map_sel_channel_ampl_cut=(TH2F*)file->Get("timing_map_sel_channel_ampl_cut");
    TH2F* timing_map_sel_fibre_ampl_cut=(TH2F*)file->Get("timing_map_sel_fibre_ampl_cut");
    TH2F* timing_map_sel_fibre_ampl_and_time_cut=(TH2F*)file->Get("timing_map_sel_fibre_ampl_and_time_cut");
    TH2F* timing_map_sel_channel_ampl_and_time_cut=(TH2F*)file->Get("timing_map_sel_channel_ampl_and_time_cut");
    TH2F* amplitude_map_sel_channel=(TH2F*)file->Get("amplitude_map_sel_channel");
    TH2F* amplitude_map_sel_channel_ampl_cut=(TH2F*)file->Get("amplitude_map_sel_channel_ampl_cut");
    TH2F* amplitude_map_sel_channel_ampl_and_time_cut=(TH2F*)file->Get("amplitude_map_sel_channel_ampl_and_time_cut");

    TH2F* timeDiffVsAmpl=(TH2F*)file->Get("timeDiffVsAmpl");
    TH2F* timeDiffVsAmpl_fibre=(TH2F*)file->Get("timeDiffVsAmpl_fibre");
    TH2F* timeDiffVsTime=(TH2F*)file->Get("timeDiffVsTime");
    TH2F* timeDiffVsTimeMcp=(TH2F*)file->Get("timeDiffVsTimeMcp");
    TH2F* timeDiffVsTime_amplCut=(TH2F*)file->Get("timeDiffVsTime_amplCut");
    TH2F* timeDiffVsTimeMcp_amplCut=(TH2F*)file->Get("timeDiffVsTimeMcp_amplCut");
    TH2F* timeDiffVsTimeMaxMinusFifty=(TH2F*)file->Get("timeDiffVsTimeMaxMinusFifty");
    TH2F* timeDiffVsTimeMaxMinusFifty_amplCut=(TH2F*)file->Get("timeDiffVsTimeMaxMinusFifty_amplCut");
    TH2F* timeDiffVsTimeMaxMinusFifty_amplCut_zoom=(TH2F*)file->Get("timeDiffVsTimeMaxMinusFifty_amplCut_zoom");

    TH2F* timing_map_sel_fibre=(TH2F*)file->Get("timing_map_sel_fibre");
    TH2F* amplitude_map_sel_fibre=(TH2F*)file->Get("amplitude_map_sel_fibre");

    TH2F* sel_channel_ampl_cut_norm=(TH2F*)file->Get("sel_channel_ampl_cut_norm");
    TH2F* sel_channel_ampl_and_time_cut_norm=(TH2F*)file->Get("sel_channel_ampl_and_time_cut_norm");

    TH1F* maxAmpl_sel_fibre=(TH1F*)file->Get("maxAmpl_sel_fibre");
    TH1F* maxAmpl_sel_channel=(TH1F*)file->Get("maxAmpl_sel_channel");



    gStyle->SetPadRightMargin(0.17);//for the palette
    TCanvas c1;
    amplitude_map_total->Draw("colz");
    c1.SaveAs(dir+"/amplitude_map.pdf");
    c1.SaveAs(dir+"/amplitude_map.png");

    c1.Clear();
    amplitude_map_fibre2->Draw("colz");
    c1.SaveAs(dir+"/amplitude_map_fibre2.pdf");
    c1.SaveAs(dir+"/amplitude_map_fibre2.png");


    c1.Clear();
    timing_map_total->Draw("colz");
    c1.SaveAs(dir+"/timing_map.pdf");
    c1.SaveAs(dir+"/timing_map.png");

    c1.Clear();
    amplitude_map_sel_channel->Draw("colz");
    c1.SaveAs(dir+"/amplitude_map_sel_channel.pdf");
    c1.SaveAs(dir+"/amplitude_map_sel_channel.png");

    c1.Clear();
    amplitude_map_sel_channel_ampl_cut->SetAxisRange(100,250,"Z");
    amplitude_map_sel_channel_ampl_cut->Draw("colz");
    c1.SaveAs(dir+"/amplitude_map_sel_channel_ampl_cut.pdf");
    c1.SaveAs(dir+"/amplitude_map_sel_channel_ampl_cut.png");


    c1.Clear();
    timing_map_sel_channel_ampl_cut->SetAxisRange(3.5,5.5,"Z");
    timing_map_sel_channel_ampl_cut->Draw("colz");
    c1.SaveAs(dir+"/timing_map_sel_channel_ampl_cut.pdf");
    c1.SaveAs(dir+"/timing_map_sel_channel_ampl_cut.png");

    c1.Clear();
    timing_map_sel_fibre_ampl_cut->SetAxisRange(3.5,5,"Z");
    timing_map_sel_fibre_ampl_cut->Draw("colz");
    c1.SaveAs(dir+"/timing_map_sel_fibre_ampl_cut.pdf");
    c1.SaveAs(dir+"/timing_map_sel_fibre_ampl_cut.png");

    c1.Clear();
    timing_map_sel_fibre_ampl_and_time_cut->SetAxisRange(3.5,5,"Z");
    timing_map_sel_fibre_ampl_and_time_cut->Draw("colz");
    c1.SaveAs(dir+"/timing_map_sel_fibre_ampl_and_time_cut.pdf");
    c1.SaveAs(dir+"/timing_map_sel_fibre_ampl_and_time_cut.png");



    c1.Clear();
    timing_map_sel_channel_ampl_and_time_cut->SetAxisRange(3.5,5.5,"Z");
    timing_map_sel_channel_ampl_and_time_cut->Draw("colz");
    c1.SaveAs(dir+"/timing_map_sel_channel_ampl_and_time_cut.pdf");
    c1.SaveAs(dir+"/timing_map_sel_channel_ampl_and_time_cut.png");


    c1.Clear();
    timing_map_sel_channel->SetAxisRange(3.5,5.5,"Z");
    timing_map_sel_channel->Draw("colz");
    c1.SaveAs(dir+"/timing_map_sel_channel.pdf");
    c1.SaveAs(dir+"/timing_map_sel_channel.png");


    c1.Clear();
    amplitude_map_sel_fibre->Draw("colz");
    c1.SaveAs(dir+"/amplitude_map_sel_fibre.pdf");
    c1.SaveAs(dir+"/amplitude_map_sel_fibre.png");

    c1.Clear();
    sel_channel_ampl_cut_norm->Draw("colz");
    c1.SaveAs(dir+"/sel_channel_ampl_cut_norm.pdf");
    c1.SaveAs(dir+"/sel_channel_ampl_cut_norm.png");

    c1.Clear();
    sel_channel_ampl_and_time_cut_norm->Draw("colz");
    c1.SaveAs(dir+"/sel_channel_ampl_and_time_cut_norm.pdf");
    c1.SaveAs(dir+"/sel_channel_ampl_and_time_cut_norm.png");



    c1.Clear();
    timing_map_sel_fibre->Draw("colz");
    c1.SaveAs(dir+"/timing_map_sel_fibre.pdf");
    c1.SaveAs(dir+"/timing_map_sel_fibre.png");

    c1.Clear();
    timeDiffVsAmpl->Draw("colz");
    c1.SaveAs(dir+"/timeDiffVsAmpl.pdf");
    c1.SaveAs(dir+"/timeDiffVsAmpl.png");

    c1.Clear();
    timeDiffVsAmpl_fibre->Draw("colz");
    c1.SaveAs(dir+"/timeDiffVsAmpl_fibre.pdf");
    c1.SaveAs(dir+"/timeDiffVsAmpl_fibre.png");


    c1.Clear();
    timeDiffVsTime->SetAxisRange(10,40,"X");
    timeDiffVsTime->Draw("colz");
    c1.SaveAs(dir+"/timeDiffVsTime.pdf");
    c1.SaveAs(dir+"/timeDiffVsTime.png");

    c1.Clear();
    timeDiffVsTime_amplCut->SetAxisRange(10,40,"X");
    timeDiffVsTime_amplCut->Draw("colz");
    c1.SaveAs(dir+"/timeDiffVsTime_amplCut.pdf");
    c1.SaveAs(dir+"/timeDiffVsTime_amplCut.png");


    c1.Clear();
    timeDiffVsTimeMcp->SetAxisRange(10,40,"X");
    timeDiffVsTimeMcp->Draw("colz");
    c1.SaveAs(dir+"/timeDiffVsTimeMcp.pdf");
    c1.SaveAs(dir+"/timeDiffVsTimeMcp.png");

    c1.Clear();
    timeDiffVsTimeMcp_amplCut->SetAxisRange(10,40,"X");
    timeDiffVsTimeMcp_amplCut->Draw("colz");
    c1.SaveAs(dir+"/timeDiffVsTimeMcp_amplCut.pdf");
    c1.SaveAs(dir+"/timeDiffVsTimeMcp_amplCut.png");



    c1.Clear();
    timeDiffVsTimeMaxMinusFifty->SetAxisRange(-2,15,"X");
    timeDiffVsTimeMaxMinusFifty->Draw("colz");
    c1.SaveAs(dir+"/timeDiffVsTimeMaxMinusFifty.pdf");
    c1.SaveAs(dir+"/timeDiffVsTimeMaxMinusFifty.png");


    c1.Clear();
    timeDiffVsTimeMaxMinusFifty_amplCut->SetAxisRange(-2,15,"X");
    timeDiffVsTimeMaxMinusFifty_amplCut->Draw("colz");
    c1.SaveAs(dir+"/timeDiffVsTimeMaxMinusFifty_amplCut.pdf");
    c1.SaveAs(dir+"/timeDiffVsTimeMaxMinusFifty_amplCut.png");

    c1.Clear();
    timeDiffVsTimeMaxMinusFifty_amplCut_zoom->Draw("colz");
    c1.SaveAs(dir+"/timeDiffVsTimeMaxMinusFifty_amplCut_zoom.pdf");
    c1.SaveAs(dir+"/timeDiffVsTimeMaxMinusFifty_amplCut_zoom.png");



    gStyle->SetPadRightMargin(0.10);
    TCanvas c2;
    maxAmpl_sel_channel->GetXaxis()->SetRangeUser(0,800.);
    maxAmpl_sel_channel->SetLineWidth(2);
    maxAmpl_sel_fibre->SetLineWidth(2);
    maxAmpl_sel_channel->DrawNormalized();
    maxAmpl_sel_fibre->DrawNormalized("same");



    TPaveText* pave = DrawTools::getLabelTop_expOnXaxis("Electron Beam");
    pave->Draw("same");
    
    c2.SaveAs(dir+"/maxAmpl_comparison.pdf");
    c2.SaveAs(dir+"/maxAmpl_comparison.png");



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

   conf.addTagFileName= Configurator::GetInt(Configurable::getElementContent(*configurator_,"addTagFileName",configurator_->root_element));

   conf.tagFileName=Configurable::getElementContent (*configurator_, "tagFileName",configurator_->root_element) ;

  content = Configurable::getElementContent (*configurator_, "energies",configurator_->root_element) ;
  Configurator::GetVecFloat(content,conf.energies);


   return conf;

}
