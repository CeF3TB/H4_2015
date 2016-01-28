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

  std::string constDirName = "plots_maps_";
  constDirName+=theConfiguration_.setup;
  if(theConfiguration_.addTagFileName){
    constDirName+="_";
    constDirName+=theConfiguration_.tagFileName;
  }
  system(Form("mkdir -p %s", constDirName.c_str()));

  
  TString dir(constDirName);



  std::vector<int> runs;
  std::vector<float> energies;
  for (int i=0;i<theConfiguration_.runs.size();++i){
    runs.push_back(theConfiguration_.runs[i]);
    energies.push_back(theConfiguration_.energies[i]);

    TString filename= "plots_timing_SiPM/timingStudies"+theConfiguration_.setup+"_"+tag+"_";
    filename+=runs[i];
    filename+=".root";

    TFile* file = TFile::Open(filename.Data());
    
    TH2F* timing_map_total=(TH2F*)file->Get("timing_map_channel");
    TH2F* amplitude_map_total=(TH2F*)file->Get("amplitude_map_channel");
    TH2F* amplitude_map_fibre2=(TH2F*)file->Get("amplitude_map_fibre2");

    TH2F* timing_map_sel_channel=(TH2F*)file->Get("timing_map_sel_channel");
    TH2F* amplitude_map_sel_channel=(TH2F*)file->Get("amplitude_map_sel_channel");

    TH2F* timing_map_sel_fibre=(TH2F*)file->Get("timing_map_sel_fibre");
    TH2F* amplitude_map_sel_fibre=(TH2F*)file->Get("amplitude_map_sel_fibre");

    TH1F* maxAmpl_sel_fibre=(TH1F*)file->Get("maxAmpl_sel_fibre");
    TH1F* maxAmpl_sel_channel=(TH1F*)file->Get("maxAmpl_sel_channel");


    TString runString;
    runString+=runs[i];
    TCanvas c1;
    amplitude_map_total->Draw("colz");
    c1.SaveAs(dir+"/amplitude_map_"+runString+".pdf");

  }

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
