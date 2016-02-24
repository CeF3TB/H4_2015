#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <sstream> 

#include "DrawTools.h"
#include "RecoTree.h"
#include "channelInfo.h"

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
#include "TProfile.h"

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
    std::cout << "./timeWalkFitter ([tag]) ([config])" << std::endl;
    exit(12345);
  }

  std::cout<<config<<std::endl;
  
  theConfiguration_=readConfiguration(config);

  float shift=10;
  bool isNino=false;

  if(theConfiguration_.setup=="Nino"){
    shift=11.5;
    isNino=true;
  }

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


  TProfile* deltaTvsAmpl=new TProfile("timeWalk","timeWalk",50,0,100,-10+shift,0+shift);
  TH2F* deltaTvsAmplHisto=new TH2F("timeWalkHisto","timeWalkHisto",50,0,100,100,-10+shift,0+shift);
  deltaTvsAmplHisto->GetXaxis()->SetTitle("Signal Amplitude [ADC Channel]");
  deltaTvsAmplHisto->GetYaxis()->SetTitle("time_{fibre}-time_{mcp} [ADC Channel]");

  TTree* recoTree=(TTree*)file->Get("recoTree");
  RecoTree t(recoTree);
  Long64_t nentries = t.fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = t.LoadTree(jentry);
      if (ientry < 0) break;
      nb = t.fChain->GetEntry(jentry);   nbytes += nb;
      if( t.mcp_time_frac50<0 ||t.mcp_max_amplitude<200)continue;
      // if(isNino && t.nino_maxAmpl<32)continue;

      float deltaTNoCorr=t.cef3_time_at_frac50->at(1)-t.mcp_time_frac50;
      if(theConfiguration_.setup=="Nino")deltaTNoCorr=t.nino_LEtime-t.mcp_time_frac50;

      if(t.nino_triggered){
	deltaTvsAmpl->Fill(t.nino_maxAmpl,deltaTNoCorr+shift);
	deltaTvsAmplHisto->Fill(t.nino_maxAmpl,deltaTNoCorr+shift);
      }

   }

   // this is the dir in which the plots will be saved:
   std::string constDirName = "timeWalkFiles_";
   constDirName+=theConfiguration_.setup;
   if(theConfiguration_.addTagFileName){
     constDirName+="_";
     constDirName+=theConfiguration_.tagFileName;
   }
   system(Form("mkdir -p %s", constDirName.c_str()));
   TString dir(constDirName);


   std::string outFileName;
   outFileName = dir+"/timeWalkFile_"+tag+".root";
   TFile* outFile = TFile::Open(outFileName.c_str(),"recreate");

   TF1 f("fit_log","[0]*log([1]*x)",0.,100.);
   f.SetParameters(-1.,20.);
   deltaTvsAmpl->Fit("fit_log","R","",14,70);


   TCanvas c1;
   c1.cd();
   deltaTvsAmplHisto->Draw("colz");
   deltaTvsAmpl->Draw("same");
   c1.SaveAs(dir+"/deltaTvsAmpl.png");
   c1.SaveAs(dir+"/deltaTvsAmpl.pdf");

   f.Write();

   deltaTvsAmpl->Write();
   deltaTvsAmplHisto->Write();
   outFile->Write();
   outFile->Close();


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
