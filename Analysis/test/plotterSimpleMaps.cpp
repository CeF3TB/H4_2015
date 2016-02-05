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
#include "TKey.h"

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

  
  
  TString filename= "plots_timingPerformance_"+theConfiguration_.setup+"/mapCreatorOutput_"+tag+".root";
  
  TFile* file = TFile::Open(filename.Data());

  std::map<TString,TH1F*> histoNames;
  std::map<TString,TH2F*> histoNames2D;
  

  TList* list = file->GetListOfKeys() ;
  if (!list) { printf("<E> No keys found in file\n") ; exit(1) ; }
  TIter next(list) ;
  TKey* key ;
  TObject* obj ;

  
  while ( key = (TKey*)next() ) {
    obj = key->ReadObj() ;
    if (    (strcmp(obj->IsA()->GetName(),"TProfile")!=0)
	    && (!obj->InheritsFrom("TH2"))
	    && (!obj->InheritsFrom("TH1")) 
	    ) {
      printf("<W> Object %s is not 1D or 2D histogram : "
	     "will not be converted\n",obj->GetName()) ;
    }
    printf("Histo name:%s title:%s\n",obj->GetName(),obj->GetTitle());
    if((strcmp(obj->IsA()->GetName(),"TH1F"))==0){
      histoNames[obj->GetName()]=(TH1F*)obj;
    }else if ((strcmp(obj->IsA()->GetName(),"TH2F"))==0){
      histoNames2D[obj->GetName()]=(TH2F*)obj;
    }
  }



  for(std::map<TString,TH2F*>::const_iterator out=histoNames2D.begin();out!=histoNames2D.end();++out){

    gStyle->SetPadRightMargin(0.19);//for the palette and ztitle
    TCanvas c1;
    TH2F* histo=(TH2F*)out->second->Clone(out->first+"_clone");
    histo->GetZaxis()->SetTitleOffset(1.4);
    if(out->first.Contains("amplitude"))histo->SetZTitle("Amplitude [ADC]");
    if(out->first.Contains("timing"))histo->SetZTitle("#DeltaT [ns]");
    histo->Draw("colz");
    c1.SaveAs(dir+"/"+out->first+".pdf");
    c1.SaveAs(dir+"/"+out->first+".png");
  }





    gStyle->SetPadRightMargin(0.10);
    TH1F* maxAmpl_sel_fibre=(TH1F*)file->Get("maxAmpl_sel_fibre");
    TH1F* maxAmpl_sel_channel=(TH1F*)file->Get("maxAmpl_sel_channel");

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
