#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <sstream> 

#include "DrawTools.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TGraphErrors.h"
#include "TVectorD.h"
#include "TF1.h"
#include "interface/Configurator.h"
#include "timingResolutionPlotsConfigurator.h"

timingResolutionPlots_Config_t readConfiguration(std::string configName);

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
    std::cout << "./timingStudies ([tag]) ([config])" << std::endl;
    exit(12345);
  }

  std::cout<<config<<std::endl;
  
  theConfiguration_=readConfiguration(config);
  


  std::vector<float> energies;
  std::vector<int> wp_channel;
  std::vector<int> wp_fibre;
  std::vector<int> runs;

  for (int i=0;i<theConfiguration_.energies.size();++i){
    energies.push_back(theConfiguration_.energies[i]);
    wp_channel.push_back(theConfiguration_.wp_channel[i]);
    wp_fibre.push_back(theConfiguration_.wp_fibre[i]);
    runs.push_back(theConfiguration_.runs[i]);
  }

  std::string constDirName = "plots_timing_";
  constDirName+=theConfiguration_.setup;
  system(Form("mkdir -p %s", constDirName.c_str()));
  TString dir(constDirName);

  TGraphErrors* resVsEnergy_channel = new TGraphErrors(0);
  TGraphErrors* resVsEnergy_fibre = new TGraphErrors(0);

  std::vector<TGraphErrors*> resVsAmplitude;


  for (int i=0;i<runs.size();++i){
    TString run;
    run.Form("%d",runs[i]); 
    resVsAmplitude.push_back(new TGraphErrors(0));


    TFile* file = TFile::Open(dir+"/timingStudiesSiPM_"+tag+"_"+run+".root");

    if( file==0 ) {
      std::cout << "ERROR! Din't find file for run" << runs[i] << std::endl;
      std::cout << "Exiting." << std::endl;
      exit(11);
    }
    
    //channel
    TVectorD* resValueTime_channel=(TVectorD*)file->Get("resValueTime_channel");
    TVectorD* resErrValueTime_channel=(TVectorD*)file->Get("resErrValueTime_channel");
    TVectorD* resValueAmplitude_channel=(TVectorD*)file->Get("resValueAmplitude_channel");
    TVectorD* resErrValueAmplitude_channel=(TVectorD*)file->Get("resErrValueAmplitude_channel");
    TVectorD* resErrRelativeValueTime_channel=(TVectorD*)file->Get("resErrRelativeValueTime_channel");
    TVectorD* nEntriesTime_channel=(TVectorD*)file->Get("nEntriesTime_channel");

    //fibre
    TVectorD* resValueTime_fibre=(TVectorD*)file->Get("resValueTime_fibre");
    TVectorD* resErrValueTime_fibre=(TVectorD*)file->Get("resErrValueTime_fibre");

    
    //    resVsEnergy->SetPoint( i, energies[i],(*resValueTime)[wp_100[i]] );
    //FIXME    resVsEnergy->SetPointError( i, 0, (*resErrValueTime)[wp_100[i]]);


    resVsEnergy_channel->SetPoint( i, energies[i], (*resValueTime_channel)[wp_channel[i]]);
    resVsEnergy_channel->SetPointError( i, 0,(*resErrValueTime_channel)[wp_channel[i]] );

    resVsEnergy_fibre->SetPoint( i, energies[i], (*resValueTime_fibre)[wp_fibre[i]]);
    resVsEnergy_fibre->SetPointError( i, 0,(*resErrValueTime_fibre)[wp_fibre[i]] );
 

    for (int j=0;j<resValueTime_channel->GetNoElements();j++){
      //      std::cout<<j<<" "<<(*resValueAmplitude_channel)[j]<<std::endl;
      if ((*resValueAmplitude_channel)[j]<10)continue;
      resVsAmplitude[i]->SetPoint(j,(20+j*10),(*resValueAmplitude_channel)[j]);
      resVsAmplitude[i]->SetPointError(j,0,(*resErrValueAmplitude_channel)[j]);
    }

    
  }
  

  TFile * outFile = new TFile(dir+"/plotsTimingStudiesSiPM_"+tag+".root","recreate");

  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );

  float yup=1.1*(resVsAmplitude[2]->GetY())[0];
  if((resVsAmplitude[2]->GetY())[0]<10) yup = 1.1*(resVsAmplitude[1]->GetY())[0];
  float  xlow=(resVsAmplitude[2]->GetX())[0]-10;
  if(xlow<10)xlow = (resVsAmplitude[1]->GetX())[0]-10;
  float xup=(resVsAmplitude[2]->GetX())[resVsAmplitude[2]->GetN()-1]+10; 
  if(xup<10)xup = (resVsAmplitude[2]->GetX())[resVsAmplitude[2]->GetN()-1]+10; 
  TH2D* h2_axes_2 = new TH2D( "axes_2", "", 100, xlow,xup , 110, 0., yup);

  //  std::cout<<(resVsEnergy_channel->GetY())[0]+(resVsEnergy_channel->GetErrorY(0))<<std::endl;

  TH2D* h2_axes_3 = new TH2D( "axes_1", "", 100, -0.0, 250. , 110, 60., 1.1*((resVsEnergy_channel->GetY())[0]+(resVsEnergy_channel->GetErrorY(0))));  
  h2_axes_3->SetXTitle("Beam Energy [GeV]");
  h2_axes_3->SetYTitle("time_{Fibre}-time_{mcp} [ps]");

  TH2D* h2_axes_4 = new TH2D( "axes_1", "", 100, -0.0, 250. , 110, 60., 1.1*((resVsEnergy_fibre->GetY())[0]+(resVsEnergy_fibre->GetErrorY(0))));  
  h2_axes_4->SetXTitle("Beam Energy [GeV]");
  h2_axes_4->SetYTitle("time_{Fibre}-time_{mcp} [ps]");

  
  //resolution vs energy
  //channel
  TF1* f_ene= new TF1("fun_ene","[1]/x+[0]",(resVsEnergy_channel->GetX())[0]-10,(resVsEnergy_channel->GetX())[resVsEnergy_channel->GetN()-1] +10);

  h2_axes_3->Draw("");
  resVsEnergy_channel->Fit("fun_ene");

  TLegend* lego_2 = new TLegend(0.47, 0.7, 0.8, 0.92);
  lego_2->SetTextSize(0.038);
  lego_2->AddEntry(  (TObject*)0 ,"f(x) = p0 + p1/x", "");
  lego_2->AddEntry(  (TObject*)0 ,Form("p0 = %.0f #pm %.0f", f_ene->GetParameter(0), f_ene->GetParError(0) ), "");
  lego_2->AddEntry(  (TObject*)0 ,Form("p1 = %.0f #pm %.0f", f_ene->GetParameter(1), f_ene->GetParError(1) ), "");
  lego_2->SetFillColor(0);
  lego_2->Draw("same");


  resVsEnergy_channel->SetMarkerStyle(20);
  resVsEnergy_channel->SetMarkerSize(1.6);
  resVsEnergy_channel->SetMarkerColor(kBlue);
  resVsEnergy_channel->Draw("p same");
  resVsEnergy_channel->SetName("resVsEnergy_channel");
  resVsEnergy_channel->Write();

  c1->SaveAs(dir+"/timingResolutionVsEnergySiPM_channel.png");
  c1->SaveAs(dir+"/timingResolutionVsEnergySiPM_channel.pdf");  
  //fibre
  h2_axes_4->Draw("");
  resVsEnergy_fibre->Fit("fun_ene");

  lego_2->Clear();
  lego_2->SetTextSize(0.038);
  lego_2->AddEntry(  (TObject*)0 ,"f(x) = p0 + p1/x", "");
  lego_2->AddEntry(  (TObject*)0 ,Form("p0 = %.0f #pm %.0f", f_ene->GetParameter(0), f_ene->GetParError(0) ), "");
  lego_2->AddEntry(  (TObject*)0 ,Form("p1 = %.0f #pm %.0f", f_ene->GetParameter(1), f_ene->GetParError(1) ), "");
  lego_2->SetFillColor(0);
  lego_2->Draw("same");


  resVsEnergy_fibre->SetMarkerStyle(20);
  resVsEnergy_fibre->SetMarkerSize(1.6);
  resVsEnergy_fibre->SetMarkerColor(kBlue);
  resVsEnergy_fibre->Draw("p same");
  resVsEnergy_fibre->SetName("resVsEnergy_fibre");
  resVsEnergy_fibre->Write();

  c1->SaveAs(dir+"/timingResolutionVsEnergySiPM_fibre.png");
  c1->SaveAs(dir+"/timingResolutionVsEnergySiPM_fibre.pdf");  

  for(int i=0;i<runs.size();++i){ 
    TString ene;
    ene.Form("%.0f",energies[i]); 
    std::string energy(Form("%.0f", energies[i]));

    h2_axes_2->SetYTitle("#sigma_{t} [ps]");
    h2_axes_2->GetXaxis()->SetTitle("Amplitude [ADC]");
    h2_axes_2->Draw(""); 

    TF1* f= new TF1("fun","[1]/x+[0]",(resVsAmplitude[i]->GetX())[0]-10,(resVsAmplitude[i]->GetX())[resVsAmplitude[i]->GetN()-1] +10);
    resVsAmplitude[i]->Fit("fun");

    resVsAmplitude[i]->SetName("resVsAmplitude_"+ene);
    resVsAmplitude[i]->SetMarkerStyle(20);
    resVsAmplitude[i]->SetMarkerSize(1.6);
    resVsAmplitude[i]->SetMarkerColor(kBlue);
    resVsAmplitude[i]->Draw("p same");


    TPaveText* pave = DrawTools::getLabelTop_expOnXaxis(energy+" GeV Electron Beam");
    pave->Draw("same");



    TLegend* lego = new TLegend(0.47, 0.7, 0.8, 0.92);
    lego->SetTextSize(0.038);
    lego->AddEntry(  (TObject*)0 ,"f(x) = p0 + p1/x", "");
    lego->AddEntry(  (TObject*)0 ,Form("p0 = %.0f #pm %.0f", f->GetParameter(0), f->GetParError(0) ), "");
    lego->AddEntry(  (TObject*)0 ,Form("p1 = %.0f #pm %.0f", f->GetParameter(1), f->GetParError(1) ), "");
    lego->SetFillColor(0);
    lego->Draw("same");

    c1->SaveAs(dir+"/timingResolutionVsAmplitudeSiPM_"+ene+"GeV.png");
    c1->SaveAs(dir+"/timingResolutionVsAmplitudeSiPM_"+ene+"GeV.pdf");  

 

    resVsAmplitude[i]->Write();
  }

  outFile->Write();
  outFile->Close();

  return 0;

}


timingResolutionPlots_Config_t readConfiguration(std::string configName){

  Configurator* configurator_ = new Configurator();
  std::string fileName = Form ("./config_timing/%s.xml",configName.c_str());
  configurator_->xmlFileName=fileName.c_str();
  configurator_->Init();

  timingResolutionPlots_Config_t conf;
 
  conf.setup=Configurable::getElementContent (*configurator_, "setup",configurator_->root_element) ;

  string content = Configurable::getElementContent (*configurator_, "energies",configurator_->root_element) ;
  Configurator::GetVecFloat(content,conf.energies);

  content = Configurable::getElementContent (*configurator_, "wp_channel",configurator_->root_element) ;
  Configurator::GetVecInt(content,conf.wp_channel);

  content = Configurable::getElementContent (*configurator_, "wp_fibre",configurator_->root_element) ;
  Configurator::GetVecInt(content,conf.wp_fibre);
  
  content = Configurable::getElementContent (*configurator_, "runs",configurator_->root_element) ;
  Configurator::GetVecInt(content,conf.runs);

   return conf;

}
