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

  shift=0;//FIXME
  TProfile* deltaTvsAmpl=new TProfile("timeWalk","timeWalk",50,0,100,-10+shift,0+shift);
  TH2F* deltaTvsAmplHisto=new TH2F("timeWalkHisto","timeWalkHisto",50,0,100,100,-10+shift,0+shift);

  TH2F* deltaTvsAmplHisto_channel=new TH2F("timeWalkHisto_channel","timeWalkHisto_channel",50,0,100,100,-10+shift,0+shift);
  TH2F* deltaTvsAmplHisto_fibre=new TH2F("timeWalkHisto_fibre","timeWalkHisto_fibre",50,0,100,100,-10+shift,0+shift);


  TProfile* deltaTvsAmpl_channel=new TProfile("timeWalk_channel","timeWalk_channel",50,0,100,-10+shift,0+shift);
  TProfile* deltaTvsAmpl_fibre=new TProfile("timeWalk_fibre","timeWalk_fibre",50,0,100,-10+shift,0+shift);


  std::vector<TH1F*> deltaTBinAmpl_channel;
  std::vector<TH1F*> deltaTBinAmpl_fibre;

  float amplCut=10;
  int shiftBin=0;



  for (int i=1;i<deltaTvsAmplHisto->GetNbinsX();++i){//i=0 is underflow
    TString icut;
    icut.Form("%d",i); 
    if(deltaTvsAmplHisto->GetXaxis()->GetBinLowEdge(i)>amplCut){
      deltaTBinAmpl_channel.push_back(new TH1F("deltaTBinAmpl_channel"+icut,"deltaTBinAmpl_channel"+icut,200,-10,0));
      deltaTBinAmpl_fibre.push_back(new TH1F("deltaTBinAmpl_fibre"+icut,"deltaTBinAmpl_fibre"+icut,200,-10,0));
    }else{
      shiftBin++;
    }
  }




  

  deltaTvsAmplHisto->GetXaxis()->SetTitle("Signal Amplitude [ADC Channel]");
  deltaTvsAmplHisto->GetYaxis()->SetTitle("time_{fibre}-time_{mcp} [ADC Channel]");

  deltaTvsAmplHisto_channel->GetXaxis()->SetTitle("Signal Amplitude [ADC Channel]");
  deltaTvsAmplHisto_channel->GetYaxis()->SetTitle("time_{fibre}-time_{mcp} [ADC Channel]");

  deltaTvsAmplHisto_fibre->GetXaxis()->SetTitle("Signal Amplitude [ADC Channel]");
  deltaTvsAmplHisto_fibre->GetYaxis()->SetTitle("time_{fibre}-time_{mcp} [ADC Channel]");


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

      if(t.nino_triggered && t.cef3_maxAmpl->at(1)>700){
	deltaTvsAmpl->Fill(t.nino_maxAmpl,deltaTNoCorr+shift);
	deltaTvsAmplHisto->Fill(t.nino_maxAmpl,deltaTNoCorr+shift);
	
	for(int i=1;i<deltaTBinAmpl_channel.size();++i) {//i=0 is underflow
	   
	      if(TopologicalSelectionHelper::passesChannelTopologicalSelection(t,isNino)){
		deltaTvsAmplHisto_channel->Fill(t.nino_maxAmpl,deltaTNoCorr+shift);
		if(t.nino_maxAmpl<amplCut)break;
		if(t.nino_maxAmpl>deltaTvsAmplHisto->GetXaxis()->GetBinLowEdge(i) && t.nino_maxAmpl<deltaTvsAmplHisto->GetXaxis()->GetBinUpEdge(i)) {
		  deltaTBinAmpl_channel[i-shiftBin]->Fill(deltaTNoCorr);
		  break;
		}
	      }else if(TopologicalSelectionHelper::passesFibreTopologicalSelection(t,isNino)){
		deltaTvsAmplHisto_fibre->Fill(t.nino_maxAmpl,deltaTNoCorr+shift);
		if(t.nino_maxAmpl<amplCut)break;
		if(t.nino_maxAmpl>deltaTvsAmplHisto->GetXaxis()->GetBinLowEdge(i) && t.nino_maxAmpl<deltaTvsAmplHisto->GetXaxis()->GetBinUpEdge(i)) {
		  deltaTBinAmpl_fibre[i-shiftBin]->Fill(deltaTNoCorr);
		  break;
		}
	
	      }


	   
	  }
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



   TVectorD mean_channel(deltaTBinAmpl_channel.size());
   TVectorD sigma_channel(deltaTBinAmpl_channel.size());
   TVectorD ampl_channel(deltaTBinAmpl_channel.size());
   TVectorD amplErr_channel(deltaTBinAmpl_channel.size());

   TVectorD mean_fibre(deltaTBinAmpl_fibre.size());
   TVectorD sigma_fibre(deltaTBinAmpl_fibre.size());
   TVectorD ampl_fibre(deltaTBinAmpl_fibre.size());
   TVectorD amplErr_fibre(deltaTBinAmpl_fibre.size());


   for(int i=0;i<deltaTBinAmpl_channel.size();++i) {
     


     TString icut;
     icut.Form("%d",i); 
     //channel
     if(deltaTBinAmpl_channel[i]->GetEntries()>25){     
       if(i==0){
	 ampl_channel[i]=amplCut+deltaTvsAmplHisto->GetXaxis()->GetBinWidth(0)/2.;
       }else{
	 ampl_channel[i]=ampl_channel[i-1]+deltaTvsAmplHisto->GetXaxis()->GetBinWidth(0);
       }
       
       amplErr_channel[i]=0;
       
       TCanvas c1;
       c1.cd();
       TF1 f_fit("f_fit","gaus",deltaTBinAmpl_channel[i]->GetMean()-deltaTBinAmpl_channel[i]->GetRMS()*2,deltaTBinAmpl_channel[i]->GetMean()+deltaTBinAmpl_channel[i]->GetRMS()/2.);
       deltaTBinAmpl_channel[i]->Fit("f_fit","R","",deltaTBinAmpl_channel[i]->GetMean()-deltaTBinAmpl_channel[i]->GetRMS()*2,deltaTBinAmpl_channel[i]->GetMean()+deltaTBinAmpl_channel[i]->GetRMS()/2.);
       deltaTBinAmpl_channel[i]->Draw("");
       mean_channel[i]=f_fit.GetParameter(1);
       sigma_channel[i]=f_fit.GetParError(1);
       
       c1.SaveAs(dir+"/deltaTBinAmpl_channel"+icut+".png");
       c1.SaveAs(dir+"/deltaTBinAmpl_channel"+icut+".pdf");
     }

     //fibre
     if(deltaTBinAmpl_fibre[i]->GetEntries()>25){     
       if(i==0){
	 ampl_fibre[i]=amplCut+deltaTvsAmplHisto->GetXaxis()->GetBinWidth(0)/2.;
       }else{
	 ampl_fibre[i]=ampl_fibre[i-1]+deltaTvsAmplHisto->GetXaxis()->GetBinWidth(0);
       }

       amplErr_fibre[i]=0;

       TCanvas c1;
       c1.cd();
       TF1 f_fit_fibre("f_fit_fibre","gaus",deltaTBinAmpl_fibre[i]->GetMean()-deltaTBinAmpl_fibre[i]->GetRMS()*2,deltaTBinAmpl_fibre[i]->GetMean()+deltaTBinAmpl_fibre[i]->GetRMS()/2.);
       deltaTBinAmpl_fibre[i]->Fit("f_fit_fibre","R","",deltaTBinAmpl_fibre[i]->GetMean()-deltaTBinAmpl_fibre[i]->GetRMS()*2,deltaTBinAmpl_fibre[i]->GetMean()+deltaTBinAmpl_fibre[i]->GetRMS()/2.);
       deltaTBinAmpl_fibre[i]->Draw("");
       mean_fibre[i]=f_fit_fibre.GetParameter(1);
       sigma_fibre[i]=f_fit_fibre.GetParError(1);

       c1.SaveAs(dir+"/deltaTBinAmpl_fibre"+icut+".png");
       c1.SaveAs(dir+"/deltaTBinAmpl_fibre"+icut+".pdf");
     }
     
   }
   

   TGraphErrors* graph_fit_channel=new TGraphErrors(ampl_channel,mean_channel,amplErr_channel,sigma_channel);
   graph_fit_channel->SetMarkerColor(kGreen+2);
   graph_fit_channel->SetLineColor(kGreen+2);


   TGraphErrors* graph_fit_fibre=new TGraphErrors(ampl_fibre,mean_fibre,amplErr_fibre,sigma_fibre);
   graph_fit_fibre->SetMarkerColor(kViolet);
   graph_fit_fibre->SetLineColor(kViolet);


   TF1 f_channel("fit_log_channel","[0]*log([1]*x)",0.,100.);
   //   TF1 f("fit_log","[0]+[1]/x",14,100);
   f_channel.SetParameters(-1.,20.);
   f_channel.SetLineColor(kGreen+2);
   //   deltaTvsAmpl->Fit("fit_log","R","",14,70);
   //   deltaTvsAmpl->Fit("fit_log","R","",20,100);
   graph_fit_channel->Fit("fit_log_channel","R","",14,60);


   TF1 f_fibre("fit_log_fibre","[0]*log([1]*x)",0.,100.);
   //   TF1 f("fit_log","[0]+[1]/x",14,100);
   f_fibre.SetParameters(-1.,20.);
   f_fibre.SetLineColor(kViolet);  
   //   deltaTvsAmpl->Fit("fit_log","R","",14,70);
   //   deltaTvsAmpl->Fit("fit_log","R","",20,100);
   graph_fit_fibre->Fit("fit_log_fibre","R","",14,60);


   gStyle->SetPadRightMargin(0.17);//for the palette

   TCanvas c1;
   c1.cd();
   deltaTvsAmplHisto->Draw("colz");
   //   deltaTvsAmpl->Draw("same");
   graph_fit_channel->Draw("sameep");
   graph_fit_fibre->Draw("sameep");
   c1.SaveAs(dir+"/deltaTvsAmpl.png");
   c1.SaveAs(dir+"/deltaTvsAmpl.pdf");


   c1.Clear();
   c1.cd();
   deltaTvsAmplHisto_channel->Draw("colz");
   //   deltaTvsAmpl->Draw("same");
   graph_fit_channel->Draw("sameep");
   c1.SaveAs(dir+"/deltaTvsAmpl_channel.png");
   c1.SaveAs(dir+"/deltaTvsAmpl_channel.pdf");

   c1.Clear();
   c1.cd();
   deltaTvsAmplHisto_fibre->Draw("colz");
   //   deltaTvsAmpl->Draw("same");
   graph_fit_fibre->Draw("sameep");
   c1.SaveAs(dir+"/deltaTvsAmpl_fibre.png");
   c1.SaveAs(dir+"/deltaTvsAmpl_channel.pdf");


   f_channel.Write();
   f_fibre.Write();

   graph_fit_channel->Write();
   graph_fit_fibre->Write();
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
