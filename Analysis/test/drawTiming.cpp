#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>


#include "DrawTools.h"
#include "RecoTree.h"
#include "channelInfo.h"

#include "TVectorD.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TF1.h"
#include "TLegend.h"

int main( int argc, char* argv[] ) {

  DrawTools::setStyle();

   if( argc<2 ) {
     std::cout << "ERROR. You need to specify the name of the run you want to process." << std::endl;  
     exit(1);
   }

   std::string runName = "";
   std::string tag = "V00";

   if( argc>1 ) {
     std::string runName_str(argv[1]);
     runName = runName_str;
     if( argc>2 ) {
       std::string tag_str(argv[2]);
       tag = tag_str;
     }
   }else{
     std::cout << "Usage:" << std::endl;
     std::cout << "./drawTiming [runName] ([tag])" << std::endl;
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

   if (t.fChain == 0) return 0;

   std::vector<TH1F*> histos;
 
   std::string   outfileName= "Timing_" + runName +"_"+ tag +".root";
   TFile* outFile = TFile::Open(outfileName.c_str(),"RECREATE");

   TH1F* maxAmpl_fit_cher = new TH1F("maxAmpl_cher","maxAmpl_cher",400,0,4001);
   TH1F* maxAmpl_fit_time_cher1 = new TH1F("maxAmpl_fit_time_cher1","maxAmpl_fit_time_cher1",100,50,100);
   TH1F* maxAmpl_fit_time_cher2 = new TH1F("maxAmpl_fit_time_cher2","maxAmpl_fit_time_cher2",100,0,200);
   TH2F* maxAmpl_fit_time_cher1VsMcp150 = new TH2F("cher1_vs_mcp","cher1_vs_mcp",100,17.,32.,100,65.,85.);
   TH2F* maxAmpl_fit_time_cher1VsMcp150_zoom = new TH2F("cher1_vs_mcp_zoom","cher1_vs_mcp_zoom",100,17.,32.,50,74.,85.);
   TH2F* maxAmpl_fit_time_cher2VsMcp150 = new TH2F("cher2_vs_mcp","cher2_vs_mcp",100,0,100,100,0,200);
   TH2F* deltaT1VsRatio= new TH2F("deltaT1_vs_ratio","deltaT1_vs_ratio",100,0,10,100,0,100);;
   TH2F* deltaT2VsRatio= new TH2F("deltaT2_vs_ratio","deltaT2_vs_ratio",100,0,10,100,0,200);;
  
   Long64_t nentries = t.fChain->GetEntries();
   std::cout<<"nentries:"<<nentries<<std::endl;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      Long64_t ientry = t.LoadTree(jentry);
      if (ientry < 0) break;
      nb = t.fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      maxAmpl_fit_cher->Fill(t.cef3_maxAmpl_fit_cher->at(1));
      if(t.cef3_maxAmpl_fit_cher_status->at(1)!=0)continue;
      maxAmpl_fit_time_cher1->Fill(t.cef3_maxAmpl_fit_time_cher1->at(1));
      maxAmpl_fit_time_cher2->Fill(t.cef3_maxAmpl_fit_time_cher2->at(1));
      maxAmpl_fit_time_cher1VsMcp150->Fill(t.mcp_time_at_thresh, t.cef3_maxAmpl_fit_time_cher1->at(1));
      maxAmpl_fit_time_cher1VsMcp150_zoom->Fill(t.mcp_time_at_thresh, t.cef3_maxAmpl_fit_time_cher1->at(1));
      maxAmpl_fit_time_cher2VsMcp150->Fill(t.mcp_time_at_thresh, t.cef3_maxAmpl_fit_time_cher2->at(1));
      deltaT1VsRatio->Fill(t.cef3_maxAmpl_fit_cher->at(1)/t.cef3_maxAmpl_fit->at(1),t.cef3_maxAmpl_fit_time_cher1->at(1)-t.mcp_time_at_thresh);
      deltaT2VsRatio->Fill(t.cef3_maxAmpl_fit_cher->at(1)/t.cef3_maxAmpl_fit->at(1),t.cef3_maxAmpl_fit_time_cher2->at(1)-t.mcp_time_at_thresh);
   }



   maxAmpl_fit_cher->Write();
   maxAmpl_fit_time_cher1->Write();
   maxAmpl_fit_time_cher2->Write();
   maxAmpl_fit_time_cher1VsMcp150->Write();
   maxAmpl_fit_time_cher1VsMcp150_zoom->Write();
   maxAmpl_fit_time_cher2VsMcp150->Write();
   deltaT1VsRatio->Write();
   deltaT2VsRatio->Write();
   outFile->Write();
   outFile->Close();

   return 0;
}
