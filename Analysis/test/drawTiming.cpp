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
   
   Long64_t nentries = t.fChain->GetEntries();
   std::cout<<"nentries:"<<nentries<<std::endl;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      Long64_t ientry = t.LoadTree(jentry);
     std::cout<<jentry<<std::endl;
      if (ientry < 0) break;
      nb = t.fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;


   }
   return 0;
}
