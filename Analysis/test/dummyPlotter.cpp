#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <sstream> 

#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TProfile.h"

#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
#include "TChain.h"
#include "TVectorD.h"
#include "TH2.h"
#include "TGraph.h"
#include "TPaveText.h"

#include "channelInfo.h"
#include "interface/HodoCluster.h"
#include "interface/TagHelper.h"
#include "interface/EnergyCalibration.h"
#include "interface/AlignmentOfficer.h"

#include "CommonTools/interface/RunHelper.h"
#include "interface/Waveform.h"
#include "interface/WaveformFit.h"
#include "TCanvas.h"
#include "interface/Configurator.h"
#include "interface/CeF3Configurator.h"
#include "interface/CeF3InputTree.h"
#include "interface/WFTree.h"
#include "interface/Event.h"

#include "DrawTools.h"


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
   }else {

     std::cout << "Usage:" << std::endl;
     std::cout << "./makeAnalysisTree [runName] ([tag]) ([config])" << std::endl;
     exit(12345);

   }

   std::string fileName = "./analysisTrees_"+tag+"/Reco_" + runName + ".root";

   TFile* file = TFile::Open(fileName.c_str());
   if( file==0 ) {
     std::cout << "ERROR! Din't find file " << fileName << std::endl;
     std::cout << "Exiting." << std::endl;
     exit(11);
   }
   TTree* tree = (TTree*)file->Get("recoTree");

   TCanvas c1("dummy");
   tree->Draw("(cef3_time_at_frac50[1]-mcp_time_frac50):cef3_maxAmpl[1]>>h","cef3_time_at_thresh[1]>0 && cef3_time_at_thresh[1]<50 && mcp_max_amplitude>300 && cef3_maxAmpl[1]>10 && TMath::Abs(cef3_time_at_frac50[1]-mcp_time_frac50)<20","colz");
   tree->Draw("(cef3_time_at_frac50[1]-mcp_time_frac50):cef3_maxAmpl[1]>>hprof","cef3_time_at_thresh[1]>0 && cef3_time_at_thresh[1]<50 && mcp_max_amplitude>300 && cef3_maxAmpl[1]>10 && TMath::Abs(cef3_time_at_frac50[1]-mcp_time_frac50)<20","profXsamep");
   c1.SaveAs("stocazzo.png");


}
