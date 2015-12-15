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

int main( int argc, char* argv[] ) {

  DrawTools::setStyle();

   if( argc<2 ) {
     std::cout << "ERROR. You need to specify the name of the run you want to process." << std::endl;  
     exit(1);
   }

   std::string runName = "";
   std::string tag = "V00";

   if( argc>1 ) {
     std::string tag_str(argv[1]);
     tag = tag_str;
   } else {
     std::cout << "Usage:" << std::endl;
     std::cout << "./timingStudies [tag]" << std::endl;
     exit(12345);
   }


  std::vector<float> energies;
  energies.push_back(20);
  energies.push_back(50);
  energies.push_back(100);
  energies.push_back(150);
  energies.push_back(200);

  std::vector<int> wp;
  wp.push_back(3);
  wp.push_back(5);
  wp.push_back(20);
  wp.push_back(23);
  wp.push_back(24);
  
  std::vector<int> runsSiPM;
  runsSiPM.push_back(4055);
  runsSiPM.push_back(4051);
  runsSiPM.push_back(4052);
  runsSiPM.push_back(4053);
  runsSiPM.push_back(4054);
  
  std::string constDirName = "plots_timing_SiPM";
  TString dir(constDirName);

  TGraphErrors* resVsEnergy = new TGraphErrors(0);
  std::vector<TGraphErrors*> resVsAmplitude;


  for (int i=0;i<runsSiPM.size();++i){
    TString run;
    run.Form("%d",runsSiPM[i]); 
    std::cout<<run<<std::endl;    


    TFile* file = TFile::Open(dir+"/timingStudiesSiPM_"+tag+"_"+run+".root");

    if( file==0 ) {
      std::cout << "ERROR! Din't find file for run" << runsSiPM[i] << std::endl;
      std::cout << "Exiting." << std::endl;
      exit(11);
    }
    
    TVectorD* resValueTime=(TVectorD*)file->Get("resValueTime");
    TVectorD* resErrValueTime=(TVectorD*)file->Get("resErrValueTime");
    TVectorD* resErrRelativeValueTime=(TVectorD*)file->Get("resErrRelativeValueTime");
    TVectorD* nEntriesTime=(TVectorD*)file->Get("nEntriesTime");
    
    resVsEnergy->SetPoint( i, energies[i],(*resValueTime)[wp[i]] );
    resVsEnergy->SetPointError( i, 0, (*resErrRelativeValueTime)[wp[i]]);
    
  }
  

  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  TH2D* h2_axes_1 = new TH2D( "axes_1", "", 100, -0.0, 250. , 110, 0., 1.1*(resVsEnergy->GetY())[0]);
  h2_axes_1->SetXTitle("Beam Energy [GeV]");
  h2_axes_1->SetYTitle("time_{Fibre}-time_{mcp} [ns]");
  
  
  h2_axes_1->Draw("");
  resVsEnergy->SetMarkerStyle(20);
  resVsEnergy->SetMarkerSize(1.6);
  resVsEnergy->SetMarkerColor(kBlue);
  resVsEnergy->Draw("p same");
  
  c1->SaveAs(dir+"/timingResolutionSiPM.png");
  c1->SaveAs(dir+"/timingResolutionSiPM.pdf");  

  return 0;

}
