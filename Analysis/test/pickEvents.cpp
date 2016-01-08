#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <sstream> 
#include <fstream>
#include "DrawTools.h"
#include "RecoTree.h"
#include "channelInfo.h"


int main( int argc, char* argv[] ) {

  DrawTools::setStyle();



  if( argc<2 ) {
    std::cout << "ERROR. You need to specify the name of the run you want to process." << std::endl;  
    exit(1);
  }

  std::string runName = "";
  std::string tag = "V00";

  int startSample=170;
  int endSample=700;


  if( argc>1 ) {

    std::string runName_str(argv[1]);
    runName = runName_str;
     
    if( argc>2 ) {
      std::string tag_str(argv[2]);
      tag = tag_str;
      if(argc>3){
	startSample = atoi(argv[3]); 
	if(argc>4){
	  endSample = atoi(argv[4]); 
	}
      }

    }

  } else {

    std::cout << "Usage:" << std::endl;
    std::cout << "./drawCherenkov [runName] ([tag])" << std::endl;
    exit(12345);

  }

  TString runNumberString(runName);

  std::string fileName;
  fileName = "analysisTrees_"+tag+"/Reco_"+runName+".root";
  if(argc>3){
    std::ostringstream convertStart, convertEnd;
    convertStart<<startSample;
    //    convertEnd<<endSample;
    fileName = "analysisTrees_"+tag+"/Reco_"+runName+".root";
    //    fileName=  "analysisTrees_"+tag+ "/Reco_" + runName + "_start_"+convertStart.str()+"_end_"+convertEnd.str()+".root";
    fileName=  "analysisTrees_"+tag+ "/Reco_" + runName + "_optStart_"+convertStart.str()+".root";
  }
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
  nentries=1000;
  Long64_t nbytes = 0, nb = 0;

  ofstream myfile;
  myfile.open("pickEvents.dat");

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = t.LoadTree(jentry);
    if (ientry < 0) break;
    nb = t.fChain->GetEntry(jentry);   nbytes += nb;

    if((t.nClusters_hodoX1==1||t.pos_2FibClust_hodoX1>-999) && (t.nClusters_hodoX2==1||t.pos_2FibClust_hodoX2>-999) && (t.nClusters_hodoY1==1||t.pos_2FibClust_hodoY1>-999) && (t.nClusters_hodoY1==1||t.pos_2FibClust_hodoY1>-999)){//exactly one cluster, or, if there are multiple clusters, exactly one 2-fiber cluster
      if(TMath::Abs(0.5* (t.cluster_pos_corr_hodoX1+t.cluster_pos_corr_hodoX2))< 3 && TMath::Abs( 0.5* (t.cluster_pos_corr_hodoY1+t.cluster_pos_corr_hodoY2))< 3 && (t.wc_x_corr-t.cluster_pos_corr_hodoX2)<4 && (t.wc_y_corr-t.cluster_pos_corr_hodoY2)< 4  && TMath::Abs( (t.cluster_pos_corr_hodoX1-t.cluster_pos_corr_hodoX2))<1.5 &&TMath::Abs( (t.cluster_pos_corr_hodoY1-t.cluster_pos_corr_hodoY2))<1.5){

	//hadron contamination on run 4001 and 4025 imposes a cut on maxAmpl FIXME!!! no hardcoded
	if((t.run==4001 && t.cef3_maxAmpl->at(2)<650) || (t.run == 4026 && t.cef3_maxAmpl->at(2)<750) || (t.run == 4063 && t.cef3_maxAmpl->at(2)<200) || (t.run == 4102 && t.cef3_maxAmpl->at(2)<300) || (t.run ==4493 && t.cef3_maxAmpl->at(2)<1400) || (t.run ==4492 && t.cef3_maxAmpl->at(2)<1500 ) || (t.run==4064 && t.cef3_maxAmpl->at(2)<200)  || (t.run==4065 && t.cef3_maxAmpl->at(2)<220) || (t.run==4103 && t.cef3_maxAmpl->at(2)<300)|| (t.run==4104 && t.cef3_maxAmpl->at(2)<300))continue;
	//error on hv for fibre 2 on first spills?
	if(t.run==4063 && t.spill<20)continue;
	if(t.spill>=2)continue;
	std::cout<<"evtNumber != "<<t.event<<" && ";
	myfile<<t.event<<" "<<t.spill<<std::endl;
      }

    }


  }

  myfile.close();

      return 0;
}
