#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <sstream> 

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

void assignValues( std::vector<float> &target, std::vector<float> source, unsigned int startPos , int skipChannel=-1, bool isOctober2015EarlyRun=false);
void assignValuesBool( std::vector<bool> &target, std::vector<bool> source, unsigned int startPos );
void computeCherenkov(std::vector<float> &cher,std::vector<float> wls);
void computeCherenkovWithFit(std::vector<float> &cher,std::vector<float> &chInt, float waveProfileInt[4],std::vector<float> cef3_maxAmpl_fit,float waveCherInt[4],std::vector<float> cef3_maxAmpl_fit_cher);
void doHodoReconstructionBool( std::vector<bool> values, int &nClusters, int *nFibres, float *pos, float fibreWidth, int clusterMaxFibres, float Cut );
std::vector<HodoCluster*> getHodoClustersBool( std::vector<bool> hodo, float fibreWidth, int nClusterMax, float Cut );

void doHodoReconstruction( std::vector<float> values, int &nClusters, int *nFibres, float *pos, float fibreWidth, int clusterMaxFibres, float Cut );
std::vector<HodoCluster*> getHodoClusters( std::vector<float> hodo, float fibreWidth, int nClusterMax, float Cut );
std::vector<HodoCluster*> getHodoClustersMultipleCuts( std::vector<float> hodo, float fibreWidth, int nClusterMax, std::vector<float> Cut );
void copyArray( int n, float *source, float *target );
void doHodoReconstructionMultipleCuts( std::vector<float> values, int &nClusters, int *nFibres, float *pos, float fibreWidth, int clusterMaxFibres, std::vector<float> Cut );
float timeSampleUnit(int drs4Freq);
CeF3_Config_t readConfiguration(std::string configName);



int main( int argc, char* argv[] ) {

  DrawTools::setStyle();


   if( argc<2 ) {
     std::cout << "ERROR. You need to specify the name of the run you want to process." << std::endl;  
     exit(1);
   }

   std::string runName = "";
   std::string tag = "V00";
   std::string config = "June2015Config";


   if( argc>1 ) {
     std::string runName_str(argv[1]);
     runName = runName_str;
     if( argc>2 ) {
       std::string tag_str(argv[2]);
       tag = tag_str;
       if(argc>3){
	 std::string config_str(argv[3]);
	 config=config_str;
       }
     }
   } else {

     std::cout << "Usage:" << std::endl;
     std::cout << "./makeAnalysisTree [runName] ([tag]) ([config])" << std::endl;
     exit(12345);

   }

   theConfiguration_=readConfiguration(config);

   //   std::string fileName = "data/Corr04_12/run_" + runName + ".root";
   std::string fileName = "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/micheli/rawData/output_run" + runName + ".root";
   char* hostName = getenv ("HOSTNAME");
   TString hostString(hostName);
   bool isLxplus = hostString.Contains("lxplus");
   if(isLxplus) fileName = "xroot://t3se01.psi.ch:1094//store/user/micheli/rawData/output_run" + runName + ".root"; 
   std::cout<<"Reading file "<<fileName<<std::endl;

   TFile* file = TFile::Open(fileName.c_str());
   if( file==0 ) {
     std::cout << "ERROR! Din't find file " << fileName << std::endl;
     std::cout << "Exiting." << std::endl;
     exit(11);
   }
   TTree* tree = (TTree*)file->Get("outputTree");

   //   TChain * chain = new TChain("tree","");
   //   tree = chain;

   //setting all the branches of inputTree
   CeF3InputTree* inputTree = new CeF3InputTree(tree);


  ///AND FINALLY THE RECO TREE/////////  
   std::string outdir = "analysisTrees_" + tag;
   system( Form("mkdir -p %s", outdir.c_str()) );

   std::string outfileName;
   outfileName= outdir + "/Reco_" + runName + ".root";
   if(theConfiguration_.addTagFileName){
     outfileName= outdir + "/Reco_" + runName + "_" + theConfiguration_.tagFileName + ".root";   
   }

   TFile* outfile = TFile::Open( outfileName.c_str(), "RECREATE" );

   TTree* outTree = new TTree("recoTree", "recoTree");

   TFile* waveformFile;
//FIXME moved to config   //   if(runName=="3076") waveformFile = TFile::Open("outWaveFormUtil_2586.root");
   waveformFile = TFile::Open(theConfiguration_.waveformFile.c_str());

   TFile* cherFile;
   bool CHER_RUN=false;
   if(theConfiguration_.cherRun==1){
     CHER_RUN=true;
     cherFile = TFile::Open(theConfiguration_.cherFile.c_str());
   }

   bool MCP_RUN=false;
   if(theConfiguration_.mcpRun==1){
     MCP_RUN=true;
   }
   if(runName=="3076") {
     gROOT->ProcessLine("gErrorIgnoreLevel=kWarning");// BE CAREFUL!this suppress warnings in root!!!Used because of a bug in minuit when doing minimization
   }


   TFile* noiseFile;
   noiseFile = TFile::Open("outWaveFormUtil_2539.root");

   //FIXME move somewhere else
   //and now naming the stuff for the recoTree//
   outTree->Branch( "run", &inputTree->runNumber, "run/i" );
   outTree->Branch( "spill", &inputTree->spillNumber, "spill/i" );
   outTree->Branch( "event", &inputTree->evtNumber, "event/i" );


   //Original
   std::vector<float> cef3( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3", &cef3 ); //for backwards compatibility
   std::vector<float> cef3_corr( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_corr", &cef3_corr );


   std::vector<float> cef3_maxAmpl( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_maxAmpl", &cef3_maxAmpl );
   std::vector<float> cef3_maxAmpl_time( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_maxAmpl_time", &cef3_maxAmpl_time );
   std::vector<float> cef3_time_at_frac50( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_time_at_frac50", &cef3_time_at_frac50 );
   std::vector<float> cef3_time_at_thresh( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_time_at_thresh", &cef3_time_at_thresh );


   std::vector<float> cef3_maxAmpl_cher( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_maxAmpl_cher", &cef3_maxAmpl_cher );
   std::vector<float> cef3_maxAmpl_wls( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_maxAmpl_wls", &cef3_maxAmpl_wls );


   std::vector<float> cef3_maxAmpl_fit( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_maxAmpl_fit", &cef3_maxAmpl_fit );
   std::vector<float> cef3_maxAmpl_fit_cher( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_maxAmpl_fit_cher", &cef3_maxAmpl_fit_cher );
   std::vector<float> cef3_maxAmpl_fit_cher_status( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_maxAmpl_fit_cher_status", &cef3_maxAmpl_fit_cher_status );
   std::vector<float> cef3_maxAmpl_fit_time_cher1( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_maxAmpl_fit_time_cher1", &cef3_maxAmpl_fit_time_cher1 );
   std::vector<float> cef3_maxAmpl_fit_time_cher2( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_maxAmpl_fit_time_cher2", &cef3_maxAmpl_fit_time_cher2 );
   std::vector<float> cef3_maxAmpl_fit_time_cher3( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_maxAmpl_fit_time_cher3", &cef3_maxAmpl_fit_time_cher3 );
   std::vector<float> cef3_maxAmpl_fit_corr( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_maxAmpl_fit_corr", &cef3_maxAmpl_fit_corr );
   std::vector<float> cef3_chaInt( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_chaInt", &cef3_chaInt );
   std::vector<float> cef3_chaInt_cher( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_chaInt_cher", &cef3_chaInt_cher );
   std::vector<float> cef3_chaInt_wls( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_chaInt_wls", &cef3_chaInt_wls );


   std::vector<float> bgo( BGO_CHANNELS, -1. );
   outTree->Branch( "bgo", &bgo );

   //Original Intercalibrated
   std::vector<float> cef3_maxAmpl_corr( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_maxAmpl_corr", &cef3_maxAmpl_corr );
   std::vector<float> cef3_chaInt_corr( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_chaInt_corr", &cef3_chaInt_corr );

   std::vector<float> bgo_corr( BGO_CHANNELS, -1. );
   outTree->Branch( "bgo_corr", &bgo_corr );


   float xTable;
   outTree->Branch( "xTable", &xTable, "xTable/F");
   float yTable;
   outTree->Branch( "yTable", &yTable, "yTable/F");

   float beamEnergy;
   outTree->Branch( "beamEnergy", &beamEnergy, "xbeamEnergy/F");

   float angle; //aka BeamTilt..
   outTree->Branch("angle", &angle, "angle/F");

   float HVCeF3;
   outTree->Branch( "HVCeF3", &HVCeF3, "HVCeF3/F");
   float HVBGO;
   outTree->Branch( "HVBGO", &HVBGO, "HVBGO/F");

   float xBeam;
   outTree->Branch( "xBeam", &xBeam, "xBeam/F");
   float yBeam;
   outTree->Branch( "yBeam", &yBeam, "yBeam/F");

   int nClusters_hodoX1;
   outTree->Branch( "nClusters_hodoX1", &nClusters_hodoX1, "nClusters_hodoX1/I" );
   int nFibres_hodoX1[HODOX1_CHANNELS];
   outTree->Branch( "nFibres_hodoX1", nFibres_hodoX1, "nFibres_hodoX1[nClusters_hodoX1]/I" );
   float pos_hodoX1[HODOX1_CHANNELS];
   outTree->Branch( "pos_hodoX1", pos_hodoX1, "pos_hodoX1[nClusters_hodoX1]/F" );
   float pos_corr_hodoX1[HODOX1_CHANNELS];
   outTree->Branch( "pos_corr_hodoX1", pos_corr_hodoX1, "pos_corr_hodoX1[nClusters_hodoX1]/F" );

   int nClusters_hodoY1;
   outTree->Branch( "nClusters_hodoY1", &nClusters_hodoY1, "nClusters_hodoY1/I" );
   int nFibres_hodoY1[HODOY1_CHANNELS];
   outTree->Branch( "nFibres_hodoY1", nFibres_hodoY1, "nFibres_hodoY1[nClusters_hodoY1]/I" );
   float pos_hodoY1[HODOY1_CHANNELS];
   outTree->Branch( "pos_hodoY1", pos_hodoY1, "pos_hodoY1[nClusters_hodoY1]/F" );
   float pos_corr_hodoY1[HODOY1_CHANNELS];
   outTree->Branch( "pos_corr_hodoY1", pos_corr_hodoY1, "pos_corr_hodoY1[nClusters_hodoY1]/F" );

   int nClusters_hodoX2;
   outTree->Branch( "nClusters_hodoX2", &nClusters_hodoX2, "nClusters_hodoX2/I" );
   int nFibres_hodoX2[HODOX2_CHANNELS];
   outTree->Branch( "nFibres_hodoX2", nFibres_hodoX2, "nFibres_hodoX2[nClusters_hodoX2]/I" );
   float pos_hodoX2[HODOX2_CHANNELS];
   outTree->Branch( "pos_hodoX2", pos_hodoX2, "pos_hodoX2[nClusters_hodoX2]/F" );
   float pos_corr_hodoX2[HODOX2_CHANNELS];
   outTree->Branch( "pos_corr_hodoX2", pos_corr_hodoX2, "pos_corr_hodoX2[nClusters_hodoX2]/F" );

   int nClusters_hodoY2;
   outTree->Branch( "nClusters_hodoY2", &nClusters_hodoY2, "nClusters_hodoY2/I" );
   int nFibres_hodoY2[HODOY2_CHANNELS];
   outTree->Branch( "nFibres_hodoY2", nFibres_hodoY2, "nFibres_hodoY2[nClusters_hodoY2]/I" );
   float pos_hodoY2[HODOY2_CHANNELS];
   outTree->Branch( "pos_hodoY2", pos_hodoY2, "pos_hodoY2[nClusters_hodoY2]/F" );
   float pos_corr_hodoY2[HODOY2_CHANNELS];
   outTree->Branch( "pos_corr_hodoY2", pos_corr_hodoY2, "pos_corr_hodoY2[nClusters_hodoY2]/F" );

   int nClusters_hodoSmallX;
   outTree->Branch( "nClusters_hodoSmallX", &nClusters_hodoSmallX, "nClusters_hodoSmallX/I" );
   int nFibres_hodoSmallX[HODOSMALLX_CHANNELS];
   outTree->Branch( "nFibres_hodoSmallX", nFibres_hodoSmallX, "nFibres_hodoSmallX[nClusters_hodoSmallX]/I" );
   float pos_hodoSmallX[HODOSMALLX_CHANNELS];
   outTree->Branch( "pos_hodoSmallX", pos_hodoSmallX, "pos_hodoSmallX[nClusters_hodoSmallX]/F" );
   float pos_corr_hodoSmallX[HODOSMALLX_CHANNELS];
   outTree->Branch( "pos_corr_hodoSmallX", pos_corr_hodoSmallX, "pos_corr_hodoSmallX[nClusters_hodoSmallX]/F" );

   int nClusters_hodoSmallY;
   outTree->Branch( "nClusters_hodoSmallY", &nClusters_hodoSmallY, "nClusters_hodoSmallY/I" );
   int nFibres_hodoSmallY[HODOSMALLY_CHANNELS];
   outTree->Branch( "nFibres_hodoSmallY", nFibres_hodoSmallY, "nFibres_hodoSmallY[nClusters_hodoSmallY]/I" );
   float pos_hodoSmallY[HODOSMALLY_CHANNELS];
   outTree->Branch( "pos_hodoSmallY", pos_hodoSmallY, "pos_hodoSmallY[nClusters_hodoSmallY]/F" );
   float pos_corr_hodoSmallY[HODOSMALLY_CHANNELS];
   outTree->Branch( "pos_corr_hodoSmallY", pos_corr_hodoSmallY, "pos_corr_hodoSmallY[nClusters_hodoSmallY]/F" );



   
   std::vector<int> nTDCHits( 4, -1. );
   outTree->Branch( "nTDCHits", &nTDCHits );




   float pos_2FibClust_hodoX1;
   outTree->Branch( "pos_2FibClust_hodoX1", &pos_2FibClust_hodoX1, "pos_2FibClust_hodoX1/F" );
   float pos_2FibClust_corr_hodoX1;
   outTree->Branch( "pos_2FibClust_corr_hodoX1", &pos_2FibClust_corr_hodoX1, "pos_2FibClust_corr_hodoX1/F" );

   float pos_2FibClust_hodoY1;
   outTree->Branch( "pos_2FibClust_hodoY1", &pos_2FibClust_hodoY1, "pos_2FibClust_hodoY1/F" );
   float pos_2FibClust_corr_hodoY1;
   outTree->Branch( "pos_2FibClust_corr_hodoY1", &pos_2FibClust_corr_hodoY1, "pos_2FibClust_corr_hodoY1/F" );

   float pos_2FibClust_hodoX2;
   outTree->Branch( "pos_2FibClust_hodoX2", &pos_2FibClust_hodoX2, "pos_2FibClust_hodoX2/F" );
   float pos_2FibClust_corr_hodoX2;
   outTree->Branch( "pos_2FibClust_corr_hodoX2", &pos_2FibClust_corr_hodoX2, "pos_2FibClust_corr_hodoX2/F" );

   float pos_2FibClust_hodoY2;
   outTree->Branch( "pos_2FibClust_hodoY2", &pos_2FibClust_hodoY2, "pos_2FibClust_hodoY2/F" );
   float pos_2FibClust_corr_hodoY2;
   outTree->Branch( "pos_2FibClust_corr_hodoY2", &pos_2FibClust_corr_hodoY2, "pos_2FibClust_corr_hodoY2/F" );


   float cluster_pos_hodoX1;
   outTree->Branch( "cluster_pos_hodoX1", &cluster_pos_hodoX1, "cluster_pos_hodoX1/F" );
   float cluster_pos_corr_hodoX1;
   outTree->Branch( "cluster_pos_corr_hodoX1", &cluster_pos_corr_hodoX1, "cluster_pos_corr_hodoX1/F" );
   float cluster_pos_hodoX2;
   outTree->Branch( "cluster_pos_hodoX2", &cluster_pos_hodoX2, "cluster_pos_hodoX2/F" );
   float cluster_pos_corr_hodoX2;
   outTree->Branch( "cluster_pos_corr_hodoX2", &cluster_pos_corr_hodoX2, "cluster_pos_corr_hodoX2/F" );

   float cluster_pos_hodoY1;
   outTree->Branch( "cluster_pos_hodoY1", &cluster_pos_hodoY1, "cluster_pos_hodoY1/F" );
   float cluster_pos_corr_hodoY1;
   outTree->Branch( "cluster_pos_corr_hodoY1", &cluster_pos_corr_hodoY1, "cluster_pos_corr_hodoY1/F" );
   float cluster_pos_hodoY2;
   outTree->Branch( "cluster_pos_hodoY2", &cluster_pos_hodoY2, "cluster_pos_hodoY2/F" );
   float cluster_pos_corr_hodoY2;
   outTree->Branch( "cluster_pos_corr_hodoY2", &cluster_pos_corr_hodoY2, "cluster_pos_corr_hodoY2/F" );



   float wc_x;
   outTree->Branch( "wc_x", &wc_x, "wc_x/F");
   float wc_y;
   outTree->Branch( "wc_y", &wc_y, "wc_y/F");
   float wc_x_corr;
   outTree->Branch( "wc_x_corr", &wc_x_corr, "wc_x_corr/F");
   float wc_y_corr;
   outTree->Branch( "wc_y_corr", &wc_y_corr, "wc_y_corr/F");

   float mcp_time_frac50;
   outTree->Branch( "mcp_time_frac50", &mcp_time_frac50, "mcp_time_frac50/F");

   float mcp_time_at_thresh;
   outTree->Branch( "mcp_time_at_thresh", &mcp_time_at_thresh, "mcp_time_at_thresh/F");

   float mcp_max_amplitude;
   outTree->Branch( "mcp_max_amplitude", &mcp_max_amplitude, "mcp_max_amplitude/F");

   int nino_triggered;
   outTree->Branch( "nino_triggered", &nino_triggered, "nino_triggered/I");

   float nino_LEtime;
   outTree->Branch( "nino_LEtime", &nino_LEtime, "nino_LEtime/F");

   float nino_LEchi2;
   outTree->Branch( "nino_LEchi2", &nino_LEchi2, "nino_LEchi2/F");

   float nino_maxAmpl;
   outTree->Branch( "nino_maxAmpl", &nino_maxAmpl, "nino_maxAmpl/F");

   float nino_chInt;
   outTree->Branch( "nino_chInt", &nino_chInt, "nino_chInt/F");
        
   float nino_timeOverThresh;
   outTree->Branch( "nino_timeOverThresh", &nino_timeOverThresh, "nino_timeOverThresh/F");


   int nentries = tree->GetEntries();
   if(nentries>theConfiguration_.limitNEntries && theConfiguration_.limitNEntries>0)nentries=theConfiguration_.limitNEntries;
   //   nentries=5000;
   RunHelper::getBeamPosition( runName, xBeam, yBeam );

   if(nentries>0)tree->GetEntry(0);     
   bool isOctober2015EarlyRun =  (inputTree->runNumber > 3900. && inputTree->runNumber<4200);//in first runs of october 2015 there was one broken channel (ch 2)
   bool isOctober2015LateRun  = (inputTree->runNumber>4490 && inputTree->runNumber<4550);

   bool isOctober2015PMTRun=false;
   bool isOctober2015EarlyMAPDRun=false;
   bool isOctober2015EarlySiPMRun=false;
   if(isOctober2015EarlyRun){
     isOctober2015PMTRun=inputTree->runNumber>=4063 && inputTree->runNumber<=4104;
     isOctober2015EarlyMAPDRun=inputTree->runNumber>=3971 && inputTree->runNumber<=4033;
     isOctober2015EarlySiPMRun=inputTree->runNumber>=4051 && inputTree->runNumber<=4055;
   }

   TH1F* signalOverNoise=new TH1F("signalOverNoise","signalOverNoise",200,0,40);

   std::cout <<"nentries:"<<nentries << std::endl;
   std::cout<< "using waveform average shape from file "<<waveformFile->GetName()<<std::endl;
   //get reference waveform for fits
   std::vector<TProfile*>  waveProfile;
   float waveProfileInt[4];
   std::vector<TProfile*>  waveCher;
   float waveCherInt[4];
   float maxTimeCher[4]; 


   int counterChannel=-1;
   for (unsigned int iCh=0; iCh<CEF3_CHANNELS; iCh++) {
     int iChannel=iCh;
     TString name="prof";
     name+=iChannel;
     TString istring;
     istring+=iChannel;
     if(inputTree->digi_frequency==1)  waveProfile.push_back(new TProfile(name,name,1024,0,0.4e-6));
     else waveProfile.push_back(new TProfile(name,name,1024,0,0.2e-6));
     TGraph* graph = (TGraph*)waveformFile->Get("waveform_"+istring);
     double X[graph->GetN()],Y[graph->GetN()];
     float maxInWLSRegion=0;
     for(int iP=0;iP<graph->GetN();iP++){
       graph->GetPoint(iP,X[iP],Y[iP]);
       waveProfile.at(iChannel)->Fill(X[iP],Y[iP]);
       if(maxInWLSRegion<Y[iP] && iP > theConfiguration_.startSample) maxInWLSRegion=Y[iP];
      }
     //     TH1F* histo = (TH1F*)waveformFile->Get("waveform_histo_"+istring);
     //     for(int i=0;i<=histo->GetNbinsX();++i) waveProfile.at(iCh)->Fill(histo->GetBinCenter(i), histo->GetBinContent(i));
     //     waveProfile.at(iChannel)->Scale(1./waveProfile.at(iChannel)->GetMaximum());
     if(iChannel!=2)    waveProfile.at(iChannel)->Scale(1./maxInWLSRegion);//in quartz fibres max could be cherenkov
     else waveProfile.at(iChannel)->Scale(1./waveProfile.at(iChannel)->GetMaximum());

     waveProfileInt[iChannel]=waveProfile.at(iChannel)->Integral(theConfiguration_.startSample,theConfiguration_.endSample);

     if(CHER_RUN){
       name+="cher";
       if(inputTree->digi_frequency==1)       waveCher.push_back(new TProfile(name,name,1024,0,0.4e-6));
       else        waveCher.push_back(new TProfile(name,name,1024,0,0.2e-6));
       TGraph* graph2 = (TGraph*)cherFile->Get("waveform_"+istring);
       double X2[graph2->GetN()],Y2[graph2->GetN()];
       for(int iP=0;iP<graph2->GetN();iP++){
	 graph2->GetPoint(iP,X2[iP],Y2[iP]);
	 waveCher.at(iChannel)->Fill(X2[iP],Y2[iP]);
       }
       waveCher.at(iChannel)->Scale(1./waveCher.at(iChannel)->GetMaximum());
       maxTimeCher[iChannel]=waveProfile.at(iChannel)->GetBinCenter(waveProfile.at(iChannel)->GetMaximumBin());
     }
   }//end loop get waveforms

   //noise histo for pedestal runs
   //   TH1F* noiseHisto[4][1024];
   TH1F* noiseHisto[4];
   TH2F* noiseHistoVsTime[4];
   float noiseValueVsChannel[4][1024];


   //noise
   std::vector<TProfile*>  waveNoiseProfile;

   if(runName=="2539"){//pedestal run
     for (unsigned int iCh=0; iCh<CEF3_CHANNELS; iCh++) {
       TString name="noiseProfile_";
       name+=iCh;
       TString istring;
       istring+=iCh;
       waveNoiseProfile.push_back(new TProfile(name,name,1024,0,1024));
       TGraph* graph = (TGraph*)noiseFile->Get("waveform_"+istring);
       double X[graph->GetN()],Y[graph->GetN()];
       for(int iP=0;iP<graph->GetN();iP++){
	 graph->GetPoint(iP,X[iP],Y[iP]);
	 waveNoiseProfile.at(iCh)->Fill(iP,Y[iP]);
       }
       noiseHisto[iCh]= new TH1F("histoNoise_"+istring,"histoNoise_"+istring,2000,-10,10);
       noiseHistoVsTime[iCh]= new TH2F("histoNoiseVsTime_"+istring,"histoNoiseVsTime_"+istring,1024,0,1024,2000,-10,10);
       for(int jj=0;jj<1024;jj++){
	 noiseValueVsChannel[iCh][jj]=0.;
	 //noiseHisto[iCh][j]=new TH1F("histoNoiseVsTime_"+istring+"_"+j,"histoNoiseVsTime_"+istring+"_"+j,);
       }
     }
   }

   
   WFTree outWFTree(theConfiguration_.cef3Channels.size(), 1024);

   float  timeMinusTimeAtMax[NDIGISAMPLES];
   float  mean[NFIBERS][NDIGISAMPLES];
   float  mean_tight[NFIBERS][NDIGISAMPLES];
   float  time[NDIGISAMPLES];
   
   float meanTimeAtMax[NFIBERS];
   TGraph* meanWaveGraphsForPlots[NFIBERS];  
   TGraph* meanWaveGraphs[NFIBERS];  
   TGraph* meanWaveGraphsForPlots_tight[NFIBERS];  
   TGraph* meanWaveGraphs_tight[NFIBERS];  
   
   if(theConfiguration_.fillWFtree){

     for (int i=0;i<NFIBERS;++i){
       meanTimeAtMax[i]=0;
       for (int j=0;j<NDIGISAMPLES;++j){
	 mean[i][j]=0;
	 mean_tight[i][j]=0;
       }
     }
     

   }


   int WFoutTreeEntries=0;
   int WFoutTreeTightEntries=0;

   for(int iEntry=0; iEntry<nentries; ++iEntry ) {
     
    tree->GetEntry( iEntry );     
    if( iEntry %  1000 == 0 ) std::cout << "Entry: " << iEntry << " / " << nentries << std::endl;


    //set the tag for calibration
    std::string theBeamEnergy = Form("%.0f",inputTree->BeamEnergy);
    TagHelper tagHelper(tag,theBeamEnergy);
    EnergyCalibration cef3Calib(tagHelper.getCeF3FileName());
    EnergyCalibration bgoCalib(tagHelper.getBGOFileName());
    AlignmentOfficer alignOfficer(tagHelper.getAlignmentFileName());
    




//     computeCherenkov(cef3_chaInt_cher,cef3_chaInt_wls);

// if(CHER_RUN)     computeCherenkovWithFit(cef3_chaInt_cher, cef3_chaInt,waveProfileInt,cef3_maxAmpl_fit,waveCherInt,cef3_maxAmpl_fit_cher);



     //     assignValues( bgo_corr, *BGOvalues, 0 );
     //     bgoCalib.applyCalibration(bgo_corr);

     std::vector<bool> hodoX1_values(HODOX1_CHANNELS, -1.);
     std::vector<bool> hodoY1_values(HODOY1_CHANNELS, -1.);
     assignValuesBool( hodoX1_values, *inputTree->HODOX1, 0. );
     assignValuesBool( hodoY1_values, *inputTree->HODOY1, 0. );


     std::vector<bool> hodoX2_values(HODOX2_CHANNELS, -1.);
     std::vector<bool> hodoY2_values(HODOY2_CHANNELS, -1.);
     assignValuesBool( hodoX2_values, *inputTree->HODOX2, 0 );
     assignValuesBool( hodoY2_values, *inputTree->HODOY2, 0 );


     std::vector<float> hodoSmallX_values(theConfiguration_.pedMeanX.size(), -1.);
     std::vector<float> hodoSmallY_values(theConfiguration_.pedMeanY.size(), -1.);
     assignValues( hodoSmallY_values, *inputTree->HODOSMALLvalues, 0 );
     assignValues( hodoSmallX_values, *inputTree->HODOSMALLvalues, theConfiguration_.pedMeanX.size() );



     // hodo cluster reconstruction
     int clusterMaxFibres = 4;
     doHodoReconstructionBool( hodoX1_values    , nClusters_hodoX1    , nFibres_hodoX1    , pos_hodoX1    , 0.5, clusterMaxFibres, 0. );
     doHodoReconstructionBool( hodoY1_values    , nClusters_hodoY1    , nFibres_hodoY1    , pos_hodoY1    , 0.5, clusterMaxFibres, 0. );
     doHodoReconstructionBool( hodoX2_values    , nClusters_hodoX2    , nFibres_hodoX2    , pos_hodoX2    , 0.5, clusterMaxFibres , 0.);
     doHodoReconstructionBool( hodoY2_values    , nClusters_hodoY2    , nFibres_hodoY2    , pos_hodoY2    , 0.5, clusterMaxFibres, 0. );

     std::vector<float> pedMeanX,pedMeanY, pedSigmaX, pedSigmaY,cutValuesX,cutValuesY;
     for(int i=0;i<theConfiguration_.pedMeanY.size();++i){
       pedMeanY.push_back(theConfiguration_.pedMeanY[i]);
       pedSigmaY.push_back(theConfiguration_.pedSigmaY[i]);
       pedMeanX.push_back(theConfiguration_.pedMeanX[i]);
       pedSigmaX.push_back(theConfiguration_.pedSigmaX[i]);

       cutValuesX.push_back(pedMeanX[i]+5*pedSigmaX[i]);
       cutValuesY.push_back(pedMeanY[i]+5*pedSigmaY[i]);
     }
     doHodoReconstructionMultipleCuts( hodoSmallX_values, nClusters_hodoSmallX, nFibres_hodoSmallX, pos_hodoSmallX, 1.0, 1, cutValuesX);
     doHodoReconstructionMultipleCuts( hodoSmallY_values, nClusters_hodoSmallY, nFibres_hodoSmallY, pos_hodoSmallY, 1.0, 1, cutValuesY);

     copyArray( nClusters_hodoX1, pos_hodoX1, pos_corr_hodoX1 );
     copyArray( nClusters_hodoY1, pos_hodoY1, pos_corr_hodoY1 );
     copyArray( nClusters_hodoX2, pos_hodoX2, pos_corr_hodoX2 );
     copyArray( nClusters_hodoY2, pos_hodoY2, pos_corr_hodoY2 );
     copyArray( nClusters_hodoSmallX, pos_hodoSmallX, pos_corr_hodoSmallX );
     copyArray( nClusters_hodoSmallY, pos_hodoSmallY, pos_corr_hodoSmallY );

     alignOfficer.fix("hodoX1", nClusters_hodoX1, pos_corr_hodoX1);
     alignOfficer.fix("hodoY1", nClusters_hodoY1, pos_corr_hodoY1);
     alignOfficer.fix("hodoX2", nClusters_hodoX2, pos_corr_hodoX2);
     alignOfficer.fix("hodoY2", nClusters_hodoY2, pos_corr_hodoY2);


     xTable = inputTree->TableX;
     yTable = inputTree->TableY;
     HVCeF3 = inputTree->CeF3HV;
     HVBGO =  inputTree->BGOHV;
     beamEnergy = inputTree->BeamEnergy;
     angle = inputTree->BeamTilt;

     nTDCHits[0] = inputTree->nTdcHits->at(0);
     nTDCHits[1] = inputTree->nTdcHits->at(1);
     nTDCHits[2] = inputTree->nTdcHits->at(2);
     nTDCHits[3] = inputTree->nTdcHits->at(3);

     wc_x = inputTree->TDCreco->at(0);
     wc_y = inputTree->TDCreco->at(1);

     if( inputTree->runNumber>=170 ) wc_y = -wc_y;

     wc_x_corr = wc_x + alignOfficer.getOffset("wc_x");
     wc_y_corr = wc_y + alignOfficer.getOffset("wc_y");

 



     int posOf2FibClustX1=0;
     int nrOf2FibreClustersX1 = 0;
     for( int i=0; i<nClusters_hodoX1; ++i ) {
       if( nFibres_hodoX1[i]==2  )  {  
	 ++nrOf2FibreClustersX1;
	 posOf2FibClustX1 = i ; } 
     }
     if( nrOf2FibreClustersX1 == 1){
     pos_2FibClust_hodoX1 =   pos_hodoX1[ posOf2FibClustX1] ;
     }else{
       pos_2FibClust_hodoX1 =  -999 ;}

     if(nClusters_hodoX1==1){
       cluster_pos_hodoX1 = pos_hodoX1[0];
     }else if(nrOf2FibreClustersX1==1){
       cluster_pos_hodoX1 = pos_hodoX1[ posOf2FibClustX1];
     }else{ cluster_pos_hodoX1 = -999;}


     int posOf2FibClustX2=0;
     int nrOf2FibreClustersX2 = 0;
     for( int i=0; i<nClusters_hodoX2; ++i ) {
       if( nFibres_hodoX2[i]==2  )  {  
	 ++nrOf2FibreClustersX2;
	 posOf2FibClustX2 = i ; } 
     }
     if( nrOf2FibreClustersX2 == 1){
     pos_2FibClust_hodoX2 =   pos_hodoX2[ posOf2FibClustX2] ;
     }else{
       pos_2FibClust_hodoX2 =  -999 ;}

     if(nClusters_hodoX2==1){
       cluster_pos_hodoX2 = pos_hodoX2[0];
     }else if(nrOf2FibreClustersX2==1){
       cluster_pos_hodoX2 = pos_hodoX2[ posOf2FibClustX2];
     }else{ cluster_pos_hodoX2 = -999;}

     int posOf2FibClustY1=0;
     int nrOf2FibreClustersY1 = 0;
     for( int i=0; i<nClusters_hodoY1; ++i ) {
       if( nFibres_hodoY1[i]==2  )  {  
	 ++nrOf2FibreClustersY1;
	 posOf2FibClustY1 = i ; } 
     }
     if( nrOf2FibreClustersY1 == 1){
     pos_2FibClust_hodoY1 =   pos_hodoY1[ posOf2FibClustY1] ;
     }else{
       pos_2FibClust_hodoY1 =  -999 ;}
     
     if(nClusters_hodoY1==1){
       cluster_pos_hodoY1 = pos_hodoY1[0];
     }else if(nrOf2FibreClustersY1==1){
       cluster_pos_hodoY1 = pos_hodoY1[ posOf2FibClustY1];
     }else{ cluster_pos_hodoY1 = -999;}
   

     int posOf2FibClustY2=0;
     int nrOf2FibreClustersY2 = 0;
     for( int i=0; i<nClusters_hodoY2; ++i ) {
       if( nFibres_hodoY2[i]==2  )  {  
	 ++nrOf2FibreClustersY2;
	 posOf2FibClustY2 = i ; } 
     }
     if( nrOf2FibreClustersY2 == 1){
     pos_2FibClust_hodoY2 =   pos_hodoY2[ posOf2FibClustY2] ;
     }else{
       pos_2FibClust_hodoY2 =  -999 ;}

     if(nClusters_hodoY2==1){
       cluster_pos_hodoY2 = pos_hodoY2[0];
     }else if(nrOf2FibreClustersY2==1){
       cluster_pos_hodoY2 = pos_hodoY2[ posOf2FibClustY2];
     }else{ cluster_pos_hodoY2 = -999;}


     pos_2FibClust_corr_hodoX1 = pos_2FibClust_hodoX1 + alignOfficer.getOffset("hodoX1");
     pos_2FibClust_corr_hodoY1 = pos_2FibClust_hodoY1 + alignOfficer.getOffset("hodoY1");   
     pos_2FibClust_corr_hodoX2 = pos_2FibClust_hodoX2 + alignOfficer.getOffset("hodoX2");   
     pos_2FibClust_corr_hodoY2 = pos_2FibClust_hodoY2 + alignOfficer.getOffset("hodoY2");


     cluster_pos_corr_hodoX1 = cluster_pos_hodoX1 + alignOfficer.getOffset("hodoX1");
     cluster_pos_corr_hodoY1 = cluster_pos_hodoY1 + alignOfficer.getOffset("hodoY1");   
     cluster_pos_corr_hodoX2 = cluster_pos_hodoX2 + alignOfficer.getOffset("hodoX2");   
     cluster_pos_corr_hodoY2 = cluster_pos_hodoY2 + alignOfficer.getOffset("hodoY2");




    //waveform creation
    std::vector<Waveform*> waveform;
    waveform.clear();
    for (unsigned int i=0; i<CEF3_CHANNELS+1*isOctober2015EarlyRun+MCP_RUN; i++) {
      waveform.push_back(new Waveform());
      waveform.at(i)->clear();
    
    }

    Waveform* waveform_ninoAnalog;
    if(isOctober2015LateRun){waveform_ninoAnalog=new Waveform();waveform_ninoAnalog->clear();}
    
    float timeOfTheEvent=inputTree->digi_time_at_1000_bare_noise_sub->at(theConfiguration_.triggerChannel);//synchronizing time of events with time of trigger
    float shiftTime=0;
    shiftTime = theConfiguration_.meanTriggerTime-timeOfTheEvent;


    int shiftSample=round(shiftTime/(1e9*timeSampleUnit(inputTree->digi_frequency)));
    shiftSample=-shiftSample;
    //FIX ME move to config    if(runName=="2539")shiftSample=0;

    bool fillWFoutTree=false;
    bool fillWFoutTreeTight=false;

    int ch=-1;
    for (int i=0;i<1024*(CEF3_CHANNELS+1*isOctober2015EarlyRun);++i){
      //      if(digi_value_ch->at(i) > 5)continue; //just to avoid not useful channels
      int iChannel=inputTree->digi_value_ch->at(i);
     //      std::cout<<i<<" "<<inputTree->digi_value_ch->at(i)<<" ";

      
      bool skipChannel=true;
      for (int j=0;j<theConfiguration_.cef3Channels.size();++j){
	if(iChannel==theConfiguration_.cef3Channels[j]){
	  skipChannel=false;
	  if(i==1024*iChannel)ch++;
	  break;
	}
      }
      if(skipChannel) continue;
      
      

	
     bool doWeWantToShift=theConfiguration_.syncChannels[ch];


     //     if(i==1024*iChannel)      std::cout<<doWeWantToShift<<"<-do we iChannel->"<<iChannel<<" ich"<<ch<< " shift:"<<shiftSample<<std::endl;

      if(inputTree->digi_max_amplitude->at(inputTree->digi_value_ch->at(i))>10000 || inputTree->digi_max_amplitude->at(inputTree->digi_value_ch->at(i))<0)continue;
      int iSample=i;
      if(i+shiftSample*doWeWantToShift>1023*inputTree->digi_value_ch->at(i) && i+shiftSample*doWeWantToShift<(1023+(1024*inputTree->digi_value_ch->at(i)))){
	iSample=i+shiftSample*doWeWantToShift;
	//		std::cout<<inputTree->digi_value_ch->at(i)<<" "<<i<<" "<<shiftSample*doWeWantToShift<<std::endl;
      }
      waveform.at(ch)->addTimeAndSample((i-1024*inputTree->digi_value_ch->at(i))*timeSampleUnit(inputTree->digi_frequency),inputTree->digi_value_bare_noise_sub->at(iSample));

      if(i==0 && (iEntry % theConfiguration_.prescaleWFtree == 0) ){
	if((nClusters_hodoX1==1||pos_2FibClust_hodoX1>-999) && (nClusters_hodoX2==1||pos_2FibClust_hodoX2>-999) && (nClusters_hodoY1==1||pos_2FibClust_hodoY1>-999) && (nClusters_hodoY1==1||pos_2FibClust_hodoY1>-999)){//exactly one cluster, or, if there are multiple clusters, exactly one 2-fiber cluster
	  

	  if(TMath::Abs(0.5* (cluster_pos_corr_hodoX1+cluster_pos_corr_hodoX2))< 3 && TMath::Abs( 0.5* (cluster_pos_corr_hodoY1+cluster_pos_corr_hodoY2))< 3 && (wc_x_corr-cluster_pos_corr_hodoX2)<4 && (wc_y_corr-cluster_pos_corr_hodoY2)< 4  && TMath::Abs( (cluster_pos_corr_hodoX1-cluster_pos_corr_hodoX2))<1.5 &&TMath::Abs( (cluster_pos_corr_hodoY1-cluster_pos_corr_hodoY2))<1.5){
	  WFoutTreeEntries++; 
	  fillWFoutTree=true;


	    if(TMath::Abs(0.5* (cluster_pos_corr_hodoX1+cluster_pos_corr_hodoX2))< 1 && TMath::Abs( 0.5* (cluster_pos_corr_hodoY1+cluster_pos_corr_hodoY2))< 1){//let's see the waveforms for events in 1x1 strip
	      WFoutTreeTightEntries++; 
	      fillWFoutTreeTight=true;
	    }    
	  }
	}
      }




      if(theConfiguration_.fillWFtree && fillWFoutTree  ){
	outWFTree.WF_ch[(i-1024*inputTree->digi_value_ch->at(i))+1024*ch] = ch;
	outWFTree.WF_time[(i-1024*inputTree->digi_value_ch->at(i))+1024*ch] = (i-1024*inputTree->digi_value_ch->at(i))*timeSampleUnit(inputTree->digi_frequency);
	outWFTree.WF_val[(i-1024*inputTree->digi_value_ch->at(i))+1024*ch] = inputTree->digi_value_bare_noise_sub->at(iSample);
	outWFTree.index[(i-1024*inputTree->digi_value_ch->at(i))+1024*ch] = WFoutTreeEntries;
	outWFTree.evtNumber[(i-1024*inputTree->digi_value_ch->at(i))+1024*ch] = iEntry;

	//average waveform
	mean[ch][i-1024*inputTree->digi_value_ch->at(i)]+=(float)(inputTree->digi_value_bare_noise_sub->at(iSample));
	if(fillWFoutTreeTight) mean_tight[ch][i-1024*inputTree->digi_value_ch->at(i)]+=(float)(inputTree->digi_value_bare_noise_sub->at(iSample));
   
	  
	if(i<1024)time[(i-1024*inputTree->digi_value_ch->at(i))]=inputTree->digi_value_time->at(i);
	if(inputTree->digi_time_at_max_bare_noise_sub->at(inputTree->digi_value_ch->at(i))>0 && inputTree->digi_time_at_max_bare_noise_sub->at(inputTree->digi_value_ch->at(i))<400 && i==1024*iChannel){
	  meanTimeAtMax[ch]+=inputTree->digi_time_at_max_bare_noise_sub->at(inputTree->digi_value_ch->at(i)); 
	}
      }

      //           std::cout<<"i:"<<inputTree->digi_value_ch->at(i)<<" time:"<<(i-1024*inputTree->digi_value_ch->at(i))*timeSampleUnit(inputTree->digi_frequency)<<" value"<<inputTree->digi_value_bare_noise_sub->at(iSample)<<std::endl;
      if(runName=="2539"){
	//	std::cout<<i<<" "<<inputTree->digi_value_bare_noise_sub->at(iSample)<<" "<<(waveNoiseProfile.at(inputTree->digi_value_ch->at(i)))->GetBinContent(iSample-1024*inputTree->digi_value_ch->at(i))<<std::endl;
	if(i>theConfiguration_.startSample+1023*inputTree->digi_value_ch->at(i) && i<theConfiguration_.endSample+1023*inputTree->digi_value_ch->at(i))	noiseHisto[ch]->Fill(inputTree->digi_value_bare_noise_sub->at(iSample)-waveNoiseProfile.at(ch)->GetBinContent(i-1024*inputTree->digi_value_ch->at(iSample)));
	noiseValueVsChannel[ch][i-1024*ch]+=(inputTree->digi_value_bare_noise_sub->at(iSample)-(waveNoiseProfile.at(ch))->GetBinContent(i-1024*inputTree->digi_value_ch->at(iSample)));
      }

    }


    //mcp info    
    if(MCP_RUN){
      timeOfTheEvent=inputTree->digi_time_at_1000_bare_noise_sub->at(theConfiguration_.mcpTriggerChannel);//synchronizing time of events with time of trigger
      shiftTime = theConfiguration_.meanTriggerTime-timeOfTheEvent;
      shiftSample=round(shiftTime/(1e9*timeSampleUnit(inputTree->digi_frequency)));
      shiftSample=-shiftSample;

      for (int i=1024*CEF3_CHANNELS;i<1024*(theConfiguration_.mcpChannel+1);++i){
	if(inputTree->digi_value_ch->at(i)!=theConfiguration_.mcpChannel)continue;
	if(inputTree->digi_max_amplitude->at(inputTree->digi_value_ch->at(i))>10000 || inputTree->digi_max_amplitude->at(inputTree->digi_value_ch->at(i))<0)continue;

	int iSample=i;
	if(i+shiftSample*theConfiguration_.syncMcp>1023*inputTree->digi_value_ch->at(i) && i+shiftSample*theConfiguration_.syncMcp<(1023+(1024*inputTree->digi_value_ch->at(i)))){
	  iSample=i+shiftSample*theConfiguration_.syncMcp;
	}
	waveform.at(CEF3_CHANNELS)->addTimeAndSample((i-1024*inputTree->digi_value_ch->at(i))*timeSampleUnit(inputTree->digi_frequency),inputTree->digi_value_bare_noise_sub->at(iSample));
      }

      

      Waveform::max_amplitude_informations wave_max_bare = waveform.at(CEF3_CHANNELS)->max_amplitude(4,900,5);
      std::vector<float> crossingTimes;
      if(wave_max_bare.time_at_max>0)    crossingTimes  = waveform.at(CEF3_CHANNELS)->time_at_threshold((const float)10.*timeSampleUnit(inputTree->digi_frequency), 100e-9,theConfiguration_.mcpTimingThresh,4);
      if(crossingTimes.size()>0)    mcp_time_at_thresh=crossingTimes[0]*1.e9;
      else mcp_time_at_thresh=-999;
      
      if(wave_max_bare.time_at_max>0)    mcp_time_frac50=waveform.at(CEF3_CHANNELS)->time_at_frac(0.,(float)100.e-9,0.5,wave_max_bare,5)*1.e9;
      else mcp_time_frac50=-999;

      mcp_max_amplitude=wave_max_bare.max_amplitude;
    }



    ///now we have all the needed waveforms, compute some quantities
    int counterChannel=-1;
    for (unsigned int i=0; i<CEF3_CHANNELS; i++) {      
	int iChannel=i;


	Waveform::max_amplitude_informations wave_max_cher = waveform.at(iChannel)->max_amplitude(30,theConfiguration_.startSample,5);
	if(wave_max_cher.max_amplitude>0)cef3_maxAmpl_cher[iChannel] = wave_max_cher.max_amplitude;

	Waveform::max_amplitude_informations wave_max_wls = waveform.at(iChannel)->max_amplitude(theConfiguration_.startSample,900,5);
	if(wave_max_wls.max_amplitude>0)cef3_maxAmpl_wls[iChannel] = wave_max_wls.max_amplitude;


	//WaveformFit using mean Waveform
	ROOT::Math::Minimizer* minimizer;
	int sampleIntegral=theConfiguration_.startSample;
	if(iChannel==2)sampleIntegral=50;
	Waveform::max_amplitude_informations wave_max = waveform.at(iChannel)->max_amplitude(sampleIntegral,900,5);
	Waveform::baseline_informations wave_pedestal = waveform.at(iChannel)->baseline(5,34);

	if(i==1 && !isOctober2015LateRun)	if(wave_max.max_amplitude>0 && wave_pedestal.rms>0)	signalOverNoise->Fill(wave_max.max_amplitude/wave_pedestal.rms);

	//fit for NINO in october 2015 runs
	if(isOctober2015LateRun && iChannel==1){
	  std::pair<float,float> timeInfo = WaveformFit::GetTimeLE(waveform.at(iChannel),wave_pedestal,300,1,1,80,120,timeSampleUnit(inputTree->digi_frequency));//window without sync
	  nino_LEtime=timeInfo.first*1.e9;
	  nino_LEchi2=timeInfo.second;

	  float timeThr = WaveformFit::GetTimeAboveThr(waveform.at(iChannel),600,80,120,timeSampleUnit(inputTree->digi_frequency));//check if there is a signal in nino (signal >650 adc)
	  if(timeThr>0) nino_triggered=1;
	  else nino_triggered=0;

	  nino_maxAmpl=inputTree->digi_max_amplitude_bare_noise_sub->at(6);//FIXME move to config!

	  if(wave_max.max_amplitude>0 && inputTree->digi_pedestal_bare_noise_sub_rms->at(6)>0)	signalOverNoise->Fill(nino_maxAmpl/inputTree->digi_pedestal_bare_noise_sub_rms->at(6));//FIXME move to config!

	  if(nino_LEtime>0){
	    //	    std::cout<<(int)(timeInfo.first/timeSampleUnit(inputTree->digi_frequency))-5<<std::endl; 
	    //	    nino_chInt=waveform.at(iChannel)->charge_integrated(80,120);
	    //	    nino_chInt=waveform.at(iChannel)->charge_integrated((int)(timeInfo.first/timeSampleUnit(inputTree->digi_frequency))-5,(int)(timeInfo.first/timeSampleUnit(inputTree->digi_frequency))+25);
	    std::pair<float,float> integralInfo = WaveformFit::GetSignalIntegral(waveform.at(iChannel),200,70);
	    nino_chInt= integralInfo.first;
	    nino_timeOverThresh=integralInfo.second*timeSampleUnit(inputTree->digi_frequency)*1.e9;
	  } else {
	    nino_chInt=-999;
	    nino_timeOverThresh=-999;
	  }
	}
	
	
	//WaveformFit on Cher using mean Waveform
	ROOT::Math::Minimizer* minimizerCher;
	
	//	if(isOctober2015LateRun && iChannel!=2) continue;
	if(wave_max.max_amplitude<0) continue;
	if(iChannel!=2) {
	    WaveformFit::fitWaveformSimple(waveform.at(iChannel),waveProfile.at(iChannel),200,200,wave_max,wave_pedestal,minimizer, true, theConfiguration_.startSample, theConfiguration_.endSample);
	    if(CHER_RUN){
	      int cherStart=150;
	      if(inputTree->digi_frequency==0)cherStart=125;
	      if(inputTree->digi_time_at_max_noise_sub->at(iChannel)<-22) WaveformFit::fitWaveformSimple(waveform.at(iChannel),waveCher.at(iChannel),200,200,wave_max,wave_pedestal,minimizerCher, true, cherStart, theConfiguration_.startSample);
	      else  WaveformFit::fitWaveformSimplePlusTime(waveform.at(iChannel),waveCher.at(iChannel),50,100,wave_max,wave_pedestal,minimizerCher, true,cherStart, theConfiguration_.startSample);
	    }
	}else{ 
	    WaveformFit::fitWaveformSimple(waveform.at(iChannel),waveProfile.at(iChannel),200,200,wave_max,wave_pedestal,minimizer);
	}
      
      const double* par=minimizer->X();
      cef3_maxAmpl_fit[iChannel]=par[0];
      



      if(CHER_RUN){
      const double* parCher=minimizerCher->X();
	cef3_maxAmpl_fit_cher[iChannel]=parCher[0];
	cef3_maxAmpl_fit_time_cher1[iChannel]=(maxTimeCher[iChannel]+parCher[1])*1.e9;
	cef3_maxAmpl_fit_time_cher2[iChannel]=waveform.at(iChannel)->time_at_frac(0.,(float)100.e-9,0.5,wave_max,7)*1.e9;
	cef3_maxAmpl_fit_cher_status[iChannel]=minimizerCher->Status();	
      }

      //timing var
      Waveform::max_amplitude_informations wave_max_bare = waveform.at(iChannel)->max_amplitude(4,900,5);
      std::vector<float> crossingTimes;
      if(wave_max_bare.time_at_max>0)    crossingTimes  = waveform.at(iChannel)->time_at_threshold((const float)10.*timeSampleUnit(inputTree->digi_frequency), 100e-9,theConfiguration_.timingThresh[iChannel],2);
      if(crossingTimes.size()>0)    cef3_time_at_thresh[iChannel]=crossingTimes[0]*1.e9;
      else cef3_time_at_thresh[iChannel]=-999;



      if(wave_max_bare.time_at_max>0) cef3_time_at_frac50[iChannel]=waveform.at(iChannel)->time_at_frac(0.,(float)100.e-9,0.5,wave_max_bare,5)*1.e9;
      else cef3_time_at_frac50[iChannel]=-999;


    }
    
    //BACKWARDS COMPATIBILITY///
    assignValues( cef3, *inputTree->digi_charge_integrated_bare_noise_sub, CEF3_START_CHANNEL,2,isOctober2015EarlyRun);     
    assignValues( cef3_corr, *inputTree->digi_charge_integrated_bare_noise_sub, CEF3_START_CHANNEL ,2,isOctober2015EarlyRun);  
    
    assignValues( cef3_maxAmpl_fit_corr, cef3_maxAmpl_fit, CEF3_START_CHANNEL );  
    cef3Calib.applyCalibration(cef3_maxAmpl_fit_corr);
    assignValues( cef3_maxAmpl, *inputTree->digi_max_amplitude_bare_noise_sub, CEF3_START_CHANNEL,2,isOctober2015EarlyRun);
    assignValues( cef3_maxAmpl_time, *inputTree->digi_time_at_max_noise_sub, CEF3_START_CHANNEL,2,isOctober2015EarlyRun);
    assignValues( cef3_chaInt, *inputTree->digi_charge_integrated_bare_noise_sub, CEF3_START_CHANNEL,2,isOctober2015EarlyRun);



    if(theConfiguration_.fillWFtree && fillWFoutTree )	{
      outWFTree.maxAmpl0=inputTree->digi_max_amplitude_bare_noise_sub->at(0);
      outWFTree.maxAmpl1=inputTree->digi_max_amplitude_bare_noise_sub->at(1);
      outWFTree.maxAmpl2=inputTree->digi_max_amplitude_bare_noise_sub->at(2);
      outWFTree.maxAmpl3=inputTree->digi_max_amplitude_bare_noise_sub->at(3);
      outWFTree.deltaT=cef3_time_at_frac50[1] -mcp_time_frac50;
      outWFTree.tree_->AutoSave("FlushBaskets");
      outWFTree.Fill();

    }



     outTree->Fill();
    
   
   } // for entries, end of loop



   outfile->cd();

      if(theConfiguration_.fillWFtree)outWFTree.Write("wf", "wf_tree");
   //  outWFTree.Write("wf", "wf_tree");
   if(runName=="2539"){//pedestal run
     for (unsigned int iCh=0; iCh<CEF3_CHANNELS; iCh++) {
       for(int jj=0;jj<1024;jj++){
	 noiseHistoVsTime[iCh]->Fill(jj,noiseValueVsChannel[iCh][jj]/nentries);
       }
       noiseHistoVsTime[iCh]->Write();
       noiseHisto[iCh]->Write();
     }
   }

  //  for(int j=0;j<1024;j++)  std::cout<<"j"<<j<<" "<<noiseValueVsChannel[1][j]<<std::endl;

   if(theConfiguration_.fillWFtree){

     float max[NFIBERS]={0};
     float max_tight[NFIBERS]={0};
     for (int i=0;i<NFIBERS;++i){
       meanTimeAtMax[i]/=WFoutTreeEntries;
       for (int j=0;j<NDIGISAMPLES;++j){
	 mean[i][j]/=WFoutTreeEntries;
	 mean_tight[i][j]/=WFoutTreeTightEntries;

	 //	 if(i==0)	 std::cout<<mean[i][j]<<std::endl;
	 if(j>800) continue;//let's avoid weard shapes at the end of the pulse
	 if(max[i]<mean[i][j] ){
	   max[i]=mean[i][j];
	 }
	 if(max_tight[i]<mean_tight[i][j]){
	   max_tight[i]=mean_tight[i][j];
	 }
       }
     }
     
     for (int i=0;i<NFIBERS;++i){
       meanWaveGraphs[i]=new TGraph(1024, time, mean[i]);
       meanWaveGraphs_tight[i]=new TGraph(1024, time, mean_tight[i]);
       TString fiber;
       fiber.Form("%d",i);
       meanWaveGraphs[i]->SetName("waveform_"+fiber);
       meanWaveGraphs_tight[i]->SetName("waveform_tight_"+fiber);

       for (int j=0;j<NDIGISAMPLES;++j){  
	 timeMinusTimeAtMax[j]=(1e9*(time[j])-meanTimeAtMax[i]);
	 mean[i][j]*=0.25; //be careful!this is to do meanWaveGraphsForPlots in milliVolt, change it if you need to reuse it
	 mean_tight[i][j]*=0.25; //be careful!this is to do meanWaveGraphsForPlots in milliVolt, change it if you need to reuse it
       }
       meanWaveGraphsForPlots[i]=new TGraph(1024, timeMinusTimeAtMax, mean[i]);
       meanWaveGraphsForPlots[i]->SetName("waveform_plots_"+fiber);
     

       meanWaveGraphs[i]->Write();
       meanWaveGraphsForPlots[i]->Write();

       meanWaveGraphsForPlots_tight[i]=new TGraph(1024, timeMinusTimeAtMax, mean_tight[i]);
       meanWaveGraphsForPlots_tight[i]->SetName("waveform_plots_tight_"+fiber);
     

       meanWaveGraphs_tight[i]->Write();
       meanWaveGraphsForPlots_tight[i]->Write();

     
       TCanvas c1("c_"+fiber);
//       TH2F* axis=new TH2F("axis","axis",100,-80-20*isOctober2015EarlyMAPDRun,280,100,0,max[i]/4*1.10);
//       TH2F* axis2=new TH2F("axis","axis",100,-80-20*isOctober2015EarlyMAPDRun,280,100,0,max_tight[i]/4*1.10);
       TH2F* axis=new TH2F("axis","axis",100,-100,280,100,0,max[i]/4*1.10);
       TH2F* axis2=new TH2F("axis","axis",100,-100,280,100,0,max_tight[i]/4*1.10);



       std::string energy(Form("%.0f", inputTree->BeamEnergy));
       TPaveText* pave = DrawTools::getLabelTop_expOnXaxis(energy+" GeV Electron Beam");
       
       //    TH2F* axis=new TH2F("axis","axis",100,-10,190,100,0,28);
       axis->GetXaxis()->SetTitle("Time [ns]");
       axis->GetYaxis()->SetTitle("Signal Amplitude [mV]");
       axis2->GetXaxis()->SetTitle("Time [ns]");
       axis2->GetYaxis()->SetTitle("Signal Amplitude [mV]");



       axis->Draw();
       meanWaveGraphsForPlots[i]->Draw("plsame");
       pave->Draw("same");
       c1.Write("pulseShape"+fiber);
       std::string name= Form("plots/pulseShape"+fiber+"_%d",inputTree->runNumber);
       if(theConfiguration_.addTagFileName)name= Form("plots/pulseShape"+fiber+"_"+theConfiguration_.tagFileName+"_%d",inputTree->runNumber);
       c1.SaveAs(Form("%s.png",name.c_str()));
       c1.SaveAs(Form("%s.pdf",name.c_str()));

       c1.Clear();
       axis2->Draw();
       meanWaveGraphsForPlots_tight[i]->Draw("plsame");
       pave->Draw("same");
       c1.Write("pulseShape_tight"+fiber);
       name= Form("plots/pulseShape_tight"+fiber+"_%d",inputTree->runNumber);
       if(theConfiguration_.addTagFileName)       name= Form("plots/pulseShape_tight"+fiber+"_"+theConfiguration_.tagFileName+"_%d",inputTree->runNumber);
       c1.SaveAs(Form("%s.png",name.c_str()));
       c1.SaveAs(Form("%s.pdf",name.c_str()));

       
     }
   }

   signalOverNoise->Write();


   outTree->Write();
   outfile->Close();

   std::cout << "-> Analysis Tree saved in: " << outfile->GetName() << std::endl;
   return 0;

}
  


void assignValues( std::vector<float> &target, std::vector<float> source, unsigned int startPos , int skipChannel, bool isOctober2015EarlyRun) {
  for( unsigned i=0; i<target.size()+1*(isOctober2015EarlyRun && skipChannel>-1); ++i ) {
    int iCh=i;
    if(i==skipChannel && isOctober2015EarlyRun){
      continue;
    }else if (i>skipChannel && isOctober2015EarlyRun){
      iCh--;
    }
    target[iCh] = source[startPos+i];
  }
}


void assignValuesBool( std::vector<bool> &target, std::vector<bool> source, unsigned int startPos ) {

  for( unsigned i=0; i<target.size(); ++i ) 
    target[i] = source[startPos+i];

}

void computeCherenkov(std::vector<float> &cher,std::vector<float> wls){
  TVectorD integrals(4);//values obtained from run2778
  integrals[0]=0.079;
  integrals[1]=0.039;
  integrals[2]=0.20;
  integrals[3]=0.042;

  for (int i=0;i<cher.size();++i){
    float chargeWlsUnderCher=(wls[i]*integrals[i])/(1-integrals[i]);//amount of wls charge under cherenkov peak. estimated by fit
    cher[i]=cher[i]-chargeWlsUnderCher;
  }
  
}

void computeCherenkovWithFit(std::vector<float> &cher,std::vector<float> &chInt, float waveProfileInt[4],std::vector<float> cef3_maxAmpl_fit,float waveCherInt[4],std::vector<float> cef3_maxAmpl_fit_cher){

  TVectorD integrals(4);//values obtained from run2778
  integrals[0]=0.079;
  integrals[1]=0.039;
  integrals[2]=0.20;
  integrals[3]=0.042;

  for (int i=0;i<cher.size();++i){
    float chWls=cef3_maxAmpl_fit[i]*waveProfileInt[i];
    //float cherPlusWlsFast=chInt[i]-chWls;
    float cherPlusWlsFast=cef3_maxAmpl_fit_cher[i]*waveCherInt[i];
    float chargeWlsUnderCher=(chWls*integrals[i])/(1-integrals[i]);//amount of wls charge under cherenkov peak. estimated by fit
    cher[i]=cherPlusWlsFast-chargeWlsUnderCher;
//if(TMath::Abs(cher[i])<100){    std::cout<<"chWls:"<<chWls<<" cherPlusWlsFast"<<cherPlusWlsFast<<" cef3_maxAmpl_fit_cher:"<<cef3_maxAmpl_fit_cher[i]<<" int:"<<waveProfileInt[i]<<std::endl;
//      std::cout<<"total "<<chInt[i]<<" cher[i]:"<<cher[i]<<" chargeWlsUnderCher"<<chargeWlsUnderCher<<std::endl;
//    }
  }


}

std::vector<HodoCluster*> getHodoClusters( std::vector<float> hodo, float fibreWidth, int nClusterMax, float Cut ) {

  std::vector<HodoCluster*> clusters;

  HodoCluster* currentCluster = new HodoCluster( hodo.size(), fibreWidth );

  for( unsigned i=0; i<hodo.size(); ++i ) {

    if( hodo[i] > Cut) { // hit

      if( currentCluster->getSize() < nClusterMax ) {

        currentCluster->addFibre( i );

      } else {

        clusters.push_back( currentCluster ); // store old one
        currentCluster = new HodoCluster( hodo.size(), fibreWidth );   // create a new one
        currentCluster->addFibre( i );        // get that fibre!

      }

    } else { // as soon as you find a hole
      
      if( currentCluster->getSize() > 0 ) {
     
        clusters.push_back( currentCluster ); // store old one
        currentCluster = new HodoCluster( hodo.size(), fibreWidth );   // create a new one

      }

    }


  } // for fibres


  if( currentCluster->getSize()>0 )
    clusters.push_back( currentCluster ); // store last cluster


  return clusters;

}

std::vector<HodoCluster*> getHodoClustersMultipleCuts( std::vector<float> hodo, float fibreWidth, int nClusterMax, std::vector<float> Cut ) {

  std::vector<HodoCluster*> clusters;

  HodoCluster* currentCluster = new HodoCluster( hodo.size(), fibreWidth );

  for( unsigned i=0; i<hodo.size(); ++i ) {

    if( hodo[i] > Cut[i]) { // hit

      if( currentCluster->getSize() < nClusterMax ) {

        currentCluster->addFibre( i );

      } else {

        clusters.push_back( currentCluster ); // store old one
        currentCluster = new HodoCluster( hodo.size(), fibreWidth );   // create a new one
        currentCluster->addFibre( i );        // get that fibre!

      }

    } else { // as soon as you find a hole
      
      if( currentCluster->getSize() > 0 ) {
     
        clusters.push_back( currentCluster ); // store old one
        currentCluster = new HodoCluster( hodo.size(), fibreWidth );   // create a new one

      }

    }


  } // for fibres


  if( currentCluster->getSize()>0 )
    clusters.push_back( currentCluster ); // store last cluster


  return clusters;

}




void doHodoReconstruction( std::vector<float> values, int &nClusters, int *nFibres, float *pos, float fibreWidth, int clusterMaxFibres, float Cut ) {

  std::vector<HodoCluster*> clusters = getHodoClusters( values, fibreWidth, clusterMaxFibres, Cut );

  nClusters = clusters.size();
  for( unsigned i=0; i<clusters.size(); ++i ) {
    nFibres[i] = clusters[i]->getSize();
    pos[i] = clusters[i]->getPosition();
  }

}

void doHodoReconstructionMultipleCuts( std::vector<float> values, int &nClusters, int *nFibres, float *pos, float fibreWidth, int clusterMaxFibres, std::vector<float> Cut ) {

  std::vector<HodoCluster*> clusters = getHodoClustersMultipleCuts( values, fibreWidth, clusterMaxFibres, Cut );

  nClusters = clusters.size();
  for( unsigned i=0; i<clusters.size(); ++i ) {
    nFibres[i] = clusters[i]->getSize();
    pos[i] = clusters[i]->getPosition();
  }

}


void copyArray( int n, float *source, float *target ) {

  for( unsigned i=0; i<n; ++i ) 
    target[i] = source[i];

}













std::vector<HodoCluster*> getHodoClustersBool( std::vector<bool> hodo, float fibreWidth, int nClusterMax, float Cut ) {

  std::vector<HodoCluster*> clusters;

  HodoCluster* currentCluster = new HodoCluster( hodo.size(), fibreWidth );

  for( unsigned i=0; i<hodo.size(); ++i ) {

    if( hodo[i] > Cut) { // hit

      if( currentCluster->getSize() < nClusterMax ) {

        currentCluster->addFibre( i );

      } else {

        clusters.push_back( currentCluster ); // store old one
        currentCluster = new HodoCluster( hodo.size(), fibreWidth );   // create a new one
        currentCluster->addFibre( i );        // get that fibre!

      }

    } else { // as soon as you find a hole
      
      if( currentCluster->getSize() > 0 ) {
     
        clusters.push_back( currentCluster ); // store old one
        currentCluster = new HodoCluster( hodo.size(), fibreWidth );   // create a new one

      }

    }


  } // for fibres


  if( currentCluster->getSize()>0 )
    clusters.push_back( currentCluster ); // store last cluster


  return clusters;

}


void doHodoReconstructionBool( std::vector<bool> values, int &nClusters, int *nFibres, float *pos, float fibreWidth, int clusterMaxFibres, float Cut ) {

  std::vector<HodoCluster*> clusters = getHodoClustersBool( values, fibreWidth, clusterMaxFibres, Cut );

  nClusters = clusters.size();
  for( unsigned i=0; i<clusters.size(); ++i ) {
    nFibres[i] = clusters[i]->getSize();
    pos[i] = clusters[i]->getPosition();
  }

}

float timeSampleUnit(int drs4Freq)
{
  if (drs4Freq == 0)
    return 0.2E-9;
  else if (drs4Freq == 1)
    return 0.4E-9;
  else if (drs4Freq == 2)
    return 1.E-9;
  return -999.;
}

//  LocalWords:  waveProfileInt computeCherenkovWithFit


CeF3_Config_t readConfiguration(std::string configName){
   //config file
  Configurator* configurator_ = new Configurator();
  std::string fileName = Form ("./config/%s.xml",configName.c_str());
  configurator_->xmlFileName=fileName.c_str();
  configurator_->Init();

  CeF3_Config_t conf;
  string content = Configurable::getElementContent (*configurator_, "cef3Channels",configurator_->root_element) ;
   Configurator::GetVecInt(content,conf.cef3Channels);

   content=Configurable::getElementContent (*configurator_, "syncChannels",configurator_->root_element) ;
   Configurator::GetVecInt(content,conf.syncChannels);

   content=Configurable::getElementContent (*configurator_, "pedMeanY",configurator_->root_element) ;
   Configurator::GetVecFloat(content,conf.pedMeanY);

   content=Configurable::getElementContent (*configurator_, "pedSigmaY",configurator_->root_element) ;
   Configurator::GetVecFloat(content,conf.pedSigmaY);

   content=Configurable::getElementContent (*configurator_, "pedMeanX",configurator_->root_element) ;
   Configurator::GetVecFloat(content,conf.pedMeanX);

   content=Configurable::getElementContent (*configurator_, "pedSigmaX",configurator_->root_element) ;
   Configurator::GetVecFloat(content,conf.pedSigmaX);


   content=Configurable::getElementContent (*configurator_, "timingThresh",configurator_->root_element) ;
   Configurator::GetVecInt(content,conf.timingThresh);


   conf.triggerChannel= Configurator::GetInt(Configurable::getElementContent(*configurator_,"triggerChannel",configurator_->root_element));


   conf.meanTriggerTime= Configurator::GetDouble(Configurable::getElementContent(*configurator_,"meanTriggerTime",configurator_->root_element));

   conf.mcpRun= Configurator::GetInt(Configurable::getElementContent(*configurator_,"mcpRun",configurator_->root_element));
   
   conf.mcpChannel= Configurator::GetInt(Configurable::getElementContent(*configurator_,"mcpChannel",configurator_->root_element));

   conf.mcpTriggerChannel= Configurator::GetInt(Configurable::getElementContent(*configurator_,"mcpTriggerChannel",configurator_->root_element));

   conf.mcpMeanTriggerTime= Configurator::GetDouble(Configurable::getElementContent(*configurator_,"mcpMeanTriggerTime",configurator_->root_element));

   conf.syncMcp= Configurator::GetInt(Configurable::getElementContent(*configurator_,"syncMcp",configurator_->root_element));

   conf.mcpTimingThresh= Configurator::GetInt(Configurable::getElementContent(*configurator_,"mcpTimingThresh",configurator_->root_element));

   conf.waveformFile=Configurable::getElementContent (*configurator_, "waveformFile",configurator_->root_element) ;

   conf.startSample= Configurator::GetInt(Configurable::getElementContent(*configurator_,"startSample",configurator_->root_element));

   conf.endSample= Configurator::GetInt(Configurable::getElementContent(*configurator_,"endSample",configurator_->root_element));

   conf.cherRun= Configurator::GetInt(Configurable::getElementContent(*configurator_,"cherRun",configurator_->root_element));

   conf.cherFile=Configurable::getElementContent (*configurator_, "cherFile",configurator_->root_element) ;

   conf.addTagFileName= Configurator::GetInt(Configurable::getElementContent(*configurator_,"addTagFileName",configurator_->root_element));

   conf.tagFileName=Configurable::getElementContent (*configurator_, "tagFileName",configurator_->root_element) ;

   conf.fillWFtree = Configurator::GetInt(Configurable::getElementContent(*configurator_,"fillWFtree",configurator_->root_element));

   conf.prescaleWFtree = Configurator::GetInt(Configurable::getElementContent(*configurator_,"prescaleWFtree",configurator_->root_element));

   conf.limitNEntries= Configurator::GetInt(Configurable::getElementContent(*configurator_,"limitNEntries",configurator_->root_element));

   return conf;

}
