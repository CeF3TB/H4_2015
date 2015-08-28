#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cmath>

#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
#include "TChain.h"
#include "TVectorD.h"

#include "channelInfo.h"
#include "interface/HodoCluster.h"
#include "interface/TagHelper.h"
#include "interface/EnergyCalibration.h"
#include "interface/AlignmentOfficer.h"

#include "CommonTools/interface/RunHelper.h"
#include "interface/Waveform.h"
#include "interface/WaveformFit.h"
#include "TCanvas.h"

void assignValues( std::vector<float> &target, std::vector<float> source, unsigned int startPos );
void assignValuesBool( std::vector<bool> &target, std::vector<bool> source, unsigned int startPos );
void computeCherenkov(std::vector<float> &cher,std::vector<float> wls);
void doHodoReconstructionBool( std::vector<bool> values, int &nClusters, int *nFibres, float *pos, float fibreWidth, int clusterMaxFibres, float Cut );
std::vector<HodoCluster*> getHodoClustersBool( std::vector<bool> hodo, float fibreWidth, int nClusterMax, float Cut );

void doHodoReconstruction( std::vector<float> values, int &nClusters, int *nFibres, float *pos, float fibreWidth, int clusterMaxFibres, float Cut );
std::vector<HodoCluster*> getHodoClusters( std::vector<float> hodo, float fibreWidth, int nClusterMax, float Cut );
std::vector<HodoCluster*> getHodoClustersMultipleCuts( std::vector<float> hodo, float fibreWidth, int nClusterMax, std::vector<float> Cut );
void copyArray( int n, float *source, float *target );
void doHodoReconstructionMultipleCuts( std::vector<float> values, int &nClusters, int *nFibres, float *pos, float fibreWidth, int clusterMaxFibres, std::vector<float> Cut );
float timeSampleUnit(int drs4Freq);




int main( int argc, char* argv[] ) {


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

   } else {

     std::cout << "Usage:" << std::endl;
     std::cout << "./makeAnalysisTree [runName] ([tag])" << std::endl;
     exit(12345);

   }

   //   std::string fileName = "data/Corr04_12/run_" + runName + ".root";
   std::string fileName = "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/micheli/rawData/output_run" + runName + ".root";
   TFile* file = TFile::Open(fileName.c_str());
   if( file==0 ) {
     std::cout << "ERROR! Din't find file " << fileName << std::endl;
     std::cout << "Exiting." << std::endl;
     exit(11);
   }
   TTree* tree = (TTree*)file->Get("outputTree");

   //   TChain * chain = new TChain("tree","");
   //   tree = chain;
   Int_t           fCurrent;
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
 

   // Declaration of leaf types
   UInt_t          runNumber;
   UInt_t          spillNumber;
   UInt_t          evtNumber;
   UInt_t          digi_frequency;
   std::vector<int>     *digi_value_ch;
   std::vector<float>   *digi_value_time;
   std::vector<float>   *digi_value_bare_noise_sub;
   std::vector<float>   *BGOvalues;
   std::vector<float>   *SCINTvalues;
   std::vector<float>   *TDCreco;
   std::vector<float>   *digi_charge_integrated;
   std::vector<float>   *digi_charge_integrated_bare_noise_sub;
   std::vector<float>   *digi_max_amplitude;
   std::vector<float>   *digi_max_amplitude_bare_noise_sub;
   std::vector<float>   *digi_charge_integrated_bare_noise_sub_fast;
   std::vector<float>   *digi_charge_integrated_bare_noise_sub_slow;
   std::vector<float>   *digi_pedestal;
   std::vector<float>   *digi_pedestal_rms;
   std::vector<float>   *digi_time_at_frac30;
   std::vector<float>   *digi_time_at_frac50_bare_noise_sub;
   std::vector<float>   *digi_time_at_1000_bare_noise_sub;
   std::vector<float>   *digi_time_at_max;
   std::vector<bool>    *HODOX1;
   std::vector<bool>    *HODOX2;
   std::vector<bool>    *HODOY1;
   std::vector<bool>    *HODOY2;
   std::vector<float>   *HODOSMALLvalues;
   Float_t         TableX;
   Float_t         TableY;
   Float_t         CeF3HV;
   Float_t         BGOHV;
   Float_t         BeamEnergy;
   Float_t         BeamTilt;
   Int_t           IsPhysics;
   std::vector<int>     *nTdcHits;
   std::vector<float>   *digi_charge_integrated_sub;
   std::vector<float>   *digi_max_amplitude_sub;
   std::vector<float>   *digi_pedestal_sub;
   std::vector<float>   *digi_pedestal_rms_sub;
   std::vector<float>   *digi_charge_integrated_corr1;
   std::vector<float>   *digi_max_amplitude_corr1;
   std::vector<float>   *digi_charge_integrated_corr2;
   std::vector<float>   *digi_max_amplitude_corr2;

   std::vector<float>   *digi_max_amplitude_bare;
   std::vector<float>   *digi_charge_integrated_bare;
   std::vector<float>   *digi_charge_integrated_frac10;
   std::vector<float>   *digi_charge_integrated_frac30;
   std::vector<float>   *digi_charge_integrated_frac50;


   // List of branches
   TBranch        *b_digi_value_ch;   //!
   TBranch        *b_digi_value_time;   //!
   TBranch        *b_digi_value_bare_noise_sub;   //!
   TBranch        *b_runNumber;   //!
   TBranch        *b_digi_frequency;   //!
   TBranch        *b_spillNumber;   //!
   TBranch        *b_evtNumber;   //!
   TBranch        *b_BGOvalues;   //!
   TBranch        *b_SCINTvalues;   //!
   TBranch        *b_TDCreco;   //!
   TBranch        *b_digi_charge_integrated;   //!
   TBranch        *b_digi_charge_integrated_bare_noise_sub;   //!
   TBranch        *b_digi_charge_integrated_bare_noise_sub_fast;   //!
   TBranch        *b_digi_charge_integrated_bare_noise_sub_slow;   //!
   TBranch        *b_digi_max_amplitude;   //!
   TBranch        *b_digi_max_amplitude_bare_noise_sub;   //!
   TBranch        *b_digi_pedestal;   //!
   TBranch        *b_digi_pedestal_rms;   //!
   TBranch        *b_digi_time_at_frac30;   //!
   TBranch        *b_digi_time_at_frac50_bare_noise_sub;   //!
   TBranch        *b_digi_time_at_1000_bare_noise_sub;   //!
   TBranch        *b_digi_time_at_max;   //!
   TBranch        *b_HODOX1;   //!
   TBranch        *b_HODOX2;   //!
   TBranch        *b_HODOY1;   //!
   TBranch        *b_HODOY2;   //!
   TBranch        *b_HODOSMALLvalues;   //!
   TBranch        *b_TableX;   //!
   TBranch        *b_TableY;   //!
   TBranch        *b_CeF3HV;   //!
   TBranch        *b_BGOHV;   //!
   TBranch        *b_BeamEnergy;   //!
   TBranch        *b_BeamTilt;   //!
   TBranch        *b_IsPhysics;   //!
   TBranch        *b_nTdcHits;   //!
   TBranch        *b_digi_charge_integrated_sub;   //!
   TBranch        *b_digi_max_amplitude_sub;   //!
   TBranch        *b_digi_pedestal_sub;   //!
   TBranch        *b_digi_pedestal_rms_sub;   //!
   TBranch        *b_digi_charge_integrated_corr1;   //!
   TBranch        *b_digi_max_amplitude_corr1;   //!
   TBranch        *b_digi_charge_integrated_corr2;   //!
   TBranch        *b_digi_max_amplitude_corr2;   //!

   TBranch *b_digi_max_amplitude_bare;
   TBranch *b_digi_charge_integrated_bare;
   TBranch *b_digi_charge_integrated_frac10;
   TBranch *b_digi_charge_integrated_frac30;
   TBranch *b_digi_charge_integrated_frac50;

   // Set object pointer
   digi_value_ch = 0;
   digi_value_time = 0;
   digi_value_bare_noise_sub = 0;

   BGOvalues = 0;
   SCINTvalues = 0;
   TDCreco = 0;
   digi_charge_integrated = 0;
   digi_max_amplitude = 0;
   digi_pedestal = 0;
   digi_pedestal_rms = 0;
   digi_time_at_frac30 = 0;
   digi_time_at_frac50_bare_noise_sub = 0;
   digi_time_at_1000_bare_noise_sub = 0;
   digi_time_at_max = 0;
   HODOX1 = 0;
   HODOX2 = 0;
   HODOY1 = 0;
   HODOY2 = 0;
   HODOSMALLvalues = 0;
   nTdcHits = 0;
   digi_charge_integrated_sub = 0;
   digi_charge_integrated_bare_noise_sub = 0;
   digi_charge_integrated_bare_noise_sub_fast = 0;
   digi_charge_integrated_bare_noise_sub_slow = 0;
   digi_max_amplitude_sub = 0;
   digi_max_amplitude_bare_noise_sub = 0;
   digi_pedestal_sub = 0;
   digi_pedestal_rms_sub = 0;
   digi_charge_integrated_corr1 = 0;
   digi_max_amplitude_corr1 = 0;
   digi_charge_integrated_corr2 = 0;
   digi_max_amplitude_corr2 = 0;

   digi_max_amplitude_bare = 0;
   digi_charge_integrated_bare = 0;
   digi_charge_integrated_frac10 = 0;
   digi_charge_integrated_frac30 = 0;
   digi_charge_integrated_frac50 = 0;
   

   // Set branch addresses and branch pointers
 fChain = tree;

   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("digi_frequency", &digi_frequency, &b_digi_frequency);
   fChain->SetBranchAddress("spillNumber", &spillNumber, &b_spillNumber);
   fChain->SetBranchAddress("evtNumber", &evtNumber, &b_evtNumber);
   fChain->SetBranchAddress("digi_value_ch", &digi_value_ch, &b_digi_value_ch);
   fChain->SetBranchAddress("digi_value_time", &digi_value_time, &b_digi_value_time);
   fChain->SetBranchAddress("digi_value_bare_noise_sub", &digi_value_bare_noise_sub, &b_digi_value_bare_noise_sub);
   fChain->SetBranchAddress("BGOvalues", &BGOvalues, &b_BGOvalues);
   fChain->SetBranchAddress("SCINTvalues", &SCINTvalues, &b_SCINTvalues);
   fChain->SetBranchAddress("TDCreco", &TDCreco, &b_TDCreco);
   fChain->SetBranchAddress("digi_charge_integrated", &digi_charge_integrated, &b_digi_charge_integrated);
   fChain->SetBranchAddress("digi_max_amplitude", &digi_max_amplitude, &b_digi_max_amplitude);
   fChain->SetBranchAddress("digi_pedestal", &digi_pedestal, &b_digi_pedestal);
   fChain->SetBranchAddress("digi_pedestal_rms", &digi_pedestal_rms, &b_digi_pedestal_rms);
   fChain->SetBranchAddress("digi_time_at_max", &digi_time_at_max, &b_digi_time_at_max);
   fChain->SetBranchAddress("digi_time_at_frac50_bare_noise_sub", &digi_time_at_frac50_bare_noise_sub, &b_digi_time_at_frac50_bare_noise_sub);
   fChain->SetBranchAddress("digi_time_at_1000_bare_noise_sub", &digi_time_at_1000_bare_noise_sub, &b_digi_time_at_1000_bare_noise_sub);
   fChain->SetBranchAddress("HODOX1", &HODOX1, &b_HODOX1);
   fChain->SetBranchAddress("HODOX2", &HODOX2, &b_HODOX2);
   fChain->SetBranchAddress("HODOY1", &HODOY1, &b_HODOY1);
   fChain->SetBranchAddress("HODOY2", &HODOY2, &b_HODOY2);
   fChain->SetBranchAddress("HODOSMALLvalues", &HODOSMALLvalues, &b_HODOSMALLvalues);
   fChain->SetBranchAddress("TableX", &TableX, &b_TableX);
   fChain->SetBranchAddress("TableY", &TableY, &b_TableY);
   fChain->SetBranchAddress("CeF3HV", &CeF3HV, &b_CeF3HV);
   fChain->SetBranchAddress("BGOHV", &BGOHV, &b_BGOHV);
   fChain->SetBranchAddress("BeamEnergy", &BeamEnergy, &b_BeamEnergy);
   fChain->SetBranchAddress("BeamTilt", &BeamTilt, &b_BeamTilt);
   fChain->SetBranchAddress("IsPhysics", &IsPhysics, &b_IsPhysics);
   fChain->SetBranchAddress("nTdcHits", &nTdcHits, &b_nTdcHits);
   fChain->SetBranchAddress("digi_pedestal_sub", &digi_pedestal_sub, &b_digi_pedestal_sub);
   fChain->SetBranchAddress("digi_pedestal_rms_sub", &digi_pedestal_rms_sub, &b_digi_pedestal_rms_sub);


   fChain->SetBranchAddress("digi_max_amplitude_bare", &digi_max_amplitude_bare, &b_digi_max_amplitude_bare);
   fChain->SetBranchAddress("digi_charge_integrated_bare", &digi_charge_integrated_bare, &b_digi_charge_integrated_bare);
   fChain->SetBranchAddress("digi_charge_integrated_bare_noise_sub", &digi_charge_integrated_bare_noise_sub, &b_digi_charge_integrated_bare_noise_sub);
   fChain->SetBranchAddress("digi_max_amplitude_bare_noise_sub", &digi_max_amplitude_bare_noise_sub, &b_digi_max_amplitude_bare_noise_sub);
   fChain->SetBranchAddress("digi_charge_integrated_bare_noise_sub_fast", &digi_charge_integrated_bare_noise_sub_fast, &b_digi_charge_integrated_bare_noise_sub_fast);
   fChain->SetBranchAddress("digi_charge_integrated_bare_noise_sub_slow", &digi_charge_integrated_bare_noise_sub_slow, &b_digi_charge_integrated_bare_noise_sub_slow);





  ///AND FINALLY THE RECO TREE/////////  
   std::string outdir = "analysisTrees_" + tag;
   system( Form("mkdir -p %s", outdir.c_str()) );

   std::string outfileName = outdir + "/Reco_" + runName + ".root";
   TFile* outfile = TFile::Open( outfileName.c_str(), "RECREATE" );

   TTree* outTree = new TTree("recoTree", "recoTree");

   TFile* waveformFile = TFile::Open("outWaveFormUtil_2778.root");


   //and now naming the stuff for the recoTree//
   outTree->Branch( "run", &runNumber, "run/i" );
   outTree->Branch( "spill", &spillNumber, "spill/i" );
   outTree->Branch( "event", &evtNumber, "event/i" );


   //Original
   std::vector<float> cef3( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3", &cef3 ); //for backwards compatibility
   std::vector<float> cef3_corr( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_corr", &cef3_corr );


   std::vector<float> cef3_maxAmpl( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_maxAmpl", &cef3_maxAmpl );
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


   int nentries = tree->GetEntries();
   nentries=1;

   RunHelper::getBeamPosition( runName, xBeam, yBeam );

   std::cout << nentries << std::endl;
 
   //get reference waveform for fits
   std::vector<TProfile*>  waveProfile;

   for (unsigned int iCh=0; iCh<CEF3_CHANNELS; iCh++) {
     TString name="prof";
     name+=iCh;
     TString istring;
     istring+=iCh;

     waveProfile.push_back(new TProfile(name,name,1024,0,0.4e-6));
     TH1F* histo = (TH1F*)waveformFile->Get("waveform_histo_"+istring);
     for(int i=0;i<=histo->GetNbinsX();++i) waveProfile.at(iCh)->Fill(histo->GetBinCenter(i), histo->GetBinContent(i));
     waveProfile.at(iCh)->Scale(1./waveProfile.at(iCh)->GetMaximum());
     waveProfile.at(iCh)->Print();
   }
 


   TCanvas dummycanvas;
   waveProfile.at(0)->Draw();
   dummycanvas.SaveAs("stocazzo.png");

  for(int  iEntry=0; iEntry<nentries; ++iEntry ) {
    
    tree->GetEntry( iEntry );     
    if( iEntry %  1000 == 0 ) std::cout << "Entry: " << iEntry << " / " << nentries << std::endl;

    //waveform creation
    std::vector<Waveform*> waveform;
    waveform.clear();
    for (unsigned int i=0; i<CEF3_CHANNELS; i++) {
      waveform.push_back(new Waveform());
      waveform.at(i)->clear();
    
    }
    

    float timeOfTheEvent=digi_time_at_1000_bare_noise_sub->at(8);//synchronizing time of events with time of trigger
    float shiftTime=190.3-timeOfTheEvent;//mean fitted on trigger run 2778
    int shiftSample=round(shiftTime/(1e9*timeSampleUnit(digi_frequency)));
    shiftSample=-shiftSample;
    for (int i=0;i<1024*4;++i){
      if(digi_value_ch->at(i) > 3)continue;
      if(digi_max_amplitude->at(digi_value_ch->at(i))>10000 || digi_max_amplitude->at(digi_value_ch->at(i))<0)continue;
      int iSample=i;
      if(i+shiftSample>1023*digi_value_ch->at(i) && i+shiftSample<(1023+(1024*digi_value_ch->at(i)))){
	iSample=i+shiftSample;
      }
      waveform.at(digi_value_ch->at(i))->addTimeAndSample((i-1024*digi_value_ch->at(i))*timeSampleUnit(digi_frequency),digi_value_bare_noise_sub->at(iSample));
      //      std::cout<<"i:"<<digi_value_ch->at(i)<<" time:"<<(i-1024*digi_value_ch->at(i))*timeSampleUnit(digi_frequency)<<" value"<<digi_value_bare_noise_sub->at(iSample)<<std::endl;
    }

    float finalFastSample=230;
    if(digi_frequency==1)finalFastSample=172;//no direct realtion between 5Gs and 2.5Gs since offset is different
    
    std::vector<float> charge_slow;
    std::vector<float> charge_fast;
    for (unsigned int i=0; i<CEF3_CHANNELS; i++) {      
      //      if(iEntry>100)std::cout<<shiftSample<<std::endl;
//      if(shiftSample<100 && shiftSample>-100)finalFastSample+=shiftSample;
//      else shiftSample = 0;
      charge_fast.push_back(waveform.at(i)->charge_integrated(4,finalFastSample));
      charge_slow.push_back(waveform.at(i)->charge_integrated(finalFastSample,900));
      
      //WaveformFit using mean Waveform
      ROOT::Math::Minimizer* minimizer;
      Waveform::max_amplitude_informations wave_max = waveform.at(i)->max_amplitude(50,900,5);
      Waveform::baseline_informations wave_pedestal = waveform.at(i)->baseline(5,34);
      WaveformFit::fitWaveform(waveform.at(i),waveProfile.at(i),200,200,wave_max,wave_pedestal,minimizer);
      exit(9);
    }




    std::string theBeamEnergy = Form("%.0f",BeamEnergy);
    if( runNumber > 272 && runNumber < 298){
      theBeamEnergy = "273"; //For the long position scan prior to tdc adjustment
    }else if( runNumber  < 273){
      theBeamEnergy = "259"; //For the first runs, before adjustments
    }
    //   std::cout << "The used constant file has label = "<< theBeamEnergy  << std::endl;
    
    
    //set the tag for calibration
    TagHelper tagHelper(tag,theBeamEnergy);
    EnergyCalibration cef3Calib(tagHelper.getCeF3FileName());
    EnergyCalibration bgoCalib(tagHelper.getBGOFileName());
    AlignmentOfficer alignOfficer(tagHelper.getAlignmentFileName());
    





     //BACKWARDS COMPATIBILITY///
     assignValues( cef3, *digi_charge_integrated_bare_noise_sub, CEF3_START_CHANNEL);      

     assignValues( cef3_corr, *digi_charge_integrated_bare_noise_sub, CEF3_START_CHANNEL );  
     cef3Calib.applyCalibration(cef3_corr);

     assignValues( cef3_maxAmpl, *digi_max_amplitude_bare_noise_sub, CEF3_START_CHANNEL);
     assignValues( cef3_chaInt, *digi_charge_integrated_bare_noise_sub, CEF3_START_CHANNEL);
     assignValues( cef3_chaInt_wls, charge_slow, CEF3_START_CHANNEL);
     assignValues( cef3_chaInt_cher, charge_fast, CEF3_START_CHANNEL);
     computeCherenkov(cef3_chaInt_cher,cef3_chaInt_wls);

     charge_slow.clear();
     charge_fast.clear();

     //     assignValues( bgo_corr, *BGOvalues, 0 );
     //     bgoCalib.applyCalibration(bgo_corr);

     std::vector<bool> hodoX1_values(HODOX1_CHANNELS, -1.);
     std::vector<bool> hodoY1_values(HODOY1_CHANNELS, -1.);
     assignValuesBool( hodoX1_values, *HODOX1, 0. );
     assignValuesBool( hodoY1_values, *HODOY1, 0. );


     std::vector<bool> hodoX2_values(HODOX2_CHANNELS, -1.);
     std::vector<bool> hodoY2_values(HODOY2_CHANNELS, -1.);
     assignValuesBool( hodoX2_values, *HODOX2, 0 );
     assignValuesBool( hodoY2_values, *HODOY2, 0 );


     std::vector<float> hodoSmallX_values(HODOSMALLX_CHANNELS, -1.);
     std::vector<float> hodoSmallY_values(HODOSMALLY_CHANNELS, -1.);
     assignValues( hodoSmallY_values, *HODOSMALLvalues, 0 );
     assignValues( hodoSmallX_values, *HODOSMALLvalues, 4 );



     // hodo cluster reconstruction
     int clusterMaxFibres = 4;
     doHodoReconstructionBool( hodoX1_values    , nClusters_hodoX1    , nFibres_hodoX1    , pos_hodoX1    , 0.5, clusterMaxFibres, 0. );
     doHodoReconstructionBool( hodoY1_values    , nClusters_hodoY1    , nFibres_hodoY1    , pos_hodoY1    , 0.5, clusterMaxFibres, 0. );
     doHodoReconstructionBool( hodoX2_values    , nClusters_hodoX2    , nFibres_hodoX2    , pos_hodoX2    , 0.5, clusterMaxFibres , 0.);
     doHodoReconstructionBool( hodoY2_values    , nClusters_hodoY2    , nFibres_hodoY2    , pos_hodoY2    , 0.5, clusterMaxFibres, 0. );

     std::vector<float> pedMeanX,pedMeanY, pedSigmaX, pedSigmaY,cutValuesX,cutValuesY;
     //values obtained fitting pedestals from run 2778
     pedMeanY.push_back(155.15);
     pedMeanY.push_back(141.30);
     pedMeanY.push_back(152.90);
     pedMeanY.push_back(152.13);

     pedMeanX.push_back(146.89);     
     pedMeanX.push_back(139.89);
     pedMeanX.push_back(151.73);
     pedMeanX.push_back(53.27);

     pedSigmaY.push_back(1.20);
     pedSigmaY.push_back(1.25);
     pedSigmaY.push_back(1.60);
     pedSigmaY.push_back(1.20);

     pedSigmaX.push_back(1.28);
     pedSigmaX.push_back(1.32);
     pedSigmaX.push_back(1.31);
     pedSigmaX.push_back(0.83);

     for(int i=0;i<pedMeanX.size();++i){
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


     xTable = TableX;
     yTable = TableY;
     HVCeF3 = CeF3HV;
     HVBGO = BGOHV;
     beamEnergy = BeamEnergy;
     angle = BeamTilt;

     nTDCHits[0] = nTdcHits->at(0);
     nTDCHits[1] = nTdcHits->at(1);
     nTDCHits[2] = nTdcHits->at(2);
     nTDCHits[3] = nTdcHits->at(3);

     wc_x = TDCreco->at(0);
     wc_y = TDCreco->at(1);

     if( runNumber>=170 ) wc_y = -wc_y;

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





     outTree->Fill();
    
   
   } // for entries

 
   outfile->cd();
   outTree->Write();
   outfile->Close();

   std::cout << "-> Analysis Tree saved in: " << outfile->GetName() << std::endl;
   return 0;

}
  


void assignValues( std::vector<float> &target, std::vector<float> source, unsigned int startPos ) {

  for( unsigned i=0; i<target.size(); ++i ) 
    target[i] = source[startPos+i];

}


void assignValuesBool( std::vector<bool> &target, std::vector<bool> source, unsigned int startPos ) {

  for( unsigned i=0; i<target.size(); ++i ) 
    target[i] = source[startPos+i];

}

void computeCherenkov(std::vector<float> &cher,std::vector<float> wls){
  TVectorD integrals(4);//values obtained from run2778
//sample 172
  integrals[0]=0.0787133;
  integrals[1]=0.0388933;
  integrals[2]=0.210458;
  integrals[3]=0.042154;

  for (int i=0;i<cher.size();++i){
    float chargeWlsUnderCher=(wls[i]*integrals[i])/(1-integrals[i]);//amount of wls charge under cherenkov peak. estimated by fit
    cher[i]=cher[i]-chargeWlsUnderCher;
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
