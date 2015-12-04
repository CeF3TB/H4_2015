//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Sep 25 17:00:23 2015 by ROOT version 5.34/18
// from TTree recoTree/recoTree
// found on file: analysisTrees_V00/Reco_3076.root
//////////////////////////////////////////////////////////

#ifndef RecoTree_h
#define RecoTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class RecoTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          spill;
   UInt_t          event;
   std::vector<float>   *cef3;
   std::vector<float>   *cef3_corr;
   std::vector<float>   *cef3_maxAmpl;
   std::vector<float>   *cef3_maxAmpl_cher;
   std::vector<float>   *cef3_maxAmpl_wls;
   std::vector<float>   *cef3_maxAmpl_time;
   std::vector<float>   *cef3_maxAmpl_fit;
   std::vector<float>   *cef3_maxAmpl_fit_cher;
   std::vector<float>   *cef3_maxAmpl_fit_cher_status;
   std::vector<float>   *cef3_maxAmpl_fit_time_cher1;
   std::vector<float>   *cef3_maxAmpl_fit_time_cher2;
   std::vector<float>   *cef3_maxAmpl_fit_time_cher3;
   std::vector<float>   *cef3_maxAmpl_fit_corr;
   std::vector<float>   *cef3_chaInt;
   std::vector<float>   *cef3_chaInt_cher;
   std::vector<float>   *cef3_chaInt_wls;
   std::vector<float>   *bgo;
   std::vector<float>   *cef3_maxAmpl_corr;
   std::vector<float>   *cef3_chaInt_corr;
   std::vector<float>   *bgo_corr;
   Float_t         xTable;
   Float_t         yTable;
   Float_t         beamEnergy;
   Float_t         angle;
   Float_t         HVCeF3;
   Float_t         HVBGO;
   Float_t         xBeam;
   Float_t         yBeam;
   Int_t           nClusters_hodoX1;
   Int_t           nFibres_hodoX1[30];   //[nClusters_hodoX1]
   Float_t         pos_hodoX1[30];   //[nClusters_hodoX1]
   Float_t         pos_corr_hodoX1[30];   //[nClusters_hodoX1]
   Int_t           nClusters_hodoY1;
   Int_t           nFibres_hodoY1[30];   //[nClusters_hodoY1]
   Float_t         pos_hodoY1[30];   //[nClusters_hodoY1]
   Float_t         pos_corr_hodoY1[30];   //[nClusters_hodoY1]
   Int_t           nClusters_hodoX2;
   Int_t           nFibres_hodoX2[30];   //[nClusters_hodoX2]
   Float_t         pos_hodoX2[30];   //[nClusters_hodoX2]
   Float_t         pos_corr_hodoX2[30];   //[nClusters_hodoX2]
   Int_t           nClusters_hodoY2;
   Int_t           nFibres_hodoY2[30];   //[nClusters_hodoY2]
   Float_t         pos_hodoY2[30];   //[nClusters_hodoY2]
   Float_t         pos_corr_hodoY2[30];   //[nClusters_hodoY2]
   Int_t           nClusters_hodoSmallX;
   Int_t           nFibres_hodoSmallX[4];   //[nClusters_hodoSmallX]
   Float_t         pos_hodoSmallX[4];   //[nClusters_hodoSmallX]
   Float_t         pos_corr_hodoSmallX[4];   //[nClusters_hodoSmallX]
   Int_t           nClusters_hodoSmallY;
   Int_t           nFibres_hodoSmallY[4];   //[nClusters_hodoSmallY]
   Float_t         pos_hodoSmallY[4];   //[nClusters_hodoSmallY]
   Float_t         pos_corr_hodoSmallY[4];   //[nClusters_hodoSmallY]
   std::vector<int>     *nTDCHits;
   Float_t         pos_2FibClust_hodoX1;
   Float_t         pos_2FibClust_corr_hodoX1;
   Float_t         pos_2FibClust_hodoY1;
   Float_t         pos_2FibClust_corr_hodoY1;
   Float_t         pos_2FibClust_hodoX2;
   Float_t         pos_2FibClust_corr_hodoX2;
   Float_t         pos_2FibClust_hodoY2;
   Float_t         pos_2FibClust_corr_hodoY2;
   Float_t         cluster_pos_hodoX1;
   Float_t         cluster_pos_corr_hodoX1;
   Float_t         cluster_pos_hodoX2;
   Float_t         cluster_pos_corr_hodoX2;
   Float_t         cluster_pos_hodoY1;
   Float_t         cluster_pos_corr_hodoY1;
   Float_t         cluster_pos_hodoY2;
   Float_t         cluster_pos_corr_hodoY2;
   Float_t         wc_x;
   Float_t         wc_y;
   Float_t         wc_x_corr;
   Float_t         wc_y_corr;
   Float_t         mcp_time_frac50;
   Float_t         mcp_time_at_150;
   Float_t         mcp_max_amplitude;
   Float_t         nino_LEtime;
   Float_t         nino_LEchi2;
   Float_t         nino_maxAmpl;





   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_spill;   //!
   TBranch        *b_event;   //!
   TBranch        *b_cef3;   //!
   TBranch        *b_cef3_corr;   //!
   TBranch        *b_cef3_maxAmpl;   //!
   TBranch        *b_cef3_maxAmpl_cher;   //!
   TBranch        *b_cef3_maxAmpl_wls;   //!
   TBranch        *b_cef3_maxAmpl_time;   //!
   TBranch        *b_cef3_maxAmpl_fit;   //!
   TBranch        *b_cef3_maxAmpl_fit_cher;   //!
   TBranch        *b_cef3_maxAmpl_fit_cher_status;   //!
   TBranch        *b_cef3_maxAmpl_fit_time_cher1;   //!
   TBranch        *b_cef3_maxAmpl_fit_time_cher2;   //!
   TBranch        *b_cef3_maxAmpl_fit_time_cher3;   //!
   TBranch        *b_cef3_maxAmpl_fit_corr;   //!
   TBranch        *b_cef3_chaInt;   //!
   TBranch        *b_cef3_chaInt_cher;   //!
   TBranch        *b_cef3_chaInt_wls;   //!
   TBranch        *b_bgo;   //!
   TBranch        *b_cef3_maxAmpl_corr;   //!
   TBranch        *b_cef3_chaInt_corr;   //!
   TBranch        *b_bgo_corr;   //!
   TBranch        *b_xTable;   //!
   TBranch        *b_yTable;   //!
   TBranch        *b_xbeamEnergy;   //!
   TBranch        *b_angle;   //!
   TBranch        *b_HVCeF3;   //!
   TBranch        *b_HVBGO;   //!
   TBranch        *b_xBeam;   //!
   TBranch        *b_yBeam;   //!
   TBranch        *b_nClusters_hodoX1;   //!
   TBranch        *b_nFibres_hodoX1;   //!
   TBranch        *b_pos_hodoX1;   //!
   TBranch        *b_pos_corr_hodoX1;   //!
   TBranch        *b_nClusters_hodoY1;   //!
   TBranch        *b_nFibres_hodoY1;   //!
   TBranch        *b_pos_hodoY1;   //!
   TBranch        *b_pos_corr_hodoY1;   //!
   TBranch        *b_nClusters_hodoX2;   //!
   TBranch        *b_nFibres_hodoX2;   //!
   TBranch        *b_pos_hodoX2;   //!
   TBranch        *b_pos_corr_hodoX2;   //!
   TBranch        *b_nClusters_hodoY2;   //!
   TBranch        *b_nFibres_hodoY2;   //!
   TBranch        *b_pos_hodoY2;   //!
   TBranch        *b_pos_corr_hodoY2;   //!
   TBranch        *b_nClusters_hodoSmallX;   //!
   TBranch        *b_nFibres_hodoSmallX;   //!
   TBranch        *b_pos_hodoSmallX;   //!
   TBranch        *b_pos_corr_hodoSmallX;   //!
   TBranch        *b_nClusters_hodoSmallY;   //!
   TBranch        *b_nFibres_hodoSmallY;   //!
   TBranch        *b_pos_hodoSmallY;   //!
   TBranch        *b_pos_corr_hodoSmallY;   //!
   TBranch        *b_nTDCHits;   //!
   TBranch        *b_pos_2FibClust_hodoX1;   //!
   TBranch        *b_pos_2FibClust_corr_hodoX1;   //!
   TBranch        *b_pos_2FibClust_hodoY1;   //!
   TBranch        *b_pos_2FibClust_corr_hodoY1;   //!
   TBranch        *b_pos_2FibClust_hodoX2;   //!
   TBranch        *b_pos_2FibClust_corr_hodoX2;   //!
   TBranch        *b_pos_2FibClust_hodoY2;   //!
   TBranch        *b_pos_2FibClust_corr_hodoY2;   //!
   TBranch        *b_cluster_pos_hodoX1;   //!
   TBranch        *b_cluster_pos_corr_hodoX1;   //!
   TBranch        *b_cluster_pos_hodoX2;   //!
   TBranch        *b_cluster_pos_corr_hodoX2;   //!
   TBranch        *b_cluster_pos_hodoY1;   //!
   TBranch        *b_cluster_pos_corr_hodoY1;   //!
   TBranch        *b_cluster_pos_hodoY2;   //!
   TBranch        *b_cluster_pos_corr_hodoY2;   //!
   TBranch        *b_wc_x;   //!
   TBranch        *b_wc_y;   //!
   TBranch        *b_wc_x_corr;   //!
   TBranch        *b_wc_y_corr;   //!
   TBranch        *b_mcp_time_frac50;   //!
   TBranch        *b_mcp_time_at_150;   //!
   TBranch        *b_mcp_max_amplitude;
   TBranch        *b_nino_LEtime;
   TBranch        *b_nino_LEchi2;
   TBranch        *b_nino_maxAmpl;



   RecoTree(TTree *tree=0);
   virtual ~RecoTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef RecoTree_cxx
RecoTree::RecoTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("analysisTrees_V00/Reco_3076.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("analysisTrees_V00/Reco_3076.root");
      }
      f->GetObject("recoTree",tree);

   }
   Init(tree);
}

RecoTree::~RecoTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t RecoTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t RecoTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void RecoTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   cef3 = 0;
   cef3_corr = 0;
   cef3_maxAmpl = 0;
   cef3_maxAmpl_cher = 0;
   cef3_maxAmpl_wls = 0;
   cef3_maxAmpl_time = 0;
   cef3_maxAmpl_fit = 0;
   cef3_maxAmpl_fit_cher = 0;
   cef3_maxAmpl_fit_cher_status = 0;
   cef3_maxAmpl_fit_time_cher1 = 0;
   cef3_maxAmpl_fit_time_cher2 = 0;
   cef3_maxAmpl_fit_time_cher3 = 0;
   cef3_maxAmpl_fit_corr = 0;
   cef3_chaInt = 0;
   cef3_chaInt_cher = 0;
   cef3_chaInt_wls = 0;
   bgo = 0;
   cef3_maxAmpl_corr = 0;
   cef3_chaInt_corr = 0;
   bgo_corr = 0;
   nTDCHits = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("spill", &spill, &b_spill);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("cef3", &cef3, &b_cef3);
   fChain->SetBranchAddress("cef3_corr", &cef3_corr, &b_cef3_corr);
   fChain->SetBranchAddress("cef3_maxAmpl", &cef3_maxAmpl, &b_cef3_maxAmpl);
   fChain->SetBranchAddress("cef3_maxAmpl_cher", &cef3_maxAmpl_cher, &b_cef3_maxAmpl_cher);
   fChain->SetBranchAddress("cef3_maxAmpl_wls", &cef3_maxAmpl_wls, &b_cef3_maxAmpl_wls);
   fChain->SetBranchAddress("cef3_maxAmpl_time", &cef3_maxAmpl_time, &b_cef3_maxAmpl_time);
   fChain->SetBranchAddress("cef3_maxAmpl_fit", &cef3_maxAmpl_fit, &b_cef3_maxAmpl_fit);
   fChain->SetBranchAddress("cef3_maxAmpl_fit_cher", &cef3_maxAmpl_fit_cher, &b_cef3_maxAmpl_fit_cher);
   fChain->SetBranchAddress("cef3_maxAmpl_fit_cher_status", &cef3_maxAmpl_fit_cher_status, &b_cef3_maxAmpl_fit_cher_status);
   fChain->SetBranchAddress("cef3_maxAmpl_fit_time_cher1", &cef3_maxAmpl_fit_time_cher1, &b_cef3_maxAmpl_fit_time_cher1);
   fChain->SetBranchAddress("cef3_maxAmpl_fit_time_cher2", &cef3_maxAmpl_fit_time_cher2, &b_cef3_maxAmpl_fit_time_cher2);
   fChain->SetBranchAddress("cef3_maxAmpl_fit_time_cher3", &cef3_maxAmpl_fit_time_cher3, &b_cef3_maxAmpl_fit_time_cher3);
   fChain->SetBranchAddress("cef3_maxAmpl_fit_corr", &cef3_maxAmpl_fit_corr, &b_cef3_maxAmpl_fit_corr);
   fChain->SetBranchAddress("cef3_chaInt", &cef3_chaInt, &b_cef3_chaInt);
   fChain->SetBranchAddress("cef3_chaInt_cher", &cef3_chaInt_cher, &b_cef3_chaInt_cher);
   fChain->SetBranchAddress("cef3_chaInt_wls", &cef3_chaInt_wls, &b_cef3_chaInt_wls);
   fChain->SetBranchAddress("bgo", &bgo, &b_bgo);
   fChain->SetBranchAddress("cef3_maxAmpl_corr", &cef3_maxAmpl_corr, &b_cef3_maxAmpl_corr);
   fChain->SetBranchAddress("cef3_chaInt_corr", &cef3_chaInt_corr, &b_cef3_chaInt_corr);
   fChain->SetBranchAddress("bgo_corr", &bgo_corr, &b_bgo_corr);
   fChain->SetBranchAddress("xTable", &xTable, &b_xTable);
   fChain->SetBranchAddress("yTable", &yTable, &b_yTable);
   fChain->SetBranchAddress("beamEnergy", &beamEnergy, &b_xbeamEnergy);
   fChain->SetBranchAddress("angle", &angle, &b_angle);
   fChain->SetBranchAddress("HVCeF3", &HVCeF3, &b_HVCeF3);
   fChain->SetBranchAddress("HVBGO", &HVBGO, &b_HVBGO);
   fChain->SetBranchAddress("xBeam", &xBeam, &b_xBeam);
   fChain->SetBranchAddress("yBeam", &yBeam, &b_yBeam);
   fChain->SetBranchAddress("nClusters_hodoX1", &nClusters_hodoX1, &b_nClusters_hodoX1);
   fChain->SetBranchAddress("nFibres_hodoX1", nFibres_hodoX1, &b_nFibres_hodoX1);
   fChain->SetBranchAddress("pos_hodoX1", pos_hodoX1, &b_pos_hodoX1);
   fChain->SetBranchAddress("pos_corr_hodoX1", pos_corr_hodoX1, &b_pos_corr_hodoX1);
   fChain->SetBranchAddress("nClusters_hodoY1", &nClusters_hodoY1, &b_nClusters_hodoY1);
   fChain->SetBranchAddress("nFibres_hodoY1", nFibres_hodoY1, &b_nFibres_hodoY1);
   fChain->SetBranchAddress("pos_hodoY1", pos_hodoY1, &b_pos_hodoY1);
   fChain->SetBranchAddress("pos_corr_hodoY1", pos_corr_hodoY1, &b_pos_corr_hodoY1);
   fChain->SetBranchAddress("nClusters_hodoX2", &nClusters_hodoX2, &b_nClusters_hodoX2);
   fChain->SetBranchAddress("nFibres_hodoX2", nFibres_hodoX2, &b_nFibres_hodoX2);
   fChain->SetBranchAddress("pos_hodoX2", pos_hodoX2, &b_pos_hodoX2);
   fChain->SetBranchAddress("pos_corr_hodoX2", pos_corr_hodoX2, &b_pos_corr_hodoX2);
   fChain->SetBranchAddress("nClusters_hodoY2", &nClusters_hodoY2, &b_nClusters_hodoY2);
   fChain->SetBranchAddress("nFibres_hodoY2", nFibres_hodoY2, &b_nFibres_hodoY2);
   fChain->SetBranchAddress("pos_hodoY2", pos_hodoY2, &b_pos_hodoY2);
   fChain->SetBranchAddress("pos_corr_hodoY2", pos_corr_hodoY2, &b_pos_corr_hodoY2);
   fChain->SetBranchAddress("nClusters_hodoSmallX", &nClusters_hodoSmallX, &b_nClusters_hodoSmallX);
   fChain->SetBranchAddress("nFibres_hodoSmallX", nFibres_hodoSmallX, &b_nFibres_hodoSmallX);
   fChain->SetBranchAddress("pos_hodoSmallX", pos_hodoSmallX, &b_pos_hodoSmallX);
   fChain->SetBranchAddress("pos_corr_hodoSmallX", pos_corr_hodoSmallX, &b_pos_corr_hodoSmallX);
   fChain->SetBranchAddress("nClusters_hodoSmallY", &nClusters_hodoSmallY, &b_nClusters_hodoSmallY);
   fChain->SetBranchAddress("nFibres_hodoSmallY", nFibres_hodoSmallY, &b_nFibres_hodoSmallY);
   fChain->SetBranchAddress("pos_hodoSmallY", pos_hodoSmallY, &b_pos_hodoSmallY);
   fChain->SetBranchAddress("pos_corr_hodoSmallY", pos_corr_hodoSmallY, &b_pos_corr_hodoSmallY);
   fChain->SetBranchAddress("nTDCHits", &nTDCHits, &b_nTDCHits);
   fChain->SetBranchAddress("pos_2FibClust_hodoX1", &pos_2FibClust_hodoX1, &b_pos_2FibClust_hodoX1);
   fChain->SetBranchAddress("pos_2FibClust_corr_hodoX1", &pos_2FibClust_corr_hodoX1, &b_pos_2FibClust_corr_hodoX1);
   fChain->SetBranchAddress("pos_2FibClust_hodoY1", &pos_2FibClust_hodoY1, &b_pos_2FibClust_hodoY1);
   fChain->SetBranchAddress("pos_2FibClust_corr_hodoY1", &pos_2FibClust_corr_hodoY1, &b_pos_2FibClust_corr_hodoY1);
   fChain->SetBranchAddress("pos_2FibClust_hodoX2", &pos_2FibClust_hodoX2, &b_pos_2FibClust_hodoX2);
   fChain->SetBranchAddress("pos_2FibClust_corr_hodoX2", &pos_2FibClust_corr_hodoX2, &b_pos_2FibClust_corr_hodoX2);
   fChain->SetBranchAddress("pos_2FibClust_hodoY2", &pos_2FibClust_hodoY2, &b_pos_2FibClust_hodoY2);
   fChain->SetBranchAddress("pos_2FibClust_corr_hodoY2", &pos_2FibClust_corr_hodoY2, &b_pos_2FibClust_corr_hodoY2);
   fChain->SetBranchAddress("cluster_pos_hodoX1", &cluster_pos_hodoX1, &b_cluster_pos_hodoX1);
   fChain->SetBranchAddress("cluster_pos_corr_hodoX1", &cluster_pos_corr_hodoX1, &b_cluster_pos_corr_hodoX1);
   fChain->SetBranchAddress("cluster_pos_hodoX2", &cluster_pos_hodoX2, &b_cluster_pos_hodoX2);
   fChain->SetBranchAddress("cluster_pos_corr_hodoX2", &cluster_pos_corr_hodoX2, &b_cluster_pos_corr_hodoX2);
   fChain->SetBranchAddress("cluster_pos_hodoY1", &cluster_pos_hodoY1, &b_cluster_pos_hodoY1);
   fChain->SetBranchAddress("cluster_pos_corr_hodoY1", &cluster_pos_corr_hodoY1, &b_cluster_pos_corr_hodoY1);
   fChain->SetBranchAddress("cluster_pos_hodoY2", &cluster_pos_hodoY2, &b_cluster_pos_hodoY2);
   fChain->SetBranchAddress("cluster_pos_corr_hodoY2", &cluster_pos_corr_hodoY2, &b_cluster_pos_corr_hodoY2);
   fChain->SetBranchAddress("wc_x", &wc_x, &b_wc_x);
   fChain->SetBranchAddress("wc_y", &wc_y, &b_wc_y);
   fChain->SetBranchAddress("wc_x_corr", &wc_x_corr, &b_wc_x_corr);
   fChain->SetBranchAddress("wc_y_corr", &wc_y_corr, &b_wc_y_corr);
   fChain->SetBranchAddress("mcp_time_frac50", &mcp_time_frac50, &b_mcp_time_frac50);
   fChain->SetBranchAddress("mcp_time_at_150", &mcp_time_at_150, &b_mcp_time_at_150);
   fChain->SetBranchAddress( "mcp_max_amplitude", &mcp_max_amplitude, &b_mcp_max_amplitude);
   fChain->SetBranchAddress( "nino_LEtime", &nino_LEtime, &b_nino_LEtime);
   fChain->SetBranchAddress( "nino_LEchi2", &nino_LEchi2, &b_nino_LEchi2);
   fChain->SetBranchAddress( "nino_maxAmpl", &nino_maxAmpl, &b_nino_maxAmpl);

   Notify();
}

Bool_t RecoTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void RecoTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t RecoTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef RecoTree_cxx
