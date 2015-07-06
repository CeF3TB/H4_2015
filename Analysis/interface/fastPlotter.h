//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jul  3 17:24:29 2015 by ROOT version 5.34/18
// from TTree outputTree/outputTree
// found on file: output2015.root
//////////////////////////////////////////////////////////

#ifndef fastPlotter_h
#define fastPlotter_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <cstdlib>

// Header file for the classes stored in the TTree if any.
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.
using namespace std;

class fastPlotter {
public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  
   // Declaration of leaf types
   UInt_t          runNumber;
   UInt_t          spillNumber;
   UInt_t          evtNumber;
   vector<float>   *BGOvalues;
   vector<float>   *SCINTvalues;
   vector<float>   *HODOSMALLvalues;
   vector<float>   *TDCreco;
   vector<float>   *digi_charge_integrated;
   vector<float>   *digi_max_amplitude;
   vector<float>   *digi_pedestal;
   vector<float>   *digi_pedestal_rms;
   vector<float>   *digi_pedestal_noise_sub;
   vector<float>   *digi_pedestal_noise_sub_rms;
   vector<float>   *digi_pedestal_bare;
   vector<float>   *digi_pedestal_bare_rms;
   vector<float>   *digi_pedestal_bare_noise_sub;
   vector<float>   *digi_pedestal_bare_noise_sub_rms;
   vector<float>   *digi_time_at_max;
   vector<bool>    *HODOX1;
   vector<bool>    *HODOX2;
   vector<bool>    *HODOY1;
   vector<bool>    *HODOY2;
   Float_t         TableX;
   Float_t         TableY;
   Float_t         CeF3HV;
   Float_t         BGOHV;
   Float_t         BeamEnergy;
   Float_t         BeamTilt;
   Int_t           IsPhysics;
   vector<int>     *nTdcHits;
   vector<float>   *digi_max_amplitude_bare;
   vector<float>   *digi_time_at_max_bare;
   vector<float>   *digi_charge_integrated_bare;
   vector<float>   *digi_max_amplitude_bare_noise_sub;
   vector<float>   *digi_time_at_max_bare_noise_sub;
   vector<float>   *digi_charge_integrated_bare_noise_sub;
   vector<float>   *digi_max_amplitude_noise_sub;
   vector<float>   *digi_time_at_max_noise_sub;
   vector<float>   *digi_charge_integrated_noise_sub;
   vector<float>   *digi_value;
   vector<float>   *digi_value_bare;
   vector<int>     *digi_value_ch;
   vector<float>   *digi_value_time;
   vector<float>   *digi_value_noise_sub;
   vector<float>   *digi_value_bare_noise_sub;

   // List of branches
   TBranch        *b_runNumber;   //!
   TBranch        *b_spillNumber;   //!
   TBranch        *b_evtNumber;   //!
   TBranch        *b_BGOvalues;   //!
   TBranch        *b_SCINTvalues;   //!
   TBranch        *b_HODOSMALLvalues;   //!
   TBranch        *b_TDCreco;   //!
   TBranch        *b_digi_charge_integrated;   //!
   TBranch        *b_digi_max_amplitude;   //!
   TBranch        *b_digi_pedestal;   //!
   TBranch        *b_digi_pedestal_rms;   //!
   TBranch        *b_digi_pedestal_noise_sub;   //!
   TBranch        *b_digi_pedestal_noise_sub_rms;   //!
   TBranch        *b_digi_pedestal_bare;   //!
   TBranch        *b_digi_pedestal_bare_rms;   //!
   TBranch        *b_digi_pedestal_bare_noise_sub;   //!
   TBranch        *b_digi_pedestal_bare_noise_sub_rms;   //!
   TBranch        *b_digi_time_at_max;   //!
   TBranch        *b_HODOX1;   //!
   TBranch        *b_HODOX2;   //!
   TBranch        *b_HODOY1;   //!
   TBranch        *b_HODOY2;   //!
   TBranch        *b_TableX;   //!
   TBranch        *b_TableY;   //!
   TBranch        *b_CeF3HV;   //!
   TBranch        *b_BGOHV;   //!
   TBranch        *b_BeamEnergy;   //!
   TBranch        *b_BeamTilt;   //!
   TBranch        *b_IsPhysics;   //!
   TBranch        *b_nTdcHits;   //!
   TBranch        *b_digi_max_amplitude_bare;   //!
   TBranch        *b_digi_time_at_max_bare;   //!
   TBranch        *b_digi_charge_integrated_bare;   //!
   TBranch        *b_digi_max_amplitude_bare_noise_sub;   //!
   TBranch        *b_digi_time_at_max_bare_noise_sub;   //!
   TBranch        *b_digi_charge_integrated_bare_noise_sub;   //!
   TBranch        *b_digi_max_amplitude_noise_sub;   //!
   TBranch        *b_digi_time_at_max_noise_sub;   //!
   TBranch        *b_digi_charge_integrated_noise_sub;   //!
   TBranch        *b_digi_value;   //!
   TBranch        *b_digi_value_bare;   //!
   TBranch        *b_digi_value_ch;   //!
   TBranch        *b_digi_value_time;   //!
   TBranch        *b_digi_value_noise_sub;   //!
   TBranch        *b_digi_value_bare_noise_sub;   //!

   fastPlotter(TTree *tree=0);
   virtual ~fastPlotter();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef fastPlotter_cxx
fastPlotter::fastPlotter(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("output2015.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("output2015.root");
      }
      f->GetObject("outputTree",tree);

   }
   Init(tree);
}

fastPlotter::~fastPlotter()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t fastPlotter::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t fastPlotter::LoadTree(Long64_t entry)
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

void fastPlotter::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   BGOvalues = 0;
   SCINTvalues = 0;
   HODOSMALLvalues = 0;
   TDCreco = 0;
   digi_charge_integrated = 0;
   digi_max_amplitude = 0;
   digi_pedestal = 0;
   digi_pedestal_rms = 0;
   digi_pedestal_noise_sub = 0;
   digi_pedestal_noise_sub_rms = 0;
   digi_pedestal_bare = 0;
   digi_pedestal_bare_rms = 0;
   digi_pedestal_bare_noise_sub = 0;
   digi_pedestal_bare_noise_sub_rms = 0;
   digi_time_at_max = 0;
   HODOX1 = 0;
   HODOX2 = 0;
   HODOY1 = 0;
   HODOY2 = 0;
   nTdcHits = 0;
   digi_max_amplitude_bare = 0;
   digi_time_at_max_bare = 0;
   digi_charge_integrated_bare = 0;
   digi_max_amplitude_bare_noise_sub = 0;
   digi_time_at_max_bare_noise_sub = 0;
   digi_charge_integrated_bare_noise_sub = 0;
   digi_max_amplitude_noise_sub = 0;
   digi_time_at_max_noise_sub = 0;
   digi_charge_integrated_noise_sub = 0;
   digi_value = 0;
   digi_value_bare = 0;
   digi_value_ch = 0;
   digi_value_time = 0;
   digi_value_noise_sub = 0;
   digi_value_bare_noise_sub = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("spillNumber", &spillNumber, &b_spillNumber);
   fChain->SetBranchAddress("evtNumber", &evtNumber, &b_evtNumber);
   fChain->SetBranchAddress("BGOvalues", &BGOvalues, &b_BGOvalues);
   fChain->SetBranchAddress("SCINTvalues", &SCINTvalues, &b_SCINTvalues);
   fChain->SetBranchAddress("HODOSMALLvalues", &HODOSMALLvalues, &b_HODOSMALLvalues);
   fChain->SetBranchAddress("TDCreco", &TDCreco, &b_TDCreco);
   fChain->SetBranchAddress("digi_charge_integrated", &digi_charge_integrated, &b_digi_charge_integrated);
   fChain->SetBranchAddress("digi_max_amplitude", &digi_max_amplitude, &b_digi_max_amplitude);
   fChain->SetBranchAddress("digi_pedestal", &digi_pedestal, &b_digi_pedestal);
   fChain->SetBranchAddress("digi_pedestal_rms", &digi_pedestal_rms, &b_digi_pedestal_rms);
   fChain->SetBranchAddress("digi_pedestal_noise_sub", &digi_pedestal_noise_sub, &b_digi_pedestal_noise_sub);
   fChain->SetBranchAddress("digi_pedestal_noise_sub_rms", &digi_pedestal_noise_sub_rms, &b_digi_pedestal_noise_sub_rms);
   fChain->SetBranchAddress("digi_pedestal_bare", &digi_pedestal_bare, &b_digi_pedestal_bare);
   fChain->SetBranchAddress("digi_pedestal_bare_rms", &digi_pedestal_bare_rms, &b_digi_pedestal_bare_rms);
   fChain->SetBranchAddress("digi_pedestal_bare_noise_sub", &digi_pedestal_bare_noise_sub, &b_digi_pedestal_bare_noise_sub);
   fChain->SetBranchAddress("digi_pedestal_bare_noise_sub_rms", &digi_pedestal_bare_noise_sub_rms, &b_digi_pedestal_bare_noise_sub_rms);
   fChain->SetBranchAddress("digi_time_at_max", &digi_time_at_max, &b_digi_time_at_max);
   fChain->SetBranchAddress("HODOX1", &HODOX1, &b_HODOX1);
   fChain->SetBranchAddress("HODOX2", &HODOX2, &b_HODOX2);
   fChain->SetBranchAddress("HODOY1", &HODOY1, &b_HODOY1);
   fChain->SetBranchAddress("HODOY2", &HODOY2, &b_HODOY2);
   fChain->SetBranchAddress("TableX", &TableX, &b_TableX);
   fChain->SetBranchAddress("TableY", &TableY, &b_TableY);
   fChain->SetBranchAddress("CeF3HV", &CeF3HV, &b_CeF3HV);
   fChain->SetBranchAddress("BGOHV", &BGOHV, &b_BGOHV);
   fChain->SetBranchAddress("BeamEnergy", &BeamEnergy, &b_BeamEnergy);
   fChain->SetBranchAddress("BeamTilt", &BeamTilt, &b_BeamTilt);
   fChain->SetBranchAddress("IsPhysics", &IsPhysics, &b_IsPhysics);
   fChain->SetBranchAddress("nTdcHits", &nTdcHits, &b_nTdcHits);
   fChain->SetBranchAddress("digi_max_amplitude_bare", &digi_max_amplitude_bare, &b_digi_max_amplitude_bare);
   fChain->SetBranchAddress("digi_time_at_max_bare", &digi_time_at_max_bare, &b_digi_time_at_max_bare);
   fChain->SetBranchAddress("digi_charge_integrated_bare", &digi_charge_integrated_bare, &b_digi_charge_integrated_bare);
   fChain->SetBranchAddress("digi_max_amplitude_bare_noise_sub", &digi_max_amplitude_bare_noise_sub, &b_digi_max_amplitude_bare_noise_sub);
   fChain->SetBranchAddress("digi_time_at_max_bare_noise_sub", &digi_time_at_max_bare_noise_sub, &b_digi_time_at_max_bare_noise_sub);
   fChain->SetBranchAddress("digi_charge_integrated_bare_noise_sub", &digi_charge_integrated_bare_noise_sub, &b_digi_charge_integrated_bare_noise_sub);
   fChain->SetBranchAddress("digi_max_amplitude_noise_sub", &digi_max_amplitude_noise_sub, &b_digi_max_amplitude_noise_sub);
   fChain->SetBranchAddress("digi_time_at_max_noise_sub", &digi_time_at_max_noise_sub, &b_digi_time_at_max_noise_sub);
   fChain->SetBranchAddress("digi_charge_integrated_noise_sub", &digi_charge_integrated_noise_sub, &b_digi_charge_integrated_noise_sub);
   fChain->SetBranchAddress("digi_value", &digi_value, &b_digi_value);
   fChain->SetBranchAddress("digi_value_bare", &digi_value_bare, &b_digi_value_bare);
   fChain->SetBranchAddress("digi_value_ch", &digi_value_ch, &b_digi_value_ch);
   fChain->SetBranchAddress("digi_value_time", &digi_value_time, &b_digi_value_time);
   fChain->SetBranchAddress("digi_value_noise_sub", &digi_value_noise_sub, &b_digi_value_noise_sub);
   fChain->SetBranchAddress("digi_value_bare_noise_sub", &digi_value_bare_noise_sub, &b_digi_value_bare_noise_sub);
   Notify();
}

Bool_t fastPlotter::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void fastPlotter::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t fastPlotter::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef fastPlotter_cxx
