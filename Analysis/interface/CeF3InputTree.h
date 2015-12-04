#ifndef Cef3InputTree_h
#define  Cef3InputTree_h

#include "TROOT.h"
#include <vector> 
#include <algorithm> 
#include <iostream> 
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"

class CeF3InputTree
{
 public:


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
   std::vector<float>   *digi_time_at_max_noise_sub;
   std::vector<float>   *digi_time_at_max_bare_noise_sub;
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
   TBranch        *b_digi_time_at_max_noise_sub;   //!
   TBranch        *b_digi_time_at_max_bare_noise_sub;   //!
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


   CeF3InputTree( TTree* tree );
   ~CeF3InputTree();
   void Init();
   void SetBranches();

};

#endif
