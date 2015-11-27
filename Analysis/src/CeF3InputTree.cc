#include "../interface/CeF3InputTree.h"

CeF3InputTree::CeF3InputTree(TTree* tree){

   fChain = tree;
   fChain->SetMakeClass(1);
   Init();
   SetBranches();
}



void CeF3InputTree::Init(){
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
   digi_time_at_max_noise_sub = 0;
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

}

void CeF3InputTree::SetBranches(){
   // Set branch addresses and branch pointers
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
   fChain->SetBranchAddress("digi_time_at_max_noise_sub", &digi_time_at_max_noise_sub, &b_digi_time_at_max_noise_sub);
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
}
