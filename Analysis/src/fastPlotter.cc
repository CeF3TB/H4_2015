#define fastPlotter_cxx
#include "fastPlotter.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#define nFibers 4
#define emptyChannelNo 4

void fastPlotter::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L fastPlotter.C
//      Root > fastPlotter t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();
   std::cout<<"nentries:"<<nentries<<std::endl;

   TFile* outFile=TFile::Open("outputFastPlotter.root","recreate");
   TH1F* max_amplitude[nFibers];
   TH1F* max_amplitude_bare[nFibers];
   TH1F* max_amplitude_bare_noiseSub[nFibers];

   for(int i =0; i<nFibers;++i){
     TString fiber;
     fiber.Form("%d",i);
     max_amplitude[i]=new TH1F("max_amplitude_"+fiber,"max_amplitude_"+fiber,100,0,1000);  
     max_amplitude_bare[i]=new TH1F("max_amplitude_bare_"+fiber,"max_amplitude_bare_"+fiber,100,0,1000);  
     max_amplitude_bare_noiseSub[i]=new TH1F("max_amplitude_bare_noiseSub_"+fiber,"max_amplitude_bare_noiseSub_"+fiber,100,0,1000);  
   }


   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      for(int i =0; i<nFibers;++i){
	max_amplitude[i]->Fill(digi_max_amplitude->at(i));
	max_amplitude_bare[i]->Fill(digi_max_amplitude_bare->at(i));
	max_amplitude_bare_noiseSub[i]->Fill(digi_max_amplitude_bare_noise_sub->at(i));
      }
   }



   for(int i =0; i<nFibers;++i){  
     TCanvas c1;

     float max=max_amplitude_bare[i]->GetMaximum();
     if(max_amplitude_bare_noiseSub[i]->GetMaximum()>max){
       max=max_amplitude_bare_noiseSub[i]->GetMaximum();
     }
     if(max_amplitude[i]->GetMaximum()>max){
       max=max_amplitude[i]->GetMaximum();
     }
     max*=1.1;
     max_amplitude_bare[i]->GetYaxis()->SetRangeUser(0.,max);
     max_amplitude_bare[i]->SetLineColor(kBlack);
     max_amplitude_bare_noiseSub[i]->SetLineColor(kBlue);
     max_amplitude[i]->SetLineColor(kRed);

     max_amplitude_bare[i]->Draw();
     max_amplitude_bare_noiseSub[i]->Draw("same");
     max_amplitude[i]->Draw("same");

     TString fiber;
     fiber.Form("%d",i);
     c1.SaveAs("plots/maxAmplitude_"+fiber+".png");

   }
   outFile->Write();
   outFile->Close();
}

