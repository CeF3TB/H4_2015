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

   TH1F* charge_integrated[nFibers];
   TH1F* charge_integrated_bare[nFibers];
   TH1F* charge_integrated_bare_noiseSub[nFibers];

   for(int i =0; i<nFibers;++i){
     TString fiber;
     fiber.Form("%d",i);
     max_amplitude[i]=new TH1F("max_amplitude_"+fiber,"max_amplitude_"+fiber,100,0,1000);  
     max_amplitude_bare[i]=new TH1F("max_amplitude_bare_"+fiber,"max_amplitude_bare_"+fiber,100,0,1000);  
     max_amplitude_bare_noiseSub[i]=new TH1F("max_amplitude_bare_noiseSub_"+fiber,"max_amplitude_bare_noiseSub_"+fiber,100,0,1000);  
     
     charge_integrated[i]=new TH1F("charge_integrated_"+fiber,"charge_integrated_"+fiber,100,0,200000);  
     charge_integrated_bare[i]=new TH1F("charge_integrated_bare_"+fiber,"charge_integrated_bare_"+fiber,100,0,200000);  
     charge_integrated_bare_noiseSub[i]=new TH1F("charge_integrated_bare_noiseSub_"+fiber,"charge_integrated_bare_noiseSub_"+fiber,100,0,200000);  
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

	charge_integrated[i]->Fill(digi_charge_integrated->at(i));
	charge_integrated_bare[i]->Fill(digi_charge_integrated_bare->at(i));
	charge_integrated_bare_noiseSub[i]->Fill(digi_charge_integrated_bare_noise_sub->at(i));

      }
   }

   TPaveText * pave= new TPaveText(0.8,0.7,0.95,0.85,"NDC");
   pave->SetFillColor(kWhite);
   pave->SetTextSize(0.030);
   pave->SetTextAlign(kHAlignRight);
   pave->SetTextFont(62);

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


     TString meanbare;
     meanbare.Form("%.1f",max_amplitude_bare[i]->GetMean());
     double rms=max_amplitude_bare[i]->GetRMS();
     TString rmsBare=Form("%.1f",rms);
     pave->AddText("Bare Mean:"+meanbare+" rms:"+rmsBare);


     TString meanbare_noiseSub;
     meanbare_noiseSub.Form("%.1f",max_amplitude_bare_noiseSub[i]->GetMean());
     rms=max_amplitude_bare_noiseSub[i]->GetRMS();
     TString rmsBare_NoiseSub=Form("%.1f",rms);
     pave->AddText("NoiseSub Mean:"+meanbare_noiseSub+" rms:"+rmsBare_NoiseSub);
     pave->SetAllWith("NoiseSub Mean:"+meanbare_noiseSub+" rms:"+rmsBare_NoiseSub,"color",kBlue);

     TString mean;
     mean.Form("%.1f",max_amplitude[i]->GetMean());
     rms=max_amplitude[i]->GetRMS();
     TString rms_fft=Form("%.1f",rms);
     pave->AddText("fft Mean:"+mean+" rms:"+rms_fft);
     pave->SetAllWith("fft Mean:"+mean+" rms:"+rms_fft,"color",kRed);

     max_amplitude_bare[i]->Draw();
     max_amplitude_bare_noiseSub[i]->Draw("same");
     max_amplitude[i]->Draw("same");
     pave->Draw("same");
     gPad->RedrawAxis();

     TString fiber;
     fiber.Form("%d",i);

     c1.SaveAs("plots/max_amplitude_"+fiber+".png");


     pave->Clear();
     c1.Clear();

     meanbare.Form("%3.f",charge_integrated_bare[i]->GetMean());
     rms=charge_integrated_bare[i]->GetRMS();
     rmsBare=Form("%3.f",rms);
     pave->AddText("Bare Mean:"+meanbare+" rms:"+rmsBare);


     meanbare_noiseSub.Form("%3.f",charge_integrated_bare_noiseSub[i]->GetMean());
     rms=charge_integrated_bare_noiseSub[i]->GetRMS();
     rmsBare_NoiseSub=Form("%3.f",rms);
     pave->AddText("NoiseSub Mean:"+meanbare_noiseSub+" rms:"+rmsBare_NoiseSub);
     pave->SetAllWith("NoiseSub Mean:"+meanbare_noiseSub+" rms:"+rmsBare_NoiseSub,"color",kBlue);

     mean.Form("%3.f",charge_integrated[i]->GetMean());
     rms=charge_integrated[i]->GetRMS();
     rms_fft=Form("%3.f",rms);
     pave->AddText("fft Mean:"+mean+" rms:"+rms_fft);
     pave->SetAllWith("fft Mean:"+mean+" rms:"+rms_fft,"color",kRed);



     max=charge_integrated_bare[i]->GetMaximum();
     if(charge_integrated_bare_noiseSub[i]->GetMaximum()>max){
       max=charge_integrated_bare_noiseSub[i]->GetMaximum();
     }
     if(charge_integrated[i]->GetMaximum()>max){
       max=charge_integrated[i]->GetMaximum();
     }
     max*=1.1;
     charge_integrated_bare[i]->GetYaxis()->SetRangeUser(0.,max);
     charge_integrated_bare[i]->SetLineColor(kBlack);
     charge_integrated_bare_noiseSub[i]->SetLineColor(kBlue);
     charge_integrated[i]->SetLineColor(kRed);

     charge_integrated_bare[i]->Draw();
     charge_integrated_bare_noiseSub[i]->Draw("same");
     charge_integrated[i]->Draw("same");
     pave->Draw("same");
     gPad->RedrawAxis();

     c1.SaveAs("plots/charge_integrated_"+fiber+".png");
     pave->Clear();


   }
   outFile->Write();
   outFile->Close();
}

