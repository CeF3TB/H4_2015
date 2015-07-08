#define LinearityLooper_cxx
#include "LinearityLooper.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>



void LinearityLooper::Loop(){

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    std::cout<<jentry<<" ";
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;      
  }
}
