#ifndef CherenkovAnalyzer_h
#define CherenkovAnalyzer_h

#include <Event.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <cstdlib>
#include <TF1.h>
#include <TGraph.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.
using namespace std;

class CherenkovAnalyzer: public Event{
 public:
  TGraph* meanWaveGraphs[4];
  TH1F* meanWaveHistos[4];

  CherenkovAnalyzer(TTree *tree=0);
  void     Loop();
  ~CherenkovAnalyzer();
};

#endif

#ifdef CherenkovAnalyzer_cxx

CherenkovAnalyzer::CherenkovAnalyzer(TTree *tree) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  fChain=0;
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("output2015.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("output2015.root");
      }
      f->GetObject("outputTree",tree);

   }
   Init(tree);
}


CherenkovAnalyzer::~CherenkovAnalyzer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}


#endif // #ifdef CherenkovAnalyzer_cxx
