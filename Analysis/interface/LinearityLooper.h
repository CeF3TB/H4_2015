#ifndef LinearityLooper_h
#define LinearityLooper_h

#include <Event.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <cstdlib>

// Header file for the classes stored in the TTree if any.
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.
using namespace std;

class LinearityLooper: public Event{
 public:
  LinearityLooper(TTree *tree=0);
  void     Loop();
  ~LinearityLooper();
};

#endif

#ifdef LinearityLooper_cxx

LinearityLooper::LinearityLooper(TTree *tree) 
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


LinearityLooper::~LinearityLooper()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}


#endif // #ifdef LinearityLooper_cxx
