#ifndef __WF_TREE__
#define __WF_TREE__

#include <string>
#include <vector>

#include "TTree.h"
#include "TString.h"

using namespace std;

typedef unsigned long int uint32;
typedef unsigned long long int uint64;
 
//****************************************************************************************

class WFTree
{
 public: 
  //---ctors---
  WFTree(int nCh, int nSamples);
  //---dtor---
  ~WFTree();

  //---utils---
  void Fill() {tree_->Fill();};
  void Write(const char* name="wf_tree", const char* title="wf_tree")
  {tree_->SetTitle(title); tree_->Write(name);};

  TTree* tree_; 
  TString suffix_;

  int*   index;
  int*   evtNumber;
  int    WF_samples;
  int*   WF_ch; 
  float* WF_time;
  float* WF_val;
  float  maxAmpl0;
  float  maxAmpl1;
  float  maxAmpl2;
  float  maxAmpl3;
  float  deltaT;
};

WFTree::WFTree(int nCh, int nSamples)
{
  tree_ = new TTree();

  //---set total number of WF samples 
  WF_samples = nSamples*nCh;
  WF_ch = new int[WF_samples];
  WF_time = new float[WF_samples];
  WF_val = new float[WF_samples];
  index = new int[WF_samples];
  evtNumber = new int[WF_samples];
  //---global branches
  tree_->Branch(Form("WF_samples%s",suffix_.Data()), &WF_samples, Form("WF_samples%s/I",suffix_.Data()));
  tree_->Branch(Form("WF_ch%s",suffix_.Data()), WF_ch,            Form("WF_ch%s[WF_samples%s]/I",suffix_.Data(),suffix_.Data()));
  tree_->Branch(Form("index%s",suffix_.Data()), index,            Form("index%s[WF_samples%s]/I",suffix_.Data(),suffix_.Data()));
  tree_->Branch(Form("evtNumber%s",suffix_.Data()), evtNumber,            Form("evtNumber%s[WF_samples%s]/I",suffix_.Data(),suffix_.Data()));
  tree_->Branch(Form("WF_time%s",suffix_.Data()), WF_time,        Form("WF_time%s[WF_samples%s]/F",suffix_.Data(),suffix_.Data()));
  tree_->Branch(Form("WF_val%s",suffix_.Data()), WF_val,          Form("WF_val%s[WF_samples%s]/F",suffix_.Data(),suffix_.Data()));
  tree_->Branch(Form("maxAmpl0%s",suffix_.Data()), &maxAmpl0, Form("maxAmpl0%s/F",suffix_.Data()));
  tree_->Branch(Form("maxAmpl1%s",suffix_.Data()), &maxAmpl1, Form("maxAmpl1%s/F",suffix_.Data()));
  tree_->Branch(Form("maxAmpl2%s",suffix_.Data()), &maxAmpl2, Form("maxAmpl2%s/F",suffix_.Data()));
  tree_->Branch(Form("maxAmpl3%s",suffix_.Data()), &maxAmpl3, Form("maxAmpl3%s/F",suffix_.Data()));
  tree_->Branch(Form("deltaT%s",suffix_.Data()), &deltaT, Form("deltaT%s/F",suffix_.Data()));
}

WFTree::~WFTree()
{
  delete[] WF_ch; 
  delete[] WF_time;
  delete[] WF_val;

  delete tree_;
}

#endif
