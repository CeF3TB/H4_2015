#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>


#include "DrawTools.h"
#include "RecoTree.h"
#include "channelInfo.h"

#include "TVectorD.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TF1.h"
#include "TLegend.h"

int main( int argc, char* argv[] ) {


  DrawTools::setStyle();
  std::vector<float> energies;
  std::vector<float> reso;
  std::vector<float> resoErr;
  std::vector<int> runs;

  energies.push_back(20);
  energies.push_back(50);
  energies.push_back(100);
  energies.push_back(150);
  energies.push_back(200);
  energies.push_back(250);

  //HV 1450
  runs.push_back(2851);
  runs.push_back(2801);
  runs.push_back(2778);
  runs.push_back(2822);
  runs.push_back(2872);
  runs.push_back(2894);


//HV 1050
//  runs.push_back(2864);
//  runs.push_back(2814);
//  runs.push_back(2795);
//  runs.push_back(2849);
//  runs.push_back(2885);
//  runs.push_back(2907);


  TGraphErrors* gr_reso_vs_energy[5];
  for(int i=0;i<5;++i) gr_reso_vs_energy[i]= new TGraphErrors(0);


  for (int i=0;i<runs.size();++i){
    TString run;
    run.Form("%d",runs[i]); 
    
    
    TFile *inputFile=TFile::Open("CherenkovPlots_"+run+".root");
    TVectorD* res=(TVectorD*)inputFile->Get("resValue");
    TVectorD* resErr=(TVectorD*)inputFile->Get("resErrValue");

    for(int j=0;j<5;++j){
      gr_reso_vs_energy[j]->SetPoint( i, energies[i], (*res)[j] );
      gr_reso_vs_energy[j]->SetPointError( i, 0, (*resErr)[j] );
    }
  }

  for(int i=0;i<5;++i){
    TString fibre;
    fibre.Form("%d",i); 

    TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
    c1->cd();
    TH2D* h2_axes2 = new TH2D( "axes", "", 100, -0.0, 305. , 10, 0.0, 25 );
    h2_axes2->SetXTitle("Beam Energy [GeV]");
    h2_axes2->SetYTitle("Energy Resolution [%]");
    h2_axes2->Draw("");

    gr_reso_vs_energy[i]->SetMarkerStyle(20);
    gr_reso_vs_energy[i]->SetMarkerSize(1.6);
    gr_reso_vs_energy[i]->SetMarkerColor(kBlue);
    gr_reso_vs_energy[i]->Draw("p same");

    //        TF1 *fun= new TF1("fun","sqrt([0]*[0]/(x*x)+[1]*[1])",1, 250.+15.);
    //TF1 *fun= new TF1("fun","sqrt([0]*[0]/x+[1]*[1])",1, 250.+15.);
    TF1 *fun= new TF1("fun",  "sqrt([0]*[0]/x+[1]*[1]+ [2]*[2]/(x*x))",1, 250+5.);
    fun->SetParameter(1, 4.);
    fun->SetParameter(0, 20.);
    gr_reso_vs_energy[i]->Fit(fun,"RN");
    fun->SetLineWidth(1.);
    fun->SetLineColor(kBlue+2);
    fun->Draw("L same");

    TLegend* leg_neat = new TLegend(0.42, 0.92-0.06*5 , 0.9, 0.92);
    leg_neat->SetTextSize(0.038);
    string ene="Data fibre ";
    //  ene+=Form("%.0f",energies[i]);
    //    ene+=" GeV";
    ene+=fibre;
    leg_neat->AddEntry(gr_reso_vs_energy[i],ene.c_str(),"p");
    leg_neat->AddEntry((TObject*)0 ,Form("S =  %.2f\n%s / #sqrt{E [GeV]}",fun->GetParameter(0),"%" ),"");
    leg_neat->AddEntry( (TObject*)0 ,Form("C =  %.2f\n%s",(fun->GetParameter(1)) ,"%" ),"");
    leg_neat->SetFillColor(0);
    leg_neat->Draw("same");
    std::cout<<"parameters: S="<<fun->GetParameter(0)<<" C="<<fun->GetParameter(1)<<" N="<<fun->GetParameter(2)<<std::endl;
    c1->SaveAs("plots_reso/reso_HV1450_"+fibre+".png");
  }
  return 0;
}
