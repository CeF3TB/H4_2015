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
  std::vector<int> runs1450;
  std::vector<int> runs1250;
  std::vector<int> runs1050;
  std::vector<int> runs850;


  energies.push_back(20);
  energies.push_back(50);
  energies.push_back(100);
  energies.push_back(150);
  energies.push_back(200);
  energies.push_back(250);

  runs1450.push_back(2851);
  runs1450.push_back(2801);
  runs1450.push_back(2778);
  runs1450.push_back(2822);
  runs1450.push_back(2872);
  runs1450.push_back(2894);

  runs1250.push_back(2862);
  runs1250.push_back(2812);
  runs1250.push_back(2791);
  runs1250.push_back(2845);
  runs1250.push_back(2883);
  runs1250.push_back(2905);

  runs1050.push_back(2864);
  runs1050.push_back(2814);
  runs1050.push_back(2793);
  runs1050.push_back(2847);
  runs1050.push_back(2885);
  runs1050.push_back(2907);

  runs850.push_back(2864+2);
  runs850.push_back(2814+2);
  runs850.push_back(2793+2);
  runs850.push_back(2847+2);
  runs850.push_back(2885+2);
  runs850.push_back(2907+2);


  TGraphErrors* gr_reso_vs_energy1450[5];
  TGraphErrors* gr_reso_vs_energy1250[5];
  TGraphErrors* gr_reso_vs_energy1050[5];
  TGraphErrors* gr_reso_vs_energy850[5];
  for(int i=0;i<5;++i){
      gr_reso_vs_energy1450[i]= new TGraphErrors(0);
      gr_reso_vs_energy1250[i]= new TGraphErrors(0);
      gr_reso_vs_energy1050[i]= new TGraphErrors(0);
      gr_reso_vs_energy850[i]= new TGraphErrors(0);
  }

  //HV 1450
  float max1450[5]={0};
  for (int i=0;i<runs1450.size();++i){
    TString run;
    run.Form("%d",runs1450[i]); 
    TFile *inputFile=TFile::Open("CherenkovPlots_"+run+".root");
    TVectorD* res=(TVectorD*)inputFile->Get("meanValue");
    TVectorD* resErr=(TVectorD*)inputFile->Get("meanErrValue");
    

    for(int j=0;j<5;++j){
      gr_reso_vs_energy1450[j]->SetPoint( i, energies[i], (*res)[j] );
      gr_reso_vs_energy1450[j]->SetPointError( i, 0, (*resErr)[j] );
      max1450[j]=(*res)[j];
    }

  }
 for(int i=0;i<5;++i){
    TString fibre;
    fibre.Form("%d",i); 

    TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
    c1->cd();
    TH2D* h2_axes2 = new TH2D( "axes", "", 100, -0.0, 305. , 1000, 0.0, 1.1*max1450[i]);
    h2_axes2->SetXTitle("Beam Energy [GeV]");
    h2_axes2->SetYTitle("CeF_{3} Response [ADC]");
    h2_axes2->Draw("");

    gr_reso_vs_energy1450[i]->SetMarkerStyle(20);
    gr_reso_vs_energy1450[i]->SetMarkerSize(1.6);
    gr_reso_vs_energy1450[i]->SetMarkerColor(kBlue);
    gr_reso_vs_energy1450[i]->Draw("p same");

        TF1 *fun= new TF1("fun","[0]*x+[1]",1, 250.+15.);
    //    TF1 *fun= new TF1("fun","sqrt([0]*[0]/x+[1]*[1])",1, 250.+15.);
    //    TF1 *fun= new TF1("fun",  "sqrt([0]*[0]/x+[1]*[1]+ [2]*[2]/(x*x))",1, 250+5.);
    fun->SetParameter(0, 1.);
    fun->SetParameter(1, 0.);
    gr_reso_vs_energy1450[i]->Fit(fun,"RN");
    fun->SetLineWidth(1.);
    fun->SetLineColor(kBlue+2);
    fun->Draw("L same");

    TLegend* leg_neat = new TLegend(0.42, 0.92-0.06*5 , 0.9, 0.92);
    leg_neat->SetTextSize(0.038);
    std::string ene="Data fibre ";
    //  ene+=Form("%.0f",energies[i]);
    //    ene+=" GeV";
    ene+=fibre;

    c1->SaveAs("plots_HVScan/mean_HV1450_"+fibre+".png");
  }

  //HV 1250
  float max1250[5]={0};
  for (int i=0;i<runs1250.size();++i){
    TString run;
    run.Form("%d",runs1250[i]); 
    TFile *inputFile=TFile::Open("CherenkovPlots_"+run+".root");
    TVectorD* res=(TVectorD*)inputFile->Get("meanValue");
    TVectorD* resErr=(TVectorD*)inputFile->Get("meanErrValue");
    

    for(int j=0;j<5;++j){
      gr_reso_vs_energy1250[j]->SetPoint( i, energies[i], (*res)[j] );
      gr_reso_vs_energy1250[j]->SetPointError( i, 0, (*resErr)[j] );
      max1250[j]=(*res)[j];
    }

  }
 for(int i=0;i<5;++i){
    TString fibre;
    fibre.Form("%d",i); 

    TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
    c1->cd();
    TH2D* h2_axes2 = new TH2D( "axes", "", 100, -0.0, 305. , 1000, 0.0, 1.1*max1250[i]);
    h2_axes2->SetXTitle("Beam Energy [GeV]");
    h2_axes2->SetYTitle("CeF_{3} Response [ADC]");
    h2_axes2->Draw("");

    gr_reso_vs_energy1250[i]->SetMarkerStyle(20);
    gr_reso_vs_energy1250[i]->SetMarkerSize(1.6);
    gr_reso_vs_energy1250[i]->SetMarkerColor(kBlue);
    gr_reso_vs_energy1250[i]->Draw("p same");

        TF1 *fun= new TF1("fun","[0]*x+[1]",1, 250.+15.);
    //    TF1 *fun= new TF1("fun","sqrt([0]*[0]/x+[1]*[1])",1, 250.+15.);
    //    TF1 *fun= new TF1("fun",  "sqrt([0]*[0]/x+[1]*[1]+ [2]*[2]/(x*x))",1, 250+5.);
    fun->SetParameter(0, 1.);
    fun->SetParameter(1, 0.);
    gr_reso_vs_energy1250[i]->Fit(fun,"RN");
    fun->SetLineWidth(1.);
    fun->SetLineColor(kBlue+2);
    fun->Draw("L same");

    TLegend* leg_neat = new TLegend(0.42, 0.92-0.06*5 , 0.9, 0.92);
    leg_neat->SetTextSize(0.038);
    std::string ene="Data fibre ";
    //  ene+=Form("%.0f",energies[i]);
    //    ene+=" GeV";
    ene+=fibre;

    c1->SaveAs("plots_HVScan/mean_HV1250_"+fibre+".png");
  }

  //HV 1050
  float max1050[5]={0};
  for (int i=0;i<runs1050.size();++i){
    TString run;
    run.Form("%d",runs1050[i]); 
    TFile *inputFile=TFile::Open("CherenkovPlots_"+run+".root");
    TVectorD* res=(TVectorD*)inputFile->Get("meanValue");
    TVectorD* resErr=(TVectorD*)inputFile->Get("meanErrValue");
    

    for(int j=0;j<5;++j){
      gr_reso_vs_energy1050[j]->SetPoint( i, energies[i], (*res)[j] );
      gr_reso_vs_energy1050[j]->SetPointError( i, 0, (*resErr)[j] );
      max1050[j]=(*res)[j];
    }

  }
 for(int i=0;i<5;++i){
    TString fibre;
    fibre.Form("%d",i); 

    TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
    c1->cd();
    TH2D* h2_axes2 = new TH2D( "axes", "", 100, -0.0, 305. , 1000, 0.0, 1.1*max1050[i]);
    h2_axes2->SetXTitle("Beam Energy [GeV]");
    h2_axes2->SetYTitle("CeF_{3} Response [ADC]");
    h2_axes2->Draw("");

    gr_reso_vs_energy1050[i]->SetMarkerStyle(20);
    gr_reso_vs_energy1050[i]->SetMarkerSize(1.6);
    gr_reso_vs_energy1050[i]->SetMarkerColor(kBlue);
    gr_reso_vs_energy1050[i]->Draw("p same");

        TF1 *fun= new TF1("fun","[0]*x+[1]",1, 250.+15.);
    //    TF1 *fun= new TF1("fun","sqrt([0]*[0]/x+[1]*[1])",1, 250.+15.);
    //    TF1 *fun= new TF1("fun",  "sqrt([0]*[0]/x+[1]*[1]+ [2]*[2]/(x*x))",1, 250+5.);
    fun->SetParameter(0, 1.);
    fun->SetParameter(1, 0.);
    gr_reso_vs_energy1050[i]->Fit(fun,"RN");
    fun->SetLineWidth(1.);
    fun->SetLineColor(kBlue+2);
    fun->Draw("L same");

    TLegend* leg_neat = new TLegend(0.42, 0.92-0.06*5 , 0.9, 0.92);
    leg_neat->SetTextSize(0.038);
    std::string ene="Data fibre ";
    //  ene+=Form("%.0f",energies[i]);
    //    ene+=" GeV";
    ene+=fibre;

    c1->SaveAs("plots_HVScan/mean_HV1050_"+fibre+".png");
  }
  //HV 850
  float max850[5]={0};
  for (int i=0;i<runs850.size();++i){
    TString run;
    run.Form("%d",runs850[i]); 
    TFile *inputFile=TFile::Open("CherenkovPlots_"+run+".root");
    TVectorD* res=(TVectorD*)inputFile->Get("meanValue");
    TVectorD* resErr=(TVectorD*)inputFile->Get("meanErrValue");
    

    for(int j=0;j<5;++j){
      gr_reso_vs_energy850[j]->SetPoint( i, energies[i], (*res)[j] );
      gr_reso_vs_energy850[j]->SetPointError( i, 0, (*resErr)[j] );
      max850[j]=(*res)[j];
    }

  }
 for(int i=0;i<5;++i){
    TString fibre;
    fibre.Form("%d",i); 

    TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
    c1->cd();
    TH2D* h2_axes2 = new TH2D( "axes", "", 100, -0.0, 305. , 1000, 0.0, 1.1*max850[i]);
    h2_axes2->SetXTitle("Beam Energy [GeV]");
    h2_axes2->SetYTitle("CeF_{3} Response [ADC]");
    h2_axes2->Draw("");

    gr_reso_vs_energy850[i]->SetMarkerStyle(20);
    gr_reso_vs_energy850[i]->SetMarkerSize(1.6);
    gr_reso_vs_energy850[i]->SetMarkerColor(kBlue);
    gr_reso_vs_energy850[i]->Draw("p same");

        TF1 *fun= new TF1("fun","[0]*x+[1]",1, 250.+15.);
    //    TF1 *fun= new TF1("fun","sqrt([0]*[0]/x+[1]*[1])",1, 250.+15.);
    //    TF1 *fun= new TF1("fun",  "sqrt([0]*[0]/x+[1]*[1]+ [2]*[2]/(x*x))",1, 250+5.);
    fun->SetParameter(0, 1.);
    fun->SetParameter(1, 0.);
    gr_reso_vs_energy850[i]->Fit(fun,"RN");
    fun->SetLineWidth(1.);
    fun->SetLineColor(kBlue+2);
    fun->Draw("L same");

    TLegend* leg_neat = new TLegend(0.42, 0.92-0.06*5 , 0.9, 0.92);
    leg_neat->SetTextSize(0.038);
    std::string ene="Data fibre ";
    //  ene+=Form("%.0f",energies[i]);
    //    ene+=" GeV";
    ene+=fibre;

    c1->SaveAs("plots_HVScan/mean_HV850_"+fibre+".png");
  }



 //total
 for(int i=0;i<5;++i){
    TString fibre;
    fibre.Form("%d",i); 

    TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
    c1->cd();
    TH2D* h2_axes2 = new TH2D( "axes", "", 100, -0.0, 305. , 1000, 0.0, 1.1*max1450[i]);
    h2_axes2->SetXTitle("Beam Energy [GeV]");
    h2_axes2->SetYTitle("CeF_{3} Response [ADC]");
    h2_axes2->Draw("");

    gr_reso_vs_energy1450[i]->SetMarkerColor(kBlack);
    gr_reso_vs_energy1450[i]->Draw("p same");
    gr_reso_vs_energy1250[i]->SetMarkerColor(kRed);
    gr_reso_vs_energy1250[i]->Draw("p same");
    gr_reso_vs_energy1050[i]->SetMarkerColor(kBlue);
    gr_reso_vs_energy1050[i]->Draw("p same");
    gr_reso_vs_energy850[i]->SetMarkerColor(kMagenta);
    gr_reso_vs_energy850[i]->Draw("p same");


    TLegend* leg_neat = new TLegend(0.8, 0.92-0.06*5 , 0.9, 0.92);
    leg_neat->SetTextSize(0.038);
    leg_neat->AddEntry(gr_reso_vs_energy1450[i],"HV 1450","p" );
    leg_neat->AddEntry(gr_reso_vs_energy1250[i],"HV 1250","p" );
    leg_neat->AddEntry(gr_reso_vs_energy1050[i],"HV 1050","p" );
    leg_neat->AddEntry(gr_reso_vs_energy850[i],"HV 850","p" );
    leg_neat->SetFillColor(0);
    //    leg_neat->Draw("same");

    c1->SaveAs("plots_HVScan/mean_Total_"+fibre+".png");
  }



  return 0;

}
