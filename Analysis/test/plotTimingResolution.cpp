#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <sstream> 

#include "DrawTools.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TGraphErrors.h"
#include "TVectorD.h"
#include "TF1.h"

int main( int argc, char* argv[] ) {

  DrawTools::setStyle();

   if( argc<2 ) {
     std::cout << "ERROR. You need to specify the name of the run you want to process." << std::endl;  
     exit(1);
   }

   std::string runName = "";
   std::string tag = "V00";

   if( argc>1 ) {
     std::string tag_str(argv[1]);
     tag = tag_str;
   } else {
     std::cout << "Usage:" << std::endl;
     std::cout << "./timingStudies [tag]" << std::endl;
     exit(12345);
   }


  std::vector<float> energies;
  energies.push_back(20);
  energies.push_back(50);
  energies.push_back(100);
  energies.push_back(150);
  energies.push_back(200);

  //5% stat error
  std::vector<int> wp;
//  wp.push_back(3);
//  wp.push_back(5);
//  wp.push_back(20);
//  wp.push_back(23);
//  wp.push_back(24);    

  wp.push_back(13);
  wp.push_back(13);
  wp.push_back(13);
  wp.push_back(13);
  wp.push_back(13);    



  std::vector<int> wp_100; 
  wp_100.push_back(5);	 
  wp_100.push_back(5);	 
  wp_100.push_back(5);	 
  wp_100.push_back(5);	 
  wp_100.push_back(5);    
  
  std::vector<int> runsSiPM;
  runsSiPM.push_back(4055);
  runsSiPM.push_back(4051);
  runsSiPM.push_back(4052);
  runsSiPM.push_back(4053);
  runsSiPM.push_back(4054);
  
  std::string constDirName = "plots_timing_SiPM";
  TString dir(constDirName);

  TGraphErrors* resVsEnergy = new TGraphErrors(0);
  TGraphErrors* resVsEnergy_opt = new TGraphErrors(0);

  std::vector<TGraphErrors*> resVsAmplitude;


  for (int i=0;i<runsSiPM.size();++i){
    TString run;
    run.Form("%d",runsSiPM[i]); 
    resVsAmplitude.push_back(new TGraphErrors(0));


    TFile* file = TFile::Open(dir+"/timingStudiesSiPM_"+tag+"_"+run+".root");

    if( file==0 ) {
      std::cout << "ERROR! Din't find file for run" << runsSiPM[i] << std::endl;
      std::cout << "Exiting." << std::endl;
      exit(11);
    }
    
    TVectorD* resValueTime=(TVectorD*)file->Get("resValueTime");
    TVectorD* resErrValueTime=(TVectorD*)file->Get("resErrValueTime");
    TVectorD* resValueAmplitude=(TVectorD*)file->Get("resValueAmplitude");
    TVectorD* resErrValueAmplitude=(TVectorD*)file->Get("resErrValueAmplitude");
    TVectorD* resErrRelativeValueTime=(TVectorD*)file->Get("resErrRelativeValueTime");
    TVectorD* nEntriesTime=(TVectorD*)file->Get("nEntriesTime");
    
    resVsEnergy->SetPoint( i, energies[i],(*resValueTime)[wp_100[i]] );
    resVsEnergy->SetPointError( i, 0, (*resErrValueTime)[wp_100[i]]);

    if(i>1) {
      resVsEnergy_opt->SetPoint( i-2, energies[i],(*resValueTime)[wp[i]] );
      resVsEnergy_opt->SetPointError( i-2, 0, (*resErrValueTime)[wp[i]]);
    }


    for (int j=0;j<resValueTime->GetNoElements();j++){
      if(j<2)continue;
      resVsAmplitude[i]->SetPoint(j-2,(j>0)*(20+(j>1)*j*7.5),(*resValueAmplitude)[j]);
      resVsAmplitude[i]->SetPointError(j-2,0,(*resErrValueAmplitude)[j]);
    }

    
  }
  

  TFile * outFile = new TFile(dir+"/plotsTimingStudiesSiPM_"+tag+".root","recreate");

  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  TH2D* h2_axes_1 = new TH2D( "axes_1", "", 100, -0.0, 250. , 110, 70., 1.1*((resVsEnergy->GetY())[0]+(resVsEnergy->GetErrorY(0))));

  TH2D* h2_axes_2 = new TH2D( "axes_2", "", 100, (resVsAmplitude[2]->GetX())[0]-10,(resVsAmplitude[2]->GetX())[resVsAmplitude[2]->GetN()-1] +10 , 110, 0., 1.1*(resVsAmplitude[2]->GetY())[0]);
  h2_axes_1->SetXTitle("Beam Energy [GeV]");
  h2_axes_1->SetYTitle("time_{Fibre}-time_{mcp} [ns]");

  TH2D* h2_axes_3 = new TH2D( "axes_1", "", 100, 80, 250. , 110, 80., 1.1*((resVsEnergy_opt->GetY())[0]+(resVsEnergy_opt->GetErrorY(0))));  
  
  h2_axes_1->Draw("");
  
  TF1* f_ene= new TF1("fun_ene","[1]/x+[0]",(resVsEnergy->GetX())[0]-10,(resVsEnergy->GetX())[resVsEnergy->GetN()-1] +10);
  resVsEnergy->Fit("fun_ene");

  TLegend* lego_1 = new TLegend(0.47, 0.7, 0.8, 0.92);
  lego_1->SetTextSize(0.038);
  lego_1->AddEntry(  (TObject*)0 ,"f(x) = p0 + p1/x", "");
  lego_1->AddEntry(  (TObject*)0 ,Form("p0 = %.0f #pm %.0f", f_ene->GetParameter(0), f_ene->GetParError(0) ), "");
  lego_1->AddEntry(  (TObject*)0 ,Form("p1 = %.0f #pm %.0f", f_ene->GetParameter(1), f_ene->GetParError(1) ), "");
  lego_1->SetFillColor(0);
  lego_1->Draw("same");
  


  resVsEnergy->SetMarkerStyle(20);
  resVsEnergy->SetMarkerSize(1.6);
  resVsEnergy->SetMarkerColor(kBlue);
  resVsEnergy->Draw("p same");
  resVsEnergy->SetName("resVsEnergy");
  resVsEnergy->Write();

  c1->SaveAs(dir+"/timingResolutionVsEnergySiPM.png");
  c1->SaveAs(dir+"/timingResolutionVsEnergySiPM.pdf");  


  h2_axes_3->Draw("");
  resVsEnergy_opt->Fit("fun_ene");

  TLegend* lego_2 = new TLegend(0.47, 0.7, 0.8, 0.92);
  lego_2->SetTextSize(0.038);
  lego_2->AddEntry(  (TObject*)0 ,"f(x) = p0 + p1/x", "");
  lego_2->AddEntry(  (TObject*)0 ,Form("p0 = %.0f #pm %.0f", f_ene->GetParameter(0), f_ene->GetParError(0) ), "");
  lego_2->AddEntry(  (TObject*)0 ,Form("p1 = %.0f #pm %.0f", f_ene->GetParameter(1), f_ene->GetParError(1) ), "");
  lego_2->SetFillColor(0);
  lego_2->Draw("same");


  resVsEnergy_opt->SetMarkerStyle(20);
  resVsEnergy_opt->SetMarkerSize(1.6);
  resVsEnergy_opt->SetMarkerColor(kBlue);
  resVsEnergy_opt->Draw("p same");
  resVsEnergy_opt->SetName("resVsEnergy_opt");
  resVsEnergy_opt->Write();

  c1->SaveAs(dir+"/timingResolutionVsEnergySiPM_opt.png");
  c1->SaveAs(dir+"/timingResolutionVsEnergySiPM_opt.pdf");  


  for(int i=0;i<runsSiPM.size();++i){ 
    TString ene;
    ene.Form("%.0f",energies[i]); 
    std::string energy(Form("%.0f", energies[i]));

    h2_axes_2->SetYTitle("#sigma_{t} [ps]");
    h2_axes_2->GetXaxis()->SetTitle("Amplitude [ADC]");
    h2_axes_2->Draw(""); 

    TF1* f= new TF1("fun","[1]/x+[0]",(resVsAmplitude[i]->GetX())[0]-10,(resVsAmplitude[i]->GetX())[resVsAmplitude[i]->GetN()-1] +10);
    resVsAmplitude[i]->Fit("fun");

    resVsAmplitude[i]->SetName("resVsAmplitude_"+ene);
    resVsAmplitude[i]->SetMarkerStyle(20);
    resVsAmplitude[i]->SetMarkerSize(1.6);
    resVsAmplitude[i]->SetMarkerColor(kBlue);
    resVsAmplitude[i]->Draw("p same");


    TPaveText* pave = DrawTools::getLabelTop_expOnXaxis(energy+" GeV Electron Beam");
    pave->Draw("same");



    TLegend* lego = new TLegend(0.47, 0.7, 0.8, 0.92);
    lego->SetTextSize(0.038);
    lego->AddEntry(  (TObject*)0 ,"f(x) = p0 + p1/x", "");
    lego->AddEntry(  (TObject*)0 ,Form("p0 = %.0f #pm %.0f", f->GetParameter(0), f->GetParError(0) ), "");
    lego->AddEntry(  (TObject*)0 ,Form("p1 = %.0f #pm %.0f", f->GetParameter(1), f->GetParError(1) ), "");
    lego->SetFillColor(0);
    lego->Draw("same");

    c1->SaveAs(dir+"/timingResolutionVsAmplitudeSiPM_"+ene+"GeV.png");
    c1->SaveAs(dir+"/timingResolutionVsAmplitudeSiPM_"+ene+"GeV.pdf");  

 

    resVsAmplitude[i]->Write();
  }

  outFile->Write();
  outFile->Close();

  return 0;

}
