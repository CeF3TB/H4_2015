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
  
  std::string tag = "V00";
  if( argc>1 ) {
    std::string tag_str(argv[1]);
    tag = tag_str;
  }

  TGraphErrors* gr_reso_vs_energy[5];
  TGraphErrors* gr_resoMaxAmplFit_vs_energy[5];
  for(int i=0;i<5;++i) {
    gr_reso_vs_energy[i]= new TGraphErrors(0);
    gr_resoMaxAmplFit_vs_energy[i]= new TGraphErrors(0);
  }

  std::string outdir = "plots_reso_" + tag;
  system( Form("mkdir -p %s", outdir.c_str()) );

  float sigma_noise[5];//sigma of maxAmpl in pedestal events
//  sigma_noise[0]=2.47;
//  sigma_noise[1]=2.54;
//  sigma_noise[2]=2.19;
//  sigma_noise[3]=2.15;
  sigma_noise[0]=1.78;
  sigma_noise[1]=1.24;
  sigma_noise[2]=1.66;
  sigma_noise[3]=1.05;
  sigma_noise[4]=sqrt(2.49*2.49*sigma_noise[0]*sigma_noise[0]+0.81*0.81*sigma_noise[1]*sigma_noise[1]+16.8*16.8*sigma_noise[2]*sigma_noise[2]+0.73*0.73*sigma_noise[3]*sigma_noise[3]);//weighted mean of sigmas with intercalib of V03 //FIXME, don't do it hardcoded
  bool subtract_noise=true;

  for (int i=0;i<runs.size();++i){
    TString run;
    run.Form("%d",runs[i]); 
    
       

    TFile *inputFile=TFile::Open("CherenkovPlots_"+run+"_"+tag.c_str()+".root");
    TVectorD* res=(TVectorD*)inputFile->Get("resValue");
    TVectorD* resErr=(TVectorD*)inputFile->Get("resErrValue");

    TVectorD* meanValuemaxAmpl_fit=(TVectorD*)inputFile->Get("meanValuemaxAmpl_fit");

    TVectorD* resMaxAmplFit=(TVectorD*)inputFile->Get("resValuemaxAmpl_fit");
    TVectorD* resErrMaxAmplFit=(TVectorD*)inputFile->Get("resErrValuemaxAmpl_fit");

    for(int j=0;j<5;++j){
      gr_reso_vs_energy[j]->SetPoint( i, energies[i], (*res)[j] );
      gr_reso_vs_energy[j]->SetPointError( i, 0, (*resErr)[j] );

      if(!subtract_noise){
	gr_resoMaxAmplFit_vs_energy[j]->SetPoint( i, energies[i], (*resMaxAmplFit)[j] );
	gr_resoMaxAmplFit_vs_energy[j]->SetPointError( i, 0, (*resErrMaxAmplFit)[j] );
      }else{
	std::cout<<"run:"<<runs[i]<<" res_original:"<<(*resMaxAmplFit)[j]<<" sigma_noise:"<<sigma_noise[j]<<" final:"<<sqrt((*resMaxAmplFit)[j]*(*resMaxAmplFit)[j]-100*sigma_noise[j]/(*meanValuemaxAmpl_fit)[j]*100*sigma_noise[j]/(*meanValuemaxAmpl_fit)[j])<<" mean:"<<(*meanValuemaxAmpl_fit)[j]<<std::endl;
	gr_resoMaxAmplFit_vs_energy[j]->SetPoint( i, energies[i], sqrt((*resMaxAmplFit)[j]*(*resMaxAmplFit)[j]-100*sigma_noise[j]/(*meanValuemaxAmpl_fit)[j]*100*sigma_noise[j]/(*meanValuemaxAmpl_fit)[j]) );
	gr_resoMaxAmplFit_vs_energy[j]->SetPointError( i, 0, (*resErrMaxAmplFit)[j] );
      }

    }
  }

  for(int i=0;i<5;++i){
    TString fibre;
    fibre.Form("%d",i); 

    if(gr_reso_vs_energy[i]){
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

      //              TF1 *fun= new TF1("fun","sqrt([0]*[0]/(x*x)+[1]*[1])",1, 250.+15.);
      //            TF1 *fun= new TF1("fun","sqrt([0]*[0]/x+[1]*[1])",1, 250.+15.);
            TF1 *fun= new TF1("fun",  "sqrt([0]*[0]/x+[1]*[1]+ [2]*[2]/(x*x))",1, 250+5.);
      fun->SetParameter(1, 4.);
      fun->SetParameter(0, 20.);
      fun->SetParameter(2, 20.);
      gr_reso_vs_energy[i]->Fit(fun,"RN");
      fun->SetLineWidth(1.);
      fun->SetLineColor(kBlue+2);
      fun->Draw("L same");

      TLegend* leg_neat = new TLegend(0.4, 0.92-0.06*5 , 0.8, 0.92);
      leg_neat->SetTextSize(0.038);
      string ene="Data fibre ";
      //  ene+=Form("%.0f",energies[i]);
      //    ene+=" GeV";
      ene+=fibre;
      if(i!=4)    leg_neat->AddEntry(gr_reso_vs_energy[i],ene.c_str(),"p");
      else leg_neat->AddEntry(gr_reso_vs_energy[i],"Data CeF_{3}","p");
      //      leg_neat->AddEntry((TObject*)0 ,Form("S =  %.2f\n%s #pm %.2f / #sqrt{E [GeV]}",fun->GetParameter(0),"%",fun->GetParError(0) ),"");
      leg_neat->AddEntry((TObject*)0 ,Form("S =  %.2f\n%s #pm %.2f",fun->GetParameter(0),"%",fun->GetParError(0) ),"");
      leg_neat->AddEntry( (TObject*)0 ,Form("C =  %.2f\n%s  #pm %.2f",(fun->GetParameter(1)) ,"%",fun->GetParError(1) ),"");
      leg_neat->AddEntry( (TObject*)0 ,Form("N =  %.2f\n%s #pm %.2f",(fun->GetParameter(2))/100 ," GeV" ,fun->GetParError(2)/100),"");
      leg_neat->SetFillColor(0);
      leg_neat->Draw("same");
      std::cout<<"parameters: S="<<fun->GetParameter(0)<<" C="<<fun->GetParameter(1)<<" N="<<fun->GetParameter(2)<<std::endl;


      c1->SaveAs(outdir+"/reso_HV1450_"+fibre+".png");
      c1->SaveAs(outdir+"/reso_HV1450_"+fibre+".pdf");
    }

    if(gr_resoMaxAmplFit_vs_energy[i]){
      TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
      c1->cd();
      TH2D* h2_axes2 = new TH2D( "axes", "", 100, -0.0, 305. , 10, 0.0, 25 );
      h2_axes2->SetXTitle("Beam Energy [GeV]");
      h2_axes2->SetYTitle("Energy Resolution [%]");
      h2_axes2->Draw("");
      
      gr_resoMaxAmplFit_vs_energy[i]->SetMarkerStyle(20);
      gr_resoMaxAmplFit_vs_energy[i]->SetMarkerSize(1.6);
      gr_resoMaxAmplFit_vs_energy[i]->SetMarkerColor(kBlue);
      gr_resoMaxAmplFit_vs_energy[i]->Draw("p same");

      //              TF1 *fun= new TF1("fun","sqrt([0]*[0]/(x*x)+[1]*[1])",1, 250.+15.);
                  TF1 *fun= new TF1("fun","sqrt([0]*[0]/x+[1]*[1])",1, 250.+15.);
      //      TF1 *fun= new TF1("fun",  "sqrt([0]*[0]/x+[1]*[1]+ [2]*[2]/(x*x))",1, 250+5.);
      fun->SetParameter(1, 4.);
      fun->SetParameter(0, 20.);
      fun->SetParameter(2, 20.);
      gr_resoMaxAmplFit_vs_energy[i]->Fit(fun,"RN");
      fun->SetLineWidth(1.);
      fun->SetLineColor(kBlue+2);
      fun->Draw("L same");


      TLegend* leg_neat = new TLegend(0.42, 0.92-0.06*5 , 0.85, 0.92);
      leg_neat->SetTextSize(0.038);
      string ene="Data fibre ";
      //  ene+=Form("%.0f",energies[i]);
      //    ene+=" GeV";
      ene+=fibre;
      if(i!=4)    leg_neat->AddEntry(gr_resoMaxAmplFit_vs_energy[i],ene.c_str(),"p");
      else leg_neat->AddEntry(gr_resoMaxAmplFit_vs_energy[i],"Data CeF_{3}","p");
      //leg_neat->AddEntry((TObject*)0 ,Form("S =  %.2f\n%s / #sqrt{E [GeV]}",fun->GetParameter(0),"%" ),"");
      leg_neat->AddEntry((TObject*)0 ,Form("S =  %.2f\n%s #pm %.2f",fun->GetParameter(0),"%",fun->GetParError(0) ),"");
      leg_neat->AddEntry( (TObject*)0 ,Form("C =  %.2f\n%s  #pm %.2f",(fun->GetParameter(1)) ,"%",fun->GetParError(1) ),"");
      //      leg_neat->AddEntry( (TObject*)0 ,Form("N =  %.2f\n%s #pm %.2f",(fun->GetParameter(2))/100 ," GeV" ,fun->GetParError(2)/100),"");

      leg_neat->SetFillColor(0);
      leg_neat->Draw("same");
      std::cout<<"parameters: S="<<fun->GetParameter(0)<<" C="<<fun->GetParameter(1)<<" N="<<fun->GetParameter(2)<<std::endl;
      c1->SaveAs(outdir+"/reso_maxAmplFit_HV1450_"+fibre+".png");
      c1->SaveAs(outdir+"/reso_maxAmplFit_HV1450_"+fibre+".pdf");
    }

  }
  return 0;
}
