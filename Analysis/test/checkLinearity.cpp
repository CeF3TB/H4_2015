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

float getRatioError( float num, float denom, float numErr, float denomErr ) {

  return sqrt( numErr*numErr/(denom*denom) + denomErr*denomErr*num*num/(denom*denom*denom*denom) );

}


int main( int argc, char* argv[] ) {



  std::string tag = "V10";

  if(argc>1){
    std::string tag_str(argv[1]);
    tag=tag_str;
  }

  DrawTools::setStyle();

  std::vector<float> energies;
  std::vector<int> runsMAPDnoDiff;
  std::vector<int> runsPMT;
  std::vector<int> runsMAPDDiff;

  //  energies.push_back(20);
  energies.push_back(50);
  energies.push_back(100);
  energies.push_back(150);
  energies.push_back(200);

  //  runsMAPDnoDiff.push_back(0);//FIXME! do we have a run at 20 GeV?
  runsMAPDnoDiff.push_back(3971);
  runsMAPDnoDiff.push_back(3990);
  runsMAPDnoDiff.push_back(4001);
  runsMAPDnoDiff.push_back(4026);

  runsPMT.push_back(4089);
  runsPMT.push_back(4076);
  runsPMT.push_back(4063);
  runsPMT.push_back(4102);

  runsMAPDDiff.push_back(4497);
  runsMAPDDiff.push_back(4496);
  runsMAPDDiff.push_back(4493);
  runsMAPDDiff.push_back(4492);


  std::vector<TGraphErrors*> gr_mean_vs_energyMAPDnoDiff;
  std::vector<TGraphErrors*> gr_mean_vs_energyPMT;
  std::vector<TGraphErrors*> gr_mean_vs_energyMAPDDiff;

  std::vector<TGraphErrors*> gr_mean_vs_energyMAPDnoDiff_fit;
  std::vector<TGraphErrors*> gr_mean_vs_energyPMT_fit;
  std::vector<TGraphErrors*> gr_mean_vs_energyMAPDDiff_fit;


  for(int i=0;i<energies.size();++i){
      gr_mean_vs_energyMAPDnoDiff.push_back( new TGraphErrors(0));
      gr_mean_vs_energyPMT.push_back( new TGraphErrors(0));
      gr_mean_vs_energyMAPDDiff.push_back( new TGraphErrors(0));

      gr_mean_vs_energyMAPDnoDiff_fit.push_back( new TGraphErrors(0));
      gr_mean_vs_energyPMT_fit.push_back( new TGraphErrors(0));
      gr_mean_vs_energyMAPDDiff_fit.push_back( new TGraphErrors(0));

   }

  // MAPDnoDiff
  float meanFirst[4];
  float meanFirstErr[4];
  float energyFirst[4];

  float meanFirst_fit[4];
  float meanFirstErr_fit[4];


  for (int i=0;i<runsMAPDnoDiff.size();++i){
    TString run;
    run.Form("%d",runsMAPDnoDiff[i]); 
    TFile *inputFile=TFile::Open("CherenkovPlots_"+run+"_"+tag.c_str()+".root");
    TVectorD* mean=(TVectorD*)inputFile->Get("meanValuemaxAmpl");
    TVectorD* meanErr=(TVectorD*)inputFile->Get("meanErrValuemaxAmpl");

    TVectorD* mean_fit=(TVectorD*)inputFile->Get("meanValuemaxAmpl_fit");
    TVectorD* meanErr_fit=(TVectorD*)inputFile->Get("meanErrValuemaxAmpl_fit");


    for(int j=0;j<4;++j){
      if(i==0){
	meanFirst[j]=(*mean)[j];
	meanFirstErr[j]=(*meanErr)[j];
	energyFirst[j]=energies[i];

	meanFirst_fit[j]=(*mean_fit)[j];
	meanFirstErr_fit[j]=(*meanErr_fit)[j];
      }


      //      std::cout<<j<<" "<<(*mean)[j]<< " "<<meanFirst[j]<< " "<<energyFirst[j]<<" "<<energies[i]<<" "<<((*mean)[j]/meanFirst[j])*(energyFirst[j]/energies[i])<<std::endl;
      gr_mean_vs_energyMAPDnoDiff[j]->SetPoint( i, energies[i], ((*mean)[j]/meanFirst[j])*(energyFirst[j]/energies[i]));
      gr_mean_vs_energyMAPDnoDiff[j]->SetPointError( i, 0, (energyFirst[j]/energies[i])*getRatioError((*mean)[j],meanFirst[j],(*meanErr)[j],meanFirstErr[j]));
      gr_mean_vs_energyMAPDnoDiff_fit[j]->SetPoint( i, energies[i], ((*mean_fit)[j]/meanFirst_fit[j])*(energyFirst[j]/energies[i]));
      gr_mean_vs_energyMAPDnoDiff_fit[j]->SetPointError( i, 0, (energyFirst[j]/energies[i])*getRatioError((*mean_fit)[j],meanFirst_fit[j],(*meanErr)[j],meanFirstErr_fit[j]));

    }

  }


  //PMT
  for (int i=0;i<runsPMT.size();++i){
    TString run;
    run.Form("%d",runsPMT[i]); 
    //    std::cout<<"CherenkovPlots_"+run+"_"+tag.c_str()+".root"<<std::endl;
    TFile *inputFile=TFile::Open("CherenkovPlots_"+run+"_"+tag.c_str()+".root");
    TVectorD* mean=(TVectorD*)inputFile->Get("meanValuemaxAmpl");
    TVectorD* meanErr=(TVectorD*)inputFile->Get("meanErrValuemaxAmpl");

    TVectorD* mean_fit=(TVectorD*)inputFile->Get("meanValuemaxAmpl_fit");
    TVectorD* meanErr_fit=(TVectorD*)inputFile->Get("meanErrValuemaxAmpl_fit");


    for(int j=0;j<4;++j){
      if(i==0){
	meanFirst[j]=(*mean)[j];
	meanFirstErr[j]=(*meanErr)[j];
	energyFirst[j]=energies[i];

	meanFirst_fit[j]=(*mean_fit)[j];
	meanFirstErr_fit[j]=(*meanErr_fit)[j];
      }


      //      if(j==2)std::cout<<"-----------PMT:"<<j<<" "<<(*mean)[j]<< " "<<meanFirst[j]<< " "<<energyFirst[j]<<" "<<energies[i]<<" "<<((*mean)[j]/meanFirst[j])*(energyFirst[j]/energies[i])<<std::endl;
      gr_mean_vs_energyPMT[j]->SetPoint( i, energies[i], ((*mean)[j]/meanFirst[j])*(energyFirst[j]/energies[i]));
      gr_mean_vs_energyPMT[j]->SetPointError( i, 0, (energyFirst[j]/energies[i])*getRatioError((*mean)[j],meanFirst[j],(*meanErr)[j],meanFirstErr[j]));

      gr_mean_vs_energyPMT_fit[j]->SetPoint( i, energies[i], ((*mean_fit)[j]/meanFirst_fit[j])*(energyFirst[j]/energies[i]));
      gr_mean_vs_energyPMT_fit[j]->SetPointError( i, 0, (energyFirst[j]/energies[i])*getRatioError((*mean_fit)[j],meanFirst_fit[j],(*meanErr)[j],meanFirstErr_fit[j]));

    }

  }

  //MAPDDiff
  for (int i=0;i<runsMAPDDiff.size();++i){
    TString run;
    run.Form("%d",runsMAPDDiff[i]); 
    //    std::cout<<"CherenkovPlots_"+run+"_"+tag.c_str()+".root"<<std::endl;
    TFile *inputFile=TFile::Open("CherenkovPlots_"+run+"_"+tag.c_str()+".root");
    TVectorD* mean=(TVectorD*)inputFile->Get("meanValuemaxAmpl");
    TVectorD* meanErr=(TVectorD*)inputFile->Get("meanErrValuemaxAmpl");

    TVectorD* mean_fit=(TVectorD*)inputFile->Get("meanValuemaxAmpl_fit");
    TVectorD* meanErr_fit=(TVectorD*)inputFile->Get("meanErrValuemaxAmpl_fit");

    for(int j=0;j<4;++j){
      if(i==0){
	meanFirst[j]=(*mean)[j];
	meanFirstErr[j]=(*meanErr)[j];
	energyFirst[j]=energies[i];

	meanFirst_fit[j]=(*mean_fit)[j];
	meanFirstErr_fit[j]=(*meanErr_fit)[j];
     }


      //            std::cout<<j<<" "<<(*mean)[j]<< " "<<meanFirst[j]<< " "<<energyFirst[j]<<" "<<energies[i]<<" "<<((*mean)[j]/meanFirst[j])*(energyFirst[j]/energies[i])<<std::endl;
      gr_mean_vs_energyMAPDDiff[j]->SetPoint( i, energies[i], ((*mean)[j]/meanFirst[j])*(energyFirst[j]/energies[i]));
      gr_mean_vs_energyMAPDDiff[j]->SetPointError( i, 0, (energyFirst[j]/energies[i])*getRatioError((*mean)[j],meanFirst[j],(*meanErr)[j],meanFirstErr[j]));

      gr_mean_vs_energyMAPDDiff_fit[j]->SetPoint( i, energies[i], ((*mean_fit)[j]/meanFirst_fit[j])*(energyFirst[j]/energies[i]));
      gr_mean_vs_energyMAPDDiff_fit[j]->SetPointError( i, 0, (energyFirst[j]/energies[i])*getRatioError((*mean_fit)[j],meanFirst_fit[j],(*meanErr)[j],meanFirstErr_fit[j]));


    }

  }


  for(int i=0;i<energies.size();++i){

    TLegend* leg_neat = new TLegend(0.64, 0.75 ,0.9, 0.9);
        leg_neat->SetTextSize(0.036);
    leg_neat->AddEntry(gr_mean_vs_energyMAPDnoDiff[i],"MAPD","p");
    leg_neat->AddEntry(gr_mean_vs_energyPMT[i],"PMT","p");
    leg_neat->AddEntry(gr_mean_vs_energyMAPDDiff[i],"MAPD+Diff","p");
    leg_neat->SetFillColor(0);

    TString fibre;
    fibre.Form("%d",i); 

    TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
    c1->cd();
    TH2D* h2_axes2 = new TH2D( "axes", "", 100, -0.0, 350. , 110, 0., 1.1);
    h2_axes2->SetXTitle("Beam Energy [GeV]");
    h2_axes2->SetYTitle("CeF_{3} Normalized Response [ADC]");
    h2_axes2->Draw("");

    gr_mean_vs_energyMAPDnoDiff[i]->SetMarkerStyle(20);
    gr_mean_vs_energyMAPDnoDiff[i]->SetMarkerSize(1.6);
    gr_mean_vs_energyMAPDnoDiff[i]->SetMarkerColor(kBlack);
    gr_mean_vs_energyMAPDnoDiff[i]->Draw("p same");

    gr_mean_vs_energyPMT[i]->SetMarkerStyle(20);
    gr_mean_vs_energyPMT[i]->SetMarkerSize(1.6);
    gr_mean_vs_energyPMT[i]->SetMarkerColor(kBlue);
    gr_mean_vs_energyPMT[i]->Draw("p same");

    gr_mean_vs_energyMAPDDiff[i]->SetMarkerStyle(20);
    gr_mean_vs_energyMAPDDiff[i]->SetMarkerSize(1.6);
    gr_mean_vs_energyMAPDDiff[i]->SetMarkerColor(kRed);
    gr_mean_vs_energyMAPDDiff[i]->Draw("p same");

    //    std::string ene="Data fibre ";
    //  ene+=Form("%.0f",energies[i]);
    //    ene+=" GeV";
    //    ene+=fibre;

    leg_neat->Draw("same");

    std::string dirName = "plots_linearity";
    system(Form("mkdir -p %s", dirName.c_str()));

    std::string outFile=dirName+"/mean_linearity_";
    c1->SaveAs(outFile.c_str()+fibre+".png");
    c1->SaveAs(outFile.c_str()+fibre+".pdf");

    //fit
    c1->Clear();

    h2_axes2->Draw("");
    gr_mean_vs_energyMAPDnoDiff[i]->SetMarkerStyle(20);
    gr_mean_vs_energyMAPDnoDiff[i]->SetMarkerSize(1.6);
    gr_mean_vs_energyMAPDnoDiff[i]->SetMarkerColor(kBlack);
    gr_mean_vs_energyMAPDnoDiff[i]->Draw("p same");

    gr_mean_vs_energyPMT[i]->SetMarkerStyle(20);
    gr_mean_vs_energyPMT[i]->SetMarkerSize(1.6);
    gr_mean_vs_energyPMT[i]->SetMarkerColor(kBlue);
    gr_mean_vs_energyPMT[i]->Draw("p same");

    gr_mean_vs_energyMAPDDiff[i]->SetMarkerStyle(20);
    gr_mean_vs_energyMAPDDiff[i]->SetMarkerSize(1.6);
    gr_mean_vs_energyMAPDDiff[i]->SetMarkerColor(kRed);
    gr_mean_vs_energyMAPDDiff[i]->Draw("p same");

    leg_neat->Draw("same");
    outFile=dirName+"/mean_linearity_fit_";
    c1->SaveAs(outFile.c_str()+fibre+".png");
    c1->SaveAs(outFile.c_str()+fibre+".pdf");


    c1->Clear();
    h2_axes2->Draw("");
    gr_mean_vs_energyMAPDnoDiff[i]->Draw("p same");
    outFile=dirName+"/mean_linearity_MAPDnoDiff_";
    c1->SaveAs(outFile.c_str()+fibre+".png");
    c1->SaveAs(outFile.c_str()+fibre+".pdf");

    c1->Clear();
    h2_axes2->Draw("");
    gr_mean_vs_energyPMT[i]->Draw("p same");
    outFile=dirName+"/mean_linearity_PMT_";
    c1->SaveAs(outFile.c_str()+fibre+".png");
    c1->SaveAs(outFile.c_str()+fibre+".pdf");

    c1->Clear();
    h2_axes2->Draw("");
    gr_mean_vs_energyMAPDDiff[i]->Draw("p same");
    outFile=dirName+"/mean_linearity_MAPDDiff_";
    c1->SaveAs(outFile.c_str()+fibre+".png");
    c1->SaveAs(outFile.c_str()+fibre+".pdf");

    //fit
    c1->Clear();
    h2_axes2->Draw("");
    gr_mean_vs_energyMAPDnoDiff_fit[i]->Draw("p same");
    outFile=dirName+"/mean_linearity_fit_MAPDnoDiff_";
    c1->SaveAs(outFile.c_str()+fibre+".png");
    c1->SaveAs(outFile.c_str()+fibre+".pdf");

    c1->Clear();
    h2_axes2->Draw("");
    gr_mean_vs_energyPMT_fit[i]->Draw("p same");
    outFile=dirName+"/mean_linearity_fit_PMT_";
    c1->SaveAs(outFile.c_str()+fibre+".png");
    c1->SaveAs(outFile.c_str()+fibre+".pdf");

    c1->Clear();
    h2_axes2->Draw("");
    gr_mean_vs_energyMAPDDiff_fit[i]->Draw("p same");
    outFile=dirName+"/mean_linearity_fit_MAPDDiff_";
    c1->SaveAs(outFile.c_str()+fibre+".png");
    c1->SaveAs(outFile.c_str()+fibre+".pdf");



  }

  return 0;

}
