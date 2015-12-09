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
  std::vector<int> runsPMTHV1250;
  std::vector<int> runsPMTHV1050;
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

  runsPMTHV1250.push_back(4090);
  runsPMTHV1250.push_back(4077);
  runsPMTHV1250.push_back(4064);
  runsPMTHV1250.push_back(4103);

  runsPMTHV1050.push_back(4091);
  runsPMTHV1050.push_back(4078);
  runsPMTHV1050.push_back(4065);
  runsPMTHV1050.push_back(4104);//??there is also 50k in 4107


  runsMAPDDiff.push_back(4497);
  runsMAPDDiff.push_back(4496);
  runsMAPDDiff.push_back(4493);
  runsMAPDDiff.push_back(4492);


  std::vector<TGraphErrors*> gr_mean_vs_energyMAPDnoDiff;
  std::vector<TGraphErrors*> gr_mean_vs_energyPMT;
  std::vector<TGraphErrors*> gr_mean_vs_energyPMTHV1250;
  std::vector<TGraphErrors*> gr_mean_vs_energyPMTHV1050;
  std::vector<TGraphErrors*> gr_mean_vs_energyMAPDDiff;

  std::vector<TGraphErrors*> gr_mean_vs_energyMAPDnoDiff_fit;
  std::vector<TGraphErrors*> gr_mean_vs_energyPMT_fit;
  std::vector<TGraphErrors*> gr_mean_vs_energyPMTHV1250_fit;
  std::vector<TGraphErrors*> gr_mean_vs_energyPMTHV1050_fit;
  std::vector<TGraphErrors*> gr_mean_vs_energyMAPDDiff_fit;

  //pol1

  std::vector<TGraphErrors*> gr_mean_vs_energyMAPDnoDiff_pol1;
  std::vector<TGraphErrors*> gr_mean_vs_energyPMT_pol1;
  std::vector<TGraphErrors*> gr_mean_vs_energyPMTHV1250_pol1;
  std::vector<TGraphErrors*> gr_mean_vs_energyPMTHV1050_pol1;
  std::vector<TGraphErrors*> gr_mean_vs_energyMAPDDiff_pol1;

  std::vector<TGraphErrors*> gr_mean_vs_energyMAPDnoDiff_fit_pol1;
  std::vector<TGraphErrors*> gr_mean_vs_energyPMT_fit_pol1;
  std::vector<TGraphErrors*> gr_mean_vs_energyPMTHV1250_fit_pol1;
  std::vector<TGraphErrors*> gr_mean_vs_energyPMTHV1050_fit_pol1;
  std::vector<TGraphErrors*> gr_mean_vs_energyMAPDDiff_fit_pol1;



  for(int i=0;i<energies.size();++i){
      gr_mean_vs_energyMAPDnoDiff.push_back( new TGraphErrors(0));
      gr_mean_vs_energyPMT.push_back( new TGraphErrors(0));
      gr_mean_vs_energyPMTHV1250.push_back( new TGraphErrors(0));
      gr_mean_vs_energyPMTHV1050.push_back( new TGraphErrors(0));
      gr_mean_vs_energyMAPDDiff.push_back( new TGraphErrors(0));

      gr_mean_vs_energyMAPDnoDiff_pol1.push_back( new TGraphErrors(0));
      gr_mean_vs_energyPMT_pol1.push_back( new TGraphErrors(0));
      gr_mean_vs_energyPMTHV1250_pol1.push_back( new TGraphErrors(0));
      gr_mean_vs_energyPMTHV1050_pol1.push_back( new TGraphErrors(0));
      gr_mean_vs_energyMAPDDiff_pol1.push_back( new TGraphErrors(0));


      gr_mean_vs_energyMAPDnoDiff_fit.push_back( new TGraphErrors(0));
      gr_mean_vs_energyPMT_fit.push_back( new TGraphErrors(0));
      gr_mean_vs_energyPMTHV1250_fit.push_back( new TGraphErrors(0));
      gr_mean_vs_energyPMTHV1050_fit.push_back( new TGraphErrors(0));
      gr_mean_vs_energyMAPDDiff_fit.push_back( new TGraphErrors(0));

      gr_mean_vs_energyMAPDnoDiff_fit_pol1.push_back( new TGraphErrors(0));
      gr_mean_vs_energyPMT_fit_pol1.push_back( new TGraphErrors(0));
      gr_mean_vs_energyPMTHV1250_fit_pol1.push_back( new TGraphErrors(0));
      gr_mean_vs_energyPMTHV1050_fit_pol1.push_back( new TGraphErrors(0));
      gr_mean_vs_energyMAPDDiff_fit_pol1.push_back( new TGraphErrors(0));

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

      gr_mean_vs_energyMAPDnoDiff_pol1[j]->SetPoint( i, energies[i], (*mean)[j]);
      gr_mean_vs_energyMAPDnoDiff_pol1[j]->SetPointError( i, 0, (*meanErr)[j]);




      gr_mean_vs_energyMAPDnoDiff_fit[j]->SetPoint( i, energies[i], ((*mean_fit)[j]/meanFirst_fit[j])*(energyFirst[j]/energies[i]));
      gr_mean_vs_energyMAPDnoDiff_fit[j]->SetPointError( i, 0, (energyFirst[j]/energies[i])*getRatioError((*mean_fit)[j],meanFirst_fit[j],(*meanErr_fit)[j],meanFirstErr_fit[j]));

      gr_mean_vs_energyMAPDnoDiff_fit_pol1[j]->SetPoint( i, energies[i], (*mean_fit)[j]);
      gr_mean_vs_energyMAPDnoDiff_fit_pol1[j]->SetPointError( i, 0, (*meanErr_fit)[j]);


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


      if(j==0)std::cout<<"-----------PMT:"<<j<<" "<<(*mean)[j]<< " "<<meanFirst[j]<< " "<<energyFirst[j]<<" "<<energies[i]<<" "<<((*mean)[j]/meanFirst[j])*(energyFirst[j]/energies[i])<<std::endl;
      gr_mean_vs_energyPMT[j]->SetPoint( i, energies[i], ((*mean)[j]/meanFirst[j])*(energyFirst[j]/energies[i]));
      gr_mean_vs_energyPMT[j]->SetPointError( i, 0, (energyFirst[j]/energies[i])*getRatioError((*mean)[j],meanFirst[j],(*meanErr)[j],meanFirstErr[j]));


      gr_mean_vs_energyPMT_pol1[j]->SetPoint( i, energies[i], (*mean)[j]);
      gr_mean_vs_energyPMT_pol1[j]->SetPointError( i, 0, (*meanErr)[j]);


      gr_mean_vs_energyPMT_fit[j]->SetPoint( i, energies[i], ((*mean_fit)[j]/meanFirst_fit[j])*(energyFirst[j]/energies[i]));
      gr_mean_vs_energyPMT_fit[j]->SetPointError( i, 0, (energyFirst[j]/energies[i])*getRatioError((*mean_fit)[j],meanFirst_fit[j],(*meanErr_fit)[j],meanFirstErr_fit[j]));


      gr_mean_vs_energyPMT_fit_pol1[j]->SetPoint( i, energies[i], (*mean_fit)[j]);
      gr_mean_vs_energyPMT_fit_pol1[j]->SetPointError( i, 0, (*meanErr_fit)[j]);


    }

  }

  //PMTHV1250
  for (int i=0;i<runsPMTHV1250.size();++i){
    TString run;
    run.Form("%d",runsPMTHV1250[i]); 
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


      //      if(j==2)std::cout<<"-----------PMTHV1250:"<<j<<" "<<(*mean)[j]<< " "<<meanFirst[j]<< " "<<energyFirst[j]<<" "<<energies[i]<<" "<<((*mean)[j]/meanFirst[j])*(energyFirst[j]/energies[i])<<std::endl;
      gr_mean_vs_energyPMTHV1250[j]->SetPoint( i, energies[i], ((*mean)[j]/meanFirst[j])*(energyFirst[j]/energies[i]));
      gr_mean_vs_energyPMTHV1250[j]->SetPointError( i, 0, (energyFirst[j]/energies[i])*getRatioError((*mean)[j],meanFirst[j],(*meanErr)[j],meanFirstErr[j]));


      gr_mean_vs_energyPMTHV1250_pol1[j]->SetPoint( i, energies[i], (*mean)[j]);
      gr_mean_vs_energyPMTHV1250_pol1[j]->SetPointError( i, 0, (*meanErr)[j]);


      gr_mean_vs_energyPMTHV1250_fit[j]->SetPoint( i, energies[i], ((*mean_fit)[j]/meanFirst_fit[j])*(energyFirst[j]/energies[i]));
      gr_mean_vs_energyPMTHV1250_fit[j]->SetPointError( i, 0, (energyFirst[j]/energies[i])*getRatioError((*mean_fit)[j],meanFirst_fit[j],(*meanErr_fit)[j],meanFirstErr_fit[j]));

      gr_mean_vs_energyPMTHV1250_fit_pol1[j]->SetPoint( i, energies[i], (*mean_fit)[j]);
      gr_mean_vs_energyPMTHV1250_fit_pol1[j]->SetPointError( i, 0, (*meanErr_fit)[j]);


    }

  }


  //PMTHV1050
  for (int i=0;i<runsPMTHV1050.size();++i){
    TString run;
    run.Form("%d",runsPMTHV1050[i]); 
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


      //      if(j==2)std::cout<<"-----------PMTHV1050:"<<j<<" "<<(*mean)[j]<< " "<<meanFirst[j]<< " "<<energyFirst[j]<<" "<<energies[i]<<" "<<((*mean)[j]/meanFirst[j])*(energyFirst[j]/energies[i])<<std::endl;
      gr_mean_vs_energyPMTHV1050[j]->SetPoint( i, energies[i], ((*mean)[j]/meanFirst[j])*(energyFirst[j]/energies[i]));
      gr_mean_vs_energyPMTHV1050[j]->SetPointError( i, 0, (energyFirst[j]/energies[i])*getRatioError((*mean)[j],meanFirst[j],(*meanErr)[j],meanFirstErr[j]));

      gr_mean_vs_energyPMTHV1050_pol1[j]->SetPoint( i, energies[i], (*mean)[j]);
      gr_mean_vs_energyPMTHV1050_pol1[j]->SetPointError( i, 0, (*meanErr)[j]);



      gr_mean_vs_energyPMTHV1050_fit[j]->SetPoint( i, energies[i], ((*mean_fit)[j]/meanFirst_fit[j])*(energyFirst[j]/energies[i]));
      gr_mean_vs_energyPMTHV1050_fit[j]->SetPointError( i, 0, (energyFirst[j]/energies[i])*getRatioError((*mean_fit)[j],meanFirst_fit[j],(*meanErr_fit)[j],meanFirstErr_fit[j]));

      gr_mean_vs_energyPMTHV1050_fit_pol1[j]->SetPoint( i, energies[i], (*mean_fit)[j]);
      gr_mean_vs_energyPMTHV1050_fit_pol1[j]->SetPointError( i, 0, (*meanErr_fit)[j]);


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
      gr_mean_vs_energyMAPDDiff[j]->SetPointError( i, 0, (energyFirst[j]/energies[i])*getRatioError((*mean)[j],meanFirst[j],(*meanErr_fit)[j],meanFirstErr[j]));


      gr_mean_vs_energyMAPDDiff_pol1[j]->SetPoint( i, energies[i], (*mean)[j]);
      gr_mean_vs_energyMAPDDiff_pol1[j]->SetPointError( i, 0, (*meanErr)[j]);



      gr_mean_vs_energyMAPDDiff_fit[j]->SetPoint( i, energies[i], ((*mean_fit)[j]/meanFirst_fit[j])*(energyFirst[j]/energies[i]));
      gr_mean_vs_energyMAPDDiff_fit[j]->SetPointError( i, 0, (energyFirst[j]/energies[i])*getRatioError((*mean_fit)[j],meanFirst_fit[j],(*meanErr)[j],meanFirstErr_fit[j]));

      gr_mean_vs_energyMAPDDiff_fit_pol1[j]->SetPoint( i, energies[i], (*mean_fit)[j]);
      gr_mean_vs_energyMAPDDiff_fit_pol1[j]->SetPointError( i, 0, (*meanErr_fit)[j]);



    }

  }


  for(int i=0;i<energies.size();++i){

    TLegend* leg_neat = new TLegend(0.64, 0.75 ,0.85, 0.9);
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
    h2_axes2->SetYTitle("Fibre Normalized Response [ADC]");

    TH2D* h2_axes_zoom = new TH2D( "axes_zoom", "", 100, -0.0, 350. , 110, 0.8, 1.1);
    h2_axes_zoom->SetXTitle("Beam Energy [GeV]");
    h2_axes_zoom->SetYTitle("Fibre Normalized Response [ADC]");


    TH2D* h2_axes_zoom_pmt = new TH2D( "axes_zoom_pmt", "", 100, -0.0, 350. , 110, 0.9, 1.1);
    h2_axes_zoom_pmt->SetXTitle("Beam Energy [GeV]");
    h2_axes_zoom_pmt->SetYTitle("Fibre Normalized Response [ADC]");

    TH2D* h2_axes_pol1 = new TH2D( "axes_pol1", "", 100, -0.0, 250. , 110, 0., 1.1*(gr_mean_vs_energyPMT_pol1[i]->GetY())[energies.size()-1]);
    h2_axes_pol1->SetXTitle("Beam Energy [GeV]");
    h2_axes_pol1->SetYTitle("Fibre Response [ADC]");

    TH2D* h2_axes_pol1_2 = new TH2D( "axes_pol1_2", "", 100, -0.0, 250. , 110, 0., 1.1*(gr_mean_vs_energyMAPDDiff_pol1[i]->GetY())[energies.size()-1]);
    h2_axes_pol1_2->SetXTitle("Beam Energy [GeV]");
    h2_axes_pol1_2->SetYTitle("Fibre Response [ADC]");


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

    std::string dirName = "plots_linearity";
    system(Form("mkdir -p %s", dirName.c_str()));

    std::string outFile=dirName+"/mean_linearity_";
    c1->SaveAs(outFile.c_str()+fibre+".png");
    c1->SaveAs(outFile.c_str()+fibre+".pdf");
    //zoom
    c1->Clear();
    h2_axes_zoom->Draw("");
    gr_mean_vs_energyMAPDnoDiff[i]->Draw("p same");
    gr_mean_vs_energyPMT[i]->Draw("p same");
    gr_mean_vs_energyMAPDDiff[i]->Draw("p same");
    leg_neat->Draw("same");
    outFile=dirName+"/mean_linearity_zoom_";
    c1->SaveAs(outFile.c_str()+fibre+".png");
    c1->SaveAs(outFile.c_str()+fibre+".pdf");

    //pol1
    c1->Clear();

    h2_axes_pol1->Draw();

    gr_mean_vs_energyPMT_pol1[i]->SetMarkerStyle(20);
    gr_mean_vs_energyPMT_pol1[i]->SetMarkerSize(1.6);
    gr_mean_vs_energyPMT_pol1[i]->SetMarkerColor(kBlue);
    gr_mean_vs_energyPMT_pol1[i]->Draw("p same");


    TF1 *fun= new TF1("fun","pol1",1, 250.+15.);
    fun->FixParameter(0,0.);
    fun->FixParameter(0,0.);
    fun->SetLineWidth(1);
    fun->SetLineColor(kRed);


    TF1 *fun2= new TF1("fun","pol1",0, 250.+15.);
    fun2->FixParameter(0,0.);
    fun2->FixParameter(1,gr_mean_vs_energyPMT_pol1[i]->GetY()[0]/energies[0]);
    fun2->SetLineWidth(1);
    fun2->SetLineColor(kBlack);

    gr_mean_vs_energyPMT_pol1[i]->Fit(fun,"RN");
    fun->Draw("same");
    fun2->Draw("same");
    outFile=dirName+"/mean_linearity_PMT_pol1_";
    c1->SaveAs(outFile.c_str()+fibre+".png"); 
    c1->SaveAs(outFile.c_str()+fibre+".pdf"); 

    //pol1 fit
    c1->Clear();

    h2_axes_pol1->Draw();

    fun2->FixParameter(0,0.);
    fun2->FixParameter(1,gr_mean_vs_energyPMT_fit_pol1[i]->GetY()[0]/energies[0]);


    gr_mean_vs_energyPMT_fit_pol1[i]->SetMarkerStyle(20);
    gr_mean_vs_energyPMT_fit_pol1[i]->SetMarkerSize(1.6);
    gr_mean_vs_energyPMT_fit_pol1[i]->SetMarkerColor(kBlue);
    gr_mean_vs_energyPMT_fit_pol1[i]->Draw("p same");


    gr_mean_vs_energyPMT_fit_pol1[i]->Fit(fun,"RN");
    fun->Draw("same");
    fun2->Draw("same");
    outFile=dirName+"/mean_linearity_PMT_fit_pol1_";
    c1->SaveAs(outFile.c_str()+fibre+".png"); 
    c1->SaveAs(outFile.c_str()+fibre+".pdf"); 



    //pol1 MAPD
    c1->Clear();

    h2_axes_pol1_2->Draw();

    fun2->FixParameter(0,0.);
    fun2->FixParameter(1,gr_mean_vs_energyMAPDDiff_pol1[i]->GetY()[0]/energies[0]);


    gr_mean_vs_energyMAPDDiff_pol1[i]->SetMarkerStyle(20);
    gr_mean_vs_energyMAPDDiff_pol1[i]->SetMarkerSize(1.6);
    gr_mean_vs_energyMAPDDiff_pol1[i]->SetMarkerColor(kBlue);
    gr_mean_vs_energyMAPDDiff_pol1[i]->Draw("p same");


    fun->FixParameter(0,0.);
    gr_mean_vs_energyMAPDDiff_pol1[i]->Fit(fun,"RN");
    fun->Draw("same");
    fun2->Draw("same");
    outFile=dirName+"/mean_linearity_MAPDDiff_pol1_";
    c1->SaveAs(outFile.c_str()+fibre+".png"); 
    c1->SaveAs(outFile.c_str()+fibre+".pdf"); 

    //pol1 fit MAPD
    c1->Clear();

    h2_axes_pol1_2->Draw();

    fun2->FixParameter(0,0.);
    fun2->FixParameter(1,gr_mean_vs_energyMAPDDiff_fit_pol1[i]->GetY()[0]/energies[0]);


    gr_mean_vs_energyMAPDDiff_fit_pol1[i]->SetMarkerStyle(20);
    gr_mean_vs_energyMAPDDiff_fit_pol1[i]->SetMarkerSize(1.6);
    gr_mean_vs_energyMAPDDiff_fit_pol1[i]->SetMarkerColor(kBlue);
    gr_mean_vs_energyMAPDDiff_fit_pol1[i]->Draw("p same");


    gr_mean_vs_energyMAPDDiff_fit_pol1[i]->Fit(fun,"RN");
    fun->Draw("same");
    fun2->Draw("same");
    outFile=dirName+"/mean_linearity_MAPDDiff_fit_pol1_";
    c1->SaveAs(outFile.c_str()+fibre+".png"); 
    c1->SaveAs(outFile.c_str()+fibre+".pdf"); 




    //fit
    c1->Clear();

    h2_axes2->Draw("");
    gr_mean_vs_energyMAPDnoDiff_fit[i]->SetMarkerStyle(20);
    gr_mean_vs_energyMAPDnoDiff_fit[i]->SetMarkerSize(1.6);
    gr_mean_vs_energyMAPDnoDiff_fit[i]->SetMarkerColor(kBlack);
    gr_mean_vs_energyMAPDnoDiff_fit[i]->Draw("p same");

    gr_mean_vs_energyPMT_fit[i]->SetMarkerStyle(20);
    gr_mean_vs_energyPMT_fit[i]->SetMarkerSize(1.6);
    gr_mean_vs_energyPMT_fit[i]->SetMarkerColor(kBlue);
    gr_mean_vs_energyPMT_fit[i]->Draw("p same");

    gr_mean_vs_energyMAPDDiff_fit[i]->SetMarkerStyle(20);
    gr_mean_vs_energyMAPDDiff_fit[i]->SetMarkerSize(1.6);
    gr_mean_vs_energyMAPDDiff_fit[i]->SetMarkerColor(kRed);
    gr_mean_vs_energyMAPDDiff_fit[i]->Draw("p same");

    leg_neat->Draw("same");
    outFile=dirName+"/mean_linearity_fit_";
    c1->SaveAs(outFile.c_str()+fibre+".png");
    c1->SaveAs(outFile.c_str()+fibre+".pdf");

    //zoom
    c1->Clear();
    h2_axes_zoom->Draw("");
    gr_mean_vs_energyMAPDnoDiff_fit[i]->Draw("p same");
    gr_mean_vs_energyPMT_fit[i]->Draw("p same");
    gr_mean_vs_energyMAPDDiff_fit[i]->Draw("p same");
    leg_neat->Draw("same");
    outFile=dirName+"/mean_linearity_fit_zoom_";
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

    //pmt comparison
    c1->Clear();
    h2_axes_zoom_pmt->Draw();
    gr_mean_vs_energyPMT[i]->Draw("p same");

    gr_mean_vs_energyPMTHV1250[i]->SetMarkerStyle(20);
    gr_mean_vs_energyPMTHV1250[i]->SetMarkerSize(1.6);
    gr_mean_vs_energyPMTHV1250[i]->SetMarkerColor(kBlue+4);
    gr_mean_vs_energyPMTHV1250[i]->Draw("p same");

    gr_mean_vs_energyPMTHV1050[i]->SetMarkerStyle(20);
    gr_mean_vs_energyPMTHV1050[i]->SetMarkerSize(1.6);
    gr_mean_vs_energyPMTHV1050[i]->SetMarkerColor(kBlue-6);
    gr_mean_vs_energyPMTHV1050[i]->Draw("p same");

    TLegend* leg_pmt = new TLegend(0.64, 0.75 ,0.9, 0.9);
    leg_pmt->SetTextSize(0.036);
    leg_pmt->AddEntry(gr_mean_vs_energyPMT[i],"PMT high","p");
    leg_pmt->AddEntry(gr_mean_vs_energyPMTHV1250[i],"PMT medium","p");
    leg_pmt->AddEntry(gr_mean_vs_energyPMTHV1050[i],"PMT small","p");
    leg_pmt->SetFillColor(0);

    leg_pmt->Draw("same");
    outFile=dirName+"/mean_linearity_PMT_comparison_";
    c1->SaveAs(outFile.c_str()+fibre+".png");
    c1->SaveAs(outFile.c_str()+fibre+".pdf");

    //pmt comparison fit
    c1->Clear();
    h2_axes_zoom_pmt->Draw();
    gr_mean_vs_energyPMT_fit[i]->Draw("p same");

    gr_mean_vs_energyPMTHV1250_fit[i]->SetMarkerStyle(20);
    gr_mean_vs_energyPMTHV1250_fit[i]->SetMarkerSize(1.6);
    gr_mean_vs_energyPMTHV1250_fit[i]->SetMarkerColor(kBlue+4);
    gr_mean_vs_energyPMTHV1250_fit[i]->Draw("p same");

    gr_mean_vs_energyPMTHV1050_fit[i]->SetMarkerStyle(20);
    gr_mean_vs_energyPMTHV1050_fit[i]->SetMarkerSize(1.6);
    gr_mean_vs_energyPMTHV1050_fit[i]->SetMarkerColor(kBlue-6);
    gr_mean_vs_energyPMTHV1050_fit[i]->Draw("p same");

    leg_pmt->Draw("same");
    outFile=dirName+"/mean_linearity_fit_PMT_comparison_";
    c1->SaveAs(outFile.c_str()+fibre+".png");
    c1->SaveAs(outFile.c_str()+fibre+".pdf");


  }

  return 0;

}
