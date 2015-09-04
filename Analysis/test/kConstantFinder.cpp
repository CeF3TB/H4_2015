
#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <sstream> 

#include "DrawTools.h"
#include "RecoTree.h"
#include "channelInfo.h"

#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooFFTConvPdf.h"
#include "RooPlot.h"
#include "RooCBShape.h"
#include "RooCruijff.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TVectorD.h"
#include "TGraphErrors.h"

float getRatioError( float num, float denom, float numErr, float denomErr ) {

  return sqrt( numErr*numErr/(denom*denom) + denomErr*denomErr*num*num/(denom*denom*denom*denom) );

}


int main( int argc, char* argv[] ) {

  DrawTools::setStyle();



   if( argc<2 ) {
     std::cout << "ERROR. You need to specify the name of the run you want to process." << std::endl;  
     exit(1);
   }

   std::string runName = "";
   std::string tag = "V00";

   if( argc>1 ) {
     std::string runName_str(argv[1]);
     runName = runName_str;
     if( argc>2 ) {
       std::string tag_str(argv[2]);
       tag = tag_str;
     }
   }else {

     std::cout << "Usage:" << std::endl;
     std::cout << "./kConstantFinder [runName] ([tag])" << std::endl;
     exit(12345);

   }

  TString runNumberString(runName);

  std::string fileName;
  fileName = "analysisTrees_"+tag+"/Reco_"+runName+".root";
  TFile* file = TFile::Open(fileName.c_str());

   TTree* recoTree=(TTree*)file->Get("recoTree");
   RecoTree t(recoTree);
   
   if (t.fChain == 0) return 0;

   TFile* outFile;
   outFile = TFile::Open("kConstantFinder_"+runNumberString+"_"+tag+".root","recreate");


   int NSTEPS=500;
   TVectorD resValue(NSTEPS);
   TVectorD resErrValue(NSTEPS);
   TH1F*   totalHistos_tight_maxAmpl_fit_tot[NSTEPS];  
   TVectorD cValue(NSTEPS);
   TVectorD cValueErr(NSTEPS);

   for(int i= 0; i< NSTEPS;i++){
     std::cout<<"c="<<(i*0.1)<<std::endl;
     cValue[i]=i*0.1;
     cValueErr[i]=0;

     TString step;
     step.Form("%d",i); 
    
     totalHistos_tight_maxAmpl_fit_tot[i]= new TH1F ("histo_step_"+step,"",200*4,0,4000*4);   

     Long64_t nentries = t.fChain->GetEntries();
     Long64_t nbytes = 0, nb = 0;
     for (Long64_t jentry=0; jentry<nentries;jentry++) {
       Long64_t ientry = t.LoadTree(jentry);
       if (ientry < 0) break;
      nb = t.fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      //      if(jentry%1000==0)std::cout<<"Entry:"<<jentry<<std::endl;

      if((t.nClusters_hodoX1==1||t.pos_2FibClust_hodoX1>-999) && (t.nClusters_hodoX2==1||t.pos_2FibClust_hodoX2>-999) && (t.nClusters_hodoY1==1||t.pos_2FibClust_hodoY1>-999) && (t.nClusters_hodoY1==1||t.pos_2FibClust_hodoY1>-999)){//exactly one cluster, or, if there are multiple clusters, exactly one 2-fiber cluster
	if(TMath::Abs(0.5* (t.cluster_pos_corr_hodoX1+t.cluster_pos_corr_hodoX2))< 3 && TMath::Abs( 0.5* (t.cluster_pos_corr_hodoY1+t.cluster_pos_corr_hodoY2))< 3 && (t.wc_x_corr-t.cluster_pos_corr_hodoX2)<4 && (t.wc_y_corr-t.cluster_pos_corr_hodoY2)< 4  && TMath::Abs( (t.cluster_pos_corr_hodoX1-t.cluster_pos_corr_hodoX2))<1.5 &&TMath::Abs( (t.cluster_pos_corr_hodoY1-t.cluster_pos_corr_hodoY2))<1.5){
	    totalHistos_tight_maxAmpl_fit_tot[i]->Fill(t.cef3_maxAmpl_fit_corr->at(0)+t.cef3_maxAmpl_fit_corr->at(1)+(i*0.1)*t.cef3_maxAmpl_fit_corr->at(2)+t.cef3_maxAmpl_fit_corr->at(3));
	}
      }
     }
     std::cout<<"reso:"<<totalHistos_tight_maxAmpl_fit_tot[i]->GetRMS()/totalHistos_tight_maxAmpl_fit_tot[i]->GetMean()<<std::endl;
     totalHistos_tight_maxAmpl_fit_tot[i]->Write();

     //fitting with Cruijff
    TH1F* histo;
    histo=totalHistos_tight_maxAmpl_fit_tot[i];
  
    double peakpos = histo->GetMean();
    double sigma = histo->GetRMS();
  
    double fitmin;
    double fitmax;
  
  
    fitmin = peakpos-5*sigma;
    fitmax = peakpos+5*sigma;
        
    RooRealVar x("x","MaxAmpl", fitmin, fitmax);
    RooDataHist data("data","dataset with x",x,RooFit::Import(*histo) );

    RooPlot* frame;
    frame = x.frame("Title");
    data.plotOn(frame);  //this will show histogram data points on canvas                                                                                                     
    std::cout<<"#######################mean"<<peakpos<<" sigma"<<sigma<<std::endl;

    RooRealVar meanr("meanr","Mean",peakpos,peakpos-2*sigma, peakpos+2*sigma);
    RooRealVar widthL("widthL","#sigmaL",sigma , 0, 5*sigma);
    RooRealVar widthR("widthR","#sigmaR",sigma , 0, 5*sigma);
    RooRealVar alphaL("alphaL","#alpha",1.08615e-02 , 0., 1.);
    RooRealVar alphaR("alphaR","#alpha",1.08615e-02 , 0., 1.);


    RooCruijff fit_fct("fit_fct","fit_fct",x,meanr,widthL,widthR,alphaL,alphaR); int ndf = 5;
    fit_fct.fitTo(data);
    fit_fct.plotOn(frame,RooFit::LineColor(4));//this will show fit overlay on canvas       
    
    double mean = meanr.getVal();
    double meanErr = meanr.getError();
    double rms = (widthL.getVal()+widthR.getVal())/2;
    double rmsErr = 0.5*sqrt(widthL.getError()*widthL.getError()+widthR.getError()*widthR.getError());
    double reso = 100.* rms/mean; //in percent                          
    double resoErr = 100.* getRatioError( rms, mean, meanErr, rmsErr );
    
    
    TCanvas* cans = new TCanvas("cans", "un canvas", 600,600);
    cans->cd();
    frame->Draw();
    TLegend* lego = new TLegend(0.6, 0.7, 0.9, 0.92);
    lego->SetTextSize(0.038);
    lego->AddEntry(  (TObject*)0 ,Form("#mu = %.0f #pm %.0f", meanr.getVal(), meanr.getError() ), "");
    lego->AddEntry(  (TObject*)0 ,Form("#sigma = %.0f #pm %.0f ", rms,  rmsErr), "");
    lego->AddEntry(  (TObject*)0 ,Form("#chi^{2} = %.2f / %d ", frame->chiSquare(ndf) , ndf ), "");
    lego->AddEntry(  (TObject*)0 ,Form("#sigma/#mu = %.1f #pm %.1f %s ", reso , resoErr ,"%"), "");
    lego->SetFillColor(0);
    lego->Draw("same");

    if(resoErr<0.2){//avoid points where fit didn't converge
      resValue[i]=reso;
      resErrValue[i]=resoErr;
    }

     cans->Write("c_"+step);

   }

   TGraphErrors* resGraph = new TGraphErrors(cValue,resValue,cValueErr,resErrValue);
   TCanvas c1;
   resGraph->SetMarkerColor(kRed);
   resGraph->SetLineColor(kRed);
   resGraph->SetFillColor(kRed);
   resGraph->SetMarkerStyle(kRed);
   resGraph->GetXaxis()->SetTitle("C_{2}");
   resGraph->GetYaxis()->SetTitle("resolution [%]");
   resGraph->GetYaxis()->SetRangeUser(2.3,4.8);
   resGraph->Draw("AP");
   c1.SaveAs("resScan.png");
   c1.SaveAs("resScan.pdf");
   resGraph->Write();

   resValue.Write("resValue");
   resErrValue.Write("resValueErr");
   outFile->Write();
   outFile->Close();

}
