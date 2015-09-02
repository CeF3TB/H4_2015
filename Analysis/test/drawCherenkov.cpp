
#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>


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

   } else {

     std::cout << "Usage:" << std::endl;
     std::cout << "./drawCherenkov [runName] ([tag])" << std::endl;
     exit(12345);

   }

  TString runNumberString(runName);

  std::string fileName = "analysisTrees_"+tag+"/Reco_"+runName+".root";
  TFile* file = TFile::Open(fileName.c_str());
  if( file==0 ) {
    std::cout << "ERROR! Din't find file " << fileName << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }
  


  TH1F* totalHistos[4];
  TH1F* totalHistosGainCorr[4];
  TH1F* cherHistos[4];
  TH1F* wlsHistos[4];

  TH1F* totalHistos_tight[4];
  TH1F* totalHistosGainCorr_tight[4];
  TH1F* cherHistos_tight[4];
  TH1F* wlsHistos_tight[4];
  TH1F* wlsHistos_tight_tot;

  TH1F* totalHistos_tight_maxAmpl[4];
  TH1F* totalHistos_tight_maxAmpl_fit[4];

  TH1F* totalHistos_tight_maxAmpl_tot;
  TH1F* totalHistos_tight_maxAmpl_fit_tot;

  for(int i=0;i<CEF3_CHANNELS;++i){
    TString fibre;
    fibre.Form("%d",i); 
    totalHistosGainCorr[i] = new TH1F ("totalHistoGainCorr_"+fibre,"",200,0,0.5);
    totalHistos[i] = new TH1F ("totalHisto_"+fibre,"",200,0,300000);
    cherHistos[i] = new TH1F ("cherHisto_"+fibre,"",200,0,55e3);
    wlsHistos[i] = new TH1F ("wlsHisto_"+fibre,"",200,0,500000);
    

    totalHistosGainCorr_tight[i] = new TH1F ("totalHistoGainCorr_tight_"+fibre,"",200,0,0.5);
    totalHistos_tight[i] = new TH1F ("totalHisto_tight_"+fibre,"",400,0,300000);
    cherHistos_tight[i] = new TH1F ("cherHisto_tight_"+fibre,"",200,0,55e3);
    wlsHistos_tight[i] = new TH1F ("wlsHisto_tight_"+fibre,"",200,0,500000);

    totalHistos_tight_maxAmpl[i] = new TH1F ("totalHisto_tight_maxAmpl"+fibre,"",700,0,4000);
    totalHistos_tight_maxAmpl_fit[i] = new TH1F ("totalHisto_tight_maxAmpl_fit"+fibre,"",700,0,4000);
    
  }

  wlsHistos_tight_tot= new TH1F ("wlsHisto_tight_total","",200*4,0,500000*4);
  totalHistos_tight_maxAmpl_tot= new TH1F ("totalHisto_tight_maxAmpl_total","",200*4,0,4000*4);
  totalHistos_tight_maxAmpl_fit_tot= new TH1F ("totalHisto_tight_maxAmpl_fit_total","",200*4,0,4000*4);

  TTree* recoTree=(TTree*)file->Get("recoTree");
  RecoTree t(recoTree);

   if (t.fChain == 0) return 0;
   
   Long64_t nentries = t.fChain->GetEntries();

   float gainR1450=1.5e+6;
   float gainR5380=7e+4;
   
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = t.LoadTree(jentry);
      if (ientry < 0) break;
      nb = t.fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      for(int i=0;i<CEF3_CHANNELS;++i){
	totalHistos[i]->Fill(t.cef3_chaInt->at(i));
	if(i==2)	totalHistosGainCorr[i]->Fill(t.cef3_chaInt->at(i)/gainR5380);
	else totalHistosGainCorr[i]->Fill(t.cef3_chaInt->at(i)/gainR1450);
	cherHistos[i]->Fill(t.cef3_chaInt_cher->at(i));
	wlsHistos[i]->Fill(t.cef3_chaInt_wls->at(i));

	if((t.nClusters_hodoX1==1||t.pos_2FibClust_hodoX1>-999) && (t.nClusters_hodoX2==1||t.pos_2FibClust_hodoX2>-999) && (t.nClusters_hodoY1==1||t.pos_2FibClust_hodoY1>-999) && (t.nClusters_hodoY1==1||t.pos_2FibClust_hodoY1>-999)){//exactly one cluster, or, if there are multiple clusters, exactly one 2-fiber cluster
	  if(TMath::Abs(0.5* (t.cluster_pos_corr_hodoX1+t.cluster_pos_corr_hodoX2))< 3 && TMath::Abs( 0.5* (t.cluster_pos_corr_hodoY1+t.cluster_pos_corr_hodoY2))< 3 && (t.wc_x_corr-t.cluster_pos_corr_hodoX2)<4 && (t.wc_y_corr-t.cluster_pos_corr_hodoY2)< 4  && TMath::Abs( (t.cluster_pos_corr_hodoX1-t.cluster_pos_corr_hodoX2))<1.5 &&TMath::Abs( (t.cluster_pos_corr_hodoY1-t.cluster_pos_corr_hodoY2))<1.5){
	  totalHistos_tight[i]->Fill(t.cef3_chaInt->at(i)); 
	  if(i==2)	totalHistosGainCorr_tight[i]->Fill(t.cef3_chaInt->at(i)/gainR5380);
	  else totalHistosGainCorr_tight[i]->Fill(t.cef3_chaInt->at(i)/gainR1450);
	  cherHistos_tight[i]->Fill(t.cef3_chaInt_cher->at(i));
	  wlsHistos_tight[i]->Fill(t.cef3_chaInt_wls->at(i));
	  if(i==0){
	    wlsHistos_tight_tot->Fill(t.cef3_chaInt_wls->at(0)+t.cef3_chaInt_wls->at(1)+t.cef3_chaInt->at(2)+t.cef3_chaInt_wls->at(3));
	    totalHistos_tight_maxAmpl_tot->Fill(t.cef3_maxAmpl->at(0)+t.cef3_maxAmpl->at(1)+t.cef3_maxAmpl->at(2)+t.cef3_maxAmpl->at(3));
	    //	    totalHistos_tight_maxAmpl_fit_tot->Fill(t.cef3_maxAmpl_fit->at(0)+t.cef3_maxAmpl_fit->at(1)+t.cef3_maxAmpl_fit->at(2)+t.cef3_maxAmpl_fit->at(3));
	    totalHistos_tight_maxAmpl_fit_tot->Fill(t.cef3_maxAmpl_fit_corr->at(0)+t.cef3_maxAmpl_fit_corr->at(1)+t.cef3_maxAmpl_fit_corr->at(2)+t.cef3_maxAmpl_fit_corr->at(3));
	  }
	  totalHistos_tight_maxAmpl[i]->Fill(t.cef3_maxAmpl->at(i)); 
	  totalHistos_tight_maxAmpl_fit[i]->Fill(t.cef3_maxAmpl_fit_corr->at(i)); 
	  }
	}
      }
   }




  TGaxis::SetMaxDigits(3);



  TFile* outFile = TFile::Open("CherenkovPlots_"+runNumberString+"_"+tag+".root","recreate");

  //plots of ch int
  TPaveText * pave= new TPaveText(0.75,0.65,0.85,0.85,"NDC");
  pave->SetFillColor(kWhite);
  pave->SetTextSize(0.030);
  pave->SetTextAlign(kHAlignLeft);
  pave->SetTextFont(62);

  TCanvas c1;
  c1.SetLogy();
  totalHistos[1]->SetLineColor(kBlue);
  totalHistos[2]->SetLineColor(kRed);
  totalHistos[3]->SetLineColor(kViolet);

  float max=-1;
  for(int i=0;i<4;++i){
      if(max<totalHistos[i]->GetMaximum()){
      max=totalHistos[i]->GetMaximum();
    }
  }


  for(int i=0;i<4;++i){
    totalHistos[i]->GetYaxis()->SetRangeUser(1,max*1.1);
    totalHistos[i]->GetXaxis()->SetTitle("Charge Integrated");
    totalHistos[i]->GetYaxis()->SetTitle("Events");
    if (i==0)    totalHistos[i]->Draw();
    else totalHistos[i]->Draw("same");
    TString fibre;
    fibre.Form("%d",i); 
    pave->AddText("Fibre "+fibre);  
    pave->SetAllWith("Fibre "+fibre,"color",totalHistos[i]->GetLineColor());
  }
  pave->Draw("same");
  c1.SaveAs("plots_drawCherenkov/chargeIntTotal_"+runNumberString+".png");
  c1.SaveAs("plots_drawCherenkov/chargeIntTotal_"+runNumberString+".pdf");
  c1.Write("chargeIntTotal");


  c1.Clear();
  totalHistosGainCorr[1]->SetLineColor(kBlue);
  totalHistosGainCorr[2]->SetLineColor(kRed);
  totalHistosGainCorr[3]->SetLineColor(kViolet);
  max=-1;
  for(int i=0;i<4;++i){
    if(max<totalHistosGainCorr[i]->GetMaximum()){
      max=totalHistosGainCorr[i]->GetMaximum();
    }
  }
  
  for(int i=0;i<4;++i){
    totalHistosGainCorr[i]->GetYaxis()->SetRangeUser(1,max*1.1);
    totalHistosGainCorr[i]->GetXaxis()->SetTitle("Charge Integrated");
    totalHistosGainCorr[i]->GetYaxis()->SetTitle("Events");
    if (i==0)    totalHistosGainCorr[i]->Draw();
    else totalHistosGainCorr[i]->Draw("same");
  }
  pave->Draw("same");
  c1.SaveAs("plots_drawCherenkov/chargeIntTotalGainCorr_"+runNumberString+".png");
  c1.SaveAs("plots_drawCherenkov/chargeIntTotalGainCorr_"+runNumberString+".pdf");
  c1.Write("chargeIntTotalGainCorr");


  c1.Clear();
  cherHistos[1]->SetLineColor(kBlue);
  cherHistos[2]->SetLineColor(kRed);
  cherHistos[3]->SetLineColor(kViolet);

  max=-1;
  for(int i=0;i<4;++i){
    if(max<cherHistos[i]->GetMaximum()){
      max=cherHistos[i]->GetMaximum();
    }
  }


  for(int i=0;i<4;++i){
    cherHistos[i]->GetYaxis()->SetRangeUser(1,max*1.1);
    cherHistos[i]->GetXaxis()->SetTitle("Charge Integrated");
    cherHistos[i]->GetYaxis()->SetTitle("Events");
    if (i==0)    cherHistos[i]->Draw();
    else cherHistos[i]->Draw("same");
  }
  pave->Draw("same");
  c1.SaveAs("plots_drawCherenkov/chargeIntCher_"+runNumberString+".png");
  c1.SaveAs("plots_drawCherenkov/chargeIntCher_"+runNumberString+".pdf");
  c1.Write("chargeIntCher");

  c1.Clear();
  wlsHistos[1]->SetLineColor(kBlue);
  wlsHistos[2]->SetLineColor(kRed);
  wlsHistos[3]->SetLineColor(kViolet);
  max=-1;
  for(int i=0;i<4;++i){
    if(max<wlsHistos[i]->GetMaximum()){
      max=wlsHistos[i]->GetMaximum();
    }
  }

  for(int i=0;i<4;++i){
    wlsHistos[i]->GetYaxis()->SetRangeUser(1,max*1.1);
    wlsHistos[i]->GetXaxis()->SetTitle("Charge Integrated");
    wlsHistos[i]->GetYaxis()->SetTitle("Events");
    if (i==0)    wlsHistos[i]->Draw();
    else wlsHistos[i]->Draw("same");
  }
  pave->Draw("same");
  c1.SaveAs("plots_drawCherenkov/chargeIntWls_"+runNumberString+".png");
  c1.SaveAs("plots_drawCherenkov/chargeIntWls_"+runNumberString+".pdf");
  c1.Write("chargeIntWls");


  for(int i=0;i<4;++i){
    totalHistos[i]->Write();
    totalHistosGainCorr[i]->Write();
    cherHistos[i] ->Write();
    wlsHistos[i] ->Write();
  }


  //tight sel
  c1.Clear();
  c1.SetLogy();
  totalHistos_tight[1]->SetLineColor(kBlue);
  totalHistos_tight[2]->SetLineColor(kRed);
  totalHistos_tight[3]->SetLineColor(kViolet);

  totalHistos_tight_maxAmpl[1]->SetLineColor(kBlue);
  totalHistos_tight_maxAmpl[2]->SetLineColor(kRed);
  totalHistos_tight_maxAmpl[3]->SetLineColor(kViolet);

  totalHistos_tight_maxAmpl_fit[1]->SetLineColor(kBlue);
  totalHistos_tight_maxAmpl_fit[2]->SetLineColor(kRed);
  totalHistos_tight_maxAmpl_fit[3]->SetLineColor(kViolet);


  max=-1;
  float maxAmpl=-1;
  for(int i=0;i<4;++i){
    if(max<totalHistos_tight[i]->GetMaximum()){
      max=totalHistos_tight[i]->GetMaximum();
    }
    if(maxAmpl<totalHistos_tight_maxAmpl[i]->GetMaximum()){
      maxAmpl= totalHistos_tight_maxAmpl[i]->GetMaximum();  
    }
  }


  for(int i=0;i<4;++i){
    totalHistos_tight[i]->GetYaxis()->SetRangeUser(1,max*1.1);
    totalHistos_tight[i]->GetXaxis()->SetTitle("Charge Integrated");
    totalHistos_tight[i]->GetYaxis()->SetTitle("Events");
    if (i==0)    totalHistos_tight[i]->Draw();
    else totalHistos_tight[i]->Draw("same");
    TString fibre;
  }
  pave->Draw("same");
  c1.SaveAs("plots_drawCherenkov/chargeIntTotal_tight_"+runNumberString+".png");
  c1.SaveAs("plots_drawCherenkov/chargeIntTotal_tight_"+runNumberString+".pdf");
  c1.Write("chargeIntTotal_tight");


  c1.Clear();

  //maxAmpl
  for(int i=0;i<4;++i){
    totalHistos_tight_maxAmpl[i]->GetXaxis()->SetRangeUser(0,2000);
    totalHistos_tight_maxAmpl[i]->GetYaxis()->SetRangeUser(1,max*1.1);
    totalHistos_tight_maxAmpl[i]->GetXaxis()->SetTitle("Max Amplitude");
    totalHistos_tight_maxAmpl[i]->GetYaxis()->SetTitle("Events");
    if (i==0)    totalHistos_tight_maxAmpl[i]->Draw();
    else totalHistos_tight_maxAmpl[i]->Draw("same");
    TString fibre;
  }
  pave->Draw("same");
  c1.SaveAs("plots_drawCherenkov/maxAmpl_tight_"+runNumberString+".png");
  c1.SaveAs("plots_drawCherenkov/maxAmpl_tight_"+runNumberString+".pdf");
  c1.Write("maxAmpl_tight");

  c1.Clear();

  //maxAmpl_fit
  for(int i=0;i<4;++i){
    totalHistos_tight_maxAmpl_fit[i]->GetXaxis()->SetRangeUser(0,2000);
    totalHistos_tight_maxAmpl_fit[i]->GetYaxis()->SetRangeUser(1,max*1.1);
    totalHistos_tight_maxAmpl_fit[i]->GetXaxis()->SetTitle("Max Amplitude");
    totalHistos_tight_maxAmpl_fit[i]->GetYaxis()->SetTitle("Events");
    if (i==0)    totalHistos_tight_maxAmpl_fit[i]->Draw();
    else totalHistos_tight_maxAmpl_fit[i]->Draw("same");
    TString fibre;
  }
  pave->Draw("same");
  c1.SaveAs("plots_drawCherenkov/maxAmpl_fit_tight_"+runNumberString+".png");
  c1.SaveAs("plots_drawCherenkov/maxAmpl_fit_tight_"+runNumberString+".pdf");
  c1.Write("maxAmpl_fit_tight");

  c1.Clear();



  totalHistosGainCorr_tight[1]->SetLineColor(kBlue);
  totalHistosGainCorr_tight[2]->SetLineColor(kRed);
  totalHistosGainCorr_tight[3]->SetLineColor(kViolet);
  max=-1;
  for(int i=0;i<4;++i){
    if(max<totalHistosGainCorr_tight[i]->GetMaximum()){
      max=totalHistosGainCorr_tight[i]->GetMaximum();
    }
  }
  
  for(int i=0;i<4;++i){
    totalHistosGainCorr_tight[i]->GetYaxis()->SetRangeUser(1,max*1.1);
    totalHistosGainCorr_tight[i]->GetXaxis()->SetTitle("Charge Integrated");
    totalHistosGainCorr_tight[i]->GetYaxis()->SetTitle("Events");
    if (i==0)    totalHistosGainCorr_tight[i]->Draw();
    else totalHistosGainCorr_tight[i]->Draw("same");
  }
  pave->Draw("same");
  c1.SaveAs("plots_drawCherenkov/chargeIntTotalGainCorr_tight_"+runNumberString+".png");
  c1.SaveAs("plots_drawCherenkov/chargeIntTotalGainCorr_tight_"+runNumberString+".pdf");
  c1.Write("chargeIntTotal_tight");


  c1.Clear();
  cherHistos_tight[1]->SetLineColor(kBlue);
  cherHistos_tight[2]->SetLineColor(kRed);
  cherHistos_tight[3]->SetLineColor(kViolet);

  max=-1;
  for(int i=0;i<4;++i){
    if(max<cherHistos_tight[i]->GetMaximum()){
      max=cherHistos_tight[i]->GetMaximum();
    }
  }


  for(int i=0;i<4;++i){
    cherHistos_tight[i]->GetYaxis()->SetRangeUser(1,max*1.1);
    cherHistos_tight[i]->GetXaxis()->SetTitle("Charge Integrated");
    cherHistos_tight[i]->GetYaxis()->SetTitle("Events");
    if (i==0)    cherHistos_tight[i]->Draw();
    else cherHistos_tight[i]->Draw("same");
  }
  pave->Draw("same");
  c1.SaveAs("plots_drawCherenkov/chargeIntCher_tight_"+runNumberString+".png");
  c1.SaveAs("plots_drawCherenkov/chargeIntCher_tight_"+runNumberString+".pdf");
  c1.Write("chargeIntCher_tight");

  c1.Clear();
  wlsHistos_tight[1]->SetLineColor(kBlue);
  wlsHistos_tight[2]->SetLineColor(kRed);
  wlsHistos_tight[3]->SetLineColor(kViolet);
  max=-1;
  for(int i=0;i<4;++i){
    if(max<wlsHistos_tight[i]->GetMaximum()){
      max=wlsHistos_tight[i]->GetMaximum();
    }
  }

  for(int i=0;i<4;++i){
    wlsHistos_tight[i]->GetYaxis()->SetRangeUser(1,max*1.1);
    wlsHistos_tight[i]->GetXaxis()->SetTitle("Charge Integrated");
    wlsHistos_tight[i]->GetYaxis()->SetTitle("Events");
    if (i==0)    wlsHistos_tight[i]->Draw();
    else wlsHistos_tight[i]->Draw("same");
  }
  pave->Draw("same");
  c1.SaveAs("plots_drawCherenkov/chargeIntWls_tight_"+runNumberString+".png");
  c1.SaveAs("plots_drawCherenkov/chargeIntWls_tight_"+runNumberString+".pdf");
  c1.Write("chargeIntWls_tight");


  for(int i=0;i<4;++i){
    totalHistos_tight[i]->Write();
    totalHistosGainCorr_tight[i]->Write();
    cherHistos_tight[i] ->Write();
    wlsHistos_tight[i] ->Write();
    totalHistos_tight_maxAmpl[i]->Write();
    totalHistos_tight_maxAmpl_fit[i]->Write();
  }

  //get the resolution of chInt and maxAmpl
  TVectorD meanValue(5);
  TVectorD meanErrValue(5);
  TVectorD widthValue(5);
  TVectorD widthErrValue(5);
  TVectorD resValue(5);
  TVectorD resErrValue(5);

    
  for(int i=0;i<4;++i){
    TH1F* histo;
    if(i!=2) histo=wlsHistos_tight[i];
    else histo=totalHistos_tight[i];
    //else histo=wlsHistos_tight[i];

    double peakpos = histo->GetMean();
    double sigma = histo->GetRMS();
    
    double fitmin;
    double fitmax;
    
    
    fitmin = peakpos-5*sigma;
    fitmax = peakpos+5*sigma;
        
    RooRealVar x("x","ChInt", fitmin, fitmax);
    RooDataHist data("data","dataset with x",x,RooFit::Import(*histo) );

    RooPlot* frame;
    frame = x.frame("Title");
    data.plotOn(frame);  //this will show histogram data points on canvas                                                                                                         
    //  data.statOn(frame);  //this will display hist stat on canvas                                                                                                                  
    //    RooRealVar meanr("meanr","Mean",peakpos-1.5*sigma,peakpos-3*sigma, peakpos+2*sigma);
    RooRealVar meanr("meanr","Mean",peakpos-sigma,peakpos-3*sigma, peakpos+2*sigma);
    //    if(i==0)meanr.setConstant(kTRUE);
    RooRealVar width("width","#sigma",sigma , 150.0, 5.*sigma);
    RooRealVar A("A","Dist",2., 0.0, 7.0);
    RooRealVar N("N","Deg",5, 0.0, 10);
    
    RooRealVar widthL("widthL","#sigmaL",sigma , 2, 5*sigma);
    RooRealVar widthR("widthR","#sigmaR",sigma , 2, 5*sigma);
    RooRealVar alphaL("alphaL","#alpha",5.08615e-02 , 0., 1.);
    RooRealVar alphaR("alphaR","#alpha",0, 0., 1.);
    int ndf;

    //  meanr.setRange( 30000. , 1000000.);
    //  width.setRange(500, 22000);
    bool fitWithCB=false;
    if(!fitWithCB){
      RooCruijff fit_fct("fit_fct","fit_fct",x,meanr,widthL,widthR,alphaL,alphaR); ndf = 5;
      fit_fct.fitTo(data);
      fit_fct.plotOn(frame,RooFit::LineColor(4));//this will show fit overlay on canvas
    }else{
      RooCBShape fit_fct("fit_fct","fit_fct",x,meanr,width,A,N); ndf = 4;
      fit_fct.fitTo(data);
      fit_fct.plotOn(frame,RooFit::LineColor(4));//this will show fit overlay on canvas
    }    
      
    TString fibre;
    fibre.Form("%d",i); 

    TH1F* fittedHisto=(TH1F*)data.createHistogram("histo_fit_"+fibre,x);
    fittedHisto->Write();
        // fit_fct.paramOn(frame); //this will display the fit parameters on canvas               
    
    double mean = meanr.getVal();
    double meanErr = meanr.getError();
    double rms,rmsErr;
    if(!fitWithCB){
      rms = (widthL.getVal()+widthR.getVal())/2;
      rmsErr = 0.5*sqrt(widthL.getError()*widthL.getError()+widthR.getError()*widthR.getError());
    }else{
      rms = width.getVal();
      rmsErr = width.getError();
    }
    double reso = 100.* rms/mean; //in percent                          
    double resoErr = 100.* getRatioError( rms, mean, meanErr, rmsErr );
    
    
    TCanvas* cans = new TCanvas("cans", "un canvas", 600,600);
    cans->cd();
    frame->Draw();
    TLegend* lego = new TLegend(0.6, 0.7, 0.9, 0.92);
    lego->SetTextSize(0.038);
    lego->AddEntry(  (TObject*)0 ,Form("#mu = %.0f #pm %.0f", meanr.getVal(), meanr.getError() ), "");
    lego->AddEntry(  (TObject*)0 ,Form("#sigma = %.0f #pm %.0f ", rms, rmsErr), "");
    lego->AddEntry(  (TObject*)0 ,Form("#chi^{2} = %.2f / %d ", frame->chiSquare(ndf) , ndf ), "");
    lego->AddEntry(  (TObject*)0 ,Form("#sigma/#mu = %.1f #pm %.1f %s ", reso , resoErr ,"%"), "");
    lego->SetFillColor(0);
    lego->Draw("same");

    meanValue[i]=meanr.getVal();
    meanErrValue[i]=meanr.getError();

    widthValue[i]=width.getVal();
    widthErrValue[i]=width.getError();

    resValue[i]=reso;
    resErrValue[i]=resoErr;
    cans->SaveAs("plots_drawCherenkov/CBFit_"+runNumberString+"_fibre_"+fibre+".png");
    cans->SaveAs("plots_drawCherenkov/CBFit_"+runNumberString+"_fibre_"+fibre+".pdf");
    cans->Write();
  }

  //maxAmplFit
  TVectorD meanValuemaxAmpl(5);
  TVectorD meanErrValuemaxAmpl(5);
  TVectorD widthValuemaxAmpl(5);
  TVectorD widthErrValuemaxAmpl(5);
  TVectorD resValuemaxAmpl(5);
  TVectorD resErrValuemaxAmpl(5);

  TVectorD meanValuemaxAmpl_fit(5);
  TVectorD meanErrValuemaxAmpl_fit(5);
  TVectorD widthValuemaxAmpl_fit(5);
  TVectorD widthErrValuemaxAmpl_fit(5);
  TVectorD resValuemaxAmpl_fit(5);
  TVectorD resErrValuemaxAmpl_fit(5);

    
  for(int i=0;i<4;++i){
    TH1F* histo;
    histo=totalHistos_tight_maxAmpl[i];

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
    //    std::cout<<"#######################mean"<<peakpos<<" sigma"<<sigma<<std::endl;
    RooRealVar meanr("meanr","Mean",peakpos,peakpos-2*sigma, peakpos+2*sigma);
    //    if(i==0)meanr.setConstant(kTRUE);
    
      // meanr.setRange( 30000. , 1000000.);
      // width.setRange(500, 22000);
    //    RooCBShape fit_fct("fit_fct","fit_fct",x,meanr,width,A,N); int ndf = 4;
    

    RooRealVar widthL("widthL","#sigmaL",sigma , 0, 5*sigma);
    RooRealVar widthR("widthR","#sigmaR",sigma , 0, 5*sigma);
    RooRealVar alphaL("alphaL","#alpha",1.08615e-02 , 0., 1.);
    RooRealVar alphaR("alphaR","#alpha",1.08615e-02 , 0., 1.);


    RooCruijff fit_fct("fit_fct","fit_fct",x,meanr,widthL,widthR,alphaL,alphaR); int ndf = 5;
    fit_fct.fitTo(data);
    fit_fct.plotOn(frame,RooFit::LineColor(4));//this will show fit overlay on canvas       

    TString fibre;
    fibre.Form("%d",i); 

    TH1F* fittedHisto=(TH1F*)data.createHistogram("histo_fit_maxAmpl"+fibre,x);
    fittedHisto->Write();
        // fit_fct.paramOn(frame); //this will display the fit parameters on canvas               
    
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

    meanValuemaxAmpl[i]=meanr.getVal();
    meanErrValuemaxAmpl[i]=meanr.getError();

    widthValuemaxAmpl[i]=rms;
    widthErrValuemaxAmpl[i]=rmsErr;

    resValuemaxAmpl[i]=reso;
    resErrValuemaxAmpl[i]=resoErr;

    cans->SaveAs("plots_drawCherenkov/CBFit_maxAmpl_"+runNumberString+"_fibre_"+fibre+".png");
    cans->SaveAs("plots_drawCherenkov/CBFit_maxAmpl_"+runNumberString+"_fibre_"+fibre+".pdf");
    cans->Write();
    //        if (i>1)exit(9);
  }
  //end max ampl

  //maxAmpl_fit
  for(int i=0;i<4;++i){
    TH1F* histo;
    histo=totalHistos_tight_maxAmpl_fit[i];

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
    RooRealVar meanr("meanr","Mean",peakpos-sigma,peakpos-3*sigma, peakpos+3*sigma);
    RooRealVar widthL("widthL","#sigmaL",sigma , 2, 5*sigma);
    RooRealVar widthR("widthR","#sigmaR",sigma , 2, 5*sigma);
    RooRealVar alphaL("alphaL","#alpha",5.08615e-02 , 0., 1.);
    RooRealVar alphaR("alphaR","#alpha",0, 0., 1.);
    int ndf;

    RooRealVar width("width","#sigma",sigma, 0, 5.*sigma);
    RooRealVar A("A","Dist",2., 0.0, 7.0);
    RooRealVar N("N","Deg",5, 0.0, 10);

    bool fitWithCB=false;
    if(!fitWithCB){
      RooCruijff fit_fct("fit_fct","fit_fct",x,meanr,widthL,widthR,alphaL,alphaR); ndf = 5;
      fit_fct.fitTo(data);
      fit_fct.plotOn(frame,RooFit::LineColor(4));//this will show fit overlay on canvas  
    }else{
      RooCBShape fit_fct("fit_fct","fit_fct",x,meanr,width,A,N); ndf = 4;
      fit_fct.fitTo(data);
      fit_fct.plotOn(frame,RooFit::LineColor(4));//this will show fit overlay on canvas  
    }    
    TString fibre;
    fibre.Form("%d",i); 

    TH1F* fittedHisto=(TH1F*)data.createHistogram("histo_fit_maxAmpl_fit"+fibre,x);
    fittedHisto->Write();
        // fit_fct.paramOn(frame); //this will display the fit parameters on canvas               
    
    double mean = meanr.getVal();
    double meanErr = meanr.getError();
    double rms,rmsErr;
    if(!fitWithCB){
      rms = (widthL.getVal()+widthR.getVal())/2;
      rmsErr = 0.5*sqrt(widthL.getError()*widthL.getError()+widthR.getError()*widthR.getError());
    }else{
      rms = width.getVal();
      rmsErr = width.getError();
    }
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

    meanValuemaxAmpl_fit[i]=meanr.getVal();
    meanErrValuemaxAmpl_fit[i]=meanr.getError();

    widthValuemaxAmpl_fit[i]=rms;
    widthErrValuemaxAmpl_fit[i]=rmsErr;

    resValuemaxAmpl_fit[i]=reso;
    resErrValuemaxAmpl_fit[i]=resoErr;

    cans->SaveAs("plots_drawCherenkov/CBFit_maxAmpl_fit_"+runNumberString+"_fibre_"+fibre+".png");
    cans->SaveAs("plots_drawCherenkov/CBFit_maxAmpl_fit_"+runNumberString+"_fibre_"+fibre+".pdf");
    cans->Write();
  }
  //end max ampl
  

  //total res
  if(wlsHistos_tight_tot){
    TH1F* histo;
    histo=wlsHistos_tight_tot;
    
    
    double peakpos = histo->GetMean();
    double sigma = histo->GetRMS();
  
    double fitmin;
    double fitmax;
  
  
    fitmin = peakpos-5*sigma;
    fitmax = peakpos+5*sigma;
        
    RooRealVar x("x","ChInt", fitmin, fitmax);
    RooDataHist data("data","dataset with x",x,RooFit::Import(*histo) );
  
    RooPlot* frame;
    frame = x.frame("Title");
    data.plotOn(frame);  //this will show histogram data points on canvas                                                                                                         
    //  data.statOn(frame);  //this will display hist stat on canvas                                                                                                                  
    RooRealVar meanr("meanr","Mean",peakpos,peakpos-2*sigma, peakpos+2*sigma);
    RooRealVar widthL("widthL","#sigmaL",sigma , 0, 5*sigma);
    RooRealVar widthR("widthR","#sigmaR",sigma , 0, 5*sigma);
    RooRealVar alphaL("alphaL","#alpha",1.08615e-02 , 0., 1.);
    RooRealVar alphaR("alphaR","#alpha",1.08615e-02 , 0., 1.);


    RooCruijff fit_fct("fit_fct","fit_fct",x,meanr,widthL,widthR,alphaL,alphaR); int ndf = 5;
    fit_fct.fitTo(data);
    fit_fct.plotOn(frame,RooFit::LineColor(4));//this will show fit overlay on canvas       

  
    TH1F* fittedHisto=(TH1F*)data.createHistogram("histo_fit_total",x);
    fittedHisto->Write();
    // fit_fct.paramOn(frame); //this will display the fit parameters on canvas               
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
  
    meanValue[4]=meanr.getVal();
    meanErrValue[4]=meanr.getError();
  
    widthValue[4]=rms;
    widthErrValue[4]=rmsErr;
  
    resValue[4]=reso;
    resErrValue[4]=resoErr;
    
  
    cans->SaveAs("plots_drawCherenkov/CBFit_"+runNumberString+"_total.png");
    cans->SaveAs("plots_drawCherenkov/CBFit_"+runNumberString+"_total.pdf");
    cans->Write();
  }

  //total res max ampl
  if(totalHistos_tight_maxAmpl_tot){
    TH1F* histo;
    histo=totalHistos_tight_maxAmpl_tot;
  
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
    //    if(i==0)meanr.setConstant(kTRUE);
    
      // meanr.setRange( 30000. , 1000000.);
      // width.setRange(500, 22000);
    //    RooCBShape fit_fct("fit_fct","fit_fct",x,meanr,width,A,N); int ndf = 4;
    

    RooRealVar widthL("widthL","#sigmaL",sigma , 0, 5*sigma);
    RooRealVar widthR("widthR","#sigmaR",sigma , 0, 5*sigma);
    RooRealVar alphaL("alphaL","#alpha",1.08615e-02 , 0., 1.);
    RooRealVar alphaR("alphaR","#alpha",1.08615e-02 , 0., 1.);


    RooCruijff fit_fct("fit_fct","fit_fct",x,meanr,widthL,widthR,alphaL,alphaR); int ndf = 5;
    fit_fct.fitTo(data);
    fit_fct.plotOn(frame,RooFit::LineColor(4));//this will show fit overlay on canvas       
  
  
    TH1F* fittedHisto=(TH1F*)data.createHistogram("histo_fit_maxAmpl_tot",x);
    fittedHisto->Write();
    // fit_fct.paramOn(frame); //this will display the fit parameters on canvas               
    
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

    meanValuemaxAmpl[4]=meanr.getVal();
    meanErrValuemaxAmpl[4]=meanr.getError();
  
    widthValuemaxAmpl[4]=rms;
    widthErrValuemaxAmpl[4]=rmsErr;
  
    resValuemaxAmpl[4]=reso;
    resErrValuemaxAmpl[4]=resoErr;

    cans->SaveAs("plots_drawCherenkov/CBFit_maxAmpl_"+runNumberString+"_total.png");
    cans->SaveAs("plots_drawCherenkov/CBFit_maxAmpl_"+runNumberString+"_total.pdf");
    cans->Write();
  }
  //end total res max Ampl  
  //total res maxAmpl_fit
  //total res max ampl
  if(totalHistos_tight_maxAmpl_fit_tot){
    TH1F* histo;
    histo=totalHistos_tight_maxAmpl_fit_tot;
  
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
    //    if(i==0)meanr.setConstant(kTRUE);
    
      // meanr.setRange( 30000. , 1000000.);
      // width.setRange(500, 22000);
    //    RooCBShape fit_fct("fit_fct","fit_fct",x,meanr,width,A,N); int ndf = 4;
    

    RooRealVar widthL("widthL","#sigmaL",sigma , 0, 5*sigma);
    RooRealVar widthR("widthR","#sigmaR",sigma , 0, 5*sigma);
    RooRealVar alphaL("alphaL","#alpha",1.08615e-02 , 0., 1.);
    RooRealVar alphaR("alphaR","#alpha",1.08615e-02 , 0., 1.);


    RooCruijff fit_fct("fit_fct","fit_fct",x,meanr,widthL,widthR,alphaL,alphaR); int ndf = 5;
    fit_fct.fitTo(data);
    fit_fct.plotOn(frame,RooFit::LineColor(4));//this will show fit overlay on canvas       
  
  
    TH1F* fittedHisto=(TH1F*)data.createHistogram("histo_fit_maxAmpl_fit_tot",x);
    fittedHisto->Write();
    // fit_fct.paramOn(frame); //this will display the fit parameters on canvas               
    
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

    meanValuemaxAmpl_fit[4]=meanr.getVal();
    meanErrValuemaxAmpl_fit[4]=meanr.getError();
  
    widthValuemaxAmpl_fit[4]=rms;
    widthErrValuemaxAmpl_fit[4]=rmsErr;
  
    resValuemaxAmpl_fit[4]=reso;
    resErrValuemaxAmpl_fit[4]=resoErr;

    cans->SaveAs("plots_drawCherenkov/CBFit_maxAmpl_fit_"+runNumberString+"_total.png");
    cans->SaveAs("plots_drawCherenkov/CBFit_maxAmpl_fit_"+runNumberString+"_total.pdf");
    cans->Write();
  }
  //end total res maxAmpl_fit


  meanValue.Write("meanValue");
  meanErrValue.Write("meanErrValue");
  
  widthValue.Write("widthValue");
  widthErrValue.Write("widthErrValue");

  resValue.Write("resValue");
  resErrValue.Write("resErrValue");
  
  //maxAmpl
  meanValuemaxAmpl.Write("meanValuemaxAmpl");
  meanErrValuemaxAmpl.Write("meanErrValuemaxAmpl");
  
  widthValuemaxAmpl.Write("widthValuemaxAmpl");
  widthErrValuemaxAmpl.Write("widthErrValuemaxAmpl");

  resValuemaxAmpl.Write("resValuemaxAmpl");
  resErrValuemaxAmpl.Write("resErrValuemaxAmpl");
 
  //maxAmpl_fit
  meanValuemaxAmpl_fit.Write("meanValuemaxAmpl_fit");
  meanErrValuemaxAmpl_fit.Write("meanErrValuemaxAmpl_fit");
  
  widthValuemaxAmpl_fit.Write("widthValuemaxAmpl_fit");
  widthErrValuemaxAmpl_fit.Write("widthErrValuemaxAmpl_fit");

  resValuemaxAmpl_fit.Write("resValuemaxAmpl_fit");
  resErrValuemaxAmpl_fit.Write("resErrValuemaxAmpl_fit");
   

  outFile->Write();
  outFile->Close();



  return 0;

}
