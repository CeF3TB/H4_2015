#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <sstream> 

#include "DrawTools.h"
#include "RecoTree.h"
#include "channelInfo.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TVectorD.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TProfile.h"

#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooFFTConvPdf.h"
#include "RooPlot.h"
#include "RooCBShape.h"
#include "RooCruijff.h"
#include "RooGenericPdf.h"
#include "RooChi2Var.h"
#include "RooMinuit.h"


#include "interface/Configurator.h"
#include "plotOverallPerformancesConfigurator.h"
#include "interface/TopologicalSelectionHelper.h"

plotOverallPerformances_Config_t readConfiguration(std::string configName);

int main( int argc, char* argv[] ) {

  DrawTools::setStyle();
  
  
  std::string tag = "V00";
  std::string config = "SiPM2015Config";
  
  if( argc>1 ) {
    std::string tag_str(argv[1]);
    tag = tag_str;
    if( argc>2 ) {
      std::string config_str(argv[2]);
      config=config_str;
    }
  } else {
    std::cout << "Usage:" << std::endl;
    std::cout << "./timeWalkFitter ([tag]) ([config])" << std::endl;
    exit(12345);
  }

  std::cout<<config<<std::endl;
  
  theConfiguration_=readConfiguration(config);

  float shift=10;
  bool isNino=false;

  if(theConfiguration_.setup=="Nino"){
    shift=11.5;
    isNino=true;
  }

   std::vector<int> runs;
   std::vector<int> energies;

  TString haddFile = "analysisTrees_"+tag+"/"+theConfiguration_.setup+"_runs.root";

  if(theConfiguration_.addTagFileName){
    haddFile = "analysisTrees_"+tag+"/"+theConfiguration_.setup+"_"+theConfiguration_.tagFileName+"_runs.root";
  }

  TString haddString = "hadd analysisTrees_"+tag+"/"+theConfiguration_.setup+"_runs.root ";

  if(theConfiguration_.addTagFileName){
    haddString = "hadd analysisTrees_"+tag+"/"+theConfiguration_.setup+"_"+theConfiguration_.tagFileName+"_runs.root ";
  }

  for (int i=0;i<theConfiguration_.runs.size();++i){
    TString fileNameString = "analysisTrees_"+tag+"/Reco_";
    runs.push_back(theConfiguration_.runs[i]);
    energies.push_back(theConfiguration_.energies[i]);
    fileNameString+=runs[i];
    haddString+=fileNameString;
    haddString+= ".root ";
  }

  system(haddString.Data());

  TFile* file = TFile::Open(haddFile.Data());

  shift=0;//FIXME
  TProfile* deltaTvsAmpl=new TProfile("timeWalk","timeWalk",50,0,100,-10.,0.);
  TH2F* deltaTvsAmplHisto=new TH2F("timeWalkHisto","timeWalkHisto",50,0,100,100,-10.,0.);

  TH2F* deltaTvsAmplHisto_channel=new TH2F("timeWalkHisto_channel","timeWalkHisto_channel",50,0,100,100,-10.,0.);
  TH2F* deltaTvsAmplHisto_fibre=new TH2F("timeWalkHisto_fibre","timeWalkHisto_fibre",50,0,100,100,-10.,0.);


  TProfile* deltaTvsAmpl_channel=new TProfile("timeWalk_channel","timeWalk_channel",50,0,100,-10.,0.);
  TProfile* deltaTvsAmpl_fibre=new TProfile("timeWalk_fibre","timeWalk_fibre",50,0,100,-10.,0.);


  std::vector<TH1F*> deltaTBinAmpl_channel;
  std::vector<TH1F*> deltaTBinAmpl_fibre;

  float amplCut=20;
  int shiftBin=0;



  for (int i=1;i<deltaTvsAmplHisto->GetNbinsX();++i){//i=0 is underflow
    TString icut;
    icut.Form("%d",i); 
    if(deltaTvsAmplHisto->GetXaxis()->GetBinLowEdge(i)>amplCut){
      deltaTBinAmpl_channel.push_back(new TH1F("deltaTBinAmpl_channel"+icut,"deltaTBinAmpl_channel"+icut,200,-10.,0.));
      deltaTBinAmpl_fibre.push_back(new TH1F("deltaTBinAmpl_fibre"+icut,"deltaTBinAmpl_fibre"+icut,200,-10.,0.));
    }else{
      shiftBin++;
    }
  }




  

  deltaTvsAmplHisto->GetXaxis()->SetTitle("Signal Amplitude [ADC Channel]");
  deltaTvsAmplHisto->GetYaxis()->SetTitle("time_{fibre}-time_{mcp} [ADC Channel]");

  deltaTvsAmplHisto_channel->GetXaxis()->SetTitle("Signal Amplitude [ADC Channel]");
  deltaTvsAmplHisto_channel->GetYaxis()->SetTitle("time_{fibre}-time_{mcp} [ADC Channel]");

  deltaTvsAmplHisto_fibre->GetXaxis()->SetTitle("Signal Amplitude [ADC Channel]");
  deltaTvsAmplHisto_fibre->GetYaxis()->SetTitle("time_{fibre}-time_{mcp} [ADC Channel]");


  TTree* recoTree=(TTree*)file->Get("recoTree");
  RecoTree t(recoTree);
  Long64_t nentries = t.fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = t.LoadTree(jentry);
      if (ientry < 0) break;
      nb = t.fChain->GetEntry(jentry);   nbytes += nb;
      if( t.mcp_time_frac50<0 ||t.mcp_max_amplitude<200)continue;
      // if(isNino && t.nino_maxAmpl<32)continue;

      float deltaTNoCorr=t.cef3_time_at_frac50->at(1)-t.mcp_time_frac50;
      if(theConfiguration_.setup=="Nino")deltaTNoCorr=t.nino_LEtime-t.mcp_time_frac50;

      bool notFilled=true;

      if(t.nino_triggered && t.cef3_maxAmpl->at(1)>700){
	deltaTvsAmpl->Fill(t.nino_maxAmpl,deltaTNoCorr+shift);
	deltaTvsAmplHisto->Fill(t.nino_maxAmpl,deltaTNoCorr+shift);
	
	for(int i=1;i<deltaTBinAmpl_channel.size();++i) {//i=0 is underflow
	   
	      if(TopologicalSelectionHelper::passesChannelTopologicalSelection(t,isNino)){
		if(notFilled){
		  deltaTvsAmplHisto_channel->Fill(t.nino_maxAmpl,deltaTNoCorr+shift);
		  notFilled=false;
		}
		if(t.nino_maxAmpl<amplCut)break;
		if(t.nino_maxAmpl>deltaTvsAmplHisto->GetXaxis()->GetBinLowEdge(i) && t.nino_maxAmpl<deltaTvsAmplHisto->GetXaxis()->GetBinUpEdge(i)) {
		  deltaTBinAmpl_channel[i-shiftBin]->Fill(deltaTNoCorr);
		  break;
		}
	      }else if(TopologicalSelectionHelper::passesFibreTopologicalSelection(t,isNino)){
		if(notFilled){
		  deltaTvsAmplHisto_fibre->Fill(t.nino_maxAmpl,deltaTNoCorr+shift);
		  notFilled=false;
		}
		if(t.nino_maxAmpl<amplCut)break;
		if(t.nino_maxAmpl>deltaTvsAmplHisto->GetXaxis()->GetBinLowEdge(i) && t.nino_maxAmpl<deltaTvsAmplHisto->GetXaxis()->GetBinUpEdge(i)) {
		  deltaTBinAmpl_fibre[i-shiftBin]->Fill(deltaTNoCorr);
		  break;
		}
	
	      }


	   
	  }
      }
	
      }

   // this is the dir in which the plots will be saved:
   std::string constDirName = "timeWalkFiles_";
   constDirName+=theConfiguration_.setup;
   if(theConfiguration_.addTagFileName){
     constDirName+="_";
     constDirName+=theConfiguration_.tagFileName;
   }
   system(Form("mkdir -p %s", constDirName.c_str()));
   TString dir(constDirName);


   std::string outFileName;
   outFileName = dir+"/timeWalkFile_"+tag+".root";
   TFile* outFile = TFile::Open(outFileName.c_str(),"recreate");



   TVectorD mean_channel(deltaTBinAmpl_channel.size());
   TVectorD meanErr_channel(deltaTBinAmpl_channel.size());
   TVectorD sigma_channel(deltaTBinAmpl_channel.size());
   TVectorD sigmaErr_channel(deltaTBinAmpl_channel.size());
   TVectorD ampl_channel(deltaTBinAmpl_channel.size());
   TVectorD amplErr_channel(deltaTBinAmpl_channel.size());

   TVectorD mean_fibre(deltaTBinAmpl_fibre.size());
   TVectorD meanErr_fibre(deltaTBinAmpl_fibre.size());
   TVectorD sigma_fibre(deltaTBinAmpl_fibre.size());
   TVectorD sigmaErr_fibre(deltaTBinAmpl_fibre.size());
   TVectorD ampl_fibre(deltaTBinAmpl_fibre.size());
   TVectorD amplErr_fibre(deltaTBinAmpl_fibre.size());

   int dimCounter_channel=0;
   int dimCounter_fibre=0;

   for(int i=0;i<deltaTBinAmpl_channel.size();++i) {
     


     TString icut;
     icut.Form("%d",i); 
     //channel

     if(deltaTBinAmpl_channel[i]->GetEntries()>25){     
       dimCounter_channel++;
       if(i==0){
	 ampl_channel[i]=amplCut+deltaTvsAmplHisto->GetXaxis()->GetBinWidth(0)/2.;
       }else{
	 ampl_channel[i]=ampl_channel[i-1]+deltaTvsAmplHisto->GetXaxis()->GetBinWidth(0);
       }

       std::cout<<i<<" "<<ampl_channel[i]<<" kkkkkkkkkkkkkkkkkkkkkkkkk"<<std::endl;
       
       amplErr_channel[i]=0;
       
       TCanvas c1;
       c1.cd();
       
       TH1F* histo;
       histo=deltaTBinAmpl_channel[i];
       double peakpos = histo->GetBinCenter(histo->GetMaximumBin());
       double sigma = histo->GetRMS();

       double fitmin;
       double fitmax;
  
  
       fitmin = peakpos-4*sigma;
       fitmax = peakpos+4*sigma;
  
       RooRealVar x("x","deltaT", fitmin, fitmax);
       RooDataHist data("data","dataset with x",x,RooFit::Import(*histo) );

       RooRealVar meanr("meanr","Mean",peakpos+2*sigma,peakpos-3*sigma, peakpos+3*sigma);
       RooRealVar widthL("widthL","#sigmaL",sigma , 0, 5*sigma);
       RooRealVar widthR("widthR","#sigmaR",sigma , 0, 5*sigma);
       RooRealVar alphaL("alphaL","#alpha",5.08615e-02 , 0., 1.);
       RooRealVar alphaR("alphaR","#alpha",5.08615e-02, 0., 1.);
       int ndf;
     
       RooPlot* frame;

       RooCruijff fit_fct("fit_fct","fit_fct",x,meanr,widthL,widthR,alphaL,alphaR); ndf = 5;
       fit_fct.fitTo(data);
    
       std::string ytitle = Form("Events");
       frame = x.frame("Title");
       frame->SetXTitle("time_{fibre}-time_{mcp} [ns]");
       frame->SetYTitle(ytitle.c_str());
     
       data.plotOn(frame);  //this will show histogram data points on canvas 
       fit_fct.plotOn(frame);//this will show fit overlay on canvas  

       mean_channel[i]=meanr.getVal();
       meanErr_channel[i]=meanr.getError();
    
       double rms,rmsErr;
       rms = (widthL.getVal()+widthR.getVal())/2;
       rmsErr = 0.5*sqrt(widthL.getError()*widthL.getError()+widthR.getError()*widthR.getError());

       sigma_channel[i]=rms*1.e3; 	
       sigmaErr_channel[i]=rmsErr*1.e3; 


       frame->Draw();
       TLegend* lego = new TLegend(0.57, 0.8, 0.89, 0.9);;
       lego->SetTextSize(0.036);
       lego->SetTextAlign(32); // align right

       TPaveText* pave = DrawTools::getLabelTop_expOnXaxis("Electron Beam");
       pave->Draw("same");


 
       c1.SaveAs(dir+"/deltaTBinAmpl_channel"+icut+".png");
       c1.SaveAs(dir+"/deltaTBinAmpl_channel"+icut+".pdf");
     }

     //fibre
     if(deltaTBinAmpl_fibre[i]->GetEntries()>25){     
       dimCounter_fibre++;
       if(i==0){
	 ampl_fibre[i]=amplCut+deltaTvsAmplHisto->GetXaxis()->GetBinWidth(0)/2.;
       }else{
	 ampl_fibre[i]=ampl_fibre[i-1]+deltaTvsAmplHisto->GetXaxis()->GetBinWidth(0);
       }

       amplErr_fibre[i]=0;

       TCanvas c1;
       c1.cd();
       
       TH1F* histo;
       histo=deltaTBinAmpl_fibre[i];
       double peakpos = histo->GetBinCenter(histo->GetMaximumBin());
       double sigma = histo->GetRMS();

       double fitmin;
       double fitmax;
  
  
       fitmin = peakpos-4*sigma;
       fitmax = peakpos+4*sigma;
  
       RooRealVar x("x","deltaT", fitmin, fitmax);
       RooDataHist data("data","dataset with x",x,RooFit::Import(*histo) );

       RooRealVar meanr("meanr","Mean",peakpos+2*sigma,peakpos-3*sigma, peakpos+3*sigma);
       RooRealVar widthL("widthL","#sigmaL",sigma , 0, 5*sigma);
       RooRealVar widthR("widthR","#sigmaR",sigma , 0, 5*sigma);
       RooRealVar alphaL("alphaL","#alpha",5.08615e-02 , 0., 1.);
       RooRealVar alphaR("alphaR","#alpha",5.08615e-02, 0., 1.);
       int ndf;
     
       RooPlot* frame;

       RooCruijff fit_fct("fit_fct","fit_fct",x,meanr,widthL,widthR,alphaL,alphaR); ndf = 5;
       fit_fct.fitTo(data);
    
       std::string ytitle = Form("Events");
       frame = x.frame("Title");
       frame->SetXTitle("time_{fibre}-time_{mcp} [ns]");
       frame->SetYTitle(ytitle.c_str());
     
       data.plotOn(frame);  //this will show histogram data points on canvas 
       fit_fct.plotOn(frame);//this will show fit overlay on canvas  

       mean_fibre[i]=meanr.getVal();
       meanErr_fibre[i]=meanr.getError();
    
       double rms,rmsErr;
       rms = (widthL.getVal()+widthR.getVal())/2;
       rmsErr = 0.5*sqrt(widthL.getError()*widthL.getError()+widthR.getError()*widthR.getError());

       sigma_fibre[i]=rms*1.e3; 	
       sigmaErr_fibre[i]=rmsErr*1.e3; 


       frame->Draw();
       TLegend* lego = new TLegend(0.57, 0.8, 0.89, 0.9);;
       lego->SetTextSize(0.036);
       lego->SetTextAlign(32); // align right

       TPaveText* pave = DrawTools::getLabelTop_expOnXaxis("Electron Beam");
       pave->Draw("same");



       c1.SaveAs(dir+"/deltaTBinAmpl_fibre"+icut+".png");
       c1.SaveAs(dir+"/deltaTBinAmpl_fibre"+icut+".pdf");
     }
     
   }
   

   ampl_channel.ResizeTo(dimCounter_channel);
   mean_channel.ResizeTo(dimCounter_channel);
   amplErr_channel.ResizeTo(dimCounter_channel);
   meanErr_channel.ResizeTo(dimCounter_channel);
   sigma_channel.ResizeTo(dimCounter_channel);
   sigmaErr_channel.ResizeTo(dimCounter_channel);

   ampl_fibre.ResizeTo(dimCounter_fibre);
   mean_fibre.ResizeTo(dimCounter_fibre);
   amplErr_fibre.ResizeTo(dimCounter_fibre);
   meanErr_fibre.ResizeTo(dimCounter_fibre);
   sigma_fibre.ResizeTo(dimCounter_fibre);
   sigmaErr_fibre.ResizeTo(dimCounter_fibre);


   TGraphErrors* graph_fit_channel=new TGraphErrors(ampl_channel,mean_channel,amplErr_channel,meanErr_channel);
   graph_fit_channel->SetMarkerColor(kGreen+2);
   //   graph_fit_channel->SetMarkerSize(2);
   graph_fit_channel->SetMarkerStyle(34);
   graph_fit_channel->SetLineColor(kGreen+2);


   TGraphErrors* graph_fit_fibre=new TGraphErrors(ampl_fibre,mean_fibre,amplErr_fibre,meanErr_fibre);
   graph_fit_fibre->SetMarkerColor(kViolet);
   graph_fit_fibre->SetLineColor(kViolet);
   graph_fit_fibre->SetMarkerStyle(34);

   TF1 f_channel("fit_log_channel","[0]*log([1]*x)",0.,100.);
   //   TF1 f("fit_log","[0]+[1]/x",14,100);
   f_channel.SetParameters(-1.,20.);
   f_channel.SetLineColor(kGreen+2);
   //   deltaTvsAmpl->Fit("fit_log","R","",14,70);
   //   deltaTvsAmpl->Fit("fit_log","R","",20,100);
   graph_fit_channel->Fit("fit_log_channel","R","",20,60);


   TF1 f_fibre("fit_log_fibre","[0]*log([1]*x)",0.,100.);
   //   TF1 f("fit_log","[0]+[1]/x",14,100);
      f_fibre.SetParameters(-4.98330e-01 ,4.06173e+04);
   //   f_fibre.FixParameter(0,-4.98330e-01);
   //   f_fibre.FixParameter(1, 4.06173e+04);
   f_fibre.SetLineColor(kViolet);  
   //   deltaTvsAmpl->Fit("fit_log","R","",14,70);
   //   deltaTvsAmpl->Fit("fit_log","R","",20,100);
   graph_fit_fibre->Fit("fit_log_fibre","R","",20,50);


   TGraphErrors*   resVsAmplitude_fibre = new TGraphErrors(0);
   TGraphErrors*   resVsAmplitude_channel = new TGraphErrors(0);
   for(int i=0;i<sigma_channel.GetNoElements();i++){
     if(sigma_channel[i]!=0){ 
       resVsAmplitude_channel->SetPoint(i,amplCut+(i+1)*deltaTvsAmplHisto->GetXaxis()->GetBinWidth(0)-deltaTvsAmplHisto->GetXaxis()->GetBinWidth(0)/2,sigma_channel[i]); 
       resVsAmplitude_channel->SetPointError(i,0,sigmaErr_channel[i]);
     }
   }
   resVsAmplitude_channel->SetName("resVsAmplitude_channel");
   //fibre
   for(int i=0;i<sigma_fibre.GetNoElements();i++){
     if(sigma_fibre[i]!=0){ 
       resVsAmplitude_fibre->SetPoint(i,amplCut+(i+1)*deltaTvsAmplHisto->GetXaxis()->GetBinWidth(0)-deltaTvsAmplHisto->GetXaxis()->GetBinWidth(0)/2,sigma_fibre[i]); 
       resVsAmplitude_fibre->SetPointError(i,0,sigmaErr_fibre[i]);
     }
   }
   resVsAmplitude_fibre->SetName("resVsAmplitude_fibre");


   TCanvas* c_res = new TCanvas( "c_res", "", 600, 600 );
   float yup=1.1*(resVsAmplitude_channel->GetY()[0]+resVsAmplitude_channel->GetEY()[0]);
   float xup=(resVsAmplitude_channel->GetX())[resVsAmplitude_channel->GetN()-1]+5; 
   TH2D* h2_axes = new TH2D( "axes", "", 100,resVsAmplitude_channel->GetX()[0]-5 ,xup , 110, 0., yup);
   
   
   //channel   
   resVsAmplitude_channel->SetName("resVsAmplitude_channel");
   resVsAmplitude_channel->SetMarkerStyle(20);
   resVsAmplitude_channel->SetMarkerSize(1.6);
   resVsAmplitude_channel->SetMarkerColor(kBlue);
   
   
   resVsAmplitude_channel->Write();
   h2_axes->SetYTitle("#sigma_{t} [ps]");
   h2_axes->GetXaxis()->SetTitle("Amplitude [ADC]");
   h2_axes->Draw(""); 
   
   TF1* f= new TF1("fun","sqrt([1]*[1]/(x*x)+[0]*[0])",(resVsAmplitude_channel->GetX())[0]-2,(resVsAmplitude_channel->GetX())[resVsAmplitude_channel->GetN()-1] +5);
   f->SetParLimits(0,50,300);
   f->SetParameter(0,100);
   f->SetParameter(1,1.12135e+04);
   
   //     f->SetParLimits(1,10,1.19257e+02);
   //   resVsAmplitude_channel->Fit("fun","","",(resVsAmplitude_channel->GetX())[0]-2,(resVsAmplitude_channel->GetX())[resVsAmplitude_channel->GetN()-4]);
   
   resVsAmplitude_channel->Draw("p same");
   
   TPaveText* pave = DrawTools::getLabelTop_expOnXaxis("Electron Beam");
   pave->Draw("same");
   
   
   
   TLegend* lego = new TLegend(0.47, 0.7, 0.8, 0.92);
   lego->SetTextSize(0.038);
   //    lego->AddEntry(  (TObject*)0 ,"f(x) = p0 + p1/x", "");
   lego->AddEntry(  (TObject*)0 ,Form("C = %.0f #pm %.0f ps", f->GetParameter(0), f->GetParError(0) ), "");
   //    lego->AddEntry(  (TObject*)0 ,Form("N = %.0f #pm %.0f", f->GetParameter(1), f->GetParError(1) ), "");
   lego->SetFillColor(0);
   //   lego->Draw("same");
   
   c_res->SaveAs(dir+"/timingResolutionVsAmplitudeSiPM_channel.png");
   c_res->SaveAs(dir+"/timingResolutionVsAmplitudeSiPM_channel.pdf");  

   //fibre
   //channel   
   c_res->Clear();
   yup=1.1*(resVsAmplitude_fibre->GetY()[0]+resVsAmplitude_fibre->GetEY()[0]);
   xup=(resVsAmplitude_fibre->GetX())[resVsAmplitude_fibre->GetN()-1]+5; 
   TH2D* h2_axes_2 = new TH2D( "axes", "", 100,resVsAmplitude_fibre->GetX()[0]-5 ,xup , 110, 0., yup);

   resVsAmplitude_fibre->SetName("resVsAmplitude_fibre");
   resVsAmplitude_fibre->SetMarkerStyle(20);
   resVsAmplitude_fibre->SetMarkerSize(1.6);
   resVsAmplitude_fibre->SetMarkerColor(kBlue);
   
   
   resVsAmplitude_fibre->Write();
   h2_axes_2->SetYTitle("#sigma_{t} [ps]");
   h2_axes_2->GetXaxis()->SetTitle("Amplitude [ADC]");
   h2_axes_2->Draw(""); 
   
   
   resVsAmplitude_fibre->Draw("p same");
   
   pave->Draw("same");
   
   
   
   c_res->SaveAs(dir+"/timingResolutionVsAmplitudeSiPM_fibre.png");
   c_res->SaveAs(dir+"/timingResolutionVsAmplitudeSiPM_fibre.pdf");  
     
     

   gStyle->SetPadRightMargin(0.17);//for the palette

   TCanvas c1;
   c1.cd();
   deltaTvsAmplHisto->Draw("colz");
   //   deltaTvsAmpl->Draw("same");
   graph_fit_channel->Draw("sameep");
   graph_fit_fibre->Draw("sameep");
   c1.SaveAs(dir+"/deltaTvsAmpl.png");
   c1.SaveAs(dir+"/deltaTvsAmpl.pdf");


   c1.Clear();
   c1.cd();
   deltaTvsAmplHisto_channel->Draw("colz");
   //   deltaTvsAmpl->Draw("same");
   graph_fit_channel->Draw("sameep");
   c1.SaveAs(dir+"/deltaTvsAmpl_channel.png");
   c1.SaveAs(dir+"/deltaTvsAmpl_channel.pdf");

   c1.Clear();
   c1.cd();
   deltaTvsAmplHisto_fibre->Draw("colz");
   //   deltaTvsAmpl->Draw("same");
   graph_fit_fibre->Draw("sameep");
   c1.SaveAs(dir+"/deltaTvsAmpl_fibre.png");
   c1.SaveAs(dir+"/deltaTvsAmpl_fibre.pdf");


   f_channel.Write();
   f_fibre.Write();

   graph_fit_channel->Write();
   graph_fit_fibre->Write();
   deltaTvsAmpl->Write();
   deltaTvsAmplHisto->Write();




   outFile->Write();
   outFile->Close();




}


plotOverallPerformances_Config_t readConfiguration(std::string configName){

  Configurator* configurator_ = new Configurator();
  std::string fileName = Form ("./config_timing/plots_performances_%s.xml",configName.c_str());
  configurator_->xmlFileName=fileName.c_str();
  configurator_->Init();

  plotOverallPerformances_Config_t conf;
 
  conf.setup=Configurable::getElementContent (*configurator_, "setup",configurator_->root_element) ;

  string  content = Configurable::getElementContent (*configurator_, "runs",configurator_->root_element) ;
  Configurator::GetVecInt(content,conf.runs);

   conf.nMaxAmplCuts= Configurator::GetInt(Configurable::getElementContent(*configurator_,"nMaxAmplCuts",configurator_->root_element));

   conf.stepAmplFibre= Configurator::GetInt(Configurable::getElementContent(*configurator_,"stepAmplFibre",configurator_->root_element));

   conf.stepAmplChannel= Configurator::GetInt(Configurable::getElementContent(*configurator_,"stepAmplChannel",configurator_->root_element));

   conf.rangeXLow= Configurator::GetDouble(Configurable::getElementContent(*configurator_,"rangeXLow",configurator_->root_element));

   conf.rangeXUp= Configurator::GetDouble(Configurable::getElementContent(*configurator_,"rangeXUp",configurator_->root_element));

   conf.channel2CutFibre= Configurator::GetDouble(Configurable::getElementContent(*configurator_,"channel2CutFibre",configurator_->root_element));
   conf.channel2CutChannel= Configurator::GetDouble(Configurable::getElementContent(*configurator_,"channel2CutChannel",configurator_->root_element));
   conf.channel1CutChannel= Configurator::GetDouble(Configurable::getElementContent(*configurator_,"channel1CutChannel",configurator_->root_element));

   conf.startCutFibre= Configurator::GetDouble(Configurable::getElementContent(*configurator_,"startCutFibre",configurator_->root_element));

   conf.startCutChannel= Configurator::GetDouble(Configurable::getElementContent(*configurator_,"startCutChannel",configurator_->root_element));

   conf.amplCut= Configurator::GetDouble(Configurable::getElementContent(*configurator_,"amplCut",configurator_->root_element));

   conf.addTagFileName= Configurator::GetInt(Configurable::getElementContent(*configurator_,"addTagFileName",configurator_->root_element));

   conf.tagFileName=Configurable::getElementContent (*configurator_, "tagFileName",configurator_->root_element) ;

  content = Configurable::getElementContent (*configurator_, "energies",configurator_->root_element) ;
  Configurator::GetVecFloat(content,conf.energies);


   return conf;

}
