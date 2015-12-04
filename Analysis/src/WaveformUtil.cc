#define WaveformUtil_cxx
#include "WaveformUtil.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TRandom3.h>

#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooGenericPdf.h"
#include "RooFFTConvPdf.h"
#include "RooPlot.h"
#include "RooTFnBinding.h"
#include "RooTFnPdfBinding.h"
#include  "RooNumConvPdf.h"
#include "TVectorD.h"
#include "TGaxis.h"
#include "DrawTools.h"

#include "interface/channelInfo.h"
#include "interface/HodoCluster.h"
#include "interface/TagHelper.h"
#include "interface/EnergyCalibration.h"
#include "interface/AlignmentOfficer.h"



void assignValues( std::vector<float> &target, std::vector<float> source, unsigned int startPos );
void assignValuesBool( std::vector<bool> &target, std::vector<bool> source, unsigned int startPos );
void doHodoReconstructionBool( std::vector<bool> values, int &nClusters, int *nFibres, float *pos, float fibreWidth, int clusterMaxFibres, float Cut );
std::vector<HodoCluster*> getHodoClustersBool( std::vector<bool> hodo, float fibreWidth, int nClusterMax, float Cut );
void doHodoReconstruction( std::vector<float> values, int &nClusters, int *nFibres, float *pos, float fibreWidth, int clusterMaxFibres, float Cut );
std::vector<HodoCluster*> getHodoClusters( std::vector<float> hodo, float fibreWidth, int nClusterMax, float Cut );
void copyArray( int n, float *source, float *target );
double funcCRRC(double *x, double *par);
bool pickEvents(int evtNumber);

void WaveformUtil::Loop(){

  DrawTools::setStyle();

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();
  if(nentries>50000)  nentries=5000;
  //  nentries=241;
  std::cout<<"nentries"<<nentries<<std::endl;

  float  timeMinusTimeAtMax[NDIGISAMPLES];
  float  mean[NFIBERS][NDIGISAMPLES];
  float  time[NDIGISAMPLES];

  float meanTimeAtMax[NFIBERS];
  TGraph* meanWaveGraphsForPlots[NFIBERS];  


  for (int i=0;i<NFIBERS;++i){
    meanTimeAtMax[i]=0;
    for (int j=0;j<NDIGISAMPLES;++j){
      mean[i][j]=0;
    }
  }

  TH1F* shiftSampleHisto=new TH1F("shiftSample","shiftSample",30,-15.5,14.5);

  TString runNumberString;
  int digiFreq;
  Long64_t nbytes = 0, nb = 0;

  bool isOctober2015Run=false;

  int nEffEntries=0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;


    // if (Cut(ientry) < 0) continue;      

    if(jentry==0){
      isOctober2015Run=  (runNumber > 3900. && runNumber<4200);
      runNumberString.Form("%d",runNumber);
      digiFreq=digi_frequency;
    }
    if(jentry%1000 == 0)std::cout<<"Processing entry:"<<jentry<<std::endl;




    //    if(passesHodoSelection()==false)continue;
      nEffEntries++;

    //    if(TMath::Abs(TDCreco->at(0))>2 || TMath::Abs(TDCreco->at(1))>2)continue;    //FIXME

    //    float timeOfTheEvent=digi_time_at_frac50_bare_noise_sub->at(8);//synchronizing time of events with time of trigger
    float timeOfTheEvent=digi_time_at_1000_bare_noise_sub->at(8);//synchronizing time of events with time of trigger
    //    float shiftTime=190.3-timeOfTheEvent;//mean fitted on trigger run 2778
    float shiftTime=0;
    if (digiFreq==1){
      shiftTime=190.2-timeOfTheEvent;//mean fitted on trigger run 2778
      if(isOctober2015Run)shiftTime=139.2-timeOfTheEvent;
      if(runNumber>4490  && runNumber<4550)shiftTime=156.6-timeOfTheEvent;//last day of beam test configuration
    }
    if (digiFreq==0)shiftTime=161.9-timeOfTheEvent;//mean fitted on trigger run 3076
    //    float shiftTime=138.1-timeOfTheEvent;//mean fitted on trigger run 329 for 2014 data
    int shiftSample=round(shiftTime/(1e9*timeSampleUnit(digiFreq)));
    shiftSample=-shiftSample;
    shiftSampleHisto->Fill(shiftSample);     

    //    if(!pickEvents(evtNumber))continue;

    for (int i=0;i<(CEF3_CHANNELS+1*isOctober2015Run);++i){
      int iChannel=i;
      if(isOctober2015Run){
	if (i==2)continue;//channel2 was broken in October test
	if (i>2)iChannel--;//channel2 was broken in October test
      }
      if(digi_time_at_max_bare_noise_sub->at(i)>0 && digi_time_at_max_bare_noise_sub->at(i)<400)meanTimeAtMax[iChannel]+=digi_time_at_max_bare_noise_sub->at(i); 
    }
    for (int i=0;i<1024*(CEF3_CHANNELS+1*isOctober2015Run);++i){
      int iChannel=digi_value_ch->at(i);
      if(isOctober2015Run){
	if(digi_value_ch->at(i)==2) continue;
	if(digi_value_ch->at(i)>2)iChannel--;//channel2 was broken in October test
      }

      if(digi_max_amplitude->at(digi_value_ch->at(i))>10000 || digi_max_amplitude->at(digi_value_ch->at(i))<0)continue;
      int iSample=i;
      if(i+shiftSample>1023*digi_value_ch->at(i) && i+shiftSample<(1023+(1024*digi_value_ch->at(i)))){
	iSample=i+shiftSample;
      }
      //      if(digi_value_ch->at(i)==1 )      std::cout<<"i:"<<i<<" isample:"<<iSample<<"digivalue:"<<digi_value_bare_noise_sub->at(i)<<" digivalue new:"<<digi_value_bare_noise_sub->at(iSample)<<" channel:"<<digi_value_ch->at(iSample)<<std::endl;

      mean[iChannel][i-1024*digi_value_ch->at(i)]+=(float)(digi_value_bare_noise_sub->at(iSample));
      if(i<1024)time[i]=digi_value_time->at(i);

   }
  }

  for (int i=0;i<NFIBERS;++i){
    for (int j=0;j<NDIGISAMPLES;++j){
      mean[i][j]/=nEffEntries;
    }
  }
  //  exit(0);

  TFile* outFile = TFile::Open("outWaveFormUtil_"+runNumberString+".root","recreate");
  int lowRange[4],highRange[4];
  TVectorD integrals(4);
  shiftSampleHisto->Write();
  lowRange[0]=190;
  lowRange[1]=170;
  lowRange[2]=130;
  lowRange[3]=170;
  highRange[0]=500;
  highRange[1]=720;
  highRange[2]=600;
  highRange[3]=680;

  for (int i=0;i<NFIBERS;++i){
    meanWaveGraphs[i]=new TGraph(1024, time, mean[i]);
    TString fiber;
    fiber.Form("%d",i);
    meanWaveGraphs[i]->SetName("waveform_"+fiber);

    meanWaveHistos[i]=new TH1F("waveform_histo_"+fiber,"",1024,0,time[1023]);
    //    meanWaveHistosForPlots[i]=new TH1F("waveform_histo_forplots_"+fiber,"",1024,-300*1e9*time[1],1e9*time[1023-300]);
    meanTimeAtMax[i]/=nEffEntries;
    //std::cout<<meanTimeAtMax[i]<< " unit:"<<timeSampleUnit(digiFreq)<<std::endl;
    if(runNumberString=="4054" || runNumberString=="4053")meanTimeAtMax[1]=24;//FIXME
    meanWaveHistosForPlots[i]=new TH1F("waveform_histo_forplots_"+fiber,"",1024,-meanTimeAtMax[i],1e9*time[1023-(int)(meanTimeAtMax[i]*1.e-9/timeSampleUnit(digiFreq))]);

    for (int j=0;j<NDIGISAMPLES;++j){  
      timeMinusTimeAtMax[j]=(1e9*(time[j])-meanTimeAtMax[i]);
      for(int k=0;k<mean[i][j];k++){
	meanWaveHistos[i]->Fill(time[j]);
	meanWaveHistosForPlots[i]->Fill(1e9*(time[j])-meanTimeAtMax[i]);
      }
      mean[i][j]*=0.25; //be careful!this is to do meanWaveGraphsForPlots in milliVolt, change it if you need to reuse it
    }

    meanWaveGraphsForPlots[i]=new TGraph(1024, timeMinusTimeAtMax, mean[i]);
    meanWaveGraphsForPlots[i]->SetName("waveform_plots_"+fiber);

    meanWaveHistos[i]->Scale(nEffEntries);
    meanWaveHistos[i]->Sumw2();
    meanWaveHistos[i]->Scale(1./nEffEntries);//this trick is just to get correct errors
    meanWaveHistosForPlots[i]->Scale(nEffEntries);
    meanWaveHistosForPlots[i]->Sumw2();
    meanWaveHistosForPlots[i]->Scale(1./nEffEntries);//this trick is just to get correct errors
    meanWaveHistosForPlots[i]->Scale(0.25);//we need millivolts

    meanWaveHistos[i]->Write();
    meanWaveHistosForPlots[i]->Write();
    meanWaveGraphs[i]->Write();
    meanWaveGraphsForPlots[i]->Write();

    TCanvas can;
    gStyle->SetPadRightMargin(0.15);
    //    meanWaveGraphs[i]->GetYaxis()->SetRangeUser(0.,meanWaveGraphs[i]->GetMaximum()*1.10);
    meanWaveGraphs[i]->GetXaxis()->SetTitle("t [s]");

    meanWaveGraphs[i]->Draw("APL");
    can.SaveAs("plots/meanWave_"+fiber+"_"+runNumberString+".png");
    can.SaveAs("plots/meanWave_"+fiber+"_"+runNumberString+".pdf");

//    TCanvas fitcanvas;
//    TF1* f1 = new TF1( "func", funcCRRC, lowRange[i]*timeSampleUnit(digiFreq), highRange[i]*timeSampleUnit(digiFreq)*3, 4 );
//    //    f1->SetParameter( 0, 1.25475e+07);
//    //    f1->SetParameter( 1, 6.71014e+08);
//    meanWaveHistos[i]->Fit(f1,"RN+");
//    meanWaveHistos[i]->Draw();
//    std::cout<<lowRange[i]*timeSampleUnit(digiFreq)<<" "<< highRange[i]*timeSampleUnit(digiFreq)<<std::endl;
//    f1->SetLineColor(kRed);
//    f1->Draw("same");
//    fitcanvas.SaveAs("plots/fitCRRC_"+fiber+".png");
    //    exit(0);

    RooRealVar t("t","t", time[0], time[1023]);
    RooDataHist data("data","dataset with t",t,RooFit::Import(*meanWaveHistos[i]) );
    //    RooDataSet data("data","dataset with t",t);



    // Construct landau(t,ml,sl) ;
    RooRealVar ml("ml","mean landau",0.05e-6,0,0.12e-6) ;
    RooRealVar sl("sl","sigma landau",0.05e-6,0.,0.12e-6) ;
    RooLandau landau("lx","lx",t,ml,sl) ;
  
    // Construct gauss(t,mg,sg)
    RooRealVar mg("mg","mg",1.06706e-07,0.,0.12e-6) ;
    RooRealVar sg("sg","sg",3.91767e-08,0.,0.12e-6) ;
    RooGaussian gauss("gauss","gauss",t,mg,sg) ;

    // Construct landau (x) gauss
    RooFFTConvPdf lxg("lxg","landau (X) gauss",t,landau,gauss) ;

    // Fit gxlx to data
    //    lxg.fitTo(data,RooFit::Range(time[lowRange[i]],time[highRange[i]])) ;

    // Plot data, landau pdf, landau (X) gauss pdf
    RooPlot* frame = t.frame(RooFit::Title("fiber_"+fiber)) ;
    t.setBins(10);
    data.plotOn(frame,RooFit::DataError(RooAbsData::SumW2)) ;

    //        lxg.plotOn(frame) ;



    gStyle->SetPadRightMargin(0.05);
    TGaxis::SetMaxDigits(3);

    TCanvas c1("c_"+fiber);
    frame->Draw();
    c1.SaveAs("plots/c_"+fiber+"_"+runNumberString+".png");
    c1.SaveAs("plots/c_"+fiber+"_"+runNumberString+".pdf");
    c1.Write();

    //try crrc function convoluted with decay of cef3 and wls
    /*
    //    TF1 *fcefAndWls = new TF1("f1","1.e10*(exp(-1.e9*x/33-exp(-1.e9*x/7)))",time[0],time[1023]);
    TF1 *fcefAndWls = new TF1("f1","1.e10*(exp(-1.e9*x/33-exp(-1.e9*x/7)))",time[0],time[1023]);


    // Create pdf from tf1
    //RooAbsReal* rfa1 = RooFit::bindFunction(fcefAndWls,t) ;
    RooAbsPdf* rfa1 = RooFit::bindPdf(fcefAndWls,t) ;

    RooRealVar tau("tau","tau",1.22451e+7,1.e+5,pow(10,30)) ;
    RooRealVar timeshift("timeshift","timeshift",0.03e-6,0,0.05e-6) ;

      RooGenericPdf crrc("crrc","(@0-0.05e-6)*@1*exp(-(@0-0.05e-6)*@1)",RooArgSet(t,tau)) ;
    //    RooGenericPdf crrc("crrc","(@0-@2)*@1*exp(-(@0-@2)*@1)",RooArgSet(t,tau,timeshift)) ;


    // Construct landau(t,ml,sl) ;
    //    RooRealVar ml2("ml","mean landau",0.,0,0.12e-6) ;
    //    RooRealVar sl2("sl","sigma landau",33e-9,0.,0.12e-6) ;
    //    RooLandau landau2("lx","lx",t,ml2,sl2) ;
    //    ml2.setConstant(kTRUE);
    //    sl2.setConstant(kTRUE);

    RooNumConvPdf excrrc("excrrc","exp (X) crrc",t,*rfa1,crrc) ;
    RooRealVar center("center","center",time[600],0,0.4e-6);
    RooRealVar width("width","width",time[599],0,0.4e-6);
    std::cout<<"####################center:"<<time[600]<<" width:"<<time[599]<<std::endl; 
    excrrc.setConvolutionWindow(center,width,1);

    excrrc.fitTo(data,RooFit::Range(0,time[highRange[i]]*2)) ;
    //    crrc.fitTo(data,RooFit::Range(time[0]+0.05e-6,time[1023])) ;
    //    exit(9);
    RooPlot* frame2 = t.frame(RooFit::Title("fiber_"+fiber)) ;
    data.plotOn(frame2,RooFit::DataError(RooAbsData::SumW2)) ;

    excrrc.plotOn(frame2) ;
    crrc.plotOn(frame2,RooFit::LineColor(kRed)) ;
    rfa1->plotOn(frame2,RooFit::LineColor(kViolet)) ;

    TCanvas c2("c_crrc_"+fiber);
    frame2->Draw();
    c2.Write();

    //END OF CRRC FUNCTION
    */

    TH1* histopdf=lxg.createHistogram("t",1024);
    histopdf->SetName("fitted_histo_"+fiber);
    float finalFastSample=230;
    if(digiFreq==1)finalFastSample=172;
    integrals[i]=histopdf->Integral(4,finalFastSample)/histopdf->Integral(4,900);
    histopdf->Write();

    c1.Clear();

    //    meanWaveHistos[i]->Scale(1./meanWaveHistos[i]->Integral());
    std::string energy(Form("%.0f", BeamEnergy));
    TPaveText* pave = DrawTools::getLabelTop(energy+" GeV Electron Beam");


    std::cout<<meanTimeAtMax[i]<<std::endl;
    //      meanWaveGraphsForPlots[i]->GetXaxis()->SetRangeUser(-meanTimeAtMax[i],1e9*time[1023]-meanTimeAtMax[i]));
    TH2F* axis=new TH2F("axis","axis",100,-80,290,100,0,meanWaveHistosForPlots[i]->GetMaximum()*1.10);
    //    TH2F* axis=new TH2F("axis","axis",100,-10,190,100,0,28);
    axis->GetXaxis()->SetTitle("Time [ns]");
    axis->GetYaxis()->SetTitle("Signal Amplitude [mV]");
    axis->Draw();
    meanWaveGraphsForPlots[i]->Draw("plsame");
    pave->Draw("same");
    c1.Write("pulseShape"+fiber);
    c1.SaveAs("plots/pulseShape"+fiber+"_"+runNumberString+".png");
    c1.SaveAs("plots/pulseShape"+fiber+"_"+runNumberString+".pdf");
  }

  integrals.Write("integrals_fast");


  outFile->Write();
  outFile->Close();

}



float WaveformUtil::timeSampleUnit(int drs4Freq)
{
  if (drs4Freq == 0)
    return 0.2E-9;
  else if (drs4Freq == 1)
    return 0.4E-9;
  else if (drs4Freq == 2)
    return 1.E-9;
  return -999.;
}


bool WaveformUtil::passesHodoSelection(){
  int nClusters_hodoX1;
  int nFibres_hodoX1[HODOX1_CHANNELS]; 
  float pos_hodoX1[HODOX1_CHANNELS];   
  float pos_corr_hodoX1[HODOX1_CHANNELS];

  int nClusters_hodoY1;
  int nFibres_hodoY1[HODOY1_CHANNELS]; 
  float pos_hodoY1[HODOY1_CHANNELS];   
  float pos_corr_hodoY1[HODOX1_CHANNELS];

  int nClusters_hodoX2;
  int nFibres_hodoX2[HODOX2_CHANNELS]; 
  float pos_hodoX2[HODOX2_CHANNELS];   
  float pos_corr_hodoX2[HODOX2_CHANNELS];   

  int nClusters_hodoY2;
  int nFibres_hodoY2[HODOY2_CHANNELS]; 
  float pos_hodoY2[HODOY2_CHANNELS];   
  float pos_corr_hodoY2[HODOX2_CHANNELS];   
  
  std::vector<bool> hodoX1_values(HODOX1_CHANNELS, -1.);
  std::vector<bool> hodoY1_values(HODOY1_CHANNELS, -1.);
  assignValuesBool( hodoX1_values, *HODOX1, 0. );
  assignValuesBool( hodoY1_values, *HODOY1, 0. );
  
  
  std::vector<bool> hodoX2_values(HODOX2_CHANNELS, -1.);
  std::vector<bool> hodoY2_values(HODOY2_CHANNELS, -1.);
  assignValuesBool( hodoX2_values, *HODOX2, 0 );
  assignValuesBool( hodoY2_values, *HODOY2, 0 );
  
  
     // hodo cluster reconstruction
  int clusterMaxFibres = 4;
  doHodoReconstructionBool( hodoX1_values    , nClusters_hodoX1    , nFibres_hodoX1    , pos_hodoX1    , 0.5, clusterMaxFibres, 0. );
  doHodoReconstructionBool( hodoY1_values    , nClusters_hodoY1    , nFibres_hodoY1    , pos_hodoY1    , 0.5, clusterMaxFibres, 0. );
  doHodoReconstructionBool( hodoX2_values    , nClusters_hodoX2    , nFibres_hodoX2    , pos_hodoX2    , 0.5, clusterMaxFibres , 0.);
  doHodoReconstructionBool( hodoY2_values    , nClusters_hodoY2    , nFibres_hodoY2    , pos_hodoY2    , 0.5, clusterMaxFibres, 0. );
  
  copyArray( nClusters_hodoX1, pos_hodoX1, pos_corr_hodoX1 );
  copyArray( nClusters_hodoY1, pos_hodoY1, pos_corr_hodoY1 );
  copyArray( nClusters_hodoX2, pos_hodoX2, pos_corr_hodoX2 );
  copyArray( nClusters_hodoY2, pos_hodoY2, pos_corr_hodoY2 );
  

  std::string theBeamEnergy = Form("%.0f",100);
  TagHelper tagHelper("V10",theBeamEnergy);
  AlignmentOfficer alignOfficer(tagHelper.getAlignmentFileName());

  alignOfficer.fix("hodoX1", nClusters_hodoX1, pos_corr_hodoX1);
  alignOfficer.fix("hodoY1", nClusters_hodoY1, pos_corr_hodoY1);
  alignOfficer.fix("hodoX2", nClusters_hodoX2, pos_corr_hodoX2);
  alignOfficer.fix("hodoY2", nClusters_hodoY2, pos_corr_hodoY2);
     

  
  if(nClusters_hodoX1!=1 || nClusters_hodoX2!=1 || nClusters_hodoY1!=1 || nClusters_hodoY2!=1) {return false;
  }else{
  return true;
  }
}



void assignValues( std::vector<float> &target, std::vector<float> source, unsigned int startPos ) {

  for( unsigned i=0; i<target.size(); ++i ) 
    target[i] = source[startPos+i];

}

double funcCRRC(double *x, double *par)
{
//  double ctime = (x[0] - par[0]) / par[1];
//  double f = 0.;
//  if(ctime > 0.)
//    f = ctime * exp( 1. - ctime );

  double  f=  par[1]*x[0]*exp(-x[0]*par[0]);

  return f;
} 

bool pickEvents(int evtNumber){

  if(evtNumber == 10 || evtNumber == 17 || evtNumber == 18 || evtNumber == 26 || evtNumber == 33 || evtNumber == 36 || evtNumber == 43 || evtNumber == 52 || evtNumber == 59 || evtNumber == 66 || evtNumber == 68 || evtNumber == 70 || evtNumber == 92 || evtNumber == 95 || evtNumber == 97 || evtNumber == 99 || evtNumber == 101 || evtNumber == 109 || evtNumber == 115 || evtNumber == 117 || evtNumber == 122 || evtNumber == 123 || evtNumber == 136 || evtNumber == 138 || evtNumber == 147 || evtNumber == 148 || evtNumber == 154 || evtNumber == 155 || evtNumber == 159 || evtNumber == 162 || evtNumber == 165 || evtNumber == 168 || evtNumber == 169 || evtNumber == 175 || evtNumber == 176 || evtNumber == 179 || evtNumber == 185 || evtNumber == 187 || evtNumber == 188 || evtNumber == 189 || evtNumber == 195 || evtNumber == 197 || evtNumber == 199 || evtNumber == 202 || evtNumber == 203 || evtNumber == 217 || evtNumber == 220 || evtNumber == 230 || evtNumber == 231 || evtNumber == 235 || evtNumber == 245 || evtNumber == 249 || evtNumber == 256 || evtNumber == 265 || evtNumber == 267 || evtNumber == 268 || evtNumber == 270 || evtNumber == 272 || evtNumber == 277 || evtNumber == 278 || evtNumber == 284 || evtNumber == 286 || evtNumber == 289 || evtNumber == 293 || evtNumber == 295 || evtNumber == 296 || evtNumber == 300 || evtNumber == 304 || evtNumber == 306 || evtNumber == 307 || evtNumber == 310 || evtNumber == 311 || evtNumber == 312 || evtNumber == 316 || evtNumber == 328 || evtNumber == 332 || evtNumber == 333 || evtNumber == 338 || evtNumber == 341 || evtNumber == 344 || evtNumber == 346 || evtNumber == 347 || evtNumber == 353 || evtNumber == 354 || evtNumber == 363 || evtNumber == 370 || evtNumber == 374 || evtNumber == 375 || evtNumber == 381 || evtNumber == 384 || evtNumber == 389 || evtNumber == 390 || evtNumber == 393 || evtNumber == 394 || evtNumber == 398 || evtNumber == 399 || evtNumber == 404 || evtNumber == 408 || evtNumber == 409 || evtNumber == 410 || evtNumber == 413 || evtNumber == 415 || evtNumber == 417 || evtNumber == 418 || evtNumber == 419 || evtNumber == 422 || evtNumber == 424 || evtNumber == 428 || evtNumber == 430 || evtNumber == 432 || evtNumber == 438 || evtNumber == 439 || evtNumber == 442 || evtNumber == 445 || evtNumber == 447 || evtNumber == 455 || evtNumber == 456 || evtNumber == 457 || evtNumber == 464 || evtNumber == 465 || evtNumber == 467 || evtNumber == 472 || evtNumber == 476 || evtNumber == 477 || evtNumber == 481 || evtNumber == 482 || evtNumber == 490 || evtNumber == 491 || evtNumber == 493 || evtNumber == 500 || evtNumber == 504 || evtNumber == 510 || evtNumber == 514 || evtNumber == 518 || evtNumber == 525 || evtNumber == 530 || evtNumber == 535 || evtNumber == 536 || evtNumber == 540 || evtNumber == 546 || evtNumber == 547 || evtNumber == 551 || evtNumber == 554 || evtNumber == 555 || evtNumber == 557 || evtNumber == 562 || evtNumber == 568 || evtNumber == 569 || evtNumber == 573 || evtNumber == 581 || evtNumber == 591 || evtNumber == 595 || evtNumber == 601 || evtNumber == 602 || evtNumber == 606 || evtNumber == 608 || evtNumber == 609 || evtNumber == 613 || evtNumber == 616 || evtNumber == 622 || evtNumber == 633 || evtNumber == 635 || evtNumber == 636 || evtNumber == 638 || evtNumber == 641 || evtNumber == 642 || evtNumber == 644 || evtNumber == 651 || evtNumber == 653 || evtNumber == 664 || evtNumber == 665 || evtNumber == 678 || evtNumber == 681 || evtNumber == 686 || evtNumber == 687 || evtNumber == 690 || evtNumber == 691 || evtNumber == 693 || evtNumber == 698 || evtNumber == 700 || evtNumber == 702 || evtNumber == 703 || evtNumber == 708 || evtNumber == 721 || evtNumber == 730 || evtNumber == 731 || evtNumber == 732 || evtNumber == 737 || evtNumber == 741 || evtNumber == 744 || evtNumber == 748 || evtNumber == 749 || evtNumber == 755 || evtNumber == 760 || evtNumber == 762 || evtNumber == 767 || evtNumber == 769 || evtNumber == 771 || evtNumber == 780 || evtNumber == 784 || evtNumber == 785 || evtNumber == 789 || evtNumber == 796 || evtNumber == 813 || evtNumber == 815 || evtNumber == 818 || evtNumber == 825 || evtNumber == 828 || evtNumber == 843 || evtNumber == 848 || evtNumber == 854 || evtNumber == 862 || evtNumber == 866 || evtNumber == 867 || evtNumber == 868 || evtNumber == 879 || evtNumber == 880 || evtNumber == 886 || evtNumber == 888 || evtNumber == 891 || evtNumber == 894 || evtNumber == 897 || evtNumber == 898 || evtNumber == 903 || evtNumber == 905 || evtNumber == 908 || evtNumber == 910 || evtNumber == 918 || evtNumber == 919 || evtNumber == 922 || evtNumber == 933 || evtNumber == 937 || evtNumber == 938 || evtNumber == 951 || evtNumber == 957 || evtNumber == 958 || evtNumber == 968 || evtNumber == 979 || evtNumber == 980 || evtNumber == 994 || evtNumber == 998) return true;
  return false;
}

void assignValuesBool( std::vector<bool> &target, std::vector<bool> source, unsigned int startPos ) {

  for( unsigned i=0; i<target.size(); ++i ) 
    target[i] = source[startPos+i];

}




std::vector<HodoCluster*> getHodoClusters( std::vector<float> hodo, float fibreWidth, int nClusterMax, float Cut ) {

  std::vector<HodoCluster*> clusters;

  HodoCluster* currentCluster = new HodoCluster( hodo.size(), fibreWidth );

  for( unsigned i=0; i<hodo.size(); ++i ) {

    if( hodo[i] > Cut) { // hit

      if( currentCluster->getSize() < nClusterMax ) {

        currentCluster->addFibre( i );

      } else {

        clusters.push_back( currentCluster ); // store old one
        currentCluster = new HodoCluster( hodo.size(), fibreWidth );   // create a new one
        currentCluster->addFibre( i );        // get that fibre!

      }

    } else { // as soon as you find a hole
      
      if( currentCluster->getSize() > 0 ) {
     
        clusters.push_back( currentCluster ); // store old one
        currentCluster = new HodoCluster( hodo.size(), fibreWidth );   // create a new one

      }

    }


  } // for fibres


  if( currentCluster->getSize()>0 )
    clusters.push_back( currentCluster ); // store last cluster


  return clusters;

}




void doHodoReconstruction( std::vector<float> values, int &nClusters, int *nFibres, float *pos, float fibreWidth, int clusterMaxFibres, float Cut ) {

  std::vector<HodoCluster*> clusters = getHodoClusters( values, fibreWidth, clusterMaxFibres, Cut );

  nClusters = clusters.size();
  for( unsigned i=0; i<clusters.size(); ++i ) {
    nFibres[i] = clusters[i]->getSize();
    pos[i] = clusters[i]->getPosition();
  }

}


void copyArray( int n, float *source, float *target ) {

  for( unsigned i=0; i<n; ++i ) 
    target[i] = source[i];

}













std::vector<HodoCluster*> getHodoClustersBool( std::vector<bool> hodo, float fibreWidth, int nClusterMax, float Cut ) {

  std::vector<HodoCluster*> clusters;

  HodoCluster* currentCluster = new HodoCluster( hodo.size(), fibreWidth );

  for( unsigned i=0; i<hodo.size(); ++i ) {

    if( hodo[i] > Cut) { // hit

      if( currentCluster->getSize() < nClusterMax ) {

        currentCluster->addFibre( i );

      } else {

        clusters.push_back( currentCluster ); // store old one
        currentCluster = new HodoCluster( hodo.size(), fibreWidth );   // create a new one
        currentCluster->addFibre( i );        // get that fibre!

      }

    } else { // as soon as you find a hole
      
      if( currentCluster->getSize() > 0 ) {
     
        clusters.push_back( currentCluster ); // store old one
        currentCluster = new HodoCluster( hodo.size(), fibreWidth );   // create a new one

      }

    }


  } // for fibres


  if( currentCluster->getSize()>0 )
    clusters.push_back( currentCluster ); // store last cluster


  return clusters;

}


void doHodoReconstructionBool( std::vector<bool> values, int &nClusters, int *nFibres, float *pos, float fibreWidth, int clusterMaxFibres, float Cut ) {

  std::vector<HodoCluster*> clusters = getHodoClustersBool( values, fibreWidth, clusterMaxFibres, Cut );

  nClusters = clusters.size();
  for( unsigned i=0; i<clusters.size(); ++i ) {
    nFibres[i] = clusters[i]->getSize();
    pos[i] = clusters[i]->getPosition();
  }

}
