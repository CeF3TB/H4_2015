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



void assignValues( std::vector<float> &target, std::vector<float> source, unsigned int startPos );
void assignValuesBool( std::vector<bool> &target, std::vector<bool> source, unsigned int startPos );
void doHodoReconstructionBool( std::vector<bool> values, int &nClusters, int *nFibres, float *pos, float fibreWidth, int clusterMaxFibres, float Cut );
std::vector<HodoCluster*> getHodoClustersBool( std::vector<bool> hodo, float fibreWidth, int nClusterMax, float Cut );
void doHodoReconstruction( std::vector<float> values, int &nClusters, int *nFibres, float *pos, float fibreWidth, int clusterMaxFibres, float Cut );
std::vector<HodoCluster*> getHodoClusters( std::vector<float> hodo, float fibreWidth, int nClusterMax, float Cut );
void copyArray( int n, float *source, float *target );
double funcCRRC(double *x, double *par);


void WaveformUtil::Loop(){

  DrawTools::setStyle();

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();
  std::cout<<"nentries"<<nentries<<std::endl;
  nentries=100;
  float  mean[NFIBERS][NDIGISAMPLES];
  float  time[NDIGISAMPLES];
  float meanTimeAtMax[NFIBERS];
  
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


    //    float timeOfTheEvent=digi_time_at_frac50_bare_noise_sub->at(8);//synchronizing time of events with time of trigger
    float timeOfTheEvent=digi_time_at_1000_bare_noise_sub->at(8);//synchronizing time of events with time of trigger
    //    float shiftTime=190.3-timeOfTheEvent;//mean fitted on trigger run 2778
    float shiftTime=0;
    if (digiFreq==1)shiftTime=190.2-timeOfTheEvent;//mean fitted on trigger run 2778
    if (digiFreq==0)shiftTime=161.9-timeOfTheEvent;//mean fitted on trigger run 3076
    //    float shiftTime=138.1-timeOfTheEvent;//mean fitted on trigger run 329 for 2014 data
    int shiftSample=round(shiftTime/(1e9*timeSampleUnit(digiFreq)));
    shiftSample=-shiftSample;
    shiftSampleHisto->Fill(shiftSample);     
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
      mean[iChannel][i-1024*digi_value_ch->at(i)]+=(float)(digi_value_bare_noise_sub->at(iSample)/nentries);
      if(i<1024)time[i]=digi_value_time->at(i);
      
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
    meanTimeAtMax[i]/=nentries;
    //std::cout<<meanTimeAtMax[i]<< " unit:"<<timeSampleUnit(digiFreq)<<std::endl;
    meanWaveHistosForPlots[i]=new TH1F("waveform_histo_forplots_"+fiber,"",1024,-meanTimeAtMax[i],1e9*time[1023-(int)(meanTimeAtMax[i]*1.e-9/timeSampleUnit(digiFreq))]);

    for (int j=0;j<NDIGISAMPLES;++j){  
      for(int k=0;k<mean[i][j];k++){

	meanWaveHistos[i]->Fill(time[j]);
	meanWaveHistosForPlots[i]->Fill(1e9*(time[j])-meanTimeAtMax[i]);
      }
    }


    meanWaveHistos[i]->Scale(nentries);
    meanWaveHistos[i]->Sumw2();
    meanWaveHistos[i]->Scale(1./nentries);//this trick is just to get correct errors
    meanWaveHistosForPlots[i]->Scale(nentries);
    meanWaveHistosForPlots[i]->Sumw2();
    meanWaveHistosForPlots[i]->Scale(1./nentries);//this trick is just to get correct errors
    meanWaveHistosForPlots[i]->Scale(0.25);//we need millivolts

    meanWaveHistos[i]->Write();
    meanWaveHistosForPlots[i]->Write();
    meanWaveGraphs[i]->Write();

    TCanvas can;
    gStyle->SetPadRightMargin(0.15);
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
    TPaveText* pave = DrawTools::getLabelTop("100 GeV Electron Beam");
    meanWaveHistosForPlots[i]->GetXaxis()->SetTitle("Time [ns]");
    meanWaveHistosForPlots[i]->GetYaxis()->SetTitle("Signal Amplitude [mV]");
    std::cout<<meanTimeAtMax[i]<<std::endl;
    //      meanWaveHistosForPlots[i]->GetXaxis()->SetRangeUser(-meanTimeAtMax[i],1e9*time[1023]-meanTimeAtMax[i]));
    
    meanWaveHistosForPlots[i]->Draw("histl");
    pave->Draw("same");
    c1.Write("pulseShape"+fiber);
    c1.SaveAs("plots/pulseShape"+fiber+".png");
    c1.SaveAs("plots/pulseShape"+fiber+".pdf");
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

  int nClusters_hodoY1;
  int nFibres_hodoY1[HODOY1_CHANNELS]; 
  float pos_hodoY1[HODOY1_CHANNELS];   

  int nClusters_hodoX2;
  int nFibres_hodoX2[HODOX2_CHANNELS]; 
  float pos_hodoX2[HODOX2_CHANNELS];   


  int nClusters_hodoY2;
  int nFibres_hodoY2[HODOY2_CHANNELS]; 
  float pos_hodoY2[HODOY2_CHANNELS];   

  
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
