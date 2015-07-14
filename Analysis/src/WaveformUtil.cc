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
#include "RooFFTConvPdf.h"
#include "RooPlot.h"
#include "TVectorD.h"

#include "interface/channelInfo.h"
#include "interface/HodoCluster.h"



void assignValues( std::vector<float> &target, std::vector<float> source, unsigned int startPos );
void assignValuesBool( std::vector<bool> &target, std::vector<bool> source, unsigned int startPos );
void doHodoReconstructionBool( std::vector<bool> values, int &nClusters, int *nFibres, float *pos, float fibreWidth, int clusterMaxFibres, float Cut );
std::vector<HodoCluster*> getHodoClustersBool( std::vector<bool> hodo, float fibreWidth, int nClusterMax, float Cut );
void doHodoReconstruction( std::vector<float> values, int &nClusters, int *nFibres, float *pos, float fibreWidth, int clusterMaxFibres, float Cut );
std::vector<HodoCluster*> getHodoClusters( std::vector<float> hodo, float fibreWidth, int nClusterMax, float Cut );
void copyArray( int n, float *source, float *target );



void WaveformUtil::Loop(){

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();
  std::cout<<"nentries"<<nentries<<std::endl;
  //    nentries=7000;
  float  mean[NFIBERS][NDIGISAMPLES];
  float  time[NDIGISAMPLES];

  for (int i=0;i<NFIBERS;++i){
    for (int j=0;j<NDIGISAMPLES;++j){
      mean[i][j]=0;
    }
  }


  TString runNumberString;
  int digiFreq;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;      
    if(jentry==0){
      runNumberString.Form("%d",runNumber);
      digiFreq=digi_frequency;
    }
    if(jentry%1000 == 0)std::cout<<"Processing entry:"<<jentry<<std::endl;
    if(passesHodoSelection()==false)continue;

    for (int i=0;i<1024*4;++i){
      //	std::cout<<"channel:"<<digi_value_ch->at(i)<<digi_value->at(i)<<" "<<digi_value_time->at(i)<<std::endl;
      if(digi_value_ch->at(i) > 3)continue;
      if(digi_max_amplitude->at(digi_value_ch->at(i))>10000 || digi_max_amplitude->at(digi_value_ch->at(i))<0)continue;
      mean[digi_value_ch->at(i)][i-1024*digi_value_ch->at(i)]+=(float)(digi_value->at(i)/nentries);
      if(i<1024)time[i]=digi_value_time->at(i);
    }
  }

  TFile* outFile = TFile::Open("outWaveFormUtil_"+runNumberString+".root","recreate");
  int lowRange[4],highRange[4];
  TVectorD integrals(4);

  lowRange[0]=130;
  lowRange[1]=170;
  lowRange[2]=145;
  lowRange[3]=170;
  highRange[0]=900;
  highRange[1]=720;
  highRange[2]=500;
  highRange[3]=720;

  for (int i=0;i<NFIBERS;++i){
    meanWaveGraphs[i]=new TGraph(1024, time, mean[i]);
    TString fiber;
    fiber.Form("%d",i);
    meanWaveGraphs[i]->SetName("waveform_"+fiber);

    meanWaveHistos[i]=new TH1F("waveform_histo_"+fiber,"",1024,0,time[1023]);
    for (int j=0;j<NDIGISAMPLES;++j){  
      for(int k=0;k<mean[i][j];k++){
	meanWaveHistos[i]->Fill(time[j]);
      }
    }
    meanWaveHistos[i]->Scale(nentries);
    meanWaveHistos[i]->Sumw2();
    meanWaveHistos[i]->Scale(1./nentries);//this trick is just to get correct errors
    meanWaveHistos[i]->Write();
    meanWaveGraphs[i]->Write();


    RooRealVar t("t","time", time[0], time[1023]);
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
    lxg.fitTo(data,RooFit::Range(time[lowRange[i]],time[highRange[i]])) ;

    // Plot data, landau pdf, landau (X) gauss pdf
    RooPlot* frame = t.frame(RooFit::Title("fiber_"+fiber)) ;
    t.setBins(10);
    data.plotOn(frame,RooFit::DataError(RooAbsData::SumW2)) ;

    lxg.plotOn(frame) ;

    TCanvas c1("c_"+fiber);
    frame->Draw();
    c1.Write();

    TH1* histopdf=lxg.createHistogram("t",1024);
    histopdf->SetName("fitted_histo_"+fiber);
    float finalFastSample=230;
    if(digiFreq==1)finalFastSample=185;
    integrals[i]=histopdf->Integral(4,finalFastSample);//the pdf is already normalized to 1
    histopdf->Write();
  }

  integrals.Write("integrals_fast");


  outFile->Write();
  outFile->Close();

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
