#include "interface/TopologicalSelectionHelper.h"



bool TopologicalSelectionHelper::passesFibreTopologicalSelection(RecoTree &t,bool isNinoRun){
  float  X=(t.cluster_pos_corr_hodoX1+t.cluster_pos_corr_hodoX2)/2.;
  float  Y=(t.cluster_pos_corr_hodoY1+t.cluster_pos_corr_hodoY2)/2.;

  if(!isNinoRun){
    if(Y>-3 && Y<0.5 && X<-2.5 && X>-6 && t.nClusters_hodoX1>0 && t.cluster_pos_corr_hodoX1>-100 && t.cluster_pos_corr_hodoY1>-100&& t.cluster_pos_corr_hodoX2>-100 && t.cluster_pos_corr_hodoY2>-100 && t.wc_x_corr>-20){
      //    if(0.66*X+4.16+Y<0){//line passing through (-5.5,-0.5) and (-2.5,2.5)
      if(X+5.5+1.5*Y<0){//line passing through (-5.5,-0.) and (-2.5,-2)
	return true; 
      }
    }
  }else {
    //    if(Y>1.8 && Y<3.5 && X<0.3 && X>-5 && t.nClusters_hodoX1>0 && t.cluster_pos_corr_hodoX1>-100 && t.cluster_pos_corr_hodoY1>-100&& t.cluster_pos_corr_hodoX2>-100 && t.cluster_pos_corr_hodoY2>-100 && t.wc_x_corr>-20){
    //      if(Y+0.85*X-0.655<0){//line passing through (-3.7,3.8) and (0.3,0.4)
    //    if(Y>1 && Y<2.8 && X<-1.5 && X>-4.5 && t.nClusters_hodoX1>0 && t.cluster_pos_corr_hodoX1>-100 && t.cluster_pos_corr_hodoY1>-100&& t.cluster_pos_corr_hodoX2>-100 && t.cluster_pos_corr_hodoY2>-100 && t.wc_x_corr>-20){
    //      if(Y+2*X+3 < 0){//line passing through (-3,3) and (-2,1)

    if(Y>1 && Y<4 && X<-0.9 && X>-3.9 && t.nClusters_hodoX1>0 && t.cluster_pos_corr_hodoX1>-100 && t.cluster_pos_corr_hodoY1>-100&& t.cluster_pos_corr_hodoX2>-100 && t.cluster_pos_corr_hodoY2>-100 && t.wc_x_corr>-20){
      if(3*Y+2*X-4.2 < 0){//line passing through (-3.9,4) and (-0.9,2)
	return true;
      }
    }
  }

  return false;

}



bool TopologicalSelectionHelper::passesChannelTopologicalSelection(RecoTree &t,bool isNinoRun){
  float  X=(t.cluster_pos_corr_hodoX1+t.cluster_pos_corr_hodoX2)/2.;
  float  Y=(t.cluster_pos_corr_hodoY1+t.cluster_pos_corr_hodoY2)/2.;
  if(!isNinoRun){
    if(Y>-2 && Y<5.5 && X<5.5 && X>-6 && t.nClusters_hodoX1>0 && t.cluster_pos_corr_hodoX1>-100 && t.cluster_pos_corr_hodoY1>-100&& t.cluster_pos_corr_hodoX2>-100 && t.cluster_pos_corr_hodoY2>-100 && t.wc_x_corr>-20){
      //    if(0.66*X+4.16+Y>0){//line passing through (-5.5,-0.5) and (-2.5,2.5)
      if(X+5.5+1.5*Y>0){//line passing through (-5.5,-0.) and (-2.5,-2)
	return true; 
      }
    }
  }else{
    //    if(Y>1.2 && Y<5 && X<5 && X>-3.7 && t.nClusters_hodoX1>0 && t.cluster_pos_corr_hodoX1>-100 && t.cluster_pos_corr_hodoY1>-100&& t.cluster_pos_corr_hodoX2>-100 && t.cluster_pos_corr_hodoY2>-100 && t.wc_x_corr>-20){
    //    if(Y>1 && Y<5 && X<5 && X>-4 && t.nClusters_hodoX1>0 && t.cluster_pos_corr_hodoX1>-100 && t.cluster_pos_corr_hodoY1>-100&& t.cluster_pos_corr_hodoX2>-100 && t.cluster_pos_corr_hodoY2>-100 && t.wc_x_corr>-20){
      //      if(Y+0.85*X-0.655 > 0){//line passing through (-3.7,3.8) and (0.3,0.4)
    //      if(Y+2*X+3 > 0){//line passing through (-3,3) and (-2,1)
    if(Y>2 && Y<5 && X<5 && X>-3.9 && t.nClusters_hodoX1>0 && t.cluster_pos_corr_hodoX1>-100 && t.cluster_pos_corr_hodoY1>-100&& t.cluster_pos_corr_hodoX2>-100 && t.cluster_pos_corr_hodoY2>-100 && t.wc_x_corr>-20){
      if(3*Y+2*X-4.2 > 0){//line passing through (-3.9,4) and (-0.9,2)
	return true;
      }
    }

  } 


  return false;
}


