#include "interface/WaveformFit.h"
#include "Math/Interpolator.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TMath.h"

#include "assert.h"
#include <iostream>
using namespace std;

#define MAX_INTERPOLATOR_POINTS 10000

namespace WaveformFit
{
  TProfile* refWave;
  TProfile* fitWave;
  Waveform* waveData;
  
  //   ROOT::Math::Minimizer* min =0;

  ROOT::Math::Interpolator* y_interpolator=0;
  ROOT::Math::Interpolator* y_err_interpolator=0;

  //   int nPointsUsed;
  //   int nSamplesToInterpolate;

  int xMin,xMax;
  float sampleRMS;
  
  double chi2(const double *par )
  {
    double chisq = 0;
    double delta = 0;
    for (int i=xMin;i<=xMax;++i)
      {
	//	delta=(refWave->GetBinContent(i)-(y_interpolator->Eval(refWave->GetBinCenter(i)+par[1])+par[0]))/(TMath::Sqrt(refWave->GetBinError(i)*refWave->GetBinError(i)+y_err_interpolator->Eval(refWave->GetBinCenter(i)+par[1])*y_err_interpolator->Eval(refWave->GetBinCenter(i)+par[1])));
	delta= (refWave->GetBinContent(i)-(y_interpolator->Eval(refWave->GetBinCenter(i)+par[1])+par[0])); //using no error for the moment in the Chi2 definition
	chisq += delta*delta;
      }
    return chisq;
  }

  double chi2Wave(const double *par )
  {
    double chisq = 0;
    double delta = 0;
    for (int i=xMin;i<=xMax;++i)
      {
	//	delta=(refWave->GetBinContent(i)-(y_interpolator->Eval(refWave->GetBinCenter(i)+par[1])+par[0]))/(TMath::Sqrt(refWave->GetBinError(i)*refWave->GetBinError(i)+y_err_interpolator->Eval(refWave->GetBinCenter(i)+par[1])*y_err_interpolator->Eval(refWave->GetBinCenter(i)+par[1])));
	//FIXME!!!!! Qui i BOUNDARIES non vanno bene
	//	fitWave->Print();
	//	std::cout<<fitWave->GetXaxis()->GetXmin()<<"<- min"<<fitWave->GetXaxis()->GetXmax()<<std::endl;
	if ( (*waveData)._times[i]*1.e9-par[1]<=fitWave->GetXaxis()->GetXmin()*1.e9+0.25 || (*waveData)._times[i]*1.e9-par[1]>fitWave->GetXaxis()->GetXmax()*1.e9-0.25)
	  {
	    std::cout << "OUT OF BOUNDS " << (*waveData)._times[i]*1.e9-par[1] << "," << fitWave->GetXaxis()->GetXmin() << "," << fitWave->GetXaxis()->GetXmax()*1e9 << std::endl;
	    //	    std::cout<<" BOUNDS"<<(fitWave->GetXaxis()->GetXmin()*1e9)<< " "<<(fitWave->GetXaxis()->GetXmax()*1e9)<<std::endl;
	    chisq += 9999;
	  }
	else
	  {
	    //	    delta= ((*waveData)._samples[i]-(y_interpolator->Eval((*waveData)._times[i]*1.e9-par[1])*par[0]))/sampleRMS; //fit par[0]*ref_shape(t-par[1]) par[0]=amplitude, par[1]=DeltaT
	    delta= ((*waveData)._samples[i]-(fitWave->GetBinContent(fitWave->FindBin((*waveData)._times[i]*1.e9-par[1])))*par[0])/sampleRMS; //fit par[0]*ref_shape(t-par[1]) par[0]=amplitude, par[1]=DeltaT
	    chisq += delta*delta;
	  }
      }
    return chisq;
  }

  double chi2WaveSimple(const double *par )
  {
    double chisq = 0;
    double delta = 0;
    for (int i=xMin;i<=xMax;++i)
      {
	//	delta= ((*waveData)._samples[i]-(fitWave->GetBinContent(fitWave->FindBin((*waveData)._times[i])))*par[0])/sampleRMS; //fit par[0]*ref_shape(t) par[0]=amplitude
	//	if(xMax>998)	std::cout<<i<<"<-i "<<(*waveData)._times[i]<<"<-eval"<<std::endl;
	delta= ((*waveData)._samples[i]-(y_interpolator->Eval((*waveData)._times[i]))*par[0])/sampleRMS; //fit par[0]*ref_shape(t) par[0]=amplitude
	    chisq += delta*delta;
      }
    return chisq;
  }

  double chi2WaveSimplePlusTime(const double *par )
  {
    double chisq = 0;
    double delta = 0;
    for (int i=xMin;i<=xMax;++i)
      {
	//	delta= ((*waveData)._samples[i]-(fitWave->GetBinContent(fitWave->FindBin((*waveData)._times[i])))*par[0])/sampleRMS; //fit par[0]*ref_shape(t) par[0]=amplitude
	//	if(xMax>998)	std::cout<<i<<"<-i "<<(*waveData)._times[i]<<"<-eval"<<std::endl;
	delta= ((*waveData)._samples[i]-(y_interpolator->Eval((*waveData)._times[i]-par[1]))*par[0])/sampleRMS; //fit par[0]*ref_shape(t-par[1]) par[0]=amplitude
	    chisq += delta*delta;
      }
    return chisq;
  }



  void residuals(const double *par )
  {
    std::cout << "++++++++++++++++++++++++++++++++++" << std::endl;
    for (int i=xMin;i<=xMax;++i)
      {
	float res = ((*waveData)._samples[i]-(y_interpolator->Eval((*waveData)._times[i]*1.e9-par[1])*par[0]))/sampleRMS; //fit par[0]*ref_shape(t-par[1]) par[0]=amplitude, par[1]=DeltaT
	std::cout << "===> " << (i-xMin) << " " << (*waveData)._samples[i] << " " << res << std::endl;
      }
  }

  void alignWaveform(TProfile* ref_profile, TProfile* fit_profile, ROOT::Math::Minimizer* &minimizer)
  {    
    xMin=1;
    xMax=990;

    refWave=ref_profile;
    fitWave=fit_profile;

    minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

    minimizer->SetMaxFunctionCalls(1000000);
    minimizer->SetMaxIterations(100000);
    minimizer->SetTolerance(1e-3);
    minimizer->SetPrintLevel(2);

    if (!y_interpolator)
      y_interpolator=new ROOT::Math::Interpolator(MAX_INTERPOLATOR_POINTS, ROOT::Math::Interpolation::kCSPLINE);
    if (!y_err_interpolator)
      y_err_interpolator=new ROOT::Math::Interpolator(MAX_INTERPOLATOR_POINTS, ROOT::Math::Interpolation::kCSPLINE);

    std::vector<double> x,y,y_err;

    for (int i=1;i<=fit_profile->GetNbinsX();++i)
      {
	x.push_back(fit_profile->GetBinCenter(i));
	y.push_back(fit_profile->GetBinContent(i));
	y_err.push_back(fit_profile->GetBinError(i));
      }

    y_interpolator->SetData(x,y);
    y_err_interpolator->SetData(x,y_err);

    ROOT::Math::Functor f(&chi2,2);
    minimizer->SetFunction(f);

    minimizer->SetVariable(0,"deltaV",0.,1e-2);
    minimizer->SetVariable(1,"deltaT",0.,1e-1);
    minimizer->Minimize();
   
    const double* par=minimizer->X();
    std::cout << "+++++ FIT RESULT: " << par[0] << "," << par[1] << std::endl;
    delete y_interpolator;
    delete y_err_interpolator;
    y_interpolator=0;
    y_err_interpolator=0;
    //     minimizer=min;
  }

  void fitWaveform(Waveform* wave, TProfile* amplitudeProfile, int nSamplesBeforeMax, int nSamplesAfterMax, const Waveform::max_amplitude_informations& max, const Waveform::baseline_informations& waveRms, ROOT::Math::Minimizer* &minimizer)
  {
    xMin=max.sample_at_max-nSamplesBeforeMax;
    xMax=max.sample_at_max+nSamplesAfterMax;

    waveData=wave;
    fitWave=amplitudeProfile;

    sampleRMS=waveRms.rms;

    minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

    minimizer->SetMaxFunctionCalls(100000);
    minimizer->SetMaxIterations(100);
    minimizer->SetTolerance(1e-3);
    minimizer->SetPrintLevel(0);


    if (!y_interpolator)
      y_interpolator=new ROOT::Math::Interpolator(MAX_INTERPOLATOR_POINTS, ROOT::Math::Interpolation::kCSPLINE);

    std::vector<double> x,y;
    for (int i=1;i<=amplitudeProfile->GetNbinsX();++i)
      {
	x.push_back(amplitudeProfile->GetBinCenter(i));
	y.push_back(amplitudeProfile->GetBinContent(i));
	//	std::cout<<"i="<<i<<" "<<amplitudeProfile->GetBinCenter(i)<<" "<<amplitudeProfile->GetBinContent(i)<<std::endl;
      }

    y_interpolator->SetData(x,y);

    ROOT::Math::Functor f(&chi2Wave,2);
    minimizer->SetFunction(f);

    //    std::cout << "INIITIAL VALUES " << max.max_amplitude << "," << max.time_at_max*1.e9 << "," << max.sample_at_max << std::endl;

    minimizer->SetLimitedVariable(0,"amplitude",max.max_amplitude,1e-2,std::max(0.,(double)max.max_amplitude-500),std::min(4000.,(double)max.max_amplitude+500.));
    minimizer->SetLimitedVariable(1,"deltaT",max.time_at_max*1.e9,1e-3,max.time_at_max*1.e9-0.5,max.time_at_max*1.e9+0.5);

    minimizer->Minimize();

     if (minimizer->Status()==0)
       {
 	const double* par=minimizer->X();

 	 std::cout << "+++++ FIT RESULT: " << par[0] << "," << par[1] << std::endl;
	 // 	residuals(par);
       }
    
    delete y_interpolator;
    y_interpolator=0;
  }

  void fitWaveformSimple(Waveform* wave, TProfile* amplitudeProfile, int nSamplesBeforeMax, int nSamplesAfterMax, const Waveform::max_amplitude_informations& max, const Waveform::baseline_informations& waveRms, ROOT::Math::Minimizer* &minimizer, bool fixedWindow, int startSample, int endSample)
  {
    if(!fixedWindow){
      xMin=std::max(10,max.sample_at_max-nSamplesBeforeMax);
      xMax=std::min(max.sample_at_max+nSamplesAfterMax,999);
    }else{
      //      xMin=170;
      //      xMax=700;
      xMin=startSample;
      xMax=endSample;
    }
    waveData=wave;
    fitWave=amplitudeProfile;

    sampleRMS=waveRms.rms;

    minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

    minimizer->SetMaxFunctionCalls(100000);
    minimizer->SetMaxIterations(100);
    minimizer->SetTolerance(1e-3);
    minimizer->SetPrintLevel(0);

    if (!y_interpolator)
      y_interpolator=new ROOT::Math::Interpolator(MAX_INTERPOLATOR_POINTS, ROOT::Math::Interpolation::kCSPLINE);

    std::vector<double> x,y;
    for (int i=1;i<=amplitudeProfile->GetNbinsX();++i)
      {
	x.push_back(amplitudeProfile->GetBinCenter(i));
	y.push_back(amplitudeProfile->GetBinContent(i));
	//	std::cout<<"i="<<i<<" "<<amplitudeProfile->GetBinCenter(i)<<" "<<amplitudeProfile->GetBinContent(i)<<std::endl;
      }

    y_interpolator->SetData(x,y);


    ROOT::Math::Functor f(&chi2WaveSimple,1);
    minimizer->SetFunction(f);

    //    std::cout << "INIITIAL VALUES " << max.max_amplitude << "," << max.time_at_max*1.e9 << "," << max.sample_at_max << " RMS:"<<sampleRMS<<std::endl;

    //    minimizer->SetLimitedVariable(0,"amplitude",max.max_amplitude,1e-2,std::max(0.,(double)max.max_amplitude-500),std::min(4000.,(double)max.max_amplitude+4000.));
    minimizer->SetLimitedVariable(0,"amplitude",max.max_amplitude,1e-2,-500,4000);

    minimizer->Minimize();

    if (minimizer->Status()==0||minimizer->Status()==3)
       {
	 //	 	const double* par=minimizer->X();
	 
	 // 	 std::cout << "+++++ FIT RESULT: " << par[0] << std::endl;
	 // 	residuals(par);
       }

     delete y_interpolator;
     y_interpolator=0;
    

  }

  void fitWaveformSimplePlusTime(Waveform* wave, TProfile* amplitudeProfile, int nSamplesBeforeMax, int nSamplesAfterMax, const Waveform::max_amplitude_informations& max, const Waveform::baseline_informations& waveRms, ROOT::Math::Minimizer* &minimizer, bool fixedWindow, int startSample, int endSample)
  {
    if(!fixedWindow){
      xMin=std::max(1,max.sample_at_max-nSamplesBeforeMax);
      xMax=std::min(max.sample_at_max+nSamplesAfterMax,999);
    }else{
      //      xMin=170;
      //      xMax=700;
      xMin=startSample;
      xMax=endSample;
    }
    waveData=wave;
    fitWave=amplitudeProfile;

    //    std::cout<<"xMin:"<<xMin<<" xMAx:"<<xMax<<std::endl;

    sampleRMS=waveRms.rms;

    minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

    minimizer->SetMaxFunctionCalls(100000);
    minimizer->SetMaxIterations(1000);
    minimizer->SetTolerance(1e-3);
    minimizer->SetPrintLevel(0);

    if (!y_interpolator)
      y_interpolator=new ROOT::Math::Interpolator(MAX_INTERPOLATOR_POINTS, ROOT::Math::Interpolation::kCSPLINE);

    std::vector<double> x,y;
    for (int i=1;i<=amplitudeProfile->GetNbinsX();++i)
      {
	x.push_back(amplitudeProfile->GetBinCenter(i));
	y.push_back(amplitudeProfile->GetBinContent(i));
	//	std::cout<<"i="<<i<<" "<<amplitudeProfile->GetBinCenter(i)<<" "<<amplitudeProfile->GetBinContent(i)<<std::endl;
      }

    y_interpolator->SetData(x,y);


    ROOT::Math::Functor f(&chi2WaveSimplePlusTime,2);
    minimizer->SetFunction(f);

    //    std::cout << "INIITIAL VALUES " << max.max_amplitude << "," << max.time_at_max*1.e9 << "," << max.sample_at_max << " RMS:"<<sampleRMS<<std::endl;

    //    minimizer->SetLimitedVariable(0,"amplitude",max.max_amplitude,1e-2,std::max(0.,(double)max.max_amplitude-500),std::min(4000.,(double)max.max_amplitude+4000.));
    minimizer->SetLimitedVariable(0,"amplitude",max.max_amplitude,1e-2,-500,4000);
    minimizer->SetLimitedVariable(1,"deltaT",0.e-9,1e-13,-10.e-9,10.e-9);

    minimizer->Minimize();

    if (minimizer->Status()==0||minimizer->Status()==3)
       {
	 //	 	 	const double* par=minimizer->X();
	 
	 //	 std::cout << "+++++ FIT RESULT: " << par[0] << ","<<par[1]<< std::endl;
	 // 	residuals(par);
       }

     delete y_interpolator;
     y_interpolator=0;
    

  }


  //----------Get leading edge time at a given threshold and in a given range---------------
  std::pair<float, float> GetTimeLE(Waveform* wave,const Waveform::baseline_informations& waveRms,float thr, int nmFitSamples, int npFitSamples, int min, int max, float timeUnit){

    float leThr=-1;
    float leSample=-1;
    float leTime=-100;
    float chi2le=-100;


    if(thr != leThr || leSample != -1)
      {
        //---find first sample above thr
        leThr = thr;
	if(max>(*wave)._samples.size()) return make_pair(-100,-100);
        for(int iSample=min; iSample<max; ++iSample)
	  {
	    //	    std::cout<<iSample<<" "<<(*wave)._samples[iSample]<<std::endl;
            if((*wave)._samples[iSample] > leThr) 
	      {
                leSample = iSample;
                break;
	      }
	  }
        //---interpolate -- A+Bx = amp
	if(leSample<1) return make_pair(-100,-100);
	float A=0, B=0;
        chi2le = LinearInterpolation(wave,waveRms, A, B, leSample-nmFitSamples, leSample+npFitSamples, timeUnit);
	if(B!=0)        leTime = (leThr - A) / B;
      }

    return make_pair(leTime, chi2le);

  }


  //----------Get the signal integral around the the max-------------------------------------
  std::pair<float,float> GetSignalIntegral(Waveform* wave,int thr, int min=-1)
  {
    //---compute position of the max
    if(min==-1)
      min = 1;

    //---find pulse borders
    int begin=-1, end=-1;
    for(int iSample=min; iSample<(*wave)._samples.size(); ++iSample)
      {
        if(begin==-1 && (*wave)._samples[iSample]>thr)
	  begin = iSample;
        else if(begin!=-1 && (*wave)._samples[iSample]<thr)
	  {
            end = iSample;
            break;
	  }
      }            
            
    //---compute integral

    if(end<0) return make_pair(-10.,-10.);
    float integral=0;
    for(int iSample=begin; iSample<end; ++iSample)
      {
        //---if signal window goes out of bound return a bad value
        if(iSample > (*wave)._samples.size() || iSample < 0)
	  return make_pair(-10.,-10.); ;        
        integral += (*wave)._samples[iSample];
      }

    return make_pair(integral,end-begin);
  }




  //----------Linear interpolation util-----------------------------------------------------
  float LinearInterpolation(Waveform* wave,const Waveform::baseline_informations& waveRms, float& A, float& B, const int& min, const int& max, float timeUnit)
  {
    //fuck ROOT, i will use the analytical formula for linear fits since it's faster
    //---definitions---
    float xx= 0.;
    float xy= 0.;
    float Sx = 0.;
    float Sy = 0.;
    float Sxx = 0.;
    float Sxy = 0.;

    //---compute sums
    int usedSamples=0;
    for(int iSample=min; iSample<=max; ++iSample)
      {
        if(iSample<0 || iSample>=(*wave)._samples.size()) 
	  continue;
        xx = iSample*iSample*timeUnit*timeUnit;
        xy = iSample*timeUnit*(*wave)._samples[iSample];
        Sx = Sx + (iSample)*timeUnit;
        Sy = Sy + (*wave)._samples[iSample];
        Sxx = Sxx + xx;
        Sxy = Sxy + xy;
        ++usedSamples;
      }
    
    float Delta = usedSamples*Sxx - Sx*Sx;
    A = (Sxx*Sy - Sx*Sxy) / Delta;
    B = (usedSamples*Sxy - Sx*Sy) / Delta;

    //---compute chi2---
    float chi2=0;
    float sigma2 = pow(waveRms.rms, 2);
    for(int iSample=min; iSample<=max; ++iSample)
      {
        if(iSample<0 || iSample>=(*wave)._samples.size()) 
	  continue;
	//	std::cout<<iSample*timeUnit*1.e9<<" "<<(*wave)._samples[iSample]<<" B:"<<B<<std::endl;
	//	std::cout<<iSample<<" "<<chi2<<" "<<(*wave)._samples[iSample]<<" "<<A<<" "<<B<<" "<<waveRms.rms<<" "<<timeUnit<<" min:"<<min<<" max:"<<max<<std::endl;
        chi2 = chi2 + pow((*wave)._samples[iSample] - A - B*iSample*timeUnit, 2)/sigma2;
	//	std::cout<<iSample<<" "<<chi2<<std::endl;
      } 

    return chi2/(usedSamples-2);
  }


  
};

