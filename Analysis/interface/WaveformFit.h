#ifndef WAVEFORM_FIT_H
#define WAVEFORM_FIT_H

#include "TProfile.h"
#include "interface/Waveform.h"
#include "Math/Minimizer.h"
#include "interface/WaveformFit.h"

namespace WaveformFit
{
  void alignWaveform(TProfile* ref_profile, TProfile* fit_profile, ROOT::Math::Minimizer* &minimizer);
  void fitWaveform(Waveform* wave, TProfile* amplitudeProfile,  int samplesBeforeMax, int samplesAfterMax, const Waveform::max_amplitude_informations& max, const Waveform::baseline_informations& rms, ROOT::Math::Minimizer* &minimizer);
  void fitWaveformSimple(Waveform* wave, TProfile* amplitudeProfile,  int samplesBeforeMax, int samplesAfterMax, const Waveform::max_amplitude_informations& max, const Waveform::baseline_informations& rms, ROOT::Math::Minimizer* &minimizer, bool fixedWindow=false, int minSample=0, int MaxSample=700);
  void fitWaveformSimplePlusTime(Waveform* wave, TProfile* amplitudeProfile,  int samplesBeforeMax, int samplesAfterMax, const Waveform::max_amplitude_informations& max, const Waveform::baseline_informations& rms, ROOT::Math::Minimizer* &minimizer, bool fixedWindow=false, int minSample=0, int MaxSample=700);
}
#endif
