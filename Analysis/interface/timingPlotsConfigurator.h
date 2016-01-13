typedef struct timingPlots_Config_t {
  int nMaxAmplCuts;
  int startCutFibre;
  int stepFibre;
  int stepAmplFibre;
  int startCutChannel;
  int stepChannel;
  int stepAmplChannel;
  float rangeXLow;
  float rangeXUp;
  std::string setup;  
  float channel2CutFibre;
  float channel1CutChannel;
  float channel2CutChannel;
  int nBinsFibre;
  int nBinsChannel;
} timingPlots_Config_t ;

timingPlots_Config_t theConfiguration_;
