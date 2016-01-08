typedef struct timingPlots_Config_t {
  int nMaxAmplCuts;
  int startCut;
  int step;
  int stepAmpl;
  float rangeXLow;
  float rangeXUp;
  std::string setup;  

} timingPlots_Config_t ;

timingPlots_Config_t theConfiguration_;
