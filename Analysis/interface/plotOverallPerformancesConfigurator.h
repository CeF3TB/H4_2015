typedef struct plotOverallPerformances_Config_t {
  int nMaxAmplCuts;
  int stepAmplFibre;
  int stepAmplChannel;
  int stepChannel;
  float rangeXLow;
  float rangeXUp;
  std::string setup;  
  float channel2CutFibre;
  float channel1CutChannel;
  float channel2CutChannel;
  std::vector<int> runs;
  float startCutFibre;
  float startCutChannel;
  float amplCut;
  int addTagFileName;
  std::string tagFileName;
  std::vector<float> energies;
} plotOverallPerformances_Config_t ;

plotOverallPerformances_Config_t theConfiguration_;
