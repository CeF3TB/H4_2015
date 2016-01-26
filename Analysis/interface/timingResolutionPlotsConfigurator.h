typedef struct timingResolutionPlots_Config_t {
  std::string setup;  
  std::vector<float> energies;
  std::vector<int> wp_channel;
  std::vector<int> wp_fibre;
  std::vector<int> runs;
} timingResolutionPlots_Config_t ;

timingResolutionPlots_Config_t theConfiguration_;
