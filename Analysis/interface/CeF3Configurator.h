typedef struct CeF3_Config_t {
  std::vector<int> cef3Channels;
  std::vector<int>  syncChannels;
  std::vector<int>  timingThresh;
  int triggerChannel;
  double meanTriggerTime;
  int mcpChannel;
  int mcpTriggerChannel;
  int mcpTimingThresh;
  double mcpMeanTriggerTime;
  int syncMcp;
  std::string waveformFile;
  int startSample;
  int endSample;
  int mcpRun;
  int cherRun;
  std::string cherFile;
  std::vector<float> pedMeanY;
  std::vector<float> pedSigmaY;
  std::vector<float> pedMeanX;
  std::vector<float> pedSigmaX;

} CeF3_Config_t ;

CeF3_Config_t theConfiguration_;
