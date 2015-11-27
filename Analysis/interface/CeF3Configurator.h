typedef struct CeF3_Config_t {
  std::vector<int> cef3Channels;
  std::vector<int>  syncChannels;
  int triggerChannel;
  int mcpChannel;
  int mcpTriggerChannel;
  int syncMcp;
  std::string waveformFile;
  int startSample;
  int endSample;
  int mcpRun;
  int cherRun;
  std::string cherFile;
} CeF3_Config_t ;

CeF3_Config_t theConfiguration_;
