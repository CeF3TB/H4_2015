#ifndef TopologicalSelectionHelper_h
#define TopologicalSelectionHelper_h
#include "RecoTree.h"

class TopologicalSelectionHelper {

 public:
  static  bool passesFibreTopologicalSelection(RecoTree &t, bool isNinoRun=false);
  static  bool passesChannelTopologicalSelection(RecoTree &t,bool isNinoRun=false);


};


#endif
