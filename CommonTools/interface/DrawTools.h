#ifndef DrawTools_h
#define DrawTools_h


#include "TStyle.h"
#include "TPaveText.h"
#include "TColor.h"


class DrawTools {

 public:

  static TStyle* setStyle();

  static TPaveText* getLabelTop( const std::string& text="H4 Test Beam 2015" );
  static TPaveText* getLabelRun( const std::string& runName, bool top=true );
  static TPaveText* getLabelTop_expOnXaxis( const std::string& text="H4 Test Beam 2014" );

};


#endif
