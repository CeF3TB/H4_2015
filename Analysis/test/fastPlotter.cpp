#include "fastPlotter.h"
#include "DrawTools.h"

int main( int argc, char* argv[] ) {

  DrawTools::setStyle();
  TFile *inputFile;
  if( argc>1 ) {
    std::string filename_str(argv[1]);
    inputFile=TFile::Open(filename_str.c_str());
    
    if( inputFile==0 ) {
      std::cout << "ERROR! Din't find file " << filename_str << std::endl;
      std::cout << "Exiting." << std::endl;
      exit(11);
    }
    
    
  }
  
  TTree* inputTree=(TTree*)inputFile->Get("outputTree");

  fastPlotter t(inputTree);
  t.Loop();

  return 0;

}
