#include "WaveformUtil.h"
#include "DrawTools.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>


std::vector<std::string> runs;
std::vector<float> beamEnergy;


int main( int argc, char* argv[] ) {

  DrawTools::setStyle();

  TFile *inputFile;

   std::string runName = "";

  if( argc>1 ) {
     std::string runName_str(argv[1]);
     runName = runName_str;
     
     std::string filename_str = "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/micheli/rawData/output_run" + runName + ".root";
     inputFile=TFile::Open(filename_str.c_str());
    
    if( inputFile==0 ) {
      std::cout << "ERROR! Din't find file " << filename_str << std::endl;
      std::cout << "Exiting." << std::endl;
      exit(11);
    }
    
    
  }


  TTree* inputTree=(TTree*)inputFile->Get("outputTree");

  WaveformUtil t(inputTree);
  t.Loop();


  return 0;

}
