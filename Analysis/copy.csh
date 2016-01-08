#!/bin/tcsh

  set dir = /shome/micheli/H4_2015/Analysis/analysisTrees_V10/
  set destinationdir = srm://t3se01.psi.ch/pnfs/psi.ch/cms/trivcat/store/user/micheli/analysis/analysisTrees_V10/

#  mkdir $i                                                                                                                                                                      
  foreach f (`ls $dir|grep Reco|grep root ` )
     set file = $dir$f
     echo "copying $file"
    gfal-copy file://$file $destinationdir
  end
#end  
