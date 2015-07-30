#!/bin/tcsh

  set dir = /shome/micheli/H4_2015/Analysis/rawData/
  set destinationdir = srm://t3se01.psi.ch/pnfs/psi.ch/cms/trivcat/store/user/micheli/rawData/

#  mkdir $i                                                                                                                                                                      
  foreach f (`ls $dir` )
     set file = $dir$f
     echo "copying $file"
    gfal-copy file://$file $destinationdir
  end
#end  
