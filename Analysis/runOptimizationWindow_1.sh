#!/bin/bash

run=$1
echo "running run $run"

for i in {85..115..5}
do 
    cp config/PMTEarlyOctober2015Config_optimization.xml config/PMTEarlyOctober2015Config_optimization_$i.xml
    sed -i "s|startSample>105|startSample>$i|" config/PMTEarlyOctober2015Config_optimization_$i.xml
    sed -i "s|optStart|optStart_$i|" config/PMTEarlyOctober2015Config_optimization_$i.xml
    ./bin/makeAnalysisTree $run V10 PMTEarlyOctober2015Config_optimization_$i; 
done