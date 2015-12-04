#!/bin/bash

run=$1
echo "running run $run"
for i in {85..115..5}
do
    ./bin/drawCherenkov $run V10 $i
    mkdir plots_drawCherenkov_opt/
    mkdir plots_drawCherenkov_opt/$run/
    mkdir plots_drawCherenkov_opt/$run/$i
    cp plots_drawCherenkov/*$run* plots_drawCherenkov_opt/$run/$i/
done

./bin/drawCherenkov $run V10