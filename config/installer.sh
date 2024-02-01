#!/usr/bin/env bash


cd ..
rm -rf EspressoPerformanceMonitor
git clone ssh://git@gitlab.cern.ch:7999/lhcb-ft/EspressoPerformanceMonitor.git
cd EspressoPerformanceMonitor
git submodule update --init --recursive

source setenv.sh
mkdir build && cd build
cmake3 ..

make -j

cd ../..

