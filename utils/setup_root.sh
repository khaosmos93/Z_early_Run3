#!/bin/bash

source /cvmfs/sft.cern.ch/lcg/views/LCG_96bpython3/x86_64-centos7-gcc9-opt/setup.sh
source /home/wunsch/workspace/root/build_own_python/bin/thisroot.sh

alias python=/home/wunsch/workspace/python/install/bin/python3.6
export PYTHONPATH=$PWD:$PYTHONPATH
