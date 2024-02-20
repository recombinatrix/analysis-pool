#!/bin/bash

# use this script to run the calculation
printf "source contactFreq.tcl \n source contactFreq.sh \n exit \n"|vmd path/to/grofile.gro path/to/trjfile.xtc -dispdev text
