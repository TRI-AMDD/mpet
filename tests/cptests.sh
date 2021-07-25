#!/bin/bash

for dir in test*
do
  refDir="20210204_124704/$dir"

  #Copy configs
  #cp $refDir/*.cfg $dir

  #Copy sim_output
  cp $refDir/sim_output/* $dir/sim_output

done
