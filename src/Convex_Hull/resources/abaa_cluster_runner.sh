#!/bin/bash

for k in ABAA/*-ABAA/
do
  cd $k
  sbatch *sh
  cd ../../
done