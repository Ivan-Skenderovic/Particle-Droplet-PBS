#!/usr/bin/env bash
set -ex

cd precursor_solution_droplet && matlab -batch 'main' 
cd ../
cd suspension_droplet && matlab -batch 'main'