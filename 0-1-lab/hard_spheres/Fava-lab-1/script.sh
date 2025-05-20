#!/bin/bash

for rho in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0; do 
# for rho in 0.1 0.2; do 
    echo "&nml density=$rho, nc=3 /" | ./initialize.x

    echo '&nml nblock=10, nstep=1000 burnout=1000 /' | ./mc_nvt_hs.x > temp.txt

    P=$(cat temp.txt | grep 'Run averages' | awk '{print $NF}')
    Rho=$(cat temp.txt | grep 'Density' | awk '{print $NF}')

    echo $Rho $P >> EOS.txt
done