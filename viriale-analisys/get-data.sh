#!/bin/bash

if [ 0 == 1 ]; then # Oxygen

    PLow=0
    # PHigh=10
    # PInc=0.2
    PHigh=1
    PInc=0.02

    # for T in 200 225 250 275 300 350 400 450 500 550 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000 2100 2200; do
    for T in 100 105 110 115 125 150 175; do

        filename="https://webbook.nist.gov/cgi/fluid.cgi?Action=Data&Wide=on&ID=C7782447&Type=IsoTherm&Digits=5&PLow=$PLow&PHigh=$PHigh&PInc=$PInc&T=$T&RefState=DEF&TUnit=K&PUnit=MPa&DUnit=mol%2Fm3&HUnit=kJ%2Fmol&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm"
        datafile=Oxygen/data/data_T$T.txt

        wget $filename -O $datafile
    done
fi


if [ True ]; then # Methane

    PLow=0

    PHigh=0.1
    PInc=0.001


    # for i in 0 1 2 3 4 5 6 7 8 9 10 11 12 13; do
        # T=$((30*$i+210))
        
    for T in 91 95 100 110 120 140 160 190 195; do

        filename="https://webbook.nist.gov/cgi/fluid.cgi?Action=Data&Wide=on&ID=C74828&Type=IsoTherm&Digits=5&PLow=$PLow&PHigh=$PHigh&PInc=$PInc&T=$T&RefState=DEF&TUnit=K&PUnit=MPa&DUnit=mol%2Fm3&HUnit=kJ%2Fmol&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm"
        datafile=Methane/data/data_T$T.txt

        wget $filename -O $datafile
    done

fi
