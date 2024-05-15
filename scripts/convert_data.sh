#!/bin/bash

res=512

for ((i=0; i<10; i++)); do
    echo $i
    basename=$(printf "output/line_test_run%i_" $i)

    parallel ~/git/afivo/tools/plot_raw_data.py {} -var sigma -project 0 -axi -inter nearest -min_pixels $res -save_npz {.}_sigmaz_$res.npz ::: "${basename}"*.silo
done

parallel --bar ~/git/afivo/tools/plot_raw_data.py {} -var phi -q -save_npz {.}_phi_$res.npz -min_pixels $res ::: output/line_test_run?_*.silo

parallel --bar ~/git/afivo/tools/plot_raw_data.py {} -var rhs -axi -q -save_npz {.}_rhs_$res.npz -min_pixels $res ::: output/line_test_run?_*.silo

parallel --bar ~/git/afivo/tools/plot_raw_data.py {} -var sigma -axi -q -save_npz {.}_sigma_$res.npz  -min_pixels $res ::: output/line_test_run?_*.silo

