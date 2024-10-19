#!/bin/bash

mcs_values='3 4'

chanBWs='80 40'  # 'CBW80 CBW40 CBW20'

NumAntsvect='2 3 4'  # '2 3 4'

base_folder='Python_code/results/'

cd ..

# SUBSTITUTE WITH THE PATH WHERE PYTHON IS INSTALLED
source myenv/bin/activate

for NumAnts in $NumAntsvect ; do

for chanBW in $chanBWs ; do

for mcs_val in $mcs_values ; do

python Python_code/figure_ber_quant_snr_standards.py $NumAnts $chanBW $mcs_val

done
done
done

deactivate
