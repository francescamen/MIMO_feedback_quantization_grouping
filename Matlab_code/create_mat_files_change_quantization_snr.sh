#!/bin/bash

quantization_values='0 1 2'

snr_values='8 10 12 14 16 18 20 22 24'

mcs_values='3 4'

standard_names='AC AX'

folder_environments='Anechoic_1 Anechoic_2 Office Classroom'

chanBWs='CBW80 CBW40 CBW20'

NumAntsvect='2 3 4'

for standard in $standard_names ; do

cd $standard

for NumAnts in $NumAntsvect ; do

for folder_environment in $folder_environments ; do

for chanBW in $chanBWs ; do

for snr_val in $snr_values ; do

for mcs_val in $mcs_values ; do

name_dir="change_quantization/${NumAnts}x${NumAnts}/${folder_environment}/SNR_${snr_val}/MCS_${mcs_val}"
mkdir -p ${name_dir}

fileorig="${standard}_emulation.m"

for quantization_val in $quantization_values ; do

filenew="${name_dir}/${standard}_${chanBW}_${NumAnts}x${NumAnts}_emulation_quantization_${quantization_val}.m"

cp $fileorig $filenew
echo $filenew

# MCS
orig_string="mcsVec      = 3;"
new_string="mcsVec       = ${mcs_val};"
echo $new_string

sed -i -e "s|$orig_string|$new_string|g" $filenew

# Quantization
orig_string="index_bf = 1;"
new_string="index_bf = ${quantization_val};"
echo $new_string

sed -i -e "s|$orig_string|$new_string|g" $filenew

# SNR
orig_string="snr           = 20;"
new_string="snr           = ${snr_val};"
echo $new_string

sed -i -e "s|$orig_string|$new_string|g" $filenew

orig_string="('../"
new_string="('../../../../../../"

sed -i -e "s|$orig_string|$new_string|g" $filenew

orig_string="addpath('utilities/')"
new_string="addpath('../../../../../utilities/')"

sed -i -e "s|$orig_string|$new_string|g" $filenew

# Bandwidth
orig_string="chanBW      = 'CBW80';"
new_string="chanBW      = '${chanBW}';"

sed -i -e "s|$orig_string|$new_string|g" $filenew

# Name folder environment
orig_string="folder_environment = 'Anechoic_2';"
new_string="folder_environment = '${folder_environment}';"

sed -i -e "s|$orig_string|$new_string|g" $filenew

# Name file environment
name_file=${folder_environment//_}

if [ "$name_file" = "Classroom" ]; then
name_file="class"
fi

name_file=${name_file,}

orig_string="name_environment = 'anechoic2';"
new_string="name_environment = '${name_file}';"

sed -i -e "s|$orig_string|$new_string|g" $filenew

# Number of antennas
orig_string="NumTxAnts = 2;"
new_string="NumTxAnts = ${NumAnts};"

sed -i -e "s|$orig_string|$new_string|g" $filenew

orig_string="NumSTSVec = 2;"
new_string="NumSTSVec = ${NumAnts};"

sed -i -e "s|$orig_string|$new_string|g" $filenew


done
done
done
done
done
done
cd ..
done
