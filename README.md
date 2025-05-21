# Evaluating the effect of MIMO feedback quantization and grouping
 
This repository contains the code of an IEEE 802.11 emulator for the evaluation of the impact of beamforming feedback compression through OFDM sub-channel grouping and angle quantization.  

The evaluation has been presented in the Wireless Communications Letter [''Evaluating the Impact of Channel Feedback Quantization and Grouping in IEEE 802.11 MIMO Wi-Fi Networks''](https://ieeexplore.ieee.org/document/10697099).

If you find the project useful and you use this code, please cite our article:
```
@article{meneghello2024evaluating,
  author={Meneghello, Francesca and Haque, Khandaker Foysal and Restuccia, Francesco},
  journal={IEEE Wireless Communications Letters}, 
  title={{Evaluating the Impact of Channel Feedback Quantization and Grouping in IEEE 802.11 MIMO Wi-Fi Networks}}, 
  year={2024},
  volume={13},
  number={12},
  pages={3419--3423}
  }
```

## How to use
Clone the repository and enter the folder with the python code:
```bash
cd <your_path>
git clone https://github.com/francescamen/MIMO_feedback_quantization_grouping.git
```

Download the input data from [http://researchdata.cab.unipd.it/id/eprint/1400](http://researchdata.cab.unipd.it/id/eprint/1400) and unzip the file.
For your convenience, you can use the ```Data``` folder inside this project to place the files but the scripts work whatever is the source folder.

### Obtain the multi-path parameters

To extract the multi-path components from the collected Channel Frequency Response (CFR) samples execute ```Python_code/aoa_toa_method_mDtrack.py``` 
The multi-path components are needed for the emulation. Note that the dataset we made available at [http://researchdata.cab.unipd.it/id/eprint/1400](http://researchdata.cab.unipd.it/id/eprint/1400) already contains also the multi-path decomposition so you can avoid executing this Python file. Instead, if you want to use our code to process a different CFR dataset you need to also execute this file.

### Create the Matlab files for the emulations with the diffrent parameters 

The Matlab files for the emulations and simulations are inside the filder ```Matlab_code/```. Two subfolders are present for the evaluations using the IEEE 802.11ac (```AC```) and the IEEE 802.11ax (```AX```) standards.  

You can find the files for emulations using the multi-path traces extracted from the CFR collected from commercial Wi-Fi devicesAX_simulation (see the details in the reference article) operating with 40 and 80 MHz channels. Matlab files for simulations considering the TG channel models in Matlab are also included to evaluate 160 MHz channels bandwidths.

The files ending with ```_grouping``` are used to evaluate the impact of OFDM sub-channel grouping.

The files ending with ```_corrupted_feedback``` are used to evaluate the impact of the channel impairments in the transmission of the feedback.

Execute the file ```Matlab_code/create_mat_change_quantization_snr.sh``` to create the Matlab files for the different simulations. We considered the following values for the evaluation parameters: 

* quantization_values='0 1 2'
* snr_values='8 10 12 14 16 18 20 22 24'
* mcs_values='3 4'
* standard_names='AC AX'
* folder_environments='Anechoic_1 Anechoic_2 Office Classroom'
* chanBWs='CBW80 CBW40 CBW20'
* NumAntsvect='2 3 4'

The outputs of the emulations will be saved in ```Matlab_code/emulation_results/```

For the evaluations changing the grouping use ```create_mat_files_change_quantization_snr_grouping_AC.sh``` or ```create_mat_files_change_quantization_snr_grouping_AX.sh```.

For the evaluations onsidering also the corruption of the feedback use ```create_mat_files_change_quantization_snr_grouping_AX_corrupted_feedback.sh```.

For the simulations create the files with ```simulation``` instead of ```emulation```.

### Analyze the results with Python

Execute the file ```Python_code/statistics_ber_snr.py``` or ```Python_code/statistics_ber_snr_grouping.py```

arguments:
```NumTxAnts=$1```
```BW=$2```
```mcs=$3```
```index_bf=$4```
```base_folder=$5```
```sub_folder_names=$6```
```name_save_dir=$7```
```standard_names=$8```

e.g.,  Python_code/statistics_ber_snr.py 2 80 3 0 Matlab_code/emulation_results/ anechoic2,class,office Python_code/results/ AC,AX

### Plot the results with Python

Execute the Python files starting with ```figure```. You can also use the bash files starting with ```figure``` that generate the plots for different sets of configuration parameters.
