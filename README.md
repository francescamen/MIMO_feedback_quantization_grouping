# Evaluating the effect of MIMO feedback quantization and grouping
 
This repository contains the code of an IEEE 802.11 emulator for the evaluation of the impact of beamforming feedback compression through OFDM sub-channel grouping and quantization.  

The evaluation has been presented in the article [''Evaluating the Impact of Channel Feedback Quantization and Grouping in IEEE 802.11 MIMO Wi-Fi Networks''](https://ieeexplore.ieee.org/document/10697099).

If you find the project useful and you use this code, please cite our article:
```
@article{meneghello2022sharp,
  author={Meneghello, Francesca and Haque, Khandaker Foysal and Restuccia, Francesco},
  journal={IEEE Wireless Communications Letters}, 
  title={{Evaluating the Impact of Channel Feedback Quantization and Grouping in IEEE 802.11 MIMO Wi-Fi Networks}}, 
  year={2024},
  volume={},
  number={},
  pages={}
  }
```

## How to use
Clone the repository and enter the folder with the python code:
```bash
cd <your_path>
git clone https://github.com/francescamen/MIMO_feedback_quantization_grouping.git
```

Download the input data from http://researchdata.cab.unipd.it/id/eprint/624 and unzip the file.
For your convenience, you can use the ```Data``` folder inside this project to place the files but the scripts work whatever is the source folder.


### Obtain the ray tracing

Execute ```Python_code/aoa_toa_method_mDtrack.py```

### Create the Matlab files for the emulation 

TBC

### Launch Matlab emulations

TBC

### Analyze the results with Python

TBC

