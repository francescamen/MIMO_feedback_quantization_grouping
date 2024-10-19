import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.io as sio
import pickle
import matplotlib.cbook as cbook
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D

mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Palatino'
#mpl.rcParams['text.usetex'] = 'true'
#mpl.rcParams['text.latex.preamble'] = [r'\usepackage{newtxmath}']
mpl.rcParams['font.size'] = 16
mpl.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.cm.Accent.colors)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('NumTxAnts', help='number of antennas', type=int)
    parser.add_argument('BW', help='bandwidth', type=int)
    parser.add_argument('mcs', help='MCS', type=int)
    parser.add_argument('index_bf', help='quantization level in {0,1,2}', type=int)
    parser.add_argument('base_folder', help='folder where the results are located')
    parser.add_argument('sub_folder_names', help='subfolders of the environments, comma separated')
    parser.add_argument('name_save_dir', help='folder where to save the results')
    parser.add_argument('standard_names', help='name of the standards, AX, AC')
    args = parser.parse_args()

    base_folder = args.base_folder  # 'Matlab_code/emulation_results/'
    sub_folder_names = args.sub_folder_names  # ['anechoic1','anechoic2','class','office']
    sub_folders = []
    for subf in sub_folder_names.split(','):
        sub_folders.append(subf)

    snr_values = [14, 16, 18, 20, 22]  # [8, 10, 12, 14, 16, 18, 20, 22, 24]
    snr_values = np.asarray(snr_values, dtype=int)

    NumTxAnts = args.NumTxAnts
    NumSTSVec = NumTxAnts
    NumRxAntsVec = NumSTSVec
    chanBW = args.BW
    mcsVec = args.mcs
    index_bf = args.index_bf
    index_grouping = 0

    end_time = 3000  # 19000
    end_iter = 4
    num_entries = end_time * end_iter

    standard_names = args.standard_names  # ['AC', 'AX']
    standard_values = []
    for stand in standard_names.split(','):
        standard_values.append(stand)

    for standard_name in standard_values:

        name_save_dir = args.name_save_dir + standard_name # 'Python_code/results/'
        if not os.path.isdir(name_save_dir):
            os.mkdir(name_save_dir)

        for sub_fold in sub_folders:

            name_save_subdir = name_save_dir + '/' + sub_fold
            if not os.path.isdir(name_save_subdir):
                os.mkdir(name_save_subdir)

            name_save_file = (name_save_subdir +
                              '/V_diff_statistics_NumTxAnts_' + str(NumTxAnts) +
                              '_NumRxAntsVec_' + str(NumTxAnts) + '_BW_' + str(chanBW) +
                              '_idx_grouping_' + str(index_grouping) + '_idx_bf_' + str(index_bf) +
                              '_mcs_' + str(mcsVec) + '.txt')
            print(name_save_file)
            # if os.path.exists(name_save_file):
            #     continue

            ber_list = []
            V_mat_list = []
            V_mat_rec_list = []
            V_error_list = []
            packetSNR_list = []
            length_list = []

            snr_values_present = []
            for snr_val in snr_values:
                name_folder = (base_folder + standard_name + '/' + sub_fold + '/NumTxAnts_' + str(NumTxAnts) +
                               '_NumRxAntsVec_' + str(NumTxAnts) + '_BW_' + str(chanBW) +
                               '_idx_grouping_' + str(index_grouping) + '_idx_bf_' + str(index_bf) +
                               '_snr_' + str(snr_val) +
                               '_mcs_' + str(mcsVec))
                try:
                    end_time = len(os.listdir(name_folder))
                except FileNotFoundError:
                    print('Folder ' + name_folder + ' does not exist.')
                    ber_list.append([])
                    length_list.append(0)
                    V_mat_list.append([])
                    V_mat_rec_list.append([])
                    V_error_list.append([])
                    packetSNR_list.append([])
                    continue
                if end_time == 0:
                    print('Folder ' + name_folder + ' is empty.')
                    ber_list.append([])
                    length_list.append(0)
                    V_mat_list.append([])
                    V_mat_rec_list.append([])
                    V_error_list.append([])
                    packetSNR_list.append([])
                    continue

                snr_values_present.append(snr_val)

                V_mat = []
                V_mat_rec = []
                V_error = []
                packetSNR = np.zeros((num_entries, 1))
                ber = np.zeros((num_entries, 1))

                index_vector = 0
                for time_idx in range(0, end_time):
                    for num_iter in range(0, end_iter):
                        name_file = (name_folder + '/results_time' + str(time_idx+1).zfill(5) + '_iteration'
                                     + str(num_iter+1).zfill(2) + '.mat')
                        # name_file = name_folder + '/results_time' + str(time_idx+1) + '_iteration' + str(num_iter+1) + '.mat'
                        try:
                            result_dict = sio.loadmat(name_file)
                            V_mat_sample = result_dict['V_mat'][0, 0]  # Nst-by-Nsts-by-Ntx
                            V_mat.append(V_mat_sample)
                            V_mat_rec_sample = result_dict['V_mat_rec'][0, 0]  # Nst-by-Nsts-by-Ntx
                            V_mat_rec.append(V_mat_rec_sample)
                            diff_V = np.mean(np.abs(V_mat_sample - V_mat_rec_sample), axis=0)
                            V_error.append(diff_V)
                            packetSNR[index_vector] = result_dict['packetSNR'][0, 0]
                            ber[index_vector] = result_dict['ber'][0, 0]
                            index_vector = index_vector + 1
                        except FileNotFoundError:
                            continue
                        except TypeError:
                            continue
                        except:
                            continue

                ber_list.append(ber[:index_vector])
                length_list.append(index_vector)
                V_mat_list.append(V_mat[:index_vector])
                V_mat_rec_list.append(V_mat_rec[:index_vector])
                V_error_list.append(V_error[:index_vector])
                packetSNR_list.append(packetSNR[:index_vector])

            if len(ber_list) < 1:
                break

            stats_ber = []
            for snr_idx, snr_val in enumerate(snr_values):
                if snr_val in snr_values_present:
                    stats_ber.append(cbook.boxplot_stats(ber_list[snr_idx], whis=(5, 95))[0])
                else:
                    stats_ber.append([])

            stats_V_error = []
            for antTx in range(NumTxAnts):
                for antRx in range(NumSTSVec):
                    stats_V_err = []
                    for snr_idx, snr_val in enumerate(snr_values):
                        if snr_val in snr_values_present:
                            length_l = length_list[snr_idx]
                            V_error_list_snr = V_error_list[snr_idx]
                            V_err_temp = np.zeros((length_l, 1))
                            for ii in range(length_l):
                                V_err_temp[ii] = V_error_list_snr[ii][antRx, antTx]
                            stats_V_err.append(cbook.boxplot_stats(V_err_temp, whis=(5, 95))[0])
                        else:
                            stats_V_err.append([])
                    stats_V_error.append(stats_V_err)

            with open(name_save_file, "wb") as fp:  # Pickling
                pickle.dump(stats_V_error, fp)

            name_save_file = (name_save_subdir +
                              '/ber_statistics_NumTxAnts_' + str(NumTxAnts) +
                              '_NumRxAntsVec_' + str(NumTxAnts) + '_BW_' + str(chanBW) +
                              '_idx_grouping_' + str(index_grouping) + '_idx_bf_' + str(index_bf) +
                              '_mcs_' + str(mcsVec) + '.txt')
            with open(name_save_file, "wb") as fp:  # Pickling
                pickle.dump(stats_ber, fp)

            name_save_file = (name_save_subdir +
                              '/length_statistics_NumTxAnts_' + str(NumTxAnts) +
                              '_NumRxAntsVec_' + str(NumTxAnts) + '_BW_' + str(chanBW) +
                              '_idx_grouping_' + str(index_grouping) + '_idx_bf_' + str(index_bf) +
                              '_mcs_' + str(mcsVec) + '.txt')
            with open(name_save_file, "wb") as fp:  # Pickling
                pickle.dump(np.asarray(length_list), fp)
