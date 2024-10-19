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
mpl.rcParams['text.usetex'] = 'true'
#mpl.rcParams['text.latex.preamble'] = [r'\usepackage{newtxmath}']
mpl.rcParams['font.size'] = 15.5
mpl.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.cm.Accent.colors)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()

    save_dir_base = 'Python_code/plots/BER_quant_snr_combined/'
    if not os.path.exists(save_dir_base):
        os.makedirs(save_dir_base)

    NumTxAnts_vals = [2, 3, 4]
    BW_vals = [40, 80]
    mcs_vals = [3, 4]
    sub_folders = ['anechoic2','class','office']  # ['anechoic1','anechoic2','class','office']
    labels_subf = ['Anechoic', 'Classroom', 'Office']

    snr_values = [14, 16, 18, 20, 22]  # [8, 10, 12, 14, 16, 18, 20, 22, 24]
    start_snr = 0
    end_snr = 6
    snr_values = np.asarray(snr_values[start_snr:end_snr], dtype=int)

    index_bf_values = [0, 1, 2]

    index_grouping_AC = [0, 2, 4]
    index_grouping_AX = [0, 4, 16]

    end_time = 2000  # 19000
    end_iter = 4
    num_entries = end_time * end_iter

    standard_values = ['AC', 'AX']
    index_grouping_standard = [index_grouping_AC, index_grouping_AX]

    colors_name = ['C4', 'C1', 'C7']
    line_styles = ['-' ,'--', ':']
    marker_style = [['<', 'v', 'o'],['>', '^', 's']]
    custom_lines = []
    legend_names = []
    for index_bf in index_bf_values:
        custom_lines.append(Line2D([0], [0], color=colors_name[index_bf], linewidth=4,
                                   linestyle=line_styles[index_bf]))
        label_name =  str(index_bf)
        if label_name == '2':
            label_name = 'inf'
        legend_names.append('CB ' + label_name)

    for mcs_idx in range(len(mcs_vals)):
        for NumTxAnts_idx in range(len(NumTxAnts_vals)):
            custom_lines.append(Line2D([0], [0], color='k', linewidth=0,
                                       marker=marker_style[mcs_idx][NumTxAnts_idx]))
            label_name =  (str(NumTxAnts_vals[NumTxAnts_idx]) + '×' + str(NumTxAnts_vals[NumTxAnts_idx]) +
                           ' MCS ' + str(mcs_vals[mcs_idx]))
            legend_names.append(label_name)

    for idx_standard, standard_name in enumerate(standard_values):
        name_save_dir = 'Python_code/results/' + standard_name
        index_grouping_values = index_grouping_standard[idx_standard]

        for index_grouping in index_grouping_values:
            for chanBW in BW_vals:
                #################################
                # BOX PLOT BER
                #################################
                fig, ax = plt.subplots(1, 3, constrained_layout=True)
                fig.set_size_inches(15, 3.2)  # (7, 2.4)

                for sub_fold_idx, sub_fold in enumerate(sub_folders):
                    sub_fold = sub_folders[sub_fold_idx]
                    name_save_subdir = name_save_dir + '/' + sub_fold

                    for NumTxAnts_idx in range(len(NumTxAnts_vals)):
                        NumTxAnts = NumTxAnts_vals[NumTxAnts_idx]
                        NumSTSVec = NumTxAnts
                        NumRxAntsVec = NumSTSVec

                        for mcs_idx in range(len(mcs_vals)):
                            mcsVec = mcs_vals[mcs_idx]

                            for index_bf in index_bf_values:

                                name_save_file = (name_save_subdir +
                                                  '/ber_statistics_NumTxAnts_' + str(NumTxAnts) +
                                                  '_NumRxAntsVec_' + str(NumTxAnts) + '_BW_' + str(chanBW) +
                                                  '_idx_grouping_' + str(index_grouping) + '_idx_bf_' + str(index_bf) +
                                                  '_mcs_' + str(mcsVec) + '.txt')

                                try:
                                    with open(name_save_file, "rb") as fp:  # Pickling
                                        stats_ber = pickle.load(fp)
                                except FileNotFoundError:
                                    print(name_save_file + ' not found -------------')
                                    continue

                                if len(stats_ber) != snr_values.shape[0]:
                                    print('dimension mismatch')
                                    break

                                name_save_file = (name_save_subdir +
                                                  '/length_statistics_NumTxAnts_' + str(NumTxAnts) +
                                                  '_NumRxAntsVec_' + str(NumTxAnts) + '_BW_' + str(chanBW) +
                                                  '_idx_grouping_' + str(index_grouping) + '_idx_bf_' + str(index_bf) +
                                                  '_mcs_' + str(mcsVec) + '.txt')
                                with open(name_save_file, "rb") as fp:  # Pickling
                                    length_stats = pickle.load(fp)

                                print(name_save_file + ' lengths ' + str(length_stats))

                                mean_ber = np.zeros_like(snr_values, dtype=float)
                                q1_ber = np.zeros_like(snr_values, dtype=float)
                                q3_ber = np.zeros_like(snr_values, dtype=float)
                                for ii in range(snr_values.shape[0]):
                                    if stats_ber[ii]:
                                        mean_ber[ii] = stats_ber[ii]['mean']
                                        q1_ber[ii] = stats_ber[ii]['q1']
                                        q3_ber[ii] = stats_ber[ii]['q3']
                                    else:
                                        mean_ber[ii] = np.NAN
                                        q1_ber[ii] = np.NAN
                                        q3_ber[ii] = np.NAN
                                q1_ber = - q1_ber + mean_ber
                                q3_ber = q3_ber - mean_ber

                                # bp = ax.errorbar(np.arange(snr_values.shape[0]), mean_verr, yerr=[q1_ber, q3_ber], capsize=5,
                                #                  marker=marker_style[NumTxAnts_idx], color=colors_name[index_bf],
                                #                  markersize=6, linewidth=1, linestyle=line_styles[mcs_idx])

                                bp = ax[sub_fold_idx].plot(np.arange(snr_values.shape[0]), mean_ber,
                                                 marker=marker_style[mcs_idx][NumTxAnts_idx], color=colors_name[index_bf],
                                                 markersize=6, linewidth=1.2, linestyle=line_styles[index_bf])

                    ax[sub_fold_idx].set_ylabel(r'BER')
                    ax[sub_fold_idx].set_yticks(np.linspace(0, 0.5, 6))
                    ax[sub_fold_idx].set_ylim([-0.01, 0.52])
                    ax[sub_fold_idx].set_xticks(np.arange(len(snr_values)))
                    ax[sub_fold_idx].set_xticklabels(snr_values)
                    ax[sub_fold_idx].set_xlabel(r'SNR [dB]')

                    ax[sub_fold_idx].set_title(labels_subf[sub_fold_idx])

                fig.legend(custom_lines, legend_names,
                           ncol=9, labelspacing=0.08, columnspacing=1, fontsize=14.8, #loc="upper center",
                           bbox_to_anchor=(0, 0.76, 1.05, 1), loc="lower center",
                           handlelength=1.5)

                for axi in ax:
                    axi.grid(True)
                    axi.label_outer()

                name_fig = (save_dir_base + 'BER_quant_snr_' + standard_name + '_BW_' +
                            str(chanBW) + '_idx_grouping_' + str(index_grouping) + '_combo.pdf')

                plt.savefig(name_fig)
                plt.close()

    a = 1
