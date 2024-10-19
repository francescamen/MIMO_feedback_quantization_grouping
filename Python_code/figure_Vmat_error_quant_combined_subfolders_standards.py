import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.io as sio
import pickle
import matplotlib.cbook as cbook
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter

mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Palatino'
mpl.rcParams['text.usetex'] = 'true'
#mpl.rcParams['text.latex.preamble'] = [r'\usepackage{newtxmath}']
mpl.rcParams['font.size'] = 16
mpl.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.cm.Accent.colors)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()

    NumTxAnts_vals = [2, 3, 4]
    labels_ants= ['2×2', '3×3', '4×4']
    antTx_plot = [0, 1, 2]
    antRx_plot = [0, 1]

    mcs_values = [3, 4]

    snr_values = [14, 16, 18, 20, 22]  # [8, 10, 12, 14, 16, 18, 20, 22, 24]
    start_snr = 0
    end_snr = 6
    snr_values = np.asarray(snr_values[start_snr:end_snr], dtype=int)

    index_bf_values = [0, 1]  # 2 has zero error

    chanBWs = [40, 80]
    standards = ['AC', 'AX']
    sub_folders = ['anechoic2','class']  # ['anechoic1','anechoic2','class','office']

    index_grouping = 0

    save_dir_base = 'Python_code/plots/V_error_bar/'
    if not os.path.exists(save_dir_base):
        os.makedirs(save_dir_base)

    end_time = 2000  # 19000
    end_iter = 4
    num_entries = end_time * end_iter

    standard_values = ['AC', 'AX']

    colors_name = ['C4', 'C1', 'C7', 'C2', 'C0', 'C6', 'C5', 'C3']
    hatch_style = ['/', '\\', '||', 'x', '///', '\\\\\\']
    line_styles = ['-' ,'--', ':']
    marker_style = [['<', '>'],['o', 's']]

    handles = []
    for antTx in antTx_plot:
        for antRx in antRx_plot:
            label_name =  ('Antenna ' + str(antTx+1) +
                           ', SS ' + str(antRx+1))
            handles.append(Patch(facecolor=colors_name[antTx*len(antRx_plot)+antRx], label=label_name,
                                 hatch=hatch_style[antTx*len(antRx_plot)+antRx], edgecolor='k'))

    labels_xaxis = []
    for index_bf in index_bf_values:
        label_name = ('CB ' + str(index_bf))
        labels_xaxis.append(label_name)

    #################################
    # BOX PLOT
    #################################
    fig, ax = plt.subplots(1, 6, constrained_layout=True)
    fig.set_size_inches(15, 2.6)  # (7, 2.4)

    barWidth = 0.8
    num_bars = len(antTx_plot) * len(antRx_plot)
    num_groups = len(index_bf_values)

    for std_idx, standard_name in enumerate(standard_values):
        name_save_dir = 'Python_code/results/' + standard_name

        sub_fold = sub_folders[std_idx]
        name_save_subdir = name_save_dir + '/' + sub_fold

        chanBW = chanBWs[std_idx]

        for NumTxAnts_idx, NumTxAnts in enumerate(NumTxAnts_vals):
            NumSTSVec = NumTxAnts
            NumRxAntsVec = NumSTSVec

            idx_plot = NumTxAnts_idx + len(NumTxAnts_vals) * std_idx

            mean_verr_max = 0
            mean_verr_min = 100
            offset_blocks = -1
            offset_blocks_list = []
            for index_bf in index_bf_values:
                offset_blocks = offset_blocks + 1
                offset_blocks_list.append(offset_blocks)

                bar_positions = (offset_blocks + np.arange(num_bars) + index_bf * num_bars)
                mean_verr = np.zeros((num_bars,), dtype=float)
                for antTx in antTx_plot:
                    if antTx >= NumTxAnts:
                        break
                    for antRx in antRx_plot:
                        if antRx >= NumRxAntsVec:
                            break
                        for mcs_idx, mcs_val in enumerate(mcs_values):
                            name_save_file = (name_save_subdir +
                                              '/V_diff_statistics_NumTxAnts_' + str(NumTxAnts) +
                                              '_NumRxAntsVec_' + str(NumTxAnts) + '_BW_' + str(chanBW) +
                                              '_idx_grouping_' + str(index_grouping) + '_idx_bf_' + str(index_bf) +
                                              '_mcs_' + str(mcs_val) + '.txt')

                            try:
                                with open(name_save_file, "rb") as fp:  # Pickling
                                    stats_V_error = pickle.load(fp)
                            except FileNotFoundError:
                                print(name_save_file + ' not found -------------')
                                continue

                            if len(stats_V_error) != NumTxAnts * NumSTSVec:
                                print('dimension mismatch')
                                break

                            name_save_file = (name_save_subdir +
                                              '/length_statistics_NumTxAnts_' + str(NumTxAnts) +
                                              '_NumRxAntsVec_' + str(NumTxAnts) + '_BW_' + str(chanBW) +
                                              '_idx_grouping_' + str(index_grouping) + '_idx_bf_' + str(index_bf) +
                                              '_mcs_' + str(mcs_val) + '.txt')
                            with open(name_save_file, "rb") as fp:  # Pickling
                                length_stats = pickle.load(fp)

                            print(name_save_file + ' lengths ' + str(length_stats))

                            for ii in range(snr_values.shape[0]):
                                # stats_V_error: 1st antenna, 1ss, 1st antenna, 2ss, 2st antenna, 1ss, 2st antenna, 2ss
                                index_mean_verr = antTx*len(antRx_plot)+antRx
                                mean_verr[index_mean_verr] = (mean_verr[index_mean_verr] +
                                                              stats_V_error[antTx * NumRxAntsVec + antRx][ii]['mean'])

                mean_verr = mean_verr / snr_values.shape[0] / len(mcs_values)
                bp = ax[idx_plot].bar(bar_positions, mean_verr, width=barWidth,
                                           color=colors_name[:num_bars], fill=True, hatch=hatch_style[:num_bars],
                                           edgecolor='black', linewidth=1)
                                           # linewidth=1.2, linestyle=line_styles[index_bf])
                mean_verr_max = max(mean_verr_max, max(mean_verr[:NumTxAnts*NumRxAntsVec]))
                mean_verr_min = min(mean_verr_min, min(mean_verr[:NumTxAnts*NumRxAntsVec]))

            # ax[NumTxAnts_idx].set_yticks(np.linspace(0, 0.5, 6))
            tolerance = 0.003

            ax[idx_plot].set_ylim([mean_verr_min-tolerance, mean_verr_max+tolerance*2])

            xticks = np.arange(0, num_groups*num_bars, num_bars) + (num_bars-1)/2 + np.asarray(offset_blocks_list)
            ax[idx_plot].set_xticks(xticks)
            ax[idx_plot].set_xticklabels(labels_xaxis)
            # ax[NumTxAnts_idx].set_xlabel(r'MIMO')

            ax[idx_plot].set_title(standard_name + ' ' + labels_ants[NumTxAnts_idx] + ' ' + str(chanBW) + ' MHz')

            ax[idx_plot].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

    fig.legend(handles=handles,
               ncol=6, labelspacing=0.08, columnspacing=1, fontsize=15.1, #loc="upper center",
               bbox_to_anchor=(-0.025, 0.72, 1.05, 1), loc="lower center",
               handlelength=1.5)

    ax[0].set_ylabel(r'reconstruction error')
    for axi in ax:
        axi.grid(True)
        # axi.label_outer()

    name_fig = (save_dir_base + 'Vmat_err_quant_idx_grouping_' + str(index_grouping) + '_combo_standard.pdf')

    plt.savefig(name_fig)
    plt.close()

    a = 1
