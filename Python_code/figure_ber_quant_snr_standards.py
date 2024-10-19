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
mpl.rcParams['font.size'] = 18
mpl.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.cm.Accent.colors)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('NumTxAnts', help='number of antennas', type=int)
    parser.add_argument('BW', help='bandwidth', type=int)
    parser.add_argument('mcs', help='MCS', type=int)
    parser.add_argument('--corrupted', help='corrupted: 0 False, 1 True', type=int, default=0)
    args = parser.parse_args()

    corrupted = args.corrupted
    save_dir_base = 'Python_code/plots/BER_quant_snr_standards/'
    if not os.path.exists(save_dir_base):
        os.makedirs(save_dir_base)

    sub_folders = ['anechoic2','class','office']  # ['anechoic1','anechoic2','class','office']
    labels_subf = ['Anechoic', 'Classroom', 'Office']

    snr_values = [14, 16, 18, 20, 22]  # [8, 10, 12, 14, 16, 18, 20, 22, 24]
    start_snr = 0
    end_snr = 6
    snr_values = np.asarray(snr_values[start_snr:end_snr], dtype=int)

    index_bf_values = [0, 1, 2]

    NumTxAnts = args.NumTxAnts
    NumSTSVec = NumTxAnts
    NumRxAntsVec = NumSTSVec
    chanBW = args.BW
    mcsVec = args.mcs
    index_grouping_AC = 0
    index_grouping_AX = 0

    end_time = 2000  # 19000
    end_iter = 4
    num_entries = end_time * end_iter

    standard_values = ['AC', 'AX']
    index_grouping_standard = [index_grouping_AC, index_grouping_AX]

    custom_lines = []
    legend_name = []
    num_bars = len(index_bf_values)
    colors_name = ['C4', 'C1', 'C7', 'C0', 'C2', 'C3', 'C5', 'C6']
    space = 0.14
    start = - space * (num_bars - 1)

    for index_bf in index_bf_values:
        custom_lines.append(Line2D([0], [0], color=colors_name[index_bf], linewidth=4))
        label_name = str(index_bf)
        if label_name == '2':
            label_name = 'inf'
        legend_name.append('CB ' + label_name)

    for sub_fold_idx, sub_fold in enumerate(sub_folders):

        #################################
        # BOX PLOT BER
        #################################
        fig, ax = plt.subplots(len(standard_values), 1, constrained_layout=True)
        fig.set_size_inches(5.4, 4.8)  # (7, 2.4)
        # Plot boxplots from our computed statistics

        for st_idx, standard_name in enumerate(standard_values):
            name_save_subdir = 'Python_code/results/' + standard_name + '/' + sub_fold

            index_grouping = index_grouping_standard[st_idx]

            for index_bf in index_bf_values:
                bar_pos = start + index_bf*space*2

                name_save_file = (name_save_subdir +
                                  '/ber_statistics_NumTxAnts_' + str(NumTxAnts) +
                                  '_NumRxAntsVec_' + str(NumTxAnts) + '_BW_' + str(chanBW) +
                                  '_idx_grouping_' + str(index_grouping) + '_idx_bf_' + str(index_bf) +
                                  '_mcs_' + str(mcsVec))

                if corrupted:
                    name_save_file = name_save_file + 'corrupted_feedback'

                name_save_file = name_save_file + '.txt'

                try:
                    with open(name_save_file, "rb") as fp:  # Pickling
                        stats_ber = pickle.load(fp)
                except FileNotFoundError:
                    print(name_save_file + ' not found -------------')
                    continue

                meanpointprops = dict(marker='D', markeredgecolor='black',
                                      markerfacecolor='firebrick')
                bp = ax[st_idx].bxp(stats_ber, positions=np.arange(snr_values.shape[0])+bar_pos,
                            showfliers=False, widths=space*1.8, showmeans=True, meanprops=meanpointprops,
                               manage_ticks=False)
                plt.setp(bp['boxes'], color='black', linewidth=1)
                plt.setp(bp['medians'], color='black', linewidth=1.5)
                plt.setp(bp['whiskers'], color='black')
                for box in bp['boxes']:
                    box_x = []
                    box_y = []
                    for j in range(5):
                        box_x.append(box.get_xdata()[j])
                        box_y.append(box.get_ydata()[j])
                    box_coords = np.column_stack([box_x, box_y])
                    ax[st_idx].add_patch(Polygon(box_coords, facecolor=colors_name[index_bf]))

            ax[st_idx].set_ylabel(r'BER')
            ax[st_idx].set_yticks(np.linspace(0, 0.5, 6))
            ax[st_idx].set_ylim([0, 0.55])
            ax[st_idx].set_xticks(np.arange(len(snr_values)))
            ax[st_idx].set_xticklabels(snr_values)
            ax[st_idx].set_xlabel(r'SNR [dB]')

            ax[st_idx].set_title('IEEE 802.11' + standard_name.lower())

            ax[st_idx].legend(custom_lines, legend_name,
                       ncol=1, labelspacing=0.1, columnspacing=0.8, fontsize=17.5,
                       bbox_to_anchor=(0.7, 0.3, 1, 0.1), loc="lower left", handlelength=0.7)

        for axi in ax:
            axi.grid(True)
            axi.label_outer()

        name_fig = (save_dir_base + 'BER_quant_snr_' + sub_fold + '_NumTxAnts_' +
                    str(NumTxAnts) + '_mcs_' + str(mcsVec) + '_BW_' + str(chanBW) + '_idx_grouping_' +
                    str(index_grouping))

        if corrupted:
            name_fig = name_fig + '_corrupted_feedback'
        name_fig = name_fig + '.pdf'

        plt.savefig(name_fig)
        plt.close()

    a = 1
