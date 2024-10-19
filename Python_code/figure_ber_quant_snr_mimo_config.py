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
mpl.rcParams['font.size'] = 16.5
mpl.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.cm.Accent.colors)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('BW', help='bandwidth', type=int)
    parser.add_argument('mcs', help='MCS', type=int)
    args = parser.parse_args()

    sub_folders = ['simulation']  # ['anechoic1','anechoic2','class','office']

    save_dir_base = 'Python_code/plots/BER_quant_snr_mimo_conf/'
    if not os.path.exists(save_dir_base):
        os.makedirs(save_dir_base)

    NumTxAnts_values = [2, 4]
    labels_ants= ['2×2 MIMO', '4×4 MIMO']

    snr_values = [14, 16, 18, 20, 22]  # [8, 10, 12, 14, 16, 18, 20, 22, 24]
    start_snr = 0
    end_snr = 6
    snr_values = np.asarray(snr_values[start_snr:end_snr], dtype=int)

    index_bf_values = [0, 1, 2]

    standard_name = 'AX'
    chanBW = args.BW
    mcsVec = args.mcs
    index_grouping = 0

    end_time = 2000  # 19000
    end_iter = 4
    num_entries = end_time * end_iter

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
        fig, ax = plt.subplots(1, len(NumTxAnts_values), constrained_layout=True)
        fig.set_size_inches(9.5, 2.4)  # (7, 2.4)
        # Plot boxplots from our computed statistics

        for ant_idx, NumTxAnts in enumerate(NumTxAnts_values):

            name_save_subdir = 'Python_code/results'

            if sub_fold == 'simulation':
                name_save_subdir = name_save_subdir + '_simulation/'

            name_save_subdir = name_save_subdir + '/' + standard_name + '/' + sub_fold

            for index_bf in index_bf_values:
                bar_pos = start + index_bf*space*2

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

                meanpointprops = dict(marker='D', markeredgecolor='black',
                                      markerfacecolor='firebrick')
                bp = ax[ant_idx].bxp(stats_ber, positions=np.arange(snr_values.shape[0]) + bar_pos,
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
                    ax[ant_idx].add_patch(Polygon(box_coords, facecolor=colors_name[index_bf]))

            ax[ant_idx].set_ylabel(r'BER')
            ax[ant_idx].set_yticks(np.linspace(0, 0.5, 6))
            ax[ant_idx].set_ylim([0, 0.55])
            ax[ant_idx].set_xticks(np.arange(len(snr_values)))
            ax[ant_idx].set_xticklabels(snr_values)
            ax[ant_idx].set_xlabel(r'SNR [dB]')

            ax[ant_idx].set_title(labels_ants[ant_idx])

            ax[ant_idx].legend(custom_lines, legend_name,
                       ncol=1, labelspacing=0.1, columnspacing=0.8, fontsize=16,
                       bbox_to_anchor=(0.7, 0.3, 1, 0.1), loc="lower left", handlelength=0.7)

        for axi in ax:
            axi.grid(True)
            axi.label_outer()

        name_fig = (save_dir_base + 'BER_quant_snr_' + sub_fold + '_mcs_' + str(mcsVec) +
                    '_BW_' + str(chanBW) + '_idx_grouping_' +
                    str(index_grouping) + '.pdf')

        plt.savefig(name_fig)
        plt.close()

    a = 1
