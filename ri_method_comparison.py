#!/usr/bin/env python3
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

#--------------------------------------------------------------------------#
def plot():

    # FDM
    fdm = [6.48, 6.66, 6.64, 6.28, 6.18]
    mean_fdm = round(np.average(fdm), 2)
    std_fdm = round(np.std(fdm), 2)
    death_fdm = round(np.average([74.13, 74.36, 74.54, 75.06, 74.43]), 1)

    # DM
    dm = [6.02, 6.12, 5.83, 5.99, 6.00]
    mean_md = round(np.average(dm), 2)
    std_md = round(np.std(dm), 2)
    death_md = round(np.average([76.56, 75.20, 76.45, 74.93, 76.77]), 1)

    # FD
    fd = [3.28, 3.24, 3.26, 3.19, 3.12]
    mean_fd = round(np.average(fd), 2)
    std_fd = round(np.std(fd), 2)
    death_fd = round(np.average([72.7, 72.5, 72.9, 72.6, 72.9]), 1)

    # D
    d = [3.17, 3.01, 2.95, 3.14, 3.12]
    mean_d = round(np.average(d), 2)
    std_d = round(np.std(d), 2)
    death_d = round(np.average([75.3, 75.2, 74.9, 74.9, 75.5]), 1)

    # 2 sample t-test
    p_values = []
    t_stat_dm_fdm, p_val_dm_fdm = stats.ttest_ind(dm, fdm, equal_var=False)
    p_values.append(p_val_dm_fdm)
    t_stat_dm_fd, p_val_dm_fd = stats.ttest_ind(dm, fd, equal_var=False)
    p_values.append(p_val_dm_fd)
    t_stat_fd_d, p_val_fd_d = stats.ttest_ind(d, fd, equal_var=False)
    p_values.append(p_val_fd_d)

    # tab containing mean and std for each group
    tab_mean = [mean_fdm, mean_md, mean_fd, mean_d]
    tab_std = [std_fdm, std_md, std_fd, std_d]
    tab_death = [death_fdm, death_md, death_fd, death_d]

    #---- figure ----#
    if not len(tab_mean) == len(tab_std) or not len(tab_mean) == len(tab_death):
        exit("error in raw data, all tab should be of same length")
    N=len(tab_mean)
    # bars x locations
    ind=np.arange(N) 
    # bars width
    width=0.35
    # bars colour
    bar_colour=['dimgrey', 'gray', 'silver', 'gainsboro']

    # bars
    plt.bar(ind, tab_mean, align='center', color=bar_colour, yerr=tab_std, error_kw=dict(elinewidth=2,ecolor='black'))

    # plot significance
    plt.ylim(0, max(tab_mean) + 0.8)
    for x, y, p in zip(ind, tab_mean, p_values):
        if p < 0.05:
            plt.plot([x, x+1], [y+0.4, y+0.4], color='black')
            if p < 0.001:
                plt.text((x+x+1)/2, y+0.45, "***", ha='center')
            elif p < 0.01:
                plt.text((x+x+1)/2, y+0.45, "**", ha='center')
            else:
                plt.text((x+x+1)/2, y+0.45, "*", ha='center')

    # x and y legend
    plt.xticks(ind, ('F-D-M', 'D-M', 'F-D', 'D'))
    plt.ylabel('Regularity Index')

    # print % of cell death above bars
#    for x, y, value in zip(ind, tab_mean, tab_death):
#        plt.text(x-0.1, y+0.2, str(value)+"%")

    # draw horizontal line at random RI value
    plt.axhline(y=1.8, color='black', linestyle='--')

    plt.show()

#--------------------------------------------------------------------------#
plot()
