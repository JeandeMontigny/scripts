#!/usr/bin/env python3
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

#TODO:  all bars in one fig?
#       ie. mechanisms conparison (D vs FD vs FDM, ...) + comparison between AB and non AB mechanisms (F vs non AB F, D vs non AB D, ...)

#--------------------------------------------------------------------------#
def plot():
    # ---- bdm simu ---- #
    #TODO: fate and movement alone
    # FDM - movement = 1.21; death = 1.216
    fdm = [6.07, 6.08, 7.28, 5.81, 6.03, 6.63, 6.26, 5.87, 6.58, 6.30]
    mean_fdm = round(np.average(fdm), 2)
    std_fdm = round(np.std(fdm), 2)
    death_fdm = round(np.average([72.5, 71.7, 70.1, 71.2, 70.8, 70.9, 71.5, 71.5, 71.1, 70.3]), 1)

    # DM - movement = 1.21; death = 1.216
    dm = [6.76, 5.83, 5.93, 6.27, 5.96, 6.00, 6.40, 5.92, 5.99, 5.86]
    mean_md = round(np.average(dm), 2)
    std_md = round(np.std(dm), 2)
    death_md = round(np.average([73.6, 74.5, 74.8, 74.7, 74.5, 74.2, 74.6, 73.9, 74.7, 73.7]), 1)

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

    # ---- non agent based ---- #
    fate_nonAB = [3.15, 3.08, 3.04, 2.88, 3.0, 3.04, 3.16, 2.97]
    mean_fate_nonAB = np.average(fate_nonAB)
    std_fate_nonAB = np.std(fate_nonAB)
    # 65%
    death_nonAB = [5.93, 5.14, 5.75, 5.39, 6.33, 5.72, 6.36, 6.06]
    mean_death_nonAB = np.average(death_nonAB)
    std_death_nonAB = np.std(death_nonAB)
    # exclusion = 30
    migration_nonAB = [16.38, 17.73, 14.48, 16.02, 15.56, 14.29, 13.36, 18.09]
    mean_migration_nonAB = np.average(migration_nonAB)
    std_migration_nonAB = np.std(migration_nonAB)

    # ---- non agent based stats ---- #
    # 2 sample t-test
    p_values = []
    t_stat_fd, p_val_fd = stats.ttest_ind(fate_nonAB, death_nonAB, equal_var=False)
    p_values.append(p_val_fd)
    t_stat_dm, p_val_dm = stats.ttest_ind(death_nonAB, migration_nonAB, equal_var=False)
    p_values.append(p_val_dm)
    # tab containing mean and std for each group
    tab_mean = [mean_migration_nonAB, mean_death_nonAB, mean_fate_nonAB]
    tab_std = [std_migration_nonAB, std_death_nonAB, std_fate_nonAB]

    # ---- non agent based figure ---- #
    N=len(tab_mean)
    ind=np.arange(N)
    width=0.35
    bar_colour=['dimgrey', 'gray', 'silver']

    f, (ax, ax2) = plt.subplots(2, 1, sharex=True)
    # plot the same data on both axes
    ax.bar(ind, tab_mean, align='center', color=bar_colour, yerr=tab_std, error_kw=dict(elinewidth=2,ecolor='black'))
    ax2.bar(ind, tab_mean, align='center', color=bar_colour, yerr=tab_std, error_kw=dict(elinewidth=2,ecolor='black'))
    # zoom-in / limit the view to different portions of the data
    ax.set_ylim(14, 18)  # outliers only
    ax2.set_ylim(0, 7)  # most of the data
    # hide the spines between ax and ax2
    ax.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax.xaxis.tick_top()
    ax.tick_params(labeltop=False)  # don't put tick labels at the top
    ax2.xaxis.tick_bottom()
    d = .01  # how big to make the diagonal lines in axes coordinates
    # arguments to pass to plot, just so we don't keep repeating them
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
    ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

    for x, y, p in zip(ind, tab_mean, p_values):
        if p < 0.05:
            ax.plot([x+0.05, x+0.95], [17.5, 17.5], color='black')
            star_x = (x+x+1)/2; star_y = 17.5
            if p < 0.001:
                ax.text(star_x, star_y, "***", ha='center')
            elif p < 0.01:
                ax.text(star_x, star_y, "**", ha='center')
            else:
                ax.text(star_x, star_y, "*", ha='center')

    plt.xticks(ind, ('Migration', 'Death', 'Fate'))
    ax.set_ylabel('Regularity Index')
    ax.set_title('Non agent based mechanisms RI')

    ax2.axhline(y=1.8, color='black', linestyle='--')

    plt.show()

    return

    # ---- stats ---- #
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

    # ---- figure ---- #
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
            star_x = (x+x+1)/2; star_y = y+0.45
            if p < 0.001:
                plt.text(star_x, star_y, "***", ha='center')
            elif p < 0.01:
                plt.text(star_x, star_y, "**", ha='center')
            else:
                plt.text(star_x, star_y, "*", ha='center')

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
