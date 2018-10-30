#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

#--------------------------------------------------------------------------#
def plot():

    # raw data
    mean_fdm = np.average([6.48, 6.66, 6.64, 6.28, 6.18])
    std_fdm = np.std([6.48, 6.66, 6.64, 6.28, 6.18])
    death_fdm = round(np.average([74.13, 74.36, 74.54, 75.06, 74.43]), 1)

    mean_fd = np.average([2.80, 3.01, 2.99,  2.89, 3.05])
    std_fd = np.std([2.80, 3.01, 2.99,  2.89, 3.05])
    death_fd = round(np.average([81.65, 81.20, 80.95, 81.40, 80.50]), 1)

    mean_d = np.average([2.73, 2.61, 2.69, 2.48, 2.55])
    std_d = np.std([2.73, 2.61, 2.69, 2.48, 2.55])
    death_d = round(np.average([83.90, 82.90, 83.45, 82.77, 83.47]), 1)

    mean_md = np.average([6.02, 6.12, 5.83, 5.99, 6.00])
    std_md = np.std([6.02, 6.12, 5.83, 5.99, 6.00])
    death_md = round(np.average([76.56, 75.20, 76.45, 74.93, 76.77]), 1)

    # tab containing mean and std for each group
    tab_mean = [mean_fdm, mean_md, mean_fd, mean_d]
    tab_std = [std_fdm, std_md, std_fd, std_d]
    tab_death = [death_fdm, death_md, death_fd, death_d]

    if not len(tab_mean) == len(tab_std) or not len(tab_mean) == len(tab_death):
        exit("error in raw data, all tab should be of same length")

    fig=plt.figure()

    N=len(tab_mean)
    ind=np.arange(N) # the x locations for the groups
    width=0.35 # the width of the bars
    bar_colour=['dimgrey', 'gray', 'silver', 'gainsboro'] # bars colour
    
    # bars
    plt.bar(ind, tab_mean, align='center', color=bar_colour, yerr=tab_std, error_kw=dict(elinewidth=2,ecolor='black'))

    # x and y legend
    plt.xticks(ind, ('F-D-M', 'M-D', 'F-D', 'D'))
    plt.ylabel('Regularity Index')

    # print % of cell death above bars
    for x, y, value in zip(ind, tab_mean, tab_death):
        plt.text(x-0.1, y+0.2, str(value)+"%")

    # draw horizontal line at random RI value
    plt.axhline(y=1.8, color='black', linestyle='--')

    plt.show()

#--------------------------------------------------------------------------#
plot()
