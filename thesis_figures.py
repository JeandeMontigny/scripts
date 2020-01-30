#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.optimize import curve_fit
import thesis_figures_data

#--------------------------------------------------------------------------#
def main():

    dev_days, rgc_nb, rgc_nb_std, rgc_ri, rgc_ri_std, surface, surface_std, corrected_rgc_nb, corrected_rgc_std, death, death_std, sac_gcl_nb, sac_gcl_nb_std, sac_gcl_ri, sac_gcl_ri_std, sac_inl_nb, sac_inl_nb_std, sac_inl_ri, sac_inl_ri_std, sac_exclusion_factor, sac_exclusion_factor_std, real_death, real_death_std, pop_real, pop_real_std, pop_simu, pop_simu_std, high_dens_death, high_dens_death_std, high_dens_death_rate, low_dens_death, low_dens_death_std, low_dens_death_rate, randon_weight, delaunay_x, delaunay_cumuls, ri_randon_weight, ri_randon_weight_std, mig_denstiy, mig_denstiy_std, mig_dist, mig_dist_std, output_fate_ave, output_fate_250, output_fate_250_std, output_fate_20, output_fate_20_std, output_death_ave, output_death_250, output_death_250_std, output_death_20, output_death_20_std, output_all_ave, output_all_250, output_all_250_std, output_all_20, output_all_20_std, fate_pop_ave, fate_ri_ave, fate_pop_std, fate_ri_std, death_pop_ave, death_ri_ave, death_pop_std, death_ri_std, all_pop_ave, all_ri_ave, all_pop_std, all_ri_std, dendrites_on_off_diam, dendrites_on_off_aniso, dendrites_on_off_branch, dendrites_on_diam, dendrites_on_aniso, dendrites_on_branch, dendrites_off_diam, dendrites_off_aniso, dendrites_off_branch, real_dendrites_all_cells, real_dendrites_on_off_cells, real_dendrites_on_cells, real_dendrites_off_cells, centre_peri_p2, centre_peri_p3, centre_peri_p4, centre_peri_p5, centre_peri_p6, centre_peri_p7, centre_peri_p8, centre_peri_p9, wave_origin_p2, wave_origin_p3, wave_origin_p4, wave_origin_p5, wave_origin_p6_7, wave_origin_p8_9, wave_origin_p10_13 = thesis_figures_data.get_data()

    # -------- front sizes -------- #
    ticks_size = 11; label_size = 12

    # -------- Figure 2 -------- #
    # plt.figure(figsize=(10, 5.12))
    # plot = plt.subplot(1, 2, 1)
    # plot.spines['right'].set_visible(False)
    # plot.spines['top'].set_visible(False)
    # plot.text(-7, 1.08, "A", size = 16)
    # for i in range(0, len(delaunay_cumuls)):
    #     plt.plot([val*10 for val in delaunay_x[i]],\
    #         [val/max(delaunay_cumuls[i]) for val in delaunay_cumuls[i]],\
    #         label = randon_weight[i])
    # plt.xlabel("Segment length", size = label_size)
    # plt.ylabel("Cumulative probability", size = label_size)
    # plt.xticks(size = ticks_size)
    # plt.yticks(size = ticks_size)
    # plt.legend()
    #
    # plot = plt.subplot(1, 2, 2)
    # plot.spines['right'].set_visible(False)
    # plot.spines['top'].set_visible(False)
    # plot.text(0, 26.5, "B", size = 16)
    # plt.errorbar(randon_weight, ri_randon_weight, ri_randon_weight_std)
    # plt.xlabel("Random weight", size = label_size)
    # plt.ylabel("Regularity index", size = label_size)
    # plt.xticks(size = ticks_size)
    # plt.yticks(size = ticks_size)
    #
    # plt.tight_layout()
    #
    # # -------- Figure 3 -------- #
    # plt.figure(figsize=(10, 12))
    #
    # plot = plt.subplot(2, 2, 1)
    # plot.spines['right'].set_visible(False)
    # plot.spines['top'].set_visible(False)
    # plot.text(0.8, 7.5e3, "A", size = 16)
    # plt.errorbar(dev_days, rgc_nb, rgc_nb_std)
    # plt.xlim(1.5, 10.5)
    # plt.xlabel("Postnatal day", size = label_size)
    # plt.ylabel("Cell density (cells/mm$^2$)", size = label_size)
    # plt.xticks(size = ticks_size)
    # plt.yticks(size = ticks_size)
    #
    # plot = plt.subplot(2, 2, 2)
    # plot.spines['right'].set_visible(False)
    # plot.spines['top'].set_visible(False)
    # plot.text(0.8, 18, "B", size = 16)
    # plt.errorbar(dev_days, surface, surface_std)
    # plt.xlim(1.5, 10.5)
    # plt.xlabel("Postnatal day", size = label_size)
    # plt.ylabel("Retinal surface (mm$^2$)", size = label_size)
    # plt.xticks(size = ticks_size)
    # plt.yticks(size = ticks_size)
    #
    # plot = plt.subplot(2, 2, 3)
    # plot.spines['right'].set_visible(False)
    # plot.spines['top'].set_visible(False)
    # plot.text(0.8, 8.2e4, "C", size = 16)
    # plt.errorbar(dev_days, corrected_rgc_nb, corrected_rgc_std)
    # plt.xlim(1.5, 10.5)
    # plt.xlabel("Postnatal day", size = label_size)
    # plt.ylabel("Estimated cell population", size = label_size)
    # plt.xticks(size = ticks_size)
    # plt.yticks(size = ticks_size)
    #
    # plot = plt.subplot(2, 2, 4)
    # plot.spines['right'].set_visible(False)
    # plot.spines['top'].set_visible(False)
    # plot.text(1.8, 135, "D", size = 16)
    # plt.errorbar(dev_days[1:], death, death_std)
    # plt.xlim(2.5, 10.5)
    # plt.xlabel("Postnatal day", size = label_size)
    # plt.ylabel("Cumulative apoptosis (%)", size = label_size)
    # plt.xticks(size = ticks_size)
    # plt.yticks(size = ticks_size)
    #
    # plt.tight_layout()
    #
    # # -------- Figure 4 -------- #
    # plt.figure(figsize=(10, 5.12))
    # plot = plt.subplot(1, 2, 1)
    # plot.spines['right'].set_visible(False)
    # plot.spines['top'].set_visible(False)
    # plot.text(2.9, 1.56e3, "A", size = 16)
    # plt.errorbar(dev_days[2:], sac_gcl_nb, sac_gcl_nb_std, label = "GCL population")
    # plt.errorbar(dev_days[2:], sac_inl_nb, sac_inl_nb_std, label = "INL population")
    # plt.xlim(3.5, 10.5)
    # plt.xlabel("Postnatal day", size = label_size)
    # plt.ylabel("Cumulative apoptosis (%)", size = label_size)
    # plt.xticks(size = ticks_size)
    # plt.yticks(size = ticks_size)
    # plt.legend()
    #
    # plot = plt.subplot(1, 2, 2)
    # plot.spines['right'].set_visible(False)
    # plot.spines['top'].set_visible(False)
    # plot.text(2.9, 5.15, "B", size = 16)
    # plt.errorbar(dev_days[2:], sac_gcl_ri, sac_gcl_ri_std, label = "GCL population")
    # plt.errorbar(dev_days[2:], sac_inl_ri, sac_inl_ri_std, label = "INL population")
    # plt.xlim(3.5, 10.5)
    # plt.xlabel("Postnatal day", size = label_size)
    # plt.ylabel("Cumulative apoptosis (%)", size = label_size)
    # plt.xticks(size = ticks_size)
    # plt.yticks(size = ticks_size)
    # plt.legend()
    #
    # plt.tight_layout()
    #
    # # -------- Figure 5 -------- #
    # plt.figure(figsize=(5.12, 5.12))
    # plot = plt.subplot(1, 1, 1)
    # plot.spines['right'].set_visible(False)
    # plot.spines['top'].set_visible(False)
    # plt.errorbar(dev_days[2:], sac_exclusion_factor, sac_exclusion_factor_std)
    # plt.xlim(3.5, 10.5)
    # plt.xlabel("Postnatal day", size = label_size)
    # plt.ylabel("Exclusion factor", size = label_size)
    # plt.xticks(size = ticks_size)
    # plt.yticks(size = ticks_size)
    #
    # plt.tight_layout()
    #
    # # -------- Figure 6 -------- #
    # plt.figure(figsize=(5.12, 5.12))
    # plot = plt.subplot(1, 1, 1)
    # plot.spines['right'].set_visible(False)
    # plot.spines['top'].set_visible(False)
    # plt.errorbar(dev_days, rgc_ri, rgc_ri_std)
    # plt.xlim(2.5, 10.5)
    # plt.xlabel("Postnatal day", size = label_size)
    # plt.ylabel("Regularity index", size = label_size)
    # plt.xticks(size = ticks_size)
    # plt.yticks(size = ticks_size)
    #
    # plt.tight_layout()
    #
    # # -------- Figure 8 -------- #
    # plt.figure(figsize=(10, 12))
    # plot = plt.subplot(3, 2, 1)
    # plot.spines['right'].set_visible(False)
    # plot.spines['top'].set_visible(False)
    # plot.text(-18, 5.2, "A", size = 16)
    # plt.errorbar([x for x in range(0, 138)], output_fate_ave,\
    #  linewidth=4, color = "black", label = "average ri")
    # plt.errorbar([x for x in range(0, 139)], output_fate_250,\
    #  yerr=output_fate_250_std, elinewidth = 0.5, capsize = 2,\
    #  label = "final density = 250")
    # plt.errorbar([x for x in range(0, 139)], output_fate_20,\
    #  yerr=output_fate_20_std, elinewidth = 0.5, capsize = 2,\
    #  label = "final density = 20")
    # plt.axvline(12, linestyle='--', color = "gray")
    # plt.ylim(1, 5)
    # plt.xlabel("Modelling time", size = label_size)
    # plt.ylabel("Regularity index", size = label_size)
    # plt.xticks(size = ticks_size)
    # plt.yticks(size = ticks_size)
    # plt.legend()
    #
    # plot = plt.subplot(3, 2, 2)
    # plot.spines['right'].set_visible(False)
    # plot.spines['top'].set_visible(False)
    # plot.text(-60, 2.92, "B", size = 16)
    # plt.errorbar(fate_pop_ave, fate_ri_ave,\
    #  xerr = fate_pop_std, yerr = fate_ri_std,\
    #  fmt = 'o', color = "black", ecolor = "gray")
    # plt.plot(np.unique(fate_pop_ave),\
    #  np.poly1d(np.polyfit(fate_pop_ave, fate_ri_ave, 1))(np.unique(fate_pop_ave)),\
    #  color = "red")
    # plt.ylim(2, 2.89)
    # plt.xlabel("Cell density", size = label_size)
    # plt.ylabel("Regularity index", size = label_size)
    # plt.xticks(size = ticks_size)
    # plt.yticks(size = ticks_size)
    #
    # plot = plt.subplot(3, 2, 3)
    # plot.spines['right'].set_visible(False)
    # plot.spines['top'].set_visible(False)
    # plot.text(-18, 5.2, "C", size = 16)
    # plt.errorbar([x for x in range(0, 138)], output_death_ave,\
    #  linewidth=4, color = "black")
    # plt.errorbar([x for x in range(0, 139)], output_death_250,\
    #  yerr=output_death_250_std, elinewidth = 0.5, capsize = 2)
    # plt.errorbar([x for x in range(0, 139)], output_death_20,\
    #  yerr=output_death_20_std, elinewidth = 0.5, capsize = 2)
    # plt.axvline(12, linestyle='--', color = "gray")
    # plt.axvline(60, linestyle='--', color = "gray")
    # plt.ylim(1, 5)
    # plt.xlabel("Modelling time", size = label_size)
    # plt.ylabel("Regularity index", size = label_size)
    # plt.xticks(size = ticks_size)
    # plt.yticks(size = ticks_size)
    #
    # plot = plt.subplot(3, 2, 4)
    # plot.spines['right'].set_visible(False)
    # plot.spines['top'].set_visible(False)
    # plot.text(-30, 4.5, "D", size = 16)
    # plt.errorbar(death_pop_ave, death_ri_ave,\
    #  xerr = death_pop_std, yerr = death_ri_std,\
    #  fmt = 'o', color = "black", ecolor = "gray")
    # plt.plot(np.unique(death_pop_ave),\
    #  np.poly1d(np.polyfit(death_pop_ave, death_ri_ave, 1))(np.unique(death_pop_ave)),\
    #  color = "red")
    # def func(x, a, b):
    #     return a * x / (b + x)
    # popt, pcov = curve_fit(func, death_pop_ave, death_ri_ave)
    # cont = np.linspace(min(death_pop_ave), max(death_pop_ave), 100)
    # fitted_data = [func(x, *popt) for x in cont]
    # plt.plot(cont, fitted_data, color = "blue")
    # plt.axvline(65, linestyle='--', color = "lightgray")
    # plt.axhline(3.4, linestyle='--', color = "lightgray")
    # plt.ylim(1.5, 4.4)
    # plt.xlabel("Cell density", size = label_size)
    # plt.ylabel("Regularity index", size = label_size)
    # plt.xticks(size = ticks_size)
    # plt.yticks(size = ticks_size)
    #
    # plot = plt.subplot(3, 2, 5)
    # plot.spines['right'].set_visible(False)
    # plot.spines['top'].set_visible(False)
    # plot.text(-18, 5.2, "E", size = 16)
    # plt.errorbar([x for x in range(0, 138)], output_all_ave,\
    #  linewidth=4, color = "black")
    # plt.errorbar([x for x in range(0, 139)], output_all_250,\
    #  yerr=output_all_250_std, elinewidth = 0.5, capsize = 2)
    # plt.errorbar([x for x in range(0, 139)], output_all_20,\
    #  yerr=output_all_20_std, elinewidth = 0.5, capsize = 2)
    # plt.axvline(12, linestyle='--', color = "gray")
    # plt.axvline(60, linestyle='--', color = "gray")
    # plt.ylim(1, 5)
    # plt.xlabel("Modelling time", size = label_size)
    # plt.ylabel("Regularity index", size = label_size)
    # plt.xticks(size = ticks_size)
    # plt.yticks(size = ticks_size)
    #
    # plot = plt.subplot(3, 2, 6)
    # plot.spines['right'].set_visible(False)
    # plot.spines['top'].set_visible(False)
    # plot.text(-22, 7.1, "F", size = 16)
    # plt.errorbar(all_pop_ave, all_ri_ave,\
    #  xerr = all_pop_std, yerr = all_ri_std,\
    #  fmt = 'o', color = "black", ecolor = "gray")
    # plt.plot(np.unique(all_pop_ave),\
    #  np.poly1d(np.polyfit(all_pop_ave, all_ri_ave, 1))(np.unique(all_pop_ave)),\
    #  color = "red")
    # plt.ylim(1.5, 6.9)
    # plt.xlabel("Cell density", size = label_size)
    # plt.ylabel("Regularity index", size = label_size)
    # plt.xticks(size = ticks_size)
    # plt.yticks(size = ticks_size)
    #
    # plt.tight_layout()
    #
    # # -------- Figure 9 -------- #
    # plt.figure(figsize=(10, 5.12))
    # plot = plt.subplot(1, 2, 1)
    # plot.spines['right'].set_visible(False)
    # plot.spines['top'].set_visible(False)
    # plot.text(-0.3, 106.5, "A", size = 16)
    # plt.errorbar([x for x in range(1, 11)], pop_real,\
    #  yerr = pop_real_std, label = "In-vitro")
    # plt.errorbar([x for x in range(1, 11)], pop_simu,\
    #  yerr = pop_simu_std, label = "In-silico")
    # plt.xlim(0.5, 10.5)
    # plt.xlabel("Postnatal day", size = label_size)
    # plt.ylabel("Cell population (%)", size = label_size)
    # plt.xticks(size = ticks_size)
    # plt.yticks(size = ticks_size)
    # plt.legend()
    #
    # plot = plt.subplot(1, 2, 2)
    # plot.spines['right'].set_visible(False)
    # plot.spines['top'].set_visible(False)
    # plot.text(-6, 11.3, "B", size = 16)
    # plt.errorbar(high_dens_death_rate, high_dens_death, high_dens_death_std, label = "High density")
    # plt.errorbar(low_dens_death_rate, low_dens_death, low_dens_death_std, label = "Low density")
    # plt.axvline(65, color='silver', linestyle='--')
    # plt.xlabel("Death rate", size = label_size)
    # plt.ylabel("Regularity index", size = label_size)
    # plt.xticks(size = ticks_size)
    # plt.yticks(size = ticks_size)
    # plt.legend()
    #
    # plt.tight_layout()
    #
    # # -------- Figure 11 -------- #
    # plt.figure(figsize=(10, 5.12))
    # plot = plt.subplot(1, 2, 1)
    # plot.spines['right'].set_visible(False)
    # plot.spines['top'].set_visible(False)
    # plot.text(-5, 1800, "A", size = 16)
    # plt.bar([6.42, 12.85, 19.28, 25.71, 32.14, 38.57, 45.0, 51.43, 57.86, 64.29], [1705, 607, 330, 174, 95, 46, 23, 10, 3, 1], width = 6)
    # plt.xlabel("Distance (µm)", size = label_size)
    # plt.ylabel("Number of cell", size = label_size)
    # plt.xticks(size = ticks_size)
    # plt.yticks(size = ticks_size)
    #
    # plot = plt.subplot(1, 2, 2)
    # plot.spines['right'].set_visible(False)
    # plot.spines['top'].set_visible(False)
    # plot.text(-20, 26.2, "B", size = 16)
    # plt.errorbar(mig_denstiy, mig_dist, xerr = mig_denstiy_std, yerr = mig_dist_std, fmt = 'o', color = "black", ecolor = "gray")
    # plt.plot([106.94, 113.36, 125.55, 132.02, 153.03, 205.31, 249.18, 345.19], [4.11, 4.34, 4.79, 5.02, 5.79, 7.68, 9.28, 12.77], color = "red")
    # plt.xlabel("Cell type density", size = label_size)
    # plt.ylabel("Distance (µm)", size = label_size)
    # plt.xticks(size = ticks_size)
    # plt.yticks(size = ticks_size)
    #
    # plt.tight_layout()

    # -------- Figure 12 -------- #
    real_cells = [real_dendrites_all_cells, real_dendrites_on_cells, real_dendrites_off_cells, real_dendrites_on_off_cells]

    dendrites_all = [[dendrites_on_diam + dendrites_off_diam + dendrites_on_off_diam], [dendrites_on_aniso + dendrites_off_aniso + dendrites_on_off_aniso], [dendrites_on_branch + dendrites_off_branch + dendrites_on_off_branch]]
    dendrites_on = [dendrites_on_diam, dendrites_on_aniso, dendrites_on_branch]
    dendrites_off = [dendrites_off_diam, dendrites_off_aniso, dendrites_off_branch]
    dendrites_on_off = [dendrites_on_off_diam, dendrites_on_off_aniso, dendrites_on_off_branch]

    tab = [dendrites_all, dendrites_on, dendrites_off, dendrites_on_off]

    colours_labels = []; initialise_label = True
    violin_labels = ["", "all", "", "on", "", "off", "", "on-off"]
    plot_labels = ["Arbour diameter (µm)", "Anisometry score", "Branching number"]
    def add_label(violin, label):
        color = violin["bodies"][0].get_facecolor().flatten()
        colours_labels.append((mpatches.Patch(color=color), label))

    plt.figure(figsize=(6.5, 10))
    for i in range(0, len(tab[0])):
        ax = plt.subplot(3, 1, i+1)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        if i == 0:
            ax.text(0.4, 275, "A", size = 16)
        if i == 1:
            ax.text(0.4, 0.95, "B", size = 16)
        if i == 2:
            ax.text(0.4, 435, "C", size = 16)
        # real data
        violin_parts = ax.violinplot([cell_type[i] for cell_type in real_cells], showmeans = True, showextrema = False)
        # keep just left half
        for b in violin_parts['bodies']:
            m = np.mean(b.get_paths()[0].vertices[:, 0])
            b.get_paths()[0].vertices[:, 0] =\
             np.clip(b.get_paths()[0].vertices[:, 0], -np.inf, m)
            # b.set_color('royalblue')
        b = violin_parts['cmeans']
        for j in range(0, len(b.get_paths())):
            m = np.mean(b.get_paths()[j].vertices[:, 0])
            b.get_paths()[j].vertices[:, 0] =\
             np.clip(b.get_paths()[j].vertices[:, 0], -np.inf, m)
            b.set_color('blue')
        if initialise_label:
            add_label(violin_parts, "In-vitro")

        # these cells
        violin_parts = ax.violinplot([cell_type[i] for cell_type in tab], showmeans = True, showextrema = False)
        # keep just right half
        for b in violin_parts['bodies']:
            m = np.mean(b.get_paths()[0].vertices[:, 0])
            b.get_paths()[0].vertices[:, 0] = np.clip(b.get_paths()[0].vertices[:, 0], m, np.inf)
            # b.set_color('orange')
        b = violin_parts['cmeans']
        for j in range(0, len(b.get_paths())):
            m = np.mean(b.get_paths()[j].vertices[:, 0])
            b.get_paths()[j].vertices[:, 0] = np.clip(b.get_paths()[j].vertices[:, 0], m, np.inf)
            b.set_color('darkorange')
        if initialise_label:
            add_label(violin_parts, "In-silico")
        initialise_label = False

        if i==0:
            plt.legend(*zip(*colours_labels), bbox_to_anchor=(0.87, 1))
        ax.set_xticklabels(violin_labels)
        ax.set_ylabel(plot_labels[i])
        # ax.set_title(violin_titles[i])

        plt.tight_layout()

    # -------- Clusters Figure 1 -------- #
    # fig, ax1 = plt.subplots()
    # ax1.spines['top'].set_visible(False)
    # ax1.set_xlabel("Postnatal day", size = label_size)
    # ax1.set_ylabel("Cluster D1/D2 ratio", size = label_size)
    # ax1.boxplot([centre_peri_p2, centre_peri_p3, centre_peri_p4, centre_peri_p5, centre_peri_p6, centre_peri_p7, centre_peri_p8, centre_peri_p9], labels=dev_days[:len(dev_days)-1], medianprops=dict(color='black'))
    # ax1.text(2.5, 0.82, "_______", horizontalalignment='center')
    # ax1.text(2.5, 0.82, "***", horizontalalignment='center')
    # ax1.text(3.5, 0.92, "_______", horizontalalignment='center')
    # ax1.text(3.5, 0.92, "***", horizontalalignment='center')
    # ax1.text(4.5, 0.99, "_______", horizontalalignment='center')
    # ax1.text(4.5, 0.99, "***", horizontalalignment='center')
    # ax2 = ax1.twinx()
    # ax2.spines['top'].set_visible(False)
    # ax2.set_ylabel("% difference", color='red', size = label_size)
    # ax2.plot(dev_days[0:len(dev_days)-2], [33.12, 63.71, 24.59, 25.76, 3.523, 4.863, 6.848], 'r--')
    # ax2.set_ylim(0, 100)
    #
    # plt.tight_layout()
    #
    # # -------- Clusters Figure 2 -------- #
    # fig, ax = plt.subplots()
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)
    # ax.set_xlabel("Postnatal day", size = label_size)
    # ax.set_ylabel("Wave origin periphery/centre ratio", size = label_size)
    # ax.boxplot([wave_origin_p2, wave_origin_p3, wave_origin_p4, wave_origin_p5, wave_origin_p6_7, wave_origin_p8_9, wave_origin_p10_13], labels=["2", "3", "4", "5", "6-7", "8-9", "10-13"], medianprops=dict(color='black'))
    #
    # plt.tight_layout()
    #
    #
    # # -------- Pyramidal simulated vs real -------- #
    #
    # plt.figure(figsize=(6, 10))
    # plot = plt.subplot(2, 1, 1)
    # plot.spines['right'].set_visible(False)
    # plot.spines['top'].set_visible(False)
    # plot.text(-5, 1800, "A", size = 16)
    # plt.bar([0, 1], [36.61, 37.2], yerr = [13.5, 10.05], color = ['dimgrey', 'darkgray'])
    # plt.xticks([0, 1], ['Simulated', 'Real'], size = ticks_size)
    # plt.ylabel("Number of branching point", size = label_size)
    # plt.yticks(size = ticks_size)
    #
    # plot = plt.subplot(2, 1, 2)
    # plot.spines['right'].set_visible(False)
    # plot.spines['top'].set_visible(False)
    # plt.bar([0, 1], [6069.17, 5942.89], yerr = [1361.96, 1590.46], color = ['dimgrey', 'darkgray'])
    # plt.xticks([0, 1], ['Simulated', 'Real'], size = ticks_size)
    # plt.ylabel("Dendritic length (µm)", size = label_size)
    # plt.yticks(size = ticks_size)

#--------------------------------------------------------------------------#

main()
plt.show()
