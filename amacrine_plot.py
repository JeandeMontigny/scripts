#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

#--------------------------------------------------------------------------#
def main():
    gcl_ri = [2.25, 4.07, 3.74, 3.36]
    gcl_ri_std = [0.1, 0.42, 0.44, 0.50]

    inl_ri = [3.76, 4.26, 4.43, 4.25]
    inl_ri_std = [0.1, 0.72, 0.36, 0.45]

    gcl_inl_ratio = [0.89, 0.89, 0.77, 0.73]
    gcl_inl_ratio_std = [0.1, 0.14, 0.10, 0.14]

    both_ri = [1.67, 2.33, 2.36, 2.42]
    both_ri_std = [0.1, 0.19, 0.06, 0.15]

    exclusion_factor = [0.76, 0.89, 0.93, 0.9]
    exclusion_factor_std = [0.1, 0.01, 0.02, 0.02]

    label = [5, 6, 7, 8]


    plt.figure(1)
    plt.errorbar(label, gcl_ri, gcl_ri_std)
    plt.xlim(min(label)-0.25, max(label)+0.25)
    plt.ylim(0, 6)
    plt.xlabel("Postnatal day")
    plt.ylabel("RI")
    plt.title("Amacrine cells mosaic regularity in the GCL")

    plt.figure(2)
    plt.errorbar(label, inl_ri, inl_ri_std)
    plt.xlim(min(label)-0.25, max(label)+0.25)
    plt.ylim(0, 6)
    plt.xlabel("Postnatal day")
    plt.ylabel("RI")
    plt.title("Amacrine cells mosaic regularity in the INL")

    plt.figure(3)
    plt.errorbar(label, gcl_inl_ratio, gcl_inl_ratio_std)
    plt.xlim(min(label)-0.25, max(label)+0.25)
    plt.ylim(0, 1.5)
    plt.xlabel("Postnatal day")
    plt.ylabel("GCL/INL ration")
    plt.title("GCL/INL population ration")

    plt.figure(4)
    plt.errorbar(label, both_ri, both_ri_std)
    plt.xlim(min(label)-0.25, max(label)+0.25)
    plt.ylim(0, 6)
    plt.xlabel("Postnatal day")
    plt.ylabel("RI")
    plt.title("Combined amacrine cells mosaic regularity (INL + GCL)")

    plt.figure(5)
    plt.errorbar(label, exclusion_factor, exclusion_factor_std)
    plt.xlim(min(label)-0.25, max(label)+0.25)
    plt.ylim(0, 1.5)
    plt.xlabel("Postnatal day")
    plt.ylabel("Index)")
    plt.title("INL and GCL amacrine cells mosaics exclusion factor over development")

    plt.show()

#--------------------------------------------------------------------------#
main()
