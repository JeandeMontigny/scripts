#!/usr/bin/env python3
import numpy as np
import random as rand
import matplotlib.pyplot as plt

#--------------------------------------------------------------------------#
def main():
    center_radius = 2
    peri_radius = 4

    point_center = 0
    point_peri = 0

    plt.figure(figsize=(8, 8))

    for i in range(0, 85):
        a = rand.uniform(0, 1) * 2 * np.pi
        r = peri_radius * np.sqrt(rand.uniform(0, 1))

        if r < center_radius:
            point_center += 1
        else:
            point_peri += 1

        #Cartesian coordinates
        # x = r * np.cos(a)
        # y = r * np.sin(a)
        # plt.scatter(x, y, c='green')

    ratio = round(point_peri/point_center, 2)
    # print("ratio center/peri:", ratio)

    normalised_point_center = point_center / (np.pi * center_radius * center_radius)
    normalised_point_peri = point_peri / ((np.pi * peri_radius * peri_radius) - (np.pi * center_radius * center_radius))
    normalised_ratio = round(normalised_point_center/normalised_point_peri, 2)
    # print("normalised ratio center/peri:", normalised_ratio)
    # print(point_center, point_peri, ratio, normalised_ratio)

    # figure
    circle_center = plt.Circle((0,0), 2, color='r', fill=False)
    circle_peri = plt.Circle((0,0), 4, color='b', fill=False)
    plt.gcf().gca().add_artist(circle_center)
    plt.gcf().gca().add_artist(circle_peri)

    plt.gcf().gca().spines['right'].set_visible(False)
    plt.gcf().gca().spines['top'].set_visible(False)

    plt.xlim(-5, 5)
    plt.ylim(-5, 5)
    plt.show()
    return ratio, normalised_ratio


#--------------------------------------------------------------------------#
list_ratio = []
list_normalised_ratio = []

for i in range(0, 10):
    results = main()
    list_ratio.append(results[0])
    list_normalised_ratio.append(results[1])

print("average ratio:", round(np.average(list_ratio), 2), "with std:", round(np.std(list_ratio), 2))
print("average normalised ratio:", round(np.average(list_normalised_ratio), 2), "with std:", round(np.std(list_normalised_ratio), 2))
