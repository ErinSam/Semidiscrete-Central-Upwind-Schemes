"""
    Plotting variety of data
"""


from math import sqrt
from scipy.ndimage import rotate
from scipy.ndimage.filters import gaussian_filter

import glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from simple_colors import *


def plotDensity(fileName):
    """
        Function that plots the density
    """
    # Reading in files
    df = pd.read_hdf("./" + fileName, "table")

    time = fileName.split("/")
    time = time[-1].split(".h")
    time = round(float(time[0]), 6)
    print("Reading for %f" % time)

 
    # PROCESSING DATA
    row_count = df.shape[0]

    # Density 
    densityDist = np.array(df['rho'].tolist())
    densityDist = densityDist.reshape((int(sqrt(row_count)), -1))
    densityDist = rotate(densityDist, angle=90)

    # Locations
    x = np.array(df['x'].tolist())
    x = x.reshape((int(sqrt(row_count)), -1))
    x = rotate(x, angle=90)

    y = np.array(df['y'].tolist())
    y = y.reshape((int(sqrt(row_count)), -1))
    y = rotate(y, angle=90)
    

    # PLOTTING
    fig = plt.figure(figsize=(18.2, 7.9))
    #fig = plt.figure()
    fig.suptitle("Density Plots for time - " + str(time))

    # Density heatmap
    ax = fig.add_subplot(1,3, 1)
    ax = sns.heatmap(densityDist, xticklabels=False, yticklabels=False, 
                     cmap='magma', vmin=0.0,
                     cbar_kws={"orientation": "horizontal", 
                               "label": "Density"})
    ax.set_title("Distribution")
    ax.set_xlabel('x')
    ax.set_ylabel('y')

    # Density Surface
    ax = fig.add_subplot(1,3, 2, projection='3d')
    ax.plot_surface(x, y, densityDist, cmap='magma')
    ax.set_zlim(0.0,5.0)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    #cset = ax.contour(x, y, densityDist, zdir='x', offset=0.0, cmap='magma')
    #cset = ax.contour(x, y, densityDist, zdir='y', offset=1.0, cmap='magma')
    #cset = ax.contour(x, y, densityDist, zdir='z', offset=0.0, cmap='plasma')
    ax.set_title("Surface Plot")
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('Density')

    # Contour Plots
    ax = fig.add_subplot(1,3, 3)
    ax.contour(x, y, densityDist, 10, cmap='magma')
    ax.set_title("Contour Plot")
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    
    plt.tight_layout()
    plt.savefig("./images/Config11_100x100(dt1e-04)/densityMap.%f.png" %time)
    #plt.show()


if __name__ == "__main__":
    #time = float(input(blue("Enter time to plot: ", "bold")))
    files = glob.glob('./results/Config11_dt1e-04/*.h5')
    for i, fileName in enumerate(files):
        print(blue("Reading in "), fileName)
        plotDensity(fileName)

