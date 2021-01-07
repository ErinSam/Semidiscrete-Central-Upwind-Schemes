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
    fig, (ax1, ax2, ax3) = plt.subplots(1,3)

    #fig.set_size_inches(6.4,7.4)
    #smoothDist = gaussian_filter(densityDist, 1)
    # Density heatmap
    ax1 = sns.heatmap(densityDist, xticklabels=False, yticklabels=False, 
                     cmap='magma', vmin=0.0,
                     cbar_kws={"orientation": "horizontal", 
                               "label": "Density"})
    ax1.set_title("Density Distibution")

    # Density Contour Plot
    ax2 = plt.axes(projection='3d')
    ax2.plot_surface(x, y, densityDist, cmap='magma')
    ax2.set_title("Surface Plot")

    # Density Surface 
    ax3 = plt.axes(projection='3d')
    ax3.contour3D(x, y, densityDist, 55, cmap='magma')
    ax3.set_title("Density Contour Plot")
    
    plt.title("Density Plots for time - " + str(time))
    #plt.savefig("./images/Config8_100x100(dt1e-04)/densityMap.%f.png" %time)
    plt.show()


if __name__ == "__main__":
    #time = float(input(blue("Enter time to plot: ", "bold")))
    files = glob.glob('./results/Config8_dt1e-04/*.h5')
    for i, fileName in enumerate(files):
        print(blue("Reading in "), fileName)
        plotDensity(fileName)
        break       

