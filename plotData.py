"""
    Plotting variety of data
"""


from math import sqrt
from scipy.ndimage import rotate

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
    time = time[1].split(".h")
    time = round(float(time[0]), 4)
    print("Reading for %f" % time)
 
    # Processing Data 
    row_count = df.shape[0]
    densityDist = np.array(df['rho'].tolist())
    densityDist = densityDist.reshape((int(sqrt(row_count)), -1))
    densityDist = rotate(densityDist, angle=90)
    

    # Plotting
    fig, ax = plt.subplots()

    fig.set_size_inches(6.4,7.4)
    ax = sns.heatmap(densityDist, xticklabels=False, yticklabels=False, 
                     cmap='viridis',
                     cbar_kws={"orientation": "horizontal", 
                               "label": "Density"})

    plt.title("Density Distribution at time - " + str(time))
    plt.savefig("./images/densityMap.%f.png" %time)
    #plt.show()


if __name__ == "__main__":
    #time = float(input(blue("Enter time to plot: ", "bold")))
    files = glob.glob('results/*.h5')
    for i, fileName in enumerate(files):
        print(fileName)
        plotDensity(fileName)
        

