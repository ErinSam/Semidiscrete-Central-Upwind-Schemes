"""
    Plotting variety of data
"""




import pandas as pd
import numpy as np
import seaborn as sns
from math import sqrt
import matplotlib.pyplot as plt
from scipy.ndimage import rotate

from simple_colors import *


def plotDensity(time):
    """
        Function that plots the density
    """
    df = pd.read_hdf("./results/" + str(time) + ".h5", "table")
 
    row_count = df.shape[0]

    densityDist = np.array(df['rho'].tolist())
    densityDist = densityDist.reshape((int(sqrt(row_count)), -1))
    densityDist = rotate(densityDist, angle=90)
    
    fig, ax = plt.subplots()

    ax = sns.heatmap(densityDist, xticklabels=False, yticklabels=False, 
                     cmap='viridis',
                     cbar_kws={"orientation": "horizontal", 
                               "label": "Density"})

    plt.title("Density Distribution at time - " + str(time))
    # plt.savefig("./images/densityMap.%f.png" %time)
    plt.show()


if __name__ == "__main__":
    time = float(input(blue("Enter time to plot: ", "bold")))
    plotDensity(time)

