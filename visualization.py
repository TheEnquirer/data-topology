import numpy as np
import math
from sympy import *
from matplotlib import colors, pyplot as plt
from mpl_toolkits import mplot3d
import timeit
import tqdm
import networkx as nx
from p5 import *
from matplotlib import style

# style.use('dark_background')

start = timeit.default_timer()


data = [[11, 0, 0, 0, 0],
[9, 0, 0, 0, 0],
[7, 0, 0, 0, 0],
[3, 1, 0, 0, 0],
[2, 0, 0, 0, 0],
[0, 0, 0, 0, 0],
[0, 1, 0, 0, 0],
[0, 0, 0, 0, 0],
[0, 0, 0, 0, 0],
[0, 0, 0, 0, 0]]


def plot_data_2d(data):
    epsilons_yes = []
    dimensions_yes = []

    epsilons_no = []
    dimensions_no = []

    num_steps = len(data)

    for i in range(len(data)):
        for j in range(len(data[i])):
            if data[i][j] != 0:
                epsilons_yes.append(i/num_steps)
                dimensions_yes.append(j)
            else:
                epsilons_no.append(i/num_steps)
                dimensions_no.append(j)

    plt.scatter(epsilons_no, dimensions_no, c="white", s=12000/num_steps, marker='s')
    plt.scatter(epsilons_yes, dimensions_yes, c="k", s=12000/num_steps, marker='s')

    plt.yticks(np.arange(0, len(data[0]), step=1))
    plt.xticks(np.arange(0, 1.1, step=1/num_steps))

    plt.title('Barcode')
    plt.xlabel('Epsilon')
    plt.ylabel('Dimension')
    plt.show()

def plot_data_3d(data):

    epsilons = []
    dimensions = []
    counts = []

    colors = ['tomato','orange','gold','yellowgreen','turquoise','deepskyblue','royalblue']
    color = []
    col1 = True

    num_steps = len(data)

    for i in range(len(data)):
        # col1 = True
        for j in range(len(data[i])):
            epsilons.append(i/num_steps)
            dimensions.append(j)
            counts.append(data[i][j])
            # if col1:
            #     color.append('tab:blue')
            # else:
            #     color.append('deepskyblue')
            # col1 = not col1
            color.append(colors[j])
    # syntax for 3-D projection
    ax = plt.axes(projection ='3d')

    # ax.scatter(epsilons, dimensions, counts, c=dimensions)
    ax.bar3d(epsilons, dimensions, 0, 1/num_steps, 1, counts, shade=True, color=color)

    ax.set_title('Barcode')
    ax.set_xlabel('Epsilon')
    ax.set_ylabel('Dimension')
    ax.set_zlabel('Count')
    plt.show()

plot_data_2d(data)