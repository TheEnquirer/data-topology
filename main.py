import numpy as np
import math
from sympy import *
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
import timeit
import tqdm
import networkx as nx
from p5 import *
start = timeit.default_timer()


### CONSTANTS ###
epsilon = 0.4 # max epsilon
num = 8 # number of datapoints
ndims = 2 # dimension of data

MAX_DIM = 4 # max dimension of simplices
# NUM_STEPS = 10 # number of different epsilon values between 0 and 1
NUM_STEPS = 10 # number of different epsilon values between 0 and 1

np.random.seed(1322)

def get_data(epsilon, num, ndims, MAX_DIM, data=None, generate_graph=False):

    # GENERATING ADJACENCY MATRIX

    # distance between two points
    def distance(p1, p2):
        sum = 0.0
        for count, value in enumerate(p1):
            sum += (value - p2[count])**2.0
        return math.sqrt(sum)

    # get actual data as np array
    if data==None:
        cloud = np.random.random((num, ndims)) # each point p = [p1, p2, ..., p_ndims]
    else:
        cloud = data

    # Adjacency Matrix
    adjacency = np.zeros((len(cloud), len(cloud)))
    i = 0
    while i < len(cloud):
        j = i + 1
        while j < len(cloud):
            if distance(cloud[i], cloud[j]) < epsilon:
                # decide whether directed connection is directed or undirected
                value = np.random.random()
                if value > 0.5:
                    adjacency[i][j] = 1
                    adjacency[j][i] = -1
                else:
                    adjacency[i][j] = -1
                    adjacency[j][i] = 1
            j += 1
        i += 1

    # GENERATING GRAPH

    if generate_graph and epsilon==0.4:
        # get graph
        graph = nx.DiGraph()
        for i in range(len(adjacency)):
            for j in range(len(adjacency[0])):
                if adjacency[i][j] == 1:
                    graph.add_edge(i, j)

        layout = nx.circular_layout(graph)
        nx.draw(graph, pos=layout, with_labels=True)
        plt.show()



    # SETUP

    # get adjacency matrix
    # adj_mat = [ [0, 1, 0, 0, -1, 1],
    #             [-1, 0, 0, -1, 1, 1],
    #             [0, 0, 0, -1, 0, 0],
    #             [0, 1, 1, 0, 0, 0],
    #             [1, -1, 0, 0, 0, 0],
    #             [-1, -1, 0, 0, 0, 0]]
    adj_mat = adjacency

    # generate list of nodes (just numbers 0 through number of nodes - 1)
    nodes = []
    for i in range(len(adj_mat)):
        nodes.append(i)

    # generate egde list and node degrees
    edges = []
    node_degrees = [0]*len(adj_mat)
    for i in nodes:
        for j in nodes:
            val = adj_mat[i][j]
            if val == -1:
                edges.append([i,j])
                node_degrees[i] += 1
                node_degrees[j] += 1


    # GENERATE HIGHER DIMENSIONAL SIMPLICES

    simplices = [nodes, edges]

    # check if subgraph is complete
    # a simplex is a "complete" subgraph, disregarding direction
    def is_complete(subgraph):
        for i in range(len(subgraph)):
            for j in range(i+1, len(subgraph)):
                if adj_mat[subgraph[i]][subgraph[j]] == 0:
                    return False
        return True

    # generate all n-simplices
    def gen_simplices(n, good_nodes):
        n_simplices = []
        for node in good_nodes:
            for simplex in simplices[n-1]:
                if node not in simplex:
                    g = [node,*simplex]
                    g.sort()
                    if is_complete(g) and g not in n_simplices:
                        n_simplices.append(g)
        return n_simplices

    # generate all simplices
    for n in range(2, MAX_DIM+1):
        good_nodes = []
        for i in nodes:
            if node_degrees[i] >= n:
                good_nodes.append(i)
        simplices.append(gen_simplices(n, good_nodes))


    # GENERATE INCIDENCE MATRICES

    # permutation calculation
    def MinSwaps(lst1, lst2):
        lst3 = [None] * len(lst1)

        for i in range(len(lst1)):
            lst3[i] = lst2.index(lst1[i])

        return calculateSwapsToSort(lst3)

    # permutation calculation
    def calculateSwapsToSort(lst):
        numOfSwaps = 0
        for index in range(len(lst)):
            if index != lst[index]:
                whereIsIndexMatchingNum = lst.index(index)
                lst[index], lst[whereIsIndexMatchingNum] = lst[whereIsIndexMatchingNum], lst[index]


                numOfSwaps +=1
            pass
        return numOfSwaps

    # get incidence between 2 simplices
    def incidence(s_large, s_small):
        # 0 if not incident
        for node in s_small:
            if node not in s_large:
                return 0
        # 1 if positive orientation (even permutation)
        # -1 if negative orientation (odd permutation)
        v = -1
        for node in s_large:
            if node not in s_small:
                v = node
        s_new = s_small.copy()
        s_new.insert(0,v)
        if MinSwaps(s_new,s_large)%2 == 0:
            return 1
        else:
            return -1

    # generate incidence matrix for n and n+1 simplices
    def n_incidences(n):
        inc_mat = []
        for i in simplices[n+1]:
            incidences = []
            for j in simplices[n]:
                incidences.append(incidence(i,j))
            inc_mat.append(incidences)
        return inc_mat

    # generate 1d incidence matrix
    def inc_mat_1d():
        inc_mat = []
        for edge in edges:
            incidences = []
            for node in nodes:
                if node == edge[1]:
                    incidences.append(1)
                elif node == edge[0]:
                    incidences.append(-1)
                else:
                    incidences.append(0)
            inc_mat.append(incidences)
        return inc_mat

    # generate all incidence matrices
    inc_mats = [inc_mat_1d()]
    for n in range(1, MAX_DIM):
        inc_mats.append(n_incidences(n))

    # ROW REDUCE

    row_reduced = []
    nullities = []
    ranks = []

    nullities.append(len(nodes)-1)
    ranks.append(0)

    for i in range(len(inc_mats)):
        m = Matrix(inc_mats[i])
        m = m.T
        rank = m.rank()
        nullity = len(m.nullspace())
        row_reduced.append(m.rref())
        nullities.append(nullity)
        ranks.append(rank)

    # CALCULATE BETTI NUMBERS
    betti_nums = []

    for i in range(MAX_DIM):
        betti_nums.append(nullities[i] - ranks[i+1])

    return [betti_nums, cloud, adjacency]

data = []
raw_points = []
adjacency = []
for e in tqdm.tqdm(range(NUM_STEPS)):
    REP = 0.3
    temp = get_data(e/NUM_STEPS, num, ndims, MAX_DIM, generate_graph=False)
    # temp = get_data(REP, num, ndims, MAX_DIM, generate_graph=False)
    data.append(temp[0])
    raw_points.append(temp[1])
    if len(adjacency) == 0:
        adjacency = temp[2]
        print(temp[2], "what?", e)

for i,v in enumerate(data):
    print(data[i])

# VISUALIZATION

def plot_data_3d(data):

    epsilons = []
    dimensions = []
    counts = []

    for i in range(len(data)):
        for j in range(len(data[i])):
            epsilons.append(i/NUM_STEPS)
            dimensions.append(j)
            counts.append(data[i][j])

    fig = plt.figure()

    # syntax for 3-D projection
    ax = plt.axes(projection ='3d')

    ax.scatter(epsilons, dimensions, counts, c=dimensions)

    # syntax for plotting
    ax.set_title('Barcode :)')
    ax.set_xlabel('Epsilon')
    ax.set_ylabel('Dimension')
    ax.set_zlabel('Count')
    plt.show()

def plot_data_2d(data):

    epsilons_yes = []
    dimensions_yes = []

    epsilons_no = []
    dimensions_no = []

    for i in range(len(data)):
        for j in range(len(data[i])):
            if data[i][j] != 0:
                epsilons_yes.append(i/NUM_STEPS)
                dimensions_yes.append(j)
            else:
                epsilons_no.append(i/NUM_STEPS)
                dimensions_no.append(j)


    plt.scatter(epsilons_no, dimensions_no, c="white", s=100, marker='s')
    plt.scatter(epsilons_yes, dimensions_yes, c="blue", s=100, marker='s')

    plt.title('Barcode :)')
    plt.xlabel('Epsilon')
    plt.ylabel('Dimension')
    # plt.show()

plot_data_2d(data)
stop = timeit.default_timer()
print('execution time: ', stop - start)

def generate_drawing(raw, adjacency):
    node_list = []
    edge_list = []
    for node in raw[0]:
        x, y = node
        node_list.append([(x*700)+100, (y*700)+100])

    for i in range(len(adjacency)):
        for j in range(len(adjacency[0])):
            if adjacency[i][j] != 0:
                edge = [node_list[i], node_list[j]]
                edge_list.append(edge)

    stroke(255)
    for i in edge_list:
        line((i[0]), (i[1]))

        fill("#666666ff")
        no_stroke()
        ellipse((i[0]), 15, 15)
        ellipse((i[1]), 15, 15)
        no_fill()
        stroke(255)

    stroke_weight(2)

    for i in node_list:
        stroke(255)
        ellipse((i[0],i[1]), 30, 30)

def setup():
        size(800, 800)

def draw():
    background("#212121ff")
    no_fill()
    no_stroke()
    generate_drawing(raw_points, adjacency)

# if __name__ == '__main__':
# run()
    # pass
