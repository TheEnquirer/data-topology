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


np.random.seed(23)
class Barcode:
    def __init__(self, random_graph=True, max_epsilon=0.4, num_datapoints=8, ndims=2, MAX_DIM=4, NUM_STEPS=10):
        self.adjacency = []
        # self.raw = []
        self.betti_nums = []
        self.all_betti = []
        self.cloud = []

        self.draw_graph = random_graph
        self.max_epsilon = max_epsilon
        self.num_datapoints = num_datapoints
        self.ndims = ndims
        self.MAX_DIM = MAX_DIM
        self.NUM_STEPS = NUM_STEPS

        self.generate_cloud()
        for e in tqdm.tqdm(range(self.NUM_STEPS)):
            self.generate_adj(e/self.NUM_STEPS)
            self.betti_nums.append(self.get_data())

    def pretty_print(self):
        for i in self.betti_nums:
            print(i)

    def distance(self, p1, p2):
        sum = 0.0
        for count, value in enumerate(p1):
            sum += (value - p2[count])**2.0
        return math.sqrt(sum)

    def generate_cloud(self):
        # get actual data as np array
        self.cloud = np.random.random((self.num_datapoints, self.ndims)) # each point p = [p1, p2, ..., p_ndims]
        return self.cloud

    def generate_adj(self, e):

        # Adjacency Matrix
        self.adjacency = np.zeros((len(self.cloud), len(self.cloud)))
        i = 0
        while i < len(self.cloud):
            j = i + 1
            while j < len(self.cloud):
                if self.distance(self.cloud[i], self.cloud[j]) < e:
                    # decide whether directed connection is directed or undirected
                    value = np.random.random()
                    if value > 0.5:
                        self.adjacency[i][j] = 1
                        self.adjacency[j][i] = -1
                    else:
                        self.adjacency[i][j] = -1
                        self.adjacency[j][i] = 1
                j += 1
            i += 1
        return self.adjacency

    def get_data(self):

        adj_mat = self.adjacency
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
                        g = [node, *simplex]
                        g.sort()
                        if is_complete(g) and g not in n_simplices:
                            n_simplices.append(g)
            return n_simplices

        # generate all simplices
        for n in range(2, self.MAX_DIM+1):
            good_nodes = []
            for i in nodes:
                if node_degrees[i] >= n:
                    good_nodes.append(i)
            simplices.append(gen_simplices(n, good_nodes))


        # GENERATE INCIDENCE MATRICES

        # permutation calculation
        def min_swaps(lst1, lst2):
            lst3 = [None] * len(lst1)

            for i in range(len(lst1)):
                lst3[i] = lst2.index(lst1[i])

            return calculate_swaps_to_sort(lst3)

        # permutation calculation
        def calculate_swaps_to_sort(lst):
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
            if min_swaps(s_new,s_large)%2 == 0:
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
        for n in range(1, self.MAX_DIM):
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
        local_betti = []
        for i in range(self.MAX_DIM):
            local_betti.append(nullities[i] - ranks[i+1])

        # return {
        #         "betti_nums": self.betti_nums,
        #         # "cloud": self.cloud,
        #         # "adjacency": self.adjacency
        #         }
        return local_betti
        # return [self.betti_nums, self.cloud, self.adjacency]
    def plot_data_2d(self):
        epsilons_yes = []
        dimensions_yes = []

        epsilons_no = []
        dimensions_no = []

        for i in range(len(self.betti_nums)):
            for j in range(len(self.betti_nums[i])):
                if self.betti_nums[i][j] != 0:
                    epsilons_yes.append(i/self.NUM_STEPS)
                    dimensions_yes.append(j)
                else:
                    epsilons_no.append(i/self.NUM_STEPS)
                    dimensions_no.append(j)

        plt.scatter(epsilons_no, dimensions_no, c="white", s=100, marker='s')
        plt.scatter(epsilons_yes, dimensions_yes, c="blue", s=100, marker='s')

        plt.title('Barcode')
        plt.xlabel('Epsilon')
        plt.ylabel('Dimension')
        plt.show()

    def plot_data_3d(self):

        epsilons = []
        dimensions = []
        counts = []

        for i in range(len(self.betti_nums)):
            for j in range(len(self.betti_nums[i])):
                epsilons.append(i/self.NUM_STEPS)
                dimensions.append(j)
                counts.append(self.betti_nums[i][j])
        # syntax for 3-D projection
        ax = plt.axes(projection ='3d')

        ax.scatter(epsilons, dimensions, counts, c=dimensions)
        # syntax for plotting
        ax.set_title('Barcode')
        ax.set_xlabel('Epsilon')
        ax.set_ylabel('Dimension')
        ax.set_zlabel('Count')
        plt.show()

barcode = Barcode(max_epsilon=0.8)
barcode.plot_data_3d()







### CONSTANTS ###
# epsilon = 0.4 # max epsilon
# num = 8 # number of datapoints
# ndims = 2 # dimension of data

# MAX_DIM = 4 # max dimension of simplices
# # NUM_STEPS = 10 # number of different epsilon values between 0 and 1
# NUM_STEPS = 10 # number of different epsilon values between 0 and 1

# np.random.seed(1322)

# data = []
# raw_points = []
# adjacency = []
# for e in tqdm.tqdm(range(NUM_STEPS)):
#     REP = 0.3
#     temp = get_data(e/NUM_STEPS, num, ndims, MAX_DIM, generate_graph=False)
#     # temp = get_data(REP, num, ndims, MAX_DIM, generate_graph=False)
#     data.append(temp[0])
#     raw_points.append(temp[1])
#     if len(adjacency) == 0:
#         adjacency = temp[2]
#         print(temp[2], "what?", e)

# for i,v in enumerate(data):
#     print(data[i])

# # VISUALIZATION

# def plot_data_3d(data):

#     epsilons = []
#     dimensions = []
#     counts = []

#     for i in range(len(data)):
#         for j in range(len(data[i])):
#             epsilons.append(i/NUM_STEPS)
#             dimensions.append(j)
#             counts.append(data[i][j])

#     fig = plt.figure()

#     # syntax for 3-D projection
#     ax = plt.axes(projection ='3d')

#     ax.scatter(epsilons, dimensions, counts, c=dimensions)

#     # syntax for plotting
#     ax.set_title('Barcode :)')
#     ax.set_xlabel('Epsilon')
#     ax.set_ylabel('Dimension')
#     ax.set_zlabel('Count')
#     plt.show()

# def plot_data_2d(data):

#     epsilons_yes = []
#     dimensions_yes = []

#     epsilons_no = []
#     dimensions_no = []

#     for i in range(len(data)):
#         for j in range(len(data[i])):
#             if data[i][j] != 0:
#                 epsilons_yes.append(i/NUM_STEPS)
#                 dimensions_yes.append(j)
#             else:
#                 epsilons_no.append(i/NUM_STEPS)
#                 dimensions_no.append(j)


#     plt.scatter(epsilons_no, dimensions_no, c="white", s=100, marker='s')
#     plt.scatter(epsilons_yes, dimensions_yes, c="blue", s=100, marker='s')

#     plt.title('Barcode :)')
#     plt.xlabel('Epsilon')
#     plt.ylabel('Dimension')
#     # plt.show()

# plot_data_2d(data)
# stop = timeit.default_timer()
# print('execution time: ', stop - start)

# def generate_drawing(raw, adjacency):
#     node_list = []
#     edge_list = []
#     for node in raw[0]:
#         x, y = node
#         node_list.append([(x*700)+100, (y*700)+100])

#     for i in range(len(adjacency)):
#         for j in range(len(adjacency[0])):
#             if adjacency[i][j] != 0:
#                 edge = [node_list[i], node_list[j]]
#                 edge_list.append(edge)

#     stroke(255)
#     for i in edge_list:
#         line((i[0]), (i[1]))

#         fill("#666666ff")
#         no_stroke()
#         ellipse((i[0]), 15, 15)
#         ellipse((i[1]), 15, 15)
#         no_fill()
#         stroke(255)

#     stroke_weight(2)

#     for i in node_list:
#         stroke(255)
#         ellipse((i[0],i[1]), 30, 30)

# def setup():
#         size(800, 800)

# def draw():
#     background("#212121ff")
#     no_fill()
#     no_stroke()
#     generate_drawing(raw_points, adjacency)

# # if __name__ == '__main__':
# # run()
#     # pass

