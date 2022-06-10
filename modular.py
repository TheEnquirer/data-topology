import numpy as np
import math
from sympy import *
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
import timeit
import tqdm
import networkx as nx
from p5 import *
# start = timeit.default_timer()

np.random.seed(7)
class Barcode:
    def __init__(self, random_graph=True, max_epsilon=0.4, num_datapoints=8, ndims=2, MAX_DIM=4, NUM_STEPS=10, drawing=False):
        self.adjacency = []
        self.betti_nums = []
        self.all_betti = []
        self.cloud = []

        self.drawing = drawing
        self.generating_manual = False
        self.mouse_was_pressed = False
        self.generating_points = []

        self.random_graph = random_graph
        self.max_epsilon = max_epsilon
        self.num_datapoints = num_datapoints
        self.ndims = ndims
        self.MAX_DIM = MAX_DIM
        self.NUM_STEPS = NUM_STEPS

        if self.random_graph:
            self.generate_random_cloud()
        else:
            self.generate_manual_cloud()


    def generate_manual_cloud(self):
        self.generating_manual = True

        if mouse_is_pressed and not self.mouse_was_pressed:
            # print("click!", mouse_x, mouse_y)
            self.cloud.append([((mouse_x-100) / 700), ((mouse_y-100) / 700)])
            self.generating_points.append([mouse_x, mouse_y])

        for i in self.generating_points:
            fill("#aac4f2")
            ellipse((i[0], i[1]), 20, 20)
            # print(i)
        # print(self.generating_points)

        if mouse_is_pressed:
            self.mouse_was_pressed = True
        else:
            self.mouse_was_pressed = False

        if len(self.cloud) == self.num_datapoints:
            self.generating_manual = False
            self.drawing = False
            self.run_gen()
            self.plot_data_3d()

    def run_gen(self):
        for e in tqdm.tqdm(range(self.NUM_STEPS)):
            self.generate_adj(e/self.NUM_STEPS)
            self.betti_nums.append(self.get_data())

        if self.random_graph == False:
            self.show_drawing()
            self.pretty_print()
            # self.plot_data_3d()

    def pretty_print(self):
        for i in self.betti_nums:
            print(i)
        print(self.cloud)

    def distance(self, p1, p2):
        sum = 0.0
        for count, value in enumerate(p1):
            sum += (value - p2[count])**2.0
        return math.sqrt(sum)

    def generate_random_cloud(self):
        # get actual data as np array
        self.cloud = np.random.random((self.num_datapoints, self.ndims)) # each point p = [p1, p2, ..., p_ndims]
        self.run_gen()
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
                        # g = [node, *simplex]
                        # g.sort()
                        # if g not in n_simplices:
                        #     if is_complete(g):
                        #         n_simplices.append(g)
                        complete = True
                        for i in range(len(simplex)):
                            if adj_mat[node][simplex[i]] == 0:
                                complete = False
                                break
                        if complete:
                            g = [node, *simplex]
                            g.sort()
                            if g not in n_simplices:
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
        # epsilons_yes = []
        # dimensions_yes = []

        # epsilons_no = []
        # dimensions_no = []

        # for i in range(len(self.betti_nums)):
        #     for j in range(len(self.betti_nums[i])):
        #         if self.betti_nums[i][j] != 0:
        #             epsilons_yes.append(i/self.NUM_STEPS)
        #             dimensions_yes.append(j)
        #         else:
        #             epsilons_no.append(i/self.NUM_STEPS)
        #             dimensions_no.append(j)

        # plt.scatter(epsilons_no, dimensions_no, c="white", s=100, marker='s')
        # plt.scatter(epsilons_yes, dimensions_yes, c="blue", s=100, marker='s')

        # plt.title('Barcode')
        # plt.xlabel('Epsilon')
        # plt.ylabel('Dimension')
        # plt.show()
        epsilons_yes = []
        dimensions_yes = []

        epsilons_no = []
        dimensions_no = []

        num_steps = len(self.betti_nums)

        for i in range(len(self.betti_nums)):
            for j in range(len(self.betti_nums[i])):
                if self.betti_nums[i][j] != 0:
                    epsilons_yes.append(i/num_steps)
                    dimensions_yes.append(j)
                else:
                    epsilons_no.append(i/num_steps)
                    dimensions_no.append(j)

        plt.scatter(epsilons_no, dimensions_no, c="white", s=12000/num_steps, marker='s')
        plt.scatter(epsilons_yes, dimensions_yes, c="k", s=12000/num_steps, marker='s')

        plt.yticks(np.arange(0, len(self.betti_nums[0]), step=1))
        plt.xticks(np.arange(0, 1.1, step=1/num_steps))

        plt.title('Barcode')
        plt.xlabel('Epsilon')
        plt.ylabel('Dimension')
        plt.show()


    def plot_data_3d(self):

        # epsilons = []
        # dimensions = []
        # counts = []

        # for i in range(len(self.betti_nums)):
        #     for j in range(len(self.betti_nums[i])):
        #         epsilons.append(i/self.NUM_STEPS)
        #         dimensions.append(j)
        #         counts.append(self.betti_nums[i][j])
        # # syntax for 3-D projection
        # ax = plt.axes(projection ='3d')

        # ax.scatter(epsilons, dimensions, counts, c=dimensions)
        # # syntax for plotting
        # ax.set_title('Barcode')
        # ax.set_xlabel('Epsilon')
        # ax.set_ylabel('Dimension')
        # ax.set_zlabel('Count')
        # plt.show()
        epsilons = []
        dimensions = []
        counts = []

        colors = ['tomato','orange','gold','yellowgreen','turquoise','deepskyblue','royalblue']
        color = []
        col1 = True

        num_steps = len(self.betti_nums)

        for i in range(len(self.betti_nums)):
            # col1 = True
            for j in range(len(self.betti_nums[i])):
                epsilons.append(i/num_steps)
                dimensions.append(j)
                counts.append(self.betti_nums[i][j])
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


    def show_drawing(self):
        self.drawing = True

    def generate_drawing(self, adj):
        node_list = []
        edge_list = []
        for node in self.cloud: # uhhh
            x, y = node
            node_list.append([(x*700)+100, (y*700)+100])

        for i in range(len(adj)):
            for j in range(len(adj[0])):
                if adj[i][j] != 0:
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



start = timeit.default_timer()
barcode = Barcode(num_datapoints=8, random_graph=False)

stop = timeit.default_timer()
print('execution time: ', stop - start)





def setup():
    size(800, 800)

ii = 0
def draw():
    global ii
    background("#212121ff")

    if barcode.drawing:
        no_fill()
        no_stroke()
        if ii < math.sqrt(2):
            ii += 0.1
        else:
            ii = 0
        barcode.generate_drawing(barcode.generate_adj(ii))

    if barcode.generating_manual:
        barcode.generate_manual_cloud()
        # for i in barcode.generating_points:
        # ellipse(i[0], i[1],


if barcode.drawing or barcode.generating_manual:
    run()
