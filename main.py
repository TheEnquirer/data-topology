import numpy as np

### CONSTANTS ###
epsilon = 0.75
num = 4
ndims = 2

# distance between two points
def distance(p1, p2):
    sum = 0.0
    for count, value in enumerate(p1):
        sum += (value - p2[count])**2.0
    return math.sqrt(sum)

cloud = np.random.random((num, ndims)) # each point p = [p1, p2, ..., p_ndims]

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

print("adjacency\n", adjacency)
