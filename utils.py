import numpy as np
import decimal
import random
import matplotlib.pyplot as plt
from time import perf_counter

def split(matrix):
    """
    Splits a given matrix into quarters.
    Input: nxn matrix
    Output: tuple containing 4 n/2 x n/2 matrices corresponding to a, b, c, d
    """
    row, col = matrix.shape
    row2, col2 = row // 2, col // 2
    return matrix[:row2, :col2], matrix[:row2, col2:], matrix[row2:, :col2], matrix[row2:, col2:]

def strassen(x, y):
    """
    Computes matrix product by divide and conquer approach, recursively.
    Input: nxn matrices x and y
    Output: nxn matrix, product of x and y
    """

    # Base case when size of matrices is 1x1
    if len(x) == 1:
        return x * y

    # Splitting the matrices into quadrants. This will be done recursively
    # until the base case is reached.
    a, b, c, d = split(x)
    e, f, g, h = split(y)

    # Computing the 7 products, recursively (p1, p2...p7)
    p1 = strassen(a, f - h)
    p2 = strassen(a + b, h)
    p3 = strassen(c + d, e)
    p4 = strassen(d, g - e)
    p5 = strassen(a + d, e + h)
    p6 = strassen(b - d, g + h)
    p7 = strassen(a - c, e + f)

    # Computing the values of the 4 quadrants of the final matrix c
    c11 = p5 + p4 - p2 + p6
    c12 = p1 + p2
    c21 = p3 + p4
    c22 = p1 + p5 - p3 - p7

    # Combining the 4 quadrants into a single matrix by stacking horizontally and vertically.
    c = np.vstack((np.hstack((c11, c12)), np.hstack((c21, c22))))

    return c

def matmult_naive(a, b):
    zip_b = zip(*b)
    zip_b = list(zip_b)
    return [[sum(ele_a * ele_b for ele_a, ele_b in zip(row_a, col_b))
             for col_b in zip_b] for row_a in a]



# A class to represent a disjoint set
class DisjointSet:
    parent = {}

    # perform MakeSet operation
    def makeSet(self, n):
        # create `n` disjoint sets (one for each vertex)
        for i in range(n):
            self.parent[i] = i

    # Find the root of the set in which element `k` belongs
    def find(self, k):
        # if `k` is root
        if self.parent[k] == k:
            return k

        # recur for the parent until we find the root
        return self.find(self.parent[k])

    # Perform Union of two subsets
    def union(self, a, b):
        # find the root of the sets in which elements `x` and `y` belongs
        x = self.find(a)
        y = self.find(b)

        self.parent[x] = y


# Function to construct MST using Kruskal’s algorithm
def runKruskalAlgorithm(edges, n):
    # stores the edges present in MST
    MST = []

    # Initialize `DisjointSet` class.
    # Create a singleton set for each element of the universe.
    ds = DisjointSet()
    ds.makeSet(n)

    index = 0

    # sort edges by increasing weight
    edges.sort(key=lambda x: x[2])

    # MST contains exactly `V-1` edges
    while len(MST) != n - 1:

        # consider the next edge with minimum weight from the graph
        (src, dest, weight) = edges[index]
        index = index + 1

        # find the root of the sets to which two endpoints
        # vertices of the next edge belongs
        x = ds.find(src)
        y = ds.find(dest)

        # if both endpoints have different parents, they belong to
        # different connected components and can be included in MST
        if x != y:
            MST.append((src, dest, weight))
            ds.union(x, y)

    return MST


if __name__ == '__main__':
    # (u, v, w) triplet represent undirected edge from
    # vertex `u` to vertex `v` having weight `w`
    edges = [
        (0, 1, 7), (1, 2, 8), (0, 3, 5), (1, 3, 9), (1, 4, 7), (2, 4, 5),
        (3, 4, 15), (3, 5, 6), (4, 5, 8), (4, 6, 9), (5, 6, 11)
    ]

    # total number of nodes in the graph (labelled from 0 to 6)
    n = 7

    # construct graph
    gg = runKruskalAlgorithm(edges, n)

    print(gg)

    import networkx as nx
    G = nx.Graph()

    for e in gg:
        G.add_edge(e[0], e[1], weight=e[2])

    import matplotlib.pyplot as plt

    pos = nx.spring_layout(G, seed=7)  # positions for all nodes - seed for reproducibility
    nx.draw_networkx_nodes(G, pos, node_size=700)
    labels = nx.get_edge_attributes(G, 'weight')
    nx.draw_networkx_edge_labels(G, pos, edge_labels=labels)
    nx.draw_networkx_labels(G, pos, font_size=20, font_family="sans-serif")

    nx.draw(G,pos)
    plt.show()


