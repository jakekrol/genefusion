#!/usr/bin/env python3
import networkx as nx
from collections import Counter
import sys

# Read the edge list from a file
with open("FI10093.el") as f:
    edges = [tuple(map(int, line.strip().split("\t"))) for line in f]

# Count occurrences of each edge
edge_counts = Counter(edges)

# Create a weighted graph
G = nx.Graph()
for (u, v), weight in edge_counts.items():
    G.add_edge(u, v, weight=weight)

# add singletons
n = 25374
for i in range(0, n+1):
    if i not in G.nodes:
        G.add_node(i)

with open("test.adjlist2", "w") as f:
    for node, neighbors in sorted(G.adjacency(), key=lambda x: int(x[0])):
        if neighbors:
            neighbor_str = ",".join(f"({nbr}:{data['weight']})" for nbr, data in neighbors.items())
            line = f"{node},{neighbor_str}" if neighbor_str else f"{node}"
            f.write(line + "\n")
        else:
            f.write(f"{node}\n")

