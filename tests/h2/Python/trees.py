from itertools import product
import numpy as np
import tskit
import pandas as pd

# Load in data
tree = tskit.load("/home/nick/Desktop/out.trees")
inds = tree.individuals()

# Filter individuals to only get the sampled ones
sample_IDs = pd.read_csv("/home/nick/Desktop/ped.csv")
sample_IDs = sample_IDs['id'].to_list()
inds = [ind for ind in inds if ind.metadata['pedigree_id'] in sample_IDs]

# Get genomes into a list of lists for relatedness measure
ind_comps = list()
for ind in inds:
    nodes = ind.nodes
    ind_comps.append(list(nodes))

# Calculate all pairwise relatednesses
i1, i2 = range(len(inds)), range(len(inds))
k = tree.genetic_relatedness(ind_comps, indexes=list(product(i1, i2)))

# Save output
np.savetxt("k_mat.csv", np.reshape(np.matrix(k), newshape=(len(inds),len(inds))), delimiter=",")
