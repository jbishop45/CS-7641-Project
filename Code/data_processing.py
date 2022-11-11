import copy
import numpy as np
import pandas as pd
import anndata as an
import scanpy as sc

import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import PercentFormatter

class SCpipeline():
	def __init__(self, data, filtered_data):
		self.data = None
		self.filtered_data = None

	def load_data(self, filename, transpose=True):
		if transpose:
			self.data = sc.read_csv(filename, first_column_names=True).T
		else:
			self.data = sc.read_csv(filename, first_column_names=True)

	def initial_viz(self, bins=60, figdims=[18,18]):
		plt.rcParams["figure.figsize"] = figdims
		fig, axs = plt.subplots(2, 2)

		axs[0,0].hist(count_depth,bins=bins)
		axs[0,0].set_xlabel("count_depth")

		axs[0,1].hist(number_of_genes,bins=bins)
		axs[0,1].set_xlabel("unique genes")

		axs[1,0].plot(range(len(count_depth)),sorted(count_depth,reverse=True))
		axs[1,0].set_yscale("log")
		axs[1,0].set_xlabel("rank")
		axs[1,0].set_ylabel("count depth - log scale")

		axs[1,1].scatter(count_depth,number_of_genes)
		axs[1,1].set_xlabel("count depth")
		axs[1,1].set_ylabel("unique genes")

	def closeup_viz(self, mindepth=1, maxdepth=800, mingenes=1, maxgenes=2500, figdims=[18,9])
		fig, axs = plt.subplots(2, 1)
		axs[0].hist(count_depth[np.where((count_depth <= maxdepth) & (count_depth >= mindepth))], bins=min([200,int(maxdepth/20)]),cumulative=False)
		axs[1].hist(number_of_genes[np.where((number_of_genes<=maxgenes) & (number_of_genes >= mingenes))],bins=min([200,int(maxgenes/20)]),cumulative=False)

	def remove_dead_cells(self, min_gene_counts = 1, min_count_depth = 300, min_genes_per_cell = 500)
		####
		# WITHIN OUR NOTEBOOK I WOULD ADVIZE MAKINGA DEEP COPY OF THE DATA AND CALLING THIS FUNCTION ON THE COPY, THIS WILL SAVE YOU FROM HAVING TO RE-IMPORT THE CSV
		####
		
		sc.pp.filter_genes(adata_filtered, min_counts = 1)
		# print(adata_filtered.X.shape)
		sc.pp.filter_cells(adata_filtered, min_counts = 300)
		# print(adata_filtered.X.shape)
		sc.pp.filter_cells(adata_filtered, min_genes = 500)
		# print(adata_filtered.X.shape)