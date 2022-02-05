#Date: 2021-10-10
#Usage: code for performing the scCODA


import seaborn as sb
import seaborn as sns
import sccoda
import scanpy as sc
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import anndata as ad
from sccoda.util import comp_ana as mod
from sccoda.util import data_visualization as viz
from sccoda.util import cell_composition_data as dat
warnings.filterwarnings("ignore")
data=sc.read_loom('/COVID/01-data/CELL/recluster_file/pbmc.loom', sparse=True)

l=[]
for s in set(data.obs["orig.ident"]):
    ad = data[data.obs["orig.ident"] == s]
    ad.uns["covariats"] = {"statas":ad.obs["statas"][0]}
    l.append(ad)
comp_data = dat.from_scanpy_list(l, "annotation", "covariats")

d = {"group":[], "cell":[], "count":[], "prop":[]}
n,m = comp_data.X.shape
for i in range(n):
    total = np.sum(comp_data.X[i,:])
    for j in range(m):
        d["group"].append(comp_data.obs.statas[i])
        d["cell"].append(comp_data.var.index[j])
        d["count"].append(comp_data.X[i,j])
        d["prop"].append(comp_data.X[i,j] / total)
df = pd.DataFrame.from_dict(d)
df.head()

g = sns.FacetGrid(df, col="cell", sharey=False, col_wrap=3, height=5)
g.map(sns.boxplot, "group", "prop", palette="Blues",order=["normal","mild","moderate","severe"]);
g = g.map(sns.swarmplot, "group", "prop", color="black", edgecolor="w",size=1).set_titles("{col_name}");
g.axes[0].set_ylabel("Proportion")
plt.savefig('/COVID/ann/Hashimoto_proportion1.pdf')
#plt.savefig('/COVID/01-data/CELL/recluster_file/loom/Hashimoto_proportion1.pdf')

data_salm = comp_data[comp_data.obs["statas"].isin(["normal", "severe"])]
model_salm = mod.CompositionalAnalysis(data_salm, formula="statas", reference_cell_type="automatic")
sim_results = model_salm.sample_hmc()
sim_results.summary_extended()
print("NS")
print(sim_results.credible_effects())
sim_results.set_fdr(est_fdr=0.3)
print("est_fdr=0.3")
print(sim_results.credible_effects())
import pickle as pkl
path = "/share/pub/qiuf/COVID/01-data/CELL/recluster_file/loom/NS"
sim_results.save(path)


data_salm = comp_data[comp_data.obs["statas"].isin(["mild", "severe"])]
model_salm = mod.CompositionalAnalysis(data_salm, formula="statas", reference_cell_type="automatic")
sim_results = model_salm.sample_hmc()
sim_results.summary_extended()
print("MildS")
print(sim_results.credible_effects())
sim_results.set_fdr(est_fdr=0.3)
print("est_fdr=0.3")
print(sim_results.credible_effects())
import pickle as pkl
path = "/COVID/01-data/CELL/recluster_file/loom/MildS"
sim_results.save(path)

data_salm = comp_data[comp_data.obs["statas"].isin(["moderate", "severe"])]
model_salm = mod.CompositionalAnalysis(data_salm, formula="statas", reference_cell_type="automatic")
sim_results = model_salm.sample_hmc()
sim_results.summary_extended()
print("ModeS")
print(sim_results.credible_effects())
sim_results.set_fdr(est_fdr=0.3)
print("est_fdr=0.3"), 
print(sim_results.credible_effects())
import pickle as pkl
path = "/COVID/01-data/CELL/recluster_file/loom/ModeS"
sim_results.save(path)
