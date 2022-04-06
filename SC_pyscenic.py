
import os
import glob
import pickle
import pandas as pd
import numpy as np
from dask.diagnostics import ProgressBar
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2
from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell
import seaborn as sns
from numba import njit, float64, int64, prange


DATA_FOLDER="/home/data/SRP130923_GSE142471_small_wound/SC/pyscenic"
RESOURCES_FOLDER="/home/data/SRP130923_GSE142471_small_wound/SC/pyscenic"
DATABASE_FOLDER = "/home/data/pyscenic/ms"  
SCHEDULER="81.70.205.254:8888/"
DATABASES_GLOB = os.path.join(DATABASE_FOLDER, "mm9-*.feather")
MOTIF_ANNOTATIONS_FNAME = os.path.join(DATABASE_FOLDER, "motifs-v9-nr.mgi-m0.001-o0.0.tbl") 
MM_TFS_FNAME = os.path.join(RESOURCES_FOLDER, 'mm_tfs.txt')    
SC_EXP_FNAME = os.path.join(RESOURCES_FOLDER, "sc_EXPforpyscenic.csv")
REGULONS_FNAME = os.path.join(DATA_FOLDER, "regulons.p")
MOTIFS_FNAME = os.path.join(DATA_FOLDER, "motifs.csv")
MODULES_FNAME = os.path.join(DATA_FOLDER, "modules.p")


ex_matrix = pd.read_csv(SC_EXP_FNAME,  header=0, index_col=0)
ex_matrix.shape

ex_matrix.head()



tf_names = load_tf_names(MM_TFS_FNAME)
db_fnames = glob.glob(DATABASES_GLOB)
def name(fname):
    return os.path.splitext(os.path.basename(fname))[0]
dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]

dbs

adjancencies = grnboost2(expression_data=ex_matrix, tf_names=tf_names, verbose=True)
ADJACENCIES_FNAME = os.path.join(DATA_FOLDER, "adjacencies.tsv")
adjancencies.to_csv(ADJACENCIES_FNAME, index=False, sep='\t') 

adjancencies.head()

modules = list(modules_from_adjacencies(adjancencies, ex_matrix))

with ProgressBar():
     df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME)  

df.head()

regulons = df2regulons(df)

df.to_csv(MOTIFS_FNAME)


with open(REGULONS_FNAME, 'wb') as f:
    pickle.dump(regulons, f)

auc_mtx = aucell(ex_matrix, regulons, num_workers=1)

sns.clustermap(auc_mtx, figsize=(12,12))
sns.savefig('plot_auc_mtx.pdf')

AUC_FNAME=os.path.join(DATA_FOLDER, "auc.tsv")
regulons = [r.rename(r.name.replace('(+)',' ('+str(len(r))+'g)')) for r in regulons]
auc_mtx.to_csv(AUC_FNAME)