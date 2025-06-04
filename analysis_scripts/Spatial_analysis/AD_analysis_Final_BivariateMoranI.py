import scanpy as sc
import pandas as pd
import numpy as np
import libpysal
from esda import Moran_BV
from libpysal.weights.spatial_lag import lag_spatial
from adjustText import adjust_text
from warnings import simplefilter
from gprofiler import GProfiler
import numpy as np
import pandas as pd
import scipy.stats as stats
from libpysal.weights import W
from libpysal.weights.spatial_lag import lag_spatial
from matplotlib import colors
from scipy import sparse
import matplotlib
# Save PDF text as text
matplotlib.rcParams['pdf.fonttype'] = 42

##Modify the Moran_BV
def _transform(w, transformation):
    """Helper to transform W or Graph"""
    # if isinstance(w, W):
    #     w.transform = transformation
    #     return w
    # else:
    #     return w.transform(transformation)
    w.transform = transformation
    return w

def _slag(w, y):
    """Helper to compute lag either for W or for Graph"""
    return lag_spatial(w, y)

class Moran_BV_modified:  # noqa: N801
    """
    Bivariate Moran's I

    Parameters
    ----------
    x : array
        x-axis variable
    y : array
        wy will be on y axis
    w : W | Graph
        spatial weights instance as W or Graph aligned with x and y
    transformation  : {'R', 'B', 'D', 'U', 'V'}
                      weights transformation, default is row-standardized "r".
                      Other options include
                      "B": binary,
                      "D": doubly-standardized,
                      "U": untransformed (general weights),
                      "V": variance-stabilizing.
    permutations    : int
                      number of random permutations for calculation of pseudo
                      p_values

    Attributes
    ----------
    zx            : array
                    original x variable standardized by mean and std
    zy            : array
                    original y variable standardized by mean and std
    w             : W | Graph
                    original w object
    permutation   : int
                    number of permutations
    I             : float
                    value of bivariate Moran's I
    sim           : array
                    (if permutations>0)
                    vector of I values for permuted samples
    p_sim         : float
                    (if permutations>0)
                    p-value based on permutations (one-sided)
                    null: spatial randomness
                    alternative: the observed I is extreme
                    it is either extremely high or extremely low
    EI_sim        : array
                    (if permutations>0)
                    average value of I from permutations
    VI_sim        : array
                    (if permutations>0)
                    variance of I from permutations
    seI_sim       : array
                    (if permutations>0)
                    standard deviation of I under permutations.
    z_sim         : array
                    (if permutations>0)
                    standardized I based on permutations
    p_z_sim       : float
                    (if permutations>0)
                    p-value based on standard normal approximation from
                    permutations

    Notes
    -----

    Inference is only based on permutations as analytical results are not too
    reliable.

    Examples
    --------
    >>> import libpysal
    >>> import numpy as np

    Set random number generator seed so we can replicate the example

    >>> np.random.seed(10)

    Open the sudden infant death dbf file and read in rates for 74 and 79
    converting each to a numpy array

    >>> f = libpysal.io.open(libpysal.examples.get_path("sids2.dbf"))
    >>> SIDR74 = np.array(f.by_col['SIDR74'])
    >>> SIDR79 = np.array(f.by_col['SIDR79'])

    Read a GAL file and construct our spatial weights object

    >>> w = libpysal.io.open(libpysal.examples.get_path("sids2.gal")).read()

    Create an instance of Moran_BV

    >>> from esda.moran import Moran_BV
    >>> mbi = Moran_BV(SIDR79,  SIDR74,  w)

    What is the bivariate Moran's I value

    >>> round(mbi.I, 3)
    0.156

    Based on 999 permutations, what is the p-value of our statistic

    >>> round(mbi.p_z_sim, 3)
    0.001


    """

    def __init__(self, x, y, w, x_mean, x_std, y_mean, y_std, transformation="r", permutations=True):
        x = np.asarray(x).flatten()
        y = np.asarray(y).flatten()
        # zy = (y - y.mean()) / y.std(ddof=1)
        # zx = (x - x.mean()) / x.std(ddof=1)
        zy = (y - y_mean) / y_std
        zx = (x - x_mean) / x_std
        self.y = y
        self.x = x
        self.zx = zx
        self.zy = zy
        n = x.shape[0]
        self.den = n - 1.0  # zx'zx = zy'zy = n-1
        w = _transform(w, transformation)
        self.w = w
        self.I = self.__calc(zy)  # noqa: E741
        if permutations:
            nrp = np.random.permutation
            sim = [self.__calc(nrp(zy)) for i in range(permutations)]
            self.sim = sim = np.array(sim)
            above = sim >= self.I
            larger = above.sum()
            if (permutations - larger) < larger:
                larger = permutations - larger
            self.p_sim = (larger + 1.0) / (permutations + 1.0)
            self.EI_sim = sim.sum() / permutations
            self.seI_sim = np.array(sim).std()
            self.VI_sim = self.seI_sim**2
            self.z_sim = (self.I - self.EI_sim) / self.seI_sim
            if self.z_sim > 0:
                self.p_z_sim = stats.norm.sf(self.z_sim)
            else:
                self.p_z_sim = stats.norm.cdf(self.z_sim)

    def __calc(self, zy):
        wzy = _slag(self.w, zy)
        self.num = (self.zx * wzy).sum()
        return self.num / self.den

    @property
    def _statistic(self):
        """More consistent hidden attribute to access ESDA statistics"""
        return self.I

    @classmethod
    def by_col(
        cls,
        df,
        x,
        y=None,
        w=None,
        inplace=False,
        pvalue="sim",
        outvals=None,
        **stat_kws,
    ):
        """
        Function to compute a Moran_BV statistic on a dataframe

        Parameters
        ----------
        df : pandas.DataFrame
            a pandas dataframe with a geometry column
        X : list of strings
            column name or list of column names to use as X values to compute
            the bivariate statistic. If no Y is provided, pairwise comparisons
            among these variates are used instead.
        Y : list of strings
            column name or list of column names to use as Y values to compute
            the bivariate statistic. if no Y is provided, pariwise comparisons
            among the X variates are used instead.
        w : W | Graph
            spatial weights instance as W or Graph aligned with the dataframe. If not
            provided, this is searched for in the dataframe's metadata
        inplace : bool
            a boolean denoting whether to operate on the dataframe inplace or to
            return a series contaning the results of the computation. If
            operating inplace, the derived columns will be named
            'column_moran_local'
        pvalue : string
            a string denoting which pvalue should be returned. Refer to the
            the Moran_BV statistic's documentation for available p-values
        outvals : list of strings
            list of arbitrary attributes to return as columns from the
            Moran_BV statistic
        **stat_kws : keyword arguments
            options to pass to the underlying statistic. For this, see the
            documentation for the Moran_BV statistic.

        Returns
        --------
        If inplace, None, and operation is conducted on dataframe
        in memory. Otherwise, returns a copy of the dataframe with
        the relevant columns attached.

        """
        return _bivariate_handler(
            df,
            x,
            y=y,
            w=w,
            inplace=inplace,
            pvalue=pvalue,
            outvals=outvals,
            swapname=cls.__name__.lower(),
            stat=cls,
            **stat_kws,
        )
# ##
if_modify_moran = True
# if_modify_moran = False
# w = libpysal.weights.KNN.from_array(spatial_coord, k=num_neighbor)
# w.transform = 'r'
# moran_bv = Moran_BV(feature1, feature1, w)
##
##load the data
adata_Xenium = sc.read_h5ad("/Users/junjietang/Desktop/AD project/Xenium_data_V1/processed_xenium.h5ad")
genelist_Xenium = adata_Xenium.var_names.tolist()
##Save genelist_Xenium
np.savetxt("/Users/junjietang/Desktop/AD project/Xenium_data_V1/genelist_Xenium.txt", genelist_Xenium, fmt="%s")
##
compartment_gene_names_subset = pd.read_csv("/Users/junjietang/Desktop/AD project/Xenium_data_V1/compartment_zscore-100k_gene_names_subset.csv").iloc[:,0].tolist()
compartment_subset = np.load("/Users/junjietang/Desktop/AD project/Xenium_data_V1/compartment_zscore-100k_subset.npy")
compartment_subset = pd.DataFrame(compartment_subset, index=adata_Xenium.obs.index, columns=compartment_gene_names_subset)
LR_Cellchat = pd.read_csv("/Users/junjietang/Desktop/AD project/Xenium_data_V1/LR_Cellchat.csv")
LR_Cellchat["if_seq"] = False
for i in range(LR_Cellchat.shape[0]):
    if LR_Cellchat.iloc[i, 0] in compartment_gene_names_subset:
        if LR_Cellchat.iloc[i, 1] in compartment_gene_names_subset:
            LR_Cellchat.iloc[i, 2] = True
LR_Cellchat = LR_Cellchat[LR_Cellchat["if_seq"] == True]
LR_geneuse = np.unique(np.hstack([LR_Cellchat['ligand'],LR_Cellchat['receptor']]))
##
num_neighbor = len(adata_Xenium.obsm['adjacency_list'][0])
##Nieghborcellindex
neighborcellindex = []
neighborcellindex_max = []
for cellindex in range(adata_Xenium.shape[0]):
    # neighborcellindex.append(np.array(adata_Xenium.obsm['adjacency_list'][cellindex,:]).tolist())
    neighborcellindex.append(np.where(adata_Xenium.obsp['adjacency_matrix'][cellindex, :].toarray() >0)[1].tolist())
    if len(neighborcellindex[cellindex]) == 0:
        neighborcellindex[cellindex] = np.nan
    else:
        neighborcellindex_max.append(np.max(neighborcellindex[cellindex]))
##Loading matrix
loading_mat = adata_Xenium.uns["M"]['"AD-2"']
loading_mat = pd.DataFrame(loading_mat, index=adata_Xenium.var_names, columns=["m_" + str(i) for i in range(loading_mat.shape[1])])
##
group1 = adata_Xenium[adata_Xenium.obs['condition'] == 'AD']
group2 = adata_Xenium[adata_Xenium.obs['condition'] == 'CT']

expr1 = group1.X.mean(axis=0).A1 if hasattr(group1.X, 'A1') else group1.X.mean(axis=0)
expr2 = group2.X.mean(axis=0).A1 if hasattr(group2.X, 'A1') else group2.X.mean(axis=0)

logfc = np.log2((expr1 + 1e-9) / (expr2 + 1e-9))
logfc = np.array(logfc).flatten()
# gene names
genes = adata_Xenium.var_names

logfc_df = pd.DataFrame({'gene': genes, 'log2FC': logfc})
logfc_df = logfc_df.sort_values('log2FC', ascending=False)
logfc_df.index = logfc_df['gene']
##
# if_choose_LR = True
# if_choose_LR = False
##
topgene_dict = {}
# tailgene_dict = {}
# if if_choose_LR:
#     num_topgene = 100
# else:
#     num_topgene = 10
topgene_dict_GP = {}
# num_topgene = 10
num_topgene = 30
num_topgene_LR = 100
# num_topgene_LR = 300
# num_tailgene = num_topgene
for m_index in range(loading_mat.shape[1]):
    loading_cur = loading_mat.iloc[:, m_index]
    topgene_cur = loading_cur.sort_values(ascending=False).index[:num_topgene_LR]
    topgene_dict["m_" + str(m_index)] = topgene_cur
    # tailgene_dict["m_" + str(m_index)] = loading_cur.sort_values(ascending=True).index[:num_tailgene]
    topgene_dict_GP["m_" + str(m_index)] = loading_cur.sort_values(ascending=False).index[:num_topgene]

topgene_metagene1_LR = topgene_dict['m_0']
topgene_metagene2_LR = topgene_dict['m_15']
topLR_list_use = []
for topgene_metagene1_cur in topgene_metagene1_LR:
    if np.sum((LR_Cellchat['ligand'] == topgene_metagene1_cur)) > 0 or np.sum((LR_Cellchat['receptor'] == topgene_metagene1_cur)) > 0:
        LR_Cellchat_cur = LR_Cellchat[
            (LR_Cellchat['ligand'] == topgene_metagene1_cur) | (LR_Cellchat['receptor'] == topgene_metagene1_cur)]
        ##remove the duplicate
        LR_Cellchat_cur = LR_Cellchat_cur.drop_duplicates(subset=['ligand', 'receptor'], keep='first')
        for i in range(LR_Cellchat_cur.shape[0]):
            topLR_list_use.append([LR_Cellchat_cur['ligand'].values[i], LR_Cellchat_cur['receptor'].values[i]])
for topgene_metagene2_cur in topgene_metagene2_LR:
    if np.sum((LR_Cellchat['ligand'] == topgene_metagene2_cur)) > 0 or np.sum((LR_Cellchat['receptor'] == topgene_metagene2_cur)) > 0:
        LR_Cellchat_cur = LR_Cellchat[
            (LR_Cellchat['ligand'] == topgene_metagene2_cur) | (LR_Cellchat['receptor'] == topgene_metagene2_cur)]
        ##remove the duplicate
        LR_Cellchat_cur = LR_Cellchat_cur.drop_duplicates(subset=['ligand', 'receptor'], keep='first')
        for i in range(LR_Cellchat_cur.shape[0]):
            topLR_list_use.append([LR_Cellchat_cur['ligand'].values[i], LR_Cellchat_cur['receptor'].values[i]])
topLR_list_df = pd.DataFrame(topLR_list_use, columns=['ligand', 'receptor'])
topLR_list_df = topLR_list_df.drop_duplicates(subset=['ligand', 'receptor'], keep='first')
##Order the rows of topLR_list_df by 'ligand'
topLR_list_df = topLR_list_df.sort_values(by='ligand')
###############################################
##Target metagene
##
if_gene_gene = True
# if_gene_gene = False
if_Hic_Hic = True
# if_Hic_Hic = False
# if_Hic_gene = True
if_Hic_gene = False
##
if_gene_Hic = True
if_gene_Hic_samegene = True
#
metagene_pair_list = [["m_0", "m_15"],
                      ["m_19", "m_15"]]

for metagene_pair_index in range(len(metagene_pair_list)):
    metagene1 = metagene_pair_list[metagene_pair_index][0]
    metagene2 = metagene_pair_list[metagene_pair_index][1]
    ##
    ##Separate metagene1 and metagene2 by "_"
    metagene1_index = int(metagene1.split("_")[1])
    metagene2_index = int(metagene2.split("_")[1])
    #
    topgene_metagene1 = topgene_dict[metagene1]
    topgene_metagene2 = topgene_dict[metagene2]
    ##Intersect with compartment_gene_names_subset
    topgene_metagene1 = list(set(topgene_metagene1).intersection(set(compartment_gene_names_subset)))
    topgene_metagene2 = list(set(topgene_metagene2).intersection(set(compartment_gene_names_subset)))
    ##
    topgene_metagene1_GP = topgene_dict_GP[metagene1]
    topgene_metagene2_GP = topgene_dict_GP[metagene2]
    ##Intersect with compartment_gene_names_subset
    topgene_metagene1_GP = list(set(topgene_metagene1_GP).intersection(set(compartment_gene_names_subset)))
    topgene_metagene2_GP = list(set(topgene_metagene2_GP).intersection(set(compartment_gene_names_subset)))
    ##
    topgene_metagene1_LR = list(set(topgene_metagene1).intersection(set(LR_geneuse)))
    topgene_metagene2_LR = list(set(topgene_metagene2).intersection(set(LR_geneuse)))
    ##Extend
    Ligand_or_receptor_match_topgene_metagene1_LR = []
    for i in range(len(topgene_metagene1_LR)):
        LR_Cellchat_cur = LR_Cellchat.iloc[np.where(
            (LR_Cellchat['ligand'] == topgene_metagene1_LR[i]) | (LR_Cellchat['receptor'] == topgene_metagene1_LR[i]))[0],
                          :]
        if LR_Cellchat_cur.iloc[0, 0] == topgene_metagene1_LR[i]:
            Ligand_or_receptor_match_topgene_metagene1_LR.extend(LR_Cellchat_cur.iloc[:, 1].tolist())
        else:
            Ligand_or_receptor_match_topgene_metagene1_LR.extend(LR_Cellchat_cur.iloc[:, 0].tolist())
    Ligand_or_receptor_match_topgene_metagene1_LR = np.unique(Ligand_or_receptor_match_topgene_metagene1_LR)
    ##
    Ligand_or_receptor_match_topgene_metagene2_LR = []
    for i in range(len(topgene_metagene2_LR)):
        LR_Cellchat_cur = LR_Cellchat.iloc[np.where(
            (LR_Cellchat['ligand'] == topgene_metagene2_LR[i]) | (LR_Cellchat['receptor'] == topgene_metagene2_LR[i]))[0],
                          :]
        if LR_Cellchat_cur.iloc[0, 0] == topgene_metagene2_LR[i]:
            Ligand_or_receptor_match_topgene_metagene2_LR.extend(LR_Cellchat_cur.iloc[:, 1].tolist())
        else:
            Ligand_or_receptor_match_topgene_metagene2_LR.extend(LR_Cellchat_cur.iloc[:, 0].tolist())
    Ligand_or_receptor_match_topgene_metagene2_LR = np.unique(Ligand_or_receptor_match_topgene_metagene2_LR)
    ##
    topgene_metagene1_LR = topgene_metagene1_LR + Ligand_or_receptor_match_topgene_metagene2_LR.tolist()
    topgene_metagene2_LR = topgene_metagene2_LR + Ligand_or_receptor_match_topgene_metagene1_LR.tolist()
    ##
    ##Obtain the top ligand-rceptor pairs
    topLR_list = []
    for topgene_metagene1_cur in topgene_metagene1_LR:
        for topgene_metagene2_cur in topgene_metagene2_LR:
            if np.sum((LR_Cellchat['ligand'] == topgene_metagene1_cur) * (LR_Cellchat['receptor'] == topgene_metagene2_cur)) > 0 or np.sum((LR_Cellchat['ligand'] == topgene_metagene2_cur) * (LR_Cellchat['receptor'] == topgene_metagene1_cur)) > 0:
                if np.sum((LR_Cellchat['ligand'] == topgene_metagene1_cur) * (LR_Cellchat['receptor'] == topgene_metagene2_cur)) > 0:
                    topLR_list.append([topgene_metagene1_cur, topgene_metagene2_cur])
                else:
                    topLR_list.append([topgene_metagene2_cur, topgene_metagene1_cur])
    topLR_list = np.array(topLR_list)
    topLR_list = pd.DataFrame(topLR_list, columns=["ligand", "receptor"])
    topLR_list.to_csv("/Users/junjietang/Desktop/AD project/Xenium_data_V1/topLR_list_" + str(metagene1) + "_" + str(metagene2) + ".csv", index=False)
    topLR_list['type'] = "LR"
    ##
    topGP_list = []
    for topgene_metagene1_cur in topgene_metagene1_GP:
        for topgene_metagene2_cur in topgene_metagene2_GP:
            topGP_list.append([topgene_metagene1_cur, topgene_metagene2_cur])
    topGP_list = np.array(topGP_list)
    topGP_list = pd.DataFrame(topGP_list, columns=["gene1", "gene2"])
    topGP_list.to_csv("/Users/junjietang/Desktop/AD project/Xenium_data_V1/topGP_list_" + str(metagene1) + "_" + str(metagene2) + ".csv", index=False)
    topGP_list['type'] = "GP"
    topLR_list_copy = topLR_list
    topLR_list_copy.columns = ["gene1", "gene2","type"]
    topGP_merge = pd.concat([topGP_list, topLR_list_copy], axis=0)
    assogene = np.unique(topGP_merge['gene1'].tolist() + topGP_merge['gene2'].tolist()).tolist()
    # ##Check the overlap between the two metagenes
    # overlap = len(set(topgene_metagene1_LR).intersection(set(topgene_metagene2_LR)))
    # print("Overlap between the two metagenes: ", overlap)
    ##
    # tailgene_metagene1 = tailgene_dict[metagene1]
    # tailgene_metagene2 = tailgene_dict[metagene2]
    # ##Intersect with compartment_gene_names_subset
    # tailgene_metagene1 = list(set(tailgene_metagene1).intersection(set(compartment_gene_names_subset)))
    # tailgene_metagene2 = list(set(tailgene_metagene2).intersection(set(compartment_gene_names_subset)))

    ##############################################################################################################################
    if if_gene_gene:
        ##Gene-to-Gene spatial correlation
        global_moranI_top1top2_sample = {}
        # global_moranI_top1tail2_sample = {}
        # global_moranI_tail1top2_sample = {}
        ##
        if if_modify_moran:
            # topgene_metagene1_mean_alltissue = np.mean(adata_Xenium[:, topgene_metagene1].X, axis=0)
            # topgene_metagene2_mean_alltissue = np.mean(adata_Xenium[:, topgene_metagene2].X, axis=0)
            # topgene_metagene1_std_alltissue = np.std(adata_Xenium[:, topgene_metagene1].X.toarray(), axis=0)
            # topgene_metagene2_std_alltissue = np.std(adata_Xenium[:, topgene_metagene2].X.toarray(), axis=0)
            topgene_assogene_mean_alltissue = np.mean(adata_Xenium[:, assogene].X, axis=0)
            topgene_assogene_std_alltissue = np.std(adata_Xenium[:, assogene].X.toarray(), axis=0)
            ##
            # topgene_metagene1_mean_alltissue = pd.DataFrame(topgene_metagene1_mean_alltissue.T, index=topgene_metagene1,
            #                                                 columns=["mean"])
            # topgene_metagene2_mean_alltissue = pd.DataFrame(topgene_metagene2_mean_alltissue.T, index=topgene_metagene2,
            #                                                 columns=["mean"])
            # topgene_metagene1_std_alltissue = pd.DataFrame(topgene_metagene1_std_alltissue.T, index=topgene_metagene1,
            #                                                columns=["std"])
            # topgene_metagene2_std_alltissue = pd.DataFrame(topgene_metagene2_std_alltissue.T, index=topgene_metagene2,
            #                                                columns=["std"])
            topgene_assogene_mean_alltissue = pd.DataFrame(topgene_assogene_mean_alltissue.T, index=assogene,
                                                            columns=["mean"])
            topgene_assogene_std_alltissue = pd.DataFrame(topgene_assogene_std_alltissue.T, index=assogene,
                                                           columns=["std"])
        ##
        sample_unique = np.unique(adata_Xenium.obs['sample_id'])
        for sample_cur in sample_unique:
            print(sample_cur)
            cellindex_cursample = np.where(adata_Xenium.obs['sample_id'] == sample_cur)[0]
            cellindex_cursample = np.where(adata_Xenium.obs['sample_id'] == sample_cur)[0]
            geneexp_cursample = adata_Xenium.X[cellindex_cursample, :].toarray()
            geneexp_cursample = pd.DataFrame(geneexp_cursample,
                                             index=adata_Xenium.obs['orig.ident'].index[cellindex_cursample],
                                             columns=adata_Xenium.var_names)
            spatial_coord_cursample = adata_Xenium.obsm['spatial'][cellindex_cursample, :]
            ##
            spatial_coord_cursample_cur = spatial_coord_cursample
            ##
            w = libpysal.weights.KNN.from_array(spatial_coord_cursample_cur, k=num_neighbor)
            w.transform = 'r'
            ##Top genes
            global_moran_I_cursample = np.zeros((len(assogene), len(assogene)))
            global_moran_I_cursample = pd.DataFrame(global_moran_I_cursample,index = assogene,
                                                    columns=assogene)
            for GP_index in range(topGP_merge.shape[0]):
                topGP_merge_cur = topGP_merge.iloc[GP_index, :]
                gene1 = topGP_merge_cur['gene1']
                gene2 = topGP_merge_cur['gene2']
                ##
                geneexp_cursample_cur = geneexp_cursample.loc[:,
                                        [gene1, gene2]].values
                ##
                moran_bv = Moran_BV_modified(x=geneexp_cursample_cur[:, 0], y=geneexp_cursample_cur[:, 1],
                                             w=w,
                                             x_mean=float(
                                                 topgene_assogene_mean_alltissue.loc[
                                                     gene1]),
                                             x_std=float(
                                                 topgene_assogene_std_alltissue.loc[
                                                     gene1]),
                                             y_mean=float(
                                                 topgene_assogene_mean_alltissue.loc[
                                                     gene2]),
                                             y_std=float(
                                                 topgene_assogene_std_alltissue.loc[
                                                     gene2]))
                ##
                global_moran_I_cursample.loc[gene1, gene2] = moran_bv.I
            global_moranI_top1top2_sample[sample_cur] = global_moran_I_cursample

        ##Condition-specific
        sample_unique_AD = ["AD-1", "AD-2"]
        sample_unique_CT = ["CT-1", "CT-2"]
        global_moranI_top1top2_AD = np.zeros((len(assogene), len(assogene)))
        for sample_cur in sample_unique_AD:
            global_moranI_top1top2_AD += global_moranI_top1top2_sample[sample_cur].values
        global_moranI_top1top2_AD /= len(sample_unique_AD)
        global_moranI_top1top2_AD = pd.DataFrame(global_moranI_top1top2_AD, index=assogene,
                                                 columns=assogene)
        ##Turn the global_moranI_top1top2 into a long dataframe
        global_moranI_top1top2_AD_long = global_moranI_top1top2_AD.stack().reset_index()
        global_moranI_top1top2_AD_long.columns = ["Gene1", "Gene2", "global_moranI"]
        global_moranI_top1top2_AD_long['Group'] = "Top1Top2_AD"
        ##
        global_moranI_top1top2_AD_long = global_moranI_top1top2_AD_long.loc[~(global_moranI_top1top2_AD_long['global_moranI'] == 0),:]
        ##
        global_moranI_top1top2_CT = np.zeros((len(assogene), len(assogene)))
        for sample_cur in sample_unique_CT:
            global_moranI_top1top2_CT += global_moranI_top1top2_sample[sample_cur].values
        global_moranI_top1top2_CT /= len(sample_unique_CT)
        global_moranI_top1top2_CT = pd.DataFrame(global_moranI_top1top2_CT, index=assogene,
                                                 columns=assogene)
        ##Turn the global_moranI_top1top2 into a long dataframe
        global_moranI_top1top2_CT_long = global_moranI_top1top2_CT.stack().reset_index()
        global_moranI_top1top2_CT_long.columns = ["Gene1", "Gene2", "global_moranI"]
        global_moranI_top1top2_CT_long['Group'] = "Top1Top2_CT"
        global_moranI_top1top2_CT_long = global_moranI_top1top2_CT_long.loc[~(global_moranI_top1top2_CT_long['global_moranI'] == 0),:]
        ##
        global_moranI_all = pd.concat([global_moranI_top1top2_AD_long, global_moranI_top1top2_CT_long], axis=0)
        ##
        choose_index = []
        for i in range(global_moranI_all.shape[0]):
            gene1 = global_moranI_all.iloc[i, 0]
            gene2 = global_moranI_all.iloc[i, 1]
            if np.sum((LR_Cellchat['ligand'] == gene2) * (LR_Cellchat['receptor'] == gene1)) > 0 or np.sum(
                    (LR_Cellchat['ligand'] == gene1) * (LR_Cellchat['receptor'] == gene2)) > 0:
                choose_index.append(i)
        # global_moranI_all = global_moranI_all.iloc[choose_index, :]
        global_moranI_all['type'] = "GP"
        global_moranI_all['type'].iloc[choose_index] = "LR"

        ##Keep the global_moranI_all['Group'] == "Top1Top2_AD" or "Top1Top2_CT"
        global_moranI_all = global_moranI_all[
            (global_moranI_all['Group'] == "Top1Top2_AD") | (global_moranI_all['Group'] == "Top1Top2_CT")]
        ##change the "Top1Top2_AD" to "AD" and "Top1Top2_CT" to "CT"
        global_moranI_all['Group'] = global_moranI_all['Group'].replace({"Top1Top2_AD": "AD", "Top1Top2_CT": "CT"})
        ##save the global_moranI_all to csv
        global_moranI_all.to_csv(
            "/Users/junjietang/Desktop/AD project/Xenium_data_V1/global_moranI_all_" + str(metagene1) + "_" + str(
                metagene2) + "_genexpression.csv", index=False)
        ##Show the density plot of global_moranI_all['global_moranI'], grouped by global_moranI_all['Group']
        import seaborn as sns
        import matplotlib.pyplot as plt
        import matplotlib
        ##
        matplotlib.use('TkAgg')
        plt.rcParams.update({
            'axes.titlesize': 20,
            'axes.labelsize': 18,
            'xtick.labelsize': 16,
            'ytick.labelsize': 16,
            'legend.title_fontsize': 16,
            'legend.fontsize': 16
        })

        plt.figure(figsize=(7, 8))
        ax = sns.kdeplot(
            data=global_moranI_all,
            x='global_moranI',
            hue='Group',
            fill=True,
            common_norm=False
        )

        # ax.set_xlabel("Bi-variate Moran's I of gene expression \nbetween associated ligand-receptor pair", fontsize=18)
        ax.set_xlabel("Bi-variate Moran's I of gene expression \nbetween associated gene pairs", fontsize=18)
        ax.set_ylabel("Density", fontsize=18)
        leg = ax.get_legend()
        leg.set_title("Group")
        plt.setp(leg.get_title(), fontsize=22)
        plt.setp(leg.get_texts(), fontsize=22)

        plt.tight_layout()

        # if if_choose_LR:
        #     # plt.savefig("/Users/junjietang/Desktop/AD project/Xenium_data_V1/global_moranI_density_LR_condition_" + str(
        #     #     metagene1) + "_" + str(metagene2) + ".png", dpi=300, bbox_inches='tight')
        #     plt.savefig("/Users/junjietang/Desktop/AD project/Xenium_data_V1/global_moranI_density_LR_condition_" + str(
        #         metagene1) + "_" + str(metagene2) + ".pdf", dpi=300, bbox_inches='tight')
        # else:
        #     # plt.savefig("/Users/junjietang/Desktop/AD project/Xenium_data_V1/global_moranI_density_condition_" + str(
        #     #     metagene1) + "_" + str(metagene2) + ".png", dpi=300, bbox_inches='tight')
        #     plt.savefig("/Users/junjietang/Desktop/AD project/Xenium_data_V1/global_moranI_density_condition_" + str(
        #         metagene1) + "_" + str(metagene2) + ".pdf", dpi=300, bbox_inches='tight')
        plt.savefig("/Users/junjietang/Desktop/AD project/Xenium_data_V1/global_moranI_density_condition_" + str(
            metagene1) + "_" + str(metagene2) + ".pdf", dpi=300, bbox_inches='tight')
        plt.close()
        ##
        global_moranI_all_AD = global_moranI_all[global_moranI_all['Group'] == "AD"]
        global_moranI_all_CT = global_moranI_all[global_moranI_all['Group'] == "CT"]

        ##Show a scatter plot bettwen global_moranI_all_AD['global_moranI'] and global_moranI_all_CT['global_moranI']
        ##Merge global_moranI_all_AD and global_moranI_all_CT
        global_moranI_all_AD_CT = pd.merge(global_moranI_all_AD, global_moranI_all_CT, on=['Gene1', 'Gene2'],
                                           suffixes=('_AD', '_CT'))
        global_moranI_all_AD_CT_LR = global_moranI_all_AD_CT.iloc[np.where(global_moranI_all_AD_CT['type_CT'] == "LR")[0],:]
        # global_moranI_all_AD_CT_LR_positive = global_moranI_all_AD_CT_LR.iloc[np.where(~((global_moranI_all_AD_CT_LR['global_moranI_AD'] < 0) * (global_moranI_all_AD_CT_LR['global_moranI_CT'] < 0)))[0], :]
        ##
        global_moranI_all_AD_CT_positive = global_moranI_all_AD_CT.loc[~(
                    (global_moranI_all_AD_CT['global_moranI_AD'] < 0) * (
                        global_moranI_all_AD_CT['global_moranI_CT'] < 0)), :]
        ##just highlight the NCAM1-NCAM1
        plt.close()
        plt.rcParams.update({
            'axes.labelsize': 16,
            'xtick.labelsize': 14,
            'ytick.labelsize': 14,
            'text.color': 'black'
        })

        x = global_moranI_all_AD_CT['global_moranI_AD']
        y = global_moranI_all_AD_CT['global_moranI_CT']
        type_cur = global_moranI_all_AD_CT['type_CT']
        unique_types = type_cur.unique()
        # palette = sns.color_palette("tab10", n_colors=len(unique_types))
        # type2color = {tp: palette[i] for i, tp in enumerate(unique_types)}
        # colors = type_cur.map(type2color)
        # custom_palette = ['#95E1D3', '#F38181']
        custom_palette = ['#3FC1C9', '#FC5185']
        type2color = {tp: custom_palette[i] for i, tp in enumerate(unique_types)}
        colors = type_cur.map(type2color)
        x_min_abs = abs(x.min())
        y_min_abs = abs(y.min())
        x_max_abs = abs(x.max())
        y_max_abs = abs(y.max())
        max_all = max(x_max_abs, y_max_abs, x_min_abs, y_min_abs)
        min_val = -max_all
        max_val = max_all
        plt.figure(figsize=(10, 8))
        plt.axhline(0, color='grey', linestyle='--', linewidth=1, zorder=0)
        plt.axvline(0, color='grey', linestyle='--', linewidth=1, zorder=0)
        plt.plot([min_val, max_val], [min_val, max_val], color='black', linestyle='--', linewidth=1, zorder=0)
        scatter = plt.scatter(x, y, c=colors, alpha=0.7, s=10)
        plt.xlabel("Bi-variate Moran's I of \ngene expression in AD group", fontsize=16)
        plt.ylabel("Bi-variate Moran's I of \ngene expression in CT group", fontsize=16)
        plt.xlim(min_val * 1.03, max_val * 1.03)
        plt.ylim(min_val * 1.03, max_val * 1.03)
        plt.gca().set_aspect('equal', adjustable='box')
        texts = []
        for LR_cur_index in range(global_moranI_all_AD_CT_positive.shape[0]):
            cur_df = global_moranI_all_AD_CT_positive.iloc[LR_cur_index, :]
            x_val = float(cur_df['global_moranI_AD'])
            y_val = float(cur_df['global_moranI_CT'])
            type_cur_highlight = cur_df['type_CT']
            LR_cur = cur_df['Gene1'] + "-" + cur_df['Gene2']
            if LR_cur == "NCAM1-NCAM1":
                plt.scatter([x_val], [y_val], color=type2color['LR'], edgecolor='yellow', s=50, linewidth=3, zorder=5)
            if LR_cur == "CNTN2-DPYSL2":
                plt.scatter([x_val], [y_val], color=type2color['GP'], edgecolor='blue', s=50, linewidth=3, zorder=5)
            if LR_cur == "NCAM1-DPYSL2":
                plt.scatter([x_val], [y_val], color=type2color['GP'], edgecolor='green', s=50, linewidth=3, zorder=5)
            # if LR_cur == "NCAM1-NCAM1":
            #     texts.append(
            #         plt.text(x_val, y_val, LR_cur,
            #                  fontsize=14, color="black", ha='left', va='bottom',
            #                  bbox=dict(boxstyle="round,pad=0.3", edgecolor="gray", facecolor="white", alpha=0.7))
            #     )
        adjust_text(texts, arrowprops=dict(arrowstyle="->", color='gray', lw=0.5))

        from matplotlib.patches import Patch

        legend_elements = [Patch(facecolor=type2color[tp], label=tp) for tp in unique_types]
        plt.legend(handles=legend_elements, title='Gene pair type', fontsize=12, title_fontsize=13, loc='best')

        filename = f"/Users/junjietang/Desktop/AD project/Xenium_data_V1/global_moranI_scatter_condition_{metagene1}_{metagene2}_NCAM1-NCAM1-highlight.pdf"
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()

    ##############################################################################################################################
    ##HIc-to-Hic spatial correlation
    if if_Hic_Hic:
        global_moranI_Hic_top1top2_sample = {}
        # global_moranI_Hic_top1tail2_sample = {}
        # global_moranI_Hic_tail1top2_sample = {}
        ##
        if if_modify_moran:
            # topgene_metagene1_mean_alltissue = np.mean(compartment_subset.loc[:, topgene_metagene1], axis=0)
            # topgene_metagene2_mean_alltissue = np.mean(compartment_subset.loc[:, topgene_metagene2], axis=0)
            # topgene_metagene1_std_alltissue = np.std(compartment_subset.loc[:, topgene_metagene1], axis=0)
            # topgene_metagene2_std_alltissue = np.std(compartment_subset.loc[:, topgene_metagene2], axis=0)
            topgene_assogene_mean_alltissue = np.mean(compartment_subset.loc[:, assogene], axis=0)
            topgene_assogene_std_alltissue = np.std(compartment_subset.loc[:, assogene], axis=0)
            ##
            # topgene_metagene1_mean_alltissue = pd.DataFrame(topgene_metagene1_mean_alltissue.T, index=topgene_metagene1,
            #                                                 columns=["mean"])
            # topgene_metagene2_mean_alltissue = pd.DataFrame(topgene_metagene2_mean_alltissue.T, index=topgene_metagene2,
            #                                                 columns=["mean"])
            # topgene_metagene1_std_alltissue = pd.DataFrame(topgene_metagene1_std_alltissue.T, index=topgene_metagene1,
            #                                                columns=["std"])
            # topgene_metagene2_std_alltissue = pd.DataFrame(topgene_metagene2_std_alltissue.T, index=topgene_metagene2,
            #                                                columns=["std"])
            topgene_assogene_mean_alltissue = pd.DataFrame(topgene_assogene_mean_alltissue.T, index=assogene,
                                                            columns=["mean"])
            topgene_assogene_std_alltissue = pd.DataFrame(topgene_assogene_std_alltissue.T, index=assogene,
                                                              columns=["std"])
        ##
        sample_unique = np.unique(adata_Xenium.obs['sample_id'])
        for sample_cur in sample_unique:
            print(sample_cur)
            cellindex_cursample = np.where(adata_Xenium.obs['sample_id'] == sample_cur)[0]
            geneexp_cursample = compartment_subset.iloc[cellindex_cursample, :]
            spatial_coord_cursample = adata_Xenium.obsm['spatial'][cellindex_cursample, :]
            ##
            spatial_coord_cursample_cur = spatial_coord_cursample
            ##
            w = libpysal.weights.KNN.from_array(spatial_coord_cursample_cur, k=num_neighbor)
            w.transform = 'r'
            ##Top genes
            global_moran_I_cursample = np.zeros((len(assogene), len(assogene)))
            global_moran_I_cursample = pd.DataFrame(global_moran_I_cursample,index = assogene,
                                                    columns=assogene)
            for GP_index in range(topGP_merge.shape[0]):
                topGP_merge_cur = topGP_merge.iloc[GP_index, :]
                gene1 = topGP_merge_cur['gene1']
                gene2 = topGP_merge_cur['gene2']
                ##
                geneexp_cursample_cur = geneexp_cursample.loc[:,
                                        [gene1, gene2]].values
                ##
                moran_bv = Moran_BV_modified(x=geneexp_cursample_cur[:, 0], y=geneexp_cursample_cur[:, 1],
                                             w=w,
                                             x_mean=float(
                                                 topgene_assogene_mean_alltissue.loc[
                                                     gene1]),
                                             x_std=float(
                                                 topgene_assogene_std_alltissue.loc[
                                                     gene1]),
                                             y_mean=float(
                                                 topgene_assogene_mean_alltissue.loc[
                                                     gene2]),
                                             y_std=float(
                                                 topgene_assogene_std_alltissue.loc[
                                                     gene2]))
                ##
                global_moran_I_cursample.loc[gene1, gene2] = moran_bv.I
            global_moranI_Hic_top1top2_sample[sample_cur] = global_moran_I_cursample
        #
        ##Condition-specific
        sample_unique_AD = ["AD-1", "AD-2"]
        sample_unique_CT = ["CT-1", "CT-2"]
        global_moranI_Hic_top1top2_AD = np.zeros((len(assogene), len(assogene)))
        for sample_cur in sample_unique_AD:
            global_moranI_Hic_top1top2_AD += global_moranI_Hic_top1top2_sample[sample_cur].values
        global_moranI_Hic_top1top2_AD /= len(sample_unique_AD)
        global_moranI_Hic_top1top2_AD = pd.DataFrame(global_moranI_Hic_top1top2_AD, index=assogene,
                                                     columns=assogene)
        ##Turn the global_moranI_Hic_top1top2 into a long dataframe
        global_moranI_Hic_top1top2_AD_long = global_moranI_Hic_top1top2_AD.stack().reset_index()
        global_moranI_Hic_top1top2_AD_long.columns = ["Gene1", "Gene2", "global_moranI_Hic"]
        global_moranI_Hic_top1top2_AD_long['Group'] = "Top1Top2_AD"
        ##
        global_moranI_Hic_top1top2_AD_long = global_moranI_Hic_top1top2_AD_long.loc[~(global_moranI_Hic_top1top2_AD_long['global_moranI_Hic'] == 0),:]
        ##
        global_moranI_Hic_top1top2_CT = np.zeros((len(assogene), len(assogene)))
        for sample_cur in sample_unique_CT:
            global_moranI_Hic_top1top2_CT += global_moranI_Hic_top1top2_sample[sample_cur].values
        global_moranI_Hic_top1top2_CT /= len(sample_unique_CT)
        global_moranI_Hic_top1top2_CT = pd.DataFrame(global_moranI_Hic_top1top2_CT, index=assogene,
                                                     columns=assogene)
        ##Turn the global_moranI_Hic_top1top2 into a long dataframe
        global_moranI_Hic_top1top2_CT_long = global_moranI_Hic_top1top2_CT.stack().reset_index()
        global_moranI_Hic_top1top2_CT_long.columns = ["Gene1", "Gene2", "global_moranI_Hic"]
        global_moranI_Hic_top1top2_CT_long['Group'] = "Top1Top2_CT"
        global_moranI_Hic_top1top2_CT_long = global_moranI_Hic_top1top2_CT_long.loc[~(global_moranI_Hic_top1top2_CT_long['global_moranI_Hic'] == 0),:]
        ##
        ##
        global_moranI_Hic_all = pd.concat([global_moranI_Hic_top1top2_AD_long, global_moranI_Hic_top1top2_CT_long],
                                          axis=0)
        ##
        ##Just keep the ligand-receptor pair
        choose_index = []
        for i in range(global_moranI_Hic_all.shape[0]):
            gene1 = global_moranI_Hic_all.iloc[i, 0]
            gene2 = global_moranI_Hic_all.iloc[i, 1]
            if np.sum((LR_Cellchat['ligand'] == gene2) * (LR_Cellchat['receptor'] == gene1)) > 0 or np.sum(
                    (LR_Cellchat['ligand'] == gene1) * (LR_Cellchat['receptor'] == gene2)) > 0:
                choose_index.append(i)
        #global_moranI_Hic_all = global_moranI_Hic_all.iloc[choose_index, :]
        global_moranI_Hic_all['type'] = "GP"
        global_moranI_Hic_all['type'].iloc[choose_index] = "LR"
        ##Keep the global_moranI_Hic_all['Group'] == "Top1Top2_AD" or "Top1Top2_CT"
        global_moranI_Hic_all = global_moranI_Hic_all[
            (global_moranI_Hic_all['Group'] == "Top1Top2_AD") | (global_moranI_Hic_all['Group'] == "Top1Top2_CT")]
        ##save global_moranI_Hic_all to csv
        global_moranI_Hic_all.to_csv(
            "/Users/junjietang/Desktop/AD project/Xenium_data_V1/global_moranI_Hic_all_" + str(metagene1) + "_" + str(
                metagene2) + "_ABcompartmentscore.csv", index=False)
        ##Show the density plot of global_moranI_Hic_all['global_moranI_Hic'], grouped by global_moranI_Hic_all['Group']
        import seaborn as sns
        import matplotlib.pyplot as plt
        import matplotlib
        ##change the "Top1Top2_AD" to "AD" and "Top1Top2_CT" to "CT"
        global_moranI_Hic_all['Group'] = global_moranI_Hic_all['Group'].replace({"Top1Top2_AD": "AD", "Top1Top2_CT": "CT"})
        matplotlib.use('TkAgg')
        plt.rcParams.update({
            'axes.titlesize': 20,
            'axes.labelsize': 18,
            'xtick.labelsize': 16,
            'ytick.labelsize': 16,
            'legend.title_fontsize': 16,
            'legend.fontsize': 16
        })

        plt.figure(figsize=(7, 8))
        ax = sns.kdeplot(
            data=global_moranI_Hic_all,
            x='global_moranI_Hic',
            hue='Group',
            fill=True,
            common_norm=False
        )

        ax.set_xlabel("Bi-variate Moran's I of A/B compartment score \nbetween associated ligand-receptor pair", fontsize=18)
        ax.set_ylabel("Density", fontsize=18)

        leg = ax.get_legend()
        leg.set_title("Group")
        plt.setp(leg.get_title(), fontsize=22)
        plt.setp(leg.get_texts(), fontsize=22)

        plt.tight_layout()
        plt.savefig("/Users/junjietang/Desktop/AD project/Xenium_data_V1/global_moranI_Hic_density_condition_" + str(
            metagene1) + "_" + str(metagene2) + ".pdf", dpi=300, bbox_inches='tight')
        plt.close()
        ##
        global_moranI_Hic_all_AD = global_moranI_Hic_all[global_moranI_Hic_all['Group'] == "AD"]
        global_moranI_Hic_all_CT = global_moranI_Hic_all[global_moranI_Hic_all['Group'] == "CT"]

        ##Show a scatter plot bettwen global_moranI_Hic_all_AD['global_moranI_Hic'] and global_moranI_Hic_all_CT['global_moranI_Hic']
        ##Merge global_moranI_Hic_all_AD and global_moranI_Hic_all_CT
        global_moranI_Hic_all_AD_CT = pd.merge(global_moranI_Hic_all_AD, global_moranI_Hic_all_CT,
                                               on=['Gene1', 'Gene2'],
                                               suffixes=('_AD', '_CT'))
        global_moranI_Hic_all_AD_CT_LR = global_moranI_Hic_all_AD_CT.loc[global_moranI_Hic_all_AD_CT['type_CT'] == "LR",:]
        ##
        global_moranI_Hic_all_AD_CT_positive = global_moranI_Hic_all_AD_CT.loc[~(
                    (global_moranI_Hic_all_AD_CT['global_moranI_Hic_AD'] < 0) * (
                        global_moranI_Hic_all_AD_CT['global_moranI_Hic_CT'] < 0)), :]
        ##

        ##Only highlight the NCAM1-NCAM1
        plt.rcParams.update({
            'axes.labelsize': 16,
            'xtick.labelsize': 14,
            'ytick.labelsize': 14,
            'text.color': 'black'
        })

        x = global_moranI_Hic_all_AD_CT['global_moranI_Hic_AD']
        y = global_moranI_Hic_all_AD_CT['global_moranI_Hic_CT']
        type_cur = global_moranI_Hic_all_AD_CT['type_CT']

        unique_types = type_cur.unique()
        # palette = sns.color_palette("tab10", n_colors=len(unique_types))
        # type2color = {tp: palette[i] for i, tp in enumerate(unique_types)}
        # colors = type_cur.map(type2color)
        # custom_palette = ['#95E1D3', '#F38181']
        custom_palette = ['#3FC1C9', '#FC5185']
        type2color = {tp: custom_palette[i] for i, tp in enumerate(unique_types)}
        colors = type_cur.map(type2color)

        x_min_abs = abs(x.min())
        y_min_abs = abs(y.min())
        x_max_abs = abs(x.max())
        y_max_abs = abs(y.max())
        max_all = max(x_max_abs, y_max_abs, x_min_abs, y_min_abs)
        min_val = (-1) * max_all
        max_val = max_all

        plt.figure(figsize=(10, 8))
        plt.axhline(0, color='grey', linestyle='--', linewidth=1,zorder=0)
        plt.axvline(0, color='grey', linestyle='--', linewidth=1,zorder=0)
        plt.plot([min_val, max_val], [min_val, max_val], color='black', linestyle='--', linewidth=1,zorder=0)

        scatter = plt.scatter(x, y, c=colors, alpha=0.7, s=10)

        plt.xlabel("Bi-variate Moran's I of \n A/B compartment score in AD group", fontsize=16)
        plt.ylabel("Bi-variate Moran's I of \n A/B compartment score in CT group", fontsize=16)

        plt.xlim(min_val * 1.03, max_val * 1.03)
        plt.ylim(min_val * 1.03, max_val * 1.03)
        plt.gca().set_aspect('equal', adjustable='box')

        texts = []
        for LR_cur_index in range(global_moranI_Hic_all_AD_CT_positive.shape[0]):
            cur_df = global_moranI_Hic_all_AD_CT_positive.iloc[LR_cur_index, :]
            x_val = float(cur_df['global_moranI_Hic_AD'])
            y_val = float(cur_df['global_moranI_Hic_CT'])
            type_cur_highlight = cur_df['type_CT']
            LR_cur = cur_df['Gene1'] + "-" + cur_df['Gene2']
            if LR_cur == "NCAM1-NCAM1":
                plt.scatter([x_val], [y_val], color=type2color['LR'], edgecolor='yellow', s=50, linewidth=3,
                            zorder=5)
            if LR_cur == "CNTN2-DPYSL2":
                plt.scatter([x_val], [y_val], color=type2color['GP'], edgecolor='blue', s=50, linewidth=3, zorder=5)
            if LR_cur == "NCAM1-DPYSL2":
                plt.scatter([x_val], [y_val], color=type2color['GP'], edgecolor='green', s=50, linewidth=3, zorder=5)
        adjust_text(texts, arrowprops=dict(arrowstyle="->", color='gray', lw=0.5))

        filename = f"/Users/junjietang/Desktop/AD project/Xenium_data_V1/global_moranI_Hic_scatter_condition_{metagene1}_{metagene2}_NCAM1-NCAM1-highlight.pdf"
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()


        ##just highlight the positive part
        ##Only highlight the NCAM1-NCAM1
        global_moranI_Hic_all_AD_CT_positivesel = global_moranI_Hic_all_AD_CT
        global_moranI_Hic_all_AD_CT_positivesel['type_CT'].loc[~((global_moranI_Hic_all_AD_CT['global_moranI_Hic_AD'] >=0) * (global_moranI_Hic_all_AD_CT['global_moranI_Hic_CT'] >=0))] = "NPGP"
        plt.rcParams.update({
            'axes.labelsize': 16,
            'xtick.labelsize': 14,
            'ytick.labelsize': 14,
            'text.color': 'black'
        })

        x = global_moranI_Hic_all_AD_CT_positivesel['global_moranI_Hic_AD']
        y = global_moranI_Hic_all_AD_CT_positivesel['global_moranI_Hic_CT']
        type_cur = global_moranI_Hic_all_AD_CT_positivesel['type_CT']

        unique_types = type_cur.unique()
        # palette = sns.color_palette("tab10", n_colors=len(unique_types))
        # type2color = {tp: palette[i] for i, tp in enumerate(unique_types)}
        # colors = type_cur.map(type2color)
        # custom_palette = ['#95E1D3', '#F38181']
        custom_palette = ['#EEEEEE', '#3FC1C9', '#FC5185']
        type2color = {tp: custom_palette[i] for i, tp in enumerate(unique_types)}
        colors = type_cur.map(type2color)

        x_min_abs = abs(x.min())
        y_min_abs = abs(y.min())
        x_max_abs = abs(x.max())
        y_max_abs = abs(y.max())
        max_all = max(x_max_abs, y_max_abs, x_min_abs, y_min_abs)
        min_val = (-1) * max_all
        max_val = max_all

        plt.figure(figsize=(10, 8))
        plt.axhline(0, color='grey', linestyle='--', linewidth=1,zorder=0)
        plt.axvline(0, color='grey', linestyle='--', linewidth=1,zorder=0)
        plt.plot([min_val, max_val], [min_val, max_val], color='black', linestyle='--', linewidth=1,zorder=0)

        scatter = plt.scatter(x, y, c=colors, alpha=0.7, s=10)

        plt.xlabel("Bi-variate Moran's I of \n A/B compartment score in AD group", fontsize=16)
        plt.ylabel("Bi-variate Moran's I of \n A/B compartment score in CT group", fontsize=16)

        plt.xlim(min_val * 1.03, max_val * 1.03)
        plt.ylim(min_val * 1.03, max_val * 1.03)
        plt.gca().set_aspect('equal', adjustable='box')

        texts = []
        for LR_cur_index in range(global_moranI_Hic_all_AD_CT_positive.shape[0]):
            cur_df = global_moranI_Hic_all_AD_CT_positive.iloc[LR_cur_index, :]
            x_val = float(cur_df['global_moranI_Hic_AD'])
            y_val = float(cur_df['global_moranI_Hic_CT'])
            type_cur_highlight = cur_df['type_CT']
            LR_cur = cur_df['Gene1'] + "-" + cur_df['Gene2']
            if LR_cur == "NCAM1-NCAM1":
                plt.scatter([x_val], [y_val], color=type2color['LR'], edgecolor='yellow', s=50, linewidth=3,
                            zorder=5)

        adjust_text(texts, arrowprops=dict(arrowstyle="->", color='gray', lw=0.5))

        filename = f"/Users/junjietang/Desktop/AD project/Xenium_data_V1/global_moranI_Hic_scatter_condition_{metagene1}_{metagene2}_NCAM1-NCAM1-highlight_positive.pdf"
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()
    ##
    global_moranI_Hic_all_AD_CT_positive = global_moranI_Hic_all_AD_CT.loc[~(
            (global_moranI_Hic_all_AD_CT['global_moranI_Hic_AD'] < 0) * (
            global_moranI_Hic_all_AD_CT['global_moranI_Hic_CT'] < 0)), :]
    global_moranI_Hic_all_AD_CT_choose = global_moranI_Hic_all_AD_CT_positive.loc[global_moranI_Hic_all_AD_CT_positive['global_moranI_Hic_AD'] * 2 < global_moranI_Hic_all_AD_CT_positive['global_moranI_Hic_CT'],:]
    global_moranI_Hic_all_AD_CT_choose["GP_name"] = global_moranI_Hic_all_AD_CT_choose['Gene1'] + "-" + global_moranI_Hic_all_AD_CT_choose['Gene2']
    ##
    global_moranI_all_AD_CT_positive = global_moranI_all_AD_CT.loc[~(
            (global_moranI_all_AD_CT['global_moranI_AD'] < 0) * (
            global_moranI_all_AD_CT['global_moranI_CT'] < 0)), :]
    global_moranI_all_AD_CT_choose = global_moranI_all_AD_CT_positive.loc[global_moranI_all_AD_CT_positive['global_moranI_AD'] * 2 < global_moranI_all_AD_CT_positive['global_moranI_CT'],:]
    global_moranI_all_AD_CT_choose["GP_name"] = global_moranI_all_AD_CT_choose['Gene1'] + "-" + global_moranI_all_AD_CT_choose['Gene2']
    ##
    GP_name_intersect = np.intersect1d(global_moranI_Hic_all_AD_CT_choose["GP_name"].values, global_moranI_all_AD_CT_choose["GP_name"].values)
    ##
    global_moranI_Hic_all_AD_CT_choose = global_moranI_Hic_all_AD_CT_choose.loc[global_moranI_Hic_all_AD_CT_choose["GP_name"].isin(GP_name_intersect),:]
    global_moranI_all_AD_CT_choose = global_moranI_all_AD_CT_choose.loc[global_moranI_all_AD_CT_choose["GP_name"].isin(GP_name_intersect),:]
    GP_use = global_moranI_Hic_all_AD_CT_choose[['Gene1', 'Gene2']].values
    ##Save GP_use
    GP_use = pd.DataFrame(GP_use, columns=["Gene1", "Gene2"])
    GP_use.to_csv(
        "/Users/junjietang/Desktop/AD project/Xenium_data_V1/GP_use_" + str(metagene1) + "_" + str(
            metagene2) + ".csv", index=False)

    if if_gene_Hic:
        ##Gene-to-Gene spatial correlation
        global_moranI_top1top2_sample = {}
        ##
        if if_modify_moran:
            topgene_assogene_mean_alltissue = np.mean(adata_Xenium[:, assogene].X, axis=0)
            topgene_assogene_std_alltissue = np.std(adata_Xenium[:, assogene].X.toarray(), axis=0)
            topgene_assogene_Hic_mean_alltissue = np.mean(compartment_subset.loc[:, assogene], axis=0)
            topgene_assogene_Hic_std_alltissue = np.std(compartment_subset.loc[:, assogene], axis=0)
            ##
            topgene_assogene_mean_alltissue = pd.DataFrame(topgene_assogene_mean_alltissue.T, index=assogene,
                                                            columns=["mean"])
            topgene_assogene_std_alltissue = pd.DataFrame(topgene_assogene_std_alltissue.T, index=assogene,
                                                           columns=["std"])
            topgene_assogene_Hic_mean_alltissue = pd.DataFrame(topgene_assogene_Hic_mean_alltissue.T, index=assogene,
                                                               columns=["mean"])
            topgene_assogene_Hic_std_alltissue = pd.DataFrame(topgene_assogene_Hic_std_alltissue.T, index=assogene,
                                                                columns=["std"])
        ##
        sample_unique = np.unique(adata_Xenium.obs['sample_id'])
        for sample_cur in sample_unique:
            print(sample_cur)
            cellindex_cursample = np.where(adata_Xenium.obs['sample_id'] == sample_cur)[0]
            geneexp_cursample = adata_Xenium.X[cellindex_cursample, :].toarray()
            geneexp_cursample = pd.DataFrame(geneexp_cursample,
                                             index=adata_Xenium.obs['orig.ident'].index[cellindex_cursample],
                                             columns=adata_Xenium.var_names)
            Hic_cursample = compartment_subset.iloc[cellindex_cursample, :]
            spatial_coord_cursample = adata_Xenium.obsm['spatial'][cellindex_cursample, :]
            ##
            spatial_coord_cursample_cur = spatial_coord_cursample
            ##
            w = libpysal.weights.KNN.from_array(spatial_coord_cursample_cur, k=num_neighbor)
            w.transform = 'r'
            ##Top genes
            global_moran_I_cursample = np.zeros((len(assogene), len(assogene)))
            global_moran_I_cursample = pd.DataFrame(global_moran_I_cursample,index = assogene,
                                                    columns=assogene)
            for GP_index in range(topGP_merge.shape[0]):
                topGP_merge_cur = topGP_merge.iloc[GP_index, :]
                ##Direction 1
                gene1 = topGP_merge_cur['gene1']
                gene2 = topGP_merge_cur['gene2']
                ##
                geneexp_cursample_cur = np.vstack([geneexp_cursample.loc[:,
                                        [gene1]].values[:,0],Hic_cursample.loc[:,gene2].values]).T
                ##
                moran_bv = Moran_BV_modified(x=geneexp_cursample_cur[:, 0], y=geneexp_cursample_cur[:, 1],
                                             w=w,
                                             x_mean=float(
                                                 topgene_assogene_mean_alltissue.loc[
                                                     gene1]),
                                             x_std=float(
                                                 topgene_assogene_std_alltissue.loc[
                                                     gene1]),
                                             y_mean=float(
                                                 topgene_assogene_Hic_mean_alltissue.loc[
                                                     gene2]),
                                             y_std=float(
                                                 topgene_assogene_Hic_std_alltissue.loc[
                                                     gene2]))
                ##
                global_moran_I_cursample.loc[gene1, gene2] = moran_bv.I
                ##Direction 2
                gene2 = topGP_merge_cur['gene1']
                gene1 = topGP_merge_cur['gene2']
                ##
                geneexp_cursample_cur = np.vstack([geneexp_cursample.loc[:,
                                        [gene1]].values[:,0],Hic_cursample.loc[:,gene2].values]).T
                ##
                moran_bv = Moran_BV_modified(x=geneexp_cursample_cur[:, 0], y=geneexp_cursample_cur[:, 1],
                                             w=w,
                                             x_mean=float(
                                                 topgene_assogene_mean_alltissue.loc[
                                                     gene1]),
                                             x_std=float(
                                                 topgene_assogene_std_alltissue.loc[
                                                     gene1]),
                                             y_mean=float(
                                                 topgene_assogene_Hic_mean_alltissue.loc[
                                                     gene2]),
                                             y_std=float(
                                                 topgene_assogene_Hic_std_alltissue.loc[
                                                     gene2]))
                ##
                global_moran_I_cursample.loc[gene1, gene2] = moran_bv.I
            global_moranI_top1top2_sample[sample_cur] = global_moran_I_cursample

        ##Condition-specific
        sample_unique_AD = ["AD-1", "AD-2"]
        sample_unique_CT = ["CT-1", "CT-2"]
        global_moranI_top1top2_AD = np.zeros((len(assogene), len(assogene)))
        for sample_cur in sample_unique_AD:
            global_moranI_top1top2_AD += global_moranI_top1top2_sample[sample_cur].values
        global_moranI_top1top2_AD /= len(sample_unique_AD)
        global_moranI_top1top2_AD = pd.DataFrame(global_moranI_top1top2_AD, index=assogene,
                                                 columns=assogene)
        ##Turn the global_moranI_top1top2 into a long dataframe
        global_moranI_top1top2_AD_long = global_moranI_top1top2_AD.stack().reset_index()
        global_moranI_top1top2_AD_long.columns = ["Gene1", "Gene2", "global_moranI"]
        global_moranI_top1top2_AD_long['Group'] = "Top1Top2_AD"
        ##
        global_moranI_top1top2_AD_long = global_moranI_top1top2_AD_long.loc[~(global_moranI_top1top2_AD_long['global_moranI'] == 0),:]
        ##
        global_moranI_top1top2_CT = np.zeros((len(assogene), len(assogene)))
        for sample_cur in sample_unique_CT:
            global_moranI_top1top2_CT += global_moranI_top1top2_sample[sample_cur].values
        global_moranI_top1top2_CT /= len(sample_unique_CT)
        global_moranI_top1top2_CT = pd.DataFrame(global_moranI_top1top2_CT, index=assogene,
                                                 columns=assogene)
        ##Turn the global_moranI_top1top2 into a long dataframe
        global_moranI_top1top2_CT_long = global_moranI_top1top2_CT.stack().reset_index()
        global_moranI_top1top2_CT_long.columns = ["Gene1", "Gene2", "global_moranI"]
        global_moranI_top1top2_CT_long['Group'] = "Top1Top2_CT"
        global_moranI_top1top2_CT_long = global_moranI_top1top2_CT_long.loc[~(global_moranI_top1top2_CT_long['global_moranI'] == 0),:]
        ##
        global_moranI_all = pd.concat([global_moranI_top1top2_AD_long, global_moranI_top1top2_CT_long], axis=0)
        ##
        choose_index = []
        for i in range(global_moranI_all.shape[0]):
            gene1 = global_moranI_all.iloc[i, 0]
            gene2 = global_moranI_all.iloc[i, 1]
            if np.sum((LR_Cellchat['ligand'] == gene2) * (LR_Cellchat['receptor'] == gene1)) > 0 or np.sum(
                    (LR_Cellchat['ligand'] == gene1) * (LR_Cellchat['receptor'] == gene2)) > 0:
                choose_index.append(i)
        # global_moranI_all = global_moranI_all.iloc[choose_index, :]
        global_moranI_all['type'] = "GP"
        global_moranI_all['type'].iloc[choose_index] = "LR"

        ##Keep the global_moranI_all['Group'] == "Top1Top2_AD" or "Top1Top2_CT"
        global_moranI_all = global_moranI_all[
            (global_moranI_all['Group'] == "Top1Top2_AD") | (global_moranI_all['Group'] == "Top1Top2_CT")]
        ##change the "Top1Top2_AD" to "AD" and "Top1Top2_CT" to "CT"
        global_moranI_all['Group'] = global_moranI_all['Group'].replace({"Top1Top2_AD": "AD", "Top1Top2_CT": "CT"})
        ##save the global_moranI_all to csv
        global_moranI_all.to_csv(
            "/Users/junjietang/Desktop/AD project/Xenium_data_V1/global_moranI_all_" + str(metagene1) + "_" + str(
                metagene2) + "_genexpressiontoHic.csv", index=False)
        ##load the global_moranI_all
        global_moranI_all = pd.read_csv(
            "/Users/junjietang/Desktop/AD project/Xenium_data_V1/global_moranI_all_" + str(metagene1) + "_" + str(
                metagene2) + "_genexpressiontoHic.csv")

        ##Show the density plot of global_moranI_all['global_moranI'], grouped by global_moranI_all['Group']
        import seaborn as sns
        import matplotlib.pyplot as plt
        import matplotlib
        ##
        matplotlib.use('TkAgg')
        plt.rcParams.update({
            'axes.titlesize': 20,
            'axes.labelsize': 18,
            'xtick.labelsize': 16,
            'ytick.labelsize': 16,
            'legend.title_fontsize': 16,
            'legend.fontsize': 16
        })

        plt.figure(figsize=(7, 8))
        ax = sns.kdeplot(
            data=global_moranI_all,
            x='global_moranI',
            hue='Group',
            fill=True,
            common_norm=False
        )

        # ax.set_xlabel("Bi-variate Moran's I of gene expression \nbetween associated ligand-receptor pair", fontsize=18)
        ax.set_xlabel("Bi-variate Moran's I of gene expression and AB score \nbetween associated gene pairs", fontsize=18)
        ax.set_ylabel("Density", fontsize=18)

        leg = ax.get_legend()
        leg.set_title("Group")
        plt.setp(leg.get_title(), fontsize=22)
        plt.setp(leg.get_texts(), fontsize=22)

        plt.tight_layout()

        plt.savefig("/Users/junjietang/Desktop/AD project/Xenium_data_V1/global_moranI_geneHic_density_condition_" + str(
            metagene1) + "_" + str(metagene2) + ".pdf", dpi=300, bbox_inches='tight')
        ##
        global_moranI_all_AD = global_moranI_all[global_moranI_all['Group'] == "AD"]
        global_moranI_all_CT = global_moranI_all[global_moranI_all['Group'] == "CT"]

        ##Show a scatter plot bettwen global_moranI_all_AD['global_moranI'] and global_moranI_all_CT['global_moranI']
        ##Merge global_moranI_all_AD and global_moranI_all_CT
        global_moranI_all_AD_CT = pd.merge(global_moranI_all_AD, global_moranI_all_CT, on=['Gene1', 'Gene2'],
                                           suffixes=('_AD', '_CT'))
        ##
        global_moranI_all_AD_CT_positive = global_moranI_all_AD_CT.loc[~(
                    (global_moranI_all_AD_CT['global_moranI_AD'] < 0) * (
                        global_moranI_all_AD_CT['global_moranI_CT'] < 0)), :]
        global_moranI_all_AD_CT_positive = global_moranI_all_AD_CT_positive.iloc[np.argsort(
            np.max(global_moranI_all_AD_CT_positive[['global_moranI_AD', 'global_moranI_CT']], axis=1))[
                                                                                 ::-1][0:20], :]
        ##
        ##
        import matplotlib.pyplot as plt
        from adjustText import adjust_text
        import seaborn as sns
        from matplotlib.cm import get_cmap

        plt.rcParams.update({
            'axes.labelsize': 16,
            'xtick.labelsize': 14,
            'ytick.labelsize': 14,
            'text.color': 'black'
        })

        x = global_moranI_all_AD_CT['global_moranI_AD']
        y = global_moranI_all_AD_CT['global_moranI_CT']
        type_cur = global_moranI_all_AD_CT['type_CT']

        unique_types = type_cur.unique()
        palette = sns.color_palette("tab10", n_colors=len(unique_types))
        type2color = {tp: palette[i] for i, tp in enumerate(unique_types)}
        colors = type_cur.map(type2color)

        x_min_abs = abs(x.min())
        y_min_abs = abs(y.min())
        x_max_abs = abs(x.max())
        y_max_abs = abs(y.max())
        max_all = max(x_max_abs, y_max_abs, x_min_abs, y_min_abs)
        min_val = -max_all
        max_val = max_all

        plt.figure(figsize=(10, 8))

        scatter = plt.scatter(x, y, c=colors, alpha=0.7, s=5)

        plt.plot([min_val, max_val], [min_val, max_val], color='red', linestyle='--', linewidth=1)

        plt.xlabel("Bi-variate Moran's I between \ngene expression and AB score in AD group", fontsize=16)
        plt.ylabel("Bi-variate Moran's I between \ngene expression and AB score in CT group", fontsize=16)
        plt.xlim(min_val * 1.1, max_val * 1.1)
        plt.ylim(min_val * 1.1, max_val * 1.1)
        plt.gca().set_aspect('equal', adjustable='box')

        texts = []
        for LR_cur_index in range(global_moranI_all_AD_CT_positive.shape[0]):
            cur_df = global_moranI_all_AD_CT_positive.iloc[LR_cur_index, :]
            x_val = float(cur_df['global_moranI_AD'])
            y_val = float(cur_df['global_moranI_CT'])
            type_cur_highlight = cur_df['type_CT']
            LR_cur = cur_df['Gene1'] + "-" + cur_df['Gene2']
            if type_cur_highlight == "LR":
                texts.append(
                    plt.text(x_val, y_val, LR_cur,
                             fontsize=14, color="red", ha='left', va='bottom',
                             bbox=dict(boxstyle="round,pad=0.3", edgecolor="gray", facecolor="white", alpha=0.7))
                )
            else:
                texts.append(
                    plt.text(x_val, y_val, LR_cur,
                             fontsize=14, color="black", ha='left', va='bottom',
                             bbox=dict(boxstyle="round,pad=0.3", edgecolor="gray", facecolor="white", alpha=0.7))
                )
        adjust_text(texts, arrowprops=dict(arrowstyle="->", color='gray', lw=0.5))

        plt.axhline(0, color='black', linestyle='--', linewidth=1)
        plt.axvline(0, color='black', linestyle='--', linewidth=1)

        from matplotlib.patches import Patch

        legend_elements = [Patch(facecolor=type2color[tp], label=tp) for tp in unique_types]
        plt.legend(handles=legend_elements, title='Cell type', fontsize=12, title_fontsize=13, loc='best')

        filename = f"/Users/junjietang/Desktop/AD project/Xenium_data_V1/global_moranI_geneHic_scatter_condition_{metagene1}_{metagene2}.pdf"
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        # plt.show()
        ##just highlight the NCAM1-NCAM1
        plt.close()
        plt.rcParams.update({
            'axes.labelsize': 16,
            'xtick.labelsize': 14,
            'ytick.labelsize': 14,
            'text.color': 'black'
        })

        x = global_moranI_all_AD_CT['global_moranI_AD']
        y = global_moranI_all_AD_CT['global_moranI_CT']
        type_cur = global_moranI_all_AD_CT['type_CT']
        unique_types = type_cur.unique()
        # palette = sns.color_palette("tab10", n_colors=len(unique_types))
        # type2color = {tp: palette[i] for i, tp in enumerate(unique_types)}
        # colors = type_cur.map(type2color)
        # custom_palette = ['#95E1D3', '#F38181']
        custom_palette = ['#3FC1C9', '#FC5185']
        type2color = {tp: custom_palette[i] for i, tp in enumerate(unique_types)}
        colors = type_cur.map(type2color)

        x_min_abs = abs(x.min())
        y_min_abs = abs(y.min())
        x_max_abs = abs(x.max())
        y_max_abs = abs(y.max())
        max_all = max(x_max_abs, y_max_abs, x_min_abs, y_min_abs)
        min_val = -max_all
        max_val = max_all
        plt.figure(figsize=(10, 8))

        plt.axhline(0, color='grey', linestyle='--', linewidth=1, zorder=0)
        plt.axvline(0, color='grey', linestyle='--', linewidth=1, zorder=0)
        plt.plot([min_val, max_val], [min_val, max_val], color='black', linestyle='--', linewidth=1, zorder=0)
        scatter = plt.scatter(x, y, c=colors, alpha=0.7, s=10)

        plt.xlabel("Bi-variate Moran's I between \ngene expression and AB score in AD group", fontsize=16)
        plt.ylabel("Bi-variate Moran's I between \ngene expression and AB score in CT group", fontsize=16)
        plt.xlim(min_val * 1.03, max_val * 1.03)
        plt.ylim(min_val * 1.03, max_val * 1.03)
        plt.gca().set_aspect('equal', adjustable='box')

        texts = []
        for LR_cur_index in range(global_moranI_all_AD_CT_positive.shape[0]):
            cur_df = global_moranI_all_AD_CT_positive.iloc[LR_cur_index, :]
            x_val = float(cur_df['global_moranI_AD'])
            y_val = float(cur_df['global_moranI_CT'])
            type_cur_highlight = cur_df['type_CT']
            LR_cur = cur_df['Gene1'] + "-" + cur_df['Gene2']
            if LR_cur == "NCAM1-NCAM1":
                plt.scatter([x_val], [y_val], color=type2color['LR'], edgecolor='yellow', s=50, linewidth=3, zorder=5)
        adjust_text(texts, arrowprops=dict(arrowstyle="->", color='gray', lw=0.5))

        from matplotlib.patches import Patch

        legend_elements = [Patch(facecolor=type2color[tp], label=tp) for tp in unique_types]
        plt.legend(handles=legend_elements, title='Gene pair type', fontsize=12, title_fontsize=13, loc='best')

        filename = f"/Users/junjietang/Desktop/AD project/Xenium_data_V1/global_moranI_geneHic_scatter_condition_{metagene1}_{metagene2}_NCAM1-NCAM1-highlight.pdf"
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        ##just show the positive part
        global_moranI_all_AD_CT_positivesel = global_moranI_all_AD_CT
        global_moranI_all_AD_CT_positivesel['type_CT'].loc[~((global_moranI_all_AD_CT['global_moranI_AD'] >=0) * (global_moranI_all_AD_CT['global_moranI_CT'] >=0))] = "NPGP"
        # .loc[(global_moranI_all_AD_CT['global_moranI_AD'] >=0) * (global_moranI_all_AD_CT['global_moranI_CT'] >=0),:]
        plt.close()
        plt.rcParams.update({
            'axes.labelsize': 16,
            'xtick.labelsize': 14,
            'ytick.labelsize': 14,
            'text.color': 'black'
        })

        x = global_moranI_all_AD_CT_positivesel['global_moranI_AD']
        y = global_moranI_all_AD_CT_positivesel['global_moranI_CT']
        type_cur = global_moranI_all_AD_CT_positivesel['type_CT']
        unique_types = type_cur.unique()
        # palette = sns.color_palette("tab10", n_colors=len(unique_types))
        # type2color = {tp: palette[i] for i, tp in enumerate(unique_types)}
        # colors = type_cur.map(type2color)
        # custom_palette = ['#95E1D3', '#F38181']
        custom_palette = ['#EEEEEE','#3FC1C9', '#FC5185']
        type2color = {tp: custom_palette[i] for i, tp in enumerate(unique_types)}
        colors = type_cur.map(type2color)
        x_min_abs = abs(x.min())
        y_min_abs = abs(y.min())
        x_max_abs = abs(x.max())
        y_max_abs = abs(y.max())
        max_all = max(x_max_abs, y_max_abs, x_min_abs, y_min_abs)
        min_val = -max_all
        max_val = max_all
        plt.figure(figsize=(10, 8))
        plt.axhline(0, color='grey', linestyle='--', linewidth=1, zorder=0)
        plt.axvline(0, color='grey', linestyle='--', linewidth=1, zorder=0)
        plt.plot([min_val, max_val], [min_val, max_val], color='black', linestyle='--', linewidth=1, zorder=0)
        scatter = plt.scatter(x, y, c=colors, alpha=0.7, s=10)
        plt.xlabel("Bi-variate Moran's I between \ngene expression and AB score in AD group", fontsize=16)
        plt.ylabel("Bi-variate Moran's I between \ngene expression and AB score in CT group", fontsize=16)
        plt.xlim(min_val * 1.03, max_val * 1.03)
        plt.ylim(min_val * 1.03, max_val * 1.03)
        plt.gca().set_aspect('equal', adjustable='box')
        texts = []
        for LR_cur_index in range(global_moranI_all_AD_CT_positive.shape[0]):
            cur_df = global_moranI_all_AD_CT_positive.iloc[LR_cur_index, :]
            x_val = float(cur_df['global_moranI_AD'])
            y_val = float(cur_df['global_moranI_CT'])
            type_cur_highlight = cur_df['type_CT']
            LR_cur = cur_df['Gene1'] + "-" + cur_df['Gene2']
            if LR_cur == "NCAM1-NCAM1":
                plt.scatter([x_val], [y_val], color=type2color['LR'], edgecolor='yellow', s=50, linewidth=3, zorder=5)
        adjust_text(texts, arrowprops=dict(arrowstyle="->", color='gray', lw=0.5))

        from matplotlib.patches import Patch

        legend_elements = [Patch(facecolor=type2color[tp], label=tp) for tp in unique_types]
        plt.legend(handles=legend_elements, title='Gene pair type', fontsize=12, title_fontsize=13, loc='best')

        filename = f"/Users/junjietang/Desktop/AD project/Xenium_data_V1/global_moranI_geneHic_scatter_condition_{metagene1}_{metagene2}_NCAM1-NCAM1-highlight_positive.pdf"
        plt.savefig(filename, dpi=300, bbox_inches='tight')


        ##load the global_moranI_all
        global_moranI_geneHic = pd.read_csv(
            "/Users/junjietang/Desktop/AD project/Xenium_data_V1/global_moranI_all_" + str(metagene1) + "_" + str(
                metagene2) + "_genexpressiontoHic.csv")
        global_moranI_genegene = pd.read_csv(
            "/Users/junjietang/Desktop/AD project/Xenium_data_V1/global_moranI_all_" + str(metagene1) + "_" + str(
                metagene2) + "_genexpression.csv")
        global_moranI_HicHic = pd.read_csv(
            "/Users/junjietang/Desktop/AD project/Xenium_data_V1/global_moranI_Hic_all_" + str(metagene1) + "_" + str(
                metagene2) + "_ABcompartmentscore.csv")
        ##Difference
        global_moranI_genegene_diff_list = []
        for GP_index in range(topGP_merge.shape[0]):
            topGP_merge_cur = topGP_merge.iloc[GP_index, :]
            gene1 = topGP_merge_cur['gene1']
            gene2 = topGP_merge_cur['gene2']
            global_moranI_genegene_cur = global_moranI_genegene.loc[(global_moranI_genegene['Gene1'] == gene1) * (global_moranI_genegene['Gene2'] == gene2),:]
            ##
            global_moranI_genegene_cur_AD = global_moranI_genegene_cur.loc[global_moranI_genegene_cur['Group'] == "AD",:]
            global_moranI_genegene_cur_CT = global_moranI_genegene_cur.loc[global_moranI_genegene_cur['Group'] == "CT",:]
            global_moranI_genegene_cur_diff = float(global_moranI_genegene_cur_AD['global_moranI']) - float(global_moranI_genegene_cur_CT['global_moranI'])
            global_moranI_genegene_cur.columns = ["Gene1", "Gene2", "global_moranI_diff", "Group","type"]
            global_moranI_genegene_cur = global_moranI_genegene_cur.iloc[0, :]
            global_moranI_genegene_cur['global_moranI_diff'] = global_moranI_genegene_cur_diff
            global_moranI_genegene_diff_list.append(global_moranI_genegene_cur)
        global_moranI_genegene_diff = pd.DataFrame(global_moranI_genegene_diff_list)

    if if_gene_Hic_samegene:
        ##Gene-to-Gene spatial correlation
        global_moranI_top1_sample = {}
        # topgene_metagene_use = topgene_metagene1
        topgene_metagene_use = topgene_metagene2
        ##
        if if_modify_moran:
            topgene_assogene_mean_alltissue = np.mean(adata_Xenium[:, topgene_metagene_use].X, axis=0)
            topgene_assogene_std_alltissue = np.std(adata_Xenium[:, topgene_metagene_use].X.toarray(), axis=0)
            topgene_assogene_Hic_mean_alltissue = np.mean(compartment_subset.loc[:, topgene_metagene_use], axis=0)
            topgene_assogene_Hic_std_alltissue = np.std(compartment_subset.loc[:, topgene_metagene_use], axis=0)
            ##
            topgene_assogene_mean_alltissue = pd.DataFrame(topgene_assogene_mean_alltissue.T, index=topgene_metagene_use,
                                                            columns=["mean"])
            topgene_assogene_std_alltissue = pd.DataFrame(topgene_assogene_std_alltissue.T, index=topgene_metagene_use,
                                                           columns=["std"])
            topgene_assogene_Hic_mean_alltissue = pd.DataFrame(topgene_assogene_Hic_mean_alltissue.T, index=topgene_metagene_use,
                                                               columns=["mean"])
            topgene_assogene_Hic_std_alltissue = pd.DataFrame(topgene_assogene_Hic_std_alltissue.T, index=topgene_metagene_use,
                                                                columns=["std"])
        ##
        sample_unique = np.unique(adata_Xenium.obs['sample_id'])
        for sample_cur in sample_unique:
            print(sample_cur)
            cellindex_cursample = np.where(adata_Xenium.obs['sample_id'] == sample_cur)[0]
            geneexp_cursample = adata_Xenium.X[cellindex_cursample, :].toarray()
            geneexp_cursample = pd.DataFrame(geneexp_cursample,
                                             index=adata_Xenium.obs['orig.ident'].index[cellindex_cursample],
                                             columns=adata_Xenium.var_names)
            Hic_cursample = compartment_subset.iloc[cellindex_cursample, :]
            spatial_coord_cursample = adata_Xenium.obsm['spatial'][cellindex_cursample, :]
            ##
            spatial_coord_cursample_cur = spatial_coord_cursample
            ##
            w = libpysal.weights.KNN.from_array(spatial_coord_cursample_cur, k=num_neighbor)
            w.transform = 'r'
            ##Top genes
            global_moran_I_cursample = np.zeros((len(topgene_metagene_use), len(topgene_metagene_use)))
            global_moran_I_cursample = pd.DataFrame(global_moran_I_cursample,index = topgene_metagene_use,
                                                    columns=topgene_metagene_use)
            for topGP_merge_cur in topgene_metagene_use:
                ##Direction 1
                gene1 = topGP_merge_cur
                gene2 = topGP_merge_cur
                ##
                geneexp_cursample_cur = np.vstack([geneexp_cursample.loc[:,
                                        [gene1]].values[:,0],Hic_cursample.loc[:,gene2].values]).T
                ##
                moran_bv = Moran_BV_modified(x=geneexp_cursample_cur[:, 0], y=geneexp_cursample_cur[:, 1],
                                             w=w,
                                             x_mean=float(
                                                 topgene_assogene_mean_alltissue.loc[
                                                     gene1]),
                                             x_std=float(
                                                 topgene_assogene_std_alltissue.loc[
                                                     gene1]),
                                             y_mean=float(
                                                 topgene_assogene_Hic_mean_alltissue.loc[
                                                     gene2]),
                                             y_std=float(
                                                 topgene_assogene_Hic_std_alltissue.loc[
                                                     gene2]))
                ##
                global_moran_I_cursample.loc[gene1, gene2] = moran_bv.I
            global_moranI_top1_sample[sample_cur] = global_moran_I_cursample

        ##Condition-specific
        sample_unique_AD = ["AD-1", "AD-2"]
        sample_unique_CT = ["CT-1", "CT-2"]
        global_moranI_top1_AD = np.zeros((len(topgene_metagene_use), len(topgene_metagene_use)))
        for sample_cur in sample_unique_AD:
            global_moranI_top1_AD += global_moranI_top1_sample[sample_cur].values
        global_moranI_top1_AD /= len(sample_unique_AD)
        global_moranI_top1_AD = pd.DataFrame(global_moranI_top1_AD, index=topgene_metagene_use,
                                                 columns=topgene_metagene_use)
        ##Turn the global_moranI_top1 into a long dataframe
        global_moranI_top1_AD_long = global_moranI_top1_AD.stack().reset_index()
        global_moranI_top1_AD_long.columns = ["Gene1", "Gene2", "global_moranI"]
        global_moranI_top1_AD_long['Group'] = "Top1Top2_AD"
        ##
        global_moranI_top1_AD_long = global_moranI_top1_AD_long.loc[~(global_moranI_top1_AD_long['global_moranI'] == 0),:]
        ##
        global_moranI_top1_CT = np.zeros((len(topgene_metagene_use), len(topgene_metagene_use)))
        for sample_cur in sample_unique_CT:
            global_moranI_top1_CT += global_moranI_top1_sample[sample_cur].values
        global_moranI_top1_CT /= len(sample_unique_CT)
        global_moranI_top1_CT = pd.DataFrame(global_moranI_top1_CT, index=topgene_metagene_use,
                                                 columns=topgene_metagene_use)
        ##Turn the global_moranI_top1 into a long dataframe
        global_moranI_top1_CT_long = global_moranI_top1_CT.stack().reset_index()
        global_moranI_top1_CT_long.columns = ["Gene1", "Gene2", "global_moranI"]
        global_moranI_top1_CT_long['Group'] = "Top1Top2_CT"
        global_moranI_top1_CT_long = global_moranI_top1_CT_long.loc[~(global_moranI_top1_CT_long['global_moranI'] == 0),:]
        ##
        global_moranI_all = pd.concat([global_moranI_top1_AD_long, global_moranI_top1_CT_long], axis=0)
        ##
        choose_index = []
        for i in range(global_moranI_all.shape[0]):
            gene1 = global_moranI_all.iloc[i, 0]
            gene2 = global_moranI_all.iloc[i, 1]
            if np.sum((LR_Cellchat['ligand'] == gene2) * (LR_Cellchat['receptor'] == gene1)) > 0 or np.sum(
                    (LR_Cellchat['ligand'] == gene1) * (LR_Cellchat['receptor'] == gene2)) > 0:
                choose_index.append(i)
        # global_moranI_all = global_moranI_all.iloc[choose_index, :]
        global_moranI_all['type'] = "GP"
        global_moranI_all['type'].iloc[choose_index] = "LR"

        ##Keep the global_moranI_all['Group'] == "Top1Top2_AD" or "Top1Top2_CT"
        global_moranI_all = global_moranI_all[
            (global_moranI_all['Group'] == "Top1Top2_AD") | (global_moranI_all['Group'] == "Top1Top2_CT")]
        ##change the "Top1Top2_AD" to "AD" and "Top1Top2_CT" to "CT"
        global_moranI_all['Group'] = global_moranI_all['Group'].replace({"Top1Top2_AD": "AD", "Top1Top2_CT": "CT"})

        ##Show the density plot of global_moranI_all['global_moranI'], grouped by global_moranI_all['Group']
        import seaborn as sns
        import matplotlib.pyplot as plt
        import matplotlib
        ##
        matplotlib.use('TkAgg')
        plt.rcParams.update({
            'axes.titlesize': 20,
            'axes.labelsize': 18,
            'xtick.labelsize': 16,
            'ytick.labelsize': 16,
            'legend.title_fontsize': 16,
            'legend.fontsize': 16
        })

        plt.figure(figsize=(7, 8))
        ax = sns.kdeplot(
            data=global_moranI_all,
            x='global_moranI',
            hue='Group',
            fill=True,
            common_norm=False
        )

        # ax.set_xlabel("Bi-variate Moran's I of gene expression \nbetween associated ligand-receptor pair", fontsize=18)
        ax.set_xlabel("Bi-variate Moran's I of gene expression and AB score \nbetween associated gene pairs", fontsize=18)
        ax.set_ylabel("Density", fontsize=18)
        leg = ax.get_legend()
        leg.set_title("Group")
        plt.setp(leg.get_title(), fontsize=22)
        plt.setp(leg.get_texts(), fontsize=22)

        plt.tight_layout()

        plt.savefig("/Users/junjietang/Desktop/AD project/Xenium_data_V1/global_moranI_geneHicsamegene_density_condition_" + str(
            metagene1) + "_" + str(metagene2) + ".pdf", dpi=300, bbox_inches='tight')
        ##
        global_moranI_all_AD = global_moranI_all[global_moranI_all['Group'] == "AD"]
        global_moranI_all_CT = global_moranI_all[global_moranI_all['Group'] == "CT"]

        ##Show a scatter plot bettwen global_moranI_all_AD['global_moranI'] and global_moranI_all_CT['global_moranI']
        ##Merge global_moranI_all_AD and global_moranI_all_CT
        global_moranI_all_AD_CT = pd.merge(global_moranI_all_AD, global_moranI_all_CT, on=['Gene1', 'Gene2'],
                                           suffixes=('_AD', '_CT'))
        ##
        ##
        global_moranI_all_AD_CT_positive = global_moranI_all_AD_CT.loc[~(
                    (global_moranI_all_AD_CT['global_moranI_AD'] < 0) * (
                        global_moranI_all_AD_CT['global_moranI_CT'] < 0)), :]
        global_moranI_all_AD_CT_positive = global_moranI_all_AD_CT_positive.iloc[np.argsort(
            np.max(global_moranI_all_AD_CT_positive[['global_moranI_AD', 'global_moranI_CT']], axis=1))[
                                                                                 ::-1][0:20], :]
        ##
        ##
        import matplotlib.pyplot as plt
        from adjustText import adjust_text
        import seaborn as sns
        from matplotlib.cm import get_cmap

        plt.rcParams.update({
            'axes.labelsize': 16,
            'xtick.labelsize': 14,
            'ytick.labelsize': 14,
            'text.color': 'black'
        })

        x = global_moranI_all_AD_CT['global_moranI_AD']
        y = global_moranI_all_AD_CT['global_moranI_CT']
        type_cur = global_moranI_all_AD_CT['type_CT']

        unique_types = type_cur.unique()
        palette = sns.color_palette("tab10", n_colors=len(unique_types))
        type2color = {tp: palette[i] for i, tp in enumerate(unique_types)}
        colors = type_cur.map(type2color)
        x_min_abs = abs(x.min())
        y_min_abs = abs(y.min())
        x_max_abs = abs(x.max())
        y_max_abs = abs(y.max())
        max_all = max(x_max_abs, y_max_abs, x_min_abs, y_min_abs)
        min_val = -max_all
        max_val = max_all

        plt.figure(figsize=(10, 8))

        scatter = plt.scatter(x, y, c=colors, alpha=0.7, s=5)

        plt.plot([min_val, max_val], [min_val, max_val], color='red', linestyle='--', linewidth=1)

        plt.xlabel("Bi-variate Moran's I between \ngene expression and AB score in AD group", fontsize=16)
        plt.ylabel("Bi-variate Moran's I between \ngene expression and AB score in CT group", fontsize=16)
        plt.xlim(min_val * 1.1, max_val * 1.1)
        plt.ylim(min_val * 1.1, max_val * 1.1)
        plt.gca().set_aspect('equal', adjustable='box')

        texts = []
        for LR_cur_index in range(global_moranI_all_AD_CT_positive.shape[0]):
            cur_df = global_moranI_all_AD_CT_positive.iloc[LR_cur_index, :]
            x_val = float(cur_df['global_moranI_AD'])
            y_val = float(cur_df['global_moranI_CT'])
            type_cur_highlight = cur_df['type_CT']
            LR_cur = cur_df['Gene1'] + "-" + cur_df['Gene2']
            if type_cur_highlight == "LR":
                texts.append(
                    plt.text(x_val, y_val, LR_cur,
                             fontsize=14, color="red", ha='left', va='bottom',
                             bbox=dict(boxstyle="round,pad=0.3", edgecolor="gray", facecolor="white", alpha=0.7))
                )
            else:
                texts.append(
                    plt.text(x_val, y_val, LR_cur,
                             fontsize=14, color="black", ha='left', va='bottom',
                             bbox=dict(boxstyle="round,pad=0.3", edgecolor="gray", facecolor="white", alpha=0.7))
                )
        adjust_text(texts, arrowprops=dict(arrowstyle="->", color='gray', lw=0.5))

        plt.axhline(0, color='black', linestyle='--', linewidth=1)
        plt.axvline(0, color='black', linestyle='--', linewidth=1)

        from matplotlib.patches import Patch

        legend_elements = [Patch(facecolor=type2color[tp], label=tp) for tp in unique_types]
        plt.legend(handles=legend_elements, title='Cell type', fontsize=12, title_fontsize=13, loc='best')

        filename = f"/Users/junjietang/Desktop/AD project/Xenium_data_V1/global_moranI_geneHicsamegene_scatter_condition_{metagene1}_{metagene2}.pdf"
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        # plt.show()
        ##just highlight the NCAM1-NCAM1
        plt.close()
        plt.rcParams.update({
            'axes.labelsize': 16,
            'xtick.labelsize': 14,
            'ytick.labelsize': 14,
            'text.color': 'black'
        })

        x = global_moranI_all_AD_CT['global_moranI_AD']
        y = global_moranI_all_AD_CT['global_moranI_CT']
        type_cur = global_moranI_all_AD_CT['type_CT']
        unique_types = type_cur.unique()
        # palette = sns.color_palette("tab10", n_colors=len(unique_types))
        # type2color = {tp: palette[i] for i, tp in enumerate(unique_types)}
        # colors = type_cur.map(type2color)
        # custom_palette = ['#95E1D3', '#F38181']
        custom_palette = ['#3FC1C9', '#FC5185']
        type2color = {tp: custom_palette[i] for i, tp in enumerate(unique_types)}
        colors = type_cur.map(type2color)
        x_min_abs = abs(x.min())
        y_min_abs = abs(y.min())
        x_max_abs = abs(x.max())
        y_max_abs = abs(y.max())
        max_all = max(x_max_abs, y_max_abs, x_min_abs, y_min_abs)
        min_val = -max_all
        max_val = max_all
        plt.figure(figsize=(10, 8))
        plt.axhline(0, color='grey', linestyle='--', linewidth=1, zorder=0)
        plt.axvline(0, color='grey', linestyle='--', linewidth=1, zorder=0)
        plt.plot([min_val, max_val], [min_val, max_val], color='black', linestyle='--', linewidth=1, zorder=0)
        scatter = plt.scatter(x, y, c=colors, alpha=0.7, s=10)

        plt.xlabel("Bi-variate Moran's I between \ngene expression and AB score in AD group", fontsize=16)
        plt.ylabel("Bi-variate Moran's I between \ngene expression and AB score in CT group", fontsize=16)
        plt.xlim(min_val * 1.03, max_val * 1.03)
        plt.ylim(min_val * 1.03, max_val * 1.03)
        plt.gca().set_aspect('equal', adjustable='box')
        texts = []
        for LR_cur_index in range(global_moranI_all_AD_CT_positive.shape[0]):
            cur_df = global_moranI_all_AD_CT_positive.iloc[LR_cur_index, :]
            x_val = float(cur_df['global_moranI_AD'])
            y_val = float(cur_df['global_moranI_CT'])
            type_cur_highlight = cur_df['type_CT']
            LR_cur = cur_df['Gene1'] + "-" + cur_df['Gene2']
            if LR_cur == "NCAM1-NCAM1":
                plt.scatter([x_val], [y_val], color=type2color['LR'], edgecolor='yellow', s=50, linewidth=3, zorder=5)
        adjust_text(texts, arrowprops=dict(arrowstyle="->", color='gray', lw=0.5))

        from matplotlib.patches import Patch

        legend_elements = [Patch(facecolor=type2color[tp], label=tp) for tp in unique_types]
        plt.legend(handles=legend_elements, title='Gene pair type', fontsize=12, title_fontsize=13, loc='best')

        filename = f"/Users/junjietang/Desktop/AD project/Xenium_data_V1/global_moranI_geneHicsamegene_scatter_condition_{metagene1}_{metagene2}_NCAM1-NCAM1-highlight.pdf"
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        ##just show the positive part
        global_moranI_all_AD_CT_positivesel = global_moranI_all_AD_CT
        global_moranI_all_AD_CT_positivesel['type_CT'].loc[~((global_moranI_all_AD_CT['global_moranI_AD'] >=0) * (global_moranI_all_AD_CT['global_moranI_CT'] >=0))] = "NPGP"
        # .loc[(global_moranI_all_AD_CT['global_moranI_AD'] >=0) * (global_moranI_all_AD_CT['global_moranI_CT'] >=0),:]
        plt.close()
        plt.rcParams.update({
            'axes.labelsize': 16,
            'xtick.labelsize': 14,
            'ytick.labelsize': 14,
            'text.color': 'black'
        })

        x = global_moranI_all_AD_CT_positivesel['global_moranI_AD']
        y = global_moranI_all_AD_CT_positivesel['global_moranI_CT']
        type_cur = global_moranI_all_AD_CT_positivesel['type_CT']
        unique_types = type_cur.unique()
        # palette = sns.color_palette("tab10", n_colors=len(unique_types))
        # type2color = {tp: palette[i] for i, tp in enumerate(unique_types)}
        # colors = type_cur.map(type2color)
        # custom_palette = ['#95E1D3', '#F38181']
        custom_palette = ['#3FC1C9', '#EEEEEE', '#FC5185']
        type2color = {tp: custom_palette[i] for i, tp in enumerate(unique_types)}
        colors = type_cur.map(type2color)
        x_min_abs = abs(x.min())
        y_min_abs = abs(y.min())
        x_max_abs = abs(x.max())
        y_max_abs = abs(y.max())
        max_all = max(x_max_abs, y_max_abs, x_min_abs, y_min_abs)
        min_val = -max_all
        max_val = max_all
        plt.figure(figsize=(10, 8))
        plt.axhline(0, color='grey', linestyle='--', linewidth=1, zorder=0)
        plt.axvline(0, color='grey', linestyle='--', linewidth=1, zorder=0)
        plt.plot([min_val, max_val], [min_val, max_val], color='black', linestyle='--', linewidth=1, zorder=0)
        scatter = plt.scatter(x, y, c=colors, alpha=0.7, s=10)

        plt.xlabel("Bi-variate Moran's I between \ngene expression and AB score in AD group", fontsize=16)
        plt.ylabel("Bi-variate Moran's I between \ngene expression and AB score in CT group", fontsize=16)
        plt.xlim(min_val * 1.03, max_val * 1.03)
        plt.ylim(min_val * 1.03, max_val * 1.03)
        plt.gca().set_aspect('equal', adjustable='box')

        texts = []
        for LR_cur_index in range(global_moranI_all_AD_CT_positive.shape[0]):
            cur_df = global_moranI_all_AD_CT_positive.iloc[LR_cur_index, :]
            x_val = float(cur_df['global_moranI_AD'])
            y_val = float(cur_df['global_moranI_CT'])
            type_cur_highlight = cur_df['type_CT']
            LR_cur = cur_df['Gene1'] + "-" + cur_df['Gene2']
            if LR_cur == "NCAM1-NCAM1":
                plt.scatter([x_val], [y_val], color=type2color['LR'], edgecolor='yellow', s=50, linewidth=3, zorder=5)
        adjust_text(texts, arrowprops=dict(arrowstyle="->", color='gray', lw=0.5))

        from matplotlib.patches import Patch

        legend_elements = [Patch(facecolor=type2color[tp], label=tp) for tp in unique_types]
        plt.legend(handles=legend_elements, title='Gene pair type', fontsize=12, title_fontsize=13, loc='best')

        filename = f"/Users/junjietang/Desktop/AD project/Xenium_data_V1/global_moranI_geneHicsamegene_scatter_condition_{metagene1}_{metagene2}_NCAM1-NCAM1-highlight_positive.pdf"
        plt.savefig(filename, dpi=300, bbox_inches='tight')


        ##load the global_moranI_all
        global_moranI_geneHic = pd.read_csv(
            "/Users/junjietang/Desktop/AD project/Xenium_data_V1/global_moranI_all_" + str(metagene1) + "_" + str(
                metagene2) + "_genexpressiontoHic.csv")
        global_moranI_genegene = pd.read_csv(
            "/Users/junjietang/Desktop/AD project/Xenium_data_V1/global_moranI_all_" + str(metagene1) + "_" + str(
                metagene2) + "_genexpression.csv")
        global_moranI_HicHic = pd.read_csv(
            "/Users/junjietang/Desktop/AD project/Xenium_data_V1/global_moranI_Hic_all_" + str(metagene1) + "_" + str(
                metagene2) + "_ABcompartmentscore.csv")
        ##Difference
        global_moranI_genegene_diff_list = []
        for GP_index in range(topGP_merge.shape[0]):
            topGP_merge_cur = topGP_merge.iloc[GP_index, :]
            gene1 = topGP_merge_cur['gene1']
            gene2 = topGP_merge_cur['gene2']
            global_moranI_genegene_cur = global_moranI_genegene.loc[(global_moranI_genegene['Gene1'] == gene1) * (global_moranI_genegene['Gene2'] == gene2),:]
            ##
            global_moranI_genegene_cur_AD = global_moranI_genegene_cur.loc[global_moranI_genegene_cur['Group'] == "AD",:]
            global_moranI_genegene_cur_CT = global_moranI_genegene_cur.loc[global_moranI_genegene_cur['Group'] == "CT",:]
            global_moranI_genegene_cur_diff = float(global_moranI_genegene_cur_AD['global_moranI']) - float(global_moranI_genegene_cur_CT['global_moranI'])
            global_moranI_genegene_cur.columns = ["Gene1", "Gene2", "global_moranI_diff", "Group","type"]
            global_moranI_genegene_cur = global_moranI_genegene_cur.iloc[0, :]
            global_moranI_genegene_cur['global_moranI_diff'] = global_moranI_genegene_cur_diff
            global_moranI_genegene_diff_list.append(global_moranI_genegene_cur)
        global_moranI_genegene_diff = pd.DataFrame(global_moranI_genegene_diff_list)