
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from itertools import permutations
from itertools import combinations
from tqdm import tqdm
from tqdm.notebook import tqdm
from functools import reduce
import trimesh
from shapely import points, contains



from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, WhiteKernel, ConstantKernel as C



from scipy.spatial import Delaunay

import numpy as np
import numpy as np
import matplotlib.pyplot as plt
import alphashape
from shapely.geometry import Point
import joblib



def gp_training_application(X, Y, y_name, X_emu, path = "/glade/work/qingyuany/repo_data/spatialtuning/", no_sen = 2, no_restart = 10):
    y = Y[y_name]
    y_norm = (y - y.mean())/y.std()

    if ~(np.isnan(y_norm).all()):
        sage1d_temp = fit_all_gp_models_1d(y_norm.values, X.values, len1d = 0.4)
        print(sage1d_temp)
        sel_para_ind = sage1d_temp[:(no_sen * 2):2].astype(int)

        #sel_para_ind = sel_para_ind.reshape(1,-1)

        kernel = C(1.0, (1e-3, 1)) * Matern(length_scale=0.001, nu=2.5, length_scale_bounds=(0.1, 3)) + WhiteKernel(noise_level=0.5, noise_level_bounds=(0.01, 0.9))
        gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer= no_restart, normalize_y=True)
        
        gp.fit(X.values[:,sel_para_ind], y_norm.values)
        
        joblib.dump(gp, path + "python_obj/" + y_name + "_gpmodel.pkl")

        y_mean_emu, y_std_emu = gp.predict(X_emu.values[:,sel_para_ind], return_std=True)

        pd.DataFrame(y_mean_emu, columns=[y_name]).to_csv(path + "emulated_ensemble" + "/gp_mean_" + y_name + ".csv")
        pd.DataFrame(y_std_emu, columns=[y_name]).to_csv(path + "emulated_ensemble" + "/gp_std_" + y_name + ".csv")

        return ([y_name, sel_para_ind])

        

def fit_gp_for_single_1d(y, X, i, length_scale=0.5):
    """Fits a GP model for a given pair of feature indices."""
    
    X_single = X[:, [i]].reshape([-1,1])  # Select the two features
    
    # Define GP kernel
    kernel =  Matern(length_scale=length_scale, 
                     length_scale_bounds="fixed", nu = 2.5)  + WhiteKernel(noise_level = 1.0)
    
    # Fit GP model
    gp = GaussianProcessRegressor(kernel=kernel, optimizer=None)
    gp.fit(X_single, y.ravel())  # Flatten y

    y_pred = gp.predict(X_single, return_std = False)
    residual_sd = (y - y_pred).std()
    
    return (i, residual_sd)



def fit_all_gp_models_1d(y, X, len1d = 0.5):
    """Fits GP models for all unique feature pairs and stores residual std."""
    residual_single_stds = {}
    
    y_std = y.std()
    
    feature_single = itertools.combinations(range(X.shape[1]), 1)
    

    for i in feature_single:
        residual_single_stds[(i)] = y_std - fit_gp_for_single_1d(y, X, i, len1d)[1]  # Store only the residual std

    
    sorted_delta = sorted(residual_single_stds.items(), key=lambda x: x[1], reverse=True)
    top_10 = dict(sorted_delta[:20])


    output = np.array([k + (v,) for k, v in top_10.items()]).flatten()
    
    return output



def compute_mi_pair(col_x, col_y):
    """
    Calculate the Mutual Information between two d series
    xxTest whether it is used laterxx
    """
    
    x = col_x.values.reshape(-1, 1)
    y = col_y.values
    mi = mutual_info_regression(x, y, discrete_features=False)
    return mi[0]



def dist_diff(df1, df2, bins=100):
    '''
    A simple way to calculate the "distribution difference" for each column of the pd dataframes using np.histogram
    Inputs:
        df1: the surviving ensemble members for each climatology (pd dataframes).
        df2: normally the original randomly sampled ensemble members (pd dataframes).
    Used to determine what parameters are varied after constraining by certain parameters
    '''
    abs_den_diff = {}
    for column in df1.columns:
        hist1, bin_edges1 = np.histogram(df1[column], bins=bins, density=True, range = (0,1))
        hist2, bin_edges2 = np.histogram(df2[column], bins=bins, density=True, range = (0,1))
    # Compute bin widths and probabilities
        abs_den_diff[column] = np.sum(abs(hist1 - hist2))

    output = pd.Series(abs_den_diff)

    return pd.Series(output)




# def filter_per_critical_para(tf_masks, sel_local_var, no_sample = 3):
#     ### ?? xxxx ####
#     output_pd = []
#     index_names = np.arange(no_sample).astype(str).tolist() + ["count"]
#     for vars1 in list(combinations(sel_local_var, no_sample)):
#         sel_vars = list(vars1)
#         sel_vars.append(tf_masks[sel_vars].all(axis = 1).sum())
#         output_pd.append(pd.Series(sel_vars, index = index_names))

#     output_pd = pd.concat(output_pd, axis = 1).transpose()
#     output_pd = output_pd.sort_values("count")
#     return output_pd




def get_top_2_info(row, n_para):
    top2 = row.nlargest(n_para)
    if n_para == 3:
        return pd.Series({
            "localvar": row.name,
            'para1': top2.index[0],
            'diff1': top2.iloc[0],
            'para2': top2.index[1],
            'diff2': top2.iloc[1],
            'para3': top2.index[2],
            'diff3': top2.iloc[2]
        })
    if n_para == 2:
        return pd.Series({
            "localvar": row.name,
            'para1': top2.index[0],
            'diff1': top2.iloc[0],
            'para2': top2.index[1],
            'diff2': top2.iloc[1]
        })




def group_para_climatology(df, tf_masks, meta, vars = np.NAN, threshold = 5000):
    
    if isinstance(vars, list):
       vars = vars
    else:
        vars = list(tf_masks.columns)


    n_para = meta.shape[0]

    paras_vars = {}

    for c in vars: 
        para_inds = meta[c]
        para_nms = sorted(list(df.columns[para_inds]))
        paras_vars.setdefault(tuple(para_nms), []).append(c)
        

    strt_paras_vars = {}
    surv_paras_vars = {}
    for k, v in paras_vars.items():

        if tf_masks[v].all(axis = 1).sum() < threshold :
            strt_paras_vars[k] = v
        else:
            surv_paras_vars[k] = v

    return paras_vars, strt_paras_vars, surv_paras_vars




def strterr_detection(paras_vars, tf_masks, n_comb = 2):
    summary_table = {}

    for para3s, vars in tqdm(list(paras_vars.items())):
        
        vars_temp = vars 
        pd_list = []
        
        for vars_comb in combinations(vars_temp, n_comb):
            pd_list.append(list(vars_comb) + [tf_masks[list(vars_comb)].all(axis = 1).sum()])
            
        pd_list = pd.DataFrame(pd_list, columns = [f"var{i+1}" for i in range(n_comb)] + ["count"])
        pd_list = pd_list.sort_values(by = "count")
        
        
        summary_table[para3s] = pd_list

    return summary_table





def range_err_detection(paras_vars, tf_masks, emu_para, meta):

    para_inds = np.sort(pd.unique(meta.values.ravel()))
    para_nm = list(emu_para.columns)
    
    para1_vars = {}
    for p in para_inds:
        for c in tf_masks.columns:
            if p in meta[c].values:
                para1_vars.setdefault(para_nm[p], []).append(c)

    para1_min_max = {}
    
    
    for p1, vars in list(para1_vars.items()):
        pts_list = []
        for var in vars:    
            temp_pts = emu_para[tf_masks[var]]
            if temp_pts.shape[0] > 40000:
                temp_pts = temp_pts.sample(40000)
            pts_list.append(temp_pts[p1].values)
            
        temp_min  = np.array([lst.min() for lst in pts_list])
        temp_max  = np.array([lst.max() for lst in pts_list])
        para1_min_max[p1] = [temp_min.max(), temp_max.min()]

    
    return para1_vars, para1_min_max


## XX Nov 6
def para1_error_localvar_detection(para1_vars_dict, p, emu_para, tf_masks, n_comb = 1):

    p = [p]
    vars = para1_vars_dict[tuple(p)]
    

    vars_min_max = []
    print(len(vars))
    for var_comb in list(combinations(vars, n_comb)):
        used_vars = [v for v in vars if v not in var_comb]
        
        pts_list = [emu_para[tf_masks[v].all(axis = 1)][p] for v in used_vars]
        pts_max = min([pts.max().values for pts in pts_list])
        pts_min = max([pts.min().values for pts in pts_list])

        vars_min_max.append(pd.Series([var_comb, pts_min, pts_max]))

    output = pd.DataFrame(vars_min_max)
    output.columns = ["excluded_vars", "min", "max"]
    output_f = output[output["min"] > output["max"]]
    return output, output_f
    

#xxx Oct 31st

def para2_error_detection(paras_vars, tf_masks, emu_para, shape_alpha = 7):
    
    paras = list(set([x for tup in list(paras_vars.keys()) for x in tup]))
    paras3 = list(paras_vars.keys())

    para2_para3 = {}
    
    for pair in combinations(paras, 2):
        pair = set(pair)
        temp_list = [t for t in paras3 if pair.issubset(t)]


        if len(temp_list) > 0:
            para2_para3[tuple(pair)] = temp_list

    pairs_hulls = {}
    
    para2_vars = {}   ## ?? xx

    for para2, para3 in tqdm(para2_para3.items()):
        pts_list = []
        vars_list = []
        for p3 in para3:
            vars_temp = paras_vars[p3]
            vars_list.extend(vars_temp)
            temp_pts = emu_para[tf_masks[vars_temp].all(axis = 1)]
            print(temp_pts.shape)
            if temp_pts.shape[0] > 20000:
                temp_pts = temp_pts.sample(20000)
            pts_list.append(temp_pts[list(para2)].values)

        hulls = [alphashape.alphashape(points, shape_alpha) for points in pts_list]
        
        intersection_hulls = reduce(lambda a, b: a.intersection(b), hulls)
        pairs_hulls[para2] = intersection_hulls
        
        para2_vars[para2] = vars_list

    print([v for v, k in list(pairs_hulls.items()) if k.is_empty])

    
    return para2_vars, pairs_hulls





def para2_error_localvar_detection(para2_localvars_dict, para2_interest, emu_para,tf_masks, no_comb = 2,  shape_alpha = 7):
    ## xx??
    para2_interest = tuple(para2_interest)
    localvars_interest = para2_localvars_dict[para2_interest]
    
    exclude_hulls ={}

    for masked_localvars in tqdm(combinations(localvars_interest, no_comb)):
        masked_localvars = list(masked_localvars)
        analyzed_localvars = [localvar for localvar in localvars_interest if localvar not in masked_localvars]

        sub_samples = emu_para[tf_masks[analyzed_localvars].all(axis = 1)][list(para2_interest)]
        if sub_samples.shape[0] > 40000:
            sub_samples = sub_samples.sample(40000)

        exclude_hulls[tuple(masked_localvars)] = alphashape.alphashape(sub_samples.values, shape_alpha)
        
    
    print([v for v, k in exclude_hulls.items() if not k.is_empty])

    return exclude_hulls

############################################################
def para3_mesh_generator(surv_para3_localvar, emu_para, tf_masks):
    para3_meshes = {}
    for k, v in tqdm(surv_para3_localvar.items()):
        temp_surv_pts = emu_para[tf_masks[v].all(axis = 1)]
        if temp_surv_pts.shape[0] > 10000:
            temp_surv_pts = temp_surv_pts.sample(10000)
            
        temp_pts_np = temp_surv_pts[list(k)].values
        temp_pts_index = alpha_shape_3D(temp_pts_np, 7)
    
        para3_meshes[k] = trimesh.Trimesh(vertices=temp_pts_np, faces=list(temp_pts_index))

    
    return para3_meshes



def sample_from_poly(no_pts, para_nm, single_min_max_range, poly_dict, existing_pts = np.NAN, monitor = True):
    ## Not thoroughly checked
    if not isinstance(existing_pts, pd.DataFrame):
        if np.isnan(existing_pts):
            random_pts = pd.DataFrame(np.random.uniform(0, 1, (no_pts, len(para_nm))), columns = para_nm)
        for k, v in single_min_max_range.items():
            random_pts[list(k)] = random_pts[list(k)] * (v[1] - v[0]) + v[0]
        
    else:
        random_pts = existing_pts
    

    print(random_pts.shape)
    for k, polyg in poly_dict.items():
        temp_pts = random_pts[list(k)].values
        temp_pts = points(temp_pts)
        
        temp_masks = contains(polyg,temp_pts)
        if monitor:
            print(f'{k}: {temp_masks.sum()}')
        if temp_masks.sum() == 0:
            print("Didn't find anything")
            return []

        random_pts = random_pts[temp_masks]

    return random_pts


def sample_from_mesh(no_pts, para_nm, single_min_max_range, mesh_dict, existing_pts = np.NAN, monitor = True):
    if not isinstance(existing_pts, pd.DataFrame):
        if np.isnan(existing_pts):
            random_pts = pd.DataFrame(np.random.uniform(0, 1, (no_pts, len(para_nm))), columns = para_nm)
        for k, v in single_min_max_range.items():
            random_pts[list(k)] = random_pts[list(k)] * (v[1] - v[0]) + v[0]
        
    else:
        random_pts = existing_pts
    
    for k, mesh in tqdm(mesh_dict.items()):
        temp_pts = random_pts[list(k)].values
        
        temp_masks =  mesh.contains(temp_pts)
        if monitor:
            print(f'{k}: {temp_masks.sum()}')
        if temp_masks.sum() == 0:
            print("Didn't find anything")
            return []

        random_pts = random_pts[temp_masks]

    return random_pts
########################################################






def dict_pairs_to_df(d, dedupe=False, sort=False):
    """
    Build a DataFrame of all 2-key combinations from dict `d`.
    
    Parameters
    ----------
    d : dict
        Keys: tuples (or any hashable).
        Values: iterables (typically lists) of elements to combine.
    dedupe : bool, default False
        If True, remove duplicate elements when combining value lists (order-preserving).
    sort : bool, default False
        If True (applied after optional dedupe), sort the combined value list.
    
    Returns
    -------
    pd.DataFrame with columns ['key1', 'key2', 'values'].
    """
    rows = []
    for k1, k2 in combinations(d.keys(), 2):
        v1 = list(d[k1])  # ensure list copy
        v2 = list(d[k2])
        combined = v1 + v2

        if dedupe:
            seen = set()
            deduped = []
            for x in combined:
                if x not in seen:
                    seen.add(x)
                    deduped.append(x)
            combined = deduped

        if sort:
            combined = sorted(combined)

        rows.append((k1, k2, combined))

    return pd.DataFrame(rows, columns=['key1', 'key2', 'localvars'])







import numpy as np
import alphashape
from shapely.geometry import Point
from scipy.spatial import Delaunay
#import plotly.graph_objects as go
import alphashape
import itertools

import itertools

def alpha_shape_3D(points, alpha):
    
    # def add_edge(edges, i, j, k):
    #     #perms = (set(p) for p in itertools.permutations([i, j, k]))
    #     #print(perms)
    #     #common_sets = perms.intersection(edges)
    #     #exists = len(common_sets) > 0
        
    #     if (i, j, k) in edges or (i, k, j) in edges or (j, i, k) in edges or (j, k, i) in edges or (k, i, j) in edges or (k, j, i) in edges:
    #         edges.remove((i, k, j))
    #         edges.remove((j, i, k))
    #         edges.remove((j, k, i))
    #         edges.remove((k, i, j))
    #         edges.remove((k, j, i))
    #         return
            
        # edges.add((i, j, k))

    def add_edge(edges, i, j, k):
        edge = tuple(sorted((i, j, k)))
        if edge in edges:    
            edges.remove(edge)
            return
        edges.add(edge)

# #######

    tetra = Delaunay(points)
    # Extract tetrahedra
    edges = set()
    for ia, ib, ic, id in tetra.simplices:
        pa = points[ia]
        pb = points[ib]
        pc = points[ic]
        pd = points[id]
        # Compute radius of circumsphere
        circumsphere_radius = circum_radius_calculation(pa, pb, pc, pd)
        if circumsphere_radius < 1.0 / alpha:
            add_edge(edges, ia, ib, ic)
            add_edge(edges, ia, ib, id)
            add_edge(edges, ia, ic, id)
            add_edge(edges, ib, ic, id)
            
    return edges




def dist_cal(pt1, pt2):
    dist = (np.sum((pt1 - pt2) **2))**0.5
    return dist




def circum_radius_calculation(pta, ptb, ptc, ptd):
    vol_mat = np.vstack([pta, ptb, ptc, ptd]).transpose()
    vol_mat = np.vstack([vol_mat, np.array([1,1,1,1])])
    det_M = np.linalg.det(vol_mat)
    V = np.abs(det_M)/6

    a = dist_cal(ptd, pta)
    b = dist_cal(ptd, ptb)
    c = dist_cal(ptd, ptc)

    A = dist_cal(ptb, ptc)
    B = dist_cal(pta, ptc)
    C = dist_cal(pta, ptb)

    q_medium = (a * A + b * B + c * C) * (a * A + b * B - c * C) * (a * A - b * B + c * C) * (-a * A + b * B + c * C)
    radius = q_medium ** 0.5 / (24 * V)
    
    return radius




def circum_center_calculation(pta, ptb, ptc, ptd):
    A = np.vstack([ptb - pta, ptc - pta, ptd - pta])
    B = np.array([np.sum(ptb ** 2) - np.sum(pta ** 2), 
              np.sum(ptc ** 2) - np.sum(pta ** 2), 
              np.sum(ptd ** 2) - np.sum(pta ** 2)]).reshape(3,1)
    B = B * 0.5
    circuncenter_xyz = np.linalg.inv(A) @ B
    return circuncenter_xyz.flatten()


def alpha_shape2d(points, alpha, only_outer=True):
    """
    Compute the alpha shape (concave hull) of a set of points.
    :param points: np.array of shape (n,2) points.
    :param alpha: alpha value.
    :param only_outer: boolean value to specify if we keep only the outer border
    or also inner edges.
    :return: set of (i,j) pairs representing edges of the alpha-shape. (i,j) are
    the indices in the points array.
    """
    assert points.shape[0] > 3, "Need at least four points"

    # def add_edge(edges, i, j):
    #     """
    #     Add an edge between the i-th and j-th points,
    #     if not in the list already
    #     """
    #     if (i, j) in edges or (j, i) in edges:
    #         # already added
    #         assert (j, i) in edges, "Can't go twice over same directed edge right?"
    #         if only_outer:
    #             # if both neighboring triangles are in shape, it's not a boundary edge
    #             edges.remove((j, i))
    #         return
    #     edges.add((i, j))

    def add_edge(edges, i, j):
        edge = tuple(sorted((i, j)))
        if edge in edges:
        
            edges.remove(edge)
            return
        edges.add(edge)

    
    tri = Delaunay(points)
    edges = set()
    # Loop over triangles:
    # ia, ib, ic = indices of corner points of the triangle
    for ia, ib, ic in tri.simplices:
        pa = points[ia]
        pb = points[ib]
        pc = points[ic]
        # Computing radius of triangle circumcircle
        # www.mathalino.com/reviewer/derivation-of-formulas/derivation-of-formula-for-radius-of-circumcircle
        a = np.sqrt((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2)
        b = np.sqrt((pb[0] - pc[0]) ** 2 + (pb[1] - pc[1]) ** 2)
        c = np.sqrt((pc[0] - pa[0]) ** 2 + (pc[1] - pa[1]) ** 2)
        s = (a + b + c) / 2.0
        area = np.sqrt(s * (s - a) * (s - b) * (s - c))
        circum_r = a * b * c / (4.0 * area)
        if circum_r < alpha:
            add_edge(edges, ia, ib)
            add_edge(edges, ib, ic)
            add_edge(edges, ic, ia)
    return edges

