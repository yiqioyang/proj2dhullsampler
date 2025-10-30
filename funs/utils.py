
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





from scipy.spatial import Delaunay

import numpy as np
import numpy as np
import matplotlib.pyplot as plt
import alphashape
from shapely.geometry import Point



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



def distribution_difference(df1, df2, bins=100):
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




def get_top_2_info(row, n_large):
    top2 = row.nlargest(n_large)
    if n_large == 3:
        return pd.Series({
            "localvar": row.name,
            'topvpara1': top2.index[0],
            'diff1': top2.iloc[0],
            'toppara2': top2.index[1],
            'diff2': top2.iloc[1],
            'toppara3': top2.index[2],
            'diff3': top2.iloc[2]
        })
    if n_large == 2:
        return pd.Series({
            "localvar": row.name,
            'topvpara1': top2.index[0],
            'diff1': top2.iloc[0],
            'toppara2': top2.index[1],
            'diff2': top2.iloc[1]
        })




def column_with_smallest_std(df, tf_masks, vars = np.NAN, group3_threshold = 5000, n_large = 3):
    
    if isinstance(vars, list):
       vars = vars
    else:
        vars = list(tf_masks.columns)
    
    output_ddiff = []
    output_info = []
    for var in tqdm(vars):
        df_s = df[tf_masks[var]]
        
        output_ddiff.append(distribution_difference(df_s, df.sample(n = 10000)))


        #temp_info = compute_pairwise_mi_parallel(df_s, para_nm, n_jobs = 64).iloc[:2,1:-1]
        #temp_info["localvar"] = var
        #output_info.append(temp_info)

        
            
    output_ddiff = pd.concat(output_ddiff, axis = 1)
    output_ddiff.columns = vars
    output_ddiff = output_ddiff.transpose().apply(get_top_2_info, args = (n_large,), axis = 1)
    
    #output_info = pd.concat(output_info, axis = 0)

    if n_large == 3:
        output_ddiff["grouped_keys"] = output_ddiff[["topvpara1", "toppara2", "toppara3"]].apply(lambda row: tuple(sorted(row)), axis = 1)
    if n_large == 2:
        output_ddiff["grouped_keys"] = output_ddiff[["topvpara1", "toppara2"]].apply(lambda row: tuple(sorted(row)), axis = 1)
        
    localvars_grouped_by_3paras = output_ddiff.groupby("grouped_keys").localvar.apply(list).to_dict()
    


    strt_para3_localvar = {}
    surv_para3_localvar = {}
    for k, v in localvars_grouped_by_3paras.items():

        print(f"{k}: {tf_masks[v].all(axis = 1).sum()}")
        if tf_masks[v].all(axis = 1).sum() < group3_threshold :
            strt_para3_localvar[k] = v
        else:
            surv_para3_localvar[k] = v

    return localvars_grouped_by_3paras, strt_para3_localvar, surv_para3_localvar

        



def strterr_localvar_detection(strt_para3_localvars, tf_masks, n_comb3 = 2):
    summary_table_3to3 = {}


    for para3s, ks in tqdm(list(strt_para3_localvars.items())):
        
    
        localvars_temp = ks    
        pd_list = []
        
        for localvars_it in tqdm(combinations(localvars_temp, n_comb3)):
            pd_list.append(list(localvars_it) + [tf_masks[list(localvars_it)].all(axis = 1).sum()])
            
        pd_list = pd.DataFrame(pd_list, columns = [f"var{i+1}" for i in range(n_comb3)] + ["cunt"])
        pd_list = pd_list.sort_values(by = "cunt")
    
        summary_table_3to3[para3s] = pd_list

    return summary_table_3to3


def para1_para3_dict_generationg(surv_para3_localvars, tf_masks, emu_para):

    modified_paras = list(set([x for tup in list(surv_para3_localvars.keys()) for x in tup]))
    set_para3 = list(surv_para3_localvars.keys())
    
    single_dict = {}
    for p in modified_paras:
        p = set([p])
        temp_list = [t for t in set_para3 if p.issubset(t)]

        if len(temp_list) > 0:
            single_dict[tuple(p)] = temp_list


    single_min_max_range = {}
    para1_localvars_dict = {}
    
    for p1, p3s in list(single_dict.items()):
        pts_list = []
        localvar_list = []
        for p3 in p3s:

            
            temp_localvars = surv_para3_localvars[p3]
            
            localvar_list.extend(temp_localvars)
            
            temp_pts = emu_para[tf_masks[temp_localvars].all(axis = 1)]
            if temp_pts.shape[0] > 40000:
                temp_pts = temp_pts.sample(40000)
            pts_list.append(temp_pts[list(p1)].values)

        para1_localvars_dict[p1] = localvar_list
        temp_min  = np.array([lst.min() for lst in pts_list])
        temp_max  = np.array([lst.max() for lst in pts_list])
        single_min_max_range[p1] = [temp_min.max(), temp_max.min()]


    print([v for v, k in single_min_max_range.items() if k[1] - k[0] < 0])
    
    
    return para1_localvars_dict, single_min_max_range




def para1_error_localvar_detection(para1_localvars_dict, para_interest, emu_para,tf_masks, no_comb = 1):
    para_interest = tuple([para_interest])
    localvars_interest = para1_localvars_dict[para_interest]

    para1_range = []

    temp_min = []
    temp_max = []
    
    var_name1 = ["var_" + str(i) for i in np.arange(1)] + ["min"] + ["max"]
    var_name2 = ["exclude_var_" + str(i) for i in np.arange(no_comb)] + ["min"] + ["max"]
    
    for localvars_a in tqdm(list(combinations(localvars_interest, 1))):
        localvars_a = list(localvars_a)
        sub_samples = emu_para[tf_masks[localvars_a].all(axis = 1)][list(para_interest)]
        if sub_samples.shape[0] > 40000:
            sub_samples = sub_samples.sample(40000)
            
        localvars_a.extend([sub_samples.min().values[0], sub_samples.max().values[0]])
        para1_range.append(localvars_a)


    para1_range = pd.DataFrame(para1_range)
    para1_range.columns = var_name1


    paran_exclude_min_max = []


    for exclude_inds in combinations(np.arange(para1_range.shape[0]), no_comb):
        exclude_inds = list(exclude_inds)
        keep_indx = np.setdiff1d(np.arange(para1_range.shape[0]), exclude_inds)
        temp_minmax = para1_range.iloc[keep_indx, :]
        
        temp_min, temp_max = temp_minmax["min"].max(), temp_minmax["max"].min()
        per_row = [para1_range.iloc[i,0] for i in exclude_inds] + [temp_min, temp_max]
        paran_exclude_min_max.append(per_row)
    
    paran_exclude_min_max = pd.DataFrame(paran_exclude_min_max)
    paran_exclude_min_max.columns = var_name2

    paran_exclude_min_max_sel = paran_exclude_min_max[paran_exclude_min_max["max"] - paran_exclude_min_max["min"] > 0]
    
    return paran_exclude_min_max_sel



def para1_para2_dict_generation(surv_para3_localvar, tf_masks, emu_para, shape_alpha = 7):

    modified_paras = list(set([x for tup in list(surv_para3_localvar.keys()) for x in tup]))
    set_para3 = list(surv_para3_localvar.keys())

    double_dict = {}
    
    for pair in combinations(modified_paras, 2):
        pair = set(pair)
        temp_list = [t for t in set_para3 if pair.issubset(t)]


        if len(temp_list) > 0:
            double_dict[tuple(pair)] = temp_list

    pairs_hulls = {}
    para2_localvars_dict = {}

    for para2, para3 in tqdm(double_dict.items()):
        pts_list = []
        localvars_list = []
        for p3 in para3:
            temp_localvars = surv_para3_localvar[p3]
            localvars_list.extend(temp_localvars)
            temp_pts = emu_para[tf_masks[temp_localvars].all(axis = 1)]
            print(temp_pts.shape)
            if temp_pts.shape[0] > 20000:
                temp_pts = temp_pts.sample(20000)
            pts_list.append(temp_pts[list(para2)].values)

        hulls = [alphashape.alphashape(points, shape_alpha) for points in pts_list]
        
        intersection_hulls = reduce(lambda a, b: a.intersection(b), hulls)
        pairs_hulls[para2] = intersection_hulls
        
        para2_localvars_dict[para2] = localvars_list

    print([v for v, k in list(pairs_hulls.items()) if k.is_empty])

    
    return para2_localvars_dict, pairs_hulls






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



# def distribution_difference(df1, df2, bins=100):
#     abs_den_diff = {}
#     for column in df1.columns:
#         hist1, bin_edges1 = np.histogram(df1[column], bins=bins, density=True, range = (0,1))
#         hist2, bin_edges2 = np.histogram(df2[column], bins=bins, density=True, range = (0,1))
#     # Compute bin widths and probabilities
#         abs_den_diff[column] = np.sum(abs(hist1 - hist2))

#     output = pd.Series(abs_den_diff)

#     return pd.Series(output)
