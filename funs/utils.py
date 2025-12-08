
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from itertools import permutations
from itertools import combinations

from tqdm import tqdm
from tqdm.notebook import tqdm
from functools import reduce
from functools import partial
import trimesh
from shapely import points, contains
from collections import Counter


from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, WhiteKernel, ConstantKernel as C
from concurrent.futures import ProcessPoolExecutor, as_completed



from scipy.spatial import Delaunay

import numpy as np
import numpy as np
import matplotlib.pyplot as plt
import alphashape
from shapely.geometry import Point
import joblib



def gp_training_application(X, Y, y_name, X_emu, path = "/glade/work/qingyuany/repo_data/spatialtuning/", n_sens_p = 2, no_restart = 10):
    y = Y[y_name]
    y_norm = (y - y.mean())/y.std()

    if ~(np.isnan(y_norm).all()):
        sage1d_temp = fit_all_gp_models_1d(y_norm.values, X.values, len1d = 0.4)
        
        sel_para_ind = sage1d_temp[:(n_sens_p * 2):2].astype(int)

        #sel_para_ind = sel_para_ind.reshape(1,-1)

        kernel = C(1.0, (1e-3, 0.8)) * Matern(length_scale=0.001, nu=2.5, length_scale_bounds=(0.1, 3)) + WhiteKernel(noise_level=0.5, noise_level_bounds=(0.01, 0.9))
        gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer= no_restart, normalize_y=True)
        
        gp.fit(X.values[:,sel_para_ind], y_norm.values)
        
        joblib.dump(gp, path + "python_obj/" + y_name + "_gpmodel.pkl")

        y_mean_emu, y_std_emu = gp.predict(X_emu.values[:,sel_para_ind], return_std=True)
        y_mean_emu = y_mean_emu.reshape(-1, 1)
        y_std_emu = y_std_emu.reshape(-1, 1)

        y_mean_std_emu = np.hstack([y_mean_emu, y_std_emu])
        y_mean_std_emu = pd.DataFrame(y_mean_std_emu, columns =[y_name + "_mean", y_name + "_std"])
        y_mean_std_emu.to_csv(path + "emulated_ensemble" + "/gp_mean_std_" + y_name + ".csv")
        #pd.DataFrame(y_mean_emu, columns=[y_name]).to_csv(path + "emulated_ensemble" + "/gp_mean_" + y_name + ".csv")
        #pd.DataFrame(y_std_emu, columns=[y_name]).to_csv(path + "emulated_ensemble" + "/gp_std_" + y_name + ".csv")

        
        return ([y_name, sel_para_ind])

        

def fit_gp_for_single_1d(y, X, i, length_scale=0.5):
    """Fits a GP model for a given pair of feature indices."""
    
    X_single = X[:, [i]].reshape([-1,1])  
    
    # Define GP kernel
    kernel =  Matern(length_scale=length_scale, 
                     length_scale_bounds="fixed", nu = 2.5)  + WhiteKernel(noise_level = 0.8)
    
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





def range_err_detection2(paras_vars, tf_masks, emu_para, meta):

    para_inds = np.sort(pd.unique(meta.values.ravel()))
    para_nm = list(emu_para.columns)
    
    para1_vars = {}
    para1_para2 = {}
    para1_min_max = {}
    
    for p in para_inds:
        for c in tf_masks.columns:
            if p in meta[c].values:
                para1_vars.setdefault(para_nm[p], []).append(c)
    
    for p in para_inds:
        for p2 in list(paras_vars.keys()):
            if para_nm[p] in p2:
                para1_para2.setdefault(para_nm[p], []).append(p2)
    
    
    
    for p1, p2s in list(para1_para2.items()):
        pts_list = []
        for p2 in p2s:    
            temp_pts = emu_para[tf_masks[paras_vars[p2]].all(axis = 1)]
            if temp_pts.shape[0] > 40000:
                temp_pts = temp_pts.sample(40000)
            pts_list.append(temp_pts[p1].values)
            
        temp_min  = np.array([lst.min() for lst in pts_list])
        temp_max  = np.array([lst.max() for lst in pts_list])
        para1_min_max[p1] = [temp_min.max(), temp_max.min()]

    
    para1_min_max = pd.DataFrame(para1_min_max).T
    para1_min_max.columns = ["min", "max"]

    print(para1_min_max[para1_min_max["min"] >=  para1_min_max["max"]])
    
    return para1_vars, para1_para2, para1_min_max






def para1_error_localvar_detection(para1_vars_dict, p, emu_para, tf_masks, n_comb = 1):
    ### take parameter p, 
    ### take combinations n_comb of the variables that are connected to p,
    ### calculate the min and max for each of the point clouds constrained by var_comb ***
    ### 
    p = [p]
    vars = para1_vars_dict[tuple(p)[0]]
    

    vars_min_max = []
    
    for var_comb in list(combinations(vars, n_comb)):
        used_vars = [v for v in vars if v not in var_comb]
        
        temp_pts = emu_para[tf_masks[used_vars].all(axis = 1)][p]
        if tf_masks[used_vars].all(axis = 1).sum() > 0:

            vars_min_max.append(pd.Series([var_comb, temp_pts.min().values, temp_pts.max().values]))

    output = pd.DataFrame(vars_min_max)
    output.columns = ["included_vars", "min", "max"]
    print(output)
    return output
    



def para2_error_detection(paras_vars, tf_masks, emu_para, meta, shape_alpha = 7):
    ## Need to check
    sel_col = []
    for k, v in paras_vars.items():
        sel_col.extend(v)

    meta = meta[sel_col]
    tf_masks = tf_masks[sel_col]
    print(meta.shape)
    
    para_inds = np.sort(pd.unique(meta.values.ravel()))
    para_nm = list(emu_para.columns)

    para2_vars = {}
    para2_vars_s = {}
    para2_area = {}
    for p2 in combinations(para_inds, 2):
        for c in tf_masks.columns:
            if set(p2).issubset(meta[c].values):
                para_nm2 = tuple(sorted([para_nm[i] for i in list(p2)]))
                para2_vars.setdefault(para_nm2, []).append(c)
    
    
    ###
    pairs_hulls = {}
    ##
    for para2, vars in tqdm(para2_vars.items()):
        
        pts_list = [emu_para[tf_masks[var]][list(para2)] for var in vars]
        pts_list =  [pts.sample(n=5000).values if pts.shape[0] > 5000 else pts.values for pts in pts_list]   
        
        hulls = [alphashape.alphashape(points, shape_alpha) for points in pts_list]
        
        intersection_hulls = reduce(lambda a, b: a.intersection(b), hulls)
        pairs_hulls[para2] = intersection_hulls
        para2_area[para2] = intersection_hulls.area
        
        
        if (intersection_hulls.area < 0.01) | (np.isnan(intersection_hulls.area)):
            para2_vars_s[para2] = vars
        
    

    return para2_vars, para2_vars_s, para2_area, pairs_hulls





def para2_error_localvar_detection(para2_vars, para2, emu_para,tf_masks, n_comb = 2,  shape_alpha = 7):
    ## xx ??
    para2 = tuple(sorted(para2))
    vars = para2_vars[para2]
    print(len(vars))    

    vars_area = []
    
    
    for include_vars in tqdm(combinations(vars, n_comb)):
        include_vars = list(include_vars)

        
        pts_list = [emu_para[tf_masks[var]][list(para2)] for var in include_vars]
        pts_list =  [pts.sample(n=5000).values if pts.shape[0] > 5000 else pts.values for pts in pts_list]   

        hulls = [alphashape.alphashape(points, shape_alpha) for points in pts_list]
        intersection_hulls = reduce(lambda a, b: a.intersection(b), hulls)

        vars_area.append(include_vars + [intersection_hulls.area])
        
    
    vars_area = pd.DataFrame(vars_area, columns = [f"var{i+1}" for i in range(n_comb)] + ["area"])
    vars_area = vars_area.sort_values(by = "area")
    
    
    return vars_area



def para2_group_summary(para2_hulls):

    para2s = list(para2_hulls.keys())
    
    ps = set([p for v, k in list(para2_hulls.items()) for p in v])
    para2s_flat = [p for v, k in list(para2_hulls.items()) for p in v]
    count = list(Counter(sorted(para2s_flat)).most_common())
    return count


def group_by_para(para2_hulls, p):
    grouped_hulls  = {}
    
    for k, v in para2_hulls.items():
        if p in list(k):
            grouped_hulls[k] = v

    return grouped_hulls



def group_by_rest(para2_hulls, list_hulls):
    hull_names = []
    rest_hulls = {}
    for l in list_hulls:
        hull_names.extend(l.keys())

    for k, v in para2_hulls.items():
        if k not in hull_names:
            rest_hulls[k] = v

    return rest_hulls


##xxx###
############################################################

def para3_mesh_generator(paras_vars, emu_para, tf_masks):
    para3_meshes = {}
    for k, v in tqdm(paras_vars.items()):
        temp_pts = emu_para[tf_masks[v].all(axis = 1)]
        if temp_pts.shape[0] > 5000:
            temp_pts = temp_pts.sample(5000)
            
        temp_pts_np = temp_pts[list(k)].values
        temp_pts_index = alpha_shape_3D(temp_pts_np, 7)
    
        para3_meshes[k] = trimesh.Trimesh(vertices=temp_pts_np, faces=list(temp_pts_index))

    
    return para3_meshes
#############################################################

def ratio_cal2d(poly_dict):
    ratio = {}
    hulls2d = {}
    for k, polyg in poly_dict.items():
        sq_area = (polyg.bounds[3] - polyg.bounds[1]) *  (polyg.bounds[2] - polyg.bounds[0])
        ratio[k] = polyg.area/sq_area

    ratio = pd.Series(ratio).sort_values()
    
    for k in list(ratio.index):
        hulls2d[k] = poly_dict[k]
    return ratio, hulls2d
    


def ratio_cal3d(hull_dict):
    ratio = {}
    hulls3d = {}
    for k, hull in hull_dict.items():
        
        sq_vol = np.prod(hull.extents)
        ratio[k] = hull.volume/sq_vol

    ratio = pd.Series(ratio).sort_values()
    
    for k in list(ratio.index):
        hulls3d[k] = hull_dict[k]
    return ratio, hulls3d
    



############################################################

def sample2d(no_pts, para_nm, single_min_max_range, poly_dict, existing_pts = np.NAN, monitor = True):
    
    if not isinstance(existing_pts, pd.DataFrame):
        if np.isnan(existing_pts):
            random_pts = pd.DataFrame(np.random.uniform(0, 1, (no_pts, len(para_nm))), columns = para_nm)
        #for k, v in single_min_max_range.items():
        for p, row in single_min_max_range.iterrows():
            random_pts[list([p])] = random_pts[list([p])] * (row['max'] - row['min']) + row['min']
        
    else:
        random_pts = existing_pts
    
    for k, polyg in poly_dict.items():
        temp_pts = random_pts[list(k)].values
        temp_pts = points(temp_pts)
        
        temp_masks = contains(polyg,temp_pts)
        if monitor:
            print(f'{k}: {temp_masks.sum()}')
        if temp_masks.sum() == 0:
            print("Didn't find anything")
            return []

        random_pts = random_pts[temp_masks].reset_index(drop=True)

    return random_pts




def sample2d_wrapper(hulls_list, para_nm, single_min_max_range, exist_pts = np.nan, n_pts = 1000000, rand_int = 1):
    
    np.random.seed()
    if isinstance(exist_pts, pd.DataFrame):
        random_pts = exist_pts
    else:
        random_pts = pd.DataFrame(np.random.uniform(0, 1, (n_pts, len(para_nm))), columns = para_nm)
        for p, row in single_min_max_range.iterrows():
            random_pts[list([p])] = random_pts[list([p])] * (row['max'] - row['min']) + row['min']
    
    for k, hull in hulls_list.items():

        temp_pts = random_pts[(sorted(list(k)))].values
        temp_pts = points(temp_pts)
        
        temp_mask = contains(hull, temp_pts)
        if temp_mask.sum() == 0:
            return []
            
        #random_pts = random_pts[temp_mask]
        random_pts = random_pts[temp_mask].reset_index(drop=True)

    return random_pts



def sample2d_perlist(hulls_list, para_nm, single_min_max_range, n_iter, exist_pts = np.nan, n_pts= 1000000, n_cpu = 30, n_need = 1000000):
    
    ################################
    pts = []
    count = 0
    for i in np.arange(0, n_iter):
        # Use all CPUs (or set max_workers)
        with ProcessPoolExecutor(max_workers= n_cpu) as ex:
            worker = partial(sample2d_wrapper, hulls_list, para_nm, single_min_max_range, exist_pts, n_pts)
            pts_list = list(ex.map(worker, np.arange(n_cpu)))
            pts_list = [x for x in pts_list if isinstance(x, pd.DataFrame)]
            
            pts_list = pd.concat(pts_list, axis = 0, ignore_index=True)
            count = count +   pts_list.shape[0]
        pts.append(pts_list)
        if count > n_need:
            break
    
    pts = pd.concat(pts, axis = 0, ignore_index= True)
    return pts
########################################








def sample3d(no_pts, para_nm, single_min_max_range, mesh_dict, existing_pts = np.NAN, monitor = True):
    if not isinstance(existing_pts, pd.DataFrame):
        if np.isnan(existing_pts):
            random_pts = pd.DataFrame(np.random.uniform(0, 1, (no_pts, len(para_nm))), columns = para_nm)
        for k, v in single_min_max_range.items():
            random_pts[list([k])] = random_pts[list([k])] * (v[1] - v[0]) + v[0]
        
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

