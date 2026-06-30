

import numpy as np
import pandas as pd
import itertools


from tqdm import tqdm
#from functools import reduce
#from functools import partial
import trimesh
#from shapely import points, contains
#from shapely.geometry import Point
#from collections import Counter

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, WhiteKernel, ConstantKernel as C
#from concurrent.futures import ProcessPoolExecutor

from scipy.spatial import Delaunay
#import alphashape

import joblib

#import plotly.graph_objects as go


def gp_training_application(X, Y, y_name, X_emu, path, n_sens_p = 2, no_restart = 10):
    y = Y[y_name]

    if y.isna().any():
        raise ValueError(
            f"Diagnostic {y_name!r} contains missing values"
        )

    y_std = y.std()

    if not np.isfinite(y_std) or y_std == 0:
        raise ValueError(
            f"Diagnostic {y_name!r} has zero or invalid variance"
        )

    y_norm = (y - y.mean())/y.std()


    sage1d_temp = fit_all_gp_models_1d(y_norm.values, X.values, len1d = 0.4)
    
    
    sel_para_ind = list(sage1d_temp.index[:n_sens_p])
    

    kernel = C(1.0, (1e-3, 0.8)) * Matern(length_scale=0.001, nu=2.5, length_scale_bounds=(0.1, 3)) + WhiteKernel(noise_level=0.5, noise_level_bounds=(0.01, 0.9))
    gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer= no_restart, normalize_y=True)
    
    gp.fit(X.values[:,sel_para_ind], y_norm.values)
    
    joblib.dump(gp, path + "/python_obj/" + y_name + "_gpmodel.pkl")

    y_mean_emu, y_std_emu = gp.predict(X_emu.values[:,sel_para_ind], return_std=True)
    y_mean_emu = y_mean_emu.reshape(-1, 1)
    y_std_emu = y_std_emu.reshape(-1, 1)

    y_mean_std_emu = np.hstack([y_mean_emu, y_std_emu])
    y_mean_std_emu = pd.DataFrame(y_mean_std_emu, columns =[y_name + "_mean", y_name + "_std"])
    y_mean_std_emu.to_csv(path + "/y_emu" + "/gp_mean_std_" + y_name + ".csv")
    
    
    return ([y_name, sel_para_ind])


        

def fit_gp_for_single_1d(y, X, i, length_scale=0.4):
    """Fits a GP model for a given pair of feature indices."""
    
    X_single = X[:, [i]].reshape([-1,1])  
    
    # Define GP kernel
    kernel =  Matern(length_scale=length_scale, 
                     length_scale_bounds="fixed", nu = 2.5)  + WhiteKernel(noise_level = 0.4)
    
    # Fit GP model
    gp = GaussianProcessRegressor(kernel=kernel, optimizer=None)
    gp.fit(X_single, y.ravel())  # Flatten y

    y_pred = gp.predict(X_single, return_std = False)
    residual_sd = (y - y_pred).std()
    
    return (i, residual_sd)



def fit_all_gp_models_1d(y, X, len1d = 0.4):
    """Fits GP models for all unique feature pairs and stores residual std."""
    residual_single_stds = {}
    
    
    
    feature_single = itertools.combinations(range(X.shape[1]), 1)
    

    for i in feature_single:
        residual_single_stds[i] = fit_gp_for_single_1d(y, X, i, len1d)[1]  # Store only the residual std

    
    sorted_delta = sorted(residual_single_stds.items(), key=lambda x: x[1])
    
    top_10 = dict(sorted_delta[:20])
    top_10 = pd.Series(top_10)
    top_10.index = top_10.index.map(lambda top_10: top_10[0])
    
    return top_10




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


















##xxx###
############################################################
