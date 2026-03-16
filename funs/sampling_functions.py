import numpy as np
import pandas as pd
import os
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import ProcessPoolExecutor, as_completed
import alphashape
from shapely import points, contains
from itertools import combinations


def sample_from_hull(X, para, h):

    minx, miny, maxx, maxy = h.bounds
    p1, p2 = para

    X = X[
        (X[p1] >= minx) & (X[p1] <= maxx) &
        (X[p2] >= miny) & (X[p2] <= maxy)
    ]

    if X.empty:
        return X
    
    pts = points(X[list(para)].values)
    tf = contains(h,pts)
    
    return X.loc[tf].reset_index(drop=True)



def _one_batch(args):
    para_l, para_nm, grouped_hulls, n_pts, seed = args

    rng = np.random.default_rng(seed)

    
    X = pd.DataFrame(
        rng.uniform(0, 1, size=(n_pts, len(para_nm))),
        columns=para_nm
    )

    for para in para_l:
        h = grouped_hulls[para]
        X = sample_from_hull(X, para, h)

        if X.shape[0] == 0:
            return None

    return X


def sample_from_hulls_n(
    para_l,
    para_nm,
    grouped_hulls,
    n_pts=100_000,
    n_threshold=100_000,
    max_workers=None,
    sample_threshold = 10**7
):
    if max_workers is None:
        max_workers = os.cpu_count() -1     
    
    
    out = []
    count = 0
    n_sampled = 0
    
    with ProcessPoolExecutor(max_workers=max_workers) as ex:
        futures = []
        MAX_IN_FLIGHT = 2 * max_workers

        while (count < n_threshold) and (n_sampled < sample_threshold):
            while (count < n_threshold) and (len(futures) < MAX_IN_FLIGHT) and (n_sampled < sample_threshold):
                seed = np.random.randint(0, 2**32 - 1)# changed
                
                futures.append(ex.submit(_one_batch,
                    (para_l, para_nm, grouped_hulls, n_pts,seed),))

                n_sampled = n_sampled + n_pts

            
            # harvest finished jobs
            for fut in as_completed(list(futures)):
                futures.remove(fut)
                X = fut.result()
                out.append(X)

                
                if X is not None:
                    count += X.shape[0]
        
                
                
        if all([x is None for x in out]):
            print(f'Find nothing from {sample_threshold} pts')
            return None
    
    out = [x for x in out if x is not None]
    out = pd.concat(out, axis=0).reset_index(drop = True)
    return out





def test_ind_vars(X_prev, X, para_nm, tf_masks, grouped_hulls, para, paras_vars, shape_alpha = 5):
    print(f'\t \tRunning test to see if {para} could be break down and lead to non-overlapping')

    vars = paras_vars[para]
    if len(vars) == 1:
        return None
    for i in reversed(range(len(vars))):
        print(f'Try {i} combinations')
        var_combs = combinations(vars, i)
        for var_comb in var_combs:
            print('Try one of the combinations')
            X_sub = X[tf_masks[list(var_comb)].all(axis = 1)]
            if X_sub.shape[0] > 5000:
                X_sub = X_sub.sample(5000)
            
            X_sub = X_sub[list(para)].values
            hull_sub = alphashape.alphashape(X_sub, shape_alpha)

            attempt = sample_from_hull(X_prev, para, hull_sub)
            
            if not attempt.empty:
                print(f'\t \t \t Found the good variable combo')
                drop_vars = [x for x in vars if x not in list(var_comb)]
                return list(var_comb), drop_vars, hull_sub
            else:
                pass 

        if i -1 == 0:
            print(f'\t \t \t No variable combinations work')
            return None






def orchestrate_test(para_seq, X, tf_masks, para_nm, grouped_hulls, paras_vars, n_pts=10000, n_threshold=10000, sample_threshold = 10**7, max_workers=None):
    para_l = []
    
    var_drop = {}
    error_sample_size_scaling = 1
    out_prev = None
    non_over_count = 0
    for p_count, p in enumerate(para_seq):
        print(f'Running {p}, the {p_count}th simulation')
        para_l.append(p)
        out = sample_from_hulls_n(para_l, para_nm, grouped_hulls, n_pts=  n_pts, n_threshold = n_threshold, sample_threshold = sample_threshold, max_workers = max_workers)
        
        if out is None:
            print("Find nothing, try to resolve it by breaking the variables into groups")
            print("First sample out_prev that needs greater sample size, which will take long")
            out_prev = sample_from_hulls_n(para_l[:-1], para_nm, grouped_hulls, n_pts=  n_pts, n_threshold = n_threshold, max_workers = max_workers, sample_threshold=sample_threshold * 120 * round(error_sample_size_scaling))
            
            print(f'The size of out_prev is {out_prev.shape[0]}')
            if (out_prev.shape[0] < 100) & (out_prev.shape[0] > 0) & (out_prev is not None):
                error_sample_size_scaling = min(100.0/out_prev.shape[0], 100)
                print(f'Increase the scheudler ratio to {error_sample_size_scaling}')
                
            if (out_prev is None):
                print("out_prev issue")
    
            check_pt = test_ind_vars(out_prev, X, para_nm, tf_masks, grouped_hulls, p, paras_vars, shape_alpha = 5)
            if check_pt is None:
                para_l.remove(p)
                
                var_drop[p] = paras_vars[p]
                del grouped_hulls[p]
                del paras_vars[p]
                print(f'\t \t \t \t {p} is causing trouble and is skipped')
                non_over_count = non_over_count + 1
                
                
                
            else:
                print(f'\t \t \t \t Modify the original hull and variable')
                print(f'\t \t \t \t Drop {check_pt[1]}')
                grouped_hulls[p] = check_pt[2]
                paras_vars[p] = check_pt[0]
                var_drop[p] = check_pt[1]
                
        else:
            print(f'There is overlap for {p}. Proceed to the next parameter pair')
            

        print("======================================================================")


        
        
    return grouped_hulls, para_l, paras_vars, var_drop, out_prev