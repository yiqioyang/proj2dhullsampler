import numpy as np
import pandas as pd
import os
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





def test_ind_vars(X_prev, X, para_nm, tf_masks, grouped_hulls, para, paras_vars, threshold_ratio_between_para_pairs, shape_alpha = 5):
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
            #xxx
            if (attempt is not None) and (not attempt.empty) and (attempt.shape[0]/ X_prev.shape[0] > threshold_ratio_between_para_pairs):
                print(f'\t \t \t \t Found the good variable combo')
                drop_vars = [x for x in vars if x not in list(var_comb)]
                return list(var_comb), drop_vars, hull_sub, attempt
            else:
                pass 

        if i -1 == 0:
            print(f'\t \t \t No variable combinations work')
            return None





def orchestrate_test(para_seq, X, tf_masks, para_nm, grouped_hulls, paras_vars, n_pts=10000, n_threshold=10000, sample_threshold = 10**7, max_workers=None, threshold_ratio_between_para_pairs = 0.02):
    para_l = []
    
    var_drop = {}
    error_sample_size_scaling = 1.0
    out_prev = None
    prev_sample_size = sample_threshold

    for p_count, p in enumerate(para_seq):
        sample_threshold = int(sample_threshold)
        prev_sample_size = int(prev_sample_size)
        print(f'Running {p}, the {p_count}th simulation')
        print(f'Current sampling size is {sample_threshold}')
        print(f'Previous sample size is {prev_sample_size}')
        
        para_l.append(p)
        out = sample_from_hulls_n(para_l, para_nm, grouped_hulls, n_pts=  n_pts, n_threshold = n_threshold, sample_threshold = sample_threshold, max_workers = max_workers)
        
        if (out is None) or (out.shape[0]/prev_sample_size < threshold_ratio_between_para_pairs):
            print("\t Find nothing or the shrink is too rapid, try to resolve it by breaking the variables into groups")
            out_prev = sample_from_hulls_n(para_l[:-1], para_nm, grouped_hulls, n_pts=  n_pts, n_threshold = n_threshold, max_workers = max_workers, sample_threshold=sample_threshold)
            if (out_prev is None):
                raise ValueError("out_prev is None")

            print(f'The size of out_prev is {out_prev.shape[0]}')
            if (out_prev.shape[0] < 100) & (out_prev.shape[0] > 0) & (out_prev is not None):
                error_sample_size_scaling = min(100.0/out_prev.shape[0], 100)
                sample_threshold = sample_threshold * error_sample_size_scaling
                out_prev = sample_from_hulls_n(para_l[:-1], para_nm, grouped_hulls, n_pts=  n_pts, n_threshold = n_threshold, max_workers = max_workers, sample_threshold=sample_threshold)
                print(f'\t \t Sample size too small, have to increase the sample size to {sample_threshold}')
                print(f'\t \t based on the scheudler ratio of {error_sample_size_scaling}')
                print(f'\t \t The size of out_prev is updated to {out_prev.shape[0]}')
                
                

            check_pt = test_ind_vars(out_prev, X, para_nm, tf_masks, grouped_hulls, p, paras_vars, threshold_ratio_between_para_pairs, shape_alpha = 5)
            if check_pt is None:
                para_l.remove(p)
                
                var_drop[p] = paras_vars[p]
                del grouped_hulls[p]
                del paras_vars[p]
                print(f'\t \t \t {p} is causing trouble and is skipped')
                prev_sample_size = out_prev.shape[0]
                
                
            else:
                print(f'\t \t \t Modify the original hull and variable')
                print(f'\t \t \t Drop {check_pt[1]}')
                grouped_hulls[p] = check_pt[2]
                paras_vars[p] = check_pt[0]
                var_drop[p] = check_pt[1]
                prev_sample_size = check_pt[3].shape[0]
                sample_threshold = sample_threshold * min(100.0/prev_sample_size, 100)
                print(f'\t \t \t After exluding some variables, the number of surviving samples is {prev_sample_size}')
                print(f' \t \t \t We also increase the sample threshold by {min(100.0/prev_sample_size, 100)}')

                
        else:
            print(f'\t There is overlap for {p}. There are {out.shape[0]} surviving ensemble members. ')
            if out.shape[0] < 100:
                error_sample_size_scaling = min(100.0/out.shape[0], 100)
                sample_threshold = sample_threshold * error_sample_size_scaling

                prev_sample_size = out.shape[0] * error_sample_size_scaling
                print(f'\t \t Since the surviving sample size is too small, we have to increase the sample threshold by {error_sample_size_scaling}')
                print(f'\t \t The ***estimated*** # of surviving pts is {prev_sample_size}')
            else:
                prev_sample_size = out.shape[0]
                print(f'\t \t The ***actual*** # of surviving pts is {prev_sample_size}')
                print(f'\t \t Proceed to the next iteration')
        print("===============================End of one iteration=======================================")


        
    return grouped_hulls, para_l, paras_vars, var_drop, out_prev