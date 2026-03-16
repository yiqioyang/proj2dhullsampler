import xarray as xr
import pandas as pd
from pathlib import Path

import alphashape
from itertools import combinations

from .sampling_functions import orchestrate_test, sample_from_hulls_n
from .aux import para_csv2nc


def meta_one_hot_shot(meta, para_nm):
    meta = meta.transpose()
    meta_one_hot = pd.DataFrame(False, index=meta.index, columns=para_nm)
    for index, row in meta.iterrows():
        for r in row.values:
            meta_one_hot.loc[index, para_nm[r]] = True

    return meta_one_hot

class EmulatedDataStorage:
    """
    Lightweight container for emulation outputs.
    No computation logic.
    """
    def __init__(self):
        pass

class HistoryMatching:
    def __init__(self, working_dir, case_name, ppe_para, threshold_level = 2.0):
        self.root = Path(working_dir) / case_name
        self.tf_masks = pd.read_csv(self.root / f'tf_masks_level_{threshold_level}.csv', index_col=0)
        
        self.meta = pd.read_csv(self.root / 'meta.csv', index_col = 0)
        self.p_emu = xr.open_dataset(self.root / 'sampled_parameters.nc').to_dataframe()
        self.var_nm = list(self.tf_masks.columns)
        self.para_nm = list(self.p_emu.columns)
        self.ppe_para = ppe_para
        self.tf_masks_raw = self.tf_masks

        self.dropped_vars = EmulatedDataStorage()
        self.n_sample = self.tf_masks.shape[0]

        self.results = EmulatedDataStorage()
        
        self.dropped_vars.nooverlap2d = []
        
    def drop_by_name(self, var_to_exclude):
        var_to_drop = []
        for v in var_to_exclude:
            var_to_drop.extend([s for s in self.var_nm if s.startswith(v)])
        
        self.tf_masks = self.tf_masks.drop(columns = var_to_drop)

        self.var_nm = list(self.tf_masks.columns)
        self.dropped_vars.by_name = var_to_drop

    def drop_by_n_survive(self, n_survive):
        survive_summary = self.tf_masks.sum(axis = 0)
        self.dropped_vars.useless = list(survive_summary[survive_summary == self.n_sample].index)
        self.dropped_vars.tight   = list(survive_summary[survive_summary < n_survive].index)

        self.tf_masks = self.tf_masks.drop(columns = self.dropped_vars.useless + self.dropped_vars.tight)
        
        self.var_nm = list(self.tf_masks.columns)

    
    def drop_by_nvar_per_pair(self, n_var_thre = 1):
        self.dropped_vars.local = []
        for k, v in list(self.paras_vars.items()):
            if len(v) <= n_var_thre:
                self.dropped_vars.local.append(v)
                del self.paras_vars[k]
    
    def update_meta(self, occurence_threshold = 2):
        self.meta = self.meta[self.var_nm]
        self.meta_onehot = meta_one_hot_shot(self.meta, self.para_nm)
        
        

    def hull_for_each(self, shape_alpha = 5):
        hull_per_var = {}
        for v in self.var_nm:
            p_ind = self.meta[v].sort_values().values
            pts = self.p_emu[self.tf_masks[v]].iloc[:,p_ind]
            if pts.shape[0] > 5000:
               pts = pts.sample(5000)

            pts = pts.values
            
            hull_per_var[v] = alphashape.alphashape(pts, shape_alpha)

        self.hull_per_var = hull_per_var

    def group_para_climatology(self):

        vars = self.var_nm
        paras_vars = {}
        paras_vars_0 = {}    
        for c in vars: 
            para_inds = self.meta[c].sort_values().values
            para_nms = [self.para_nm[ind] for ind in para_inds]
            paras_vars.setdefault(tuple(para_nms), []).append(c)



        paras_vars = dict(
                sorted(paras_vars.items(), key=lambda item: len(item[1]), reverse=True)
        )
        
        self.paras_vars = paras_vars
        
        
        for k, v in paras_vars.items():
            temp_count = self.tf_masks[v].all(axis = 1).sum()
            print(f'{k[0]:<40} and {k[1]:<40}: {temp_count:>8}')
            if temp_count < 10000:
                paras_vars_0[k] = v

        self.paras_vars_0 = paras_vars_0

        
        
    def shuffle_vars(self, n_comb = 2):
        summary_table = {}
    
        for paras, vars in self.paras_vars_0.items():
            
            vars_temp = vars 
            pd_list = []
            
            for vars_comb in combinations(vars_temp, n_comb):
                pd_list.append(list(vars_comb) + [self.tf_masks[list(vars_comb)].all(axis = 1).sum()])
                
            pd_list = pd.DataFrame(pd_list, columns = [f"var{i+1}" for i in range(n_comb)] + ["count"])
            pd_list = pd_list.sort_values(by = "count")
            
            
            summary_table[paras] = pd_list

        print(f'There are {len(self.paras_vars_0)} groups that have no overlapping within own groups')
        return summary_table

    def drop_no_overlap2d_vars(self, vars_to_drop):
        self.tf_masks = self.tf_masks.drop(columns = vars_to_drop)
        self.var_nm = list(self.tf_masks.columns)
        self.meta = self.meta[self.var_nm]
        self.meta_onehot = meta_one_hot_shot(self.meta, self.para_nm)
        self.dropped_vars.nooverlap2d.append(vars_to_drop)

        
    
    def visualize(self, para_pair):
        para_pair = tuple(para_pair)
        survive_pts = self.p_emu[self.tf_masks[self.paras_vars[para_pair]].all(axis = 1)]

        if survive_pts.shape[0] > 5000:
            survive_pts = survive_pts.sample(5000)

        return survive_pts

    def build_hulls(self, shape_alpha = 5):
        grouped_hulls = {}

        for para2, vars in self.paras_vars.items():

            tf_mask = self.tf_masks[vars].all(axis = 1)
            pts = self.p_emu[tf_mask]
            if pts.shape[0] > 5000:
                pts = pts.sample(5000)
            
            pts_np = pts[list(para2)].values
            
            grouped_hulls[para2] = alphashape.alphashape(pts_np, shape_alpha)
            
            
        self.grouped_hulls = grouped_hulls
    
    def rescale_para(self, sampled_para):
        ppe_para = self.ppe_para
        return (ppe_para.max() - ppe_para.min()) * sampled_para + ppe_para.min()
    

    def orchestrate(self, n_pts = 10000, n_threshold = 100, sample_threshold = 10**5, max_workers = 31):


        para_seq = list(self.grouped_hulls.keys())
        
        check = orchestrate_test(para_seq, self.p_emu, self.tf_masks,  
                         self.para_nm, self.grouped_hulls, self.paras_vars, n_pts, n_threshold, sample_threshold, max_workers)

        self.results.valid_hulls = check[0]
        self.results.para_l = check[1]
        self.dropped_vars.during_iteration = check[3]

    def draw(self, n_pts=50000, n_threshold=5000, sample_threshold=10**8, max_workers=32, n_max = 1000):
        valid_hulls = self.results.valid_hulls
        samples = sample_from_hulls_n(list(valid_hulls.keys()), self.para_nm, valid_hulls, n_pts, n_threshold, max_workers, sample_threshold)
        if samples.shape[0]>n_max:
            samples = samples.iloc[:n_max]
            
        self.results.unscaled_samples = samples
        self.results.realscale_samples = self.rescale_para(samples)



    def save_samples(self, n = 100):
        csv_path1 = self.root / 'full_sel_para_realscale.csv'
        nc_path1 = self.root / 'full_sel_para_realscale.nc'

        csv_path2 = self.root / 'sel_para_realscale.csv'
        nc_path2 = self.root / 'sel_para_realscale.nc'


        
        self.results.realscale_samples.to_csv(csv_path1)
        para_csv2nc(csv_path1, nc_path1, self.results.realscale_samples.shape[0])

        self.results.realscale_samples.iloc[:n,:].to_csv(csv_path2)
        para_csv2nc(csv_path2, nc_path2, n)
        
        