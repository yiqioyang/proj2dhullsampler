import xarray as xr
import pandas as pd
from pathlib import Path
import json
import alphashape
from itertools import combinations
#from joblib import Parallel, delayed

from .sampling_functions import orchestrate_test, sample_from_hulls_n
from .plotting import visualize_emulation
from .aux import para_csv2nc
from .prep_class import Prepare_Case
import glob
import matplotlib.pyplot as plt


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
    def __init__(self, working_dir, case_name):
        self.working_dir = working_dir
        self.case_name = case_name
        self.root = Path(working_dir) / case_name

        
    def create_case(self, para, tabs, ppe, obs, obs_dict, lat_bins, manul_ppe_info, n_sample):
        if self.root.exists():
            print("Directory already exists")
            return 
        else:
            print("Start creating new case")
            prep_case = Prepare_Case(self.working_dir, self.case_name, para, tabs, ppe, obs, obs_dict, lat_bins, manul_ppe_info, n_sample)
            self.prep_case = prep_case
            

    def load_case(self):
        if self.root.exists():    
            self.meta = pd.read_csv(self.root / 'meta.csv', index_col = 0)
            self.p_emu = xr.open_dataset(self.root / 'sampled_parameters.nc').to_dataframe()
            
            self.para_nm = list(self.p_emu.columns)
            self.ppe_para = pd.read_csv(self.root / 'tabs/parameters.csv', index_col = 0)
            self.data_obs = pd.read_csv(self.root / 'tabs/obs_data.csv', index_col = 0).iloc[:, 0] #%xx        
            self.data_ppe = pd.read_csv(self.root / 'tabs/ppe_data.csv', index_col = 0)
            self.var_nm = list(self.data_ppe.columns) 

            self.dropped_vars = EmulatedDataStorage()
            
            self.ppe_para_norm = self.ppe_para.copy()
            self.ppe_para_norm = (self.ppe_para_norm - self.ppe_para_norm.min())/(self.ppe_para_norm.max() - self.ppe_para_norm.min())

            self.results = EmulatedDataStorage()
            self.specifications = EmulatedDataStorage()
            
            self.dropped_vars.nooverlap2d = []
        else:
            print("No case created")


    def prepare_case(self, config):
        self.prep_case.sensitivity_emulation(n_cpus=config['n_cpus'])
        self.load_case()
        for level in config['threshold_levels']:
            self.create_mask(threshold_level=level)




    def create_mask(self, threshold_level):
###########################
        mean_paths = glob.glob(str(self.root) + "/y_emu/" + "*mean*", recursive=True)
        tf_masks = []

        for path in mean_paths:
            var_name_file = path.split("/")[-1].split("_mean_std_")[1]
            var_name = var_name_file.split(".")[0]
            emulated_mean_std = pd.read_csv(path,index_col=0)
            emulated_mean = emulated_mean_std.iloc[:,0]
            emulated_std = emulated_mean_std.iloc[:,1]            
            
            obs_temp = self.data_obs.loc[var_name]
            y_ppe = self.data_ppe[var_name]
        
            yscale = y_ppe.std()
            ymu = y_ppe.mean()            
            emulated_mean = emulated_mean * yscale + ymu
            emulated_std = emulated_std * yscale
            
            temp_tf_mask = ((emulated_mean - threshold_level * emulated_std) < obs_temp) & ((emulated_mean + threshold_level * emulated_std) > obs_temp)
            temp_tf_mask.name = var_name
            tf_masks.append(temp_tf_mask)


        tf_masks = pd.concat(tf_masks, axis = 1)
        tf_masks.to_csv(self.root / f"tf_masks_level_{threshold_level}.csv")
        
        
##########################
    def load_mask(self, threshold_level):
        self.tf_masks = pd.read_csv(self.root / f'tf_masks_level_{threshold_level}.csv', index_col=0)
        self.tf_masks_raw = self.tf_masks.copy()
        self.specifications.uncertainty_threshold = threshold_level
        

        self.var_nm = list(self.tf_masks.columns)
        self.n_sample = self.tf_masks.shape[0]

    def visualize_check(self, yname):
        y_emu_norm = pd.read_csv(self.root / f"y_emu/gp_mean_std_{yname}.csv", index_col=0)
        visualize_emulation(X_gcm_norm = self.ppe_para_norm, X_emu = self.p_emu, y_gcm = self.data_ppe[yname], y_emu_norm = y_emu_norm, 
                            para_inds = self.meta[yname], tf_mask = self.tf_masks[yname], 
                            para_nm = self.para_nm, obs = self.data_obs[yname])


    def drop_by_name(self, var_to_exclude):
        var_to_drop = []
        for v in var_to_exclude:
            var_to_drop.extend([s for s in self.var_nm if s.startswith(v)])
        
        self.tf_masks = self.tf_masks.drop(columns = var_to_drop)

        self.var_nm = list(self.tf_masks.columns)
        self.dropped_vars.by_name = var_to_drop
        self.update_meta()
        #xxx%self.specifications.drop_by_name = var_to_exclude

    def drop_by_n_survive(self, n_survive):
        survive_summary = self.tf_masks.sum(axis = 0)
        self.dropped_vars.useless = list(survive_summary[survive_summary == self.n_sample].index)
        self.dropped_vars.tight   = list(survive_summary[survive_summary < n_survive].index)

        self.tf_masks = self.tf_masks.drop(columns = self.dropped_vars.useless + self.dropped_vars.tight)
        self.specifications.n_survive = n_survive
        self.var_nm = list(self.tf_masks.columns)
        self.update_meta()
    
    def drop_by_nvar_per_pair(self, n_var_thre = 1):
        self.dropped_vars.too_few_vars = []
        for k, v in list(self.paras_vars.items()):
            if len(v) <= n_var_thre:
                self.dropped_vars.too_few_vars.extend(v)
                del self.paras_vars[k]
    

        self.tf_masks = self.tf_masks.drop(columns = self.dropped_vars.too_few_vars)
        self.var_nm = list(self.tf_masks.columns)
        self.specifications.n_var_thre_per_parapair = n_var_thre
        self.update_meta()

    def update_meta(self):
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


    def group_para_climatology(self, overlapping_threshold = 10000):

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

            if temp_count < overlapping_threshold:
                paras_vars_0[k] = v
                print(f'{k[0]:<40} and {k[1]:<40}: {temp_count:>8}')

        self.paras_vars_0 = paras_vars_0
        self.specifications.overlapping_threshold = overlapping_threshold
        
        
    def shuffle_vars(self, n_comb = 2):
        summary_table = {}
    
        for paras, vars_temp in self.paras_vars_0.items():
             
            pd_list = []
            
            for vars_comb in combinations(vars_temp, n_comb):
                if len(vars_temp) >= n_comb:
                    pd_list.append(list(vars_comb) + [self.tf_masks[list(vars_comb)].all(axis = 1).sum()])
                else:
                    raise ValueError("n_comb is too large")

            pd_list = pd.DataFrame(pd_list, columns = [f"var{i+1}" for i in range(n_comb)] + ["count"])
            pd_list = pd_list.sort_values(by = "count")
            
            
            summary_table[paras] = pd_list

        print(f'There are {len(self.paras_vars_0)} groups that have no overlapping within own groups')
        if len(self.paras_vars_0) == 0:
            print('No non-overlapping variables at this stage for the variables sharing the same 2 sensitive parameters')

        return summary_table, len(self.paras_vars_0)

    
    def remove_var2d_auto(self, overlapping_threshold, no_iter = 100):
    
        for i in range(no_iter):

            self.group_para_climatology(overlapping_threshold)
            summary2d, no_over_count = self.shuffle_vars()

            if (i == 0) & (no_over_count == 0):
                print('No need to consider non-overlapping')

            if (no_over_count > 0) & (i < no_iter -1) & (i > 0):
                c_s = pd.concat(list(summary2d.values()), axis = 0).sort_values(by='count')
                c_s = c_s[c_s['count'] < overlapping_threshold]
                no_overlap_2d_var = list(c_s[['var1', 'var2']].stack().value_counts()[:1].index) 
                #no_overlap_2d_vars = no_overlap_2d_vars.append(no_overlap_2d_var[0])
                print(f'Drop variable {no_overlap_2d_var}')
                self.drop_no_overlap2d_vars(no_overlap_2d_var)

            if (no_over_count == 0) & (i < no_iter -1) & (i > 0): 
                print('Finished dropping variables')
                return 
            
            if (no_over_count > 0) & (i == no_iter -1):
                print('Failed to resolve non-overlapping')

    


    def drop_no_overlap2d_vars(self, vars_to_drop):
        if not set(vars_to_drop).issubset(list(self.tf_masks.columns)):
            print('Warning: some variables that are proposed to drop are already dropped in previous steps')
        self.tf_masks = self.tf_masks.drop(columns = vars_to_drop, errors = "ignore")
        self.var_nm = list(self.tf_masks.columns)
        self.meta = self.meta[self.var_nm]
        self.meta_onehot = meta_one_hot_shot(self.meta, self.para_nm)
        self.dropped_vars.nooverlap2d.append(vars_to_drop[0])
        #self.specifications.drop_vars_2d = vars_to_drop
        
    
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
    

    def orchestrate(self, n_pts, n_threshold, sample_threshold, max_workers):


        para_seq = list(self.grouped_hulls.keys())
        
        check = orchestrate_test(para_seq, self.p_emu, self.tf_masks,  
                         self.para_nm, self.grouped_hulls, self.paras_vars, n_pts, n_threshold, sample_threshold, max_workers)

        self.results.valid_hulls = check[0]
        self.results.para_l = check[1]
        self.dropped_vars.during_iteration = check[3]
        self.specifications.dropped_during_orchastrate = check[3]


    def prepare_for_sampling(self, shape_alpha = 5, n_pts = 10000, n_threshold = 100, sample_threshold = 10**5, max_workers = 2):
        self.build_hulls(shape_alpha)
        print('Finish constructing the 2d polygons/convex hulls')
        self.orchestrate(n_pts, n_threshold, sample_threshold, max_workers)
        print('Finish preparing for drawing samples')


    def draw(self, n_pts=50000, n_threshold=5000, sample_threshold=10**8, max_workers=32, n_max = 1000):
        valid_hulls = self.results.valid_hulls
        samples = sample_from_hulls_n(list(valid_hulls.keys()), self.para_nm, valid_hulls, n_pts, n_threshold, max_workers, sample_threshold)
        if samples.shape[0]>n_max:
            samples = samples.iloc[:n_max]
            
        self.results.unscaled_samples = samples
        self.results.realscale_samples = self.rescale_para(samples)



    def save_samples_specifications(self, result_name, top_n = 100):

        self.result_name = result_name

        csv_path1 = self.root / ('output/' + result_name + '_all_para_realscale.csv')
        nc_path1 = self.root / ('output/' + result_name + '_all_para_realscale.nc')

        csv_path2 = self.root / ('output/' + result_name + '_topn_para_realscale.csv')
        nc_path2 = self.root / ('output/' + result_name + '_topn_para_realscale.nc')

        if not csv_path1.exists():
            self.results.realscale_samples.to_csv(csv_path1)
            para_csv2nc(csv_path1, nc_path1, self.results.realscale_samples.shape[0])

            self.results.realscale_samples.iloc[:top_n,:].to_csv(csv_path2)
            para_csv2nc(csv_path2, nc_path2, top_n)
        else:
            raise FileExistsError(f'result_name {result_name} already exists')

        self.write_specifications()


    def write_specifications(self):
        def _key_to_str(key):
            if isinstance(key, tuple):
                return "(" + ", ".join(str(k) for k in key) + ")"
            return str(key)

        def _to_jsonable(value):
            if isinstance(value, dict):
                return {_key_to_str(k): _to_jsonable(v) for k, v in value.items()}
            if isinstance(value, (list, tuple, set)):
                return [_to_jsonable(v) for v in value]
            if hasattr(value, "item"):
                # Handles numpy scalar types (e.g., np.int64, np.float64).
                try:
                    return value.item()
                except Exception:
                    pass
            return value

        spec_dict = {
            "uncertainty_threshold": getattr(self.specifications, "uncertainty_threshold", None),
            "n_survive": getattr(self.specifications, "n_survive", None),
            "n_var_thre_per_parapair": getattr(self.specifications, "n_var_thre_per_parapair", None),
            "overlapping_threshold": getattr(self.specifications, "overlapping_threshold", None),
            #"drop_vars_2d": getattr(self.specifications, "drop_vars_2d", None),
            "dropped_during_orchastrate": getattr(self.specifications, "dropped_during_orchastrate", None),
            'para_var' : self.paras_vars,
        }
        spec_dict = _to_jsonable(spec_dict)

        # Python-friendly + human-friendly canonical output
        specificatiaon_path = self.root / ('output/' + self.result_name + "_specifications.json")
        with open(specificatiaon_path, "w", encoding="utf-8") as f:
            json.dump(spec_dict, f, indent=2, sort_keys=True, ensure_ascii=False)

        dropped_vars_dict = _to_jsonable(vars(self.dropped_vars))
        dropped_vars_path = self.root / ('output/' + self.result_name + "_dropped_vars.json")
        with open(dropped_vars_path, "w", encoding="utf-8") as f:
            json.dump(dropped_vars_dict, f, indent=2, sort_keys=True, ensure_ascii=False)




    def compare_with_original(self, df_vline=None, bins=30, density=True):

        dfs = [self.ppe_para, self.results.realscale_samples.iloc[:self.ppe_para.shape[0],:]]

        cols = dfs[0].columns
        ncols = 5
        nrows = (len(cols) + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 3*nrows), squeeze=False)

        df_vline = df_vline 

        for ax, c in zip(axes.ravel(), cols):
            # histograms
            for i, df in enumerate(dfs):
                ax.hist(df[c].dropna(), bins=bins, density=density, alpha=0.4, label=f"hist{i}")
            # vlines
            if df_vline is not None:
                for j, df in enumerate(df_vline):
                    vals = df[c].dropna()
                    for k, v in enumerate(vals):
                        ax.axvline(v, alpha=0.7, lw=1.5, linestyle="--",
                                label=f"vline{j}" if k == 0 else None)
            ax.set_title(c)

        for ax in axes.ravel()[len(cols):]:
            ax.axis("off")
        axes[0,0].legend()
        fig.tight_layout()
