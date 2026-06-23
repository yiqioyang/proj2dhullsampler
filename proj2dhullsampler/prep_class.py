import xarray as xr
import pandas as pd
import glob

import numpy as np
from joblib import Parallel, delayed
from pathlib import Path
import matplotlib.pyplot as plt

from .preprocess import feature_builder
from .utils import gp_training_application

class EmulatedDataStorage:
    """
    Lightweight container for emulation outputs.
    No computation logic.
    """
    def __init__(self):
        pass



class CaseDirectory:
    def __init__(self, working_dir, case_name):
        self.root = Path(working_dir) / case_name
        self._init_dirs()

    def _init_dirs(self):
        self.root.mkdir(parents=True, exist_ok=True)
        (self.root / "y_emu").mkdir(exist_ok=True)
        (self.root / "tabs").mkdir(exist_ok=True)
        (self.root / "python_obj").mkdir(exist_ok=True)
        (self.root / "class_obj").mkdir(exist_ok=True)

    @property
    def path_y_emu(self):
        return self.root / "y_emu"

    @property
    def path_tabs(self):
        return self.root / "tabs"

    @property
    def path_python_obj(self):
        return self.root / "python_obj"

    @property
    def path_class_obj(self):
        return self.root / "class_obj"



def visualize_emulation(X_gcm_norm, X_emu, y_gcm, y_emu_norm, para_inds, tf_mask, para_nm, obs):

    y_emu_norm.iloc[:,0] = y_emu_norm.iloc[:,0] * y_gcm.std() + y_gcm.mean()
    y_emu_norm.iloc[:,1] = y_emu_norm.iloc[:,1] * y_gcm.std()
    
    xy_emu = pd.concat([X_emu, y_emu_norm], axis = 1)
    xy_emu_sub = xy_emu[tf_mask]

    xy_emu = xy_emu.sample(50000)
    if xy_emu_sub.shape[0] > 50000:
        xy_emu_sub = xy_emu_sub.sample(50000)
            


    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))  # 1 row, 2 columns
    
    xy_emu.sort_values(by = para_nm[para_inds[0]])
    xy_emu_sub.sort_values(by =para_nm[para_inds[0]])

    ax1.scatter(xy_emu.iloc[:, para_inds[0]], xy_emu.iloc[:, -2])
    ax1.scatter(xy_emu_sub.iloc[:, para_inds[0]], xy_emu_sub.iloc[:, -2])
    # ax1.plot(xy_emu.iloc[:, para_inds[0]], xy_emu.iloc[:, -2] - xy_emu.iloc[:, -1], color = 'gray')
    # ax1.plot(xy_emu.iloc[:, para_inds[0]], xy_emu.iloc[:, -2] + xy_emu.iloc[:, -1], color = 'gray')
    
    ax1.scatter(X_gcm_norm.iloc[:,para_inds[0]], y_gcm)
    ax1.axhline(obs)
    ax1.set_xlabel(para_nm[para_inds[0]])
#############################################################################
    xy_emu.sort_values(by = para_nm[para_inds[1]])
    xy_emu_sub.sort_values(by =para_nm[para_inds[1]])

    ax2.scatter(xy_emu.iloc[:, para_inds[1]], xy_emu.iloc[:, -2])
    ax2.scatter(xy_emu_sub.iloc[:, para_inds[1]], xy_emu_sub.iloc[:, -2])
    
    ax2.scatter(X_gcm_norm.iloc[:,para_inds[1]], y_gcm)
    ax2.axhline(obs)
    ax2.set_xlabel(para_nm[para_inds[1]])
    plt.show()


def meta_one_hot_shot(meta, para_nm):
    meta = meta.transpose()
    meta_one_hot = pd.DataFrame(False, index=meta.index, columns=para_nm)
    for index, row in meta.iterrows():
        for r in row.values:
            meta_one_hot.loc[index, para_nm[r]] = True

    return meta_one_hot


class Prepare_Case:
    def __init__(self, working_dir, case_name, para, tabs, ppe, obs, obs_dict, lat_bins, manul_ppe_info, n_sample = 1000000):
        
        self.wd = working_dir
        self.case_name = case_name
        self.case = CaseDirectory(working_dir, case_name)
        self.n_sample = n_sample

        #### Process ppe data
        ppe_data, obs_data = feature_builder(tabs, ppe, obs, obs_dict, lat_bins, manul_ppe_info)
        
        if para.index.equals(ppe_data.index):
            self.data_gcm = EmulatedDataStorage()
            para_norm = para.copy()
            para_norm = (para_norm - para_norm.min())/(para_norm.max() - para_norm.min())
            self.data_gcm.para = para
            self.data_gcm.para_norm = para_norm
            self.data_gcm.ppe_data = ppe_data
            self.data_gcm.obs_data = obs_data
            self.data_gcm.var_nm = list(ppe_data.columns)
            self.data_gcm.para_nm = list(para.columns)
        else:
            raise ValueError("Parameters and simulation output indices do not match")

        para.to_csv(self.case.path_tabs / 'parameters.csv')
        ppe_data.to_csv(self.case.path_tabs / 'ppe_data.csv')
        obs_data.to_csv(self.case.path_tabs / 'obs_data.csv')
        

        ##### Sample parameters
        self.sample_uniform(n_sample)

    def sample_uniform(self, n):
        samples = pd.DataFrame(np.random.rand(n, len(self.data_gcm.para_nm)),
                                columns=self.data_gcm.para_nm
                                )
        xr.Dataset.from_dataframe(samples).to_netcdf(self.case.root / "sampled_parameters.nc")


    def sensitivity_emulation(self, n_sens_p = 2, n_cpus = 15):
        
        sampled_paras = xr.open_dataset(self.case.root / "sampled_parameters.nc").to_dataframe()
        
        results = Parallel(n_jobs=n_cpus)(
                        delayed(gp_training_application)(self.data_gcm.para_norm, self.data_gcm.ppe_data, y_name, sampled_paras, path = str(self.case.root) + "/", n_sens_p=n_sens_p)
                        for y_name in list(self.data_gcm.ppe_data.columns)
                    )

        del sampled_paras

        meta_xy_dict = {pair[0]: pd.Series(pair[1]) for pair in results if pair is not None}
        meta = pd.concat(list(meta_xy_dict.values()), axis = 1)
        meta.columns = list(meta_xy_dict.keys())
        self.meta = meta
        
        meta.to_csv(self.case.root / "meta.csv", index=True)
        
        # mark
    
    def mask_generation(self, threshold_level = 2.0):
        mean_paths = glob.glob(str(self.case.root) + "/y_emu/" + "*mean*", recursive=True)
        tf_masks = []

        for path in mean_paths:
            var_name_file = path.split("/")[-1].split("_mean_std_")[1]

                    
            var_name = var_name_file.split(".")[0]
            
            emulated_mean_std = pd.read_csv(path,index_col=0)
            emulated_mean = emulated_mean_std.iloc[:,0]
            emulated_std = emulated_mean_std.iloc[:,1]            
            
            obs_temp = self.data_gcm.obs_data.loc[var_name]
            y_ppe = self.data_gcm.ppe_data[var_name]
        
            yscale = y_ppe.std()
            ymu = y_ppe.mean()            
            emulated_mean = emulated_mean * yscale + ymu
            emulated_std = emulated_std * yscale
            
            temp_tf_mask = ((emulated_mean - threshold_level * emulated_std) < obs_temp) & ((emulated_mean + threshold_level * emulated_std) > obs_temp)
            temp_tf_mask.name = var_name
            tf_masks.append(temp_tf_mask)


        tf_masks = pd.concat(tf_masks, axis = 1)
        tf_masks.to_csv(self.case.root / f"tf_masks_level_{threshold_level}.csv")
        self.threshold_level = threshold_level
        

### Add function to update parameters such that they match
    def load_mask(self, threshold_level):
        return(pd.read_csv(self.case.root / f"tf_masks_level_{threshold_level}.csv", index_col = 0))

    def load_certainy(self, yname):
        return(pd.read_csv(self.case.root / f"y_emu/gp_mean_std_{yname}.csv", index_col=0))

    def load_para_emu(self):
        return xr.open_dataset(self.case.root / "sampled_parameters.nc").to_dataframe()

    # def visualize_check(self, yname, threshold):
    #     y_emu_norm = pd.read_csv(self.case.root / f"y_emu/gp_mean_std_{yname}.csv", index_col=0)
    #     X_emu = self.load_para_emu()
    #     tf_masks = self.load_mask(threshold)[yname]

    #     visualize_emulation(X_gcm_norm = self.data_gcm.para_norm, X_emu = X_emu, y_gcm = self.data_gcm.ppe_data[yname], y_emu_norm = y_emu_norm, 
    #                         para_inds = self.meta[yname], tf_mask = tf_masks, 
    #                         para_nm = self.data_gcm.para_nm, obs = self.data_gcm.obs_data[yname])



