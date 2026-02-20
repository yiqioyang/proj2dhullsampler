import xarray as xr
import pandas as pd
import glob
import os
import math

import numpy as np
import re
from joblib import Parallel, delayed
from pathlib import Path
import matplotlib.pyplot as plt


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
    def y_emu(self):
        return self.root / "y_emu"

    @property
    def tabs(self):
        return self.root / "tabs"

    @property
    def python_obj(self):
        return self.root / "python_obj"

    @property
    def class_obj(self):
        return self.root / "class_obj"


class DataPrep:
    def __init__(self, ppe, obs, obs_dict, para, case: CaseDirectory):
        self.ppe = ppe
        self.obs = obs
        self.obs_dict = obs_dict
        self.para = para
        self.para_nm = list(para.columns)
        para_norm = para.copy()
        para_norm = (para_norm - para_norm.min())/(para_norm.max() - para_norm.min())
        self.para_norm = para_norm
        self.case = case
        if ppe is None:
            print('No PPE input')
        else:
            if np.array_equal(ppe.ppe_ind.to_numpy(), para.index.to_numpy()):
                print("Parameter and simulation indices match")
            else:
                print("Indices not matching!")

    def sample_uniform(self, n = 1000000):
        samples = pd.DataFrame(
            np.random.rand(n, len(self.para_nm)),
            columns=self.para_nm
        )
        xr.Dataset.from_dataframe(samples).to_netcdf(self.case.root / "sampled_parameters.nc")




class FeatureBuilder:
    def __init__(self, ppe, obs, obs_dict, added_data, case: CaseDirectory):
        self.ppe = ppe
        self.obs = obs
        self.obs_dict = obs_dict
        self.case = case  
        self.added_data = added_data
        
    def zonalize_obs_ppe(self, lat_bins, manul_ppe_info):

        
        if self.ppe is None:
            self.added_data[0].to_csv(self.case.tabs / "ppe_tab.csv", index=True)
            self.added_data[1].to_csv(self.case.tabs / "obs_tab.csv", index=True)
    
            self.ppe_pd = self.added_data[0]
            self.obs_pd = self.added_data[1]
            self.var_nm = list(ppe_pd.columns)

        else:
            ppe_zonal_list = []
            obs_zonal_list = []
    
            
            lab_bin_labels = np.char.add(np.char.add(lat_bins[:-1].astype(str), "to"), lat_bins[1:].astype(str))
            lab_bin_labels = np.char.add("zonal_", lab_bin_labels)
            ############################################################################################################
            for cam_nm, obs_nm in self.obs_dict.items():
            
                ppe_da = self.ppe[cam_nm]
                filter_tf = self.obs[obs_nm].notnull()           ## ## Take out the na values that are in obs from the PPE 
                ppe_da = ppe_da.where(filter_tf)
    
    
                
                zonal_ppe_temp = (ppe_da.mean(dim  = "lon", 
                                              skipna = True).groupby_bins("lat",lat_bins, labels = lab_bin_labels).mean(dim = "lat", skipna = True).to_dataframe().unstack(level = 1))
                
                zonal_ppe_temp.columns.name = None # At this point, zonal_ppe_temp has two level in the columns
                zonal_ppe_temp.columns = ["_".join(col) for col in list(zonal_ppe_temp.columns)] # The comprehension unpack the two levels
                
            
                zonal_obs_temp = self.obs[obs_nm].mean(dim = "lon", skipna = True).groupby_bins("lat",lat_bins, labels = lab_bin_labels).mean(dim = "lat", skipna = True).to_series()
                zonal_obs_temp.index = zonal_ppe_temp.columns
            
                ppe_zonal_list.append(zonal_ppe_temp)
                obs_zonal_list.append(zonal_obs_temp)
            
            ppe_zonal_pd = pd.concat(ppe_zonal_list, axis = 1)
            obs_zonal_pd = pd.concat(obs_zonal_list)
    
            nan_obs_vars = list(obs_zonal_pd.index[obs_zonal_pd.isna()])
            obs_zonal_pd = obs_zonal_pd.dropna()
            
            nan_ppe_vars = list(ppe_zonal_pd.columns[ppe_zonal_pd.isna().any()])
            ppe_zonal_pd = ppe_zonal_pd.dropna(axis = 1)
    
            if sorted(nan_obs_vars) == sorted(nan_ppe_vars):
                print("nan variables matching between obs and simulation")
            else:
                print("non variables not matching between obs and simulation")
            ############################################################################################################
            
            ppe_manual_list = []
            obs_manual_list = []
            manual_name_list = []
        
            for row_ind, row in manul_ppe_info.iterrows():
                
                temp_obs = self.obs[self.obs_dict[row.nm]].sel(lat = slice(row.lat_min, row.lat_max), lon = slice(row.lon_min, row.lon_max)).mean(dim = ["lat", "lon"]).values
                if ~np.isnan(temp_obs):
                    temp_ppe = self.ppe[row.nm].sel(lat = slice(row.lat_min, row.lat_max), lon = slice(row.lon_min, row.lon_max)).mean(dim = ["lat", "lon"]).to_dataframe()
                    
                    manual_name_list.append("_".join(row.astype(str)))
                    ppe_manual_list.append(temp_ppe)
                    obs_manual_list.append(temp_obs)
                else:
                    print("??")        
        
        
            ppe_manual_pd = pd.concat(ppe_manual_list,axis = 1)
            obs_manual_pd = pd.Series(obs_manual_list)
        
            ppe_manual_pd.columns = manual_name_list
            obs_manual_pd.index = manual_name_list
            ############################################################################################################
            if self.added_data is not None:
                if np.array_equal(self.added_data[0].index.to_numpy(), ppe_zonal_pd.index.to_numpy()):
                    print("Added data index matching")
                    ppe_pd = pd.concat([ppe_zonal_pd, ppe_manual_pd, self.added_data[0]], axis = 1)
                    obs_pd = pd.concat([obs_zonal_pd, obs_manual_pd, self.added_data[1]])
    
                else:
                    print('Added data index not matching, break')
                    return 
                    
            else:
                ppe_pd = pd.concat([ppe_zonal_pd, ppe_manual_pd], axis = 1)
                obs_pd = pd.concat([obs_zonal_pd, obs_manual_pd])
                
            #ppe_pd.to_csv(os.path.join(self.path, "tabs/", "ppe_tab.csv"), index = True)
            #obs_pd.to_csv(os.path.join(self.path, "tabs/", "obs_tab.csv"), index=True)
    
            ppe_pd.to_csv(self.case.tabs / "ppe_tab.csv", index=True)
            obs_pd.to_csv(self.case.tabs / "obs_tab.csv", index=True)
    
            print("Zonalized and manually selected obs and ppe written as csv")
            self.ppe_pd = ppe_pd
            self.obs_pd = obs_pd
            self.var_nm = list(ppe_pd.columns)


def visualize_emulation(X_gcm_norm, X_emu, y_gcm, y_emu_norm, para_inds, tf_mask, para_nm, obs, yname):

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


class Prep_Mask_Generation:
    def __init__(self, working_dir, case_name, ppe, obs, obs_dict, para, lat_bins, manul_ppe_info, added_ppe_obs = None, n_sample = 1000000):
        
        self.case = CaseDirectory(working_dir, case_name)
        self.data_gcm = DataPrep(ppe, obs, obs_dict, para, self.case)
        self.data_gcm.sample_uniform(n_sample)
        
        self.features = FeatureBuilder(ppe, obs, obs_dict, added_ppe_obs, self.case)
        self.features.zonalize_obs_ppe(lat_bins, manul_ppe_info)

        
        self.data_gcm.ppe_pd = self.features.ppe_pd
        self.data_gcm.obs_pd = self.features.obs_pd
        self.data_gcm.var_nm = self.features.var_nm
        
        self.data_emu = EmulatedDataStorage()
        self.var_excluded = EmulatedDataStorage()
    
    def load_certainy(self, yname):
        return(pd.read_csv(self.case.root / f"y_emu/gp_mean_std_{yname}.csv", index_col=0))

    def load_para_emu(self):
        return xr.open_dataset(self.case.root / "sampled_parameters.nc").to_dataframe()
    
    def sensitivity_emulation(self, n_sens_p = 2, n_cpus = 15):
        from funs.utils import gp_training_application, fit_gp_for_single_1d, fit_all_gp_models_1d
        
        para_s = xr.open_dataset(self.case.root / "sampled_parameters.nc").to_dataframe()
        
        
        results = Parallel(n_jobs=n_cpus)(
                        delayed(gp_training_application)(self.data_gcm.para_norm, self.data_gcm.ppe_pd, y_name, para_s, path = str(self.case.root) + "/", n_sens_p=n_sens_p)
                        for y_name in list(self.data_gcm.ppe_pd.columns)
                    )

        del para_s

        meta_xy_dict = {pair[0]: pd.Series(pair[1]) for pair in results if pair is not None}
        meta = pd.concat(list(meta_xy_dict.values()), axis = 1)
        meta.columns = list(meta_xy_dict.keys())
        self.meta = meta
        
        meta.to_csv(self.case.root / "meta.csv", index=True)
        
    
    def mask_generation(self, threshold_level = 2.0):
        mean_paths = glob.glob(str(self.case.root) + "/y_emu/" + "*mean*", recursive=True)
        tf_masks = []

        for path in mean_paths:
            var_name_file = path.split("/")[-1].split("_mean_std_")[1]
            if "zonal" in var_name_file:
                cli_temp = var_name_file.split("_zonal")[0]
                
            else:
                cli_temp = var_name_file.split("_")[0]  #% xxx
                    
            var_name = var_name_file.split(".")[0]
            
            emulated_mean_std = pd.read_csv(path,index_col=0)
            emulated_mean = emulated_mean_std.iloc[:,0]
            emulated_std = emulated_mean_std.iloc[:,1]            
            
            obs_temp = self.data_gcm.obs_pd.loc[var_name]
            y_ppe = self.data_gcm.ppe_pd[var_name]
        
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
        self.data_emu.tf_masks = tf_masks

    def visualize_check(self, yname):
        
        y_emu_norm = self.load_certainy(yname)
        X_emu = self.load_para_emu()
        
        visualize_emulation(X_gcm_norm = self.data_gcm.para_norm, X_emu = X_emu, y_gcm = self.data_gcm.ppe_pd[yname], y_emu_norm = y_emu_norm, 
                            para_inds = self.meta[yname], tf_mask = self.data_emu.tf_masks[yname], 
                            para_nm = self.data_gcm.para_nm, obs = self.data_gcm.obs_pd[yname], yname = yname)



