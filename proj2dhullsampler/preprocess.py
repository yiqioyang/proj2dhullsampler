import numpy as np
import pandas as pd
import xarray as xr

def zonal_process(ppe, obs, obs_dict, lat_bins):
    ppe_zonal_list = []
    obs_zonal_list = []

    lab_bin_labels = np.char.add(np.char.add(lat_bins[:-1].astype(str), "to"), lat_bins[1:].astype(str))
    lab_bin_labels = np.char.add("zonal_", lab_bin_labels)
    ############################################################################################################
    for cam_nm, obs_nm in obs_dict.items():

        ppe_da = ppe[cam_nm]
        filter_tf = obs[obs_nm].notnull()           ## ## Take out the na values that are in obs from the PPE 
        ppe_da = ppe_da.where(filter_tf)

        
        zonal_ppe_temp = (ppe_da.mean(dim  = "lon", 
                                        skipna = True).groupby_bins("lat",lat_bins, labels = lab_bin_labels).mean(dim = "lat", skipna = True).to_dataframe().unstack(level = 1))
        
        zonal_ppe_temp.columns.name = None # At this point, zonal_ppe_temp has two level in the columns
        zonal_ppe_temp.columns = ["_".join(col) for col in list(zonal_ppe_temp.columns)] # The comprehension unpack the two levels
        

        zonal_obs_temp = obs[obs_nm].mean(dim = "lon", skipna = True).groupby_bins("lat",lat_bins, labels = lab_bin_labels).mean(dim = "lat", skipna = True).to_series()
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
        print("Check on nan variables passes: nan variables matching between obs and simulation")
    else:
        print("Check on nan variables FAILED: nan variables do not match between obs and simulation")

    print(f'obs_shape: {obs_zonal_pd.shape}, and ppe shape: {ppe_zonal_pd.shape}')

    return ppe_zonal_pd, obs_zonal_pd


def local_process(ppe, obs, obs_dict, manul_ppe_info):

    ppe_manual_list = []
    obs_manual_list = []
    manual_name_list = []

    for row_ind, row in manul_ppe_info.iterrows():
        temp_obs = obs[obs_dict[row.nm]].sel(lat = slice(row.lat_min, row.lat_max), lon = slice(row.lon_min, row.lon_max)).mean(dim = ["lat", "lon"], skipna=False).values
        if ~np.isnan(temp_obs):
            temp_ppe = ppe[row.nm].sel(lat = slice(row.lat_min, row.lat_max), lon = slice(row.lon_min, row.lon_max)).mean(dim = ["lat", "lon"], skipna=False).to_dataframe()
            
            manual_name_list.append("_".join(row.astype(str)))
            ppe_manual_list.append(temp_ppe)
            obs_manual_list.append(temp_obs)
        else:
            print("There are NaN values in the selected area for the obs")        


    manual_name_list = ["local_" + name for name in manual_name_list]
    ppe_manual_pd = pd.concat(ppe_manual_list,axis = 1)
    obs_manual_pd = pd.Series(obs_manual_list)

    ppe_manual_pd.columns = manual_name_list
    obs_manual_pd.index = manual_name_list

    return ppe_manual_pd, obs_manual_pd
    ############################################################################################################




def feature_builder(tabs, ppe = None, obs = None, obs_dict = None, lat_bins = None, manul_ppe_info = None):
    if ppe is None:
        ppe_tab, obs_tab = tabs
        if obs_tab.index.equals(ppe_tab.columns):
            print("obs and ppe variable names match")
            return ppe_tab, obs_tab
        else:
            print("obs and ppe variable names do not match")
            return 

    else:
        if manul_ppe_info is not None and lat_bins is not None:
            ppe_manual_pd, obs_manual_pd = local_process(ppe, obs, obs_dict, manul_ppe_info)
            ppe_zonal_pd, obs_zonal_pd = zonal_process(ppe, obs, obs_dict, lat_bins)
                        
            if ppe_manual_pd.index.equals(ppe_zonal_pd.index):
                print("Manual and zonal indices match")

            else:
                print("Manual and zonal indices NOT match")
                return
            
            ppe_zonal_manual = pd.concat([ppe_zonal_pd, ppe_manual_pd], axis = 1)
            obs_zonal_manual = pd.concat([obs_zonal_pd, obs_manual_pd])
            
            
        if manul_ppe_info is None and lat_bins is not None:
            ppe_zonal_manual, obs_zonal_manual = zonal_process(ppe, obs, obs_dict, lat_bins)
            
        if manul_ppe_info is not None and lat_bins is None:
            ppe_zonal_manual, obs_zonal_manual = local_process(ppe, obs, obs_dict, manul_ppe_info)

        if manul_ppe_info is None and lat_bins is None:
            raise ValueError("Need to specify lat_bins or manul_ppe_info")

        if tabs is None:
            if obs_zonal_manual.index.equals(ppe_zonal_manual.columns):
                print("obs and ppe variable names match")
                return ppe_zonal_manual, obs_zonal_manual
            else:
                print("obs and ppe variable names do not match")
                return 
        else:
            ppe_tab, obs_tab = tabs
            if ppe_tab.index.equals(ppe_zonal_manual.index):
                print("Tabulated and processed data indices match")

                ppe_combine = pd.concat([ppe_tab, ppe_zonal_manual], axis = 1)
                obs_combine = pd.concat([obs_tab, obs_zonal_manual])

                if obs_combine.index.equals(ppe_combine.columns):
                    print("obs and ppe variable names match")
                    return ppe_combine, obs_combine
                else:
                    print("obs and ppe variable names do not match")
                    return 

            else:
                print("Tabulated and processed data indices do not match")
                return



