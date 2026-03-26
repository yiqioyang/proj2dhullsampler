#from pyDOE import *
import xarray as xr

import numpy as np
import pandas as pd
#from pyDOE2 import lhs # PyDOE2 is no more supported

#from netCDF4 import Dataset


def para_csv2nc(path_to_csv, output_path, sample_n):
    df = pd.read_csv(path_to_csv, index_col=0)
    random_array = df.to_numpy()
    number_sample = random_array
    print(random_array.shape)

        
    var_arr=list(df.columns)
    
    nmb_var=len(var_arr)
    
    
    #design1  = random_array
    #name_arr = var_arr
    
    samp=sample_n
    number_sample = []

    for k in range(samp):
        number_sample.append("{0:03}".format(k+1))
    print(number_sample)
    
    design1 = np.empty((samp,nmb_var))
    
    name_arr = []
    count=0
    for i in range(nmb_var):
        count=count+1
        parname=var_arr[i]
        #print(parname)
        name_arr.append(parname)
        design1[:,i]= random_array[:,i]
    

    
    print(type(design1))
    print(type(name_arr))
    ds=xr.Dataset({name_arr[0]:(["nmb_sim"], design1[:,0])})
    #ds=ds.assign_coords({"sample_nmb": number_sample}) 
       
    for n in range(nmb_var-1):
        if name_arr[n+1] == 'zmconv_num_cin':
            designh=np.rint(design1[:,n+1])
            designh = designh.astype(int)
    #       ds=xr.merge([ds,xr.Dataset({name_arr[n+1]:(["nmb_sim"], np.rint(design1[:,n+1]))})])
            ds=xr.merge([ds,xr.Dataset({name_arr[n+1]:(["nmb_sim"], designh)})])
        else:
            ds=xr.merge([ds,xr.Dataset({name_arr[n+1]:(["nmb_sim"], design1[:,n+1])})])
    
    ds=xr.merge([ds,xr.Dataset({"Sample_nmb":(["nmb_sim"], number_sample)})])
    
    
    ds.to_netcdf(path=output_path)
    return



def metric_cal_single(obs_dict, obs, ppe):
    weights = np.cos(np.deg2rad(obs.lat))
    weights = weights / weights.mean()

    ave = []
    bias = []
    rmse = []

    
    
    for k, v in obs_dict.items():
        
        ave.append((ppe[k] - obs[v]).weighted(weights).mean(dim = ["lat", "lon"]).rename(k).to_dataframe())
        bias.append((ppe[k] - obs[v]).weighted(weights).mean(dim = ["lat", "lon"]).rename(k).to_dataframe())
        rmse.append(((ppe[k] - obs[v])**2).weighted(weights).mean(dim = ["lat", "lon"]).rename(k).to_dataframe())


    ave = pd.concat(ave, axis = 1)
    bias = pd.concat(bias, axis = 1)
    rmse = pd.concat(rmse, axis = 1)

    return ave, bias, rmse


#metrics = metric_cal_single(obs_dict, obs_ds, ppev1)



    
    




    
