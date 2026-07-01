import xarray as xr
import pandas as pd


import numpy as np
from joblib import Parallel, delayed
from pathlib import Path


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
        (self.root / "output").mkdir(exist_ok=True)


    @property
    def path_tabs(self):
        return self.root / "tabs"

    



class Prepare_Case:
    def __init__(self, working_dir, case_name, para, tabs, ppe, obs, obs_dict, lat_bins, manul_ppe_info, n_sample = 1000000):
        

        para_range = para.max() - para.min()
        constant_parameters = para_range[
            para_range == 0
        ].index.tolist()

        if constant_parameters:
            raise ValueError(
                f"Constant parameters are not allowed: "
                f"{constant_parameters}"
            )
        

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
        
        

