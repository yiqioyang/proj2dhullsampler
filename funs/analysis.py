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
import alphashape
from itertools import combinations
from collections import defaultdict, deque
from concurrent.futures import ThreadPoolExecutor

from concurrent.futures import ProcessPoolExecutor, as_completed

from shapely import points, contains
import random



from prep_class import (meta_one_hot_shot,
visualize_emulation)

class EmulatedDataStorage:
    """
    Lightweight container for emulation outputs.
    No computation logic.
    """
    def __init__(self):
        pass



class Analysis:
    def __init__(self, working_dir, case_name, ppe_para, threshold_level = 2.0):
        self.root = Path(working_dir) / case_name

        self.tf_masks = pd.read_csv(self.root / f'tf_masks_level_{threshold_level}.csv', index_col=0)
        self.meta = pd.read_csv(self.root / 'meta.csv', index_col = 0)
        
        self.p_emu = xr.open_dataset(self.root / 'sampled_parameters.nc').to_dataframe()
        
        self.var_nm = list(self.tf_masks.columns)
        self.para_nm = list(self.p_emu.columns)

        self.data = EmulatedDataStorage()
        self.data.pd_ppe = pd.read_csv(self.root / 'tabs/ppe_tab.csv', index_col = 0)
        self.data.pd_obs = pd.read_csv(self.root / 'tabs/obs_tab.csv', index_col = 0)
        
        self.data.ppe_para = ppe_para
        ppe_para_norm = ppe_para.copy()
        
        self.data.ppe_para_norm = (ppe_para_norm - ppe_para_norm.min())/(ppe_para_norm.max() - ppe_para_norm.min())
        
        self.meta_onehot = meta_one_hot_shot(self.meta, self.para_nm)

    def load_certainy(self, yname):
        return(pd.read_csv(self.root / f"y_emu/gp_mean_std_{yname}.csv", index_col=0))


    
    def visualize_check(self, yname):
        
        y_emu_norm = self.load_certainy(yname)
        X_emu = self.p_emu
        
        visualize_emulation(X_gcm_norm = self.data.ppe_para_norm, X_emu = X_emu, y_gcm = self.data.pd_ppe[yname], y_emu_norm = y_emu_norm, 
                            para_inds = self.meta[yname], tf_mask = self.tf_masks[yname], 
                            para_nm = self.para_nm, obs = self.data.pd_obs.loc[yname].values, yname = yname)
