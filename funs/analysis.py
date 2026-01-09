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




class EmulatedDataStorage:
    """
    Lightweight container for emulation outputs.
    No computation logic.
    """
    def __init__(self):
        pass



class Analysis:
    def __init__(self, path, ppe_para, case, threshold_level = 2.0):
        self.root = Path(working_dir) / case_name

        self.tf_masks = pd.read_csv(self.root / f'tf_masks_level_{threshold_level}.csv', index_col=0)
        self.meta = pd.read_csv(self.root / 'meta.csv', index_col = 0)
        
        self.p_emu = xr.open_dataset(self.root / 'sampled_parameters.nc').to_dataframe()
        
        self.var_nm = list(self.tf_masks.columns)
        self.para_nm = list(self.p_emu.columns)

        self.data = EmulatedDataStorage()
        self.data.ppe_para = ppe_para
        self.data.pd_ppe = pd.read_csv(self.root / 'tabs/ppe_tab.csv', index_col = 0)
        self.data.pd_obs = pd.read_csv(self.root / 'tabs/obs_tab.csv', index_col = 0)
        
    
