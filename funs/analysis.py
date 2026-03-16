import xarray as xr
import pandas as pd
import math

import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns


from .prep_class import meta_one_hot_shot, visualize_emulation

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


    def plot_by_para(self, para, index = None,ncols=4, figsize=(12, 8), sharex=True, sharey=False):

        columns = sorted(list(self.meta_onehot.index[self.meta_onehot[para]]))

        
        n = len(columns)
        nrows = math.ceil(n / ncols)
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=sharex, sharey=sharey)
        axes = axes.ravel() if n > 1 else [axes]

        x_val = self.data.ppe_para[para]
        
        for i, col in enumerate(columns):
            axes[i].scatter(x_val, self.data.pd_ppe[col])
            axes[i].set_xlabel(para)
            if index is not None:
                axes[i].scatter(x_val.loc[index], self.data.pd_ppe[col].loc[index])
            axes[i].set_ylabel(col)
            axes[i].axhline(self.data.pd_obs.loc[col].values)
            
    
        # Hide any unused axes
        for j in range(i + 1, len(axes)):
            axes[j].axis("off")
    
        fig.tight_layout();
        return fig, axes    

    def plot_onehot(self, figsize = (12, 8), dpi = 120, xstep = 3):
    
        data = self.meta_onehot.T.astype(int)
        
        plt.figure(figsize=(12, 8), dpi=120)
        ax = sns.heatmap(
            data,
            cmap="viridis",
            cbar_kws={"label": "False=0, True=1"},
            yticklabels=True,
            xticklabels=False
        )
        
        # show every Nth x label

        x_labels = self.meta_onehot.index[::xstep]
        x_pos = np.arange(0, data.shape[1], xstep) + 0.5  # center ticks on cells
        
        ax.set_xticks(x_pos)
        ax.set_xticklabels(x_labels, rotation=90, fontsize=8)
        
        ax.set_xlabel("")
        ax.set_ylabel("")
        plt.tight_layout()
        plt.show()
        