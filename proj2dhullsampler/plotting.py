import pandas as pd

import matplotlib.pyplot as plt

def biplot(dfs, subsample_size = 5000, figsize = (30, 30)):
    """
    Purpose: visualize the sampled/PPE/CPE NORMALIZED parameters from biplots
    Inputs: 
        dfs:                A list of pd dataframes or a single pd dataframe
        subsample_size:     Randomly subset samples to this size if the number of 
                            samples (# of rows) within the pd dataframes is above this value. 
    Note:
        x and y axes are constrained to range from 0-1
                    
    """
    if isinstance(dfs, pd.DataFrame):
        dfs = [dfs]

    fig, axes = plt.subplots(nrows=dfs[0].shape[1], ncols=dfs[0].shape[1], figsize=(30, 30))
    plt.subplots_adjust(wspace=0.05, hspace=0.05)

    
    for df in dfs:
        if df.shape[0] > subsample_size:
            df = df.sample(subsample_size)
        
        for i, row_col in enumerate(df.columns):
            for j, col_col in enumerate(df.columns):
                ax = axes[i, j]
        
                if i == j:
                    # Plot histogram on the diagonal
                    ax.hist(df[row_col], bins=10, alpha=0.3, density = True)
                    
                    ax.set_xlim(0,1)
                elif i > j:
                    # Plot scatter off-diagonal
                    ax.scatter(df[col_col], df[row_col], s=5, alpha=0.2)
        
                    ax.set_xlim(0, 1)
                    ax.set_ylim(0, 1)
                    
                if i == df.shape[1] - 1:
                    ax.set_xlabel(col_col, fontsize=6)
                else:
                    ax.set_xticks([])
        
                if j == 0:
                    ax.set_ylabel(row_col, fontsize=6)
                else:
                    ax.set_yticks([])
    
    plt.suptitle("Pairwise parameters (0-1 normalized)", fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.98])
    plt.show()



def biplot_original_scale(dfs, subsample_size = 5000, figsize = (50, 50)):

    '''
    Purpose: visualize the sampled/PPE/CPE ORIGINAL-SCALED parameters from biplots
    Inputs: 
        dfs:                A list of pd dataframes or a single pd dataframe
        subsample_size:     Randomly subset samples to this size if the number of 
                            samples (# of rows) within the pd dataframes is above this value. 

                            
    '''

    if isinstance(dfs, pd.DataFrame):
        dfs = [dfs]

    fig, axes = plt.subplots(nrows=dfs[0].shape[1], ncols=dfs[0].shape[1], figsize=(30, 30))
    plt.subplots_adjust(wspace=0.05, hspace=0.05)

    for df in dfs:
        
        if df.shape[0] > subsample_size:
            df = df.sample(subsample_size)
        
        for i, row_col in enumerate(df.columns):
            for j, col_col in enumerate(df.columns):
                ax = axes[i, j]
        
                if i == j:
                    # Plot histogram on the diagonal
                    ax.hist(df[row_col], bins=10, alpha=0.3, density = True)
                    
                elif i > j:
                    # Plot scatter off-diagonal
                    ax.scatter(df[col_col], df[row_col], s=5, alpha=0.2)
                    
                if i == df.shape[1] - 1:
                    ax.set_xlabel(col_col, fontsize=6)
                    ax.set_xticks([])
                else:
                    ax.set_xticks([])
        
                if j == 0:
                    ax.set_ylabel(row_col, fontsize=6)
                    ax.set_yticks([])
                else:
                    ax.set_yticks([])
    
            
    plt.suptitle("Pairwise parameters (original scaled)", fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.98])
    plt.show()


def plot_histograms_grid_5(dfs, df_vline=None, bins=30, density=False):
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

