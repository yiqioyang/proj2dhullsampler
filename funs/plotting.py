
import matplotlib.pyplot as plt

def biplot(dfs, subsample_size = 5000, figsize = (30, 30)):
    ''''
    Purpose: visualize the sampled/PPE/CPE NORMALIZED parameters from biplots
    Inputs: 
        dfs:                A list of pd dataframes or a single pd dataframe
        subsample_size:     Randomly subset samples to this size if the number of 
                            samples (# of rows) within the pd dataframes is above this value. 
    Note:
        x and y axes are constrained to range from 0-1

                            
    ''''
    if isinstance(dfs, pd.DataFrame):
        dfs = [dfs]

    fig, axes = plt.subplots(nrows=dfs[0].shape[1], ncols=dfs[0].shape[1], figsize=(30, 30))
    plt.subplots_adjust(wspace=0.05, hspace=0.05)

    
    for df in dfs:
        if df.shape[0] > subsample_size
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

    ''''
    Purpose: visualize the sampled/PPE/CPE ORIGINAL-SCALED parameters from biplots
    Inputs: 
        dfs:                A list of pd dataframes or a single pd dataframe
        subsample_size:     Randomly subset samples to this size if the number of 
                            samples (# of rows) within the pd dataframes is above this value. 

                            
    ''''

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






