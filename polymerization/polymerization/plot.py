import matplotlib.pyplot as plt
import numpy as np

def breakpoint2scatter(
    df_breakpoint,
    outfile,
    title=None,
    facet_by_strand_config=False,
    jitter=0,
    size = 3,
    alpha = 0.7
):
    '''
    scatter breakpoint positions of left and right genes from breakpoint dataframe
    '''
    if jitter:
        df_breakpoint['left_breakpoint'] += np.random.normal(0, jitter, size=len(df_breakpoint))
        df_breakpoint['right_breakpoint'] += np.random.normal(0, jitter, size=len(df_breakpoint))
    if facet_by_strand_config:
        df_breakpoint['strand_config'] = df_breakpoint.apply(lambda row: f"{row['left_strand']},{row['right_strand']}", axis=1)
        color_map = {
        '1,1': 'blue',
        '1,-1': 'red',
        '-1,1': 'green',
        '-1,-1': 'purple'
        }
        marker_map = {
            '1,1': 'o',
            '1,-1': 'v',
            '-1,1': '^',
            '-1,-1': '<'
        }
        df_breakpoint['color'] = df_breakpoint['strand_config'].map(color_map)
        df_breakpoint['marker'] = df_breakpoint['strand_config'].map(marker_map)

        # plot
        ncols = 4
        nrows = 1
        fig, ax = plt.subplots(nrows, ncols, figsize=(10,4), sharex=True, sharey=True)
        for i in range(ncols):
            strand_config = list(color_map.keys())[i]
            strand_df = df_breakpoint[df_breakpoint['strand_config'] == strand_config]
            if not strand_df.empty:
                ax[i].scatter(
                    strand_df['left_breakpoint'],
                    strand_df['right_breakpoint'],
                    color=strand_df['color'],
                    marker=strand_df['marker'].iat[0],
                    alpha=alpha,
                    s=size,
                )
                ax[i].set_title(f"Strand config: {strand_config}")
            else:
                ax[i].set_title(f"Strand config: {strand_config} (no reads)")
            ax[i].set_xlabel('Left gene breakpoint')
            ax[i].set_ylabel('Right gene breakpoint')
        plt.tight_layout()
        if title:
            plt.suptitle(title)
            plt.subplots_adjust(top=0.88)
        plt.savefig(outfile)
    else:
        plt.figure(figsize=(6,6))
        plt.scatter(df_breakpoint['left_breakpoint'], df_breakpoint['right_breakpoint'], alpha=alpha, s=size)
        plt.xlabel('Left gene breakpoint')
        plt.ylabel('Right gene breakpoint')
        if title:
            plt.title(title)
        plt.savefig(outfile)
        
            
    
    return None