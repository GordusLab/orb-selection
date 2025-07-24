#!/usr/bin/env python3

"""
Manhattan plot visualization for genomic data.

This script creates Manhattan plots for visualizing p-values or other statistics
across chromosomes, with support for overlaying multiple datasets.
"""

import argparse
import itertools
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def manhattan_plot(pvals_df, coords_df, LOCs_col, pvals_col, alpha=0.05,
                   other_var=None, title=None, drop_dups=True, logscale=False,
                   x=None, y=None, dot_label=None, hline_val=None, colors='pastel', ax=None, fig=None,
                   overlay_params=None, label_col_name=None, counter=0):
    """
    Create a Manhattan plot for genomic data visualization.
    
    Args:
        pvals_df: DataFrame containing p-values or statistics
        coords_df: DataFrame containing genomic coordinates
        LOCs_col: Column name for gene/locus identifiers
        pvals_col: Column name for p-values or statistics
        alpha: Significance threshold (default: 0.05)
        other_var: Alternative variable to plot on y-axis
        title: Plot title
        drop_dups: Whether to drop duplicate entries (default: True)
        logscale: Use log scale for y-axis (default: False)
        x: X coordinates for additional points
        y: Y coordinates for additional points
        dot_label: Label for additional points
        hline_val: Y-value for horizontal reference line
        colors: Color palette ('pastel', 'bright', etc.)
        ax: Matplotlib axes object (for subplots)
        fig: Matplotlib figure object (for subplots)
        overlay_params: Parameters for overlaying additional data
        label_col_name: Column name for point labels
        counter: Internal counter for recursive calls
        
    Returns:
        tuple: (figure, axes, merged_dataframe)
    """

    pvals_df = pvals_df.rename(columns={pvals_col: 'P', LOCs_col: 'LOC'})

    if drop_dups:
        pvals_df = pvals_df.sort_values(by=['P'], ascending=True)
        pvals_df = pvals_df.drop_duplicates(subset=['LOC'], keep='first').dropna()

    else:
        pvals_df = pvals_df.dropna()

    if other_var is not None:    
        pvals_df = pvals_df.astype({'LOC': str, 'P': float, other_var: float})

    else:
        pvals_df = pvals_df.astype({'LOC': str, 'P': float})


    coords_df=coords_df.rename(columns={'NCBI GeneID': 'LOC',
                                   'Annotation Genomic Range Start': 'Position',
                                   'Chromosomes': 'Chromosome'})

    coords_df=coords_df.astype({'LOC': str, 'Position': int, 'Chromosome': str})

    merged_df = coords_df.merge(pvals_df, how='right', on='LOC')
    merged_df = merged_df[merged_df['Chromosome'] != 'Un']
    merged_df = merged_df.astype({'Chromosome': int})

    chrom_lens = [272330431,221850441,218885758,195363447,185519777,172884048,172099698,161310338,159483530,139172590]
    chrom_ends = list(itertools.accumulate(chrom_lens))
    chrom_starts = [0] + chrom_ends[:-1]

    values = [0.5*(chrom_ends[i] + chrom_starts[i]) for i in range(10)]

    merged_df['Position'] = merged_df.apply(
        lambda row: row['Position'] + chrom_starts[row['Chromosome']-1], axis=1)

    if other_var is not None:
        y_ax_var = other_var
    else:
        y_ax_var = '-log10(P)'

    merged_df['-log10(P)'] = -np.log10(merged_df['P'])

    if counter == 0:
        fig, ax=plt.subplots(figsize=(12, 6))
        ax.set_title(title)
        ax.set_xlabel('Chromosome')
        ax.set_ylabel(y_ax_var)
        if logscale:
            ax.set_yscale('log')
        ax.set_xticks(values)
        ax.set_xticklabels(['1','2','3','4','5','6','7','8','9','10'])
        ax.scatter(x,y, color='red', s=10, label=(dot_label if 'dot_label' in locals() else None))
        if other_var is None:
            ax.hlines(y=-np.log10(alpha), xmin=np.min(merged_df['Position']), 
                        xmax= np.max(merged_df['Position']),
                        color='gray', linestyle='--', label='p=0.05')
        else:
            if hline_val is not None:
                ax.hlines(y=hline_val, xmin=np.min(merged_df['Position']), 
                            xmax= np.max(merged_df['Position']),
                            color='gray', linestyle='--', label='p=0.05')
        size=10
    else:
        size=15
    
    sns.scatterplot(
        data=merged_df, 
        x='Position', 
        y=y_ax_var, 
        hue='Chromosome', 
        palette=colors, 
        s=size, 
        alpha=0.7,
        ax=ax
    )

    if label_col_name is not None:
        merged_df.apply(lambda x: ax.text(x['Position']+0.01, x[y_ax_var], 
                    x[label_col_name], fontsize=8, ha='right', va='bottom'), axis=1)

        # for i in range(len(merged_df)):
        #     ax.text(merged_df['Position'].iloc[i], merged_df[y_ax_var].iloc[i], 
        #             merged_df[label_col_name].iloc[i], fontsize=8, ha='right', va='bottom')
        

    # if counter == 0:
    #     # legends = []
    #     # legend = ax.legend(loc='upper right', fontsize='small', title=legend_title if legend_title else '')
    #     # legends.append(legend)
    #     ax.legend_.remove()

    counter += 1
    
    if overlay_params is not None:
        # overlay_params['pvals_df']
        manhattan_plot(overlay_params['pvals_df'], 
                       overlay_params['coords_df'],
                       LOCs_col=overlay_params['LOCs_col'], 
                       pvals_col=overlay_params['pvals_col'],
                       alpha=overlay_params.get('alpha', 0.05), 
                       other_var=overlay_params.get('other_var'),
                       title=overlay_params.get('title', ''),
                       drop_dups=overlay_params.get('drop_dups', True), 
                       logscale=overlay_params.get('logscale', False),
                       x=overlay_params.get('x'), 
                       y=overlay_params.get('y'),
                       dot_label=overlay_params.get('dot_label'), 
                       hline_val=overlay_params.get('hline_val'), 
                       colors=overlay_params.get('colors', 'bright'),
                       overlay_params=overlay_params.get('overlay_params', None), 
                       label_col_name=overlay_params.get('label_col_name', None),
                       ax=ax, 
                       fig=fig, 
                       counter=counter)
 
    if ax.get_legend() is not None:
        ax.get_legend().remove()

    fig.tight_layout()
    
    return fig, ax, merged_df


def main():
    """Main function to run Manhattan plot from command line."""
    parser = argparse.ArgumentParser(
        description='Create Manhattan plots for genomic data visualization',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python manhattan_plot.py pvals.csv coords.csv LOC p_value --title "GWAS Results"
  python manhattan_plot.py data.csv coords.csv gene_id pval --alpha 0.01 --output plot.png
        """
    )
    
    # Required arguments
    parser.add_argument('pvals_file', help='CSV file containing p-values or statistics')
    parser.add_argument('coords_file', help='CSV file containing genomic coordinates')
    parser.add_argument('locs_col', help='Column name for gene/locus identifiers')
    parser.add_argument('pvals_col', help='Column name for p-values or statistics')
    
    # Optional arguments
    parser.add_argument('--alpha', type=float, default=0.05,
                       help='Significance threshold (default: 0.05)')
    parser.add_argument('--title', type=str, default='Manhattan Plot',
                       help='Plot title')
    parser.add_argument('--output', '-o', type=str, default='manhattan_plot.png',
                       help='Output file name (default: manhattan_plot.png)')
    parser.add_argument('--colors', type=str, default='pastel',
                       choices=['pastel', 'bright', 'dark', 'colorblind'],
                       help='Color palette (default: pastel)')
    parser.add_argument('--logscale', action='store_true',
                       help='Use log scale for y-axis')
    parser.add_argument('--no-drop-dups', action='store_true',
                       help='Do not drop duplicate entries')
    parser.add_argument('--figsize', nargs=2, type=float, default=[12, 6],
                       help='Figure size as width height (default: 12 6)')
    parser.add_argument('--dpi', type=int, default=300,
                       help='Output resolution in DPI (default: 300)')
    
    args = parser.parse_args()
    
    try:
        # Load data files
        print(f"Loading p-values from: {args.pvals_file}")
        pvals_df = pd.read_csv(args.pvals_file)
        
        print(f"Loading coordinates from: {args.coords_file}")
        coords_df = pd.read_csv(args.coords_file)
        
        # Validate required columns
        if args.locs_col not in pvals_df.columns:
            raise ValueError(f"Column '{args.locs_col}' not found in p-values file")
        if args.pvals_col not in pvals_df.columns:
            raise ValueError(f"Column '{args.pvals_col}' not found in p-values file")
            
        # Check for expected coordinate columns (with flexibility)
        coord_cols = coords_df.columns.tolist()
        print(f"Available coordinate columns: {coord_cols}")
        
        # Create the plot
        print("Creating Manhattan plot...")
        fig, ax, merged_df = manhattan_plot(
            pvals_df=pvals_df,
            coords_df=coords_df,
            LOCs_col=args.locs_col,
            pvals_col=args.pvals_col,
            alpha=args.alpha,
            title=args.title,
            drop_dups=not args.no_drop_dups,
            logscale=args.logscale,
            colors=args.colors
        )
        
        # Set figure size
        fig.set_size_inches(args.figsize[0], args.figsize[1])
        
        # Save the plot
        print(f"Saving plot to: {args.output}")
        fig.savefig(args.output, dpi=args.dpi, bbox_inches='tight')
        
        print(f"Successfully created Manhattan plot!")
        print(f"- Total points plotted: {len(merged_df)}")
        print(f"- Significance threshold: {args.alpha}")
        print(f"- Output saved as: {args.output}")
        
        # Show basic statistics
        if args.pvals_col in merged_df.columns:
            sig_count = len(merged_df[merged_df[args.pvals_col] <= args.alpha])
            print(f"- Significant points (p <= {args.alpha}): {sig_count}")
        
    except FileNotFoundError as e:
        print(f"Error: File not found - {e}")
        return 1
    except ValueError as e:
        print(f"Error: {e}")
        return 1
    except Exception as e:
        print(f"Unexpected error: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())
