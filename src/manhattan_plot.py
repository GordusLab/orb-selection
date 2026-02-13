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


def manhattan_plot(data_df, chrom_col, pos_col, pvals_col, alpha=0.05,
                   other_var=None, title=None, drop_dups=True, logscale=False,
                   x=None, y=None, dot_label=None, hline_val=None, colors='bright', ax=None, fig=None,
                   overlay_params=None, label_col_name=None, counter=0):
    """
    Create a Manhattan plot for genomic data visualization.
    
    Args:
        data_df: DataFrame containing chromosome, position, and p-values
        chrom_col: Column name for chromosome
        pos_col: Column name for genomic start position
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

    merged_df = data_df.rename(columns={
        chrom_col: 'Chromosome',
        pos_col: 'Position',
        pvals_col: 'P'
    }).copy()

    if drop_dups:
        merged_df = merged_df.sort_values(by=['P'], ascending=True)
        merged_df = merged_df.drop_duplicates(subset=['Chromosome', 'Position'], keep='first').dropna()
    else:
        merged_df = merged_df.dropna()

    if other_var is not None:
        merged_df = merged_df.astype({'Chromosome': str, 'Position': int, 'P': float, other_var: float})
    else:
        merged_df = merged_df.astype({'Chromosome': str, 'Position': int, 'P': float})

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
    
    if overlay_params is not None:
        colors='pastel'

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
        manhattan_plot(overlay_params['data_df'],
                   chrom_col=overlay_params['chrom_col'],
                   pos_col=overlay_params['pos_col'],
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
    python manhattan_plot.py data.csv Chromosome Position p_value --title "GWAS Results"
    python manhattan_plot.py data.csv chrom start pval --other-var-col logFC --alpha 0.01 --output plot.png
        """
    )
    
    # Required arguments
    parser.add_argument('data_file', help='CSV file containing chromosome, position, and p-values')
    parser.add_argument('chrom_col', help='Column name for chromosome')
    parser.add_argument('pos_col', help='Column name for genomic start position')
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
    parser.add_argument('--other-var-col', type=str, default=None,
                       help='Optional column name for alternative y-axis variable')
    parser.add_argument('--label-col', type=str, default=None,
                       help='Optional column name for point labels')
    parser.add_argument('--figsize', nargs=2, type=float, default=[12, 6],
                       help='Figure size as width height (default: 12 6)')
    parser.add_argument('--dpi', type=int, default=300,
                       help='Output resolution in DPI (default: 300)')
    
    args = parser.parse_args()
    
    try:
        # Load data files
        print(f"Loading data from: {args.data_file}")
        data_df = pd.read_csv(args.data_file)

        # Validate required columns
        for col in [args.chrom_col, args.pos_col, args.pvals_col]:
            if col not in data_df.columns:
                raise ValueError(f"Column '{col}' not found in data file")
        if args.other_var_col is not None and args.other_var_col not in data_df.columns:
            raise ValueError(f"Column '{args.other_var_col}' not found in data file")
        if args.label_col is not None and args.label_col not in data_df.columns:
            raise ValueError(f"Column '{args.label_col}' not found in data file")
        
        # Create the plot
        print("Creating Manhattan plot...")
        fig, ax, merged_df = manhattan_plot(
            data_df=data_df,
            chrom_col=args.chrom_col,
            pos_col=args.pos_col,
            pvals_col=args.pvals_col,
            alpha=args.alpha,
            other_var=args.other_var_col,
            title=args.title,
            drop_dups=not args.no_drop_dups,
            logscale=args.logscale,
            colors=args.colors,
            label_col_name=args.label_col
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
