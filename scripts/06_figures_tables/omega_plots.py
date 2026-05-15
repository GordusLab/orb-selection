import numpy as np
from matplotlib import (
    pyplot as plt,
    patheffects as pe,
    ticker
)

logbins = np.geomspace(0.001, 10000, 100)

def plot_omega_distributions(
        x, 
        result, 
        top_title, 
        bottom_title,
        numeral="",
        suptitle=None,
        filename=None,
        shift_top_title=False,
        transparent=True
        ):

    """
    Plot the distributions of ω1, ω2, and ω3 for both TEST and REFERENCE groups.
    Parameters:
    - x: DataFrame containing the ω values and their proportions.
    - result: String indicating the type of result (e.g., 'all', 'intensified', 'relaxed').
    - top_title: Title for the top plot.
    - bottom_title: Title for the bottom plot.
    - filename: Path to save the figure (optional).
    - transparent: Boolean indicating whether to save with transparent background.
    Returns:
    - fig: The figure object.
    - axs: The axes object.
    """
    # Initialize lists to store the maximum y-axis limits for each group
    # and the second highest maximums
    # for the two groups.


    ymaxs = []

    ymaxs_1 = []
    ymins_1 = []

    ymaxs_2 = []

    # Calculate the maximum counts for each group and set the y-axis limits
    for i, group in enumerate(['test','ref']):
        ω3_counts = np.histogram(x[f'ω3_{group}'], bins=logbins)[0]
        ω2_counts = np.histogram(x[f'ω2_{group}'], bins=logbins)[0]

        ω3_max = np.max(ω3_counts)

        # second heighest maximum for ω3
        ω3_max2 = np.partition(ω3_counts, -2)[-2] 
        ω2_max = np.max(ω2_counts)

        maxes = np.array([ω3_max, ω3_max2, ω2_max])

        ymax = np.max(maxes)

        # second heighest maximum
        ymax2 = np.partition(maxes, -2)[-2] 

        pad = ymax/10

        ymaxs_1.append(ymax+0.75*pad)
        ymins_1.append(ymax-pad)
        
        ymaxs_2.append(ymax2+pad)
        ymaxs.append(ymax)

    if result != 'busted':
    # if any(y > 50 for y in ymaxs):
    # if ymaxs[0]>50 | ymaxs[1]>50:
        # Calculate the height ratios for the subplots
        a = ymaxs_1[0]-ymins_1[0]
        b = ymaxs_2[0]
        c = ymaxs_2[1]
        d = ymaxs_1[0]-ymins_1[1]

        hra = 4/(1+b/a)
        hrb = 4-hra
        hrd = 4/(1+c/d)
        hrc = 4-hrd

        fig, axs = plt.subplots(6,1, sharex=True, height_ratios=[hra,0.2,hrb,hrc,0.2,hrd], figsize=(6,5))

        plt.subplots_adjust(hspace=0)
        plt.xlim(0.001, 100000)
        plt.xscale('log')
        plt.rcParams['font.family'] = 'Verdana'

        axs[1].set_axis_off()
        axs[4].set_axis_off()

        axs[0].set_ylim(ymins_1[0], ymaxs_1[0])
        axs[2].set_ylim(0, ymaxs_2[0])
        axs[3].set_ylim(0, ymaxs_2[1])
        axs[3].invert_yaxis()
        axs[5].set_ylim(ymins_1[1],ymaxs_1[0])
        axs[5].invert_yaxis()   
        
        # Plot distributions for ω1, ω2, and ω3 for the TEST group
        intervals= []
        for ax in [axs[0], axs[2]]:
            ax.hist(x['ω1_test'], bins=logbins, histtype='stepfilled', 
                    color='salmon', alpha=0.17, label='ω1 distribution',
                    path_effects=[pe.Stroke(linewidth=1, foreground='darkred')])
            ax.hist(x['ω2_test'], bins=logbins, histtype='stepfilled', 
                    color='steelblue', alpha =0.17, label='ω2 distribution',
                    path_effects=[pe.Stroke(linewidth=1, foreground='darkblue')])
            ax.hist(x['ω3_test'], bins=logbins, histtype='stepfilled', 
                    color='goldenrod', alpha =0.17, label='ω3 distribution',
                    path_effects=[pe.Stroke(linewidth=1, foreground='brown')])

            locator = ticker.AutoLocator()
            ax.yaxis.set_major_locator(locator)
            intervals.append(locator()[1]-locator()[0])

        # axs[0].legend(fontsize='small')

        # Plot distributions for ω1, ω2, and ω3 for the REFERENCE group
        for ax in [axs[3], axs[5]]:
            ax.hist(x['ω1_ref'], bins=logbins, histtype='stepfilled', color='salmon', 
                    alpha=0.17, path_effects=[pe.Stroke(linewidth=1, foreground='darkred')])
            ax.hist(x['ω2_ref'], bins=logbins, histtype='stepfilled', color='steelblue', 
                    alpha =0.17, path_effects=[pe.Stroke(linewidth=1, foreground='darkblue')])
            ax.hist(x['ω3_ref'], bins=logbins, histtype='stepfilled', color='goldenrod', 
                    alpha =0.17, path_effects=[pe.Stroke(linewidth=1, foreground='brown')])

            locator = ticker.AutoLocator()
            ax.yaxis.set_major_locator(locator)
            intervals.append(locator()[1]-locator()[0])

        if any(i>=100 for i in intervals):
            [axs[j].yaxis.set_major_locator(ticker.MultipleLocator(100)) for j in [0,2,3,5]]
        else:
            [axs[j].yaxis.set_major_locator(ticker.MultipleLocator(np.max(intervals))) for j in [0,2,3,5]]

        axs[0].spines.bottom.set_visible(False)
        axs[2].spines.top.set_visible(False)
        axs[3].spines.bottom.set_visible(False)
        axs[5].spines.top.set_visible(False)

        # Add subplots to plot the means of ω1, ω2, and ω3 for both TEST and REFERENCE groups
        gs = axs[0].get_gridspec()
        ax_avgs = fig.add_subplot(gs[:3], sharex=axs[5])
        ax_avgs2 = fig.add_subplot(gs[3:], sharex=axs[5])


        ax_avgs.set_facecolor('none')
        ax_avgs.set_axis_off()

        ax_avgs = ax_avgs.twinx()
        means = [x['ω1_test'].mean(), x['ω2_test'].mean(), x['ω3_test'].mean()]
        weights = [x['ω1_test_P'].mean(), x['ω2_test_P'].mean(), x['ω3_test_P'].mean()]

        # Plot vertical lines for the means of ω1, ω2, and ω3 for the TEST group
        ax_avgs.axvline(means[0], linewidth=6, color='salmon', ymax=weights[0], 
                        path_effects=[pe.withStroke(linewidth=7, foreground='white'), pe.Normal()])
        ax_avgs.axvline(means[1], linewidth=6, color='steelblue', ymax=weights[1], 
                        path_effects=[pe.withStroke(linewidth=7, foreground='white'), pe.Normal()])
        ax_avgs.axvline(means[2], linewidth=6, color='goldenrod', ymax=weights[2], 
                        path_effects=[pe.withStroke(linewidth=7, foreground='white'), pe.Normal()])
        ax_avgs.spines.right.set_visible(False)

        ax_avgs2.set_facecolor('none')
        ax_avgs2.set_axis_off()

        ax_avgs2 = ax_avgs2.twinx()
        ax_avgs2.invert_yaxis()
        means = [x['ω1_ref'].mean(), x['ω2_ref'].mean(), x['ω3_ref'].mean()]
        weights = [x['ω1_ref_P'].mean(), x['ω2_ref_P'].mean(), x['ω3_ref_P'].mean()]

        # Plot vertical lines for the means of ω1, ω2, and ω3 for the REFERENCE group
        ax_avgs2.axvline(means[0],  linewidth=6, color='salmon', 
                        ymin=(1-weights[0]), label='mean inferred ω1',
                        path_effects=[pe.withStroke(linewidth=7.5, foreground='white'), pe.Normal()])
        ax_avgs2.axvline(means[1], linewidth=6, color='steelblue', 
                        ymin=(1-weights[1]), label='mean inferred ω2',
                        path_effects=[pe.withStroke(linewidth=7.5, foreground='white'), pe.Normal()])
        ax_avgs2.axvline(means[2], linewidth=6, color='goldenrod', 
                        ymin=(1-weights[2]), label='mean inferred ω3',
                        path_effects=[pe.withStroke(linewidth=7.5, foreground='white'), pe.Normal()])
        ax_avgs2.spines.right.set_visible(False)

        ax_avgs.axvline(1, linewidth=0.5, linestyle='dashed', color='k', alpha=0.5)
        ax_avgs2.axvline(1, linewidth=0.5, linestyle='dashed', color='k', alpha=0.5)

        if result == 'all':
            
            # labels for top vs bottom axes
            ax_avgs.set_title(top_title, x=0.2, y=0.7, fontsize=13, color='white', weight='bold', backgroundcolor='lightgray')
            ax_avgs2.set_title(bottom_title, x=0.2, y=0.1, fontsize=13, color='white', weight='bold', backgroundcolor='lightgray')

        else: 
            # labels for top vs bottom axes
            bottom_title_x = 0.185 if result in ('relaxed', 'intensified') else 0.2
            ax_avgs.set_title(top_title, x=0.75, y=0.7, fontsize=13, color='white', weight='bold', backgroundcolor='lightgray')
            ax_avgs2.set_title(bottom_title, x=bottom_title_x, y=0.075, fontsize=13, color='silver', weight='bold', backgroundcolor='white')

            if result == 'intensified':
                plt.text(0.8, 0.7, 
                    '$\it{p}$≤0.05\n$\it{k}$>1', 
                    fontsize=15, ha='left', va='center', transform=ax.transAxes, 
                    bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
                
            else:
                plt.text(0.8, 0.35, 
                    '$\it{p}$≤0.05\n$\it{k}$<1', 
                    fontsize=15, ha='left', va='center', transform=ax.transAxes, 
                    bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

        # ax_avgs2.legend(loc='lower right', fontsize='small')
        # for ax in [axs[3],axs[4],axs[5], ax_avgs2]:
        #     ax.set_facecolor('lightgray')

        for ax in [ax_avgs, ax_avgs2]:
            ax.yaxis.set_label_position('left')
            ax.yaxis.tick_left()
            for tick in ax.get_yticklabels():
                tick.set_fontweight('bold')
                tick.set_fontsize(14)
                # tick.set_alpha(0.8)

        for i in [0,2,3,5]:
            axs[i].yaxis.set_label_position('right')
            axs[i].yaxis.tick_right()
            for tick in axs[i].get_yticklabels():
                tick.set_fontsize(14)
                # tick.set_alpha(0.7)


        fig.supylabel('mean proportion of sites', x=0.01, weight='bold', fontsize=16)

        dummy_ax=fig.add_subplot(1,1,1)
        dummy_ax.set_xticks([])
        dummy_ax.set_yticks([])
        [dummy_ax.spines[side].set_visible(False) for side in ('left', 'top', 'right', 'bottom')]
        dummy_ax.patch.set_visible(False)
        dummy_ax.yaxis.set_label_position('right')
        # dummy_ax.yaxis.label.set(rotation=270)
        dummy_ax.set_ylabel('number of orthogroups', 
                            labelpad=48, rotation=270, 
                            color='whitesmoke', alpha=0.6, fontsize=16,
                            path_effects=[pe.withStroke(linewidth=1, foreground="black")])

        axs[0].tick_params(bottom=False, which='both')
        axs[3].tick_params(bottom=False, which='both')
        # axs[5].tick_params(bottom=False, which='minor')


        # Add markers to split axes 
        d = 0.1

        kwargs = dict(marker=[(-1, -d), (1, d)], markersize=8,
                linestyle="none", color='k', mec='dimgrey', mew=1, clip_on=False)
        axs[0].plot([None,0], transform=axs[0].transAxes, **kwargs)
        axs[2].plot((None,1), transform=axs[2].transAxes, **kwargs)
        axs[3].plot((None,0), transform=axs[3].transAxes, **kwargs)
        axs[5].plot((None,1), transform=axs[5].transAxes, **kwargs)

        for ax in fig.axes:
            ax.tick_params(direction='in')
            ax.tick_params(axis='x', which='major', pad=10)

        for tick in axs[5].get_xticklabels():
            tick.set_fontsize(14)

    else:
        # If the maximum counts are less than or equal to 50, use a simpler layout
        # This is what we are using for the BUSTED results

        fig, axs = plt.subplots(2,1, sharex=True, figsize=(6,5))

        plt.subplots_adjust(hspace=0)
        plt.xlim(0.001, 100000)
        plt.xscale('log')
        plt.rcParams['font.family'] = 'Verdana'

        # make the y axes the same height
        axs[0].set_ylim(0, np.max(ymaxs_1))
        axs[1].set_ylim(0, np.max(ymaxs_1))
        
        axs[1].invert_yaxis()
        # axs[0].set_facecolor('lightgray')

        # axs[0].legend(fontsize='small')

        # Plot distributions for ω1, ω2, and ω3 for the TEST group
        axs[0].hist(x['ω1_test'], bins=logbins, histtype='stepfilled', color='salmon', 
                    alpha=0.17, path_effects=[pe.Stroke(linewidth=1, foreground='darkred')])
        axs[0].hist(x['ω2_test'], bins=logbins, histtype='stepfilled', color='steelblue', 
                    alpha=0.17, path_effects=[pe.Stroke(linewidth=1, foreground='darkblue')])
        axs[0].hist(x['ω3_test'], bins=logbins, histtype='stepfilled', color='goldenrod', 
                    alpha=0.17, path_effects=[pe.Stroke(linewidth=1, foreground='brown')])                                 
        
        # Plot distributions for ω1, ω2, and ω3 for the REF group
        axs[1].hist(x['ω1_ref'], bins=logbins, histtype='stepfilled', 
                    color='salmon', alpha=0.17, label='ω1 distribution',
                    path_effects=[pe.Stroke(linewidth=1, foreground='darkred')])
        axs[1].hist(x['ω2_ref'], bins=logbins, histtype='stepfilled', 
                    color='steelblue', alpha =0.17, label='ω2 distribution',
                    path_effects=[pe.Stroke(linewidth=1, foreground='darkblue')])
        axs[1].hist(x['ω3_ref'], bins=logbins, histtype='stepfilled', 
                    color='goldenrod', alpha =0.17, label='ω3 distribution',
                    path_effects=[pe.Stroke(linewidth=1, foreground='brown')])


        ax_avgs = axs[0].twinx()
        ax_avgs.set_facecolor('none')
        # ax_avgs.set_axis_off()

        means = [x['ω1_test'].mean(), x['ω2_test'].mean(), x['ω3_test'].mean()]
        weights = [x['ω1_test_P'].mean(), x['ω2_test_P'].mean(), x['ω3_test_P'].mean()]

        # Plot vertical lines for the means of ω1, ω2, and ω3 for the TEST group
        ax_avgs.axvline(means[0], linewidth=6, color='salmon', ymax=weights[0], 
                        path_effects=[pe.withStroke(linewidth=7.5, foreground='white'), pe.Normal()])
        ax_avgs.axvline(means[1], linewidth=6, color='steelblue', ymax=weights[1], 
                        path_effects=[pe.withStroke(linewidth=7.5, foreground='white'), pe.Normal()])
        ax_avgs.axvline(means[2], linewidth=6, color='goldenrod', ymax=weights[2], 
                        path_effects=[pe.withStroke(linewidth=7.5, foreground='white'), pe.Normal()])
        ax_avgs.spines.right.set_visible(False)

        ax_avgs2 = axs[1].twinx()

        # ax_avgs2.set_axis_off()

        ax_avgs2.invert_yaxis()
        
        means = [x['ω1_ref'].mean(), x['ω2_ref'].mean(), x['ω3_ref'].mean()]
        weights = [x['ω1_ref_P'].mean(), x['ω2_ref_P'].mean(), x['ω3_ref_P'].mean()]

        # Plot vertical lines for the means of ω1, ω2, and ω3 for the REF group
        ax_avgs2.axvline(means[0],  linewidth=6, color='salmon', 
                        ymin=(1-weights[0]), label='mean inferred ω1',
                        path_effects=[pe.withStroke(linewidth=7, foreground='white'), pe.Normal()])
        ax_avgs2.axvline(means[1], linewidth=6, color='steelblue', 
                        ymin=(1-weights[1]), label='mean inferred ω2',
                        path_effects=[pe.withStroke(linewidth=7, foreground='white'), pe.Normal()])
        ax_avgs2.axvline(means[2], linewidth=6, color='goldenrod', 
                        ymin=(1-weights[2]), label='mean inferred ω3',
                        path_effects=[pe.withStroke(linewidth=7, foreground='white'), pe.Normal()])
        ax_avgs2.spines.right.set_visible(False)

        ax_avgs.axvline(1, linewidth=0.5, linestyle='dashed', color='k', alpha=0.5)
        ax_avgs2.axvline(1, linewidth=0.5, linestyle='dashed', color='k', alpha=0.5)

        # labels for top vs bottom axes
        if shift_top_title:
            ax_avgs.set_title(top_title, x=0.8, y=0.6, fontsize=13, color='white', weight='bold', backgroundcolor='lightgray')
        else:
            ax_avgs.set_title(top_title, x=0.8, y=0.7, fontsize=13, color='white', weight='bold', backgroundcolor='lightgray')
        ax_avgs2.set_title(bottom_title, x=0.185, y=0.1, fontsize=13, color='silver', weight='bold', backgroundcolor= 'white')

        # ax_avgs2.legend(loc='lower right', fontsize='small')

        for ax in [ax_avgs, ax_avgs2]:
            ax.set_facecolor('none')
            ax.yaxis.set_label_position('left')
            ax.yaxis.tick_left()
            for tick in ax.get_yticklabels():
                tick.set_fontweight('bold')
                tick.set_fontsize(14)
                # tick.set_alpha(0.8)

        for i in [0,1]:
            axs[i].yaxis.set_label_position('right')
            axs[i].yaxis.tick_right()
            for tick in axs[i].get_yticklabels():
                tick.set_fontsize(14)
                # tick.set_alpha(0.7)

        for tick in axs[1].get_xticklabels():
            tick.set_fontsize(14)
            

        # label the two different y axes
        fig.supylabel('mean proportion of sites', x=0.01, weight='bold', fontsize=16)

        dummy_ax=fig.add_subplot(1,1,1)
        dummy_ax.set_xticks([])
        dummy_ax.set_yticks([])
        [dummy_ax.spines[side].set_visible(False) for side in ('left', 'top', 'right', 'bottom')]
        dummy_ax.patch.set_visible(False)
        dummy_ax.yaxis.set_label_position('right')
        # dummy_ax.yaxis.label.set(rotation=270)
        dummy_ax.set_ylabel('number of orthogroups', 
                            labelpad=48, rotation=270, 
                            color='whitesmoke', alpha=0.6, fontsize=16,
                            path_effects=[pe.withStroke(linewidth=1, foreground="black")])

        for ax in fig.axes:
            ax.tick_params(direction='in')
            ax.tick_params(axis='x', which='major', pad=10)

        plt.text(0.775, 0.125, 
                 '$\it{p_1}$≤0.05\n$\it{p_2}$>0.05\n$\it{p_3}$≤0.05', 
                 fontsize=15, ha='left', va='center', transform=ax.transAxes, 
                 bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

    if suptitle is not None:
        y_pos = 0.96 if result == 'all' else 0.98
        fig.suptitle(f"{numeral}{suptitle}", y=y_pos, fontsize=16)

    if filename is not None:
        plt.savefig(filename, dpi=600, transparent=transparent, bbox_inches='tight')

    return fig, axs


def plot_omega_single_gene(
        df, 
        gene, 
        suptitle=None,
        subtitle=None,
        i="", 
        offset_zero=False, 
        k=False,
        filename=None,
        transparent=True,
        build_in=False
        ):

    """
    Plot the ω values for a single gene, including vertical lines for ω1, ω2, and ω3 for both TEST and REFERENCE groups.
    Parameters:
    - df: DataFrame containing the ω values and their proportions for the gene.
    - gene: The gene to plot.
    - plot_title: Title for the plot.
    - i: Index or identifier for the gene (optional).
    - offset_zero: Boolean indicating whether to offset the zero line.
    - k: Boolean indicating whether to include the k value in the plot.
    - filename: Path to save the figure (optional).
    - transparent: Boolean indicating whether to save with transparent background.
    Returns:
    - fig: The figure object.
    - ax: The axes object.
    """

    x = df.loc[gene]

    fig, ax = plt.subplots(figsize=(5.5,5))

    thresh = np.min([x['ω2_ref'], x['ω2_test']]) * 0.1

    plt.subplots_adjust(hspace=0)
    plt.xscale('symlog', linthresh=thresh)
    plt.rcParams['font.family'] = 'Verdana'

    # Plot vertical lines for ω1, ω2, and ω3 for the REFERENCE group
    if offset_zero:
        ref_path_effects = [pe.SimpleLineShadow(offset=(-2,2), alpha=0.3, foreground='salmon'), pe.Normal()]
    else:
        ref_path_effects = [pe.withStroke(linewidth=22, foreground='white'), pe.Normal()]
    ax.axvline(x['ω1_ref'], linewidth=20, color='salmon', ymax=x['ω1_ref_P'], alpha=0.17,
               path_effects=ref_path_effects)
    ax.axvline(x['ω2_ref'], linewidth=20, color='steelblue', ymax=x['ω2_ref_P'], alpha=0.17, 
                    path_effects=[pe.withStroke(linewidth=22, foreground='white'), pe.Normal()])
    ax.axvline(x['ω3_ref'], linewidth=20, color='goldenrod', ymax=x['ω3_ref_P'], alpha=0.17, 
                    path_effects=[pe.withStroke(linewidth=22, foreground='white'), pe.Normal()])
        
    # Plot vertical lines for ω1, ω2, and ω3 for the TEST group (invisible if build_in=True)
    test_alpha = 0 if build_in else 1
    test_effects = [] if build_in else [pe.withStroke(linewidth=22, foreground='white'), pe.Normal()]
    
    ax.axvline(x['ω1_test'], linewidth=20, color='salmon', ymax=x['ω1_test_P'], alpha=test_alpha,
                    path_effects=test_effects)
    ax.axvline(x['ω2_test'], linewidth=20, color='steelblue', ymax=x['ω2_test_P'], alpha=test_alpha,
                    path_effects=test_effects)
    ax.axvline(x['ω3_test'], linewidth=20, color='goldenrod', ymax=x['ω3_test_P'], alpha=test_alpha,
                    path_effects=test_effects)
    
    # mark omega = 1
    ax.axvline(1, linewidth=1, linestyle='dashed', color='k', alpha=0.5)

    for tick in ax.get_yticklabels():
        tick.set_fontweight('bold')
        tick.set_fontsize(14)

    for tick in ax.get_xticklabels():
        tick.set_fontsize(14)

    fig.supylabel('proportion of sites', x=0.01, weight='bold', fontsize=16)
   
    for ax in fig.axes:
        ax.tick_params(direction='in')

    if suptitle is not None:
        if subtitle is not None:
            fig.suptitle(f"{i}{suptitle}", y=0.96, fontsize=14, fontweight='bold')
        else:
            fig.suptitle(f"{i}{suptitle}", y=0.92, fontsize=14, fontweight='bold')
    
    if build_in==False:
        if subtitle is not None:
            plt.title(subtitle, fontsize=12, fontweight='bold')

        if k:
            plt.text(0.7, 0.9,
                f'$\\it{{k}}$={round(x["k"], 2)}\n$\\it{{p}}$={x["p_value"]:.2e}',
                fontsize=13, ha='left', va='center', transform=ax.transAxes,
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
        else:
            plt.text(0.65, 0.9,
                f'$\\it{{p_1}}$={x["test_pval"]:.2e}\n$\\it{{p_3}}$={x["shared_pval"]:.2e}',
                fontsize=13, ha='left', va='center', transform=ax.transAxes,
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

    if filename is not None:
        plt.savefig(filename, dpi=300, transparent=transparent)

    return fig, ax