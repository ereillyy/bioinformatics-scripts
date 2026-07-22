import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def normalise_data(df, fluorophores, by_group_id=True, by_date_only=False, condition='S9only',
                    method_nucleus='Mean', output_csv=None):
    """
    Normalise intensity data per group_id or per date using S9 only control.

    output_csv : str or None
        Full path to save the normalised dataframe to, or None to skip saving.
    """
    # Check that exactly one normalisation method is selected
    if sum([by_group_id, by_date_only]) != 1:
        raise ValueError("Please select exactly one normalisation method: by_group_id or by_date_only.")
    df_normalised = df.copy()
    # Choose grouping column
    group_col = 'group_id' if by_group_id else 'Date'
    for name in fluorophores.values():
        intensity_col = f'Intensity Nucleus {name} {method_nucleus}'
        norm_col = f'Normalised Intensity Nucleus {name} {method_nucleus}'
        df_normalised[norm_col] = np.nan
        for group_val, group in df_normalised.groupby(group_col):
            s9_mask = group['label'] == condition
            s9_median = group.loc[s9_mask, intensity_col].median()
            if pd.isna(s9_median) or s9_median == 0:
                print(f"Warning: {condition} median is NaN or zero for {group_col} {group_val}, cannot normalise {name}.")
                continue
            df_normalised.loc[group.index, norm_col] = group[intensity_col] / s9_median
            print(f"{group_col} {group_val}: Normalised {name} using {condition} median {s9_median:.2f}")
    # Save normalised data
    if output_csv is not None:
        dirname = os.path.dirname(output_csv)
        if dirname:
            os.makedirs(dirname, exist_ok=True)
        df_normalised.to_csv(output_csv, index=False)
        print(f"\nNormalised data saved to {output_csv}")
    return df_normalised


def filter_data(df, dates='all', cells='all', samples='all', recovery_time='all'):
    """
    Filter data based on specified dates, cell lines, and sample treatments.
    """
    df_filtered = df.copy()
    if dates != 'all':
        df_filtered = df_filtered[df_filtered['Date'].astype(str).isin(dates)]
        print(f"Filtered by dates {dates}: {df_filtered.shape[0]} nuclei remain.")
    if cells != 'all':
        df_filtered = df_filtered[df_filtered['Cells'].isin(cells)]
        print(f"Filtered by cells {cells}: {df_filtered.shape[0]} nuclei remain.")
    if samples != 'all':
        df_filtered = df_filtered[df_filtered['label'].isin(samples)]
        print(f"Filtered by samples {samples}: {df_filtered.shape[0]} nuclei remain.")
    if recovery_time != 'all':
        df_filtered = df_filtered[df_filtered['Recovery_time'].isin(recovery_time)]
        print(f"Filtered by recovery time {recovery_time}: {df_filtered.shape[0]} nuclei remain.")
    return df_filtered


def plot_box(df, plot_factor, palettes, hierarchy=['Cells', 'label', 'Date'], colour_by='label',
             downsample=True, plate_thresholds=None, y_col=None, ylim=None, method_nucleus='Mean',
             output_path=None):
    """
    Plot box plot with scatter points showing hierarchical grouping.
    IF DATE IS NOT IN HIERARCHY, THEN ALL DATES ARE COMBINED FOR EACH CELL LINE + SAMPLE COMBO.

    output_path : str or None
        Full path to save the plot to, or None to skip saving.
    """

    # set up col to plot
    y_col = f"Normalised Intensity Nucleus {plot_factor} {method_nucleus}" if y_col is None else y_col

    sns.set(style="whitegrid")
    # get all x axis order
    combos = df[hierarchy].drop_duplicates()
    if 'Cells' in hierarchy:
        combos['Cells_order'] = combos['Cells'].map({k: i for i, k in enumerate(list(palettes['cellline'].keys()))})
    if 'label' in hierarchy:
        combos['label_order'] = combos['label'].map({k: i for i, k in enumerate(list(palettes['sample'].keys()))})
    if 'Recovery_time' in hierarchy:
        combos['Recovery_time_order'] = combos['Recovery_time'].map({k: i for i, k in enumerate(list(palettes['recovery_time'].keys()))})
    if 'Date' in hierarchy:
        combos['Date_order'] = combos['Date'].astype(str).map({k: i for i, k in enumerate(sorted(combos['Date'].astype(str).unique()))})
    else:
        print(f"Date not in hierarchy: combining all dates for each Cells + label combo\nFound dates: {df['Date'].unique()}")
    combos = combos.sort_values([f"{h}_order" for h in hierarchy])
    x_order = combos[hierarchy].astype(str).agg('_'.join, axis=1).tolist()

    # create group labels
    combos['group_label'] = combos[hierarchy].astype(str).agg('_'.join, axis=1)
    df['group_label'] = df[hierarchy].astype(str).agg('_'.join, axis=1)

    # create colour palette for x_order which adapts to hierarchy
    if colour_by == 'Cells':
        palette = combos['Cells'].map(palettes['cellline']).to_list()
    elif colour_by == 'label':
        palette = combos['label'].map(palettes['sample']).to_list()
    elif colour_by == 'Date':
        # generate a color palette with seaborn
        unique_dates = combos['Date'].astype(str).unique()
        date_palette = sns.color_palette("hsv", len(unique_dates))
        date_color_map = {date: date_palette[i] for i, date in enumerate(unique_dates)}
        palette = combos['Date'].astype(str).map(date_color_map).to_list()
    else:
        raise ValueError(f"colour_by must be one of 'Cells', 'label', or 'Date', got '{colour_by}'")

    # Prepare data for jittered scatter plot
    if downsample:
        group_counts = df.groupby('group_label')[y_col].count()
        min_n = group_counts.min()
        min_group = group_counts.idxmin()
        print(f"Downsampling to {min_n} nuclei per group (min: '{min_group}')")
        df_jitter = (
            df
            .groupby('group_label', group_keys=False)
            .apply(lambda x: x.sample(n=min(len(x), min_n), random_state=0), include_groups=False)
            .reset_index(drop=True)
        )
        df_jitter['group_label'] = df_jitter[hierarchy].astype(str).agg('_'.join, axis=1)
    else:
        print("No downsampling: plotting all points")
        df_jitter = df

    # Create box plot
    fig = plt.figure(figsize=(1 + (len(x_order) * 0.25), 5))
    ax = sns.boxplot(
        x='group_label', y=y_col, data=df,
        order=x_order,
        showcaps=True, boxprops={'facecolor': 'none', 'edgecolor': 'black'},
        whiskerprops={'linewidth': 1.5, 'color': 'black'},
        medianprops={'linewidth': 2, 'color': 'black'},
        showfliers=False,
        legend=False
    )

    # Scatter (jitter) plot
    sns.stripplot(
        x='group_label', y=y_col, data=df_jitter,
        order=x_order, palette=palette,
        size=1.6, alpha=0.3, jitter=0.25, ax=ax,
        legend=False
    )

    if ylim:
        plt.ylim(ylim)
    ax.set_xlabel("", fontsize=12)
    ax.set_ylabel(y_col, fontsize=12)
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    if output_path is not None:
        dirname = os.path.dirname(output_path)
        if dirname:
            os.makedirs(dirname, exist_ok=True)
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"\nPlot saved to {output_path}")

    plt.show()

def plot_median_bars(
    medians,
    plot_factor,
    palettes,
    hierarchy=['Cells', 'label', 'Date'],
    method_nucleus='Mean',
    combine_dates=True,
    normalise_to_s9=False,
    output_path=None,
    figsize=(max(6, 0.25 * 20), 5),
    date_palette_name="hsv",
):
    """
    Barplot of median values.

    - medians: path to a medians CSV or a pd.DataFrame (one row per Date/group).
    - plot_factor: readable fluorophore name (e.g. 'EdU').
    - palettes: same palettes dict used by `plot_box` for ordering/colours.
    - hierarchy: must include 'Date'. When `combine_dates=True`, bars are
      computed per hierarchy except Date (i.e. group across dates).
    - normalise_to_s9: if True, divides each row's median by the S9only median
      for the same Date and Cells before computing means/stds.
    - combine_dates: when True, one bar per (Cells,label,... excluding Date).
    - output_path: save figure if not None.
    """
    if isinstance(medians, str):
        dfm = pd.read_csv(medians)
    else:
        dfm = medians.copy()

    if 'Date' not in hierarchy:
        raise ValueError("`hierarchy` must include 'Date' so dates can be handled.")

    value_col = f'Intensity Nucleus {plot_factor} {method_nucleus}'

    if value_col not in dfm.columns:
        raise ValueError(f"{value_col} not found in medians data.")

    # Optionally normalise medians to S9only per Date+Cells
    if normalise_to_s9:
        if 'label' not in dfm.columns:
            raise ValueError("medians DataFrame must include 'label' to normalise to S9only.")
        # build S9 lookup: S9 value per Date + Cells (and any extra grouping except label)
        s9_lookup = (
            dfm[dfm['label'] == 'S9only']
            .set_index(['Date', 'Cells'])[value_col]
            .to_dict()
        )
        def _norm_row(r):
            key = (r['Date'], r['Cells'])
            s9 = s9_lookup.get(key, np.nan)
            if pd.isna(s9) or s9 == 0:
                return np.nan
            return r[value_col] / s9
        dfm[f'Norm_{value_col}'] = dfm.apply(_norm_row, axis=1)
        dfm_used = dfm.rename(columns={f'Norm_{value_col}': value_col})
    else:
        dfm_used = dfm

    # Determine grouping when combining dates
    if combine_dates:
        base_hierarchy = [h for h in hierarchy if h != 'Date']
        missing = [h for h in base_hierarchy if h not in dfm_used.columns]
        if missing:
            raise ValueError(f"Missing grouping columns for combine: {missing}")
        dfm_used['group_base'] = dfm_used[base_hierarchy].astype(str).agg('_'.join, axis=1)
        # generate ordering consistent with plot_box
        combos = dfm_used[base_hierarchy].drop_duplicates().reset_index(drop=True)
        if 'Cells' in base_hierarchy:
            combos['Cells_order'] = combos['Cells'].map({k: i for i, k in enumerate(list(palettes['cellline'].keys()))})
        if 'label' in base_hierarchy:
            combos['label_order'] = combos['label'].map({k: i for i, k in enumerate(list(palettes['sample'].keys()))})
        if 'Date' in base_hierarchy:
            combos['Date_order'] = combos['Date'].astype(str).map({k: i for i, k in enumerate(sorted(combos['Date'].astype(str).unique()))})
        order_cols = [f"{h}_order" for h in base_hierarchy if f"{h}_order" in combos.columns]
        if order_cols:
            combos = combos.sort_values(order_cols)
        x_order = combos[base_hierarchy].astype(str).agg('_'.join, axis=1).tolist()
        combos['group_base'] = combos[base_hierarchy].astype(str).agg('_'.join, axis=1)
        # palette per bar (adapts to colour_by behaviour used in plot_box)
        # default: colour by 'label' ordering
        if 'label' in base_hierarchy:
            # map each group_base -> label (first label found)
            label_map = dfm_used.groupby('group_base')['label'].first().to_dict()
            # collect used labels, detect NA separately
            labels_used = {l for l in label_map.values() if pd.notna(l)}
            has_na_label = any(pd.isna(l) for l in label_map.values())
            palette_keys = set(palettes.get('sample', {}).keys())
            missing = sorted(labels_used - palette_keys)
            if has_na_label:
                missing.append('<NA>')
            if missing:
                raise ValueError(f"Missing colour entries in palettes['sample'] for labels: {missing}")
            # build palette in the same order as combos
            bar_palette = [palettes['sample'][label_map[gb]] for gb in combos['group_base']]
        else:
            bar_palette = ['grey'] * len(x_order)
        # compute mean and std across dates for each group_base from medians table
        stats = (
            dfm_used
            .groupby('group_base')[value_col]
            .agg(['mean', 'std', 'count'])
            .reindex(x_order)
            .reset_index()
        )
        stats.rename(columns={'mean': 'mean_value', 'std': 'std_value', 'count': 'n_dates'}, inplace=True)

        # Prepare mapping from group_base->x position
        stats['x'] = np.arange(len(stats))

        # create date colour map (rainbow in date order)
        unique_dates = sorted(dfm_used['Date'].astype(str).unique())
        date_palette = sns.color_palette(date_palette_name, len(unique_dates))
        date_color_map = {d: date_palette[i] for i, d in enumerate(unique_dates)}

        # build plot
        fig, ax = plt.subplots(figsize=figsize)
        ax.bar(stats['x'], stats['mean_value'], yerr=stats['std_value'], color=bar_palette, capsize=5, edgecolor='black')
        # overlay per-date points (one dot per date per group_base)
        jitter = 0.12
        for i, gb in enumerate(stats['group_base']):
            sub = dfm_used[dfm_used['group_base'] == gb]
            # for visual separation, offset each date by small step around x
            for j, (_, row) in enumerate(sub.sort_values('Date').iterrows()):
                dx = (j - 0.5 * (len(sub) - 1)) * (jitter / max(len(sub), 1))
                ax.scatter(i + dx, row[value_col], color=date_color_map[str(row['Date'])], zorder=5, s=30, edgecolor='white')

        ax.set_xticks(stats['x'])
        ax.set_xticklabels(stats['group_base'], rotation=45, ha='right')
        ax.set_ylabel(value_col)
        ax.set_title(f"{plot_factor} medians (bars = mean over dates; dots = per-date medians)")
        # date legend
        handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=date_color_map[d], markersize=6) for d in unique_dates]
        ax.legend(handles, unique_dates, title='Date', bbox_to_anchor=(1.02, 1), loc='upper left', frameon=False)
        plt.tight_layout()

        if output_path is not None:
            dirname = os.path.dirname(output_path)
            if dirname:
                os.makedirs(dirname, exist_ok=True)
            fig.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"Saved median barchart to {output_path}")

        return stats

    else:
        raise NotImplementedError("combine_dates=False is not implemented; pass medians grouped how you want.")