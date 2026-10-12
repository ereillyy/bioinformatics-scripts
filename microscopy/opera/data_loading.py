import pandas as pd
import matplotlib.pyplot as plt


def load_multiple_plates(plate_list, fluorophores, output_csv, plate_var,
                          method_nucleus='Mean', header_rows=8):
    """
    Load and merge multiple plates with platemap data and thresholding.

    Parameters
    ----------
    plate_list : list of tuples
        Each tuple is (platemap_csv, data_txt, plate_threshold, plate_var_val).
        plate_threshold is a dict {fluorophore_name: threshold}.
    fluorophores : dict
        Maps channel name (as it appears in the Harmony/Columbus export) to a
        readable name, e.g. {'Alexa 488': 'MSH2'}.
    output_csv : str or None
        Path to save the combined dataframe, or None to skip saving.
    plate_var : str or None
        Name of an extra per-plate column to add (e.g. 'treatment'), or None.
    method_nucleus : str
        Which intensity statistic to use (e.g. 'Mean', 'Median').
    header_rows : int
        Number of metadata rows before the real header row in the .txt export
        (default matches the standard Harmony/Columbus export format).

    Returns
    -------
    pd.DataFrame
        One row per nucleus, pooled across all plates. df['group_id'] is
        Date_Cells(_plate_var) - a physical plate can contribute more than one
        group_id when it holds multiple cell lines. Thresholds are still
        applied per physical plate (via plate_threshold), and are also carried
        through as a df['{name}_threshold'] column per fluorophore so they can
        be looked up later without needing plate_list again.
    """
    all_plates = []

    for i, (platemap_csv, data_txt, plate_threshold, plate_var_val) in enumerate(plate_list):
        pm = pd.read_csv(platemap_csv)
        plate_date = pm['Date'].iloc[0]
        print(f"\n{'='*60}")
        print(f"Loading plate {i+1}/{len(plate_list)}, ({plate_date}): {platemap_csv}")

        intensity_cols = [
            f'Nuclei Filt Selected - Intensity Nucleus {channel} {method_nucleus}'
            for channel in fluorophores.keys()
        ]
        area_col = 'Nuclei Filt Selected - Nucleus Area [µm²]'
        keep_cols = ['Row', 'Column', area_col] + intensity_cols

        df = pd.read_csv(data_txt, header=header_rows, sep='\t', usecols=lambda c: c in keep_cols)

        # Fail loudly if a fluorophore/area column doesn't actually exist in the file,
        # rather than silently ending up with fewer columns than expected
        missing = set(keep_cols) - set(df.columns)
        if missing:
            raise ValueError(f"{data_txt} is missing expected column(s): {missing}")

        # Row/Column dtype must match between platemap and data or the merge below
        # will silently drop every row that doesn't match
        for col in ['Row', 'Column']:
            pm[col] = pm[col].astype(df[col].dtype)

        # Merge and verify no wells were lost on either side
        pm_wells = set(zip(pm['Row'], pm['Column']))
        df_wells = set(zip(df['Row'], df['Column']))
        if pm_wells != df_wells:
            only_pm = sorted(pm_wells - df_wells)
            only_df = sorted(df_wells - pm_wells)
            raise ValueError(
                f"Row/Column mismatch between {platemap_csv} and {data_txt}. "
                f"Wells only in platemap: {only_pm[:5]}{'...' if len(only_pm) > 5 else ''}. "
                f"Wells only in data: {only_df[:5]}{'...' if len(only_df) > 5 else ''}."
            )
        df = pd.merge(pm, df, on=['Row', 'Column'], how='inner', validate='one_to_many')

        if plate_var is not None:
            df[plate_var] = plate_var_val
            print(f"added {plate_var} of {plate_var_val} to plate.")

        # Rename intensity columns and add status/threshold columns (using RAW thresholds
        # for this physical plate)
        for channel, name in fluorophores.items():
            old_col = f'Nuclei Filt Selected - Intensity Nucleus {channel} {method_nucleus}'
            new_col = f'Intensity Nucleus {name} {method_nucleus}'
            df[new_col] = df[old_col]
            df[f'{name}_pos'] = df[new_col] > plate_threshold[name]
            df[f'{name}_threshold'] = plate_threshold[name]

        df['Nuclear area'] = df[area_col]
        df.drop(columns=intensity_cols + [area_col], inplace=True)

        # group_id groups nuclei by the combination that actually matters for
        # visualisation/analysis (cell line, recovery time, etc) - Date_Cells(_plate_var).
        # A single physical plate can produce more than one group_id (e.g. half a
        # plate of one cell line, half of another).
        df['group_id'] = (
            df['Date'].astype(str) + "_" + df['Cells'].astype(str)
            + ("" if plate_var is None else "_" + df[plate_var].astype(str))
        )

        print(f"\nPlate loaded: {df.shape[0]} nuclei, {df.shape[1]} columns")
        print(f"  Date: {df['Date'].unique()}")
        print(f"  Cells: {df['Cells'].unique()}")
        print(f"  group_id(s): {df['group_id'].unique()}")
        print(f"  Conditions: {df['label'].unique().tolist()}")
        all_plates.append(df)

    combined_df = pd.concat(all_plates, ignore_index=True)
    print(f"\n{'='*60}")
    print(f"All plates loaded. In total: {combined_df.shape[0]} nuclei, {combined_df.shape[1]} columns")

    if output_csv is not None:
        combined_df.to_csv(output_csv, index=False)
        print(f"\nCombined raw data saved to {output_csv}")

    return combined_df


def plot_thresholds(df, fluorophores, method_nucleus='Mean', bins=100, ncols=3, group_col='group_id'):
    """
    Sense-check thresholding: one figure per fluorophore, with one subplot per
    group_id (3 per row) showing the intensity distribution and the threshold
    used to call it positive/negative, plus a final subplot with all group_ids
    overlaid in different colours.

    Thresholds are read from the df['{name}_threshold'] column added by
    load_multiple_plates, so this needs only df - no plate_list required.

    Parameters
    ----------
    df : pd.DataFrame
        Combined dataframe returned by load_multiple_plates.
    fluorophores : dict
        Same dict passed to load_multiple_plates (channel -> readable name).
    method_nucleus : str
        Same value passed to load_multiple_plates.
    bins : int
        Number of histogram bins.
    ncols : int
        Number of subplots per row.
    group_col : str
        Column to group/subplot by (default 'group_id').
    """
    groups = list(df[group_col].unique())
    n_subplots = len(groups) + 1  # one per group, plus one combined
    nrows = -(-n_subplots // ncols)  # ceiling division

    for name in fluorophores.values():
        col = f'Intensity Nucleus {name} {method_nucleus}'
        threshold_col = f'{name}_threshold'

        fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.5 * nrows))
        axes = axes.ravel()

        for ax, group in zip(axes, groups):
            group_df = df[df[group_col] == group]

            ax.hist(group_df[col].dropna(), bins=bins, color='steelblue')
            for threshold in group_df[threshold_col].unique():
                ax.axvline(threshold, color='red', linestyle='--', label=f'threshold = {threshold:g}')

            pct_pos = 100 * (group_df[col] > group_df[threshold_col]).mean()
            ax.set_title(f'{group}\n{pct_pos:.1f}% positive')
            ax.set_xlabel(col)
            ax.set_ylabel('Nuclei count')
            ax.legend(fontsize=8)

        # Combined subplot: all group_ids overlaid, one colour each
        combined_ax = axes[len(groups)]
        for group in groups:
            group_df = df[df[group_col] == group]
            combined_ax.hist(group_df[col].dropna(), bins=bins, histtype='step',
                              linewidth=1.5, label=group)
        combined_ax.set_title('All samples combined')
        combined_ax.set_xlabel(col)
        combined_ax.set_ylabel('Nuclei count')
        combined_ax.legend(fontsize=6)

        # Hide any unused axes in the grid
        for ax in axes[n_subplots:]:
            ax.axis('off')

        fig.suptitle(name)
        fig.tight_layout()

    plt.show()


import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def calculate_group_medians(
    input_csv,
    output_csv=None,
    fluorophores=None,
    method_nucleus='Mean',
    hierarchy=None,
    value_cols=None,
    use_normalised=False,
):
    """
    Read a CSV produced by load_multiple_plates (or a similarly structured table),
    calculate medians for the same groups used by the box plots, while always
    keeping Date separate.

    Parameters
    ----------
    input_csv : str or pd.DataFrame
        Path to the CSV exported from load_multiple_plates, or an in-memory dataframe.
    output_csv : str or None
        Path to save the median table to, or None to skip saving.
    fluorophores : dict or None
        Mapping from channel name to readable name (e.g. {'Alexa 488': 'EdU'}).
        Used to infer the intensity columns when value_cols is not supplied.
    method_nucleus : str
        Which intensity statistic to summarise (e.g. 'Mean', 'Median').
    hierarchy : list or None
        Grouping columns. Defaults to ['Cells', 'label', 'Date'], which matches
        the box-plot grouping while keeping Date separate.
    value_cols : list or None
        Specific columns to summarise. If None, these are inferred from fluorophores.
    use_normalised : bool
        If True, prefer Normalised Intensity columns when available.

    Returns
    -------
    pd.DataFrame
        One row per group with the median values for each selected intensity column.
    """
    if isinstance(input_csv, pd.DataFrame):
        df = input_csv.copy()
    else:
        df = pd.read_csv(input_csv)

    if hierarchy is None:
        hierarchy = ['Cells', 'label', 'Date']
    elif isinstance(hierarchy, str):
        hierarchy = [hierarchy]

    if 'Date' not in hierarchy:
        raise ValueError("Date must be included in hierarchy so groups stay date-separated.")

    missing_group_cols = [col for col in hierarchy if col not in df.columns]
    if missing_group_cols:
        raise ValueError(f"Input data is missing grouping column(s): {missing_group_cols}")

    if value_cols is None:
        if fluorophores is None:
            raise ValueError("Please provide fluorophores or value_cols.")

        if use_normalised:
            value_cols = [
                f'Normalised Intensity Nucleus {name} {method_nucleus}'
                for name in fluorophores.values()
            ]
            if not all(col in df.columns for col in value_cols):
                value_cols = [
                    f'Intensity Nucleus {name} {method_nucleus}'
                    for name in fluorophores.values()
                ]
        else:
            value_cols = [
                f'Intensity Nucleus {name} {method_nucleus}'
                for name in fluorophores.values()
            ]
            if not all(col in df.columns for col in value_cols):
                normalised_cols = [
                    f'Normalised Intensity Nucleus {name} {method_nucleus}'
                    for name in fluorophores.values()
                ]
                if all(col in df.columns for col in normalised_cols):
                    value_cols = normalised_cols

    if isinstance(value_cols, str):
        value_cols = [value_cols]

    missing_value_cols = [col for col in value_cols if col not in df.columns]
    if missing_value_cols:
        raise ValueError(f"Input data is missing value column(s): {missing_value_cols}")

    summary = (
        df.groupby(hierarchy, dropna=False)[value_cols]
          .median()
          .reset_index()
    )

    # Match the box-plot grouping label style
    summary['group_label'] = summary[hierarchy].astype(str).agg('_'.join, axis=1)

    if output_csv is not None:
        dirname = os.path.dirname(output_csv)
        if dirname:
            os.makedirs(dirname, exist_ok=True)
        summary.to_csv(output_csv, index=False)
        print(f"\nGroup medians saved to {output_csv}")

    return summary