import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb


# ===========================================================================
# Colour schemes
# ===========================================================================

SBS_COLORS = {
    'C>A': '#03BCEE', 'C>G': '#010101', 'C>T': '#E32926',
    'T>A': '#CAC9C9', 'T>C': '#A1CE63', 'T>G': '#EBC6C4',
}
SBS_MUT_ORDER = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
BASES = ['A', 'C', 'G', 'T']

DBS_COLORS = {
    'AC': '#A6CEE3',  # light blue
    'AT': '#1F78B4',  # dark blue
    'CC': '#B2DF8A',  # light green
    'CG': '#33A02C',  # dark green
    'CT': '#FB9A99',  # pink/salmon
    'GC': '#E31A1C',  # red
    'TA': '#FDBF6F',  # light orange
    'TC': '#FF7F00',  # orange
    'TG': '#CAB2D6',  # light purple
    'TT': '#6A3D9A',  # dark purple
}
DBS_REF_ORDER = ['AC', 'AT', 'CC', 'CG', 'CT', 'GC', 'TA', 'TC', 'TG', 'TT']# Alt dinucleotides per reference, in standard COSMIC DBS78 order.
DBS_ALTS = {
    'AC': ['CA', 'CG', 'CT', 'GA', 'GG', 'GT', 'TA', 'TG', 'TT'],
    'AT': ['CA', 'CC', 'CG', 'GA', 'GC', 'TA'],
    'CC': ['AA', 'AG', 'AT', 'GA', 'GG', 'GT', 'TA', 'TG', 'TT'],
    'CG': ['AT', 'GC', 'GT', 'TA', 'TC', 'TT'],
    'CT': ['AA', 'AC', 'AG', 'GA', 'GC', 'GG', 'TA', 'TC', 'TG'],
    'GC': ['AA', 'AG', 'AT', 'CA', 'CG', 'TA'],
    'TA': ['AT', 'CG', 'CT', 'GC', 'GG', 'GT'],
    'TC': ['AA', 'AG', 'AT', 'CA', 'CG', 'CT', 'GA', 'GG', 'GT'],
    'TG': ['AA', 'AC', 'AT', 'CA', 'CC', 'CT', 'GA', 'GC', 'GT'],
    'TT': ['AA', 'AC', 'AG', 'CA', 'CC', 'CG', 'GA', 'GC', 'GG'],
}

# Canonical COSMIC/SigProfilerPlotting ID83 colour scheme.
ID_COLORS = {
    'DEL_C': '#FDBE6F', 'DEL_T': '#FF8002',
    'INS_C': '#B0DD8B', 'INS_T': '#36A12E',
}
ID_DEL_REPEAT_COLORS = ['#FDCAB5', '#FC8A6A', '#F14432', '#BC191A']
ID_INS_REPEAT_COLORS = ['#D0E1F2', '#94C4DF', '#4A98C9', '#1764AB']
ID_MH_COLORS = ['#E2E2EF', '#B6B6D8', '#8683BD', '#62409B']


# ===========================================================================
# Per-type order/label/group builders
# ===========================================================================

def _sbs96_order_and_groups():
    order = []
    for mut in SBS_MUT_ORDER:
        for five in BASES:
            for three in BASES:
                order.append(f"{five}[{mut}]{three}")
    labels = [t[0] + t[2] + t[6] for t in order]
    groups = []
    for mut in SBS_MUT_ORDER:
        types = [t for t in order if f"[{mut}]" in t]
        groups.append((mut, types, SBS_COLORS[mut]))
    return order, labels, groups, None


def _dbs78_order_and_groups():
    order, labels, groups = [], [], []
    for ref in DBS_REF_ORDER:
        types = [f"{ref}>{alt}" for alt in DBS_ALTS[ref]]
        order.extend(types)
        labels.extend(DBS_ALTS[ref])
        groups.append((f"{ref}>NN", types, DBS_COLORS[ref]))
    return order, labels, groups, None

def _id83_order_and_groups():
    order, labels, groups = [], [], []

    def add_run(types_labels, glabel, color):
        types, labs = zip(*types_labels)
        order.extend(types)
        labels.extend(labs)
        groups.append((glabel, list(types), color))

    # Deletion homopolymer/repeat length is displayed 1-indexed (1..6+, COSMIC
    # convention), one more than the 0-indexed count stored in the type string.
    add_run([(f"1:Del:C:{n}", f"{n + 1}+" if n == 5 else str(n + 1)) for n in range(0, 6)],
             'C', ID_COLORS['DEL_C'])
    add_run([(f"1:Del:T:{n}", f"{n + 1}+" if n == 5 else str(n + 1)) for n in range(0, 6)],
             'T', ID_COLORS['DEL_T'])
    add_run([(f"1:Ins:C:{n}", f"{n}+" if n == 5 else str(n)) for n in range(0, 6)],
             'C', ID_COLORS['INS_C'])
    add_run([(f"1:Ins:T:{n}", f"{n}+" if n == 5 else str(n)) for n in range(0, 6)],
             'T', ID_COLORS['INS_T'])

    for i, length in enumerate([2, 3, 4, 5]):
        glabel = f"{length}{'+' if length == 5 else ''}"
        add_run([(f"{length}:Del:R:{n}", f"{n + 1}+" if n == 5 else str(n + 1)) for n in range(0, 6)],
                 glabel, ID_DEL_REPEAT_COLORS[i])
    for i, length in enumerate([2, 3, 4, 5]):
        glabel = f"{length}{'+' if length == 5 else ''}"
        add_run([(f"{length}:Ins:R:{n}", f"{n}+" if n == 5 else str(n)) for n in range(0, 6)],
                 glabel, ID_INS_REPEAT_COLORS[i])

    for i, length in enumerate([2, 3, 4, 5]):
        n_mh = length - 1 if length < 5 else 5
        glabel = f"{length}{'+' if length == 5 else ''}"
        add_run([(f"{length}:Del:M:{m}", f"{m}+" if (length == 5 and m == 5) else str(m))
                  for m in range(1, n_mh + 1)],
                 glabel, ID_MH_COLORS[i])

    # Broader category headers, each spanning a fixed number of the groups
    # above, mirroring the two-tier header used in canonical COSMIC ID83 plots.
    supergroups = [
        ("1bp Deletion", 2),
        ("1bp Insertion", 2),
        (">1bp Deletion at Repeats\n(Deletion Length)", 4),
        (">1bp Insertion at Repeats\n(Insertion Length)", 4),
        ("Microhomology\n(Deletion Length)", 4),
    ]

    return order, labels, groups, supergroups


_CATALOGUE_BUILDERS = {
    'SBS96': _sbs96_order_and_groups,
    'DBS78': _dbs78_order_and_groups,
    'ID83': _id83_order_and_groups,
}
_DEFAULT_FIGSIZE = {
    'SBS96': (9, 3.5),
    'DBS78': (9, 3.5),
    'ID83': (10.5, 3.5),
}


# ===========================================================================
# Shared helpers
# ===========================================================================

def _check_counts(df):
    """Reject input that looks like proportions or percentages (non-integer columns summing to 1 or 100) or is negative.
    Any other counts are accepted, including non-integer (context-normalised) ones."""
    for col in df.columns:
        values = df[col].values
        total = values.sum()
        if np.any(values < 0):
            raise ValueError(f"Column '{col}' has negative values; expected non-negative counts.")
        if not np.all(np.isclose(values, np.round(values))) and (np.isclose(total, 1) or np.isclose(total, 100)):
            raise ValueError(
                f"Column '{col}' looks like proportions or percentages (sums to {total:g}); "
                f"expected counts. Pass proportion=True if this data is already normalized."
            )


def _normalize_to_coverage(df, coverage_path):
    if isinstance(coverage_path, pd.DataFrame):
        print("[INFO] Using the coverage passed as a DataFrame")
        coverage_df = coverage_path
    else:
        print(f"[INFO] Reading coverage from {coverage_path}")
        coverage_df = pd.read_csv(coverage_path, index_col='sample')
    print(f"[INFO] Normalizing by coverage (per billion bases)")
    for col in df.columns:
        if col in coverage_df.index:
            coverage = coverage_df.loc[col, 'effective_coverage']
            df[col] = (df[col] / coverage) * 1e9
        else:
            print(f"[WARNING] Coverage not found for {col}")
    return df


def _normalize_to_proportion(df):
    for col in df.columns:
        total = df[col].sum()
        if not np.isclose(total, 1.0):
            df[col] = df[col] / total
    return df


def _apply_background_subtraction(df, background_subtraction):
    print(f"[INFO] Applying background subtraction")
    for sample, background in background_subtraction.items():
        if sample in df.columns and background in df.columns:
            df[sample] = df[sample] - df[background]
            print(f"[INFO] Subtracted {background} from {sample}")
        else:
            print(f"[WARNING] Cannot subtract: {sample} or {background} not found")
    background_samples = set(background_subtraction.values())
    return [c for c in df.columns if c not in background_samples]


def _shared_ylim(df, plot_cols):
    all_values = df[plot_cols].values
    data_min = min(0.0, all_values.min())
    data_max = max(0.0, all_values.max())
    margin = (data_max - data_min) * 0.05 if data_max > data_min else 1.0
    ymin = data_min - margin if data_min < 0 else 0.0
    ymax = data_max + margin
    return ymin, ymax


# ===========================================================================
# Main plotting function
# ===========================================================================

def plot_mutational_signatures(
    df,
    mut_type='SBS96',
    figsize=None,
    ylabel=None,
    background_subtraction=None,
    proportion=False,
    coverage_path=None,
    fix_y_axis=True,
    save_path=None,
    show_plot=True,
):
    """
    Plot mutational signature bar plots (one plot per sample/signature),
    for SBS96, DBS78, or ID83 catalogues.

    Parameters
    ----------
    df : pd.DataFrame
        Index = mutation type strings for the chosen `mut_type`, columns = samples.
    mut_type : {'SBS96', 'DBS78', 'ID83'}
    figsize : tuple, optional
        Defaults to a sensible per-type size if not given.
    ylabel : str, optional
        Inferred from `proportion`/background_subtraction if not given.
    background_subtraction : dict, optional
        Maps sample -> background sample. Subtracted before plotting, values
        NOT clipped at zero (bars can go negative). Background columns are
        excluded from plotting.
    proportion : bool, default False
        If False, `df` must be raw counts and `coverage_path` is required;
        normalized to mutations per billion bases using 'effective_coverage'.
        Errors if the data doesn't look like counts.
        If True, each column is normalized to sum to 1 if it doesn't already.
    coverage_path : str or DataFrame, optional
        Required when proportion=False. A path to a CSV, or a DataFrame, indexed by 'sample'
        with an 'effective_coverage' column.
    fix_y_axis : bool, default True
        Share the same y-axis range across all plotted samples.
    save_path : str, optional
        Folder path. This folder is created (if it doesn't exist), and every
        sample's plot is saved inside it as '{sample}_{foldername}.png',
        where foldername is the last component of save_path.
    show_plot : bool, default True
        Whether to display each plot.

    Returns
    -------
    dict of {column_name: (fig, ax)}
    """
    if mut_type not in _CATALOGUE_BUILDERS:
        raise ValueError(f"mut_type must be one of {list(_CATALOGUE_BUILDERS)}, got '{mut_type}'")

    order, tick_labels, groups, supergroups = _CATALOGUE_BUILDERS[mut_type]()
    n = len(order)

    missing = set(order) - set(df.index)
    if missing:
        raise ValueError(
            f"Missing {len(missing)} mutation types for {mut_type}, e.g. "
            f"{list(missing)[:5]}. Check that your data's index format matches "
            f"the expected {mut_type} type strings."
        )
    df = df.loc[order].copy()

    user_ylabel = ylabel

    if not proportion:
        if coverage_path is None:
            raise ValueError("coverage_path is required when proportion=False")
        _check_counts(df)
        df = _normalize_to_coverage(df, coverage_path)
        if user_ylabel is None:
            ylabel = 'Mutations per billion bases'
    else:
        df = _normalize_to_proportion(df)
        if user_ylabel is None:
            ylabel = 'Mutation proportion'

    if background_subtraction:
        plot_cols = _apply_background_subtraction(df, background_subtraction)
        if user_ylabel is None:
            ylabel = ylabel + ' (delta)'
    else:
        plot_cols = df.columns.tolist()

    colors = [None] * n
    for _, types, color in groups:
        for t in types:
            colors[order.index(t)] = color

    if figsize is None:
        figsize = _DEFAULT_FIGSIZE[mut_type]

    fixed_ymin = fixed_ymax = None
    if fix_y_axis:
        fixed_ymin, fixed_ymax = _shared_ylim(df, plot_cols)

    group_positions = []
    cum = 0
    for glabel, types, color in groups:
        group_positions.append((cum, len(types), color, glabel))
        cum += len(types)

    # Set up save folder, named after the last component of save_path
    if save_path:
        save_path = save_path.rstrip('/')
        folder_name = os.path.basename(save_path)
        os.makedirs(save_path, exist_ok=True)
        print(f"[INFO] Saving plots to {save_path}/")

    figs = {}
    for col in plot_cols:
        values = df[col].values
        fig, ax = plt.subplots(figsize=figsize)
        x = np.arange(n)
        ax.bar(x, values, color=colors, width=0.8, edgecolor='none')

        ax.set_xticks(x)
        ax.set_xticklabels(tick_labels, rotation=90, fontsize=6, family='monospace')
        ax.set_xlim(-0.5, n - 0.5)
        ax.set_ylabel(ylabel)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        if fix_y_axis:
            ymin, ymax = fixed_ymin, fixed_ymax
        else:
            ymin, ymax = ax.get_ylim()
            ymin = min(ymin, 0.0)
        ax.set_ylim(ymin, ymax)

        if ymin < 0:
            ax.axhline(0, color='black', linewidth=0.8, zorder=3)

        # coloured strip + short group labels above the plot. Drawn in axes coordinates (x in data, y as a fraction
        # of the axes height) so the strip is the same size whatever the y-range, including negative values.
        strip_y, strip_h = 1.03, 0.11
        label_fontsize = ax.yaxis.label.get_size()
        for start, size, color, glabel in group_positions:
            ax.add_patch(plt.Rectangle(
                (start - 0.5, strip_y), size, strip_h,
                color=color, clip_on=False, transform=ax.get_xaxis_transform()
            ))
            r, g, b = to_rgb(color)
            text_color = 'white' if 0.299 * r + 0.587 * g + 0.114 * b < 0.5 else 'black'      # white on dark colours
            ax.text(start + size / 2 - 0.5, strip_y + strip_h * 0.5, glabel,
                    ha='center', va='center', fontsize=label_fontsize, fontweight='bold', color=text_color,
                    transform=ax.get_xaxis_transform())

        # Broader category header spanning several groups (ID83 only).
        if supergroups:
            header_y = strip_y + strip_h + 0.04
            gi = 0
            for slabel, n_spanned in supergroups:
                spanned = group_positions[gi:gi + n_spanned]
                sstart = spanned[0][0]
                swidth = sum(size for _, size, _, _ in spanned)
                ax.text(sstart + swidth / 2 - 0.5, header_y, slabel,
                        ha='center', va='bottom', fontsize=10, fontweight='bold',
                        transform=ax.get_xaxis_transform())
                gi += n_spanned
        ax.set_ylim(ymin, ymax)

        # Title sits above the strip (and, for ID83, the extra header row
        # above it), not overlapping them.
        title_pad = 70 if supergroups else 48
        ax.set_title(col, fontsize=12, fontweight='bold', pad=title_pad)

        fig.tight_layout()
        fig.subplots_adjust(top=0.68 if supergroups else 0.78)

        if save_path:
            out_file = os.path.join(save_path, f"{col}_{folder_name}.png")
            fig.savefig(out_file, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
            print(f"[INFO] Saved {out_file}")

        if show_plot:
            plt.show()
        else:
            plt.close(fig)

        figs[col] = (fig, ax)

    return figs


# ===========================================================================
# Convenience wrappers
# ===========================================================================

def plot_sbs96(df, **kwargs):
    return plot_mutational_signatures(df, mut_type='SBS96', **kwargs)


def plot_dbs78(df, **kwargs):
    return plot_mutational_signatures(df, mut_type='DBS78', **kwargs)


def plot_id83(df, **kwargs):
    return plot_mutational_signatures(df, mut_type='ID83', **kwargs)