import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os


# Standard color palettes for mutation signatures
STANDARD_COLORS = {
    # SBS6 colors (standard COSMIC colors)
    'C>A': '#03BCEE',  # cyan
    'C>G': '#010101',  # black
    'C>T': '#E32926',  # red
    'T>A': '#CAC9C9',  # light gray
    'T>C': '#A1CE63',  # green
    'T>G': '#ECC6C5',  # pink
    
    # SBS96 colors - 16 shades of each SBS6 color (ordered by trinucleotide context)
    # C>A (cyan shades) - A[C>A]A through T[C>A]T
    'A[C>A]A': '#9EEAF9', 'A[C>A]C': '#8AE6F7', 'A[C>A]G': '#76E2F5', 'A[C>A]T': '#62DEF4',
    'C[C>A]A': '#4EDAF2', 'C[C>A]C': '#3AD6F0', 'C[C>A]G': '#26D2EE', 'C[C>A]T': '#12CEEC',
    'G[C>A]A': '#03BCEE', 'G[C>A]C': '#03AAD6', 'G[C>A]G': '#0398BE', 'G[C>A]T': '#0286A6',
    'T[C>A]A': '#02748E', 'T[C>A]C': '#026276', 'T[C>A]G': '#01505E', 'T[C>A]T': '#013E46',
    
    # C>G (black/gray shades) - A[C>G]A through T[C>G]T
    'A[C>G]A': '#E8E8E8', 'A[C>G]C': '#D1D1D1', 'A[C>G]G': '#BABABA', 'A[C>G]T': '#A3A3A3',
    'C[C>G]A': '#8C8C8C', 'C[C>G]C': '#757575', 'C[C>G]G': '#5E5E5E', 'C[C>G]T': '#474747',
    'G[C>G]A': '#303030', 'G[C>G]C': '#282828', 'G[C>G]G': '#202020', 'G[C>G]T': '#181818',
    'T[C>G]A': '#101010', 'T[C>G]C': '#0C0C0C', 'T[C>G]G': '#080808', 'T[C>G]T': '#040404',
    
    # C>T (red shades) - A[C>T]A through T[C>T]T
    'A[C>T]A': '#FFB3B1', 'A[C>T]C': '#FF9C9A', 'A[C>T]G': '#FF8583', 'A[C>T]T': '#FF6E6C',
    'C[C>T]A': '#FF5755', 'C[C>T]C': '#FF403E', 'C[C>T]G': '#F73632', 'C[C>T]T': '#E32926',
    'D[C>T]A': '#D02522', 'G[C>T]C': '#BD211E', 'G[C>T]G': '#AA1D1A', 'G[C>T]T': '#971916',
    'T[C>T]A': '#841512', 'T[C>T]C': '#71110E', 'T[C>T]G': '#5E0D0A', 'T[C>T]T': '#4B0906',
    
    # T>A (light gray shades) - A[T>A]A through T[T>A]T
    'A[T>A]A': '#FFFFFF', 'A[T>A]C': '#F5F5F5', 'A[T>A]G': '#EBEBEB', 'A[T>A]T': '#E1E1E1',
    'C[T>A]A': '#D7D7D7', 'C[T>A]C': '#CDCDCD', 'C[T>A]G': '#C3C3C3', 'C[T>A]T': '#B9B9B9',
    'G[T>A]A': '#CAC9C9', 'G[T>A]C': '#BEBDBD', 'G[T>A]G': '#B2B1B1', 'G[T>A]T': '#A6A5A5',
    'T[T>A]A': '#9A9999', 'T[T>A]C': '#8E8D8D', 'T[T>A]G': '#828181', 'T[T>A]T': '#767575',
    
    # T>C (green shades) - A[T>C]A through T[T>C]T
    'A[T>C]A': '#E1F5C9', 'A[T>C]C': '#D7F1BA', 'A[T>C]G': '#CDEDAB', 'A[T>C]T': '#C3E99C',
    'C[T>C]A': '#B9E58D', 'C[T>C]C': '#AFE17E', 'C[T>C]G': '#A5DD6F', 'C[T>C]T': '#9BD960',
    'G[T>C]A': '#A1CE63', 'G[T>C]C': '#91B959', 'G[T>C]G': '#81A44F', 'G[T>C]T': '#718F45',
    'T[T>C]A': '#617A3B', 'T[T>C]C': '#516531', 'T[T>C]G': '#415027', 'T[T>C]T': '#313B1D',
    
    # T>G (pink shades) - A[T>G]A through T[T>G]T
    'A[T>G]A': '#FBEEF0', 'A[T>G]C': '#F9E3E8', 'A[T>G]G': '#F7D8E0', 'A[T>G]T': '#F5CDD8',
    'C[T>G]A': '#F3C2D0', 'C[T>G]C': '#F1B7C8', 'C[T>G]G': '#EFACC0', 'C[T>G]T': '#EDA1B8',
    'G[T>G]A': '#ECC6C5', 'G[T>G]C': '#D5B2B1', 'G[T>G]G': '#BE9E9D', 'G[T>G]T': '#A78A89',
    'T[T>G]A': '#907675', 'T[T>G]C': '#796261', 'T[T>G]G': '#624E4D', 'T[T>G]T': '#4B3A39',
    
    # DBS colors
    'AC>NN': '#03BCEE',
    'AT>NN': '#010101',
    'CC>NN': '#E32926',
    'CG>NN': '#CAC9C9',
    'CT>NN': '#A1CE63',
    'GC>NN': '#ECC6C5',
    'TA>NN': '#EDBF69',
    'TC>NN': '#E770A0',
    'TG>NN': '#B0DD8B',
    'TT>NN': '#36A12E',
    
    # ID colors
    'DEL_C': '#FDBE6F',   # orange (C deletions)
    'DEL_T': '#FF8001',   # dark orange (T deletions)
    'INS_C': '#B0DD8B',   # light green (C insertions)
    'INS_T': '#36A12E',   # dark green (T insertions)
    'DEL_long': '#FC8CC3', # pink (long deletions)
    'INS_long': '#D43C96', # magenta (long insertions)
    'DEL_MH': '#8E038E',   # purple (microhomology)
    'COMPLEX': '#958F8F',  # gray (complex)
}


def get_mutation_color(mutation_type):
    """
    Get standard color for a mutation type.
    
    Args:
        mutation_type: Mutation type string (e.g., 'C>A', 'T:C>A', '1:Del:C:5')
    
    Returns:
        str: Hex color code
    """
    # Direct match (includes SBS6 and SBS96)
    if mutation_type in STANDARD_COLORS:
        return STANDARD_COLORS[mutation_type]
    
    # SBS18/24/384/1536 - extract base substitution or trinucleotide
    if '>' in mutation_type:
        # Try to match full trinucleotide context first (e.g., 'A[C>A]T')
        if '[' in mutation_type and ']' in mutation_type:
            if mutation_type in STANDARD_COLORS:
                return STANDARD_COLORS[mutation_type]
            # Extract just substitution as fallback
            sub = mutation_type.split('[')[1].split(']')[0]  # 'A[C>A]T' -> 'C>A'
        # T:C>A format
        elif ':' in mutation_type:
            sub = mutation_type.split(':')[1]  # 'T:C>A' -> 'C>A'
        else:
            sub = mutation_type
        
        if sub in STANDARD_COLORS:
            return STANDARD_COLORS[sub]
    
    # ID signatures
    if 'Del' in mutation_type or 'DEL' in mutation_type:
        if 'C' in mutation_type:
            return STANDARD_COLORS['DEL_C']
        elif 'T' in mutation_type:
            return STANDARD_COLORS['DEL_T']
        elif 'long' in mutation_type.lower():
            return STANDARD_COLORS['DEL_long']
    
    if 'Ins' in mutation_type or 'INS' in mutation_type:
        if 'C' in mutation_type:
            return STANDARD_COLORS['INS_C']
        elif 'T' in mutation_type:
            return STANDARD_COLORS['INS_T']
        elif 'long' in mutation_type.lower():
            return STANDARD_COLORS['INS_long']
    
    if 'MH' in mutation_type or 'microhomology' in mutation_type.lower():
        return STANDARD_COLORS['DEL_MH']
    
    if 'complex' in mutation_type.lower():
        return STANDARD_COLORS['COMPLEX']
    
    # Default fallback
    return '#808080'  # gray


def get_default_colors(mutation_types):
    """
    Generate color palette for mutation types using standard colors.
    Falls back to colormap if too many types don't have standard colors.
    
    Args:
        mutation_types: List of mutation type strings
        
    Returns:
        list: List of color values
    """
    colors = [get_mutation_color(mt) for mt in mutation_types]
    
    # Check if too many fell back to gray (unknown types)
    gray_count = sum(1 for c in colors if c == '#808080')
    
    # If more than 30% are gray, use colormap instead
    if len(colors) > 0 and gray_count / len(colors) > 0.3:
        print(f"[INFO] Many mutation types without standard colors, using colormap")
        if len(mutation_types) <= 20:
            cmap = plt.cm.get_cmap('tab20')
        elif len(mutation_types) <= 40:
            cmap = plt.cm.get_cmap('tab20b')
        else:
            cmap = plt.cm.get_cmap('viridis')
        return [cmap(i / len(mutation_types)) for i in range(len(mutation_types))]
    
    return colors


def plot_stacked_signatures(
    matrix_path,
    coverage_path,
    samples=None,
    mutation_types=None,
    background_subtraction=None,
    figsize=None,
    colors=None,
    save_path=None,
    show_plot=True,
    show_legend=True,
    title=None,
    ylabel=None,
    show_stack_labels=False,
    stack_label_threshold=0.02,
    stack_label_fontsize=6
):
    """
    Plot mutation signatures as stacked bar charts, normalized by coverage.
    
    Args:
        matrix_path: Path to mutation matrix file (e.g., sigpro.SBS6.all)
        coverage_path: Path to coverage summary CSV (required)
        samples: List of sample names to plot (if None, uses all samples in matrix)
        mutation_types: List of mutation types to include (if None, uses all types)
        background_subtraction: Dict mapping sample->background (e.g., {'treated': 'control'})
        figsize: Figure size tuple (if None, auto-calculated based on sample count)
        colors: List of colors for mutation types (if None, uses standard colors)
        save_path: Path to save the figure
        show_plot: Whether to display the plot
        show_legend: Whether to show legend
        title: Plot title (if None, inferred from matrix filename)
        ylabel: Y-axis label (if None, defaults to 'Mutations per billion bases')
        show_stack_labels: Whether to show labels on each stack segment
        stack_label_threshold: Minimum proportion of total height to show label (default 0.02 = 2%)
        stack_label_fontsize: Font size for stack labels (default 6)
        
    Returns:
        tuple: (fig, ax) matplotlib objects
    """
    # Read mutation matrix
    print(f"[INFO] Reading mutation matrix from {matrix_path}")
    matrix_df = pd.read_csv(matrix_path, sep='\t', index_col='MutationType')
    
    # Filter mutation types if specified
    if mutation_types is not None:
        available_types = [mt for mt in mutation_types if mt in matrix_df.index]
        missing_types = [mt for mt in mutation_types if mt not in matrix_df.index]
        
        if missing_types:
            print(f"[WARNING] Mutation types not found: {missing_types}")
        
        if not available_types:
            print(f"[ERROR] No valid mutation types found")
            return None, None
        
        matrix_df = matrix_df.loc[available_types]
        print(f"[INFO] Filtered to {len(available_types)} mutation types")
    
    # Determine samples to plot
    if samples is None:
        samples = matrix_df.columns.tolist()
    else:
        samples = [s for s in samples if s in matrix_df.columns]
        missing = [s for s in (samples if samples else []) if s not in matrix_df.columns]
        if missing:
            print(f"[WARNING] Samples not found: {missing}")
    
    if not samples:
        print(f"[ERROR] No valid samples found")
        return None, None
    
    # Filter out background samples from plotting if background subtraction is enabled
    if background_subtraction:
        background_samples = set(background_subtraction.values())
        plot_samples = [s for s in samples if s not in background_samples]
        print(f"[INFO] Background subtraction enabled, excluding {len(background_samples)} background samples")
        print(f"[INFO] Background pairs: {background_subtraction}")
    else:
        plot_samples = samples
    
    # Auto-calculate figure width based on sample count if not provided
    if figsize is None:
        width = 2 + (len(plot_samples) * 0.8)
        height = 6
        figsize = (width, height)
        print(f"[INFO] Auto-calculated figure size: {figsize}")
    
    plot_data = matrix_df[samples].copy()
    
    # Normalize by coverage (always)
    print(f"[INFO] Reading coverage from {coverage_path}")
    coverage_df = pd.read_csv(coverage_path, index_col='sample')
    
    print(f"[INFO] Normalizing by coverage (per billion bases)")
    for sample in samples:
        if sample in coverage_df.index:
            try:
                coverage = coverage_df.loc[sample, 'total_coverage']
            except:
                coverage = coverage_df.loc[sample, 'effective_coverage']
            plot_data[sample] = (plot_data[sample] / coverage) * 1e9
        else:
            print(f"[WARNING] Coverage not found for {sample}")
    
    # Apply background subtraction if specified
    if background_subtraction:
        print(f"[INFO] Applying background subtraction")
        for treated_sample, control_sample in background_subtraction.items():
            if treated_sample in plot_data.columns and control_sample in plot_data.columns:
                plot_data[treated_sample] = plot_data[treated_sample] - plot_data[control_sample]
                # Set negative values to zero (can't have negative mutations)
                plot_data[treated_sample] = plot_data[treated_sample].clip(lower=0)
                print(f"[INFO] Subtracted {control_sample} from {treated_sample}")
            else:
                print(f"[WARNING] Cannot subtract: {treated_sample} or {control_sample} not found")
    
    # Select only samples to plot (excluding backgrounds)
    plot_data = plot_data[plot_samples]
    
    mutation_type_list = plot_data.index.tolist()
    print(f"[INFO] Plotting {len(plot_samples)} samples with {len(mutation_type_list)} mutation types")
    
    # Generate colors using standard palette if not provided
    if colors is None:
        colors = get_default_colors(mutation_type_list)
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize, facecolor='white')
    
    # Create stacked bar chart
    x_pos = np.arange(len(plot_samples))
    bottom = np.zeros(len(plot_samples))
    
    for idx, mut_type in enumerate(mutation_type_list):
        values = plot_data.loc[mut_type].values
        ax.bar(x_pos, values, width=0.6, bottom=bottom,
               color=colors[idx], edgecolor='white', linewidth=0.5,
               label=mut_type, alpha=0.9)
        
        # Add labels to stack segments if requested
        if show_stack_labels:
            for i, (x, val, bot) in enumerate(zip(x_pos, values, bottom)):
                # Calculate total height for this sample
                total_height = plot_data[plot_samples[i]].sum()
                # Only label if segment is large enough (above threshold)
                if total_height > 0 and val / total_height >= stack_label_threshold:
                    # Position label in middle of segment
                    y_pos = bot + val / 2
                    # Choose text color based on background brightness
                    color_rgb = plt.matplotlib.colors.to_rgb(colors[idx])
                    brightness = 0.299 * color_rgb[0] + 0.587 * color_rgb[1] + 0.114 * color_rgb[2]
                    text_color = 'white' if brightness < 0.5 else 'black'
                    # Add text label
                    ax.text(x, y_pos, mut_type, ha='center', va='center',
                           fontsize=stack_label_fontsize, fontweight='bold',
                           color=text_color, rotation=0)
        
        bottom += values
    
    # Formatting
    if ylabel is None:
        ylabel = 'Mutations per billion bases'
        if background_subtraction:
            ylabel += ' (delta)'
    ax.set_ylabel(ylabel, fontsize=14, fontweight='bold')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(plot_samples, fontsize=12, fontweight='bold', rotation=45, ha='right')
    ax.tick_params(axis='y', labelsize=12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.tick_params(axis='x', length=0)
    ax.set_ylim(bottom=0)
    
    # Add title
    if title:
        ax.set_title(title, fontsize=16, fontweight='bold', pad=20)
    else:
        matrix_type = os.path.basename(matrix_path).split('.')[1]
        title_text = f'{matrix_type} Mutation Signatures'
        if background_subtraction:
            title_text += ' (delta)'
        ax.set_title(title_text, fontsize=16, fontweight='bold', pad=20)
    
    # Add legend
    if show_legend and len(mutation_type_list) <= 30:
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10, frameon=False)
    
    plt.tight_layout()
    
    # Save figure
    if save_path:
        print(f"[INFO] Saving plot to {save_path}")
        os.makedirs(os.path.dirname(save_path) if os.path.dirname(save_path) else '.', exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
        print(f"[SUCCESS] Plot saved")
    
    # Show plot
    if show_plot:
        plt.show()
    else:
        plt.close()
    
    return fig, ax
