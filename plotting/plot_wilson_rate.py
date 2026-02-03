import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from plotting.sample_config import get_colors_for_samples

def wilson_ci(count, denom, conf=0.95):
    """Calculate Wilson score interval for binomial proportion."""
    if denom == 0 or np.isnan(count) or np.isnan(denom):
        return np.nan, np.nan
    z = norm.ppf(1 - (1 - conf) / 2)
    phat = count / denom
    denom = float(denom)
    center = (phat + z**2/(2*denom)) / (1 + z**2/denom)
    margin = (z * np.sqrt(phat*(1-phat)/denom + z**2/(4*denom**2))) / (1 + z**2/denom)
    return center - margin, center + margin

def plot_wilson_rate(
    df, count_metric, denom_metric, samples, figsize=(8, 6), colors=None,
    save_path=None, show_plot=False, show_counts=True
):
    """
    Plot rates (per billion) with 95% Wilson confidence intervals.
    Error bars are straight lines.
    Optionally show value labels for counts.
    """
    fig, ax = plt.subplots(figsize=figsize, facecolor='white')
    counts = df.loc[count_metric, samples]
    denoms = df.loc[denom_metric, samples]
    rates, lowers, uppers = [], [], []
    for idx, (c, d) in enumerate(zip(counts, denoms)):
        try:
            c = float(c) if c not in ['', ' ', None] else np.nan
            d = float(d) if d not in ['', ' ', None] else np.nan
        except (ValueError, TypeError):
            c, d = np.nan, np.nan
        # print(f"Sample: {samples[idx]}, Count: {c}, Denom: {d}")
        if np.isnan(c) or np.isnan(d) or d == 0:
            rates.append(np.nan)
            lowers.append(np.nan)
            uppers.append(np.nan)
            print(f"  Skipped (NaN or zero denominator)")
        else:
            rate = (c / d) * 1e9  # per billion
            lo, up = wilson_ci(c, d)
            lo = lo * 1e9
            up = up * 1e9
            rates.append(rate)
            lowers.append(lo)
            uppers.append(up)
            # print(f"  Rate per billion: {rate:.3e}, Wilson CI: [{lo:.3e}, {up:.3e}]")
    if colors is None:
        colors = get_colors_for_samples(samples)
    else:
        colors = colors[:len(samples)]
    x_pos = np.arange(len(samples))
    bars = ax.bar(x_pos, rates, width=0.6, color=colors, edgecolor='black', linewidth=1.0, alpha=0.8)
    # Error bars: straight lines, no whiskers
    err_low = [r - l if not np.isnan(r) and not np.isnan(l) else 0 for r, l in zip(rates, lowers)]
    err_up = [u - r if not np.isnan(r) and not np.isnan(u) else 0 for r, u in zip(rates, uppers)]
    # print(f"Error bars (low): {err_low}")
    # print(f"Error bars (up): {err_up}")
    # Draw error bars as straight lines
    for i, (rate, lo, up) in enumerate(zip(rates, lowers, uppers)):
        if not np.isnan(rate) and not np.isnan(lo) and not np.isnan(up):
            ax.plot([i, i], [lo, up], color='grey', linewidth=2)
    # Value labels
    if show_counts:
        for i, (rate, c, d) in enumerate(zip(rates, counts, denoms)):
            if not np.isnan(rate):
                label = f'{rate:.0f}'
                ax.text(i, rate + (np.nanmax(rates) * 0.02 if np.nanmax(rates) else 0), label, ha='center', va='bottom', fontsize=8)
    ax.set_ylabel(f'Rate: {count_metric} / {denom_metric} (per 1e9)', fontsize=14, fontweight='bold')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(samples, fontsize=12, fontweight='bold', rotation=45, ha='right')
    ax.tick_params(axis='y', labelsize=12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    if rates and np.nanmin(rates) >= 0:
        ax.set_ylim(bottom=0)
    if rates and np.nanmax(rates) > 0:
        ax.set_ylim(top=np.nanmax(rates) * 1.15)
    ax.tick_params(axis='x', length=0)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    if show_plot:
        plt.show()
    plt.close()
    return

