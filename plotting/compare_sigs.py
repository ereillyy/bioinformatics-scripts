import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics.pairwise import cosine_similarity


def top_reference_matches(sig, reference_sigs, n=10):
    """
    Return the top-n most cosine-similar signatures from a reference catalogue.

    Parameters
    ----------
    sig : pd.Series
        The query signature (index = mutation types / features).
    reference_sigs : pd.DataFrame
        Reference catalogue where rows are features and columns are signatures.
    n : int
        Number of top matches to return.

    Returns
    -------
    pd.Series
        Cosine similarities, sorted descending, length <= n.
    """
    aligned = sig.reindex(reference_sigs.index)
    sims = cosine_similarity(aligned.values.reshape(1, -1), reference_sigs.T.values)[0]
    return pd.Series(sims, index=reference_sigs.columns).sort_values(ascending=False).head(n)


def plot_top_matches(
    matches,
    title=None,
    ax=None,
    no_labels=False,
    annotate="line",
    threshold=0.8,
):
    """
    Bar plot of cosine-similarity matches against a reference catalogue.

    Parameters
    ----------
    matches : pd.Series
        Output of top_reference_matches (or any similarity Series).
    title : str, optional
    ax : matplotlib Axes, optional
        Draws on a new figure if not supplied.
    no_labels : bool
        If True, suppresses the y-axis label and title.
    annotate : {'line', 'value', None}
        'line'  – draw a horizontal threshold line.
        'value' – print each bar's value above it.
        None    – no annotation.
    threshold : float
        Value for the threshold line (only used when annotate='line').

    Returns
    -------
    matplotlib Axes
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(0.3 * len(matches) + 0.5, 1.5))

    bars = ax.bar(matches.index, matches.values, color="#818181")

    if annotate == "line":
        ax.axhline(threshold, color="black", linestyle="--", linewidth=1)
    elif annotate == "value":
        ax.bar_label(bars, fmt="%.2f", padding=2, fontsize=8)

    ax.set_ylim(0, 1.15)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    # ax.set_yticks([0.1, 0.3, 0.5, 0.7, 0.9], minor=True)  # Keep in case want to use 0.9 as the threshold for a dashed line
    plt.setp(ax.get_xticklabels(), rotation=90, ha="center")

    ax.set_ylabel("Cosine similarity")
    if not no_labels:
        if title:
            ax.set_title(title)

    ax.spines[["top", "right"]].set_visible(False)
    return ax

