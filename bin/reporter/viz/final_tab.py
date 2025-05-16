import pandas as pd
import numpy as np
import panel as pn
import plotly.graph_objects as go
import plotly.express as px

from .plot_style import style_plotly_figure
from bin.reporter.utils import shannon_entropy


def plot_final_diversity(clones_df: pd.DataFrame):
    clone_sizes = clones_df['duplicate_count']
    freqs = clone_sizes / clone_sizes.sum()

    def hill_diversity(freqs, q):
        if q == 1:
            return np.exp(shannon_entropy(freqs, base=np.e))
        return np.power((freqs ** q).sum(), 1 / (1 - q))

    hill_qs = np.linspace(0, 5, 100)
    hill_values = [hill_diversity(freqs, q) for q in hill_qs]

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=hill_qs, y=hill_values, mode='lines', line=dict(color='teal'), name='Hill'))
    fig.update_layout(title="Hill Diversity Curve", xaxis_title="q", yaxis_title="Hill Diversity",
                      width=600, height=350, showlegend=False)
    return style_plotly_figure(fig)


def plot_final_distributions(clones_df: pd.DataFrame):
    clones_df = clones_df.copy()
    clones_df['clone_size'] = clones_df['duplicate_count']
    clones_df['cdr3_len'] = clones_df['cdr3_aa'].str.len()
    clones_df['log_pgen'] = np.log10(clones_df['pgen'].replace(0, np.nan)).fillna(-100)

    p1 = px.histogram(clones_df, x='log_pgen', nbins=20, title='log10(Pgen)', color_discrete_sequence=['slateblue'])
    p1.update_layout(width=400, height=350, xaxis_title='log10(pgen)', yaxis_title='Count')

    p2 = px.histogram(clones_df, x='cdr3_len', nbins=15, title='CDR3 Length', color_discrete_sequence=['darkgreen'])
    p2.update_layout(width=400, height=350, xaxis_title='Length', yaxis_title='Count')

    p3 = px.histogram(clones_df, x='clone_size', nbins=10, title='Clone Sizes', color_discrete_sequence=['darkorange'])
    p3.update_layout(width=400, height=350, xaxis_title='Clone size', yaxis_title='Frequency')

    return pn.Row(
        pn.pane.Plotly(style_plotly_figure(p1)),
        pn.pane.Plotly(style_plotly_figure(p2)),
        pn.pane.Plotly(style_plotly_figure(p3))
    )


def create_final_tab(clones_df):
    dist = plot_final_distributions(clones_df)
    hill = plot_final_diversity(clones_df)
    return pn.Column("## Final Clone Summary", dist, hill)


def create_final_per_chain_panel(df: pd.DataFrame):
    tabs = []
    for locus in sorted(df["locus"].dropna().unique()):
        subset = df[df["locus"] == locus]

        top_n = 15
        row_genes = pn.Row(
            pn.pane.Plotly(plot_top_genes(subset, "v_call", "Top V genes", top_n), width=400),
            pn.pane.Plotly(plot_top_genes(subset, "d_call", "Top D genes", top_n),
                           width=400) if "d_call" in subset.columns and subset["d_call"].notna().any() else pn.Spacer(width=0),
            pn.pane.Plotly(plot_top_genes(subset, "j_call", "Top J genes", top_n), width=400)
        )

        top_clones = subset.nlargest(10, "duplicate_count")
        fig_clones = go.Figure(go.Bar(
            y=top_clones['cdr3_aa'],
            x=top_clones['duplicate_count'],
            orientation='h',
            marker_color='indianred'))
        fig_clones.update_layout(title="Top 10 Clones by Count", xaxis_title="Count", width=600, height=350, showlegend=False)

        diversity = plot_final_diversity(subset)
        second_summary = plot_final_distributions(subset)

        row_summary = pn.Row(pn.pane.Plotly(style_plotly_figure(fig_clones), width=500), pn.pane.Plotly(diversity, width=500))

        tab = pn.Column(pn.pane.Markdown(f"### Chain: {locus}"), row_genes, second_summary, row_summary, width=1300)
        tabs.append((locus, tab))
    return pn.Tabs(*tabs)


def plot_top_genes(df, gene_column, title, top_n=15):
    counts = df[gene_column].dropna().apply(lambda x: x.split(',')[0]).value_counts().nlargest(top_n)
    fig = go.Figure(go.Bar(y=counts.index, x=counts.values, marker_color="indianred", orientation='h'))
    fig.update_layout(title=title, yaxis_title=gene_column, xaxis_title="Count", height=650, showlegend=False)
    return style_plotly_figure(fig)
