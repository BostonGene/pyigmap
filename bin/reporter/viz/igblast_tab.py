import pandas as pd
import panel as pn
import plotly.graph_objects as go
import plotly.express as px

from .plot_style import style_plotly_figure


def plot_top_genes(df, gene_column, title, top_n=15):
    counts = df[gene_column].dropna().apply(lambda x: x.split(',')[0]).value_counts().nlargest(top_n)
    fig = go.Figure(go.Bar(y=counts.index, x=counts.values, marker_color="indianred", orientation='h'))
    fig.update_layout(title=title, yaxis_title=gene_column, xaxis_title="Count", height=650, showlegend=False)
    return style_plotly_figure(fig)


def plot_cdr3_length(df):
    if "junction_aa" in df.columns:
        lengths = df["junction_aa"].dropna().apply(len)
        title = "CDR3 Length (Amino Acids)"
    elif "junction_length" in df.columns:
        lengths = df["junction_length"].dropna().astype(int)
        title = "CDR3 Length (Nucleotides)"
    else:
        return style_plotly_figure(go.Figure())
    fig = px.histogram(lengths, nbins=30, title=title)
    fig.update_layout(height=350, showlegend=False)
    return style_plotly_figure(fig)


def plot_productivity(df):
    counts = df["productive"].value_counts()
    fig = px.pie(values=counts.values, names=counts.index, title="Productivity Distribution", hole=0.3)
    fig.update_traces(textinfo='label+percent', textfont_size=12, marker=dict(line=dict(color='#000000', width=1)))
    fig.update_layout(height=350, showlegend=False)
    return style_plotly_figure(fig)


def plot_v_identity(df):
    if "v_identity" not in df.columns:
        return style_plotly_figure(go.Figure())
    values = df["v_identity"].dropna().astype(float)
    fig = px.histogram(values, nbins=40, title="V-gene Alignment Identity (%)")
    fig.update_layout(height=350, showlegend=False)
    return style_plotly_figure(fig)


def create_igblast_summary_panel(df: pd.DataFrame):
    tabs = []
    for locus in sorted(df["locus"].dropna().unique()):
        subset = df[df["locus"] == locus]
        row1 = pn.Row(
            pn.pane.Plotly(plot_top_genes(subset, "v_call", "Top V genes"), width=400),
            pn.pane.Plotly(plot_top_genes(subset, "d_call", "Top D genes"), width=400) if "d_call" in subset.columns and
                                                                                          subset["d_call"].notna().any() else pn.Spacer(width=0),
            pn.pane.Plotly(plot_top_genes(subset, "j_call", "Top J genes"), width=400)
        )

        row2 = pn.Row(
            pn.pane.Plotly(plot_cdr3_length(subset), width=400),
            pn.pane.Plotly(plot_productivity(subset), width=400),
            pn.pane.Plotly(plot_v_identity(subset), width=400)
        )

        tab = pn.Column(pn.pane.Markdown(f"### Chain: {locus}"), row1, row2, width=1300)
        tabs.append((locus, tab))
    return pn.Tabs(*tabs)


def add_igblast_tab_to_report(tabs, igblast_df):
    tabs.append(("IgBlast Summary", create_igblast_summary_panel(igblast_df)))
    return tabs
