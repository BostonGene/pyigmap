import pandas as pd
import panel as pn
import plotly.graph_objects as go
from bokeh.plotting import figure
from .plot_style import style_plotly_figure

def plot_locus_barplot(df: pd.DataFrame, locus_column: str, title: str, color: str = "steelblue"):
    counts = df[locus_column].value_counts().sort_values(ascending=False)
    categories = list(counts.index)
    p = figure(x_range=categories, height=450, width=600, title=title, toolbar_location=None, tools="")
    p.vbar(x=counts.index, top=counts.values, width=0.8, color=color)
    p.xgrid.grid_line_color = None
    p.y_range.start = 0
    p.xaxis.major_label_orientation = "vertical"
    return p, counts.max()

def plot_locus_distributions(vidjil_df: pd.DataFrame, igblast_df: pd.DataFrame):
    p1, max1 = plot_locus_barplot(vidjil_df, "locus", "Vidjil: chain distribution", color="steelblue")
    p2, max2 = plot_locus_barplot(igblast_df, "locus", "IgBlast: chain distribution", color="darkorange")
    y_max = max(max1, max2) * 1.05
    p1.y_range.end = y_max
    p2.y_range.end = y_max
    return pn.Row(p1, p2)

def build_chain_sankey(vidjil_df, igblast_df):
    merged = pd.DataFrame({
        "sequence_id": vidjil_df["sequence_id"],
        "vidjil_locus": vidjil_df["locus"]
    }).merge(
        igblast_df[["sequence_id", "locus"]], on="sequence_id", how="inner"
    ).rename(columns={"locus": "igblast_locus"})

    counts = merged.groupby(["vidjil_locus", "igblast_locus"]).size().reset_index(name="value")
    counts["source"] = "Vidjil: " + counts["vidjil_locus"]
    counts["target"] = "IgBlast: " + counts["igblast_locus"]

    all_labels = list(pd.unique(counts[["source", "target"]].values.ravel()))
    label_map = {l: i for i, l in enumerate(all_labels)}

    fig = go.Figure(data=[go.Sankey(
        node=dict(label=all_labels),
        link=dict(
            source=counts["source"].map(label_map),
            target=counts["target"].map(label_map),
            value=counts["value"]
        )
    )])

    fig.update_layout(width=1200, height=500)
    return style_plotly_figure(fig)