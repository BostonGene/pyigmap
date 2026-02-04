from collections import Counter

import numpy as np
import pandas as pd
import panel as pn
import plotly.express as px
import plotly.graph_objects as go
from bokeh.models import HoverTool, NumeralTickFormatter
from bokeh.plotting import figure


def compute_logo_matrix(sequences):
    seq_len = len(sequences[0])
    counts = {i: Counter() for i in range(seq_len)}
    for seq in sequences:
        for i, char in enumerate(seq):
            counts[i][char] += 1
    freq_df = pd.DataFrame.from_dict(
        {i: {char: count / sum(counts[i].values()) for char, count in counts[i].items()} for i in counts},
        orient='index',
    ).fillna(0)
    return freq_df


def plotly_logo(sequences):
    freq_df = compute_logo_matrix(sequences)
    fig = go.Figure()
    y_offset = np.zeros(len(freq_df))
    for letter in sorted(freq_df.columns, key=lambda x: -freq_df[x].sum()):
        heights = freq_df[letter].to_numpy()
        hovertext = [f'{letter}: {100 * h:.1f}%' for h in heights]
        writtentext = [letter if h > 0.05 else '' for h in heights]
        fig.add_trace(
            go.Bar(
                x=list(freq_df.index),
                y=heights,
                base=y_offset,
                hovertext=hovertext,
                hoverinfo='text',
                text=writtentext,
                name=letter,
                marker={'line': {'width': 0}},
                width=0.9,
            )
        )
        y_offset += heights
    fig.update_layout(
        barmode='stack',
        title='Sequence Logo (Interactive)',
        xaxis_title='Position',
        yaxis_title='Frequency',
        xaxis={'tickmode': 'linear'},
        yaxis={'range': [0, 1]},
        showlegend=False,
        height=400,
    )
    return fig


def bokeh_histogram(res, plot_by='read_num_pre_dedup', bins=20, x_max=None, y_range=None):
    data = res[plot_by]
    hist, edges = np.histogram(data, bins=np.logspace(np.log10(1), np.log10(data.max() + 1), bins))
    p = figure(
        title=f'UMI Usage ({plot_by})',
        x_axis_type='log',
        x_axis_label='Read Count',
        y_axis_label='Number of UMIs',
        y_range=y_range,  # type: ignore[arg-type]
        height=500,
        width=600,
        tools='pan,wheel_zoom,box_zoom,reset,save',  # type: ignore[arg-type]
    )
    p.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:], fill_color='skyblue', line_color='black')

    p.x_range.end = x_max  # type: ignore[union-attr]
    p.yaxis.formatter = NumeralTickFormatter(format='0')  # type: ignore[index]
    p.add_tools(HoverTool(tooltips=[('Bin', '@left - @right'), ('Count', '@top')], mode='mouse'))
    return p


def plot_umi_scatter(res):
    res = res.copy()
    res = res[res['read_num_pre_dedup'] > 0]
    res['delta'] = res['read_num_post_dedup'] / res['read_num_pre_dedup']
    top10 = res.sort_values('delta', ascending=False).head(int(len(res) * 0.1))

    max_val = max(top10['read_num_pre_dedup'].max(), top10['read_num_post_dedup'].max()) * 1.01

    fig = px.scatter(
        top10,
        x='read_num_pre_dedup',
        y='read_num_post_dedup',
        hover_name='umi',
        log_x=True,
        log_y=True,
        title='Top 10% UMI Usage Increase After Deduplication',
        labels={'read_num_pre_dedup': 'Pre-dedup Read Count', 'read_num_post_dedup': 'Post-dedup Read Count'},
    )

    fig.add_trace(
        go.Scatter(
            x=[0, max_val], y=[0, max_val], mode='lines', line={'dash': 'dash', 'color': 'gray'}, showlegend=False
        )
    )

    fig.update_layout(
        height=600,
        width=600,
        # xaxis=dict(range=[np.log10(1), np.log10(max_val)], type="log"),
        # yaxis=dict(range=[np.log10(1), np.log10(max_val)], type="log")
    )
    fig.update_traces(marker={'size': 6, 'color': 'darkorange'})
    return fig


def bokeh_weighted_usage(res, col='read_num_pre_dedup', bins=20, x_max=None, y_range=None):
    df = res.copy()
    total_reads = df[col].sum()
    df['weight'] = df[col] / total_reads
    values = df[col].values
    weights = df['weight'].values
    hist, edges = np.histogram(values, bins=np.logspace(np.log10(1), np.log10(values.max() + 1), bins), weights=weights)
    p = figure(
        title=f'Weighted UMI Usage ({col})',
        x_axis_type='log',
        x_axis_label='Read Count',
        y_axis_label='Weighted UMI Frequency',
        #
        y_range=y_range,  # type: ignore[arg-type]
        height=500,
        width=600,
        tools='pan,wheel_zoom,box_zoom,reset,save',  # type: ignore[arg-type]
    )
    p.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:], fill_color='skyblue', line_color='black')

    # p.x_range.start = x_min
    p.x_range.end = x_max  # type: ignore[union-attr]
    p.yaxis.formatter = NumeralTickFormatter(format='0.00')  # type: ignore[index]
    p.add_tools(HoverTool(tooltips=[('Bin', '@left - @right'), ('Weighted Freq', '@top')], mode='mouse'))
    return p


def plot_barcode_rank(res):
    df = res.copy()
    df = df[df['read_num_post_dedup'] > 0]
    sorted_reads = df['read_num_post_dedup'].sort_values(ascending=False).reset_index(drop=True)
    cumulative_reads = sorted_reads.cumsum()

    if len(sorted_reads) > 2000:
        log_space = np.geomspace(1, len(sorted_reads), num=2000).astype(int) - 1
        log_space = np.unique(np.clip(log_space, 0, len(sorted_reads) - 1))
        sorted_reads = sorted_reads.iloc[log_space].reset_index(drop=True)
        cumulative_reads = cumulative_reads.iloc[log_space].reset_index(drop=True)

    fig = go.Figure()

    fig.add_trace(
        go.Scatter(x=cumulative_reads, y=sorted_reads, mode='lines', line={'color': 'teal'}, name='read_counts')
    )

    fig.update_layout(
        title='UMI Rank Plot (Cumulative Reads vs UMI Read Count)',
        xaxis_title='Cumulative Reads',
        yaxis_title='Read Count per UMI',
        xaxis_type='log',
        yaxis_type='log',
        height=500,
        width=600,
    )

    return fig


def create_report(res, umi_sequences, report_file_name='panel_report.html'):
    x_max = max(res['read_num_pre_dedup'].max(), res['read_num_post_dedup'].max())
    y_max_hist = max(
        np.histogram(
            res['read_num_pre_dedup'], bins=np.logspace(np.log10(1), np.log10(res['read_num_pre_dedup'].max() + 1), 20)
        )[0].max(),
        np.histogram(
            res['read_num_post_dedup'],
            bins=np.logspace(np.log10(1), np.log10(res['read_num_post_dedup'].max() + 1), 20),
        )[0].max(),
    )
    y_max_weight = max(
        res['read_num_pre_dedup'].value_counts(normalize=True).max(),
        res['read_num_post_dedup'].value_counts(normalize=True).max(),
    )
    y_range_hist = (0, y_max_hist * 1.01)
    y_range_weight = (0, y_max_weight * 1.01)

    tab_logo = pn.Column(
        pn.pane.Markdown('## Sequence Logo'),
        pn.pane.Plotly(plotly_logo([x[:12].decode('ascii') for x in umi_sequences]), config={'responsive': True}),
    )
    tab_umi_usage = pn.Column(
        pn.pane.Markdown('## UMI Usage Pre- and Post-Deduplication'),
        pn.Row(
            pn.pane.Bokeh(bokeh_histogram(res, plot_by='read_num_pre_dedup', x_max=x_max, y_range=y_range_hist)),
            pn.pane.Bokeh(bokeh_histogram(res, plot_by='read_num_post_dedup', x_max=x_max, y_range=y_range_hist)),
        ),
    )
    tab_scatter = pn.Column(
        pn.pane.Markdown('## Top 10% UMI Usage Increase (Scatterplot)'),
        pn.pane.Plotly(plot_umi_scatter(res), config={'responsive': True}),
    )
    tab_weighted = pn.Column(
        pn.pane.Markdown('## Weighted UMI Usage by Read Count'),
        pn.Row(
            pn.pane.Bokeh(bokeh_weighted_usage(res, col='read_num_pre_dedup', x_max=x_max, y_range=y_range_weight)),
            pn.pane.Bokeh(bokeh_weighted_usage(res, col='read_num_post_dedup', x_max=x_max, y_range=y_range_weight)),
        ),
    )
    tab_rank = pn.Column(
        pn.pane.Markdown('## UMI Barcode Rank Plot'),
        pn.pane.Plotly(plot_barcode_rank(res), config={'responsive': True}),
    )

    tabs = pn.Tabs(
        ('Sequence Logo', tab_logo),
        ('UMI Usage', tab_umi_usage),
        ('Weighted UMI Usage', tab_weighted),
        ('UMI Usage Changes', tab_scatter),
        ('Read Rank Plot', tab_rank),
    )
    report = pn.template.FastListTemplate(title='UMI Analysis Report', main=[tabs])
    report.save(report_file_name)
