import tempfile

import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import pandas as pd


def read_fastq_file_chunk(fastq_file: str) -> list[int]:
    """Reads a FASTQ file chunk and returns a list of reads."""
    consensus_group_size_per_read = []
    with open(fastq_file, 'r') as f:
        while header := f.readline().strip():  # read header string
            if not header:
                break
            _ = f.readline()
            _ = f.readline()
            _ = f.readline()
            consensus_group_size = len(header.split('\t')[1].split(';')) - 1  # @0  145607;265853;563279;
            consensus_group_size_per_read.append(consensus_group_size)
    return consensus_group_size_per_read


def draw_umi_coverage_plot(reads_group_size_list: list[int], output_plot: str):
    arr = np.array(reads_group_size_list)

    df = pd.DataFrame(arr, columns=['number'])

    frequency = df['number'].value_counts(normalize=True)

    frequency_df = frequency.reset_index()
    frequency_df.columns = ['number', 'percentage']
    frequency_df['percentage'] = frequency_df['percentage'] * frequency_df['number']
    frequency_df['number'] = np.log2(frequency_df['number'])

    fig = px.histogram(frequency_df, x='number', y='percentage',
                       title="UMI coverage plot",
                       labels={"number": "reads count", "percentage": "frequency"}, nbins=100, histnorm='probability')

    fig.update_layout(xaxis=dict(title="reads count"), yaxis=dict(title="frequency (weighted by the number of reads)"),
                      title_x=0.5)

    fig.write_html(output_plot)


def draw_barcode_rank_plot(reads_group_size_list: list[int], output_plot_path: str):
    umi_counts_sorted = sorted(reads_group_size_list, reverse=True)

    # Создаем массив рангов
    ranks = np.arange(1, len(umi_counts_sorted) + 1)

    # Создаем график
    fig = go.Figure()

    # Добавляем линию графика
    fig.add_trace(go.Scatter(
        x=ranks,
        y=umi_counts_sorted,
        mode='lines+markers',
        name='UMI counts'
    ))

    fig.update_layout(
        title='Barcode Rank Plot',
        xaxis_title='Rank',
        yaxis_title='UMI counts',
        xaxis=dict(type='linear'),
        yaxis=dict(type='log')
    )

    fig.write_html(output_plot_path)
    fig.write_image("barcode_rank_plot.png")


if __name__ == '__main__':
    group_size_list = read_fastq_file_chunk('~/cR1.fastq')
    coverage_plot = tempfile.NamedTemporaryFile(suffix=".html").name
    draw_umi_coverage_plot(group_size_list, coverage_plot)
    print(coverage_plot)

    barcode_rank_plot = tempfile.NamedTemporaryFile(suffix=".html").name
    draw_barcode_rank_plot(group_size_list, barcode_rank_plot)
    print(barcode_rank_plot)
