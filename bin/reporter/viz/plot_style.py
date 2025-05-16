import plotly.graph_objects as go

def style_plotly_figure(fig: go.Figure, line_color='black', line_width=0.5) -> go.Figure:
    fig.update_layout(
        paper_bgcolor='rgba(0,0,0,0)',
        # plot_bgcolor='rgba(0,0,0,0)',
    )
    for trace in fig.data:
        if hasattr(trace, "marker") and isinstance(trace.marker, dict) and "line" in trace.marker:
            trace.marker.line.update(color=line_color, width=line_width)
        elif hasattr(trace, "marker"):
            trace.marker.line = dict(color=line_color, width=line_width)
    return fig
