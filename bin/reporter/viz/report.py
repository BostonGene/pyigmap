import panel as pn
from .umi_tab import create_umi_tab
from .igblast_tab import create_igblast_summary_panel
from .final_tab import create_final_tab, create_final_per_chain_panel
from .chain_flow import plot_locus_distributions, build_chain_sankey


def create_report(final_clones_df, igblast_df,
                  umi_res=None, umi_sequences=None,
                  vidjil_df=None,
                  report_file_name='panel_report.html'):
    tab_stats = pn.Column(
        pn.pane.Markdown("## Vidjil+IgBlast Report"),
        pn.pane.Markdown("## Chain distributions"),
        plot_locus_distributions(vidjil_df, igblast_df),
        pn.pane.Markdown("## Chain distribution flow: Vidjil → IgBlast"),
        pn.pane.Plotly(build_chain_sankey(vidjil_df, igblast_df))
    ) if vidjil_df is not None else pn.Column(
        "Vidjil tool was not launched in this run.")

    tab_igblast = create_igblast_summary_panel(igblast_df)

    tab_umi = create_umi_tab(umi_res, umi_sequences) if umi_res is not None else pn.Column(
        "UMI information was not provided.")

    tab_final = create_final_tab(final_clones_df)

    tab_final_per_chain = create_final_per_chain_panel(final_clones_df)

    tabs = pn.Tabs(
        ("UMI Info", tab_umi),
        ('Vidjil + IgBlast', tab_stats),
        ('IgBlast Stats', tab_igblast),
        ('Final Clones', tab_final),
        ('Final Per-Chain', tab_final_per_chain)
    )
    report = pn.template.FastListTemplate(
        title="Pyigmap launch report",
        main=[tabs]
    )
    report.save(report_file_name)
