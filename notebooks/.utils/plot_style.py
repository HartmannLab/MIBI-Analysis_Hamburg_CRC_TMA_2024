"""Text sizes of figures for print.

A panel prints at its layout scale: the size it is placed at in the assembled figure divided by the size it is saved
at. Text sizes in code are the print sizes divided by that scale, so that tick labels print at 5 pt and titles at 6 pt.
A plot saved at its placed size has a scale of 1. From a notebook folder such as notebooks/lineage:

    import sys
    sys.path.insert(0, "../.utils")
    from plot_style import paper_theme, paper_text_size

    layout_scale = 0.43  # Fig 2a is placed at 43% of its saved size
    gp + paper_theme(layout_scale)

From the repository root, the path is "./notebooks/.utils". R notebooks source plot_style.R, which defines the same
functions.
"""
from plotnine import element_text, theme

TICK_PT = 5  # tick labels, legend labels, facet strips and annotations
TITLE_PT = 6  # axis titles, legend titles and plot titles


def paper_text_size(scale=1.0, pt=TICK_PT):
    """Size in code of text that prints at pt, for geom_text and geom_label"""
    return pt / scale


def paper_theme(scale=1.0, tick=TICK_PT, title=TITLE_PT):
    """plotnine theme with all text at its print size"""
    small, large = tick / scale, title / scale
    return theme(axis_text=element_text(size=small), axis_text_x=element_text(size=small),
                 axis_text_y=element_text(size=small), legend_text=element_text(size=small),
                 strip_text=element_text(size=small), strip_text_x=element_text(size=small),
                 strip_text_y=element_text(size=small), plot_subtitle=element_text(size=small),
                 axis_title=element_text(size=large), axis_title_x=element_text(size=large),
                 axis_title_y=element_text(size=large), legend_title=element_text(size=large),
                 plot_title=element_text(size=large))


def paper_rc(scale=1.0, tick=TICK_PT, title=TITLE_PT):
    """matplotlib rcParams with the same sizes, for plt.rc_context"""
    small, large = tick / scale, title / scale
    return {"font.size": small, "xtick.labelsize": small, "ytick.labelsize": small, "legend.fontsize": small,
            "legend.title_fontsize": large, "axes.labelsize": large, "axes.titlesize": large,
            "figure.titlesize": large}
