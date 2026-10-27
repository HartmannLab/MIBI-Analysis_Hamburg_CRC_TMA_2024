"""Names shown in text outputs and figures.

Data files keep the names of the annotation, for example `Cancer_cell` and `CAF`. Notebooks use these helpers when
they print or plot them. From a notebook folder such as notebooks/sampling:

    import sys
    sys.path.insert(0, "../.utils")
    from plot_labels import display_name, display_types
"""
import pandas as pd

# Cell types named differently in figures
display_names = {"CAF": "Fibroblast", "Cancer_cell": "Epithelial cell"}


def display_name(cell_type):
    return display_names.get(cell_type, cell_type).replace("_", " ")


def display_types(cell_types, order):
    """Display names of cell types, as a categorical in the given order"""
    return pd.Categorical([display_name(t) for t in cell_types], categories=[display_name(t) for t in order])
