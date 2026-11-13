"""Names shown in text outputs and figures.

Data files keep the names of the annotation, for example `Cancer_cell` and `CAF`. Notebooks use these helpers when
they print or plot them. The import path follows the working directory that the data paths of the notebook assume.
From the notebook folder, as with data_dir = "../../data/":

    import sys
    sys.path.insert(0, "../.utils")
    from plot_labels import display_name, display_types, display_text

From the repository root, as with data_dir = "./data/":

    sys.path.insert(0, "./notebooks/.utils")

R notebooks source plot_labels.R, which defines the same functions.
"""
import re

import pandas as pd

# Cell types named differently in figures
display_names = {"CAF": "Fibroblast", "Cancer_cell": "Epithelial cell"}

# Cell types of the consensus annotation (lineage/annotation_consensus.ipynb) and of the groups of
# multicellular/PrepareMuData.ipynb, as found in feature names
cell_types = ["APC", "B_cell", "CAF", "CD163_Macrophage", "CD4_Tcell", "CD68_Macrophage", "CD8_Tcell", "Cancer_cell",
              "Endothelial_cell", "Monocyte", "NK_cell", "Neutrophil", "Other_immune_cell", "T_reg_cell", "Unclear",
              "CD4_lymphocyte", "Cytotoxic_lymphocyte", "Epithelial_cell", "Fibroblast", "Macrophage"]
# Longest names first, so that CD68_Macrophage is not read as CD68_ followed by Macrophage
cell_type_pattern = re.compile("|".join(re.escape(t) for t in sorted(cell_types, key=len, reverse=True)))


def display_name(cell_type):
    return display_names.get(cell_type, cell_type).replace("_", " ")


def display_types(cell_types, order):
    """Display names of cell types, as a categorical in the given order"""
    return pd.Categorical([display_name(t) for t in cell_types], categories=[display_name(t) for t in order])


def display_text(text):
    """Text with the display name of each cell type it contains, for example in cell type pairs such as CAF_Cancer_cell.

    MISTy feature names describe the features of the model and keep the annotation names.
    """
    return cell_type_pattern.sub(lambda match: display_name(match.group()), text)
