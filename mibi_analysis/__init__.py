"""
MIBI Analysis Package

A package for preprocessing, analyzing, and visualizing MIBI (Multiplexed Ion Beam Imaging) data
using MuVI factor analysis and other multicellular analysis approaches.

This package provides helper functions for:
- Preprocessing and normalizing muon objects
- Running MuVI factor analysis
- Visualizing factor scores and loadings
- Interpreting results with clinical associations
"""

from .preprocessing import (
    normalize_mudata_features,
    prepare_mudata_for_mofa,
    get_modality_stats
)

from .mofa import (
    run_mofa_analysis,
    calculate_model_r2,
    calculate_reconstruction_r2,
    extract_factor_scores,
    extract_factor_loadings
)

from .visualization import (
    get_cov_ellipse,
    plot_factor_scores,
    plot_factor_loadings,
    plot_r2_explained,
    plot_factor_comparison
)

from .interpretation import (
    test_factor_associations,
    annotate_factors,
    summarize_factor_by_group,
    identify_top_loadings,
    calculate_factor_stability
)

__version__ = "0.1.0"
__author__ = "MIBI Analysis Team"

__all__ = [
    # Preprocessing
    "normalize_mudata_features",
    "prepare_mudata_for_mofa",
    "get_modality_stats",
    # MOFA analysis
    "run_mofa_analysis", 
    "calculate_model_r2",
    "calculate_reconstruction_r2",
    "extract_factor_scores",
    "extract_factor_loadings",
    # Visualization
    "get_cov_ellipse",
    "plot_factor_scores",
    "plot_factor_loadings", 
    "plot_r2_explained",
    "plot_factor_comparison",
    # Interpretation
    "test_factor_associations",
    "annotate_factors",
    "summarize_factor_by_group",
    "identify_top_loadings",
    "calculate_factor_stability"
]