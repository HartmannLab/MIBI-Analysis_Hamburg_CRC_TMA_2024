# MIBI Analysis Package

A Python package for preprocessing, analyzing, and visualizing MIBI (Multiplexed Ion Beam Imaging) data using MuVI factor analysis and other multicellular analysis approaches, following the patterns from `notebooks/multicellular/MOFACell.ipynb`.

## Features

- **Preprocessing**: Normalize and prepare muon objects for factor analysis
- **MuVI/MOFA Integration**: Convenient wrappers for running factor analysis using muon and muvi
- **Liana Integration**: Uses liana utilities for factor score and loading extraction (as in MOFACell.ipynb)
- **Device Detection**: Automatic GPU/CPU device detection using muvi
- **Visualization**: Generate plots using plotnine for factor scores, loadings, and confidence ellipses
- **Statistical Analysis**: Test associations between factors and clinical variables
- **Interpretation**: Extract insights from factor analysis results

## Installation

### From Source (Development)

```bash
# Clone the repository
git clone https://github.com/HartmannLab/MIBI-Analysis_Hamburg_CRC_TMA_2024.git
cd MIBI-Analysis_Hamburg_CRC_TMA_2024

# Install in development mode
pip install -e .
```

### Dependencies

The package requires Python 3.8+ and the following dependencies (following MOFACell.ipynb requirements):

- numpy >= 1.20.0
- pandas >= 1.3.0
- scipy >= 1.7.0
- matplotlib >= 3.5.0
- seaborn >= 0.11.0
- plotnine >= 0.10.0
- mudata >= 0.2.0
- muon >= 0.1.0
- mofax >= 0.3.0
- **muvi >= 0.2.0** (for GPU detection and MuVI-specific functionality)
- **liana >= 1.0.0** (for factor score and loading extraction utilities)
- scikit-learn >= 1.0.0

## Quick Start

```python
import mibi_analysis as ma
import muon as mu
import liana as li
from plotnine import *

# Check device availability
device = ma.get_device()
print(f"Using device: {device}")

# Load your MuData object
features = mu.read_h5mu("path/to/your/data.h5mu")

# Preprocess data (center and scale as in MOFACell.ipynb)
ma.prepare_mudata_for_mofa(features, normalize=True, filter_features=True)

# Run MOFA analysis with same parameters as MOFACell.ipynb
model_path = ma.run_mofa_analysis(
    features, 
    n_factors=10,
    use_obs='union',
    convergence_mode='medium',
    scale_groups=False,
    scale_views=False,
    seed=1337
)

# Calculate model performance
r2_results = ma.calculate_model_r2(model_path)
print(f"Overall R²: {r2_results['mean_r2']:.3f}")

# Extract factor scores using liana utilities (as in MOFACell.ipynb)
factor_scores = ma.extract_factor_scores(
    features, 
    obsm_key='X_mofa', 
    obs_keys=['Stage'],
    use_liana=True
)

# Extract variable loadings using liana utilities  
variable_loadings = ma.extract_factor_loadings(
    features,
    varm_key='LFs',
    use_liana=True
)

# Test factor-clinical associations using Kruskal-Wallis
associations = ma.test_factor_associations(
    factor_scores, 
    clinical_variables=['Stage', 'Sex'],
    test_type='kruskal'
)

# Visualize factor scores with confidence ellipses (as in MOFACell.ipynb)
plot = ma.plot_factor_scores(
    factor_scores, 
    'Factor1', 'Factor2',  # Note: liana uses 'Factor1' not 'Factor 1'
    color_by='Stage', 
    add_ellipses=True
)
```

## Modules

### `mibi_analysis.preprocessing`

Functions for normalizing and preparing MuData objects:

- `normalize_mudata_features()`: Center and scale features across modalities
- `prepare_mudata_for_mofa()`: Complete preprocessing pipeline for MOFA
- `get_modality_stats()`: Generate summary statistics for modalities

### `mibi_analysis.mofa`

Wrapper functions for MOFA/MuVI analysis following MOFACell.ipynb patterns:

- `get_device()`: Detect GPU/CPU device using muvi (as in MOFACell.ipynb)
- `run_mofa_analysis()`: Run factor analysis with same parameters as MOFACell.ipynb
- `calculate_model_r2()`: Compute variance explained by the model
- `extract_factor_scores()`: Get factor scores using liana utilities (preferred) or model file
- `extract_factor_loadings()`: Get factor loadings using liana utilities (preferred) or model file

### `mibi_analysis.visualization`

Plotting functions for results visualization:

- `get_cov_ellipse()`: Generate confidence ellipse points
- `plot_factor_scores()`: Scatter plots of factor scores
- `plot_factor_loadings()`: Heatmaps of factor loadings
- `plot_r2_explained()`: Variance explained bar plots
- `plot_factor_comparison()`: Compare factors across groups

### `mibi_analysis.interpretation`

Statistical analysis and interpretation functions:

- `test_factor_associations()`: Test factor-clinical variable associations
- `annotate_factors()`: Add clinical metadata to factor scores
- `summarize_factor_by_group()`: Generate group summary statistics
- `identify_top_loadings()`: Find most important features per factor

## Examples

See the `examples/` directory for detailed usage examples:

- `examples/mibi_analysis_example.ipynb`: Complete workflow demonstration

## Development

### Running Tests

```bash
# Install development dependencies
pip install -e ".[dev]"

# Run tests
pytest tests/
```

### Code Formatting

```bash
# Format code
black mibi_analysis/

# Check code style
flake8 mibi_analysis/
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this package in your research, please cite:

```
[Citation information will be added upon publication]
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Support

For questions or issues, please open an issue on GitHub or contact the development team.