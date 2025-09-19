# MIBI Analysis Package

A Python package for preprocessing, analyzing, and visualizing MIBI (Multiplexed Ion Beam Imaging) data using MuVI factor analysis and other multicellular analysis approaches.

## Features

- **Preprocessing**: Normalize and prepare muon objects for factor analysis
- **MuVI/MOFA Integration**: Convenient wrappers for running factor analysis
- **Visualization**: Generate plots for factor scores, loadings, and confidence ellipses
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

The package requires Python 3.8+ and the following dependencies:

- numpy >= 1.20.0
- pandas >= 1.3.0
- scipy >= 1.7.0
- matplotlib >= 3.5.0
- seaborn >= 0.11.0
- plotnine >= 0.10.0
- mudata >= 0.2.0
- muon >= 0.1.0
- mofax >= 0.3.0
- scikit-learn >= 1.0.0

## Quick Start

```python
import mibi_analysis as ma
import mudata as mu

# Load your MuData object
mdata = mu.read_h5mu("path/to/your/data.h5mu")

# Preprocess data
ma.prepare_mudata_for_mofa(mdata, normalize=True, filter_features=True)

# Run MOFA analysis
model_path = ma.run_mofa_analysis(mdata, n_factors=10, outfile='mofa_model.h5ad')

# Calculate model performance
r2_results = ma.calculate_model_r2(model_path)
print(f"Overall R²: {r2_results['mean_r2']:.3f}")

# Extract factor scores with metadata
factor_scores = ma.extract_factor_scores(model_path, mdata=mdata, include_metadata=True)

# Test factor-clinical associations
associations = ma.test_factor_associations(
    factor_scores, 
    clinical_variables=['Stage', 'Sex'],
    test_type='auto'
)

# Visualize factor scores
plot = ma.plot_factor_scores(
    factor_scores, 
    'Factor 1', 'Factor 2', 
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

Wrapper functions for MOFA/MuVI analysis:

- `run_mofa_analysis()`: Run factor analysis with sensible defaults
- `calculate_model_r2()`: Compute variance explained by the model
- `extract_factor_scores()`: Get factor scores with metadata
- `extract_factor_loadings()`: Get factor loadings for interpretation

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