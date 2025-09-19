"""
Visualization functions for MIBI factor analysis results.

This module provides plotting functions for factor scores, loadings, and 
confidence ellipses as used in the MOFACell analysis.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from plotnine import *
from typing import Optional, Dict, List, Tuple, Union
import warnings


def get_cov_ellipse(
    cov: np.ndarray, 
    center: np.ndarray, 
    nstd: float = 2.0,
    n_points: int = 100
) -> pd.DataFrame:
    """
    Generate points for a covariance ellipse.
    
    This function creates an ellipse representing the covariance matrix centered
    at a given point, as used in the MOFACell analysis for confidence ellipses.
    
    Parameters
    ----------
    cov : np.ndarray
        2x2 covariance matrix
    center : np.ndarray
        Center point of the ellipse (x, y)
    nstd : float, default 2.0
        Number of standard deviations for the ellipse radius
    n_points : int, default 100
        Number of points to generate for the ellipse
        
    Returns
    -------
    pd.DataFrame
        DataFrame with 'x' and 'y' columns containing ellipse points
        
    Examples
    --------
    >>> import mibi_analysis as ma
    >>> import numpy as np
    >>> # Create example covariance matrix and center
    >>> cov = np.array([[1, 0.5], [0.5, 2]])
    >>> center = np.array([0, 0])
    >>> ellipse_points = ma.get_cov_ellipse(cov, center, nstd=2)
    """
    # Find and sort eigenvalues and eigenvectors
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = eigvals.argsort()[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]

    # Calculate the angle for the ellipse
    angle = np.degrees(np.arctan2(*eigvecs[:, 0][::-1]))

    # Width and height are 2 * nstd * sqrt(eigenvalue)
    width, height = 2 * nstd * np.sqrt(eigvals)

    # Generate ellipse points
    theta = np.linspace(0, 2 * np.pi, n_points)
    ellipse_x = center[0] + (width/2) * np.cos(theta) * np.cos(np.radians(angle)) - (height/2) * np.sin(theta) * np.sin(np.radians(angle))
    ellipse_y = center[1] + (width/2) * np.cos(theta) * np.sin(np.radians(angle)) + (height/2) * np.sin(theta) * np.cos(np.radians(angle))

    return pd.DataFrame({'x': ellipse_x, 'y': ellipse_y})


def plot_factor_scores(
    factor_scores: pd.DataFrame,
    x_factor: str,
    y_factor: str,
    color_by: Optional[str] = None,
    add_ellipses: bool = False,
    ellipse_nstd: float = 2.0,
    figsize: Tuple[int, int] = (8, 6),
    title: Optional[str] = None,
    **kwargs
) -> Union[object, plt.Figure]:
    """
    Plot factor scores with optional confidence ellipses.
    
    Parameters
    ----------
    factor_scores : pd.DataFrame
        DataFrame with factor scores and metadata
    x_factor : str
        Column name for x-axis factor
    y_factor : str
        Column name for y-axis factor  
    color_by : str, optional
        Column name to color points by
    add_ellipses : bool, default False
        Whether to add confidence ellipses for each group
    ellipse_nstd : float, default 2.0
        Number of standard deviations for ellipses
    figsize : tuple, default (8, 6)
        Figure size (width, height)
    title : str, optional
        Plot title
    **kwargs
        Additional arguments for scatter plot styling
        
    Returns
    -------
    plotnine plot object or matplotlib Figure
        Plot object
        
    Examples
    --------
    >>> import mibi_analysis as ma
    >>> # Basic factor score plot
    >>> plot = ma.plot_factor_scores(scores, 'Factor 1', 'Factor 2')
    >>> 
    >>> # Colored by clinical variable with ellipses
    >>> plot = ma.plot_factor_scores(
    ...     scores, 'Factor 1', 'Factor 2', 
    ...     color_by='Stage', add_ellipses=True
    ... )
    """
    # Check if required columns exist
    missing_cols = []
    for col in [x_factor, y_factor]:
        if col not in factor_scores.columns:
            missing_cols.append(col)
    if color_by and color_by not in factor_scores.columns:
        missing_cols.append(color_by)
        
    if missing_cols:
        raise ValueError(f"Missing columns in factor_scores: {missing_cols}")
    
    # Create base plot
    if color_by:
        # Filter out missing values for color variable
        plot_data = factor_scores.dropna(subset=[color_by])
        
        p = (ggplot(plot_data, aes(x=x_factor, y=y_factor, color=color_by)) +
             geom_point(size=2, alpha=0.7, **kwargs) +
             theme_classic())
    else:
        p = (ggplot(factor_scores, aes(x=x_factor, y=y_factor)) +
             geom_point(size=2, alpha=0.7, **kwargs) +
             theme_classic())
    
    # Add confidence ellipses if requested
    if add_ellipses and color_by:
        ellipses_data = pd.DataFrame()
        
        for group in plot_data[color_by].unique():
            if pd.isna(group):
                continue
                
            group_data = plot_data[plot_data[color_by] == group]
            
            if len(group_data) < 3:  # Need at least 3 points for covariance
                warnings.warn(f"Skipping ellipse for group '{group}' (insufficient data points)")
                continue
                
            x = group_data[x_factor]
            y = group_data[y_factor]
            
            # Calculate mean and covariance
            mean = np.array([x.mean(), y.mean()])
            cov = np.cov(x, y)
            
            # Generate ellipse points
            ellipse_df = get_cov_ellipse(cov, mean, nstd=ellipse_nstd)
            ellipse_df[color_by] = group
            ellipses_data = pd.concat([ellipses_data, ellipse_df], ignore_index=True)
        
        if not ellipses_data.empty:
            p = p + geom_path(
                data=ellipses_data,
                mapping=aes(x='x', y='y', color=color_by, group=color_by),
                linetype='dashed', size=1
            )
    
    # Add title if provided
    if title:
        p = p + ggtitle(title)
    
    # Add axis labels
    p = p + xlab(x_factor) + ylab(y_factor)
    
    return p


def plot_factor_loadings(
    loadings: Dict[str, pd.DataFrame],
    factors_to_plot: Optional[List[str]] = None,
    views_to_plot: Optional[List[str]] = None,
    figsize: Tuple[int, int] = (12, 8),
    cmap: str = 'RdBu_r',
    center: float = 0.0,
    title: Optional[str] = None
) -> plt.Figure:
    """
    Plot factor loadings heatmap for multiple views.
    
    Parameters
    ----------
    loadings : Dict[str, pd.DataFrame]
        Dictionary mapping view names to loading DataFrames
    factors_to_plot : List[str], optional
        List of factors to include in plot. If None, plots all factors
    views_to_plot : List[str], optional
        List of views to include in plot. If None, plots all views
    figsize : tuple, default (12, 8)
        Figure size (width, height)
    cmap : str, default 'RdBu_r'
        Colormap for heatmap
    center : float, default 0.0
        Value to center colormap at
    title : str, optional
        Plot title
        
    Returns
    -------
    plt.Figure
        Matplotlib figure object
        
    Examples
    --------
    >>> import mibi_analysis as ma
    >>> loadings = ma.extract_factor_loadings('mofa_model.h5ad')
    >>> fig = ma.plot_factor_loadings(loadings, factors_to_plot=['Factor 1', 'Factor 2'])
    """
    views_to_use = views_to_plot if views_to_plot else list(loadings.keys())
    
    # Combine loadings from all views
    combined_loadings = []
    view_labels = []
    
    for view in views_to_use:
        if view not in loadings:
            warnings.warn(f"View '{view}' not found in loadings. Skipping.")
            continue
            
        view_loadings = loadings[view]
        
        # Select factors to plot
        if factors_to_plot:
            available_factors = [f for f in factors_to_plot if f in view_loadings.columns]
            if not available_factors:
                warnings.warn(f"No requested factors found in view '{view}'. Skipping.")
                continue
            view_loadings = view_loadings[available_factors]
        
        combined_loadings.append(view_loadings.values)
        view_labels.extend([f"{view}_{i}" for i in range(len(view_loadings))])
    
    if not combined_loadings:
        raise ValueError("No valid loadings data found for plotting")
    
    # Create combined matrix
    combined_matrix = np.vstack(combined_loadings)
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot heatmap
    sns.heatmap(
        combined_matrix,
        ax=ax,
        cmap=cmap,
        center=center,
        yticklabels=view_labels,
        xticklabels=view_loadings.columns if factors_to_plot else None,
        cbar_kws={'label': 'Loading'}
    )
    
    if title:
        ax.set_title(title)
    
    ax.set_xlabel('Factors')
    ax.set_ylabel('Features (by View)')
    
    plt.tight_layout()
    return fig


def plot_r2_explained(
    r2_results: Dict[str, float],
    figsize: Tuple[int, int] = (10, 6),
    title: Optional[str] = None
) -> plt.Figure:
    """
    Plot R² (variance explained) for different modalities.
    
    Parameters
    ----------
    r2_results : Dict[str, float]
        Dictionary with R² values from calculate_model_r2()
    figsize : tuple, default (10, 6)
        Figure size (width, height)
    title : str, optional
        Plot title
        
    Returns
    -------
    plt.Figure
        Matplotlib figure object
        
    Examples
    --------
    >>> import mibi_analysis as ma
    >>> r2_results = ma.calculate_model_r2('mofa_model.h5ad')
    >>> fig = ma.plot_r2_explained(r2_results)
    """
    # Extract modality-specific R² values
    modality_r2 = {}
    for key, value in r2_results.items():
        if key.endswith('_r2') and key != 'mean_r2':
            modality_name = key.replace('_r2', '')
            modality_r2[modality_name] = value
    
    if not modality_r2:
        raise ValueError("No modality-specific R² values found in results")
    
    # Create DataFrame for plotting
    plot_data = pd.DataFrame(list(modality_r2.items()), columns=['Modality', 'R2'])
    plot_data = plot_data.sort_values('R2', ascending=True)
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create horizontal bar plot
    bars = ax.barh(plot_data['Modality'], plot_data['R2'])
    
    # Add value labels on bars
    for i, (bar, r2) in enumerate(zip(bars, plot_data['R2'])):
        ax.text(bar.get_width() + 0.005, bar.get_y() + bar.get_height()/2, 
                f'{r2:.3f}', ha='left', va='center')
    
    # Add mean R² line if available
    if 'mean_r2' in r2_results:
        ax.axvline(r2_results['mean_r2'], color='red', linestyle='--', 
                   label=f"Mean R² = {r2_results['mean_r2']:.3f}")
        ax.legend()
    
    ax.set_xlabel('R² (Variance Explained)')
    ax.set_ylabel('Modality')
    
    if title:
        ax.set_title(title)
    else:
        ax.set_title('Variance Explained by Modality')
    
    # Set x-axis limits
    ax.set_xlim(0, max(plot_data['R2']) * 1.1)
    
    plt.tight_layout()
    return fig


def plot_factor_comparison(
    factor_scores: pd.DataFrame,
    factor1: str,
    factor2: str,
    group_by: str,
    plot_type: str = 'violin',
    figsize: Tuple[int, int] = (12, 6)
) -> plt.Figure:
    """
    Create comparison plots for factors across groups.
    
    Parameters
    ----------
    factor_scores : pd.DataFrame
        DataFrame with factor scores and metadata
    factor1 : str
        First factor to plot
    factor2 : str  
        Second factor to plot
    group_by : str
        Column to group by
    plot_type : str, default 'violin'
        Type of plot: 'violin', 'box', or 'strip'
    figsize : tuple, default (12, 6)
        Figure size
        
    Returns
    -------
    plt.Figure
        Matplotlib figure with subplots
        
    Examples
    --------
    >>> import mibi_analysis as ma
    >>> fig = ma.plot_factor_comparison(scores, 'Factor 1', 'Factor 2', 'Stage')
    """
    # Filter out missing values
    plot_data = factor_scores.dropna(subset=[group_by])
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    
    # Plot first factor
    if plot_type == 'violin':
        sns.violinplot(data=plot_data, x=group_by, y=factor1, ax=ax1)
    elif plot_type == 'box':
        sns.boxplot(data=plot_data, x=group_by, y=factor1, ax=ax1)
    elif plot_type == 'strip':
        sns.stripplot(data=plot_data, x=group_by, y=factor1, ax=ax1)
    else:
        raise ValueError(f"Unknown plot_type: {plot_type}")
    
    ax1.set_title(f'{factor1} by {group_by}')
    ax1.tick_params(axis='x', rotation=45)
    
    # Plot second factor
    if plot_type == 'violin':
        sns.violinplot(data=plot_data, x=group_by, y=factor2, ax=ax2)
    elif plot_type == 'box':
        sns.boxplot(data=plot_data, x=group_by, y=factor2, ax=ax2)
    elif plot_type == 'strip':
        sns.stripplot(data=plot_data, x=group_by, y=factor2, ax=ax2)
    
    ax2.set_title(f'{factor2} by {group_by}')
    ax2.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    return fig