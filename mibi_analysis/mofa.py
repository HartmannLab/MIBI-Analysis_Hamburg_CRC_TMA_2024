"""
MOFA/MuVI wrapper functions for factor analysis of MIBI data.

This module provides convenient wrapper functions for running MOFA factor analysis
on MuData objects and extracting results.
"""

import numpy as np
import pandas as pd
import mudata as mu
import muon
import mofax as mofa
from typing import Optional, Dict, Any, Union
import warnings
import os


def run_mofa_analysis(
    mdata: mu.MuData,
    n_factors: int = 10,
    outfile: Optional[str] = None,
    convergence_mode: str = 'medium',
    scale_groups: bool = False,
    scale_views: bool = False,
    seed: int = 1337,
    use_var: Optional[str] = None,
    use_obs: str = 'union',
    **kwargs
) -> str:
    """
    Run MOFA factor analysis on a MuData object.
    
    This function provides a convenient wrapper around muon.tl.mofa with 
    sensible defaults for MIBI analysis.
    
    Parameters
    ----------
    mdata : mu.MuData
        Input MuData object with multiple modalities
    n_factors : int, default 10
        Number of factors to infer
    outfile : str, optional
        Path to save the MOFA model. If None, uses temporary file
    convergence_mode : str, default 'medium'
        Convergence criteria: 'fast', 'medium', or 'slow'
    scale_groups : bool, default False
        Whether to scale groups (samples)
    scale_views : bool, default False
        Whether to scale views (modalities)
    seed : int, default 1337
        Random seed for reproducibility
    use_var : str, optional
        Variable selection criteria
    use_obs : str, default 'union'
        Observation selection criteria
    **kwargs
        Additional arguments passed to muon.tl.mofa
        
    Returns
    -------
    str
        Path to the saved MOFA model file
        
    Examples
    --------
    >>> import mibi_analysis as ma
    >>> # Run MOFA with default parameters
    >>> model_path = ma.run_mofa_analysis(mdata, outfile='mofa_model.h5ad')
    >>> 
    >>> # Run with custom number of factors
    >>> model_path = ma.run_mofa_analysis(mdata, n_factors=15, convergence_mode='slow')
    """
    # Generate output filename if not provided
    if outfile is None:
        outfile = 'mofa_model.h5ad'
    
    # Run MOFA analysis
    muon.tl.mofa(
        mdata,
        n_factors=n_factors,
        use_obs=use_obs,
        convergence_mode=convergence_mode,
        scale_groups=scale_groups,
        scale_views=scale_views,
        seed=seed,
        outfile=outfile,
        use_var=use_var,
        **kwargs
    )
    
    return outfile


def calculate_model_r2(model_path: str, close_model: bool = True) -> Dict[str, float]:
    """
    Calculate R² (variance explained) for a MOFA model.
    
    Parameters
    ----------
    model_path : str
        Path to the saved MOFA model file
    close_model : bool, default True
        Whether to close the model after calculation
        
    Returns
    -------
    Dict[str, float]
        Dictionary with R² values including 'mean_r2' for overall model performance
        and individual modality R² values
        
    Examples
    --------
    >>> import mibi_analysis as ma
    >>> r2_results = ma.calculate_model_r2('mofa_model.h5ad')
    >>> print(f"Overall R²: {r2_results['mean_r2']:.3f}")
    """
    # Load the model
    model = mofa.mofa_model(model_path)
    
    try:
        # Calculate variance explained
        model_r2 = model.calculate_variance_explained().sort_values("R2", ascending=False)
        
        # Extract overall R² (macro average)
        mean_r2 = model_r2.R2.mean() / 100  # Convert from percentage
        
        # Create results dictionary
        results = {'mean_r2': mean_r2}
        
        # Add individual modality R² values
        for _, row in model_r2.iterrows():
            results[f"{row['view']}_r2"] = row['R2'] / 100
            
    finally:
        if close_model:
            model.close()
    
    return results


def calculate_reconstruction_r2(
    model_path: str, 
    mdata: mu.MuData, 
    close_model: bool = True
) -> Dict[str, float]:
    """
    Calculate reconstruction R² for each modality by comparing original and reconstructed data.
    
    This provides an alternative R² calculation based on correlation between
    original features and their reconstruction from factor loadings.
    
    Parameters
    ----------
    model_path : str
        Path to the saved MOFA model file
    mdata : mu.MuData
        Original MuData object used for training
    close_model : bool, default True
        Whether to close the model after calculation
        
    Returns
    -------
    Dict[str, float]
        Dictionary with reconstruction R² for each modality and macro average
        
    Examples
    --------
    >>> import mibi_analysis as ma
    >>> recon_r2 = ma.calculate_reconstruction_r2('mofa_model.h5ad', mdata)
    >>> print(f"Macro reconstruction R²: {recon_r2['macro_r2']:.3f}")
    """
    model = mofa.mofa_model(model_path)
    
    try:
        # Get factor scores
        Z = np.array(model.expectations["Z"]["group1"])
        
        r2_values = []
        results = {}
        
        for view_name in model.views:
            if view_name not in mdata.mod:
                warnings.warn(f"View '{view_name}' not found in MuData object. Skipping.")
                continue
                
            # Get factor loadings for this view
            W = model.expectations["W"][view_name]
            
            # Reconstruct data
            reconstructed = Z.T.dot(W)
            
            # Get original data
            original = mdata[view_name].X
            
            # Calculate correlation (R²)
            correlation_df = pd.DataFrame({
                'x': original.flatten(), 
                'y': reconstructed.flatten()
            }).corr()
            
            r2 = correlation_df.iloc[0, 1]
            r2_values.append(r2)
            results[f"{view_name}_reconstruction_r2"] = r2
            
        # Calculate macro average
        results['macro_reconstruction_r2'] = np.mean(r2_values) if r2_values else 0.0
        
    finally:
        if close_model:
            model.close()
            
    return results


def extract_factor_scores(
    model_path: str, 
    mdata: Optional[mu.MuData] = None,
    include_metadata: bool = True,
    close_model: bool = True
) -> pd.DataFrame:
    """
    Extract factor scores from a MOFA model and optionally merge with metadata.
    
    Parameters
    ----------
    model_path : str
        Path to the saved MOFA model file
    mdata : mu.MuData, optional
        Original MuData object to extract metadata from
    include_metadata : bool, default True
        Whether to include metadata from mdata.obs if mdata is provided
    close_model : bool, default True
        Whether to close the model after extraction
        
    Returns
    -------
    pd.DataFrame
        DataFrame with factor scores and optionally metadata
        
    Examples
    --------
    >>> import mibi_analysis as ma
    >>> # Extract factor scores only
    >>> scores = ma.extract_factor_scores('mofa_model.h5ad')
    >>> 
    >>> # Extract with metadata
    >>> scores_with_meta = ma.extract_factor_scores('mofa_model.h5ad', mdata=mdata)
    """
    model = mofa.mofa_model(model_path)
    
    try:
        # Get factor scores
        Z = np.array(model.expectations["Z"]["group1"])
        
        # Create DataFrame with factor scores
        n_factors = Z.shape[0]
        factor_names = [f'Factor {i+1}' for i in range(n_factors)]
        
        scores_df = pd.DataFrame(Z.T, columns=factor_names)
        
        # If MuData object is provided and metadata should be included
        if mdata is not None and include_metadata:
            # Ensure indices match
            if len(scores_df) == len(mdata.obs):
                scores_df.index = mdata.obs.index
                # Add metadata columns
                for col in mdata.obs.columns:
                    scores_df[col] = mdata.obs[col].values
            else:
                warnings.warn(
                    f"Mismatch in number of samples between model ({len(scores_df)}) "
                    f"and MuData ({len(mdata.obs)}). Metadata not included."
                )
        
    finally:
        if close_model:
            model.close()
            
    return scores_df


def extract_factor_loadings(
    model_path: str,
    view_names: Optional[list] = None,
    close_model: bool = True
) -> Dict[str, pd.DataFrame]:
    """
    Extract factor loadings from a MOFA model.
    
    Parameters
    ----------
    model_path : str
        Path to the saved MOFA model file
    view_names : list, optional
        List of view names to extract loadings for. If None, extracts all views
    close_model : bool, default True
        Whether to close the model after extraction
        
    Returns
    -------
    Dict[str, pd.DataFrame]
        Dictionary mapping view names to DataFrames of factor loadings
        
    Examples
    --------
    >>> import mibi_analysis as ma
    >>> loadings = ma.extract_factor_loadings('mofa_model.h5ad')
    >>> macrophage_loadings = loadings['Macrophage']
    """
    model = mofa.mofa_model(model_path)
    
    try:
        loadings_dict = {}
        
        views_to_extract = view_names if view_names is not None else model.views
        
        for view_name in views_to_extract:
            if view_name not in model.views:
                warnings.warn(f"View '{view_name}' not found in model. Skipping.")
                continue
                
            # Get loadings matrix for this view
            W = model.expectations["W"][view_name]
            
            # Create DataFrame
            n_factors = W.shape[1]
            factor_names = [f'Factor {i+1}' for i in range(n_factors)]
            
            loadings_df = pd.DataFrame(W, columns=factor_names)
            loadings_dict[view_name] = loadings_df
            
    finally:
        if close_model:
            model.close()
            
    return loadings_dict