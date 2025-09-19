"""
MOFA/MuVI wrapper functions for factor analysis of MIBI data.

This module provides convenient wrapper functions for running MOFA and MuVI factor analysis
on MuData objects and extracting results, following the patterns from MOFACell.ipynb.
"""

import numpy as np
import pandas as pd
import mudata as mu
import muon
import mofax as mofa
import muvi
import liana as li
from typing import Optional, Dict, Any, Union
import warnings
import os


def get_device() -> str:
    """
    Get the appropriate device for MuVI computations.
    
    This function attempts to get a free GPU device using muvi.get_free_gpu_idx(),
    following the same pattern as in MOFACell.ipynb. Falls back to CPU if no GPU available.
    
    Returns
    -------
    str
        Device string ('cpu' or 'cuda:X' where X is the GPU index)
        
    Examples
    --------
    >>> import mibi_analysis as ma
    >>> device = ma.get_device()
    >>> print(f"Using device: {device}")
    """
    device = "cpu"
    try:
        device = f"cuda:{muvi.get_free_gpu_idx()}"
    except Exception as e:
        warnings.warn(f"Could not get GPU device, using CPU: {e}")
    return device


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
    mdata: mu.MuData, 
    obsm_key: str = 'X_mofa',
    obs_keys: Optional[list] = None,
    use_liana: bool = True,
    model_path: Optional[str] = None,
    close_model: bool = True
) -> pd.DataFrame:
    """
    Extract factor scores from a MuData object or MOFA model.
    
    This function extracts factor scores following the patterns from MOFACell.ipynb,
    using liana utilities when available or falling back to direct model access.
    
    Parameters
    ----------
    mdata : mu.MuData
        MuData object containing factor scores in obsm
    obsm_key : str, default 'X_mofa'
        Key in mdata.obsm containing factor scores
    obs_keys : list, optional
        List of observation metadata keys to include (e.g., ['Stage', 'Sex'])
    use_liana : bool, default True
        Whether to use liana utilities for extraction (as in MOFACell.ipynb)
    model_path : str, optional
        Path to MOFA model file (used if use_liana=False)
    close_model : bool, default True
        Whether to close the model after extraction (if using model_path)
        
    Returns
    -------
    pd.DataFrame
        DataFrame with factor scores and metadata
        
    Examples
    --------
    >>> import mibi_analysis as ma
    >>> # Extract using liana (preferred method from MOFACell.ipynb)
    >>> scores = ma.extract_factor_scores(mdata, obs_keys=['Stage', 'Sex'])
    >>> 
    >>> # Extract from model file
    >>> scores = ma.extract_factor_scores(mdata, use_liana=False, model_path='mofa_model.h5ad')
    """
    if use_liana and obsm_key in mdata.obsm:
        # Use liana utilities as in MOFACell.ipynb
        if obs_keys is None:
            obs_keys = ['Stage']  # Default from notebook
            
        try:
            factor_scores = li.ut.get_factor_scores(
                mdata, 
                obsm_key=obsm_key, 
                obs_keys=obs_keys
            )
            return factor_scores
        except Exception as e:
            warnings.warn(f"Could not extract factor scores using liana: {e}")
            # Fall back to manual extraction
    
    # Manual extraction or fallback
    if obsm_key in mdata.obsm:
        # Extract from obsm
        factor_matrix = mdata.obsm[obsm_key]
        n_factors = factor_matrix.shape[1]
        factor_names = [f'Factor {i+1}' for i in range(n_factors)]
        
        scores_df = pd.DataFrame(factor_matrix, columns=factor_names, index=mdata.obs.index)
        
        # Add metadata if requested
        if obs_keys:
            available_keys = [key for key in obs_keys if key in mdata.obs.columns]
            if available_keys:
                for key in available_keys:
                    scores_df[key] = mdata.obs[key].values
        
        return scores_df
    
    elif model_path:
        # Extract from model file (original implementation)
        model = mofa.mofa_model(model_path)
        
        try:
            # Get factor scores
            Z = np.array(model.expectations["Z"]["group1"])
            
            # Create DataFrame with factor scores
            n_factors = Z.shape[0]
            factor_names = [f'Factor {i+1}' for i in range(n_factors)]
            
            scores_df = pd.DataFrame(Z.T, columns=factor_names)
            
            # Add metadata if MuData object provided
            if len(scores_df) == len(mdata.obs):
                scores_df.index = mdata.obs.index
                if obs_keys:
                    available_keys = [key for key in obs_keys if key in mdata.obs.columns]
                    for key in available_keys:
                        scores_df[key] = mdata.obs[key].values
        
        finally:
            if close_model:
                model.close()
                
        return scores_df
    
    else:
        raise ValueError(f"No factor scores found in obsm['{obsm_key}'] and no model_path provided")


def extract_factor_loadings(
    mdata: mu.MuData,
    varm_key: str = 'LFs',
    use_liana: bool = True,
    model_path: Optional[str] = None,
    view_names: Optional[list] = None,
    close_model: bool = True
) -> Union[pd.DataFrame, Dict[str, pd.DataFrame]]:
    """
    Extract factor loadings from a MuData object or MOFA model.
    
    This function extracts factor loadings following the patterns from MOFACell.ipynb,
    using liana utilities when available or falling back to direct model access.
    
    Parameters
    ----------
    mdata : mu.MuData
        MuData object containing factor loadings in varm
    varm_key : str, default 'LFs'
        Key in varm containing factor loadings ('LFs' for MOFA, 'MuVI' for MuVI)
    use_liana : bool, default True
        Whether to use liana utilities for extraction (as in MOFACell.ipynb)
    model_path : str, optional
        Path to MOFA model file (used if use_liana=False)
    view_names : list, optional
        List of view names to extract loadings for. If None, extracts all views
    close_model : bool, default True
        Whether to close the model after extraction
        
    Returns
    -------
    pd.DataFrame or Dict[str, pd.DataFrame]
        If use_liana=True, returns a single DataFrame with loadings and view information.
        Otherwise, returns a dictionary mapping view names to DataFrames of factor loadings.
        
    Examples
    --------
    >>> import mibi_analysis as ma
    >>> # Extract using liana (preferred method from MOFACell.ipynb)
    >>> loadings = ma.extract_factor_loadings(mdata, varm_key='LFs')
    >>> 
    >>> # Extract from model file 
    >>> loadings = ma.extract_factor_loadings(mdata, use_liana=False, model_path='mofa_model.h5ad')
    """
    if use_liana:
        # Use liana utilities as in MOFACell.ipynb
        try:
            variable_loadings = li.ut.get_variable_loadings(mdata, varm_key=varm_key)
            
            # Add cell type (view) to the variable loadings as in the notebook
            variable_loadings['view'] = ''
            for view in mdata.mod.keys():
                if view in mdata.varm and varm_key in mdata.varm:
                    view_mask = mdata.varm[view] if hasattr(mdata.varm[view], '__len__') else np.arange(len(mdata.mod[view].var))
                    variable_loadings.loc[view_mask, "view"] = view
            
            return variable_loadings
        except Exception as e:
            warnings.warn(f"Could not extract variable loadings using liana: {e}")
            # Fall back to manual extraction
    
    # Manual extraction or fallback
    if model_path:
        # Extract from model file (original implementation)
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
    
    else:
        raise ValueError(f"Cannot extract loadings: use_liana=False but no model_path provided")