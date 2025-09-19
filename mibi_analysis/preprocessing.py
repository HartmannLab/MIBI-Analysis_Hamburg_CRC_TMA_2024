"""
Preprocessing functions for MIBI analysis using muon objects.

This module provides functions to normalize and prepare MuData objects for 
factor analysis using MuVI/MOFA.
"""

import numpy as np
import pandas as pd
import mudata as mu
import warnings
from typing import Optional, Union, List


def normalize_mudata_features(
    mdata: mu.MuData,
    method: str = 'z_score',
    ignore_na: bool = True,
    modalities: Optional[List[str]] = None,
    copy: bool = False
) -> Optional[mu.MuData]:
    """
    Normalize features in a MuData object across modalities.
    
    This function applies normalization to each modality in the MuData object,
    centering and scaling the data as done in the MOFACell analysis.
    
    Parameters
    ----------
    mdata : mu.MuData
        The input MuData object containing multiple modalities
    method : str, default 'z_score'
        Normalization method. Currently only 'z_score' (center and scale) is supported
    ignore_na : bool, default True
        Whether to ignore NA values when computing mean and standard deviation
    modalities : List[str], optional
        List of modality names to normalize. If None, normalizes all modalities
    copy : bool, default False
        Whether to return a copy of the MuData object
        
    Returns
    -------
    Optional[mu.MuData]
        Normalized MuData object if copy=True, otherwise modifies in place and returns None
        
    Examples
    --------
    >>> import mibi_analysis as ma
    >>> # Normalize all modalities in place
    >>> ma.normalize_mudata_features(mdata)
    >>> 
    >>> # Normalize specific modalities and return copy
    >>> normalized_mdata = ma.normalize_mudata_features(
    ...     mdata, 
    ...     modalities=['Macrophage', 'T_cell'],
    ...     copy=True
    ... )
    """
    if copy:
        mdata = mdata.copy()
    
    if method != 'z_score':
        raise ValueError(f"Method '{method}' not supported. Only 'z_score' is currently available.")
    
    # Get modalities to process
    modalities_to_process = modalities if modalities is not None else list(mdata.mod.keys())
    
    for modality_name in modalities_to_process:
        if modality_name not in mdata.mod:
            warnings.warn(f"Modality '{modality_name}' not found in MuData object. Skipping.")
            continue
            
        modality_data = mdata.mod[modality_name]
        
        # Center and scale each column of the data
        if ignore_na:
            # Use nanmean and nanstd to ignore NA values
            mean_vals = np.nanmean(modality_data.X, axis=0, keepdims=True)
            std_vals = np.nanstd(modality_data.X, axis=0, keepdims=False)
        else:
            mean_vals = np.mean(modality_data.X, axis=0, keepdims=True)
            std_vals = np.std(modality_data.X, axis=0, keepdims=False)
            
        # Center the data
        modality_data.X = modality_data.X - mean_vals
        
        # Scale by standard deviation (avoid division by zero)
        std_vals = np.where(std_vals == 0, 1, std_vals)
        modality_data.X = modality_data.X / std_vals
        
        # Store normalization parameters in modality metadata
        if 'normalization' not in modality_data.uns:
            modality_data.uns['normalization'] = {}
        modality_data.uns['normalization']['method'] = method
        modality_data.uns['normalization']['mean'] = mean_vals.flatten()
        modality_data.uns['normalization']['std'] = std_vals
        
    if copy:
        return mdata


def prepare_mudata_for_mofa(
    mdata: mu.MuData,
    normalize: bool = True,
    filter_features: bool = True,
    min_variance: float = 0.0,
    remove_constant: bool = True,
    copy: bool = False
) -> Optional[mu.MuData]:
    """
    Prepare a MuData object for MOFA factor analysis.
    
    This function performs common preprocessing steps before running MOFA,
    including normalization and feature filtering.
    
    Parameters
    ----------
    mdata : mu.MuData
        Input MuData object
    normalize : bool, default True
        Whether to normalize features using z-score normalization
    filter_features : bool, default True
        Whether to filter features based on variance
    min_variance : float, default 0.0
        Minimum variance threshold for feature filtering
    remove_constant : bool, default True
        Whether to remove constant features (zero variance)
    copy : bool, default False
        Whether to return a copy of the MuData object
        
    Returns
    -------
    Optional[mu.MuData]
        Prepared MuData object if copy=True, otherwise modifies in place and returns None
        
    Examples
    --------
    >>> import mibi_analysis as ma
    >>> # Basic preparation with normalization
    >>> ma.prepare_mudata_for_mofa(mdata)
    >>> 
    >>> # Prepare with custom variance filtering
    >>> prepared_mdata = ma.prepare_mudata_for_mofa(
    ...     mdata,
    ...     min_variance=0.1,
    ...     copy=True
    ... )
    """
    if copy:
        mdata = mdata.copy()
        
    # Normalize features if requested
    if normalize:
        normalize_mudata_features(mdata, copy=False)
    
    # Filter features if requested
    if filter_features:
        for modality_name, modality_data in mdata.mod.items():
            # Calculate variance for each feature
            feature_var = np.var(modality_data.X, axis=0)
            
            # Create filter mask
            if remove_constant:
                # Remove features with zero variance or below minimum
                keep_features = feature_var > min_variance
            else:
                # Only apply minimum variance filter
                keep_features = feature_var >= min_variance
                
            # Apply filter if any features should be removed
            if not np.all(keep_features):
                n_removed = np.sum(~keep_features)
                n_total = len(keep_features)
                warnings.warn(
                    f"Removing {n_removed}/{n_total} features from modality '{modality_name}' "
                    f"due to low variance (< {min_variance})"
                )
                
                # Filter the data and variable annotations
                modality_data._inplace_subset_var(keep_features)
                
                # Store filtering information
                if 'preprocessing' not in modality_data.uns:
                    modality_data.uns['preprocessing'] = {}
                modality_data.uns['preprocessing']['n_features_removed'] = n_removed
                modality_data.uns['preprocessing']['min_variance_threshold'] = min_variance
    
    if copy:
        return mdata


def get_modality_stats(mdata: mu.MuData) -> pd.DataFrame:
    """
    Get summary statistics for each modality in a MuData object.
    
    Parameters
    ----------
    mdata : mu.MuData
        Input MuData object
        
    Returns
    -------
    pd.DataFrame
        DataFrame with statistics for each modality including number of samples,
        features, and basic data characteristics
        
    Examples
    --------
    >>> import mibi_analysis as ma
    >>> stats = ma.get_modality_stats(mdata)
    >>> print(stats)
    """
    stats_list = []
    
    for modality_name, modality_data in mdata.mod.items():
        n_obs, n_vars = modality_data.shape
        
        # Calculate basic statistics
        data_mean = np.nanmean(modality_data.X)
        data_std = np.nanstd(modality_data.X)
        data_min = np.nanmin(modality_data.X)
        data_max = np.nanmax(modality_data.X)
        n_na = np.sum(np.isnan(modality_data.X))
        na_proportion = n_na / (n_obs * n_vars)
        
        stats_list.append({
            'modality': modality_name,
            'n_samples': n_obs,
            'n_features': n_vars,
            'mean': data_mean,
            'std': data_std,
            'min': data_min,
            'max': data_max,
            'n_missing': n_na,
            'missing_proportion': na_proportion
        })
    
    return pd.DataFrame(stats_list)