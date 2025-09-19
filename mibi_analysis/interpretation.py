"""
Interpretation functions for MIBI factor analysis results.

This module provides functions for statistical testing and interpretation of 
factor scores with clinical variables.
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import kruskal
from typing import Optional, Dict, List, Tuple, Union
import warnings


def test_factor_associations(
    factor_scores: pd.DataFrame,
    clinical_variables: List[str],
    factors: Optional[List[str]] = None,
    test_type: str = 'auto',
    correction_method: str = 'bonferroni',
    alpha: float = 0.05
) -> pd.DataFrame:
    """
    Test associations between factors and clinical variables.
    
    This function performs statistical tests to identify significant associations
    between factor scores and clinical/metadata variables.
    
    Parameters
    ----------
    factor_scores : pd.DataFrame
        DataFrame with factor scores and clinical metadata
    clinical_variables : List[str]
        List of clinical variable column names to test
    factors : List[str], optional
        List of factor column names to test. If None, tests all Factor columns
    test_type : str, default 'auto'
        Statistical test type: 'auto', 'kruskal', 'anova', 'ttest', 'chi2'
    correction_method : str, default 'bonferroni'
        Multiple testing correction method
    alpha : float, default 0.05
        Significance level
        
    Returns
    -------
    pd.DataFrame
        Results DataFrame with test statistics and p-values
        
    Examples
    --------
    >>> import mibi_analysis as ma
    >>> # Test associations with clinical variables
    >>> results = ma.test_factor_associations(
    ...     factor_scores, 
    ...     clinical_variables=['Stage', 'Sex', 'Age']
    ... )
    >>> # Show significant associations
    >>> significant = results[results['p_adjusted'] < 0.05]
    """
    # Identify factor columns if not specified
    if factors is None:
        factors = [col for col in factor_scores.columns if col.startswith('Factor')]
    
    # Check that all specified columns exist
    missing_vars = [var for var in clinical_variables if var not in factor_scores.columns]
    missing_factors = [factor for factor in factors if factor not in factor_scores.columns]
    
    if missing_vars:
        raise ValueError(f"Clinical variables not found: {missing_vars}")
    if missing_factors:
        raise ValueError(f"Factors not found: {missing_factors}")
    
    results = []
    
    for clinical_var in clinical_variables:
        for factor in factors:
            # Get data for this test, removing missing values
            test_data = factor_scores[[clinical_var, factor]].dropna()
            
            if len(test_data) < 3:
                warnings.warn(f"Insufficient data for {factor} vs {clinical_var} test")
                continue
            
            # Determine appropriate test based on variable types
            if test_type == 'auto':
                # Check if clinical variable is categorical or continuous
                unique_values = test_data[clinical_var].nunique()
                total_values = len(test_data[clinical_var])
                
                if unique_values <= 10 or test_data[clinical_var].dtype == 'object':
                    # Categorical variable - use Kruskal-Wallis
                    current_test = 'kruskal'
                else:
                    # Continuous variable - use correlation
                    current_test = 'correlation'
            else:
                current_test = test_type
            
            # Perform the test
            test_result = _perform_single_test(
                test_data[factor], 
                test_data[clinical_var], 
                current_test
            )
            
            results.append({
                'factor': factor,
                'clinical_variable': clinical_var,
                'test_type': current_test,
                'statistic': test_result['statistic'],
                'p_value': test_result['p_value'],
                'n_samples': len(test_data),
                'effect_size': test_result.get('effect_size', np.nan)
            })
    
    results_df = pd.DataFrame(results)
    
    if len(results_df) == 0:
        return results_df
    
    # Apply multiple testing correction
    if correction_method == 'bonferroni':
        results_df['p_adjusted'] = results_df['p_value'] * len(results_df)
        results_df['p_adjusted'] = np.minimum(results_df['p_adjusted'], 1.0)
    elif correction_method == 'fdr':
        from statsmodels.stats.multitest import fdrcorrection
        _, p_adjusted = fdrcorrection(results_df['p_value'], alpha=alpha)
        results_df['p_adjusted'] = p_adjusted
    else:
        results_df['p_adjusted'] = results_df['p_value']
    
    # Add significance indicator
    results_df['significant'] = results_df['p_adjusted'] < alpha
    
    # Sort by adjusted p-value
    results_df = results_df.sort_values('p_adjusted')
    
    return results_df


def _perform_single_test(factor_values: pd.Series, clinical_values: pd.Series, test_type: str) -> Dict:
    """
    Perform a single statistical test between factor and clinical variable.
    """
    if test_type == 'kruskal':
        # Kruskal-Wallis test for multiple groups
        groups = clinical_values.unique()
        if len(groups) < 2:
            return {'statistic': np.nan, 'p_value': 1.0}
        
        group_data = [factor_values[clinical_values == group] for group in groups]
        group_data = [group for group in group_data if len(group) > 0]
        
        if len(group_data) < 2:
            return {'statistic': np.nan, 'p_value': 1.0}
        
        statistic, p_value = kruskal(*group_data)
        
        # Calculate eta-squared as effect size
        n_total = len(factor_values)
        k_groups = len(groups)
        eta_squared = (statistic - k_groups + 1) / (n_total - k_groups)
        
        return {
            'statistic': statistic,
            'p_value': p_value,
            'effect_size': eta_squared
        }
    
    elif test_type == 'correlation':
        # Pearson correlation for continuous variables
        correlation, p_value = stats.pearsonr(factor_values, clinical_values)
        return {
            'statistic': correlation,
            'p_value': p_value,
            'effect_size': correlation**2  # R-squared
        }
    
    elif test_type == 'anova':
        # One-way ANOVA
        groups = clinical_values.unique()
        if len(groups) < 2:
            return {'statistic': np.nan, 'p_value': 1.0}
        
        group_data = [factor_values[clinical_values == group] for group in groups]
        group_data = [group for group in group_data if len(group) > 0]
        
        if len(group_data) < 2:
            return {'statistic': np.nan, 'p_value': 1.0}
        
        statistic, p_value = stats.f_oneway(*group_data)
        
        return {
            'statistic': statistic,
            'p_value': p_value
        }
    
    elif test_type == 'ttest':
        # Two-sample t-test (assumes exactly 2 groups)
        groups = clinical_values.unique()
        if len(groups) != 2:
            warnings.warn(f"t-test requires exactly 2 groups, found {len(groups)}")
            return {'statistic': np.nan, 'p_value': 1.0}
        
        group1_data = factor_values[clinical_values == groups[0]]
        group2_data = factor_values[clinical_values == groups[1]]
        
        statistic, p_value = stats.ttest_ind(group1_data, group2_data)
        
        # Calculate Cohen's d as effect size
        pooled_std = np.sqrt(((len(group1_data)-1)*group1_data.var() + 
                             (len(group2_data)-1)*group2_data.var()) / 
                            (len(group1_data) + len(group2_data) - 2))
        cohens_d = (group1_data.mean() - group2_data.mean()) / pooled_std
        
        return {
            'statistic': statistic,
            'p_value': p_value,
            'effect_size': cohens_d
        }
    
    else:
        raise ValueError(f"Unknown test type: {test_type}")


def annotate_factors(
    factor_scores: pd.DataFrame,
    clinical_metadata: pd.DataFrame,
    on: Optional[str] = None,
    how: str = 'left'
) -> pd.DataFrame:
    """
    Add clinical metadata to factor scores DataFrame.
    
    Parameters
    ----------
    factor_scores : pd.DataFrame
        DataFrame with factor scores
    clinical_metadata : pd.DataFrame
        DataFrame with additional clinical metadata
    on : str, optional
        Column name to join on. If None, uses DataFrame indices
    how : str, default 'left'
        Type of join to perform
        
    Returns
    -------
    pd.DataFrame
        Factor scores with additional clinical metadata
        
    Examples
    --------
    >>> import mibi_analysis as ma
    >>> # Add additional clinical metadata
    >>> annotated_scores = ma.annotate_factors(factor_scores, clinical_data)
    """
    if on is None:
        # Join on index
        return factor_scores.join(clinical_metadata, how=how, rsuffix='_meta')
    else:
        # Join on specified column
        return factor_scores.merge(clinical_metadata, on=on, how=how, suffixes=('', '_meta'))


def summarize_factor_by_group(
    factor_scores: pd.DataFrame,
    factor: str,
    group_by: str,
    summary_stats: List[str] = ['mean', 'std', 'median', 'count']
) -> pd.DataFrame:
    """
    Summarize factor scores by groups.
    
    Parameters
    ----------
    factor_scores : pd.DataFrame
        DataFrame with factor scores and metadata
    factor : str
        Factor column name to summarize
    group_by : str
        Column name to group by
    summary_stats : List[str], default ['mean', 'std', 'median', 'count']
        Statistics to calculate
        
    Returns
    -------
    pd.DataFrame
        Summary statistics by group
        
    Examples
    --------
    >>> import mibi_analysis as ma
    >>> # Summarize Factor 1 by cancer stage
    >>> summary = ma.summarize_factor_by_group(scores, 'Factor 1', 'Stage')
    """
    if factor not in factor_scores.columns:
        raise ValueError(f"Factor '{factor}' not found in factor_scores")
    if group_by not in factor_scores.columns:
        raise ValueError(f"Group variable '{group_by}' not found in factor_scores")
    
    # Remove missing values
    clean_data = factor_scores[[factor, group_by]].dropna()
    
    # Calculate summary statistics
    summary = clean_data.groupby(group_by)[factor].agg(summary_stats)
    
    return summary


def identify_top_loadings(
    loadings: Dict[str, pd.DataFrame],
    factor: str,
    n_top: int = 10,
    by_absolute: bool = True
) -> Dict[str, pd.DataFrame]:
    """
    Identify top feature loadings for a specific factor across views.
    
    Parameters
    ----------
    loadings : Dict[str, pd.DataFrame]
        Dictionary of factor loadings by view
    factor : str
        Factor name to analyze
    n_top : int, default 10
        Number of top loadings to return per view
    by_absolute : bool, default True
        Whether to rank by absolute loading values
        
    Returns
    -------
    Dict[str, pd.DataFrame]
        Top loadings for each view
        
    Examples
    --------
    >>> import mibi_analysis as ma
    >>> loadings = ma.extract_factor_loadings('mofa_model.h5ad')
    >>> top_loadings = ma.identify_top_loadings(loadings, 'Factor 1', n_top=5)
    """
    top_loadings = {}
    
    for view_name, view_loadings in loadings.items():
        if factor not in view_loadings.columns:
            warnings.warn(f"Factor '{factor}' not found in view '{view_name}'. Skipping.")
            continue
        
        factor_loadings = view_loadings[factor]
        
        if by_absolute:
            # Sort by absolute loading values
            sorted_loadings = factor_loadings.abs().sort_values(ascending=False)
        else:
            # Sort by loading values (preserving sign)
            sorted_loadings = factor_loadings.sort_values(ascending=False)
        
        # Get top n features
        top_features = sorted_loadings.head(n_top)
        
        # Create result DataFrame with original loading values
        result_df = pd.DataFrame({
            'feature_index': top_features.index,
            'loading': factor_loadings[top_features.index],
            'abs_loading': factor_loadings[top_features.index].abs(),
            'rank': range(1, len(top_features) + 1)
        })
        
        top_loadings[view_name] = result_df
    
    return top_loadings


def calculate_factor_stability(
    factor_scores: pd.DataFrame,
    factors: Optional[List[str]] = None,
    bootstrap_n: int = 1000,
    sample_frac: float = 0.8
) -> pd.DataFrame:
    """
    Assess factor stability through bootstrap resampling.
    
    Parameters
    ----------
    factor_scores : pd.DataFrame
        DataFrame with factor scores
    factors : List[str], optional
        List of factors to analyze. If None, analyzes all Factor columns
    bootstrap_n : int, default 1000
        Number of bootstrap samples
    sample_frac : float, default 0.8
        Fraction of samples to use in each bootstrap
        
    Returns
    -------
    pd.DataFrame
        Stability metrics for each factor
        
    Examples
    --------
    >>> import mibi_analysis as ma
    >>> stability = ma.calculate_factor_stability(factor_scores, bootstrap_n=500)
    """
    if factors is None:
        factors = [col for col in factor_scores.columns if col.startswith('Factor')]
    
    stability_results = []
    
    for factor in factors:
        if factor not in factor_scores.columns:
            warnings.warn(f"Factor '{factor}' not found. Skipping.")
            continue
        
        factor_data = factor_scores[factor].dropna()
        original_mean = factor_data.mean()
        original_std = factor_data.std()
        
        bootstrap_means = []
        bootstrap_stds = []
        
        for _ in range(bootstrap_n):
            # Bootstrap sample
            bootstrap_sample = factor_data.sample(
                n=int(len(factor_data) * sample_frac), 
                replace=True
            )
            bootstrap_means.append(bootstrap_sample.mean())
            bootstrap_stds.append(bootstrap_sample.std())
        
        bootstrap_means = np.array(bootstrap_means)
        bootstrap_stds = np.array(bootstrap_stds)
        
        stability_results.append({
            'factor': factor,
            'original_mean': original_mean,
            'original_std': original_std,
            'bootstrap_mean_mean': np.mean(bootstrap_means),
            'bootstrap_mean_std': np.std(bootstrap_means),
            'bootstrap_std_mean': np.mean(bootstrap_stds),
            'bootstrap_std_std': np.std(bootstrap_stds),
            'mean_stability': 1 - (np.std(bootstrap_means) / abs(original_mean)) if original_mean != 0 else np.nan,
            'std_stability': 1 - (np.std(bootstrap_stds) / original_std) if original_std != 0 else np.nan
        })
    
    return pd.DataFrame(stability_results)