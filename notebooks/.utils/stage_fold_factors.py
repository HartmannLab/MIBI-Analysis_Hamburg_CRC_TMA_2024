"""MuVI factors fit within the stage prediction splits of StageXGB.ipynb.

Run with the liana-squidpy environment, which holds muvi and the squidpy scores. spatial/EncodingClinicalInformation.ipynb
calls this script, because its own environment, pam-keras3, has no muvi.

    python stage_fold_factors.py export
        Writes the spatial block of each feature set of SpatialFeatureContribution.ipynb, one core per row.
    python stage_fold_factors.py fit --worker 0 --workers 14 --threads 2
        Fits MuVI on the training cores of each split, feature set and seed handled by this worker, and projects
        every core of the split onto the factors. Each result is cached in its own file, so reruns only fill gaps.
        The fits run on the CPU: the model is small, so a GPU spends its time launching kernels, and separate
        processes cannot share it without CUDA MPS. Several CPU workers are much faster.
"""
import argparse
import logging
import os
import re

import anndata as ad
import muon as mu
import muvi
import numpy as np
import pandas as pd
import torch

logging.getLogger("muvi").setLevel(logging.WARNING)

data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../data/")
output_dir = data_dir + "spatial_features/"
splits_path = output_dir + "stage_splits.csv"
encoding_dir = output_dir + "stage_encodings/"
fold_factor_dir = output_dir + "stage_fold_factors/"

n_factors = 10
seeds = [0, 1, 2, 3, 4]
factor_order = [f"Factor {i}" for i in range(1, n_factors + 1)]

# Feature classes and feature sets, as in SpatialFeatureContribution.ipynb
spatial_classes = ["MISTy juxtaview", "MISTy paraview", "COZI enrichment", "Squidpy enrichment"]


def feature_class(name):
    if name.startswith("cozi"):
        return "COZI enrichment"
    if name.startswith("squidpy"):
        return "Squidpy enrichment"
    if name.startswith("juxtaview"):
        return "MISTy juxtaview"
    if name.startswith("paraview"):
        return "MISTy paraview"
    if name == "proportion":
        return "Composition"
    if re.search(r"(eccentricity|perimeter|area)", name):
        return "Morphology"
    if name in ["PDL1", "STING1", "MSH2", "MSH6"]:
        return "Functional"
    return "Metabolism"


cell_types = ["CD163_Macrophage", "CD68_Macrophage", "Other_immune_cell", "Endothelial_cell", "Cancer_cell",
              "CD4_Tcell", "CD8_Tcell", "NK_cell", "Monocyte", "Neutrophil", "APC", "CAF", "B_cell", "T_reg_cell"]


def target_and_predictor(name):
    '''Split <prefix>.<target>_<neighbour>, where both cell types may contain underscores'''
    pair = name.split(".")[-1]
    for target in sorted(cell_types, key=len, reverse=True):
        if pair.startswith(target + "_"):
            return target, pair[len(target) + 1:]
    return pair, ""


coarse_of = {fine: coarse for coarse, fines in {
    "Other_immune_cell": ["APC", "B_cell", "Neutrophil", "Other_immune_cell"],
    "Fibroblast": ["CAF"],
    "Macrophage": ["CD163_Macrophage", "CD68_Macrophage"],
    "CD4_lymphocyte": ["CD4_Tcell", "T_reg_cell"],
    "Epithelial_cell": ["Cancer_cell"],
    "Endothelial_cell": ["Endothelial_cell"],
    "Monocyte": ["Monocyte"],
    "Cytotoxic_lymphocyte": ["NK_cell", "CD8_Tcell"],
}.items() for fine in fines}

set_policies = {"MISTy, all features": "all", "MISTy, juxtaview features": "juxtaview",
                "MISTy, own-type and cancer neighbours": "self and cancer",
                "COZI enrichment": "cozi", "Squidpy enrichment": "squidpy", "No spatial": "none"}


def load_feature_sets():
    features = mu.read_h5mu(data_dir + "celltype_features_with_functional_markers.h5mu")
    view_names = list(features.mod.keys())
    published = {v: features.mod[v].to_df() for v in view_names}

    spatial = pd.DataFrame([{"view": v, "feature": n} for v in view_names for n in published[v].columns
                            if feature_class(n) in spatial_classes])
    spatial[["target", "predictor"]] = [target_and_predictor(n) for n in spatial.feature]
    compact = (spatial.target == spatial.predictor) | (spatial.predictor == "Cancer_cell")
    compact_spatial = set(zip(spatial.loc[compact, "view"], spatial.loc[compact, "feature"]))

    def keep_misty(view, name, policy):
        if policy == "all":
            return True
        if policy == "juxtaview":
            return feature_class(name) == "MISTy juxtaview"
        if policy == "self and cancer":
            return (view, name) in compact_spatial
        return False

    cozi = pd.read_csv(data_dir + "cozi_celltypes.csv", index_col="fov")[["interaction", "zscore"]]
    cozi.columns = ["interaction", "value"]
    cozi[["target", "neighbour"]] = [target_and_predictor("." + i) for i in cozi.interaction]
    cozi = cozi[cozi.target.isin(coarse_of) & cozi.neighbour.isin(coarse_of)]
    cozi["view"] = cozi.target.map(coarse_of)
    cozi["feature"] = "cozi." + cozi.neighbour.map(coarse_of)
    cozi_features = {v: (frame.pivot_table(index=frame.index, columns="feature", values="value", aggfunc="mean")
                         .reindex(features.obs_names))
                     for v, frame in cozi.groupby("view")}

    squidpy_long = pd.read_csv(output_dir + "squidpy_nhood_enrichment.csv")
    squidpy_long["zscore"] = squidpy_long.zscore.replace([np.inf, -np.inf], np.nan)
    squidpy_long["feature"] = "squidpy." + squidpy_long.neighbour
    squidpy_features = {v: frame.pivot_table(index="fov", columns="feature", values="zscore")
                        .reindex(features.obs_names)
                        for v, frame in squidpy_long.groupby("view")}

    def build(policy):
        out = {}
        for view in view_names:
            frame = published[view]
            kept = [c for c in frame.columns if feature_class(c) not in spatial_classes or keep_misty(view, c, policy)]
            out[view] = frame[kept]
            if policy == "cozi":
                out[view] = pd.concat([out[view], cozi_features[view]], axis=1)
            if policy == "squidpy":
                out[view] = pd.concat([out[view], squidpy_features[view]], axis=1)
        return out

    return {label: build(policy) for label, policy in set_policies.items()}


def slug(label):
    return re.sub(r"[^a-z0-9]+", "_", label.lower()).strip("_")


def export_encodings(feature_sets):
    '''The spatial features of each feature set, all views side by side, named view:feature'''
    os.makedirs(encoding_dir, exist_ok=True)
    for label, views in feature_sets.items():
        block = pd.concat([frame[[c for c in frame.columns if feature_class(c) in spatial_classes]]
                           .add_prefix(view + ":") for view, frame in views.items()], axis=1)
        if block.shape[1]:
            block.rename_axis("fov").to_csv(encoding_dir + slug(label) + ".csv")


def fit_and_project(views, train, cores, seed, device):
    '''Fit MuVI on the training cores and return the posterior mean factor scores of all cores.

    Each feature is centred and scaled with the training cores only. Features without variance in the training
    cores are dropped. The scores of every core, training cores included, come from the same projection: the
    posterior mean of z given the loadings W and the noise scale sigma of the fit, with the prior z ~ N(0, I) of MuVI
    and only the features observed in that core.'''
    modalities, scaling = {}, {}
    for view, frame in views.items():
        x = frame.loc[train].values.astype(np.float64)
        mean, std = np.nanmean(x, axis=0), np.nanstd(x, axis=0)
        kept = np.isfinite(std) & (std > 0)
        scaling[view] = (frame.columns[kept], mean[kept], std[kept])
        adata = ad.AnnData((x[:, kept] - mean[kept]) / std[kept])
        adata.obs_names, adata.var_names = train, frame.columns[kept]
        modalities[view] = adata
    model = muvi.tl.from_mdata(mu.MuData(modalities), n_factors=n_factors, nmf=False, device=device)
    model.fit(seed=seed, verbose=False)

    loadings = model.get_factor_loadings()
    sigma = {v: np.ravel(s) for v, s in zip(model.view_names, model._guide.get_sigma(as_list=True))}
    view_scale = dict(zip(model.view_names, model._model.view_scales))
    scores = np.zeros((len(cores), n_factors))
    for i, core in enumerate(cores):
        precision, projected = np.eye(n_factors), np.zeros(n_factors)
        for view in model.view_names:
            columns, mean, std = scaling[view]
            y = (views[view].loc[core, columns].values.astype(np.float64) - mean) / std
            observed = np.isfinite(y)
            weighted = loadings[view][:, observed] / sigma[view][observed] ** 2 * view_scale[view]
            precision += weighted @ loadings[view][:, observed].T
            projected += weighted @ y[observed]
        scores[i] = np.linalg.solve(precision, projected)

    fitted = pd.DataFrame(model.get_factor_scores(), index=train)
    check = min(np.corrcoef(fitted.values[:, k], scores[np.isin(cores, train), k])[0, 1] for k in range(n_factors))
    del model
    torch.cuda.empty_cache()
    return pd.DataFrame(scores, index=pd.Index(cores, name="fov"), columns=factor_order), check


def fit_all(feature_sets, worker, workers, device):
    os.makedirs(fold_factor_dir, exist_ok=True)
    splits = pd.read_csv(splits_path, index_col="fov")
    split_names = [c for c in splits.columns if c != "Stage"]
    tasks = [(label, split, seed) for label in feature_sets for split in split_names for seed in seeds]
    for i, (label, split, seed) in enumerate(tasks):
        if i % workers != worker:
            continue
        path = fold_factor_dir + f"{slug(label)}__{slug(split)}__seed{seed}.csv"
        if os.path.exists(path):
            continue
        train = splits.index[splits[split] == "train"]
        cores = splits.index[splits[split].isin(["train", "test"])]
        scores, check = fit_and_project(feature_sets[label], train, cores, seed, device)
        scores.insert(0, "role", splits.loc[cores, split].values)
        scores.insert(0, "projection check", check)
        scores.insert(0, "seed", seed)
        scores.insert(0, "split", split)
        scores.insert(0, "feature set", label)
        scores.to_csv(path + ".tmp")
        os.replace(path + ".tmp", path)
        print(label, split, seed, f"lowest correlation of fitted and projected scores: {check:.4f}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["export", "fit"])
    parser.add_argument("--worker", type=int, default=0)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--threads", type=int, default=4, help="torch threads per worker")
    parser.add_argument("--device", default="cpu")
    args = parser.parse_args()
    torch.set_num_threads(args.threads)
    feature_sets = load_feature_sets()
    if args.mode == "export":
        export_encodings(feature_sets)
    else:
        fit_all(feature_sets, args.worker, args.workers, args.device)
