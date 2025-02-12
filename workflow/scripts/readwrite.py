import yaml
import pandas as pd
import os
import json
import h5py
import numpy as np
import scipy
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed

try:
    import msgspec
except ImportError:
    pass

script_dir = os.path.dirname(os.path.abspath(__file__))
config_path = os.path.join(script_dir, "../../config/config.yml")


def config(path=config_path):
    """
    Read the configuration file and return a dictionary of config values.

    Parameters
    ----------
    path : str
        The path to the configuration file. Defaults to the value of
        `config_path` if not provided.

    Returns
    -------
    cfg : dict
        A dictionary of configuration values. All values are strings and
        have been converted to absolute paths by prepending the value of
        `cfg["base_dir"]`.
    """
    with open(path, "r") as ymlfile:
        cfg = yaml.safe_load(ymlfile)
        for k in cfg.keys():
            cfg[k] = os.path.join(cfg["base_dir"], cfg[k])
    return cfg


def soma_to_anndata(soma_uri, measurement_name, X_layer_name, return_experiment=False, **kwargs):
    """
    Export a SOMA experiment to an anndata object.

    Parameters
    ----------
    soma_uri (str): The URI (path) of the SOMA experiment.
    measurement_name (str): The measurement name (e.g., 'RNA') to extract.
    X_layer_name (str): The layer name in the X matrix to extract (e.g., 'counts').
    return_experiment (bool): Whether to return the SOMA experiment.

    Returns
    -------
    anndata.AnnData: The anndata object.
    """
    import tiledbsoma as soma
    import tiledbsoma.io

    # Open the SOMA experiment
    experiment = soma.open(soma_uri)
    # Export the SOMA experiment to anndata format
    ad = soma.io.to_anndata(
        experiment=experiment,
        measurement_name=measurement_name,
        X_layer_name=X_layer_name,
        **kwargs,
    )

    if return_experiment:
        return ad, experiment
    else:
        return ad


def xenium_specs(path):
    path = Path(path)
    with open(path / "experiment.xenium") as f:
        specs = json.load(f)
    return specs


######### Xenium readers
def xenium_samples_files(dir_segmentation_condition, segmentation=None, samples=None):
    """
    Get a dictionary of files for each sample in a Xenium segmentation run.

    Parameters
    ----------
    dir_segmentation_condition (str): The directory path of the segmentation run.
    segmentation (str): The segmentation name, e.g., 'default', '10x_5um', 'baysor'.
    samples (list): The sample names to include. If None, include all samples.

    Returns
    -------
    dict: A dictionary of files for each sample.
    """
    files = {}
    for sample_path in Path(dir_segmentation_condition).iterdir():
        for sample_path in sample_path.iterdir():
            sample_name = sample_path.stem

            if samples is not None and sample_name not in samples:
                continue
            elif "corrupted" in sample_name:
                continue
            else:
                if segmentation != "default":
                    files[sample_name] = sample_path / "normalised_results/outs"
                else:
                    files[sample_name] = sample_path

    return files


def read_xenium_sample(
    sample_name,
    path,
    cells_as_circles=False,
    cells_boundaries=False,
    nucleus_boundaries=False,
    cells_labels=False,
    nucleus_labels=False,
    transcripts=False,
    morphology_mip=False,
    morphology_focus=False,
    aligned_images=False,
    anndata_only=False,
):
    """
    Reads a xenium sample from a directory path.

    Parameters
    ----------
    sample_name (str): The sample name.
    path (str): The directory path of the segmentation run.
    cells_as_circles (bool): Whether to include cell polygons as circles or not.
    cells_boundaries (bool): Whether to include cell boundaries or not.
    nucleus_boundaries (bool): Whether to include nucleus boundaries or not.
    cells_labels (bool): Whether to include cell labels or not.
    nucleus_labels (bool): Whether to include nucleus labels or not.
    transcripts (bool): Whether to include transcript locations or not.
    morphology_mip (bool): Whether to include morphology MIP or not.
    morphology_focus (bool): Whether to include morphology focus or not.
    aligned_images (bool): Whether to include aligned images or not.
    anndata_only (bool): Whether to return only the anndata object or the full spatialdata object.

    Returns
    -------
    If anndata_only, returns a tuple of the sample name and anndata object.
    Otherwise, returns a tuple of the sample name and spatialdata object.
    """
    import spatialdata_io

    xdata = spatialdata_io.xenium(
        path,
        cells_as_circles=cells_as_circles,
        cells_boundaries=cells_boundaries,
        nucleus_boundaries=nucleus_boundaries,
        cells_labels=cells_labels,
        nucleus_labels=nucleus_labels,
        transcripts=transcripts,
        morphology_mip=morphology_mip,
        morphology_focus=morphology_focus,
        aligned_images=aligned_images,
    )

    ad = xdata["table"]
    ad.obs_names = ad.obs["cell_id"].values

    metrics_path = Path(path) / "metrics_summary.csv"
    if metrics_path.exists():
        ad.uns["metrics_summary"] = pd.read_csv(metrics_path)
    else:
        print("metrics_summary.csv not found at:", metrics_path)

    if anndata_only:
        return sample_name, ad
    return sample_name, xdata


def read_xenium_samples(
    data_dirs,
    cells_as_circles=False,
    cells_boundaries=False,
    nucleus_boundaries=False,
    cells_labels=False,
    nucleus_labels=False,
    transcripts=False,
    morphology_mip=False,
    morphology_focus=False,
    aligned_images=False,
    anndata_only=False,
    sample_name_as_key=True,
):
    """
    Reads in a dictionary of sample directories and returns a dictionary of
    AnnData objects or spatialdata objects depending on the anndata_only flag.

    Parameters
    ----------
    data_dirs : dict or list
        A dictionary of sample directories or a list of paths to sample directories.
    cells_as_circles : bool, optional
        Whether to include cell boundary data as circles, by default False
    cells_boundaries : bool, optional
        Whether to include cell boundary data, by default False
    nucleus_boundaries : bool, optional
        Whether to include nucleus boundary data, by default False
    cells_labels : bool, optional
        Whether to include cell labels, by default False
    nucleus_labels : bool, optional
        Whether to include nucleus labels, by default False
    transcripts : bool, optional
        Whether to include transcript data, by default False
    morphology_mip : bool, optional
        Whether to include morphology data at the maximum intensity projection, by default False
    morphology_focus : bool, optional
        Whether to include morphology data at the focus, by default False
    aligned_images : bool, optional
        Whether to include aligned images, by default False
    anndata_only : bool, optional
        Whether to only return an AnnData object, by default False
    sample_name_as_key: bool, optional
        Whether to use the sample name as the key in the return dictionary, otherwise returns full path

    Returns
    -------
    dict
        A dictionary of sample names mapped to AnnData objects or spatialdata objects.
    """
    if isinstance(data_dirs, list):
        sample_names = [Path(path).stem if sample_name_as_key else path for path in data_dirs]
        data_dirs = {sample_name: path for sample_name, path in zip(sample_names, data_dirs)}

    # Parallel processing
    xdatas = {}
    with ProcessPoolExecutor() as executor:
        futures = [
            executor.submit(
                read_xenium_sample,
                sample_name,
                path,
                cells_as_circles,
                cells_boundaries,
                nucleus_boundaries,
                cells_labels,
                nucleus_labels,
                transcripts,
                morphology_mip,
                morphology_focus,
                aligned_images,
                anndata_only,
            )
            for sample_name, path in data_dirs.items()
        ]

        for future in as_completed(futures):
            try:
                sample_name, result = future.result()
                xdatas[sample_name] = result
            except Exception as e:
                print(f"Error processing {e}")

    return xdatas


######### 10x writers


def write_10X_h5(adata, file):
    """Writes adata to a 10X-formatted h5 file.
    taken from https://github.com/scverse/anndata/issues/595

    Note that this function is not fully tested and may not work for all cases.
    It will not write the following keys to the h5 file compared to 10X:
    '_all_tag_keys', 'pattern', 'read', 'sequence'

    Args:
        adata (AnnData object): AnnData object to be written.
        file (str): File name to be written to. If no extension is given, '.h5' is appended.

    Raises:
        FileExistsError: If file already exists.

    Returns:
        None
    """

    if ".h5" not in file:
        file = f"{file}.h5"
    if Path(file).exists():
        raise FileExistsError(f"There already is a file `{file}`.")

    def int_max(x):
        return int(max(np.floor(len(str(int(max(x)))) / 4), 1) * 4)

    def str_max(x):
        return max([len(i) for i in x])

    if not scipy.sparse.issparse(adata.X):
        adata.X = scipy.sparse.csr_matrix(adata.X)
    if "genome" not in adata.var:
        adata.var["genome"] = "undefined"
    if "feature_types" not in adata.var:
        adata.var["feature_types"] = "Gene Expression"
    if "gene_ids" not in adata.var:
        adata.var["gene_ids"] = adata.var_names

    w = h5py.File(file, "w")
    grp = w.create_group("matrix")
    grp.create_dataset(
        "barcodes",
        data=np.array(adata.obs_names, dtype=f"|S{str_max(adata.obs_names)}"),
    )
    grp.create_dataset("data", data=np.array(adata.X.data, dtype=f"<i{int_max(adata.X.data)}"))
    ftrs = grp.create_group("features")
    # this group will lack the following keys:
    # '_all_tag_keys', 'feature_type', 'genome', 'id', 'name', 'pattern', 'read', 'sequence'
    ftrs.create_dataset(
        "feature_type",
        data=np.array(adata.var.feature_types, dtype=f"|S{str_max(adata.var.feature_types)}"),
    )
    ftrs.create_dataset(
        "genome",
        data=np.array(adata.var.genome, dtype=f"|S{str_max(adata.var.genome)}"),
    )
    ftrs.create_dataset(
        "id",
        data=np.array(adata.var.gene_ids, dtype=f"|S{str_max(adata.var.gene_ids)}"),
    )
    ftrs.create_dataset("name", data=np.array(adata.var.index, dtype=f"|S{str_max(adata.var.index)}"))
    grp.create_dataset("indices", data=np.array(adata.X.indices, dtype=f"<i{int_max(adata.X.indices)}"))
    grp.create_dataset("indptr", data=np.array(adata.X.indptr, dtype=f"<i{int_max(adata.X.indptr)}"))
    grp.create_dataset(
        "shape",
        data=np.array(list(adata.X.shape)[::-1], dtype=f"<i{int_max(adata.X.shape)}"),
    )


######### RCTD readers


def read_json(file_path):
    with open(file_path, "r") as file:
        return json.load(file)


def read_json_msgspec(file_path):
    with open(file_path, "rb") as file:
        return msgspec.json.decode(file.read())


def _rds2py_dict_to_df(r_obj_df, mode="results_df"):
    if mode == "results_df":
        r_obj_df_columns = r_obj_df["attributes"]["names"]["data"]
        r_obj_df_index = r_obj_df["attributes"]["row.names"]["data"]
        pandas_df = pd.DataFrame(
            [r_obj_df["data"][i]["data"] for i in range(len(r_obj_df["data"]))],
            index=r_obj_df_columns,
            columns=r_obj_df_index,
        ).T
    elif mode == "weights":
        r_obj_df_columns = r_obj_df["attributes"]["dimnames"]["data"][1]["data"]
        r_obj_df_index = r_obj_df["attributes"]["dimnames"]["data"][0]["data"]
        r_obj_df["data"] = r_obj_df["data"].reshape(r_obj_df["attributes"]["dim"]["data"], order="F")
        pandas_df = pd.DataFrame(r_obj_df["data"], index=r_obj_df_index, columns=r_obj_df_columns)

    return pandas_df


def read_rctd_sample(sample_name, rctd_results_path):
    """
    Reads RCTD results from a single sample and returns a dictionary containing:

    - results_df: a pandas DataFrame with columns to be added to the anndata object's obs
    - weights: a pandas Series with the weights for each cell for the given reference
    - weights_doublet: (not implemented) a pandas Series with the weights for each cell for doublets for the given reference
    - singlet_scores: (not implemented) a pandas Series with the singlet scores for each cell for the given reference

    Parameters
    ----------
    sample_name: str
        The name of the sample
    rctd_results_path : str
        The path to the sample RCTD results
    rsuffix : str, optional
        The suffix to append to the reference name when storing to the anndata objects

    Returns
    -------
    A tuple containing the sample name and the results dictionary
    """
    from rds2py import read_rds

    r_obj = read_rds(rctd_results_path)

    results = r_obj["attributes"]["results"]
    results_keys = results["attributes"]["names"]["data"]
    results_keys_idx = {k: results_keys.index(k) for k in results_keys}

    pandas_results = {}
    for k in ["results_df", "weights"]:
        pandas_results[k] = _rds2py_dict_to_df(results["data"][results_keys_idx[k]], mode=k)

    return sample_name, pandas_results


def read_rctd_samples(ads, rctd_results_paths, prefix=""):
    """
    Read RCTD results into anndata objects in parallel using ProcessPoolExecutor.

    Parameters
    ----------
    ads : dict of anndata.AnnData
        The anndata objects to be updated.
    rctd_results_paths : str
        The directory containing the RCTD results.
    prefix : str, optional
        The prefix to append to the reference name when storing to the anndata objects.

    Returns
    -------
    None
    """

    # Use ProcessPoolExecutor for CPU-bound tasks
    with ProcessPoolExecutor() as executor:
        futures = {
            executor.submit(
                read_rctd_sample,
                sample_name,
                rctd_results_paths[sample_name],
            ): sample_name
            for sample_name in ads.keys()
        }

        # Update anndata objects in the parent process
        for future in as_completed(futures):
            try:
                sample_name, results = future.result()
                if results:
                    ad = ads[sample_name]
                    ad.obs = ad.obs.join(results["results_df"].add_prefix(prefix))
                    ad.uns[f"{prefix}_weights"] = results["weights"]

            except Exception as e:
                print(f"Error processing sample {futures[future]}: {e}")


###### coexpression files readers
def read_coexpression_file(k, method, target_count, results_dir):
    """
    Worker function to read the coexpression and positivity rate parquet for a single sample.

    Parameters
    ----------
    k : tuple
        The sample name tuple (segmentation, condition, panel, sample, sample).
    method : str
        The coexpression method.
    target_count : int
        The target count of the coexpression method.
    results_dir : Path
        The directory containing the coexpression results.

    Returns
    -------
    method : str
        The coexpression method.
    target_count : int
        The target count of the coexpression method.
    cc : pd.DataFrame
        The coexpression matrix.
    pos_rate : pd.Series
        The positivity rate.
    """
    out_file_coexpr = results_dir / f"{'/'.join(k)}/coexpression_{method}_{target_count}.parquet"
    out_file_pos_rate = results_dir / f"{'/'.join(k)}/positivity_rate_{method}_{target_count}.parquet"

    cc = pd.read_parquet(out_file_coexpr)
    pos_rate = pd.read_parquet(out_file_pos_rate)[0]
    return method, target_count, cc, pos_rate


def read_coexpression_files(cc_paths, results_dir):
    """
    Reads coexpression parquet files for multiple methods and target counts in parallel using ThreadPoolExecutor.

    Parameters
    ----------
    cc_paths : list of tuples
        A list of tuples containing the key `k`, i.e., a sample name tuple (segmentation, condition, panel, sample, sample)
        and the method and target count to read.
    results_dir : str
        The directory containing the coexpression results.

    Returns
    -------
    CC : dict
        A dictionary with the coexpression matrices for each method and target count.
    pos_rate : dict
        A dictionary with the positivity rates for each method and target count.
    """
    with ThreadPoolExecutor() as executor:
        futures = [
            executor.submit(read_coexpression_file, k, method, target_count, results_dir)
            for k, method, target_count in cc_paths
        ]

        CC = {}
        pos_rate = {}
        for future in as_completed(futures):
            method, target_count, cc, pr = future.result()
            k = cc_paths[futures.index(future)][0]  # Retrieve the `k` corresponding to this future

            if k not in CC:
                CC[k] = {}
            if k not in pos_rate:
                pos_rate[k] = {}

            CC[k][method, target_count] = cc
            pos_rate[k][method, target_count] = pr
    return CC, pos_rate


def get_gene_panel_info(path):
    with open(path, "r") as f:
        gene_panel = json.load(f)["payload"]["targets"]

    gene_panel_info = pd.DataFrame(columns=["codewords"])
    for i, g in enumerate(gene_panel):
        gene_panel_info.at[i, "gene_coverage"] = g["info"]["gene_coverage"]
        gene_panel_info.at[i, "id"] = g["type"]["data"].get("id")
        gene_panel_info.at[i, "name"] = g["type"]["data"]["name"]
        gene_panel_info.at[i, "codewords"] = g["codewords"]
        gene_panel_info.at[i, "source_category"] = g["source"]["category"]
        gene_panel_info.at[i, "source_design_id"] = g["source"]["identity"]["design_id"]
        gene_panel_info.at[i, "source_name"] = g["source"]["identity"]["name"]
        gene_panel_info.at[i, "source_version"] = g["source"]["identity"].get("version")
    return gene_panel_info
