import pandas as pd
import pathlib


def download_cellmarker_markers(
    species_list=["Human"],
    canonical_url="http://bio-bigdata.hrbmu.edu.cn/CellMarker/CellMarker_download_files/file/Cell_marker_All.xlsx",
    computational_url="http://bio-bigdata.hrbmu.edu.cn/CellMarker/CellMarker_download_files/file/Cell_marker_Seq.xlsx",
    save_dir=None,
):
    """
    Function to retrieve cell markers based on a given dataset and species from the CellMarker database.

    Parameters
    ----------
    species_list: List of species to filter from the marker datasets. 'Human' and 'Mouse' are supported.
    canonical_url: URL to download the canonical marker dataset (Excel file).
    computational_url: URL to download the computational marker dataset (Excel file).

    Returns
    -------
    df_markers_canonical_found: DataFrame with canonical markers.
    df_markers_computational_found: DataFrame with computational markers.
    """

    # Process canonical markers
    df_markers_canonical = pd.read_excel(canonical_url)
    df_markers_canonical = df_markers_canonical[
        df_markers_canonical["species"].isin(species_list)
    ]
    df_markers_canonical["cellontology_id"] = df_markers_canonical[
        "cellontology_id"
    ].str.replace("_", ":")
    # df_markers_canonical_found = df_markers_canonical[df_markers_canonical['cellontology_id'].isin(u_cell_types)]

    # Process computational markers
    df_markers_computational = pd.read_excel(computational_url)
    df_markers_computational = df_markers_computational[
        df_markers_computational["species"].isin(species_list)
    ]
    df_markers_computational["cellontology_id"] = df_markers_computational[
        "cellontology_id"
    ].str.replace("_", ":")
    # df_markers_computational_found = df_markers_computational[df_markers_computational['cellontology_id'].isin(u_cell_types)]

    if save_dir is not None:
        save_dir = pathlib.Path(save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)

        df_markers_canonical.to_csv(save_dir / "cellmarker_canonical.csv")
        df_markers_computational.to_csv(save_dir / "cellmarker_computational.csv")

    return df_markers_canonical, df_markers_computational


def load_cellmarker_markers(save_dir):
    """
    Load the hubmap markers data from the save directory.

    Parameters
    ----------
    save_dir : str
        The directory containing the hubmap markers data.

    Returns:
    --------
    df_markers_canonical_found: DataFrame with canonical markers.
    df_markers_computational_found: DataFrame with computational markers.
    """
    save_dir = pathlib.Path(save_dir)
    df_markers_canonical = pd.read_csv(
        save_dir / "cellmarker_canonical.csv", index_col=0
    )
    df_markers_computational = pd.read_csv(
        save_dir / "cellmarker_computational.csv", index_col=0
    )

    return df_markers_canonical, df_markers_computational


def get_cellmarker_markers(
    df_cell_types, df_markers, custom_map={}, column="cellontology_id"
):
    """
    Get the cell type markers for a given set of cell types.

    Parameters
    ----------
    u_cell_types : list
        A list of unique cell types.
    df_markers : pd.DataFrame
        A DataFrame containing the cell type markers and references columns from cellmarker.
    mode : str
        The mode to match the cell types. If "strict", it will only match exact cell types. If not "strict", it will match any cell type that contains the given cell type.
    custom_map : dict
        A dictionary containing custom mappings from one cell type to another.
    column:
        Whether to use 'cell_name' or 'cellontology_id' as the column to match.
    Returns
    -------
    ct_descr : dict
        A dictionary containing the description of each cell type.
    ct_markers : dict
        A dictionary containing the markers for each cell type.
    ct_markers_df : pd.DataFrame
        A DataFrame containing the cell type markers.
    """
    ct_descr = {}
    ct_markers = {}
    for i, (ct_name, ct) in df_cell_types.iterrows():
        if ct in custom_map.keys():
            ct_ = custom_map[ct]
        else:
            ct_ = ct

        ix = df_markers[column].str.upper() == ct_.upper()

        if ix.sum():
            ct_descr[ct_name] = df_markers[ix]
            ct_markers[ct_name] = list(df_markers[ix]["Symbol"].dropna().unique())
            if not len(ct_markers[ct_name]):
                print(ct_name, ct, "has no markers")
        else:
            ct_descr[ct_name] = None
            ct_markers[ct_name] = []
            print(ct_name, ct, "is not found")

    ct_markers_df = []
    for cell_type, genes in ct_markers.items():
        for gene in genes:
            ct_markers_df.append([cell_type, gene])
    ct_markers_df = pd.DataFrame(ct_markers_df, columns=["Cell Type", "Gene"])

    return ct_descr, ct_markers, ct_markers_df
