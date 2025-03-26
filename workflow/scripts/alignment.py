import uuid
from typing import Callable

import numpy as np
import pandas as pd
import geopandas as gpd

import shapely
from shapely.geometry.base import BaseGeometry
from shapely.geometry import mapping, Polygon, MultiPolygon


def transform_mtx_by_coords(
    coords: np.ndarray,
    transform_mtx: np.ndarray,
) -> np.ndarray:
    """
    Apply transformation matrix to a set of coordinates.

    Parameters
    ----------
    coords : np.ndarray
        Shape (n, 2) array of coordinates to be transformed.
    transform_mtx : np.ndarray
        Shape (3, 3) transformation matrix.

    Returns
    -------
    np.ndarray
        Shape (n, 2) array of transformed coordinates.
    """
    return np.matmul(
        transform_mtx,
        np.c_[
            coords,
            np.ones(coords.shape[0]),
        ].transpose(),
    ).transpose()[:, :2]


def process_polygon(
    geom: Polygon,
    transform_mtx: np.ndarray | None = None,
    mat_mul_func: Callable | None = transform_mtx_by_coords,
    scalers: tuple[float, float] = (1 / 0.2125, 1.0),  # (um2px, px2um)
) -> Polygon:
    """
    Process a shapely Polygon by applying a transformation matrix and/or scaling.

    Parameters
    ----------
    geom : Polygon
        The input Polygon to be processed.
    transform_mtx : np.ndarray, optional
        Shape (3, 3) transformation matrix to be applied to the coordinates of the Polygon.
    mat_mul_func : Callable, optional
        Function to apply the transformation matrix to the coordinates. Defaults to
        `transform_mtx_by_coords`.
    scalers : tuple[float, float], optional
        Tuple of two scaling factors: (um2px, px2um). The first factor is applied before the
        transformation matrix, and the second factor is applied afterwards. Defaults to
        `(1 / 0.2125, 1.0)`.

    Returns
    -------
    Polygon
        The processed Polygon.
    """
    coords = np.array(geom.exterior.coords)
    if transform_mtx is None or mat_mul_func is None:
        ret: Polygon = Polygon(coords * scalers[0] * scalers[1])  # Keep original geometry
    else:
        ret = Polygon(
            mat_mul_func(
                coords * scalers[0],
                transform_mtx,
            )
            * scalers[1]
        )
    return ret


def process_multipolygon(
    geom: MultiPolygon,
    transform_mtx: np.ndarray | None = None,
    mat_mul_func: Callable | None = transform_mtx_by_coords,
    scalers: tuple[float, float] = (1 / 0.2125, 1.0),  # (um2px, px2um)
) -> MultiPolygon:
    """
    Applies the same processing as `process_polygon` to each Polygon of the input MultiPolygon.

    Parameters
    ----------
    geom : MultiPolygon
        The input MultiPolygon to be processed.
    transform_mtx : np.ndarray, optional
        Shape (3, 3) transformation matrix to be applied to the coordinates of the Polygon.
    mat_mul_func : Callable, optional
        Function to apply the transformation matrix to the coordinates. Defaults to
        `transform_mtx_by_coords`.
    scalers : tuple[float, float], optional
        Tuple of two scaling factors: (um2px, px2um). The first factor is applied before the
        transformation matrix, and the second factor is applied afterwards. Defaults to
        `(1 / 0.2125, 1.0)`.

    Returns
    -------
    MultiPolygon
        The processed MultiPolygon.
    """
    ret: list[Polygon] = []
    for i in geom.geoms:
        ret.append(
            process_polygon(
                i,
                transform_mtx,
                mat_mul_func,
                scalers,
            )
        )

    return MultiPolygon(ret)


def process_raw_boundaries(
    df: gpd.GeoDataFrame,
    *,
    object_type: str = "annotation",
    transform_mtx: np.ndarray | None = None,
    mat_mul_func: Callable | None = transform_mtx_by_coords,
    scalers: tuple[float, float] = (1 / 0.2125, 1.0),  # (um2px, px2um)
    first_n: int | None = None,
) -> list:
    """
    Processes a GeoDataFrame of geometries, applying transformations and scaling, and returns
    a list of features in GeoJSON-like format.

    Parameters
    ----------
    df : gpd.GeoDataFrame
        The GeoDataFrame containing geometries to process.
    object_type : str, optional
        The type of object being processed, used as a property in the resulting features.
        Defaults to "annotation".
    transform_mtx : np.ndarray, optional
        A 3x3 transformation matrix to apply to each geometry's coordinates. If not provided,
        no transformation is applied.
    mat_mul_func : Callable, optional
        A function that applies the transformation matrix to the coordinates. Defaults to
        `transform_mtx_by_coords`.
    scalers : tuple[float, float], optional
        Tuple of scaling factors (um2px, px2um). The first factor is applied before any
        transformation, and the second factor is applied after. Defaults to (1 / 0.2125, 1.0).
    first_n : int, optional
        If provided, only the first `n` geometries in the GeoDataFrame are processed.

    Returns
    -------
    list
        A list of dictionaries, each representing a feature with a unique ID, transformed
        geometry, and properties including the object type and feature name.
    """

    features = []

    idx = 0
    for _, row in df.iterrows():
        if first_n is not None and idx >= first_n:
            break

        geom = row["geometry"]

        if geom.is_empty:
            print("Empty shape found. Skipped.")
            continue

        if geom.geom_type == "Polygon":
            _geom: BaseGeometry = process_polygon(
                geom,
                transform_mtx,
                mat_mul_func,
                scalers,
            )
        elif geom.geom_type == "MultiPolygon":
            _geom = process_multipolygon(
                geom,
                transform_mtx,
                mat_mul_func,
                scalers,
            )
        else:
            raise ValueError(f"Error! Do not support geometry type: {geom.geom_type}")

        features.append(
            {
                "type": "Feature",
                "id": str(uuid.uuid4()),  # Generate a random uuid
                "geometry": mapping(_geom),
                "properties": {
                    "name": str(row.name),  # Use the row index or 'cell_id' for feature ID
                    "objectType": object_type,
                },
            }
        )

        idx += 1

    return features


def shape_in_regions(
    shape: BaseGeometry,
    regions: dict[str, BaseGeometry],
) -> bool:
    """
    Checks if a shape is in any of the given regions.

    Parameters
    ----------
    shape : BaseGeometry
        The shape to check.
    regions : dict[str, BaseGeometry]
        A dictionary of region names to their corresponding geometries.

    Returns
    -------
    pandas.Series
        A boolean Series with the region names as index. Each value is True if the shape is
        in the corresponding region, False otherwise. The last value is True if the shape is
        in any of the regions, False otherwise.
    """
    keys: list[str] = list(regions.keys())
    is_in: list[bool] = [shapely.contains(regions[key], shape) for key in keys]

    return pd.Series(
        is_in + [np.any(is_in)],
        index=keys + ["in_regions"],
    )


def split_cells_from_samples(
    gdf: gpd.GeoDataFrame,
    regions: dict[str, BaseGeometry],
):
    """
    Splits a GeoDataFrame of cells into samples based on given regions.

    The function takes a GeoDataFrame of cells and a dictionary of region names to their
    corresponding geometries. It then applies the shape_in_regions function to each cell
    and returns a new GeoDataFrame with the same columns as the original, but with
    additional columns indicating which regions each cell is in.

    The function is intended to be used in a pipeline for processing segmentation data
    from different sources.

    Parameters
    ----------
    gdf : GeoDataFrame
        A GeoDataFrame of cells.
    regions : dict[str, BaseGeometry]
        A dictionary of region names to their corresponding geometries.

    Returns
    -------
    GeoDataFrame
        A GeoDataFrame with the same columns as the original, but with additional columns
        indicating which regions each cell is in.
    """
    return (
        gdf.set_index(
            "id",
        )
        .apply(
            lambda x: shape_in_regions(x.geometry, regions),
            axis=1,
        )
        .reset_index(
            None,
        )
    )
