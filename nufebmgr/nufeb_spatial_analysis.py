import warnings
import polars as pl
import numpy as np
from sklearn.neighbors import KDTree

def neighbors_radius(df: pl.DataFrame, radius: float, periodicity: str,
                     xlen: float=None, ylen:float=None, zlen:float=None) -> dict[str, list[int]]:
    """
    Find all neighbours for each point within a radius.
    Points are specified as their index into the dataframe.to_numpy()
    Periodicity must be specified and currently supports "none" or "xy".

    Warns if dimensions are set when not needed

    See test_nufeb_spatial_analysis.py for example usages.

    :raises Value error if periodicity is not recognized

    :param df: a polars dataframe with columns x, y, and z f
    :param radius: Radius within which to search
    :param periodicity: Either "none" or "xy". "xy" indicates the x and y directions wrap around in a toroidal geometry
    :param xlen: x-dimension length (required for 'xy' periodicity)
    :param ylen: x-dimension length (required for 'xy' periodicity)
    :param zlen: z-dimension length (not yet required )
    :return: A dictionary of lists. Each key is the index of a point in the dataframe. The list items are the indices
    of neighbors within the search radius, accounting for periodicity. The items are sorted in order of increasing
    distance
    """
    # others may be supported in the future
    periodicities = ["none", "xy"]
    # explicit checking, combined with non-default argument is intentional. Users WILL run into issues with any
    # default
    if periodicity not in periodicities:
        raise ValueError(f'Unrecognized periodicity: {periodicity}. Must be one of {",".join(periodicities)}')

    coords = df.select(['x','y','z']).to_numpy()
    if periodicity == "none":
        if xlen is not None or ylen is not None or zlen is not None:
            warnings.warn("Either xlen, ylen, or zlen is set but is not needed for no periodicity. Are you sure you're asking for what you're expecting?", UserWarning)
        tree = KDTree(coords)
        indices, dists = tree.query_radius(coords, r=radius,  return_distance = True, sort_results=True)
        neighbor_lists = {}
        for i, index in enumerate(indices):
            ilist = index.tolist()
            ilist.remove(i)
            neighbor_lists[i] = ilist
        return neighbor_lists
    if periodicity == "xy":
        if xlen is None and ylen is None:
            raise ValueError('Periodicity of "xy" specified but xlen and ylen are not set')
        if xlen is None:
            raise ValueError('Periodicity of "xy" specified but xlen is not set')
        if ylen is None:
            raise ValueError('Periodicity of "xy" specified but ylen is not set')
        if xlen <= 0 and ylen <=0:
            raise ValueError(f'Periodicity of "xy" specified but xlen and ylen not > 0. xlen: {xlen}, ylen: {ylen}')
        if xlen <= 0:
            raise ValueError(f'Periodicity of "xy" specified but xlen is not > 0. xlen: {xlen}')
        if ylen <= 0:
            raise ValueError(f'Periodicity of "xy" specified but ylen is not > 0. ylen: {ylen}')
        offsets = np.array([
            [0, 0],
            [xlen, 0],
            [-xlen, 0],
            [0, ylen],
            [0, -ylen],
            [xlen, ylen],
            [xlen, -ylen],
            [-xlen, ylen],
            [-xlen, -ylen]
        ])
        tiled_coords_list = []
        for dx, dy in offsets:
            # Shift x and y by tile offsets, leave z unchanged
            shifted = coords.copy()
            shifted[:, 0] += dx
            shifted[:, 1] += dy
            tiled_coords_list.append(shifted)
        tiled_coords = np.vstack(tiled_coords_list)
        tree = KDTree(tiled_coords)
        indices, dists = tree.query_radius(coords, r=radius, return_distance=True, sort_results=True)
        pass
    else:
        raise ValueError(f'Unrecognized periodicity: {periodicity}. Must be one of {",".join(periodicities)}')