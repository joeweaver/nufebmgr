import warnings
import polars as pl
import numpy as np
from sklearn.neighbors import KDTree
from scipy.spatial import distance_matrix
from enum import Enum

class Periodicity(Enum):
    """
    Boundary conditions for simulation regarding how bugs and substances 'wrap around'. Currently supports 'none' and
    'x-y' wrapround.

    Attributes:
        XY: Bugs and substances wrap around the sides, but not through the top and bottom (height direction)
        NONE: Boundaries act like walls
    """
    XY = "xy"
    NONE = "none"


def neighbors_radius(df: pl.DataFrame, radius: float, periodicity: Periodicity,
                     xlen: float=None, ylen:float=None, zlen:float=None) -> dict[int, list[int]]:
    """
    Find all neighbours for each point within a radius.
    Points are specified as their index into the dataframe.to_numpy()
    Periodicity must be specified and currently supports "none" or "xy". It is intentionally not a default parameter
    to prevent erroneous assumptions on behaviour.

    Warns if dimensions are set when not needed

    See test_nufeb_spatial_analysis.py for example usages.

    :raises Value error if periodicity is not recognized

    :param df: a polars dataframe with columns x, y, and z f
    :param radius: Radius within which to search
    :param periodicity: is the simulation non-periodic or does it wrap around in xy plane
    :param xlen: x-dimension length (required for 'xy' periodicity)
    :param ylen: x-dimension length (required for 'xy' periodicity)
    :param zlen: z-dimension length (not yet required )
    :return: A dictionary of lists. Each key is an ID of a point in the dataframe. The list items are the IDs
    of neighbors within the search radius, accounting for periodicity. The items are sorted in order of increasing
    distance. In the case of matching distance, there is no guarantee of order.
    """
    _validate_periodicity_lens(df, periodicity, xlen, ylen, zlen)

    coords = df.select(['x', 'y', 'z']).to_numpy()
    ids = df.select(['id']).to_numpy()
    match periodicity:
        case Periodicity.NONE:
            tree = KDTree(coords)
            indices, dists = tree.query_radius(coords, r=radius,  return_distance = True, sort_results=True)
            neighbor_lists = {}
            for i, index in enumerate(indices):
                ilist = index.tolist()
                if i in ilist:
                    ilist.remove(i)
                neighbor_lists[int(ids[i][0])] = [int(ids[x][0]) for x in ilist]
            return neighbor_lists
        case Periodicity.XY:
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
            results = [list(dict.fromkeys((arr % coords.shape[0]).astype(int))) for arr in indices]
            neighbor_lists = {}
            for i, index in enumerate(results):
                if i in index:
                    index.remove(i)
                neighbor_lists[int(ids[i][0])] = [int(ids[x][0]) for x in index]
            return neighbor_lists
        case _:
            raise ValueError(
                f'Unrecognized periodicity: {periodicity}. Must be one of {",".join([p.value for p in Periodicity])}')

def local_population_structure(df: pl.DataFrame, radius: float, periodicity: Periodicity,
                     xlen: float=None, ylen:float=None, zlen:float=None) -> pl.DataFrame:
    """
    Determine the population structure around a bug.
    :param df: a polars dataframe with columns x, y, and z f
    :param radius: Radius within which to search
    :param periodicity: is the simulation non-periodic or does it wrap around in xy plane
    :param xlen: x-dimension length (required for 'xy' periodicity)
    :param ylen: x-dimension length (required for 'xy' periodicity)
    :param zlen: z-dimension length (not yet required )
    :return: A polars dataframe with n+1 columns. One column lists the ID of a bug. There other n columns are for each
    bug type (group).  The values in each column are the counts of bugs of that type within the radius of bug ID (or 0
    if none). There is no guarantee of order.
    """
    # TODO update this to use the periodicity enum
    neighbor_ids = neighbors_radius(df, radius, periodicity, xlen, ylen, zlen)
    rows = [(k, n) for k, vals in neighbor_ids.items() for n in vals]
    neighbor_df = pl.DataFrame(rows, schema=["id", "neighbor_id"])
    # TODO join this all up and do a lazy_eval
    neighbor_df = neighbor_df.join(
        df.select(["id", "group"]),
        left_on="neighbor_id",
        right_on="id",
        how="left"
    ).rename({"group": "neighbor_group"})

    counts = (
        neighbor_df
        .group_by(["id", "neighbor_group"])
        .count()
        .rename({"count": "n"})
    )

    wide = counts.pivot(
        values="n",
        index="id",
        columns="neighbor_group"
    ).fill_null(0).sort('id')
    return wide

def _validate_periodicity_lens(df: pl.DataFrame, periodicity: Periodicity, xlen: float=None, ylen:float=None, zlen:float=None) -> None:
    """
    Internal method which raises errors if Periodicity is unknown or the x, y, or z lengths don't make sense in that context.
    :param df: Dataframe with id, x, y, z columns representing bugs at a timepoint
    :param periodicity: is the simulation non-periodic or does it wrap around in xy plane
    :param xlen: simulation size in xdim
    :param ylen:  simulation size in ydim
    :param zlen:  simulation size in zdim
    :return: None
    """
    match periodicity:
        case Periodicity.NONE:
            if xlen is not None or ylen is not None or zlen is not None:
                warnings.warn(
                    "Either xlen, ylen, or zlen is set but is not needed for no periodicity. Are you sure you're asking for what you're expecting?",
                    UserWarning)
        case Periodicity.XY:
            if xlen is None and ylen is None:
                raise ValueError('Periodicity of "xy" specified but xlen and ylen are not set')
            if xlen is None:
                raise ValueError('Periodicity of "xy" specified but xlen is not set')
            if ylen is None:
                raise ValueError('Periodicity of "xy" specified but ylen is not set')
            if xlen <= 0 and ylen <= 0:
                raise ValueError(f'Periodicity of "xy" specified but xlen and ylen not > 0. xlen: {xlen}, ylen: {ylen}')
            if xlen <= 0:
                raise ValueError(f'Periodicity of "xy" specified but xlen is not > 0. xlen: {xlen}')
            if ylen <= 0:
                raise ValueError(f'Periodicity of "xy" specified but ylen is not > 0. ylen: {ylen}')
            max_x = df.select(['x']).max().item()
            max_y = df.select(['y']).max().item()
            if xlen < max_x and ylen < max_y:
                raise ValueError(f'xlen, ylen are {xlen}, {ylen}, lower than max values in dataset:{max_x} {max_y}')
            if xlen < max_x:
                raise ValueError(f'xlen is specified to {xlen}, lower than max x-value of points in dataset: {max_x}')
            if ylen < max_y:
                raise ValueError(f'ylen is specified to {ylen}, lower than max y-value of points in dataset: {max_y}')
        case _:
            raise ValueError(f'Unrecognized periodicity: {periodicity}. Must be one of {",".join([p.value for p in Periodicity])}')


def distance_to_each_group(df: pl.DataFrame, periodicity: Periodicity, xlen: float=None, ylen:float=None, zlen:float=None) -> pl.DataFrame:
    """
    Determines, for each bug, the distance to the nearest bug of each group in the simulation at that timepoint.
    Bugs cannot be closest to themselves UNLESS they are the only bug in that group.

    :param df: Dataframe with id, x, y, z columns representing bugs at a timepoint
    :param periodicity: is the simulation non-periodic or does it wrap around in xy plane
    :param xlen: simulation size in xdim
    :param ylen:  simulation size in ydim
    :param zlen:  simulation size in zdim
    :return: A dataframe of the form id, type-1-dist, type-1-id, type-2-dist, type-2-id, ... type-n-dist, type-n-id
    """
    _validate_periodicity_lens(df, periodicity, xlen, ylen, zlen)
    # TODO DRY out the common stuff regarding distances
    match periodicity:
        case Periodicity.NONE:
            # for bookkeeping and individual operations
            coords = df.select(['x', 'y', 'z']).to_numpy()
            ids = df.select(['id']).to_numpy().flatten()
            groups = df.select(['group']).to_numpy().flatten()

            # Get squared periodic distances
            periodic_lengths = np.array([xlen, ylen])
            ## Non-periodic differences
            diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
            ## avoiding sqrt on N*N since it's not necessary
            # D = np.sqrt((diff ** 2).sum(axis=2))
            dist_sq = (diff ** 2).sum(axis=2)

            # Filter based on group and min squared distance
            unique_types = np.unique(groups)
            N = len(ids)
            num_types = len(unique_types)
            columns = {"id": ids}
            for t in unique_types:
                mask = np.array(groups).flatten() == t
                sq_dist_masked = np.where(mask[np.newaxis, :], dist_sq, np.inf)

                if mask.sum() == 1:
                    # handle the case where a bug is the only one of its type
                    only_idx = np.where(mask)[0][0]
                    closest_idx = np.full(N, only_idx)
                    closest_sq = sq_dist_masked[:, only_idx]
                    closest_sq[only_idx] = 0.0
                else:
                    # avoid self distance for relevant mask
                    np.fill_diagonal(sq_dist_masked, np.inf)
                    # get index and calc sqrt for only nearest
                    closest_idx = sq_dist_masked.argmin(axis=1)
                    closest_sq = sq_dist_masked[np.arange(N), closest_idx]

                # update dict used to create dataframe
                columns[f"type-{t}-dist"] = np.sqrt(closest_sq)
                columns[f"type-{t}-id"] = ids[closest_idx]
            return pl.DataFrame(columns)
        case Periodicity.XY:
            # for bookkeeping and individual operations
            coords = df.select(['x', 'y', 'z']).to_numpy()
            ids = df.select(['id']).to_numpy().flatten()
            groups = df.select(['group']).to_numpy().flatten()

            # Get squared periodic distances
            periodic_lengths = np.array([xlen, ylen])
            ## Non-periodic differences
            diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
            diff[:, :, :2] = np.abs(diff[:, :, :2])
            ## adjust x-y for periodicity
            diff[:, :, :2] = np.minimum(diff[:, :, :2], periodic_lengths - diff[:, :, :2])
            ## avoiding sqrt on N*N since it's not necessary
            #D = np.sqrt((diff ** 2).sum(axis=2))
            dist_sq = (diff ** 2).sum(axis=2)

            # Filter based on group and min squared distance
            unique_types = np.unique(groups)
            N = len(ids)
            num_types = len(unique_types)
            columns = {"id": ids}
            for t in unique_types:
                mask = np.array(groups).flatten() == t
                sq_dist_masked = np.where(mask[np.newaxis, :], dist_sq, np.inf)

                if mask.sum() == 1:
                    # handle the case where a bug is the only one of its type
                    only_idx = np.where(mask)[0][0]
                    closest_idx = np.full(N, only_idx)
                    closest_sq = sq_dist_masked[:, only_idx]
                    closest_sq[only_idx] = 0.0
                else:
                    # avoid self distance for relevant mask
                    np.fill_diagonal(sq_dist_masked, np.inf)
                    # get index and calc sqrt for only nearest
                    closest_idx = sq_dist_masked.argmin(axis=1)
                    closest_sq = sq_dist_masked[np.arange(N), closest_idx]

                # update dict used to create dataframe
                columns[f"type-{t}-dist"] = np.sqrt(closest_sq)
                columns[f"type-{t}-id"] = ids[closest_idx]
            return pl.DataFrame(columns)
        case _:
            # Probably caught by call to _validate_periodicity_lens above, but playing safe
            raise ValueError(f'Unrecognized periodicity: {periodicity}. Must be one of {",".join([p.value for p in Periodicity])}')
    # tile if periodic
    # get index, distances
    # map index to id
    # join with group codes
    # filter to get min for each group
    pass