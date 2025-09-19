"""
Calculate various spatial relationships between bugs in simulation results.

This module provides functions which, lists all the neighbors within a given radius,
lists a population summary (e.g., number of neighbors in each group-type), and
lists the nearest neighbor of each type.
"""

import warnings
from enum import Enum

import numpy as np
import polars as pl
from sklearn.neighbors import KDTree


class Periodicity(Enum):
    """
    Boundary conditions for simulation regarding how bugs and substances 'wrap around'.

    Currently supports 'none' and 'x-y' wrapround.

    Attributes:
        XY: Bugs and substances wrap around the sides, but not through the top and
            bottom (height direction)
        NONE: Boundaries act like walls

    """

    XY = "xy"
    NONE = "none"


def neighbors_radius(
    df: pl.DataFrame,
    radius: float,
    periodicity: Periodicity,
    xlen: float = None,
    ylen: float = None,
    zlen: float = None,
    xbleed: float = 0.0,
    ybleed: float = 0.0,
) -> dict[int, list[int]]:
    """
    Find all neighbours for each point within a radius.

    Points are specified as their index into the dataframe.to_numpy()
    Periodicity must be specified and currently supports "none" or "xy".
    It is intentionally not a default parameter to prevent erroneous assumptions on
    behaviour.

    Warns if dimensions are set when not needed

    See test_nufeb_spatial_analysis.py for example usages.

    :raises Value error if periodicity is not recognized

    :param df: a polars dataframe with columns x, y, and z f
    :param radius: Radius within which to search
    :param periodicity: is the simulation non-periodic or does it wrap xy plane
    :param xlen: x-dimension length (required for 'xy' periodicity)
    :param ylen: x-dimension length (required for 'xy' periodicity)
    :param zlen: z-dimension length (not yet required )
    :param xbleed: allow bugs to be slightly outside bounds for len check
    :param ybleed: allow bugs to be slightly outside bounds for len check
    :return: A dictionary of lists. Each key is an ID of a point in the dataframe.
    The list items are the IDs of neighbors within the search radius, accounting for
    periodicity. The items are sorted in order of increasing distance. In the case of
    matching distance, there is no guarantee of order.
    """
    _validate_periodicity_lens(df, periodicity, xlen, ylen, zlen, xbleed, ybleed)
    if radius <= 0:
        raise (ValueError(f"Radius is {radius}, but must be greater than 0"))
    coords = df.select(["x", "y", "z"]).to_numpy()
    ids = df.select(["id"]).to_numpy()
    match periodicity:
        case Periodicity.NONE:
            tree = KDTree(coords)
            indices, dists = tree.query_radius(
                coords, r=radius, return_distance=True, sort_results=True
            )
            neighbor_lists = {}
            for i, index in enumerate(indices):
                ilist = index.tolist()
                if i in ilist:
                    ilist.remove(i)
                neighbor_lists[int(ids[i][0])] = [int(ids[x][0]) for x in ilist]
            return neighbor_lists
        case Periodicity.XY:
            offsets = np.array(
                [
                    [0, 0],
                    [xlen, 0],
                    [-xlen, 0],
                    [0, ylen],
                    [0, -ylen],
                    [xlen, ylen],
                    [xlen, -ylen],
                    [-xlen, ylen],
                    [-xlen, -ylen],
                ]
            )
            tiled_coords_list = []
            for dx, dy in offsets:
                # Shift x and y by tile offsets, leave z unchanged
                shifted = coords.copy()
                shifted[:, 0] += dx
                shifted[:, 1] += dy
                tiled_coords_list.append(shifted)
            tiled_coords = np.vstack(tiled_coords_list)
            tree = KDTree(tiled_coords)
            indices, dists = tree.query_radius(
                coords, r=radius, return_distance=True, sort_results=True
            )
            results = [
                list(dict.fromkeys((arr % coords.shape[0]).astype(int)))
                for arr in indices
            ]
            neighbor_lists = {}
            for i, index in enumerate(results):
                if i in index:
                    index.remove(i)
                neighbor_lists[int(ids[i][0])] = [int(ids[x][0]) for x in index]
            return neighbor_lists
        case _:
            raise ValueError(
                f"Unrecognized periodicity: {periodicity}. "
                f"Must be one of {','.join([p.value for p in Periodicity])}"
            )


def local_population_structure(
    df: pl.DataFrame,
    radius: float,
    periodicity: Periodicity,
    xlen: float = None,
    ylen: float = None,
    zlen: float = None,
    xbleed: float = 0.0,
    ybleed: float = 0.0,
) -> pl.DataFrame:
    """
    Determine the population structure around a bug.

    :param df: a polars dataframe with columns x, y, and z
    :param radius: Radius within which to search
    :param periodicity: is the simulation non-periodic or does it wrap around xy plane
    :param xlen: x-dimension length (required for 'xy' periodicity)
    :param ylen: x-dimension length (required for 'xy' periodicity)
    :param zlen: z-dimension length (not yet required )
    :param xbleed: allow bugs to be slightly outside bounds for len check
    :param ybleed: allow bugs to be slightly outside bounds for len check
    :return: A polars dataframe with n+1 columns. One column lists the ID of a bug.
    There other n columns are for each bug type (group).  The values in each column are
    the counts of bugs of that type within the radius of bug ID (or 0 if none). There is
    no guarantee of order. If a bug has no neighbors, it returns 0 for every group
    """
    neighbor_ids = neighbors_radius(
        df,
        radius,
        periodicity,
        xlen, ylen, zlen,
        xbleed, ybleed
    )
    rows = [(k, n) for k, vals in neighbor_ids.items() for n in vals]
    if rows != []:
        neighbor_df = pl.DataFrame(rows, schema=["id", "neighbor_id"], orient='row')
        # ensure all ids appear
        all_ids = pl.DataFrame({"id": list(neighbor_ids.keys())})
        neighbor_df = all_ids.join(neighbor_df, on="id", how="left")
    else:
        # deal with the case where nobody has neighbors
        neighbor_df = pl.DataFrame({"id": list(neighbor_ids.keys())}).with_columns(
            pl.lit(None, pl.UInt32).alias("neighbor_id")
        )
    neighbor_df_wide = (
        neighbor_df.lazy()
        .join(
            df.lazy().select(["id", "group"]),
            left_on="neighbor_id", right_on="id", how="left"
        )
        .rename({"group": "neighbor_group"})
        .group_by(["id", "neighbor_group"]).count().rename({"count": "n"})
        # before lazy, did a pivot. Pivot doesn't work well with lazy frames,
        # so replaced with group_by->agg
        #.pivot(values="n", index="id", columns="neighbor_group")
        .group_by('id')
        .agg([
            pl.when(pl.col("neighbor_group") == g)
                .then(pl.col("n"))
                .sum()
                .fill_null(0)
                .alias(str(g))
           for g in df.select("group").unique().to_series().to_list()
        ])
        .sort("id")
        .drop("null", strict=False)
    )

    # append any groups which exist but which were not in any neighbor list as all 0's
    all_groups = df["group"].unique().cast(pl.Utf8).to_list()
    for g in all_groups:
        if g not in neighbor_df_wide.columns:
            neighbor_df_wide = neighbor_df_wide.with_columns(pl.lit(0, dtype=pl.UInt32).alias(g))
    return neighbor_df_wide.collect()


def _validate_periodicity_lens(
    df: pl.DataFrame,
    periodicity: Periodicity,
    xlen: float = None,
    ylen: float = None,
    zlen: float = None,
    xbleed: float = 0.0,
    ybleed: float = 0.0,
) -> None:
    """
    Validate x, y, z lengths based on periodicity and data values.

    Raise errors if Periodicity is unknown or the x, y, or z lengths don't make sense
    in that context.

    Note that it's an internal method, as this validation is common to multiple funcs
    :param df: Dataframe with id, x, y, z columns representing bugs at a timepoint
    :param periodicity: is simulation non-periodic or does it wrap around in xy plane
    :param xlen: simulation size in xdim
    :param ylen:  simulation size in ydim
    :param zlen:  simulation size in zdim
    :param xbleed: allow bugs to be slightly outside bounds for len check
    :param ybleed: allow bugs to be slightly outside bounds for len check
    :return: None
    """
    match periodicity:
        case Periodicity.NONE:
            if zlen is not None:
                pv = periodicity.value
                warnings.warn(
                    f"zlen is set but is not needed for Periodicity:{pv}. "
                    f"Are you sure you're asking for what you're expecting?",
                    UserWarning,
                    stacklevel=2
                )
            if xlen is not None or ylen is not None:
                pv = periodicity.value
                warnings.warn(
                    f"Either xlen or ylen is set but is not needed for Periodicity:{pv}"
                    f". Are you sure you're asking for what you're expecting?",
                    UserWarning,
                    stacklevel=2
                )
        case Periodicity.XY:
            if zlen is not None:
                pv = periodicity.value
                warnings.warn(
                    f"zlen is set but is not needed for Periodicity:{pv}."
                    f"Are you sure you're asking for what you're expecting?",
                    UserWarning,
                    stacklevel=2
                )
            if xlen is None and ylen is None:
                raise ValueError(
                    'Periodicity of "xy" specified but xlen and ylen are not set'
                )
            if xlen is None:
                raise ValueError('Periodicity of "xy" specified but xlen is not set')
            if ylen is None:
                raise ValueError('Periodicity of "xy" specified but ylen is not set')
            if xlen <= 0 and ylen <= 0:
                raise ValueError(
                    f'Periodicity of "xy" specified but xlen and ylen not > 0. '
                    f'xlen: {xlen}, ylen: {ylen}'
                )
            if xlen <= 0:
                raise ValueError(
                    f'Periodicity of "xy" specified but xlen is not > 0. xlen: {xlen}'
                )
            if ylen <= 0:
                raise ValueError(
                    f'Periodicity of "xy" specified but ylen is not > 0. ylen: {ylen}'
                )

            max_x = df.select(["x"]).max().item()
            max_y = df.select(["y"]).max().item()
            # bleeds allow for situations where LAMMPS has briefly allowed
            # a bug to be slightly out of bounds
            # build partial string for better error message
            xbleedstr = ''
            if xbleed > 0:
                xbleedstr = f' with bleed of {xbleed}'
            ybleedstr = ''
            if ybleed > 0:
                ybleedstr = f' with bleed of {ybleed}'
            xybleedstr = ''
            if xbleed > 0 and ybleed > 0:
                xybleedstr = f' with bleeds of {xbleed}, {ybleed}'
            if xlen < (max_x - xbleed) and ylen < (max_y - ybleed):
                raise ValueError(
                    f"xlen, ylen are {xlen}, {ylen}, "
                    f"lower than max values in dataset:{max_x} {max_y}{xybleedstr}"
                )
            if xlen < max_x - xbleed:
                raise ValueError(
                    f"xlen is specified to {xlen}, "
                    f"lower than max x-value of points in dataset: {max_x}{xbleedstr}"
                )
            if ylen < max_y - ybleed:
                raise ValueError(
                    f"ylen is specified to {ylen}, "
                    f"lower than max y-value of points in dataset: {max_y}{ybleedstr}"
                )
        case _:
            raise ValueError(
                f"Unrecognized periodicity: {periodicity}."
                f"Must be one of {','.join([p.value for p in Periodicity])}"
            )


def distance_to_each_group(
    df: pl.DataFrame,
    periodicity: Periodicity,
    xlen: float = None,
    ylen: float = None,
    zlen: float = None,
    xbleed: float = 0.0,
    ybleed: float = 0.0,
) -> pl.DataFrame:
    """
    List distances to nearest bugs based on group.

    Determine, for each bug, the distance to the nearest bug of each group in the
    simulation at that timepoint. Bugs cannot be closest to themselves UNLESS they are
    the only bug in that group.

    :param df: Dataframe with id, x, y, z columns representing bugs at a timepoint
    :param periodicity: non-periodic or does it wrap around in xy plane
    :param xlen: simulation size in xdim
    :param ylen:  simulation size in ydim
    :param zlen:  simulation size in zdim
    :param xbleed: allow bugs to be slightly outside bounds for len check
    :param ybleed: allow bugs to be slightly outside bounds for len check
    :return: A dataframe of the form id, type-1-dist, type-1-id, type-2-dist, type-2-id,
             ... type-n-dist, type-n-id
    """
    _validate_periodicity_lens(df, periodicity, xlen, ylen, zlen, xbleed, ybleed)
    # TODO DRY out the common stuff regarding distances
    match periodicity:
        case Periodicity.NONE:
            # for bookkeeping and individual operations
            coords = df.select(["x", "y", "z"]).to_numpy()
            ids = df.select(["id"]).to_numpy().flatten()
            groups = df.select(["group"]).to_numpy().flatten()

            # Get squared periodic distances
            periodic_lengths = np.array([xlen, ylen])
            ## Non-periodic differences
            diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
            ## avoiding sqrt on N*N since it's not necessary
            # D = np.sqrt((diff ** 2).sum(axis=2))
            dist_sq = (diff**2).sum(axis=2)

            # Filter based on group and min squared distance
            unique_types = np.unique(groups)
            n_ids = len(ids)
            columns = {"id": ids}
            for t in unique_types:
                mask = np.array(groups).flatten() == t
                sq_dist_masked = np.where(mask[np.newaxis, :], dist_sq, np.inf)

                if mask.sum() == 1:
                    # handle the case where a bug is the only one of its type
                    only_idx = np.where(mask)[0][0]
                    closest_idx = np.full(n_ids, only_idx)
                    closest_sq = sq_dist_masked[:, only_idx]
                    closest_sq[only_idx] = 0.0
                else:
                    # avoid self distance for relevant mask
                    np.fill_diagonal(sq_dist_masked, np.inf)
                    # get index and calc sqrt for only nearest
                    closest_idx = sq_dist_masked.argmin(axis=1)
                    closest_sq = sq_dist_masked[np.arange(n_ids), closest_idx]

                # update dict used to create dataframe
                columns[f"type-{t}-dist"] = np.sqrt(closest_sq)
                columns[f"type-{t}-id"] = ids[closest_idx]
            return pl.DataFrame(columns)
        case Periodicity.XY:
            # for bookkeeping and individual operations
            coords = df.select(["x", "y", "z"]).to_numpy()
            ids = df.select(["id"]).to_numpy().flatten()
            groups = df.select(["group"]).to_numpy().flatten()

            # Get squared periodic distances
            periodic_lengths = np.array([xlen, ylen])
            ## Non-periodic differences
            diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
            diff[:, :, :2] = np.abs(diff[:, :, :2])
            ## adjust x-y for periodicity
            diff[:, :, :2] = np.minimum(
                diff[:, :, :2], periodic_lengths - diff[:, :, :2]
            )
            ## avoiding sqrt on N*N since it's not necessary
            # D = np.sqrt((diff ** 2).sum(axis=2))
            dist_sq = (diff**2).sum(axis=2)

            # Filter based on group and min squared distance
            unique_types = np.unique(groups)
            n_ids = len(ids)
            columns = {"id": ids}
            for t in unique_types:
                mask = np.array(groups).flatten() == t
                sq_dist_masked = np.where(mask[np.newaxis, :], dist_sq, np.inf)

                if mask.sum() == 1:
                    # handle the case where a bug is the only one of its type
                    only_idx = np.where(mask)[0][0]
                    closest_idx = np.full(n_ids, only_idx)
                    closest_sq = sq_dist_masked[:, only_idx]
                    closest_sq[only_idx] = 0.0
                else:
                    # avoid self distance for relevant mask
                    np.fill_diagonal(sq_dist_masked, np.inf)
                    # get index and calc sqrt for only nearest
                    closest_idx = sq_dist_masked.argmin(axis=1)
                    closest_sq = sq_dist_masked[np.arange(n_ids), closest_idx]

                # update dict used to create dataframe
                columns[f"type-{t}-dist"] = np.sqrt(closest_sq)
                columns[f"type-{t}-id"] = ids[closest_idx]
            return pl.DataFrame(columns)
        case _:
            # Probably caught by _validate_periodicity_lens above, but playing safe
            raise ValueError(
                f"Unrecognized periodicity: {periodicity}."
                f"Must be one of {','.join([p.value for p in Periodicity])}"
            )

