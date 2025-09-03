import pytest
import polars as pl
from polars.testing import assert_frame_equal
from nufebmgr import nufeb_spatial_analysis as nu_spa
from nufebmgr.nufeb_spatial_analysis import Periodicity

# Example dataset
df_simple = pl.DataFrame({
    "id": [1, 2, 3, 4],
    "group": [1, 2, 1, 3],
    "x": [0.0, 1.0, 3.0, 0.0],
    "y": [0.0, 1.0, 2.0, 4.0],
    "z": [0.0, 0.0, 0.0, 5.0],
})

# for easy reference, within df_simple, distances, by ID:
    # 1 to 2: 1.41421356
    # 1 to 3: 3.6055128
    # 1 to 4: 6.40312424
    # 2 to 3: 2.23606798
    # 2 to 4: 5.91607978
    # 3 to 4: 6.164414

df_neighbor_periodicity = pl.DataFrame({
    'id': [1, 2, 3, 4, 5],
    'group': [1, 1, 2, 2, 3],
    'x': [0.1, 4.9, 2.5, 0.2, 4.8],
    'y': [0.1, 0.2, 2.5, 4.8, 7.9],
    'z': [0.0, 0.0, 1.0, 0.5, 0.5]})

df_neighbor_periodicity_non_contig_id = pl.DataFrame({
    'id': [14, 12, 2, 25, 4],
    'group': [1, 1, 2, 2, 3],
    'x': [0.1, 4.9, 2.5, 0.2, 4.8],
    'y': [0.1, 0.2, 2.5, 4.8, 7.9],
    'z': [0.0, 0.0, 1.0, 0.5, 0.5]})
# from scipy.spatial import distance_matrix
# coords = df_neighbor_periodicity.select(['x','y','z']).to_numpy()
# dist_mat = distance_matrix(coords, coords)
# non-periodic distance matrix for above
# array([[0.        , 4.80104155, 3.5383612 , 4.72757866, 9.12030701],
#        [4.80104155, 0.        , 3.47131099, 6.59545298, 7.71686465],
#        [3.5383612 , 3.47131099, 0.        , 3.29089653, 5.89067059],
#        [4.72757866, 6.59545298, 3.29089653, 0.        , 5.5470713 ],
#        [9.12030701, 7.71686465, 5.89067059, 5.5470713 , 0.        ]])

# more points in site. intended to be x=6, y=5 distance units
# intentionally set up to exercise boundary conditions in xy, especially along diagonals. note all z=0
df_neighbor_periodicity_two = pl.DataFrame({
    'id': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
    'group': [1, 1, 2, 2, 3, 1, 2, 2, 1, 3, 1, 1],
    'x': [0.5, 1.0, 5.5, 0.5, 2.0, 2.0, 5.5, 2.0, 3.0, 0.5, 2.0, 5.5],
    'y': [0.5, 0.5, 0.5, 1.0, 1.0, 2.0, 2.0, 2.5, 2.5, 4.5, 4.5, 4.5],
    'z': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]})


def test_neighbors_radius_non_periodic():
    # array([[0.        , 4.80104155, 3.5383612 , 4.72757866, 9.12030701],
    #        [4.80104155, 0.        , 3.47131099, 6.59545298, 7.71686465],
    #        [3.5383612 , 3.47131099, 0.        , 3.29089653, 5.89067059],
    #        [4.72757866, 6.59545298, 3.29089653, 0.        , 5.5470713 ],
    #        [9.12030701, 7.71686465, 5.89067059, 5.5470713 , 0.        ]])
    neighbour_lists = nu_spa.neighbors_radius(df_neighbor_periodicity, radius=6, periodicity=Periodicity.NONE)
    expected = { 1: [3, 4, 2],
                 2: [3, 1],
                 3: [4, 2, 1, 5],
                 4: [3, 1, 5],
                 5: [4, 3]}
    assert neighbour_lists == expected

    neighbour_lists = nu_spa.neighbors_radius(df_neighbor_periodicity_two, radius=3, periodicity=Periodicity.NONE)
    expected = {1: [2, 4, 5, 6, 8],
                2: [1, 4, 5, 6, 8, 9],
                3: [7],
                4: [1, 2, 5, 6, 8, 9],
                5: [6, 2, 8, 4, 1, 9],
                6: [8, 5, 9, 2, 4, 1, 11, 10],
                7: [3, 12, 9],
                8: [6, 9, 5, 11, 4, 2, 1, 10],
                9: [8, 6, 5, 11, 7, 2, 4],
                10: [11, 8, 6],
                11: [10, 8, 9, 6],
                12: [7]
                }
    assert neighbour_lists == expected

    # using non-contiguous, out of numeric order id column to check if ids are correctly preserved
    # For ref:      'id': [14, 12, 2, 25, 4],
    neighbour_lists = nu_spa.neighbors_radius(df_neighbor_periodicity_non_contig_id, radius=6, periodicity=Periodicity.NONE)
    expected = {14: [2, 25, 12],
                12: [2, 14],
                2: [25, 12, 14, 4],
                25: [2, 14, 4],
                4: [25, 2]}
    assert neighbour_lists == expected


def test_neighbors_radius_periodic_xy():
    # distance matrix for x=5, y=9. Lower triangle. truncated without rounding
    #   0           1           2           3           4
    # 0 0
    # 1 0.2236067   0
    # 2 3.5382612   3.471310    0
    # 3 4.3301270   4.438468    3.2908965   0
    # 4 1.3341664   1.396424    4.3011626   3.1256999   0

    neighbour_lists = nu_spa.neighbors_radius(df_neighbor_periodicity, radius=3.5, periodicity=Periodicity.XY, xlen=5, ylen=9)
    expected = { 1: [2, 5],
                 2: [1, 5, 3],
                 3: [4, 2],
                 4: [5, 3],
                 5: [1, 2, 4]}
    assert neighbour_lists == expected

    # tighter radius
    neighbour_lists = nu_spa.neighbors_radius(df_neighbor_periodicity, radius=2, periodicity=Periodicity.XY, xlen=5, ylen=9)
    expected = { 1: [2, 5],
                 2: [1, 5],
                 3: [],
                 4: [],
                 5: [1, 2]}
    assert neighbour_lists == expected

    # using non-contiguous, out of numeric order id column to check if ids are correctly preserved
    # For ref:      'id': [14, 12, 2, 25, 4],
    neighbour_lists = nu_spa.neighbors_radius(df_neighbor_periodicity_non_contig_id, radius=2, periodicity=Periodicity.XY, xlen=5, ylen=9)
    expected = { 14: [12, 4],
                 12: [14, 4],
                 2: [],
                 25: [],
                 4: [14, 12]}
    assert neighbour_lists == expected

    # different setup, just to be extra careful
    neighbour_lists = nu_spa.neighbors_radius(df_neighbor_periodicity_two, radius=3, periodicity=Periodicity.XY, xlen=6, ylen=5)
    expected = {1: [4, 2, 3, 10, 12, 5, 11, 7, 6, 8],
                2: [1, 4, 10, 5, 11, 3, 6, 12, 7, 8, 9],
                3: [1, 12, 4, 10, 7, 2, 5, 11, 6],
                4: [1, 2, 3, 7, 10, 5, 12, 6, 8, 11, 9],
                5: [6, 2, 4, 8, 11, 1, 9, 10, 3, 7, 12],
                6: [8, 5, 9, 2, 4, 1, 11, 7, 10, 3],
                7: [4, 3, 1, 2, 12, 6, 8, 9, 5, 10],
                8: [6, 9, 5, 11, 4, 2, 1, 10, 7],
                9: [8, 6, 5, 11, 7, 2, 4],
                10: [1, 12, 2, 3, 4, 11, 5, 8, 7, 6],
                11: [2, 5, 10, 1, 8, 4, 9, 6, 12, 3],
                12: [10, 3, 1, 4, 2, 11, 7, 5]}

    assert neighbour_lists == expected


def test_neighbors_and_pop_struct_radius_edge():
    df = pl.DataFrame({
        "id": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        "group": [2, 2, 1, 3, 1, 3, 2, 2, 1, 3],
        "x": [0.0, 0.5, 0.2, 0.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        "y": [0.0, 0.0, 0.0, 0.0, 0.2, 0.8, 0.0, 0.0, 0.5, 0.0],
        "z": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.8, 0.0, 0.5],
    })
    radius = 0.2
    ## XY periodicity neighbours
    expected_periodic_xy_neighbors = {1: [4, 6, 5, 7, 3],
                2: [],
                3: [1],
                4: [1],
                5: [1],
                6: [1],
                7: [1],
                8: [],
                9: [],
                10: []}
    results_periodic_xy_neighbors = nu_spa.neighbors_radius(df, radius, Periodicity.XY, xlen=1.0, ylen=1.0)
    assert results_periodic_xy_neighbors == expected_periodic_xy_neighbors

    ## NONE periodicity neighbors
    expected_periodic_none_neighbors = {1: [7, 5, 3],
                                        2: [],
                                        3: [1],
                                        4: [],
                                        5: [1],
                                        6: [],
                                        7: [1],
                                        8: [],
                                        9: [],
                                        10: []}
    results_periodic_none_neighbors = nu_spa.neighbors_radius(df, radius, Periodicity.NONE, xlen=1.0, ylen=1.0)
    assert results_periodic_none_neighbors == expected_periodic_none_neighbors

    # XY periodicity local pop struct
    expected_XY_groups = pl.DataFrame({
               'id': pl.Series([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], dtype=pl.Int64),
               '1': pl.Series([2, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], dtype=pl.UInt32),
               '2': pl.Series([1, 0, 1, 1, 1, 1, 1, 0, 0, 0], dtype=pl.UInt32),
               '3': pl.Series([2, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype=pl.UInt32)})

    results_XY_groups = nu_spa.local_population_structure(df, radius, Periodicity.XY, xlen=1.0, ylen=1.0)
    assert_frame_equal(results_XY_groups, expected_XY_groups, check_column_order=False)

    # NONE periodicity local pop struct
    expected_none_groups = pl.DataFrame({
               'id': pl.Series([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], dtype=pl.Int64),
               '1': pl.Series([2, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], dtype=pl.UInt32),
               '2': pl.Series([1, 0, 1, 0, 1, 0, 1, 0, 0, 0], dtype=pl.UInt32),
               '3': pl.Series([0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype=pl.UInt32)})

    results_none_groups = nu_spa.local_population_structure(df, radius, Periodicity.NONE, xlen=1.0, ylen=1.0)
    assert_frame_equal(results_none_groups, expected_none_groups, check_column_order=False)


def test_local_pop_structure_one_group_never_in():
    # sometimes a group may not be in any local population structures, make sure it still shows as a column of 0's
    # this was surfaced as a coincidental bug in other testing. The root cause was that the group list was
    # originally built by looking at neighbours. Now, we also cheaply get all the groups and append empty columns if
    # necessary. This test is to guard against regression
    df = pl.DataFrame({
        "id": [1, 2, 3],
        "group": [1, 2, 3],
        "x": [0.0, 0.5, 0.2],
        "y": [0.0, 0.0, 0.0],
        "z": [0.0, 0.0, 0.0],
    })
    radius = 0.2
    expected = pl.DataFrame({
        'id': pl.Series([1, 2, 3], dtype=pl.Int64),
        '1': pl.Series([0, 0, 1], dtype=pl.UInt32),
        '2': pl.Series([0, 0, 0], dtype=pl.UInt32),
        '3': pl.Series([1, 0, 0], dtype=pl.UInt32)})

    results = nu_spa.local_population_structure(df, radius, Periodicity.NONE, xlen=1.0, ylen=1.0)
    assert_frame_equal(results, expected, check_column_order=False)


def test_point_on_boundary_periodic_none():
    df = pl.DataFrame({
        "id": [1, 2, 3, 4, 5, 6, 7, 8, 9],
        "group": [2, 2, 1, 3, 1, 2, 2, 3, 1],
        "x": [0.0, 1.0, 0.0, 1.0, 0.5, 0.2, 0.2, 0.8, 0.8],
        "y": [0.0, 0.0, 0.5, 0.5, 0.25, 0.0, 0.5, 0.0, 0.5],
        "z": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    })
    radius = 0.2

    ### Neighbors
    expected_neighbors_none = {1: [6],
                             2: [8],
                             3: [7],
                             4: [9],
                             5: [],
                             6: [1],
                             7: [3],
                             8: [2],
                             9: [4]}
    neighbors_none = nu_spa.neighbors_radius(df, radius, periodicity=Periodicity.NONE, xlen=1.0, ylen=0.5)
    assert neighbors_none == expected_neighbors_none

    ### Groups
    expected_none_groups = pl.DataFrame({
        'id': pl.Series([1, 2, 3, 4, 5, 6, 7, 8, 9], dtype=pl.Int64),
        '1': pl.Series([0, 0, 0, 1, 0, 0, 1, 0, 0], dtype=pl.UInt32),
        '2': pl.Series([1, 0, 1, 0, 0, 1, 0, 1, 0], dtype=pl.UInt32),
        '3': pl.Series([0, 1, 0, 0, 0, 0, 0, 0, 1], dtype=pl.UInt32)})
    results_none_groups = nu_spa.local_population_structure(df, radius, periodicity=Periodicity.NONE, xlen=1.0, ylen=0.5)
    assert_frame_equal(results_none_groups, expected_none_groups, check_column_order=False)

    ## nearest in each group
    expected_nearest_none = pl.DataFrame({
        'id': pl.Series([1, 2, 3, 4, 5, 6, 7, 8, 9], dtype=pl.Int64),
        'type-1-dist': pl.Series([.5, .538516, .5590169, .2, .390512, .390512, .2, .390512, .390512], dtype=pl.Float64),
        'type-1-id': pl.Series([3, 9, 5, 9, 9, 5, 3, 5, 5], dtype=pl.Int64),
        'type-2-dist': pl.Series([.2, .8, .2, .5, .390512, .2,.5, .2, .538516], dtype=pl.Float64),
        'type-2-id': pl.Series([6, 6, 7, 2, 6, 1, 6, 2, 2], dtype=pl.Int64),
        'type-3-dist': pl.Series([.8, .2, .943398, .538516, .390512, .6, .781025, .538516, .2], dtype=pl.Float64),
        'type-3-id': pl.Series([8, 8, 8, 8, 8, 8, 8, 4, 4], dtype=pl.Int64),})
    result_nearest_none = nu_spa.distance_to_each_group(df,periodicity=nu_spa.Periodicity.NONE, xlen=1.0, ylen=0.5)
    assert_frame_equal(result_nearest_none, expected_nearest_none, check_column_order=False, check_exact=False)


def test_point_on_boundary_periodic_xy():
    df = pl.DataFrame({
        "id": [1, 2, 3, 4, 5, 6, 7, 8, 9],
        "group": [2, 2, 1, 3, 1, 2, 2, 3, 1],
        "x": [0.0, 1.0, 0.0, 1.0, 0.5, 0.2, 0.2, 0.8, 0.8],
        "y": [0.0, 0.0, 0.5, 0.5, 0.25, 0.0, 0.5, 0.0, 0.5],
        "z": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    })
    radius = 0.2
    # Neighbors
    expected_neighbors_xy = {1: [3, 2, 4, 9, 8, 6, 7],
                             2: [4, 3, 1, 9, 7, 8, 6],
                             3: [2, 4, 1, 9, 8, 6, 7],
                             4: [1, 2, 3, 6, 8, 9, 7],
                             5: [],
                             6: [7, 1, 3, 2, 4],
                             7: [6, 4, 2, 3, 1],
                             8: [9, 2, 1, 4, 3],
                             9: [8, 4, 3, 1, 2]}
    neighbors_xy = nu_spa.neighbors_radius(df, radius, periodicity=Periodicity.XY, xlen=1.0, ylen=0.5)
    assert neighbors_xy == expected_neighbors_xy

    # Groups
    expected_xy_groups = pl.DataFrame({
        'id': pl.Series([1, 2, 3, 4, 5, 6, 7, 8, 9], dtype=pl.Int64),
        '1': pl.Series([2, 2, 1, 2, 0, 1, 1, 2, 1], dtype=pl.UInt32),
        '2': pl.Series([3, 3, 4, 4, 0, 3, 3, 2, 2], dtype=pl.UInt32),
        '3': pl.Series([2, 2, 2, 1, 0, 1, 1, 1, 2], dtype=pl.UInt32)})
    results_xy_groups = nu_spa.local_population_structure(df, radius, periodicity=Periodicity.XY, xlen=1.0, ylen=0.5)
    assert_frame_equal(results_xy_groups, expected_xy_groups, check_column_order=False)

    # Nearest in each group
    expected_nearest_xy = pl.DataFrame({
        'id': pl.Series([1, 2, 3, 4, 5, 6, 7, 8, 9], dtype=pl.Int64),
        'type-1-dist': pl.Series([0, 0, .2, 0, .390512, .2, .2, 0, 0.2], dtype=pl.Float64),
        'type-1-id': pl.Series([3, 3, 9, 3, 9, 3, 3, 9, 3], dtype=pl.Int64),
        'type-2-dist': pl.Series([0, 0, 0, 0, .390512, 0, 0, .2, .2], dtype=pl.Float64),
        'type-2-id': pl.Series([2, 1, 1, 1, 6, 7, 6, 1, 1], dtype=pl.Int64),
        'type-3-dist': pl.Series([0, 0, 0, .2, .390512, .2, .2, .2, 0], dtype=pl.Float64),
        'type-3-id': pl.Series([4, 4, 4,8, 8, 4, 4, 4, 8], dtype=pl.Int64), })
    result_nearest_xy = nu_spa.distance_to_each_group(df,periodicity=nu_spa.Periodicity.XY, xlen=1.0, ylen=0.5)
    assert_frame_equal(result_nearest_xy, expected_nearest_xy, check_column_order=False, check_exact=False)


def test_overlapping_points():
    df = pl.DataFrame({
        "id": [1, 2, 3, 4, 5],
        "group": [1, 2, 2, 3, 4],
        "x": [0.1, 0.5, 0.5, 0.5, 0.6],
        "y": [0.5, 0.5, 0.5, 0.5, 0.5],
        "z": [0.0, 0.0, 0.0, 0.0, 0.0],
    })
    radius = 0.2
    xlen = 1.0
    ylen = 1.0
    # Neighbors
    expected = {1: [],
                2: [4, 3, 5],
                3: [4, 2, 5],
                4: [3, 2, 5],
                5: [2, 4, 3]}
    neighbors = nu_spa.neighbors_radius(df, radius, periodicity=Periodicity.XY, xlen=xlen, ylen=ylen)
    assert neighbors == expected

    # Groups
    expected = pl.DataFrame({
        'id': pl.Series([1, 2, 3, 4, 5], dtype=pl.Int64),
        '1': pl.Series([0, 0, 0, 0, 0], dtype=pl.UInt32),
        '2': pl.Series([0, 1, 1, 2, 2], dtype=pl.UInt32),
        '3': pl.Series([0, 1, 1, 0, 1], dtype=pl.UInt32),
        '4': pl.Series([0, 1, 1, 1, 0], dtype=pl.UInt32)})
    groups = nu_spa.local_population_structure(df, radius, periodicity=Periodicity.XY, xlen=xlen, ylen=ylen)
    assert_frame_equal(groups, expected, check_column_order=False)

    # Nearest in each group
    expected = pl.DataFrame({
        'id': pl.Series([1, 2, 3, 4, 5], dtype=pl.Int64),
        'type-1-dist': pl.Series([0, 0.4, 0.4, 0.4, 0.5], dtype=pl.Float64),
        'type-1-id': pl.Series([1, 1, 1, 1, 1], dtype=pl.Int64),
        'type-2-dist': pl.Series([0.4, 0, 0, 0, 0.1], dtype=pl.Float64),
        'type-2-id': pl.Series([2, 3, 2, 2, 2], dtype=pl.Int64),
        'type-3-dist': pl.Series([0.4, 0, 0, 0, 0.1], dtype=pl.Float64),
        'type-3-id': pl.Series([4, 4, 4, 4, 4], dtype=pl.Int64),
        'type-4-dist': pl.Series([0.5, 0.1, 0.1, 0.1, 0], dtype=pl.Float64),
        'type-4-id': pl.Series([5, 5, 5, 5, 5], dtype=pl.Int64), })
    result_nearest = nu_spa.distance_to_each_group(df, periodicity=nu_spa.Periodicity.XY, xlen=xlen, ylen=ylen)
    assert_frame_equal(result_nearest, expected, check_column_order=False, check_exact=False)


def test_one_point_only():
    df = pl.DataFrame({
        "id": [3],
        "group": [2],
        "x": [0.1],
        "y": [0.3],
        "z": [0.2],
    })
    radius = 0.2
    xlen = 1.0
    ylen = 1.0
    # Neighbors
    expected = {3: []}
    neighbors = nu_spa.neighbors_radius(df, radius, periodicity=Periodicity.XY, xlen=xlen, ylen=ylen)
    assert neighbors == expected

    # Groups
    expected = pl.DataFrame({
        'id': pl.Series([3], dtype=pl.Int64),
        '2': pl.Series([0], dtype=pl.UInt32)})
    groups = nu_spa.local_population_structure(df, radius, periodicity=Periodicity.XY, xlen=xlen, ylen=ylen)
    assert_frame_equal(groups, expected, check_column_order=False)

    # Nearest in each group
    expected = pl.DataFrame({
        'id': pl.Series([3], dtype=pl.Int64),
        'type-2-dist': pl.Series([0], dtype=pl.Float64),
        'type-2-id': pl.Series([3], dtype=pl.Int64), })
    result_nearest = nu_spa.distance_to_each_group(df, periodicity=nu_spa.Periodicity.XY, xlen=xlen, ylen=ylen)
    assert_frame_equal(result_nearest, expected, check_column_order=False, check_exact=False)


def test_z_coordinate():
    df = pl.DataFrame({
        "id": [1,2,3],
        "group": [3,1,2],
        "x": [0.5,0.5,0.5],
        "y": [0.5,0.5,0.5],
        "z": [0.0,0.5,0.6],
    })
    radius = 0.2
    xlen = 1.0
    ylen = 1.0
    # Neighbors
    expected = {1: [],
                2: [3],
                3: [2]}
    neighbors = nu_spa.neighbors_radius(df, radius, periodicity=Periodicity.XY, xlen=xlen, ylen=ylen)
    assert neighbors == expected

    # Groups
    expected = pl.DataFrame({
        'id': pl.Series([1,2,3], dtype=pl.Int64),
        '1': pl.Series([0, 0, 1], dtype=pl.UInt32),
        '2': pl.Series([0, 1, 0], dtype=pl.UInt32),
        '3': pl.Series([0, 0, 0], dtype=pl.UInt32)})
    groups = nu_spa.local_population_structure(df, radius, periodicity=Periodicity.XY, xlen=xlen, ylen=ylen)
    assert_frame_equal(groups, expected, check_column_order=False)

    # Nearest in each group
    expected = pl.DataFrame({
        'id': pl.Series([1,2,3], dtype=pl.Int64),
        'type-1-dist': pl.Series([0.5, 0, 0.1], dtype=pl.Float64),
        'type-1-id': pl.Series([2, 2, 2], dtype=pl.Int64),
        'type-2-dist': pl.Series([0.6, 0.1, 0], dtype=pl.Float64),
        'type-2-id': pl.Series([3, 3, 3], dtype=pl.Int64),
        'type-3-dist': pl.Series([0, 0.5, 0.6], dtype=pl.Float64),
        'type-3-id': pl.Series([1, 1, 1], dtype=pl.Int64), })
    result_nearest = nu_spa.distance_to_each_group(df, periodicity=nu_spa.Periodicity.XY, xlen=xlen, ylen=ylen)
    assert_frame_equal(result_nearest, expected, check_column_order=False, check_exact=False)


def test_no_neighbours_in_radius():
    df = pl.DataFrame({
        "id": [1, 2, 3],
        "group": [3, 1, 2],
        "x": [0.5, 0.5, 0.5],
        "y": [0.5, 0.5, 0.5],
        "z": [0.0, 0.5, 0.6],
    })
    radius = 0.05
    xlen = 1.0
    ylen = 1.0
    # Neighbors
    expected = {1: [],
                2: [],
                3: []}
    neighbors = nu_spa.neighbors_radius(df, radius, periodicity=Periodicity.XY, xlen=xlen, ylen=ylen)
    assert neighbors == expected

    # Groups
    expected = pl.DataFrame({
        'id': pl.Series([1, 2, 3], dtype=pl.Int64),
        '1': pl.Series([0, 0, 0], dtype=pl.UInt32),
        '2': pl.Series([0, 0, 0], dtype=pl.UInt32),
        '3': pl.Series([0, 0, 0], dtype=pl.UInt32)})
    groups = nu_spa.local_population_structure(df, radius, periodicity=Periodicity.XY, xlen=xlen, ylen=ylen)
    assert_frame_equal(groups, expected, check_column_order=False)

    # Nearest in each group
    expected = pl.DataFrame({
        'id': pl.Series([1, 2, 3], dtype=pl.Int64),
        'type-1-dist': pl.Series([0.5, 0, 0.1], dtype=pl.Float64),
        'type-1-id': pl.Series([2, 2, 2], dtype=pl.Int64),
        'type-2-dist': pl.Series([0.6, 0.1, 0], dtype=pl.Float64),
        'type-2-id': pl.Series([3, 3, 3], dtype=pl.Int64),
        'type-3-dist': pl.Series([0, 0.5, 0.6], dtype=pl.Float64),
        'type-3-id': pl.Series([1, 1, 1], dtype=pl.Int64), })
    result_nearest = nu_spa.distance_to_each_group(df, periodicity=nu_spa.Periodicity.XY, xlen=xlen, ylen=ylen)
    assert_frame_equal(result_nearest, expected, check_column_order=False, check_exact=False)


@pytest.mark.parametrize(
    "func", [nu_spa.neighbors_radius, nu_spa.local_population_structure,]
)
@pytest.mark.parametrize(
    "radius, expect_error",
    [ (0, True), (-1, True), (2, False),]
)
def test_radius_lte_zero(func, radius, expect_error):
    if expect_error:
        with pytest.raises(ValueError) as excinfo:
            func(
                df_neighbor_periodicity,
                radius=radius,
                periodicity=Periodicity.NONE,
                xlen=6,
                ylen=5,
            )
        expected_msg = f"Radius is {radius}, but must be greater than 0"
        assert expected_msg in str(excinfo.value)
    else:
        result = func(
            df_neighbor_periodicity,
            radius=radius,
            periodicity=Periodicity.NONE,
            xlen=6,
            ylen=5,
        )
        assert result is not None


def test_distance_to_nearest_of_each_group_periodic():
    # values calculated using test data spreadsheet
    expected = pl.DataFrame({
        'id': pl.Series([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], dtype=pl.Int64),
        'type-1-dist': pl.Series([0.5, 0.5, 1.0, 0.5, 1, 1.118033988749895, 1.8027756377319946, 0.5, 1.118033988749895, 1, 1.4142135623730951, 1.4142135623730951], dtype=pl.Float64),
        'type-1-id': pl.Series([2, 1, 1, 1, 6, 9, 1, 6, 6, 1, 2, 1], dtype=pl.Int64),
        'type-2-dist': pl.Series([0.5, 0.7071, 1.118033988749895, 1.118033988749895, 1.5, 0.5, 1.414213562370951, 2.1213203435596424, 1.0, 1.414213562370951, 2.0, 1.0], dtype=pl.Float64),
        'type-2-id': pl.Series([4, 4, 4, 3, 4, 8, 4, 4, 8, 3, 8, 3], dtype=pl.Int64),
        'type-3-dist': pl.Series([1, 1.118033988749895, 1.414213562370951, 1.5, 2.1213203435596424, 1, 2.6926, 1.5, 1.8027756377319946, 2.1213203435596424, 1.5, 1.0], dtype=pl.Float64),
        'type-3-id': pl.Series([10, 5, 10, 5, 10, 5, 5, 5, 5, 5, 5, 10], dtype=pl.Int64),})
    result = nu_spa.distance_to_each_group(df_neighbor_periodicity_two,periodicity=nu_spa.Periodicity.XY, xlen=6, ylen=5)
    assert_frame_equal(result, expected, check_column_order=False, check_exact=False)


def test_distance_to_nearest_of_each_group_one_bug_in_a_group():
    # This is currently a repeat of the non-contig ID test. It is repeated here to show explicit intent.
    # During the non-contig testing the special case of 'only one bug in a group' was highlighted and the
    # underlying code was changed to follow the convention that a bug can be its nearest neighbor of group
    # (with dist 0) if it is the only one in that group.
    # values calculated using test data spreadsheet
    expected_non_periodic = pl.DataFrame({
        'id': pl.Series([14, 12, 2, 25, 4], dtype=pl.Int64),
        'type-1-dist': pl.Series([4.80104155366, 4.80104155366, 3.47131, 4.72757, 7.71686], dtype=pl.Float64),
        'type-1-id': pl.Series([12, 14, 12, 14, 12], dtype=pl.Int64),
        'type-2-dist': pl.Series([3.53836, 3.47131, 3.29089, 3.29089, 5.54707], dtype=pl.Float64),
        'type-2-id': pl.Series([2, 2, 25, 2, 25], dtype=pl.Int64),
        'type-3-dist': pl.Series([9.12030, 7.71686, 5.89067, 5.54707, 0], dtype=pl.Float64),
        'type-3-id': pl.Series([4, 4, 4, 4, 4], dtype=pl.Int64),})
    result = nu_spa.distance_to_each_group(df_neighbor_periodicity_non_contig_id,periodicity=nu_spa.Periodicity.NONE)
    assert_frame_equal(result, expected_non_periodic, check_column_order=False, check_exact=False)

    expected_periodic = pl.DataFrame({
        'id': pl.Series([14, 12, 2, 25, 4], dtype=pl.Int64),
        'type-1-dist': pl.Series([0.2236067, 0.2236067, 3.47131099, 3.33916, 0.591607978], dtype=pl.Float64),
        'type-1-id': pl.Series([12, 14, 12, 14, 12], dtype=pl.Int64),
        'type-2-dist': pl.Series([3.33916, 3.44963, 3.29089, 3.29089, 3.12569], dtype=pl.Float64),
        'type-2-id': pl.Series([25, 25, 25, 2, 25], dtype=pl.Int64),
        'type-3-dist': pl.Series([0.6164414, 0.591607978, 3.50713, 3.12569, 0], dtype=pl.Float64),
        'type-3-id': pl.Series([4, 4, 4, 4, 4], dtype=pl.Int64),})
    result = nu_spa.distance_to_each_group(df_neighbor_periodicity_non_contig_id,periodicity=nu_spa.Periodicity.XY, xlen=5, ylen=8)
    assert_frame_equal(result, expected_periodic, check_column_order=False, check_exact=False)


def test_distance_to_nearest_of_each_group_non_contig_coords():
    # values calculated using test data spreadsheet
    expected_non_periodic = pl.DataFrame({
        'id': pl.Series([14, 12, 2, 25, 4], dtype=pl.Int64),
        'type-1-dist': pl.Series([4.80104155366, 4.80104155366, 3.47131, 4.72757, 7.71686], dtype=pl.Float64),
        'type-1-id': pl.Series([12, 14, 12, 14, 12], dtype=pl.Int64),
        'type-2-dist': pl.Series([3.53836, 3.47131, 3.29089, 3.29089, 5.54707], dtype=pl.Float64),
        'type-2-id': pl.Series([2, 2, 25, 2, 25], dtype=pl.Int64),
        'type-3-dist': pl.Series([9.12030, 7.71686, 5.89067, 5.54707, 0], dtype=pl.Float64),
        'type-3-id': pl.Series([4, 4, 4, 4, 4], dtype=pl.Int64),})
    result = nu_spa.distance_to_each_group(df_neighbor_periodicity_non_contig_id,periodicity=nu_spa.Periodicity.NONE)
    assert_frame_equal(result, expected_non_periodic, check_column_order=False, check_exact=False)

    expected_periodic = pl.DataFrame({
        'id': pl.Series([14, 12, 2, 25, 4], dtype=pl.Int64),
        'type-1-dist': pl.Series([0.2236067, 0.2236067, 3.47131099, 3.33916, 0.591607978], dtype=pl.Float64),
        'type-1-id': pl.Series([12, 14, 12, 14, 12], dtype=pl.Int64),
        'type-2-dist': pl.Series([3.33916, 3.44963, 3.29089, 3.29089, 3.12569], dtype=pl.Float64),
        'type-2-id': pl.Series([25, 25, 25, 2, 25], dtype=pl.Int64),
        'type-3-dist': pl.Series([0.6164414, 0.591607978, 3.50713, 3.12569, 0], dtype=pl.Float64),
        'type-3-id': pl.Series([4, 4, 4, 4, 4], dtype=pl.Int64),})
    result = nu_spa.distance_to_each_group(df_neighbor_periodicity_non_contig_id,periodicity=nu_spa.Periodicity.XY, xlen=5, ylen=8)
    assert_frame_equal(result, expected_periodic, check_column_order=False, check_exact=False)


def test_distance_to_nearest_of_each_group_non_periodic():
    # values calculated using test data spreadsheet
    expected = pl.DataFrame({
        'id': pl.Series([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], dtype=pl.Int64),
        'type-1-dist': pl.Series([0.5, 0.5, 3.20156211871642, 0.5, 1.0, 1.1180339887499, 2.50, 0.5, 1.1180339887499, 1.500, 2.23606797749979, 3.20156211871642],
                                 dtype=pl.Float64),
        'type-1-id': pl.Series([2, 1, 9, 1, 6, 9, 12, 6, 6, 11, 9, 9], dtype=pl.Int64),
        'type-2-dist': pl.Series([0.5, 0.7071067811865, 1.500, 2.1213203435596, 1.50, 0.5, 1.500, 2.1213203435596, 1.0, 2.5, 2.0, 2.5],
                                 dtype=pl.Float64),
        'type-2-id': pl.Series([4, 4, 7, 8, 4, 8, 3, 4, 8, 8, 8, 7], dtype=pl.Int64),
        'type-3-dist': pl.Series([1.5811388300841898, 1.118033988749895, 3.5355339059327378, 1.500, 3.8078865529319543, 1.0, 3.640054944640259, 1.5, 1.8027767377, 3.80788655293, 1.5, 4.949747468],
                                 dtype=pl.Float64),
        'type-3-id': pl.Series([5, 5, 5, 5, 10, 5, 5, 5, 5, 5, 10, 5], dtype=pl.Int64), })
    result = nu_spa.distance_to_each_group(df_neighbor_periodicity_two,periodicity=nu_spa.Periodicity.NONE)
    assert_frame_equal(result, expected, check_column_order=False, check_exact=False)


def test_local_population_structure_periodic():
    expected_periodic = pl.DataFrame({
        'id': pl.Series([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], dtype=pl.Int64),
        '1': pl.Series([4, 5, 5, 6, 6, 4, 5, 5, 3, 5, 5, 3], dtype=pl.UInt32),
        '2': pl.Series([4, 4, 2, 3, 4, 4, 3, 2, 3, 4, 3, 3], dtype=pl.UInt32),
        '3': pl.Series([2, 2, 2, 2, 1, 2, 2, 2, 1, 1, 2, 2], dtype=pl.UInt32)})

    # For df_neighbor_periodicity_two, xlen=5, ylen=6, radius=3
    # Tallied by hand. ID, group type of neighbors, sorted, counts
    #1:  2 1 2 3 1 3 1 2 1 2     -> 1111   2222   33   -> 4 4 2
    #2:  1 2 3 3 1 2 1 1 2 2 1   -> 11111  2222   33   -> 5 4 2
    #3:  1 1 2 3 2 1 3 1 1       -> 11111  22     33   -> 5 2 2
    #4:  1 1 2 2 3 3 1 1 2 1 1   -> 111111 222    33   -> 6 3 2
    #5:  1 1 2 2 1 1 1 3 2 2 1   -> 111111 2222   3    -> 6 4 1
    #6:  2 3 1 1 2 1 1 2 3 2     -> 1111   2222   33   -> 4 4 2
    #7:  2 2 1 1 1 1 2 1 3 3     -> 11111  222    33   -> 5 3 2
    #8:  1 1 3 1 2 1 1 3 2       -> 11111  22     33   -> 5 2 2
    #9:  2 1 3 1 2 1 2           -> 111    222    3    -> 3 3 1
    #10: 1 1 1 2 2 1 3 2 2 1     -> 11111  2222   3    -> 5 4 1
    #11: 1 3 3 1 2 2 1 1 1 2     -> 11111  222    33   -> 5 3 2
    #12: 3 2 1 2 1 1 2 3         -> 111    222    33   -> 3 3 2

    local_pop = nu_spa.local_population_structure(df_neighbor_periodicity_two, radius= 3, periodicity=Periodicity.XY, xlen=6, ylen=5)
    assert_frame_equal(local_pop, expected_periodic, check_column_order=False)


def test_local_population_structure_non_periodic():
    # nb got the column ids screwed up when reading autogenerated
    expected = pl.DataFrame({
        'id': pl.Series([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], dtype=pl.Int64),
        '1': pl.Series([2, 3, 0, 4, 4, 4, 2, 5, 3, 2, 2, 0], dtype=pl.UInt32),
        '2': pl.Series([2, 2, 1, 1, 2, 2, 1, 1, 3, 1, 1, 1], dtype=pl.UInt32),
        '3': pl.Series([1, 1, 0, 1, 0, 2, 0, 2, 1, 0, 1, 0], dtype=pl.UInt32)})

    # For df_neighbor_periodicity_two, xlen=5, ylen=6, radius=3
    # Tallied by hand. ID, group type of neighbors, sorted, counts
    #1:  1 2 3 1 2          -> 11       22  3   -> 2 2 1
    #2:  1 2 3 1 2 1        -> 111      22  3   -> 3 2 1
    #3:  2                  ->          2       -> 0 1 0
    #4:  1 1 3 1 2 1        -> 1111     2   3   -> 4 1 1
    #5:  1 1 2 2 1 1        -> 1111     22      -> 4 2 0
    #6:  2 3 1 1 2 1 1 3    -> 1111     22  33  -> 4 2 2
    #7:  2 1 1              -> 11       2       -> 2 1 0
    #8:  1 1 3 1 2 1 1 3    -> 11111    2   33  -> 5 1 2
    #9:  2 1 3 1 2 1 2      -> 111      222 3   -> 3 3 1
    #10: 1 2 1              -> 11       2       -> 2 1 0
    #11: 3 2 1 1            -> 11       2   3   -> 2 1 1
    #12: 2                  ->          2       -> 0 1 0

    local_pop = nu_spa.local_population_structure(df_neighbor_periodicity_two, radius= 3, periodicity=Periodicity.NONE)
    assert_frame_equal(local_pop, expected, check_column_order=False)

###########################
# Input validation checks #
###########################

# TODO wrap all such of these into a pytest.mark.paramterize
@pytest.mark.skip(reason="periodicity validation")
def test_periodicity_validation():
    pass


@pytest.mark.parametrize(
    "kwargs",[
        {"xlen": 50},
        {"ylen": 50},
        {"zlen": 150},
        {"xlen": 50, "ylen": 70},
        {"ylen": 50, "xlen": 20},
        {"ylen": 50, "zlen": 10, "xlen": 20}],
)
def test_neighbors_radius_non_periodic_warns_if_dim_set(kwargs):
    with pytest.warns(UserWarning, match="is set but is not needed for no periodicity"):
        nu_spa.neighbors_radius(
            df_neighbor_periodicity,
            radius=6,
            periodicity=Periodicity.NONE,
            **kwargs
        )


def test_neighbors_radius_periodic_xy_dim_checks():
    with pytest.raises(ValueError) as excinfo:
        nu_spa.neighbors_radius(df_neighbor_periodicity, radius=6, periodicity=Periodicity.XY)
    assert f'Periodicity of "xy" specified but xlen and ylen are not set' in str(excinfo.value)

    with pytest.raises(ValueError) as excinfo:
        nu_spa.neighbors_radius(df_neighbor_periodicity, radius=6, periodicity=Periodicity.XY, ylen=20)
    assert f'Periodicity of "xy" specified but xlen is not set' in str(excinfo.value)

    with pytest.raises(ValueError) as excinfo:
        nu_spa.neighbors_radius(df_neighbor_periodicity, radius=6, periodicity=Periodicity.XY, xlen=50)
    assert f'Periodicity of "xy" specified but ylen is not set' in str(excinfo.value)

    with pytest.raises(ValueError) as excinfo:
        nu_spa.neighbors_radius(df_neighbor_periodicity, radius=6, periodicity=Periodicity.XY, xlen=-1, ylen=0)
    assert f'Periodicity of "xy" specified but xlen and ylen not > 0. xlen: -1, ylen: 0' in str(excinfo.value)

    with pytest.raises(ValueError) as excinfo:
        nu_spa.neighbors_radius(df_neighbor_periodicity, radius=6, periodicity=Periodicity.XY, xlen=-1, ylen=1)
    assert f'Periodicity of "xy" specified but xlen is not > 0. xlen: -1' in str(excinfo.value)

    with pytest.raises(ValueError) as excinfo:
        nu_spa.neighbors_radius(df_neighbor_periodicity, radius=6, periodicity=Periodicity.XY, xlen=1, ylen=-2)
    assert f'Periodicity of "xy" specified but ylen is not > 0. ylen: -2' in str(excinfo.value)


def test_min_xyz_len():
    #xlen or ylen not greater than max"
    max_x = 4.9
    bad_x = max_x-1
    good_x = max_x+1
    max_y = 7.9
    bad_y = max_y - 1
    good_y = max_y + 1
    with pytest.raises(ValueError) as excinfo:
        nu_spa.neighbors_radius(df_neighbor_periodicity, radius=6, periodicity=Periodicity.XY, xlen=bad_x, ylen=good_y)
    assert f'xlen is specified to {bad_x}, lower than max x-value of points in dataset: {max_x}' in str(excinfo.value)

    with pytest.raises(ValueError) as excinfo:
        nu_spa.neighbors_radius(df_neighbor_periodicity, radius=6, periodicity=Periodicity.XY, xlen=good_x, ylen=bad_y)
    assert f'ylen is specified to {bad_y}, lower than max y-value of points in dataset: {max_y}' in str(excinfo.value)

    with pytest.raises(ValueError) as excinfo:
         nu_spa.neighbors_radius(df_neighbor_periodicity, radius=6, periodicity=Periodicity.XY, xlen=bad_x, ylen=bad_y)
    assert f'xlen, ylen are {bad_x}, {bad_y}, lower than max values in dataset:{max_x} {max_y}' in str(excinfo.value)

    try:
        nu_spa.neighbors_radius(df_neighbor_periodicity, radius=6, periodicity=Periodicity.NONE, xlen=max_x, ylen=max_y)
    except Exception as e:
        pytest.fail(f'Unexpected exception {e}')


def test_neighbors_radius_periodic_xy_warn_zlen_set():
    with pytest.warns(UserWarning, match="is set but is not needed for XY periodicity"):
        neighbour_lists = nu_spa.neighbors_radius(df_neighbor_periodicity, radius=3.5, periodicity=Periodicity.XY,
                                              xlen=5, ylen=9, zlen=10)