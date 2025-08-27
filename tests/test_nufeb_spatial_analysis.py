import pytest
import polars as pl
from polars.testing import assert_frame_equal
import numpy as np
from sklearn.neighbors import KDTree
from nufebmgr import nufeb_spatial_analysis as nu_spa

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

def test_smoke():
    coords = df_simple.select(["x", "y", "z"])
    assert coords['x'][1] == 1.0
    assert coords['y'][2] == 2.0
    assert coords['z'][3] == 5.0

def test_neighbors_radius_unsupported_periodicity():
    periodicities = ["none", "xy"]
    with pytest.raises(ValueError) as excinfo:
        nu_spa.neighbors_radius(df_neighbor_periodicity, radius=2, periodicity='banana', xlen=5, ylen=8)
    assert f'Unrecognized periodicity: banana. Must be one of {",".join(periodicities)}' in str(excinfo.value)

def test_neighbors_radius_non_periodic_warns_if_dim_set():
    with pytest.warns(UserWarning, match="is set but is not needed for no periodicity"):
        nu_spa.neighbors_radius(df_neighbor_periodicity, radius=6, periodicity='none', xlen=50)

    with pytest.warns(UserWarning, match="is set but is not needed for no periodicity"):
        nu_spa.neighbors_radius(df_neighbor_periodicity, radius=6, periodicity='none', ylen=50)

    with pytest.warns(UserWarning, match="is set but is not needed for no periodicity"):
        nu_spa.neighbors_radius(df_neighbor_periodicity, radius=6, periodicity='none', zlen=150)

    with pytest.warns(UserWarning, match="is set but is not needed for no periodicity"):
        nu_spa.neighbors_radius(df_neighbor_periodicity, radius=6, periodicity='none', xlen=50, ylen=70)

    with pytest.warns(UserWarning, match="is set but is not needed for no periodicity"):
        nu_spa.neighbors_radius(df_neighbor_periodicity, radius=6, periodicity='none', ylen=50, xlen=20)

    with pytest.warns(UserWarning, match="is set but is not needed for no periodicity"):
        nu_spa.neighbors_radius(df_neighbor_periodicity, radius=6, periodicity='none', ylen=50, zlen=10, xlen=20)

    with pytest.warns(UserWarning, match="is set but is not needed for no periodicity"):
        nu_spa.neighbors_radius(df_neighbor_periodicity, radius=6, periodicity='none', xlen=50)


def test_neighbors_radius_non_periodic():
    # array([[0.        , 4.80104155, 3.5383612 , 4.72757866, 9.12030701],
    #        [4.80104155, 0.        , 3.47131099, 6.59545298, 7.71686465],
    #        [3.5383612 , 3.47131099, 0.        , 3.29089653, 5.89067059],
    #        [4.72757866, 6.59545298, 3.29089653, 0.        , 5.5470713 ],
    #        [9.12030701, 7.71686465, 5.89067059, 5.5470713 , 0.        ]])
    neighbour_lists = nu_spa.neighbors_radius(df_neighbor_periodicity, radius=6, periodicity='none')
    expected = { 1: [3, 4, 2],
                 2: [3, 1],
                 3: [4, 2, 1, 5],
                 4: [3, 1, 5],
                 5: [4, 3]}
    assert neighbour_lists == expected

    neighbour_lists = nu_spa.neighbors_radius(df_neighbor_periodicity_two, radius=3, periodicity='none')
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
    neighbour_lists = nu_spa.neighbors_radius(df_neighbor_periodicity_non_contig_id, radius=6, periodicity='none')
    expected = {14: [2, 25, 12],
                12: [2, 14],
                2: [25, 12, 14, 4],
                25: [2, 14, 4],
                4: [25, 2]}
    assert neighbour_lists == expected
    ## different setup, just to be careful

@pytest.mark.skip(reason="check dists internally for certainty")
def test_neighbours_radius_dist():
    pass

def test_neighbors_radius_periodic_xy_dim_checks():
    with pytest.raises(ValueError) as excinfo:
        nu_spa.neighbors_radius(df_neighbor_periodicity, radius=6, periodicity='xy')
    assert f'Periodicity of "xy" specified but xlen and ylen are not set' in str(excinfo.value)

    with pytest.raises(ValueError) as excinfo:
        nu_spa.neighbors_radius(df_neighbor_periodicity, radius=6, periodicity='xy', ylen=20)
    assert f'Periodicity of "xy" specified but xlen is not set' in str(excinfo.value)

    with pytest.raises(ValueError) as excinfo:
        nu_spa.neighbors_radius(df_neighbor_periodicity, radius=6, periodicity='xy', xlen=50)
    assert f'Periodicity of "xy" specified but ylen is not set' in str(excinfo.value)

    with pytest.raises(ValueError) as excinfo:
        nu_spa.neighbors_radius(df_neighbor_periodicity, radius=6, periodicity='xy', xlen=-1, ylen=0)
    assert f'Periodicity of "xy" specified but xlen and ylen not > 0. xlen: -1, ylen: 0' in str(excinfo.value)

    with pytest.raises(ValueError) as excinfo:
        nu_spa.neighbors_radius(df_neighbor_periodicity, radius=6, periodicity='xy', xlen=-1, ylen=1)
    assert f'Periodicity of "xy" specified but xlen is not > 0. xlen: -1' in str(excinfo.value)

    with pytest.raises(ValueError) as excinfo:
        nu_spa.neighbors_radius(df_neighbor_periodicity, radius=6, periodicity='xy', xlen=1, ylen=-2)
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
        nu_spa.neighbors_radius(df_neighbor_periodicity, radius=6, periodicity='xy', xlen=bad_x, ylen=good_y)
    assert f'xlen is specified to {bad_x}, lower than max x-value of points in dataset: {max_x}' in str(excinfo.value)

    with pytest.raises(ValueError) as excinfo:
        nu_spa.neighbors_radius(df_neighbor_periodicity, radius=6, periodicity='xy', xlen=good_x, ylen=bad_y)
    assert f'ylen is specified to {bad_y}, lower than max y-value of points in dataset: {max_y}' in str(excinfo.value)

    with pytest.raises(ValueError) as excinfo:
         nu_spa.neighbors_radius(df_neighbor_periodicity, radius=6, periodicity='xy', xlen=bad_x, ylen=bad_y)
    assert f'xlen, ylen are {bad_x}, {bad_y}, lower than max values in dataset:{max_x} {max_y}' in str(excinfo.value)

    try:
        nu_spa.neighbors_radius(df_neighbor_periodicity, radius=6, periodicity='none', xlen=max_x, ylen=max_y)
    except Exception as e:
        pytest.fail(f'Unexpected exception {e}')

@pytest.mark.skip(reason="periodic_xy should warn if zlen is set")
def test_neighbors_radius_periodic_xy_warn_too_low():
    max_x = 4.9
    max_y = 7.9
    neighbour_lists = nu_spa.neighbors_radius(df_neighbor_periodicity, radius=3.5, periodicity='xy', xlen=5, ylen=9)
    expected = {0: [1, 4],
                1: [0, 4, 2],
                2: [3, 1],
                3: [4, 2],
                4: [0, 1, 3]}
    assert neighbour_lists == expected
    pass

def test_neighbors_radius_periodic_xy():
    # distance matrix for x=5, y=9. Lower triangle. truncated without rounding
    #   0           1           2           3           4
    # 0 0
    # 1 0.2236067   0
    # 2 3.5382612   3.471310    0
    # 3 4.3301270   4.438468    3.2908965   0
    # 4 1.3341664   1.396424    4.3011626   3.1256999   0

    neighbour_lists = nu_spa.neighbors_radius(df_neighbor_periodicity, radius=3.5, periodicity='xy', xlen=5, ylen=9)
    expected = { 1: [2, 5],
                 2: [1, 5, 3],
                 3: [4, 2],
                 4: [5, 3],
                 5: [1, 2, 4]}
    assert neighbour_lists == expected

    # tighter radius
    neighbour_lists = nu_spa.neighbors_radius(df_neighbor_periodicity, radius=2, periodicity='xy', xlen=5, ylen=9)
    expected = { 1: [2, 5],
                 2: [1, 5],
                 3: [],
                 4: [],
                 5: [1, 2]}
    assert neighbour_lists == expected

    # using non-contiguous, out of numeric order id column to check if ids are correctly preserved
    # For ref:      'id': [14, 12, 2, 25, 4],
    neighbour_lists = nu_spa.neighbors_radius(df_neighbor_periodicity_non_contig_id, radius=2, periodicity='xy', xlen=5, ylen=9)
    expected = { 14: [12, 4],
                 12: [14, 4],
                 2: [],
                 25: [],
                 4: [14, 12]}
    assert neighbour_lists == expected

    # different setup, just to be extra careful
    neighbour_lists = nu_spa.neighbors_radius(df_neighbor_periodicity_two, radius=3, periodicity='xy', xlen=6, ylen=5)
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

@pytest.mark.skip(reason="extensive varied periodic bc checks to be sure")
def test_varied_periodic_check():
    pass

@pytest.mark.skip(reason="case where it's within radius for both periodic and non and it's closer when periodic/when not periodic")
def test_more_periodic_dumbness():
    pass
@pytest.mark.skip(reason="what about exactly on 0 or boundary or right at radius")
def test_neighbors_radius_edge():
    pass



@pytest.mark.skip(reason="n neighbors")
def test_n_neighbors():
    pass

@pytest.mark.skip(reason="nearest of type")
def test_distance_to_nearest_of_group():
    pass

@pytest.mark.skip(reason="xlen, ylen, raises, etc for population structure")
def test_local_population_structure_raises():
    pass

def test_local_population_structure_periodic():
    # nb got the column ids screwed up when reading autogenerated
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

    local_pop = nu_spa.local_population_structure(df_neighbor_periodicity_two, radius= 3, periodicity="xy", xlen=6, ylen=5)
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

    local_pop = nu_spa.local_population_structure(df_neighbor_periodicity_two, radius= 3, periodicity="none")
    assert_frame_equal(local_pop, expected, check_column_order=False)