import pytest
import polars as pl
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
# from scipy.spatial import distance_matrix
# coords = df_neighbor_periodicity.select(['x','y','z']).to_numpy()
# dist_mat = distance_matrix(coords, coords)
# non-periodic distance matrix for above
# array([[0.        , 4.80104155, 3.5383612 , 4.72757866, 9.12030701],
#        [4.80104155, 0.        , 3.47131099, 6.59545298, 7.71686465],
#        [3.5383612 , 3.47131099, 0.        , 3.29089653, 5.89067059],
#        [4.72757866, 6.59545298, 3.29089653, 0.        , 5.5470713 ],
#        [9.12030701, 7.71686465, 5.89067059, 5.5470713 , 0.        ]])

def test_smoke():
    coords = df_simple.select(["x", "y", "z"])
    assert coords['x'][1] == 1.0
    assert coords['y'][2] == 2.0
    assert coords['z'][3] == 5.0

def test_index_gaurantee():
    # just reassuring ourselves that we can map back to id and group
    # some of this also serves a reminders on using kdtree, etc

    coords = df_simple.select(["x", "y", "z"]).to_numpy()
    kd_tree = KDTree(coords, metric="euclidean")

    #get nearest neighbor of ID 3
    row_for_ID3 = 2
    n_nearest = 1
    dist, ind = kd_tree.query(coords[row_for_ID3].reshape(1, -1),n_nearest+1)
    row_for_nearest = df_simple[int(ind[0][1])]
    nearest_id = row_for_nearest['id'][0]
    assert nearest_id == 2

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
    expected = { 0: [2, 3, 1],
                 1: [2, 0],
                 2: [3, 1, 0, 4],
                 3: [2, 0, 4],
                 4: [3, 2]}
    assert neighbour_lists == expected


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


@pytest.mark.skip(reason="periodic should warn if zlen is set")
def test_dimset_xy():
    pass

@pytest.mark.skip(reason="Left off here")
def test_neighbors_radius_periodic_xy():
    neighbour_lists = nu_spa.neighbors_radius(df_neighbor_periodicity, radius=6, periodicity='xy', xlen=5, ylen=8)
    assert False

@pytest.mark.skip(reason="extensive varied periodic bc checks to be sure")
def test_varied_periodic_check():
    pass

@pytest.mark.skip(reason="case where it's within radius for both periodic and non and it's closer when periodic/when not periodic")
def test_more_periodic_dumbness():
    pass
@pytest.mark.skip(reason="what about exactly on 0 or boundary or right at radius")
def test_neighbors_radius_edge():
    pass

@pytest.mark.skip(reason="single query for neighbours")
def test_neighbors_radius_single():
    pass

@pytest.mark.skip(reason="all query for neighbours")
def test_neighbours_radius_all():
    pass

@pytest.mark.skip(reason="n neighbors")
def test_n_neighbors():
    pass

@pytest.mark.skip(reason="nearest of type")
def test_nearest_of_type():
    pass
