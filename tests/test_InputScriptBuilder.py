import pytest
from nufebmgr.InputScriptBuilder import InputScriptBuilder
from nufebmgr.SimulationBox import SimulationBox
from nufebmgr.NufebProject import Substrate

def test_grid_cell_size_picker():
    isb = InputScriptBuilder()
    s = SimulationBox() # 100 x 100 x 100 default
    grid_size = isb._pick_grid_size(s)
    assert grid_size == pytest.approx(2.0)

    s = SimulationBox(xlen=90,ylen=90,zlen=60)
    grid_size = isb._pick_grid_size(s)
    assert grid_size == pytest.approx(2.0)

    s = SimulationBox(xlen=140,ylen=126,zlen=56)
    grid_size = isb._pick_grid_size(s)
    assert grid_size == pytest.approx(2.0)

    for expected in [1.5,2.0]:
        _test_a_grid(expected)

    # with pytest.raises(ValueError):
    #     s = SimulationBox(xlen=141,ylen=126,zlen=56)
    #     grid_size = isb._pick_grid_size(s)

    s = SimulationBox(xlen=141,ylen=126,zlen=56)
    grid_size = isb._pick_grid_size(s)
    assert grid_size == pytest.approx(1)


    assert grid_size == pytest.approx(1)

    s = SimulationBox(xlen=141.5, ylen=126, zlen=56)
    with pytest.raises(ValueError) as excinfo:
        grid_size = isb._pick_grid_size(s)
    assert f'No valid grid size (option 1, 1.5, 2 microns) for a simulation of dimensions {s.dim_string()}' in str(excinfo.value)

def _test_a_grid(expected):
    isb = InputScriptBuilder()
    s = SimulationBox(xlen=expected*10, ylen=expected*9, zlen=expected*4)
    grid_size = isb._pick_grid_size(s)
    assert grid_size == pytest.approx(expected)

def test_build_grid_modify_dict():
    s = Substrate(name='test_sub', init_concentration='1e-3', bulk_concentration='1e-4',
              x_boundaries='pp',
              y_boundaries='pp',
              z_boundaries='nd')

    isb = InputScriptBuilder()
    result = isb._build_grid_modify_dict(s)

    expected = s.as_grid_modify_dict()
    assert expected == result