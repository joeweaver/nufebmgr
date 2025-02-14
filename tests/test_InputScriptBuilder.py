import pytest
from nufebmgr.InputScriptBuilder import InputScriptBuilder
from nufebmgr.SimulationBox import SimulationBox

def test_grid_cell_size_picker():
    isb = InputScriptBuilder()
    s = SimulationBox() # 100 x 100 x 100 default
    grid_size = isb._pick_grid_size(s)
    assert grid_size == pytest.approx(2e-6)

    s = SimulationBox(xlen=90,ylen=90,zlen=60)
    grid_size = isb._pick_grid_size(s)
    assert grid_size == pytest.approx(1.5e-6)

    s = SimulationBox(xlen=140,ylen=126,zlen=56)
    grid_size = isb._pick_grid_size(s)
    assert grid_size == pytest.approx(2e-6)

    for expected in [1.5,2.0,2.5]:
        _test_a_grid(expected)

    # with pytest.raises(ValueError):
    #     s = SimulationBox(xlen=141,ylen=126,zlen=56)
    #     grid_size = isb._pick_grid_size(s)

    s = SimulationBox(xlen=141,ylen=126,zlen=56)
    grid_size = isb._pick_grid_size(s)
    assert grid_size == pytest.approx(1e-6)

def _test_a_grid(expected):
    isb = InputScriptBuilder()
    s = SimulationBox(xlen=expected*10, ylen=expected*9, zlen=expected*4)
    grid_size = isb._pick_grid_size(s)
    assert grid_size == pytest.approx(expected*1e-6)
