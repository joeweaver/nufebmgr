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
    s = Substrate(name='test_sub', init_concentration='1e-3', bulk_concentration='1e-4')
    isb = InputScriptBuilder()

    boundary_scenario = "bioreactor"
    result = isb._build_grid_modify_dict(s, boundary_scenario)
    expected = {'name': 'grid_modify',
                   'action': 'set',
                   'substrate': 'test_sub',
                   'xbound': 'pp',
                   'ybound': 'pp',
                   'zbound': 'nd',
                   'init_conc': '1e-3',
                   'bulk-kw': 'bulk',
                   'bulkd_conc': '1e-4'}
    assert expected == result

    boundary_scenario = "microwell"
    result = isb._build_grid_modify_dict(s, boundary_scenario)
    expected = {'name': 'grid_modify',
                   'action': 'set',
                   'substrate': 'test_sub',
                   'xbound': 'nn',
                   'ybound': 'nn',
                   'zbound': 'nn',
                   'init_conc': '1e-3',
                   'bulk-kw': 'bulk',
                   'bulkd_conc': '1e-4'}
    assert expected == result

    boundary_scenario = "floating"
    result = isb._build_grid_modify_dict(s, boundary_scenario)
    expected = {'name': 'grid_modify',
                   'action': 'set',
                   'substrate': 'test_sub',
                   'xbound': 'dd',
                   'ybound': 'dd',
                   'zbound': 'dd',
                   'init_conc': '1e-3',
                   'bulk-kw': 'bulk',
                   'bulkd_conc': '1e-4'}
    assert expected == result

    boundary_scenario = "agar"
    result = isb._build_grid_modify_dict(s, boundary_scenario)
    expected = {'name': 'grid_modify',
                   'action': 'set',
                   'substrate': 'test_sub',
                   'xbound': 'pp',
                   'ybound': 'pp',
                   'zbound': 'dn',
                   'init_conc': '1e-3',
                   'bulk-kw': 'bulk',
                   'bulkd_conc': '1e-4'}
    assert expected == result

    boundary_scenario = "unknown"
    with pytest.raises(KeyError) as excinfo:
        isb._build_grid_modify_dict(s, boundary_scenario)
    assert f'{boundary_scenario }' in str(excinfo.value)