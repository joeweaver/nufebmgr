import pytest
from nufebmgr.InputScriptBuilder import InputScriptBuilder
from nufebmgr.SimulationBox import SimulationBox
from nufebmgr.Substrate import Substrate

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
    s = Substrate(name='test_sub', init_concentration='1e-3', bulk_concentration='1e-4', diffusion_coefficient=2.0e-9,
                  biofilm_diffusion_ratio=0.5)
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
    with pytest.raises(ValueError) as excinfo:
        isb._build_grid_modify_dict(s, boundary_scenario)
    assert f'Invalid boundary_scenario:' in str(excinfo.value)

def test_build_diffusion():
    sa = Substrate(name='test_sub_a', init_concentration='1e-3', bulk_concentration='1.5e-4',
                  diffusion_coefficient=2.0e-9,
                  biofilm_diffusion_ratio=0.5)
    sb = Substrate(name='test_sub_b', init_concentration='1.1e-3', bulk_concentration='1.4e-4',
                  diffusion_coefficient=2.1e-9,
                  biofilm_diffusion_ratio=0.7)

    isb = InputScriptBuilder()

    subs_set1={'test_sub_a':sa, 'test_sub_b':sb}

    diffusion_stuff = isb.config_vals['chemical_processes'][0]
    # should start with empty values
    assert [] == diffusion_stuff['diffusion_coefficients']
    assert [] == diffusion_stuff['diffusion_biofilm_ratios']
    isb.build_diffusion(subs_set1)
    diffusion_stuff = isb.config_vals['chemical_processes'][0]

    expected = {'title':'#------------------------------Chemical Processes------------------------------#',
                'diffusion_biofilm_ratios': [
                    {'title': 'Ratio of diffusion in biofilm as compared to water'},
                    {'name': 'fix',
                     'fix_name': f'coeff_test_sub_a',
                     'group': 'all',
                     'fix_loc': 'nufeb/diffusion_coeff',
                     'sub1': 'test_sub_a',
                     'coeff1': 0.5,
                     'comment': ''},
                    {'name': 'fix',
                     'fix_name': f'coeff_test_sub_b',
                     'group': 'all',
                     'fix_loc': 'nufeb/diffusion_coeff',
                     'sub1': 'test_sub_b',
                     'coeff1': 0.7,
                     'comment': ''}
                ],
                'diffusion_coefficients': [
                    {'title': 'Diffusion in water'},
                    {'name': 'fix',
                         'fix_name': f'diff_test_sub_a',
                         'group': 'all',
                         'fix_loc': 'nufeb/diffusion_reaction',
                         'sub1': 'test_sub_a',
                         'coeff1': 2.0e-9,
                         'comment': ''},
                    {'name': 'fix',
                     'fix_name': f'diff_test_sub_b',
                     'group': 'all',
                     'fix_loc': 'nufeb/diffusion_reaction',
                     'sub1': 'test_sub_b',
                     'coeff1': 2.1e-9,
                     'comment': ''}
                ]
                }

    assert expected == diffusion_stuff


def test_build_postphysical():
    isb = InputScriptBuilder()
    current_data = isb.config_vals['post_physical_processes'][0]
    # should start with empty values
    assert [] == current_data['boundary_layer']
    isb.build_post_physical(20)
    new_data = isb.config_vals['post_physical_processes'][0]

    expected = {'title': '#---------------------------Post-Physical Processes----------------------------#',
                'boundary_layer': [
                    {'title': 'Elastic boundary layer above biofilm surface'},
                    {'name': 'fix',
                     'fix_name': 'blayer',
                     'group': 'all',
                     'fix_loc': 'nufeb/boundary_layer',
                     'distance_m': "20e-6",
                     'comment': ''}
                ]
                }

    assert expected == new_data
