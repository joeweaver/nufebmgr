import pytest
from nufebmgr.InputScriptBuilder import InputScriptBuilder
from nufebmgr.SimulationBox import SimulationBox
from nufebmgr.Substrate import Substrate
from nufebmgr.BugDumpSpec import BugDumpSpec
from nufebmgr.ChemDumpSpec import ChemDumpSpec
from nufebmgr.HDF5DumpSpec import HDF5DumpSpec

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

def test_empty_hdf5():
    isb = InputScriptBuilder()
    # should start with empty values
    assert isb.config_vals['computation_output'][0]['hdf5_output'] == []


    dump_A = HDF5DumpSpec(dumpname="custom_dump.h5", dumpdir="hdf5", nsteps=4, dump_bugs=BugDumpSpec("all").hdf5_vars(),
                          dump_chems=ChemDumpSpec("conc").hdf5_vars())
    isb.add_hdf5_output([dump_A])

    isb.add_hdf5_output([])
    assert isb.config_vals['computation_output'][0]['hdf5_output'] == []

def test_single_hdf5():
    isb = InputScriptBuilder()

    # should start with empty values
    assert isb.config_vals['computation_output'][0]['hdf5_output'] == []

    expected_A = [
        {'linetype': 'subsection', 'title': 'HDF5 output, efficient binary format, for storing many atom properties'},
        {'linetype': 'comment', 'title': 'NOTE: requires NUFEB built with HDF5 option'},
        {'linetype': 'subsection', 'title': 'Create directory(s) for dump'},
        {'linetype': 'command', 'name': 'shell', 'command': 'mkdir hdf5'},
        {'linetype': 'subsection', 'title': 'Dump specifications'},
        {'linetype': 'command', 'name': 'dump', 'dumpname': 'du_hdf5_0', 'group': 'all', 'format': 'nufeb/hdf5',
         'every_n': '4',
         'loc': 'hdf5/custom_dump.h5', 'dumpvars': 'id type x y z vx vy vz fx fy fz radius conc'},
    ]


    dump_A = HDF5DumpSpec(dumpname="custom_dump.h5", dumpdir="hdf5", nsteps=4, dump_bugs=BugDumpSpec("all").hdf5_vars(),
                          dump_chems=ChemDumpSpec("conc").hdf5_vars())
    isb.add_hdf5_output([dump_A])
    assert isb.config_vals['computation_output'][0]['hdf5_output'] == expected_A


def test_multi_hdf5():
    isb = InputScriptBuilder()

    # should start with empty values
    assert isb.config_vals['computation_output'][0]['hdf5_output'] == []

    expected_ABC = [
        {'linetype': 'subsection', 'title': 'HDF5 output, efficient binary format, for storing many atom properties'},
        {'linetype': 'comment', 'title': 'NOTE: requires NUFEB built with HDF5 option'},
        {'linetype': 'subsection', 'title': 'Create directory(s) for dump'},
        {'linetype': 'command', 'name': 'shell', 'command': 'mkdir alt_hdf5'},
        {'linetype': 'command', 'name': 'shell', 'command': 'mkdir hdf5'},
        {'linetype': 'subsection', 'title': 'Dump specifications'},
        {'linetype': 'command', 'name': 'dump', 'dumpname': 'du_hdf5_0', 'group': 'all', 'format': 'nufeb/hdf5',
         'every_n': '4', 'loc': 'hdf5/custom_dump.h5', 'dumpvars': 'id type x y z vx vy vz fx fy fz radius conc'},
        {'linetype': 'command', 'name': 'dump', 'dumpname': 'du_hdf5_1', 'group': 'all', 'format': 'nufeb/hdf5',
         'every_n': '2', 'loc': 'hdf5/custom_dumpB.h5', 'dumpvars': 'id type x y z reac'},
        {'linetype': 'command', 'name': 'dump', 'dumpname': 'du_hdf5_2', 'group': 'all', 'format': 'nufeb/hdf5',
         'every_n': '9', 'loc': 'alt_hdf5/custom_dump.h5', 'dumpvars': 'id type x y z vx vy vz fx fy fz radius conc reac'},
    ]

    dump_A = HDF5DumpSpec(dumpname="custom_dump.h5", dumpdir="hdf5", nsteps=4, dump_bugs=BugDumpSpec("all").hdf5_vars(),
                          dump_chems=ChemDumpSpec("conc").hdf5_vars())
    dump_B = HDF5DumpSpec(dumpname="custom_dumpB.h5", dumpdir="hdf5", nsteps=2, dump_bugs=BugDumpSpec("location").hdf5_vars(),
                          dump_chems=ChemDumpSpec("reac").hdf5_vars())
    dump_C = HDF5DumpSpec(dumpname="custom_dump.h5", dumpdir="alt_hdf5", nsteps=9, dump_bugs=BugDumpSpec("all").hdf5_vars(),
                          dump_chems=ChemDumpSpec("all").hdf5_vars())

    isb.add_hdf5_output([dump_A, dump_B, dump_C])
    assert isb.config_vals['computation_output'][0]['hdf5_output'] == expected_ABC

def test_add_vol_computes():
    isb = InputScriptBuilder()

    # should start with empty values
    assert isb.config_vals['computation_output'][0]['vol_track'] == []

    isb.group_assignments = {'taxon_a': '1', 'taxon_b': '2'}
    isb.add_vol_vars()
    expected_a = [{'name': 'Compute per-taxa volumes'},
                {'command': 'compute', 'name': 'vol_taxon_a', 'group': 'taxon_a', 'fix_loc': 'nufeb/volume'},
                {'command': 'variable', 'name': 'var_vol_taxon_a', 'op': 'equal', 'to': '"c_vol_taxon_a"'},
                {'command': 'compute', 'name': 'vol_taxon_b', 'group': 'taxon_b', 'fix_loc': 'nufeb/volume'},
                {'command': 'variable', 'name': 'var_vol_taxon_b', 'op': 'equal', 'to': '"c_vol_taxon_b"'},]

    assert isb.config_vals['computation_output'][0]['vol_track'] == expected_a

    isb = InputScriptBuilder()

    # should start with empty values
    assert isb.config_vals['computation_output'][0]['vol_track'] == []

    isb.group_assignments = {'bug3': '3', 'bug4': '4', 'bugn': 'n'}
    isb.add_vol_vars()
    expected_b = [{'name': 'Compute per-taxa volumes'},
                  {'command': 'compute', 'name': 'vol_bug3', 'group': 'bug3', 'fix_loc': 'nufeb/volume'},
                  {'command': 'variable', 'name': 'var_vol_bug3', 'op': 'equal', 'to': '"c_vol_bug3"'},
                  {'command': 'compute', 'name': 'vol_bug4', 'group': 'bug4', 'fix_loc': 'nufeb/volume'},
                  {'command': 'variable', 'name': 'var_vol_bug4', 'op': 'equal', 'to': '"c_vol_bug4"'},
                  {'command': 'compute', 'name': 'vol_bugn', 'group': 'bugn', 'fix_loc': 'nufeb/volume'},
                  {'command': 'variable', 'name': 'var_vol_bugn', 'op': 'equal', 'to': '"c_vol_bugn"'},]

    assert isb.config_vals['computation_output'][0]['vol_track'] == expected_b