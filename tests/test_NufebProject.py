import pytest
from nufebmgr.NufebProject import NufebProject
from typing import Any
from pathlib import Path
from nufebmgr.ChemDumpSpec import ChemDumpSpec
from nufebmgr.BugDumpSpec import BugDumpSpec
from nufebmgr.HDF5DumpSpec import HDF5DumpSpec

HERE = Path(__file__).parent  # directory containing the test file
DATA_DIR = HERE / "data"

@pytest.fixture
def prj() -> "NufebProject":
    """
    Fixture that returns a fresh instance of NufebProject for each test.
    """
    return NufebProject()

def test_initialization():
    project = NufebProject()
    assert project is not None

def test_error_on_assign_taxa_not_all_taxa_have_entries_or_compositions():
    def setup_local(taxa_filename:str) -> NufebProject:
        prj = NufebProject()
        prj.add_taxon_by_jsonfile(taxa_filename)
        prj.layout_uniform(nbugs=20)
        prj.set_composition({'anammox_two_pathway': '33.333',
                             'denitrifier': '33.33',
                             'imperfect_denitrifier_no': '33.33'})
        prj.distribute_spatially_even()
        return prj

    with pytest.raises(ValueError) as excinfo:
        setup_local(DATA_DIR / "example_het_taxa.json")._assign_taxa()
    assert f'Setting composition with unknown taxon: anammox_two_pathway, denitrifier, imperfect_denitrifier_no' in str(excinfo.value)

    try:
        setup_local(DATA_DIR / "example_nitrogen_anaerobic_taxa.json")._assign_taxa()
    except Exception as e:
        pytest.fail(f'Unexpected exception {e}')

    with pytest.raises(ValueError) as excinfo:
        setup_local(DATA_DIR / "example_nitrogen_anaerobic_extra_het_taxa.json")._assign_taxa()
    assert f'Composition not set for taxon taxon (use 0 if not initially present): small_heterotroph' in str(excinfo.value)

    with pytest.raises(ValueError) as excinfo:
        setup_local(DATA_DIR / "example_nitrogen_anaerobic_taxa_missing_amx.json")._assign_taxa()
    assert f'Setting composition with unknown taxon: anammox_two_pathway' in str(excinfo.value)

    with pytest.raises(ValueError) as excinfo:
        setup_local(DATA_DIR / "example_nitrogen_anaerobic_taxa_missing_amx_extra_het.json")._assign_taxa()
    assert f'Setting composition with unknown taxon: anammox_two_pathway' in str(excinfo.value)


@pytest.mark.parametrize(
    "setter_name, getter_name, input_value, expected_value",
    [
        ("run_for_N_steps", "run_steps", 42, 42),
        ("use_seed", "seed", 19, 19),
        ("use_seed", "seed", None, 1701),
        ("set_biological_timestep_size_s", "biostep", 1100, 1100),
    ]
)
def test_API_contract(prj:NufebProject, setter_name: str, getter_name:str, input_value:Any, expected_value:Any) -> None:
    """
    Test the public API contract for NufebProject setter/getter pairs.

    This test verifies that:
      - The setter correctly updates the internal state - apart from side effects
      - The corresponding getter returns the expected value.
      - Default parameter behavior (if applicable) works as intended.

    Parameters
    ----------
    prj : NufebProject
        Instance of NufebProject provided by the fixture.
    setter_name : str
        Name of the setter method to call on the instance.
    getter_name : str
        Name of the getter attribute or method to check the value.
    input_value : Any
        Value to pass to the setter. Can be any type compatible with the setter.
    expected_value : Any
        The value expected to be returned by the getter after calling the setter.

    Notes
    -----
      - Side effects are not tested here and should be covered in
        separate, they will be tested in dedicated explicit tests
      - The goal here is to avoid boilerplate code.
      - We are testing simple getter/setters as this is part of the Public facing API and we want to know if *anything*
      about that breaks. In the future, 'simple' getter/setters may no longer be such.
    """
    #prj = NufebProject()
    setter = getattr(prj, setter_name)
    if input_value is None:
        setter()
    else:
        setter(input_value)
    getter_value = getattr(prj, getter_name)

    assert getter_value == expected_value

@pytest.mark.parametrize(
    "hours, biostep, expected",
    [
        (2, 900, 8),
        (3, None, 12),
        (2, 1800, 4),
        (2.5, 900, 10),
        (1.25, 900, 5),
        (1.24, 900, 4),
        (1.26, 900, 5),
        (1.1, 1, int(60*60*1.1)),
    ]
)
def test_run_for_N_hours(prj:NufebProject, hours: float, biostep:int, expected:int) -> None:
    if biostep is None:
        prj.run_for_N_hours(hours)
    else:
        prj.run_for_N_hours(hours,biostep)
    assert prj.run_steps == expected

def test_default_single_custom_hdf5(prj:NufebProject,):
    bds = BugDumpSpec()
    cds = ChemDumpSpec()
    prj.clear_hdf5_output()
    assert prj.hdf5_dump_specs == []

    prj.add_custom_hdf5_output()
    assert len(prj.hdf5_dump_specs) == 1
    ho = prj.hdf5_dump_specs[0]
    expected = HDF5DumpSpec(dumpname ="dump.h5", dumpdir ="hdf5", nsteps =1, dump_bugs =BugDumpSpec().hdf5_vars(), dump_chems =ChemDumpSpec().hdf5_vars())
    assert ho == expected

def test_single_custom_hdf5(prj: NufebProject, ):
    bds = BugDumpSpec()
    cds = ChemDumpSpec()
    prj.clear_hdf5_output()
    assert prj.hdf5_dump_specs == []

    prj.add_custom_hdf5_output(dumpname="custom_dump.h5", dumpdir="alt_hdf5", nsteps=10, dump_bugs=BugDumpSpec("location"),
                            dump_chems=ChemDumpSpec("reac"))
    assert len(prj.hdf5_dump_specs) == 1
    expected = HDF5DumpSpec(dumpname="custom_dump.h5", dumpdir="alt_hdf5", nsteps=10, dump_bugs=BugDumpSpec("location").hdf5_vars(),
                            dump_chems=ChemDumpSpec("reac").hdf5_vars())
    assert prj.hdf5_dump_specs[0] == expected

def test_multiple_custom_hdf5(prj: NufebProject, ):
    bds = BugDumpSpec()
    cds = ChemDumpSpec()
    prj.clear_hdf5_output()
    assert prj.hdf5_dump_specs == []

    prj.add_custom_hdf5_output(dumpname="custom_dump.h5", dumpdir="alt_hdf5", nsteps=1, dump_bugs=BugDumpSpec("location"),
                            dump_chems=ChemDumpSpec("reac"))
    prj.add_custom_hdf5_output(dumpname="custom_dumpb.h5", dumpdir="altb_hdf5", nsteps=1, dump_bugs=BugDumpSpec("all"),
                            dump_chems=ChemDumpSpec("conc"))
    assert len(prj.hdf5_dump_specs) == 2
    expected0 = HDF5DumpSpec(dumpname="custom_dump.h5", dumpdir="alt_hdf5", nsteps=1, dump_bugs=BugDumpSpec("location").hdf5_vars(),
                            dump_chems=ChemDumpSpec("reac").hdf5_vars())
    expected1 = HDF5DumpSpec(dumpname="custom_dumpb.h5", dumpdir="altb_hdf5", nsteps=1, dump_bugs=BugDumpSpec("all").hdf5_vars(),
                            dump_chems=ChemDumpSpec("conc").hdf5_vars())
    assert prj.hdf5_dump_specs[0] == expected0
    assert prj.hdf5_dump_specs[1] == expected1

def test_multiple_custom_no_file_clobber_hdf5(prj: NufebProject ):
    bds = BugDumpSpec()
    cds = ChemDumpSpec()
    prj.clear_hdf5_output()
    assert prj.hdf5_dump_specs == []

    prj.add_custom_hdf5_output(dumpname="custom_dump.h5", dumpdir="hdf5", nsteps=1, dump_bugs=BugDumpSpec("location"),
                            dump_chems=ChemDumpSpec("reac"))

    with pytest.raises(ValueError) as excinfo:
        prj.add_custom_hdf5_output(dumpname="custom_dump.h5", dumpdir="hdf5", nsteps=1,
                                   dump_bugs=BugDumpSpec("all"),
                                   dump_chems=ChemDumpSpec("conc"))
    assert f'Specified multiple HDF5 dump files named custom_dump.h5 in directory hdf5' in str(excinfo.value)

    # should not raise an error if different directories
    try:
        prj.add_custom_hdf5_output(dumpname="custom_dump.h5", dumpdir="alt_hdf5", nsteps=1,
                                   dump_bugs=BugDumpSpec("all"),
                                   dump_chems=ChemDumpSpec("conc"))
    except Exception as e:
        pytest.fail(f'Unexpected exception {e}')

def test_assign_groups():
    def setup_local(taxa_filename:str) -> NufebProject:
        prj = NufebProject()
        prj.add_taxon_by_jsonfile(taxa_filename)
        prj.layout_uniform(nbugs=20)
        prj.set_composition({'anammox_two_pathway': '33.333',
                             'denitrifier': '33.33',
                             'imperfect_denitrifier_no': '33.33'})
        prj.distribute_spatially_even()
        return prj

    #error if generate called without assignment
    prj = setup_local(DATA_DIR / "example_nitrogen_anaerobic_taxa.json")
    prj._assign_taxa()
    with pytest.raises(ValueError) as excinfo:
        prj.generate_case()
    assert f'Taxa Group IDs are not defined, please assign using set_taxa_groups' in str(excinfo.value)

    #error if group asigned to non-existent taxon
    prj = setup_local(DATA_DIR / "example_nitrogen_anaerobic_taxa.json")
    prj._assign_taxa()
    prj.set_taxa_groups({'anammox_two_pathway':'4', 'denitrifier':'6', 'imperfect_denitrifier_no':'8',
                         'this_one_dne':'9'})
    with pytest.raises(ValueError) as excinfo:
        prj.generate_case()
    assert f'Group ID assigned to "this_one_dne", but they do no appear in taxa list.' in str(excinfo.value)

    #error if taxon has no group assignment
    prj = setup_local(DATA_DIR / "example_nitrogen_anaerobic_taxa.json")
    prj._assign_taxa()
    prj.set_taxa_groups({'anammox_two_pathway':'4'})
    with pytest.raises(ValueError) as excinfo:
        prj.generate_case()
    assert f'No group ID assigned to "denitrifier, imperfect_denitrifier_no"' in str(excinfo.value)

    #happy path 1
    prj = setup_local(DATA_DIR / "example_nitrogen_anaerobic_taxa.json")
    prj._assign_taxa()
    prj.set_taxa_groups({'anammox_two_pathway': '4', 'denitrifier': '6', 'imperfect_denitrifier_no': '8'})
    try:
        prj.generate_case()
    except Exception as e:
        pytest.fail(f'Unexpected exception {e}')

    #happy path 2
    prj = setup_local(DATA_DIR / "example_nitrogen_anaerobic_taxa.json")
    prj._assign_taxa()
    prj.set_taxa_groups({'anammox_two_pathway': '7', 'denitrifier': '6', 'imperfect_denitrifier_no': '8'})
    try:
        prj.generate_case()
    except Exception as e:
        pytest.fail(f'Unexpected exception {e}')

    #should account for lysis groups
    prj = setup_local(DATA_DIR / "example_nitrogen_anaerobic_taxa.json")
    prj.add_lysis_group_by_json('vuln_intoxicated',
                                {'name': 'vuln_intoxicated', 'releases': 'sub', 'rate': '2e-3', 'percent': '0.2'})
    prj._assign_taxa()
    prj.set_taxa_groups({'anammox_two_pathway': '7', 'denitrifier': '6', 'imperfect_denitrifier_no': '8'})
    with pytest.raises(ValueError) as excinfo:
        prj.generate_case()
    assert f'No group ID assigned to "vuln_intoxicated"' in str(excinfo.value)

    #happy path 3 with lysis group
    prj = setup_local(DATA_DIR / "example_nitrogen_anaerobic_taxa.json")
    prj._assign_taxa()
    prj.set_taxa_groups({'anammox_two_pathway': '7', 'denitrifier': '6', 'imperfect_denitrifier_no': '8', 'vuln_intoxicated':'7'})
    # adding after the fact intentionally here to exercise order-independence
    prj.add_lysis_group_by_json('vuln_intoxicated',
                                {'name': 'vuln_intoxicated', 'releases': 'sub', 'rate': '2e-3', 'percent': '0.2'})
    prj.set_track_abs()
    try:
        prj.generate_case()
    except Exception as e:
        pytest.fail(f'Unexpected exception {e}')