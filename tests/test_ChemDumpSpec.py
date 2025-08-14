import pytest
from nufebmgr.ChemDumpSpec import ChemDumpSpec

def test_default_hdf5_vars():
    cds = ChemDumpSpec()
    to_dump = cds.hdf5_vars()
    assert to_dump == ['conc', 'reac']

def test_all_hdf5_vars():
    cds = ChemDumpSpec("all")
    to_dump = cds.hdf5_vars()
    assert to_dump == ['conc', 'reac']

def test_conc_hdf5_vars():
    cds = ChemDumpSpec("conc")
    to_dump = cds.hdf5_vars()
    assert to_dump == ['conc']

def test_reac_hdf5_vars():
    cds = ChemDumpSpec("reac")
    to_dump = cds.hdf5_vars()
    assert to_dump == ['reac']

def test_custom_hdf5_vars():
    bds = ChemDumpSpec("custom")
    to_dump = bds.hdf5_vars()
    assert to_dump == []

    bds = ChemDumpSpec("custom", [])
    to_dump = bds.hdf5_vars()
    assert to_dump ==  []

    bds = ChemDumpSpec("custom", ['x', 'y', 'z'])
    to_dump = bds.hdf5_vars()
    assert to_dump == ['x', 'y', 'z']

    # explicitly allow anything the user wants in custom
    bds = ChemDumpSpec("custom", ['type', 'id', 'banana'])
    to_dump = bds.hdf5_vars()
    assert to_dump == ['type', 'id', 'banana']

def test_unknown_spec_hdf5_vars():
    valids = ChemDumpSpec._VALID_H5
    with pytest.raises(ValueError) as excinfo:
        bds = ChemDumpSpec("banana")
        to_dump = bds.hdf5_vars()
    assert f'Unknown output specification "banana". Valid values are: {", ".join(valids)}' in str(excinfo.value)

def test_eq():
    bds1 = ChemDumpSpec("custom", ['banana'])
    bds2 = ChemDumpSpec("custom", [ 'banana'])
    bds3 = ChemDumpSpec("custom", ['orange'])

    assert bds1 == bds2
    assert bds1 != bds3
    assert (bds1 == "not ChemDumpSpec") is False