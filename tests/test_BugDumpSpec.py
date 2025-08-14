import pytest
from nufebmgr.BugDumpSpec import BugDumpSpec

def test_default_hdf5_vars():
    bds = BugDumpSpec()
    to_dump = bds.hdf5_vars()
    assert to_dump == ['id', 'type', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'fx', 'fy', 'fz', 'radius']

def test_all_hdf5_vars():
    bds = BugDumpSpec("all")
    to_dump = bds.hdf5_vars()
    assert to_dump == ['id', 'type', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'fx', 'fy', 'fz', 'radius']

def test_location_hdf5_vars():
    bds = BugDumpSpec("location")
    to_dump = bds.hdf5_vars()
    assert to_dump == ['id', 'type', 'x', 'y', 'z']

def test_location_hdf5_vars():
    bds = BugDumpSpec("location_radius")
    to_dump = bds.hdf5_vars()
    assert to_dump == ['id', 'type', 'x', 'y', 'z', 'radius']

def test_custom_hdf5_vars():
    bds = BugDumpSpec("custom")
    to_dump = bds.hdf5_vars()
    assert to_dump == []

    bds = BugDumpSpec("custom", [])
    to_dump = bds.hdf5_vars()
    assert to_dump ==  []

    bds = BugDumpSpec("custom", ['x', 'y', 'z'])
    to_dump = bds.hdf5_vars()
    assert to_dump == ['x', 'y', 'z']

    bds = BugDumpSpec("custom", ['type', 'id', 'fx'])
    to_dump = bds.hdf5_vars()
    assert to_dump == ['type', 'id', 'fx']

    # explicitly allow anything the user wants in custom
    bds = BugDumpSpec("custom", ['type', 'id', 'banana'])
    to_dump = bds.hdf5_vars()
    assert to_dump == ['type', 'id', 'banana']

def test_unknown_spec_hdf5_vars():
    valids = BugDumpSpec._VALID_H5
    with pytest.raises(ValueError) as excinfo:
        bds = BugDumpSpec("banana")
        to_dump = bds.hdf5_vars()
    assert f'Unknown output specification "banana". Valid values are: {", ".join(valids)}' in str(excinfo.value)


def test_eq():
    bds1 = BugDumpSpec("custom", ['type', 'id', 'banana'])
    bds2 = BugDumpSpec("custom", ['type', 'id', 'banana'])
    bds3 = BugDumpSpec("custom", ['type', 'id', 'orange'])

    assert bds1 == bds2
    assert bds1 != bds3
    assert (bds1 == "not BugDumpSpec") is False
