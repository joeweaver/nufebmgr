import pytest
import polars as pl
from nufebmgr.DumpTools import DumpFile
import pyarrow.parquet as pq
from pyarrow import fs
import polars as pl

def test_num_timesteps():

    with DumpFile("data/dump_1.h5") as dump:
        assert 20 == dump.num_timesteps()

    with DumpFile("data/dump_2.h5") as dump:
        assert 25 == dump.num_timesteps()

def test_timesteps():
    with DumpFile("data/dump_1.h5") as dump:
        assert list(range(0,20)) == dump.timesteps()

    with DumpFile("data/dump_2.h5") as dump:
        assert list(range(0, 25)) == dump.timesteps()

def test_population_abs():
    expected = pl.read_csv("data/expected_population_abs_dump_1.csv")
    with DumpFile("data/dump_1.h5") as dump:
        result = dump.population_abs()
        assert expected.equals(result)









