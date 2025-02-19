import pytest
import polars as pl
from nufebmgr.DumpTools import DumpFile
import pyarrow.parquet as pq
from pyarrow import fs
import numpy as np
import numpy.testing as npt
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

def test_births():
    expected_blank = np.array([], dtype=int)
    expected_births5_dump1 = np.loadtxt("data/birth5_dump1.txt", dtype=int)
    with DumpFile("data/dump_1.h5") as dump:
        result = dump.births(2)
        npt.assert_equal(expected_blank, result)

        result = dump.births(5)
        npt.assert_equal(expected_births5_dump1,result)

def test_deaths():
    expected_blank = np.array([], dtype=int)
    expected_births5_dump1 = np.loadtxt("data/birth5_dump1.txt", dtype=int)
    with DumpFile("data/t6ss_1.h5") as dump:
        result = dump.deaths(2)
        npt.assert_equal(expected_blank, result)

        result = dump.deaths(12)
        npt.assert_equal(np.array([14], dtype=int),result)

        result = dump.deaths(13)
        npt.assert_equal(expected_blank,result)

        result = dump.deaths(16)
        npt.assert_equal(np.array([20, 125, 152], dtype=int),result)








