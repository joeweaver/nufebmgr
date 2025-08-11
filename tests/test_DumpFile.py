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
    expected = pl.read_csv("data/expected_births_dump_1.csv")
    expected_just_taxa1 = expected.filter(pl.col("type").is_in([1]))
    expected_just_taxa23 = expected.filter(pl.col("type").is_in([2,3]))
    with DumpFile("data/dump_1.h5") as dump:
        result = dump.births()
        assert expected.equals(result['all'])

        result = dump.births(groups={'taxa1':[1]})
        assert expected_just_taxa1.equals(result['taxa1'])
        assert list(result.keys()) == ['taxa1']

        result = dump.births(groups={'taxa23':[2,3]})
        assert expected_just_taxa23.equals(result['taxa23'])
        assert list(result.keys()) == ['taxa23']

        result = dump.births(groups={'taxa1':[1], 'taxa23':[2,3]})
        assert expected_just_taxa1.equals(result['taxa1'])
        assert expected_just_taxa23.equals(result['taxa23'])
        assert list(result.keys()) == ['taxa1','taxa23']

def test_deaths():
    expected = pl.read_csv("data/expected_deaths_t6ss_1.csv")
    expected_just_taxa4 = expected.filter(pl.col("type").is_in([4]))
    expected_just_taxa5 = expected.filter(pl.col("type").is_in([5]))
    #TODO set upa  system which has deaths for more than 1 type and use that for test data
    with DumpFile("data/t6ss_1.h5") as dump:
        result = dump.deaths()
        assert expected.equals(result['all'])

        result = dump.deaths(groups={'taxa4':[4]})
        assert expected_just_taxa4.equals(result['taxa4'])
        assert list(result.keys()) == ['taxa4']

        result = dump.deaths(groups={'taxa5':[5]})
        assert expected_just_taxa5.equals(result['taxa5'])
        assert list(result.keys()) == ['taxa5']

        result = dump.deaths(groups={'taxa4':[4], 'taxa5':[5]})
        assert expected_just_taxa4.equals(result['taxa4'])
        assert expected_just_taxa5.equals(result['taxa5'])
        assert list(result.keys()) == ['taxa4','taxa5']

def test_births_at_time():
    expected_blank = np.array([], dtype=int)
    expected_births5_dump1 = np.loadtxt("data/birth5_dump1.txt", dtype=int)
    with DumpFile("data/dump_1.h5") as dump:
        result = dump.births_at_time(2)
        npt.assert_equal(expected_blank, result)

        result = dump.births_at_time(5)
        npt.assert_equal(expected_births5_dump1,result)

def test_deaths_at_time():
    expected_blank = np.array([], dtype=int)
    expected_births5_dump1 = np.loadtxt("data/birth5_dump1.txt", dtype=int)
    with DumpFile("data/t6ss_1.h5") as dump:
        result = dump.deaths_at_time(2)
        npt.assert_equal(expected_blank, result)

        result = dump.deaths_at_time(12)
        npt.assert_equal(np.array([14], dtype=int),result)

        result = dump.deaths_at_time(13)
        npt.assert_equal(expected_blank,result)

        result = dump.deaths_at_time(16)
        npt.assert_equal(np.array([20, 125, 152], dtype=int),result)

# def test_curtis_numbers():
#     #TODO set upa  system which has deaths for more than 1 type and use that for test data
#     with DumpFile("data/t6ss_1.h5") as dump:
#         dump.curtis_numbers()
#
#     with DumpFile("/home/joe/professional/software/nufeb/cases/T6SS/force_biofilm/case/rand1701-160wide/hdf5/dump.h5") as dump:
#         dump.curtis_numbers()







