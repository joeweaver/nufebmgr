import h5py
import polars as pl
import numpy as np
from typing import List


class DumpFile:
    """
    A context manager for dealing with NUFEB HDF5 dumpfiles.

    The focus is high level abstractions for managing output, e.g. 'Get the abundances of all types over time'

    Although high level abstractions are preferable, it should also be possible to 'export' data as tables for custom
    processing.

    Eventually we may want to have something like polars reading via pyarrow loading hdf5, but right now we're doing
    things a bit more simply.

    Attributes:
        dumpfile_name (str): The path to the hdf5 file
        dumpfile (Optional[IO]): The the actual HDF5 file
    """

    def __init__(self,dumpfile_name: str):
        """
        Initialize the DumpFile with the given filename

        :param dumpfile_name (str): Path to the dump file
        """
        self.dumpfile_name = dumpfile_name
        self.dumpfile = None

    def __enter__(self) -> "Dumpfile":
        """
        Open the dumpfile and provide context.

        :return: Instance of itself with an open dump file
        """
        self.dumpfile = h5py.File(self.dumpfile_name, 'r')
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """
        Closes the dump file on exit context
           """
        if self.dumpfile:
            self.dumpfile.close()

    def num_timesteps(self) -> int:
        return max(self._id_list())+1

    def timesteps(self) -> List[int]:
        return sorted(self._id_list())

    def _id_list(self) -> List[int]:
        return list(map(int, self.dumpfile['/id']))
    def population_abs(self):
        return pl.DataFrame([self.unique_types_at_time(time) for time in self.timesteps()])

    def fields_at_time(self,field:str, t:int):
        try:
            return list(self.dumpfile[f'{field}/{t}'])
        except KeyError:
            print(f"Received a key error trying to read field: {field} at time {t} from {self.dumpfile_name}")

    def types_at_time(self,t:int) -> list:
        return self.fields_at_time("type", t)
    def _count_uniques(self,x:list) ->dict:
        uniques, counts = np.unique(x, return_counts=True)
        unique_str = [str(unique) for unique in uniques]
        return dict(zip(unique_str, counts))

    def unique_types_at_time(self, t:int):
        return self._count_uniques(self.types_at_time(t))

    def births(self, timestep:int) -> List[int]:
        id_now = self.fields_at_time('id', timestep)
        try:
            id_past = self.fields_at_time('id', timestep-1)
            return np.setdiff1d(id_now, id_past)
        except KeyError:
            print(f'Trying to infer births at time {timestep}. It appears data for the immediate previous {timestep-1} does not exist.')


