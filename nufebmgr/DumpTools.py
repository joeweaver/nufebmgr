import h5py
import polars
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

        self.df = pl.DataFrame()
        timesteps = np.array(self.dumpfile['/id'])
        for timestep in timesteps:
            ids = np.array(self.dumpfile[f'/id/{timestep}'])
            types = np.array(self.dumpfile[f'/type/{timestep}'])
            times = np.full(ids.shape, int(timestep), dtype=int)
            tdf = polars.DataFrame(np.stack([times,ids,types],axis=1),
                                      schema=['timestep','id','type'])
            self.df = self.df.vstack(tdf)

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

    # TODO could DRY between births() and deaths()
    # TODO also, look into just ingesting the h5 into polars on construction and doing all this with dataframe operations
    #      main argument against that is that the straightforward way would be memory inefficient
    def births_at_time(self, timestep:int) -> List[int]:
        id_now = self.fields_at_time('id', timestep)
        try:
            id_past = self.fields_at_time('id', timestep-1)
            return np.setdiff1d(id_now, id_past)
        except KeyError:
            print(f'Trying to infer births at time {timestep}. It appears data for the immediate previous {timestep-1} does not exist.')


    def deaths_at_time(self, timestep:int) -> List[int]:
        id_now = self.fields_at_time('id', timestep)
        try:
            id_past = self.fields_at_time('id', timestep-1)
            return np.setdiff1d(id_past, id_now)
        except KeyError:
            print(f'Trying to infer births at time {timestep}. It appears data for the immediate previous {timestep-1} does not exist.')

    def births(self, groups:dict=None, as_df=False) -> dict:
        # TODO raise error if dump doesn't consist of consecutive timessteps a la
        # (df['time'].unique().sort().diff().drop_null() == 1).all()
        new_ids_per_timestep = self.df.join(self.df,
                                             left_on=[pl.col("timestep") - 1, "id"],
                                             right_on=["timestep", "id"],
                                             how="anti").filter(pl.col("timestep") > self.df['timestep'].min())
        births = {}
        if groups is None:
            if as_df:
                return new_ids_per_timestep
            else:
                births['all'] = new_ids_per_timestep
                return births

        for group_name,group_types in groups.items():
           births[group_name] = new_ids_per_timestep.filter(pl.col("type").is_in(group_types))
        return births

    def deaths(self, groups:dict=None, as_df=False) -> dict:
        # TODO raise error if dump doesn't consist of consecutive timessteps a la
        # (df['time'].unique().sort().diff().drop_null() == 1).all()

        deaths_per_timestep = (self.df.join(self.df,
                                            left_on=[pl.col("timestep") + 1, "id"],
                                            right_on=["timestep", "id"],  # N
                                            how="anti")
                                            .with_columns((pl.col("timestep") + 1).alias("timestep"))
                                            .filter(pl.col("timestep") < self.df['timestep'].max() + 1))

        deaths = {}
        if groups is None:
            if as_df:
                return deaths_per_timestep
            else:
                deaths['all'] = deaths_per_timestep
                return deaths

        for group_name,group_types in groups.items():
           deaths[group_name] = deaths_per_timestep.filter(pl.col("type").is_in(group_types))
        return deaths

    # def curtis_numbers(self):
    #     births_agg =pl.DataFrame()
    #     for timestep in self.timesteps():
    #         print(timestep)
    #         ids=self.births_at_time(timestep)
    #         bt=self.df.filter(pl.col('timestep')==timestep).filter(pl.col('id').is_in(ids))
    #         births_agg = births_agg.vstack(bt.group_by(pl.col('timestep', 'type')).agg(pl.count().alias('births')))
    #     deaths_agg =pl.DataFrame()
    #     for timestep in self.timesteps():
    #         print(timestep)
    #         ids=self.deaths_at_time(timestep)
    #         bt=self.df.filter(pl.col('timestep')==timestep-1).filter(pl.col('id').is_in(ids))
    #         deaths_agg = deaths_agg.vstack(bt.group_by(pl.col('timestep', 'type')).agg(pl.count().alias('deaths')))
    #
    #     bds = births_agg.join(deaths_agg, on=["timestep", "type"], how='full', coalesce=True).fill_null(0)
    #     taxon_mapping = pl.DataFrame({
    #         "type": [1, 2, 3, 4, 5, 6],  # Known types
    #         "group": ["A", "B", "C", "A", "B", "C"]
    #     })
    #     bds=bds.join(taxon_mapping, on="type", how="left")
    #     aggregated = (
    #         bds.group_by(["timestep", "group"])
    #             .agg([
    #                     pl.col("births").sum(),
    #                     pl.col("deaths").sum()
    #             ])
    #     )
    #     cns = aggregated.with_columns(
    #         pl.when(pl.col("births") > 0)
    #                     .then(pl.col("deaths") / pl.col("births"))
    #                     .otherwise(None)
    #                     .alias("Curtis_number"))
    #     cns.write_csv('curtis.csv')
    #     return



