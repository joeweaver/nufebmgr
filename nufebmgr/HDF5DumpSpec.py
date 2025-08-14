from dataclasses import dataclass
from typing import List
@dataclass
class HDF5DumpSpec:
    """
    Specification for a single HDF5 output dump.

    Attributes
    ----------
    dumpname : str
        The name of the HDF5 file to create.
    dumpdir : str
        Directory where the HDF5 file will be written.
    nsteps : int
        Number of simulation steps between consecutive dumps.
    dump_bugs : List[str]
        List of bug-related data fields to include in the dump.
    dump_chems : List[str]
        List of chemical-related data fields to include in the dump.
    """
    dumpname: str
    dumpdir: str
    nsteps: int
    dump_bugs: List[str]
    dump_chems: List[str]