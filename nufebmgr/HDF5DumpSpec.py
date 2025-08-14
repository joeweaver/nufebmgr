from dataclasses import dataclass
from typing import List
@dataclass
class HDF5DumpSpec:
    dumpname: str
    dumpdir: str
    nsteps: int
    dump_bugs: List[str]
    dump_chems: List[str]