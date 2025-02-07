'''

Encapsulate the nitty-gritty details of assigning taxa to a list of positions

'''

import numpy as np
import pandas as pd
from dataclasses import dataclass, asdict
from typing import Literal
from .BugPos import BugPos

class TaxaAssignmentManager:
    def __init__(self,bug_locs):
        self.bug_locs = bug_locs
        pass

    # assign as evenly spaced vertical or horizontal strips
    def even_strips(self, taxa_names, cut_dir: Literal['x', 'y']):
        allowed_dirs = {"x", "y"}
        if cut_dir not in allowed_dirs:
            raise ValueError(f"Invalid cut_dir: {cut_dir}. Must be one of {allowed_dirs}.")


        df_locs = pd.DataFrame([asdict(bug_loc) for bug_loc in self.bug_locs])


        df_locs['taxon_name'] = pd.cut(df_locs[cut_dir], bins=len(taxa_names), labels=taxa_names)

        strips =[]
        for row in df_locs.iterrows():
            strips.append(BugPos(x=row[1]["x"],y=row[1]["y"],taxon_name=row[1]["taxon_name"]))


        return strips