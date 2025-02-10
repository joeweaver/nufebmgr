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
    def proportional_strips(self, taxa_names, proportions, cut_len, cut_dir: Literal['x', 'y'], noise=0):
        allowed_dirs = {"x", "y"}
        if cut_dir not in allowed_dirs:
            raise ValueError(f"Invalid cut_dir: {cut_dir}. Must be one of {allowed_dirs}.")

        df_locs = pd.DataFrame([asdict(bug_loc) for bug_loc in self.bug_locs])

        # Extract labels and proportions (order is preserved in Python 3.7+)
        labels = list(proportions.keys())
        bin_proportions = list(map(float,proportions.values()))
        bin_proportions = bin_proportions/np.sum(bin_proportions)

        quantiles = [0.0]
        csum = 0
        for proportion in bin_proportions:
            csum = csum + proportion
            quantiles.append(csum*cut_len)

        # Convert proportions to cumulative quantiles
        #quantiles =  [0,0.2*120e-6,0.5*120e-6,120e-6]#bin_proportions / np.cumsum(bin_proportions)
        #quantiles = bin_proportions / np.cumsum(bin_proportions)
        # Use qcut to bin the data
        #df_locs['taxon_name'] = pd.qcut(df_locs[cut_dir], q=quantiles, labels=labels)

        df_locs['taxon_name'] = pd.cut(df_locs[cut_dir], bins=quantiles, labels=taxa_names)

        if noise != 0:
            df_locs = self._shuffle_assignments(df_locs, taxa_names, noise)
        strips =[]
        for row in df_locs.iterrows():
            strips.append(BugPos(x=row[1]["x"],y=row[1]["y"],taxon_name=row[1]["taxon_name"]))

        return strips

    def even_strips(self, taxa_names, cut_dir: Literal['x', 'y'], noise=0):
        allowed_dirs = {"x", "y"}
        if cut_dir not in allowed_dirs:
            raise ValueError(f"Invalid cut_dir: {cut_dir}. Must be one of {allowed_dirs}.")

        df_locs = pd.DataFrame([asdict(bug_loc) for bug_loc in self.bug_locs])


        df_locs['taxon_name'] = pd.cut(df_locs[cut_dir], bins=len(taxa_names), labels=taxa_names)

        if noise != 0:
            df_locs = self._shuffle_assignments(df_locs, taxa_names, noise)
        strips =[]
        for row in df_locs.iterrows():
            strips.append(BugPos(x=row[1]["x"],y=row[1]["y"],taxon_name=row[1]["taxon_name"]))

        return strips

    def _shuffle_assignments(self,locs, taxa_names, percent_shuffle):
        # Shuffle the DataFrame just a little bit
        shuffle_fraction = percent_shuffle/100

        # get a random subset to reassign
        to_reassign = locs.sample(frac=shuffle_fraction)
        nonshuffled = locs.drop(to_reassign.index)

        to_reassign_exploded = to_reassign.explode('taxon_name')

        total_counts = to_reassign_exploded['taxon_name'].value_counts()
        absolute_abundances = total_counts.to_dict()
        relative_abundances = (total_counts / total_counts.sum()).to_dict()

        reassignments = [nonshuffled]
        potentials = to_reassign
        for taxon in taxa_names:
            reassigned = potentials.sample(absolute_abundances[taxon])
            reassigned['taxon_name'] = taxon
            reassignments.append(reassigned)
            potentials = potentials.drop(reassigned.index)

        return pd.concat(reassignments)
