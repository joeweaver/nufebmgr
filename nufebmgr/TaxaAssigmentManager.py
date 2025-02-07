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

        # choose new A's
        potA = to_reassign
        a_reassigned = potA.sample(absolute_abundances['basic_het']) #frac=relative_abundances['basic_het'])
        a_reassigned['taxon_name'] = 'basic_het'


        # choose new B's
        potB = to_reassign.drop(a_reassigned.index)

        b_reassigned = potB.sample(absolute_abundances['small_het'])
        b_reassigned['taxon_name'] = 'small_het'

        potC = to_reassign.drop(a_reassigned.index).drop(b_reassigned.index)
        pass
        n_rows,_ = potC.shape
        if n_rows != absolute_abundances['slow_het']:
            raise "Error"
        c_reassigned = potC.sample(absolute_abundances['slow_het'])  # count_c_species / potC.species.count())
        c_reassigned['taxon_name'] = 'slow_het'
        #  sns.relplot(data=c_reassigned, x='x', y='y', hue='species', aspect=1.61).set(title="new C's")
        #  plt.xlim(0,width)
        # plt.ylim(0,height)
        #  plt.show()

        shuffled_df2 = pd.concat([a_reassigned, b_reassigned, c_reassigned, nonshuffled])
        return shuffled_df2
