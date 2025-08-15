from jinja2 import Template
from pathlib import Path
from .SimulationBox import SimulationBox
from .Substrate import Substrate
from typing import Any, Dict, List
from .default_inputscript import DEFAULT_INPUTSCRIPT
from .inputscript_template import TEMPLATE_STR
from .HDF5DumpSpec import HDF5DumpSpec
import copy

class InputScriptBuilder:
    def __init__(self):
        self.config_vals = copy.deepcopy(DEFAULT_INPUTSCRIPT)
        self.group_assignments = {}

    def clear_bug_groups(self, keep_dead=True):
        if(keep_dead):
            self.config_vals['microbes_and_groups'][0]['bug_groups'] = [{"name": "group", "group_name": "dead", 'param': "empty", 'comment': "# Dead cells"}]
        else:
            self.config_vals['microbes_and_groups'][0]['bug_groups'] = []

    def build_bug_groups(self, group_assignemnts, active_taxa, lysis_groups, keep_dead=True):
        self.group_assignments = group_assignemnts
        self.clear_bug_groups(keep_dead)
        all_groups = {**active_taxa, **lysis_groups}
        for k, v in group_assignemnts.items():
            entry = {'name': 'group', 'group_name': k, 'param': 'type', 'group_num': v}
            if 'description' in all_groups[k]:
                entry['comment'] = f'# {all_groups[k]["description"]}'
            self.config_vals['microbes_and_groups'][0]['bug_groups'].append(entry)

    def build_substrate_grid(self, substrates, simbox, boundary_scenario, forced_size=None):
        if forced_size is not None:
            if simbox.xlen % forced_size == simbox.ylen % forced_size == simbox.zlen % forced_size == 0:
                grid_size = f'{forced_size}e-6'
            else:
                raise ValueError(f'Grid size was explicity set to {forced_size}, this does not fit a simulation of dimensions {simbox.dim_string()}')
        else:
            grid_size = f'{self._pick_grid_size(simbox)}e-6'
        contents = self.config_vals['mesh_grid_and_substrates'][0]['content']
        new_contents = []
        for content in contents:
            if content['name'] != 'grid_modify':
                if content['name'] != 'grid_style':
                    new_contents.append(content)
                else:
                    grid_style_dict = {'name': 'grid_style', 'loc': 'nufeb/chemostat', 'nsubs': len(substrates)}
                    for i,substrate in enumerate(substrates):
                        grid_style_dict[f's{i}'] = substrate
                    grid_style_dict['grid_cell'] = f'{grid_size}'
                    new_contents.append(grid_style_dict)

        for substrate in substrates.values():
            new_contents.append(self._build_grid_modify_dict(substrate, boundary_scenario))
        self.config_vals['mesh_grid_and_substrates'][0]['content'] = new_contents


    def limit_biofilm_height(self, max_height):
         self.config_vals['system_settings'][0]['content'].append({
                                                                 'name': 'region',
                                                                 'region_id': 'biofilm_height_limiter',
                                                                 'shape': 'block',
                                                                 'dims': f'INF INF INF INF {max_height}e-6 INF',
                                                                 'comment': '# Limit biofilm max height'
                                                                 })

         self.config_vals['biological_processes'][0]['death'].append({
                                                                 'name': 'fix',
                                                                 'fix_name': 'rem_tall',
                                                                 'fig_group': 'all',
                                                                 'fix_loc': 'evaporate',
                                                                 'comment': '# Limit biofilm max height',
                                                                 'p1': '1',
                                                                 'p2': '1000000',
                                                                 'p3': 'biofilm_height_limiter', 'seed': '1701',
                                                                 'comment': ''
                                                                 })


    def track_percent_biomass(self, s: SimulationBox):
        track_dicts = [{'name':'# Compute the volume of all bacteria '},
                       {'name':'compute', 'vname': 'biomass_vol', 'group': 'all', 'loc':'nufeb/volume'},
                       {'name': '# Convert to a percent of simulation volume '},
                       {'name': 'variable', 'vname':'sim_vol', 'op':'equal', 'calc':f'"{s.volume()}"'},
                       {'name': 'variable', 'vname': 'biomass_pct', 'op': 'equal', 'calc': '"c_biomass_vol/v_sim_vol"'}
                       ]

        self.config_vals['computation_output'][0]['percent_biomass']=track_dicts
        pass

    def end_on_biomass(self,percent):
        float_percent = percent/100
        halt_dict = {'name':'fix', 'id': 'halt_vol', 'group':'all', 'cmd':'halt',
                     'N-check': '1', 'comparison': f'v_biomass_pct > {float_percent}',
                     'action':'error', 'error-type':'soft', 'comment':'# end at percent biomass vol'}
        self.config_vals['run'][0]['content'].insert(1,halt_dict)
        pass

        # "run": [
        #     {"title": "#----Run----#",
        #      "content": [{'name': 'run with timestep for physical (pairdt) and chemical (diffdt) processes'},
        #                  {"name": "run_style", "run_type": "nufeb", 'p1': 'diffdt', 'v1': '1e-4',
        #                   'p2': 'difftol', 'v2': '1e-6', 'p3': 'pairdt', 'v3': '1e-2',
        #                   'p4': 'pairtol', 'v4': '1', 'p5': 'screen', 'v5': 'no',
        #                   'p6': 'pairmax', 'v6': '1000', 'p7': 'diffmax', 'v7': '5000',
        #                   'comment': ''},
        #                  {"name": "timestep", "val": "900",
        #                   'comment': '# define biological timestep (900s = 15 min)'},
        #                  {"name": "run", "val": "86400",
        #                   'comment': '# run duration (24h)'}
        #                  ]
        #      }
        # ]
    def _pick_grid_size(self,s:SimulationBox):
        possible_sizes = [2.0,1.5]
        for size in possible_sizes:
            if s.xlen%size==s.ylen%size==s.zlen%size==0:
                return size

        # we'd prefer slightly larger, but return a grid size of 1 if needed
        size = 1.0
        if s.xlen % size == s.ylen % size == s.zlen % size == 0:
            return size
        raise ValueError(f'No valid grid size (option 1, 1.5, 2 microns) for a simulation of dimensions {s.dim_string()}')


    def clear_growth_strategy(self):
        self.config_vals['biological_processes'][0]['growth'] = []

    def _build_grid_modify_dict(self, s: Substrate, boundary_scenario: str) -> dict[str, Any]:
        """Build the 'grid_modify' dictionary for a given substrate and boundary scenario.

            This is an internal helper method and not part of the public API.

            Boundary scenarios:
                - 'bioreactor': open top, periodic sides, no-flux bottom
                - 'microwell': fully closed
                - 'floating': fully open
                - 'agar': open bottom, periodic sides, no-flux top

            Args:
                s: Substrate dataclass
                boundary_scenario: One of {'bioreactor', 'microwell', 'floating', 'agar'}.

            Returns:
                A dictionary of grid modification parameters, including boundaries,
                substrate name, concentrations, and molecular weight if provided.

            Raises:
                ValueError: If `boundary_scenario` is not one of the allowed values.
            """
        #define boundary conditions for different scenarios
        #bioreactor has open top, microwell is fully closed, floating is fully open, agar has open bottom
        bs={'bioreactor': {'x_boundaries': 'pp', 'y_boundaries': 'pp', 'z_boundaries': 'nd'},
            'microwell': {'x_boundaries': 'nn', 'y_boundaries': 'nn', 'z_boundaries': 'nn'},
            'floating': {'x_boundaries': 'dd', 'y_boundaries': 'dd', 'z_boundaries': 'dd'},
            'agar': {'x_boundaries': 'pp', 'y_boundaries': 'pp', 'z_boundaries': 'dn'}}

        if boundary_scenario not in bs:
            raise ValueError(f"Invalid boundary_scenario: {boundary_scenario!r}. "
                             f"Must be one of {set(bs)}.")

        gm_dict = {'name': 'grid_modify',
                   'action': 'set',
                   'substrate': s.name,
                   'xbound': bs[boundary_scenario]['x_boundaries'],
                   'ybound': bs[boundary_scenario]['y_boundaries'],
                   'zbound': bs[boundary_scenario]['z_boundaries'],
                   'init_conc': s.init_concentration,
                   'bulk-kw': 'bulk',
                   'bulkd_conc': s.bulk_concentration}

        if s.molecular_weight:
            gm_dict['mw-kw'] = 'mw'
            gm_dict['mw'] = s.molecular_weight

        return gm_dict


    def build_growth_strategy(self, active_taxa):
        self.clear_growth_strategy()
        self.config_vals['biological_processes'][0]['growth'].append({'name':'Growth Strategies'})
        for k,v in self.group_assignments.items():
            if k not in active_taxa:
                continue
            if 'growth_strategy' not in active_taxa[k]:
                raise KeyError(f"Growth strategy is not defined for taxon: {k}")
            # TODO need to abstract this better. Overarching goal is that new fixes and taxa jsons shouldn't necessitate touching this
            if active_taxa[k]['growth_strategy']['name'] == 'growth_het':
                entry ={}
                entry["name"] = 'fix'
                entry["ID"] = f'growth_{k}'
                entry['group-ID'] = k
                entry['fix_loc'] = "nufeb/growth/het"
                for kj,vj in active_taxa[k]['growth_strategy'].items():
                    if kj != 'name':
                        entry[kj] = vj

                if 'description' in active_taxa[k]:
                    entry['comment'] = active_taxa[k]['description']
                else:
                    entry['comment'] = ''
                self.config_vals['biological_processes'][0]['growth'].append(entry)
            elif active_taxa[k]['growth_strategy']['name'] == 'growth_denit':
                entry = {}
                entry["name"] = '\nfix'
                entry["ID"] = f'growth_{k}'
                entry['group-ID'] = k
                entry['fix_loc'] = "nufeb/growth/denit &\n"
                # TODO essentially preprocessing to make life easier in jinja, look into better ways
                growth_strat = active_taxa[k]["growth_strategy"]
                entry[growth_strat['sub-ID']] = self._composite_kvp(growth_strat, growth_strat['sub-ID'],
                                                                    ['K_s1','K_s2','K_s3','K_s4','K_s5'])
                entry[growth_strat['o2-ID']] = self._composite_kvp(growth_strat, growth_strat['o2-ID'],
                                                                    ['K_oh1', 'K_oh2', 'K_oh3', 'K_oh4', 'K_oh5'])
                entry[growth_strat['no3-ID']] = self._composite_kvp(growth_strat, growth_strat['no3-ID'],  ['K_no3'])
                entry[growth_strat['no2-ID']] = self._composite_kvp(growth_strat, growth_strat['no2-ID'], ['K_no2'])
                entry[growth_strat['n2o-ID']] = self._composite_kvp(growth_strat, growth_strat['n2o-ID'], ['K_n2o'])
                entry[growth_strat['no-ID']] = self._composite_kvp(growth_strat,growth_strat['no-ID'],
                                                                   ['K_no', 'K_I3no', 'K_I4no', 'K_I5no'])
                entry[growth_strat['nh-ID']] = self._composite_kvp(growth_strat, growth_strat['nh-ID'], ['inxb'])
                # the label bit is a bit ridiculous, need to fix when streamlining how fixes, taxon library,
                # and templating work together
                entry['growth_label'] = '\tgrowth'
                entry['growth'] = growth_strat['mu_max']
                entry['yield_label'] = 'yield'
                entry['yield'] = growth_strat['yield']
                entry['decay_label'] = 'decay'
                entry['decay'] = str(growth_strat['decay']) + ' &\n'
                entry['eta_g2_label'] = '\teta_g2'
                entry['eta_g2'] = growth_strat['eta_g2']
                entry['eta_g3_label'] = 'eta_g3'
                entry['eta_g3'] = growth_strat['eta_g3']
                entry['eta_g4_label'] = 'eta_g4'
                entry['eta_g4'] = growth_strat['eta_g4']
                entry['eta_g5_label'] = 'eta_g5'
                entry['eta_g5'] = growth_strat['eta_g5']
                entry['eta_Y_label'] = 'eta_Y'
                entry['eta_Y'] = growth_strat['eta_Y']
                self.config_vals['biological_processes'][0]['growth'].append(entry)
            elif active_taxa[k]['growth_strategy']['name'] == 'growth_imperf_denit_no':
                entry = {}
                entry["name"] = '\nfix'
                entry["ID"] = 'growth_imperf_denit_no'
                entry['group-ID'] = k
                entry['fix_loc'] = "nufeb/growth/imperf_denit_nitric_oxide &\n"
                # TODO essentially preprocessing to make life easier in jinja, look into better ways
                growth_strat = active_taxa[k]["growth_strategy"]
                entry[growth_strat['sub-ID']] = self._composite_kvp(growth_strat, growth_strat['sub-ID'],
                                                                    ['K_s1','K_s2','K_s3'])
                entry[growth_strat['o2-ID']] = self._composite_kvp(growth_strat, growth_strat['o2-ID'],
                                                                    ['K_oh1', 'K_oh2', 'K_oh3'])
                entry[growth_strat['no3-ID']] = self._composite_kvp(growth_strat, growth_strat['no3-ID'],  ['K_no3'])
                entry[growth_strat['no2-ID']] = self._composite_kvp(growth_strat, growth_strat['no2-ID'], ['K_no2'])
                entry[growth_strat['no-ID']] = self._composite_kvp(growth_strat,growth_strat['no-ID'], ['K_I3no'])
                entry[growth_strat['nh-ID']] = self._composite_kvp(growth_strat, growth_strat['nh-ID'], ['inxb'])
                # the label bit is a bit ridiculous, need to fix when streamlining how fixes, taxon library,
                # and templating work together
                entry['growth_label'] = '\tgrowth'
                entry['growth'] = growth_strat['mu_max']
                entry['yield_label'] = 'yield'
                entry['yield'] = growth_strat['yield']
                entry['decay_label'] = 'decay'
                entry['decay'] = str(growth_strat['decay']) + ' &\n'
                entry['eta_g2_label'] = '\teta_g2'
                entry['eta_g2'] = growth_strat['eta_g2']
                entry['eta_g3_label'] = 'eta_g3'
                entry['eta_g3'] = growth_strat['eta_g3']
                entry['eta_Y_label'] = 'eta_Y'
                entry['eta_Y'] = growth_strat['eta_Y']
                self.config_vals['biological_processes'][0]['growth'].append(entry)
            elif active_taxa[k]['growth_strategy']['name'] == 'growth_anammox_two_pathway':
                entry = {}
                entry["name"] = '\nfix'
                entry["ID"] = 'growth_anammox_two_pathway'
                entry['group-ID'] = k
                entry['fix_loc'] = "nufeb/growth/anammox_two_pathway &\n"
                # TODO essentially preprocessing to make life easier in jinja, look into better ways
                growth_strat = active_taxa[k]["growth_strategy"]
                entry[growth_strat['o2-ID']] = self._composite_kvp(growth_strat, growth_strat['o2-ID'], ['K_oh_an'])
                entry[growth_strat['no2-ID']] = self._composite_kvp(growth_strat, growth_strat['no2-ID'], ['K_no2_an'])
                entry[growth_strat['no-ID']] = self._composite_kvp(growth_strat, growth_strat['no-ID'], ['K_no_an'])
                entry[growth_strat['nh-ID']] = self._composite_kvp(growth_strat, growth_strat['nh-ID'],
                                                                   ['K_nh_an', 'inxb_an'])
                entry['no3-ID'] = '\tno3 &\n'
                # the label bit is a bit ridiculous, need to fix when streamlining how fixes, taxon library,
                # and templating work together
                entry['growth_label'] = '\tgrowth'
                entry['growth'] = growth_strat['mu_max']
                entry['yield_label'] = 'yield'
                entry['yield'] = growth_strat['yield']
                entry['decay_label'] = 'decay'
                entry['decay'] = str(growth_strat['decay']) + ' &\n'
                entry['eta_I_AN_label'] = '\teta_I_an'
                entry['eta_I_AN'] = growth_strat['eta_I_an']
                entry['eta_S_AN_label'] = 'eta_S_an'
                entry['eta_S_an'] = growth_strat['eta_S_an']
                self.config_vals['biological_processes'][0]['growth'].append(entry)
            else:
                raise KeyError(f"Taxon {k} has unrecognized growth strategy: {active_taxa[k]['growth_strategy']['name'] }")

    def _composite_kvp(self,db: Dict[str,str], main_key:str, value_keys: List[str]) -> str:
        """
        Utility function to take a growth strategy dictionary where multiple parameters follow a named entry and
        combine them into one flat string, prepended by a tab and ending with a LAMMPS line continuation and newline.

        This is a somewhat stopgap measure as the original jinja template just takes key value pairs. It also adds a
        bit of pretty print formatting with tabs and line continuations, which is nice even for single pairs.

        Example: Assuming a growth strategy dict contains something like a sub-ID:'sub' followed by multiple K_S
        constants.  {'sub_id':sub, 'ks1':0.2 'ks2':03}
        If the fix does not have named parameters, it might 'want' to be expressed as '\t sub 0.2 0.3 &\n'.
        This function helps with that formatting.
        Here, db would be active_taxa[k]["growth_strategy"]
        and value_keys would be ['ks1', 'ks2']

        :param db: The growth strategy dict. e.g active_taxa[k]["growth_strategy"]
        :param main_key: The label to use at the start
        :param value_keys: The keys with values to be concatenated
        :return: The formatted string
        """
        return (f'\t{main_key} ' + ' '.join(str(db[key]) for key in value_keys) +' &\n')

    def clear_division(self):
        self.config_vals['biological_processes'][0]['division'] = []

    def build_division(self, active_taxa, seed):
        self.clear_division()
        self.config_vals['biological_processes'][0]['division'].append({'name':'Division'})

        for k,v in self.group_assignments.items():
            if k not in active_taxa: # TODO should this be an error?
                continue
            if 'division_strategy' not in active_taxa[k]:
                raise KeyError(f"Division strategy is not defined for taxon: {k}")
            if active_taxa[k]['division_strategy']['name'] == 'divide_coccus':
                entry ={}
                entry["name"] = 'fix'
                entry["ID"] = f'division_{k}'
                entry['group-ID'] = k
                entry['fix_loc'] = "nufeb/division/coccus"
                for kj,vj in active_taxa[k]['division_strategy'].items():
                    if kj not in ['name','eps_dens', 'version']:
                        entry[kj] = vj

                entry['seed']= seed
                if  'eps_dens' in active_taxa[k]['division_strategy']:
                    entry['eps_dens_name'] = 'eps_dens'
                    entry['eps_dens_val'] = entry['eps_dens']

                if 'description' in active_taxa[k]:
                    entry['comment'] = active_taxa[k]['description']
                else:
                    entry['comment'] = ''
                self.config_vals['biological_processes'][0]['division'].append(entry)
            else:
                raise KeyError(f"Taxon {k} has unrecognized division strategy: {active_taxa[k]['division_strategy']['name'] }")

    def add_hdf5_output(self, hdf5_dump_specs: List[HDF5DumpSpec]):
        if hdf5_dump_specs == []:
            self.config_vals['computation_output'][0]['hdf5_output'] = []
            return

        self.config_vals['computation_output'][0]['hdf5_output'] = [
            {'linetype': 'subsection',
             'title': 'HDF5 output, efficient binary format, for storing many atom properties'},
            {'linetype': 'comment', 'title': 'NOTE: requires NUFEB built with HDF5 option'},
            {'linetype': 'subsection', 'title': 'Create directory(s) for dump'}]

        # create a list of unique dumpdirs
        dump_dirs = list({Path(hds.dumpdir) for hds in hdf5_dump_specs})
        for dump_dir in sorted(dump_dirs):
            self.config_vals['computation_output'][0]['hdf5_output'].append(
              {'linetype': 'command', 'name': 'shell', 'command': f'mkdir {dump_dir}'})

        self.config_vals['computation_output'][0]['hdf5_output'].append({'linetype': 'subsection', 'title': 'Dump specifications'})

        for i, hds in enumerate(hdf5_dump_specs):
            save_loc = Path(hds.dumpdir) / hds.dumpname
            dump_vars = ' '.join(hds.dump_bugs + hds.dump_chems)
            self.config_vals['computation_output'][0]['hdf5_output'].append(
              {'linetype': 'command', 'name': 'dump', 'dumpname': f'du_hdf5_{i}', 'group': 'all',
               'format': 'nufeb/hdf5', 'every_n': f'{hds.nsteps}',
                  'loc': f'{save_loc}', 'dumpvars': f'{dump_vars}'})


    def add_vtk_output(self):
        self.config_vals['computation_output'][0]['vtk_output'] = [
            {'name': 'VTK output, useful for paraview visualizations'},
            {'name': 'requires NUFEB built with VTK option'},
            {'name': 'shell', 'command': 'mkdir vtk', 'comment': '#Create directory for dump'},
            {'name': 'dump', 'dumpname': 'du1', 'group': 'all', 'format': 'vtk', 'p1': '1',
             'loc': 'vtk/dump*.vtu',
             'd1': 'id', 'd2': 'type', 'd3': 'diameter', 'comment': ''},
            {'name': 'dump', 'dumpname': 'du2', 'group': 'all', 'format': 'grid/vtk', 'p1': '10',
             'loc': 'vtk/dump_%_*.vti',
             'd1': 'con', 'd2': 'rea', 'd3': 'den', 'd4': 'gro', 'comment': ''},
        ]

    def add_thermo_output(self,track_abs,timestep):
        self.config_vals['computation_output'][0]['thermo_output'] = []
        self.config_vals['computation_output'][0]['thermo_output'].append({'name': 'Output to screen'})
        thermo_style = {'name': 'thermo_style',
                          'args': 'custom',
                          'custom_key1': 'step',
                          'custom_key2': 'cpu',
                          'col0':'atoms'}
        if track_abs:
            for k,v in self.group_assignments.items():
                thermo_style[f'abs_var_{k}'] = f'v_n_{k}'
                thermo_style[f'rel_var_{k}'] = f'v_ra_{k}'
        self.config_vals['computation_output'][0]['thermo_output'].append(thermo_style)
        self.config_vals['computation_output'][0]['thermo_output'].append({'name':'thermo', 'step':timestep})

    def enable_csv_output(self,tracking_abs,tracking_biomass_pct):
        csv_vars = ['current_step']
        csv_header=['step']
        if tracking_abs:
            csv_vars.append('biomass_pct')
            csv_header.append('percent_vol_biomass')
        if tracking_abs:
            for k, v in self.group_assignments.items():
                # abs abundance
                csv_vars.append(f'n_{k}')
                csv_header.append(f'{k}_abundance')
                # rel abundance
                csv_vars.append(f'ra_{k}')
                csv_header.append(f'{k}_relative_abundance')

        var_strings=[]
        for csv_var in csv_vars:
            var_strings.append(f'${{{csv_var}}}')
        var_string = ','.join(var_strings)
        header_string = ','.join(csv_header)
        csv_dict = [{'name':'# Save some results as CSV '},
                    {'name': 'variable', 'vname':'current_step', 'op':'equal', 'c':'step'},
                       {'name':'fix', 'vname': 'volcsv', 'group': 'all', 'loc':'print',
                        'time':'1', 'vars':f'"{var_string}"', 'screen': 'screen no',
                        'v1':'file', 'file': 'output.csv', 't':'title', 'header':f'"{header_string}"'}
                       ]

        self.config_vals['computation_output'][0]['csv_output']=csv_dict

            # fix
            # volcsv
            # all
            # print
            # 100
            # "${step},${simvol},${fillpct},${bug1_relab},${bug2_relab}" &
            # screen
            # no &
            # file
            # "cell_rel_volumes.csv" &
            # title
            # "step,simulation_volume,fill_percent,bug1_relab,bug2_relab"

    def add_abs_vars(self):
        self.config_vals['computation_output'][0]['ab_track'] = []
        self.config_vals['computation_output'][0]['ab_track'].append({'name': 'Tracking abundances'})
        self.config_vals['computation_output'][0]['ab_track'].append({'name':'variable',
                                                                      'varname':'n_all',
                                                                      'op': 'equal',
                                                                      'expression':'"count(all)"'})
        for k,v in self.group_assignments.items():
            # abs abundance
            entry = {}
            entry['name'] = 'variable'
            entry['varname'] = f'n_{k}'
            entry['op'] = 'equal'
            entry['expression'] = f'"count({k})"'
            self.config_vals['computation_output'][0]['ab_track'].append(entry)
            # rel abundance
            entry = {}
            entry['name'] = 'variable'
            entry['varname'] = f'ra_{k}'
            entry['op'] = 'equal'
            entry['expression'] = f'"v_n_{k}/v_n_all"'

            self.config_vals['computation_output'][0]['ab_track'].append(entry)

    def build_lysis(self, lysis_groups):
        if not lysis_groups:
            return
        self.config_vals['biological_processes'][0]['lysis'] = []
        self.config_vals['biological_processes'][0]['lysis'].append({'name': 'Lysis, requires make yes-T6ss build'})

        for k, v in self.group_assignments.items():
            if k in lysis_groups:
                entry = {}
                entry["name"] = 'fix'
                entry["ID"] = 'lysis_' + k
                entry['group-ID'] = k
                entry['fix_loc'] = "nufeb/T6SS/lysis"
                for kj, vj in lysis_groups[k].items():
                    if kj != 'name':
                        entry[kj] = vj
                self.config_vals['biological_processes'][0]['lysis'].append(entry)

    def build_run(self, runtime, biostep_s):
        content = self.config_vals['run'][0]['content']
        bio_timestep=biostep_s
        for i,item in enumerate(content):
            if item['name']=='run':
                run_hours = (runtime*bio_timestep)/(60*60)
                new_item = {"name": "run", "val": f'{runtime}', 'comment': f'# run duration ({run_hours} H)'}
                content[i]=new_item
            if item['name']=='timestep':
                new_item = {"name": "timestep", "val": f'{bio_timestep}', 'comment': f'# seconds per timestep)'}
                content[i]=new_item

    def build_post_physical(self, elastic_bl:float) -> None:
        """Build the dictionary used by the jinja template to populate the postphysical effects section.
        Currently only adds a blayer fix.

        Note, this only makes sense right now for the bioreactor boundary scenario and uses the implicit 'z is up' convention.

        Calling the function clears anything already in the relevant data structure.
        i.e., it clobbers
            self.config_vals['post_physical_processes'][0]['boundary_layer']

        Args:
            elastic_bl: distance (in microns) above the biofilm to place the boundary layer.

        Returns:
            None

        Side Effects:
            Updates:
                self.config_vals['post_physical_processes'][0]['boundary_layer']
        """
        self.config_vals['post_physical_processes'][0]['boundary_layer'] = []
        self.config_vals['post_physical_processes'][0]['boundary_layer'].append({'title': 'Elastic boundary layer above biofilm surface'})
        fix_dict = {'name': 'fix',
                     'fix_name': 'blayer',
                     'group': 'all',
                     'fix_loc': 'nufeb/boundary_layer',
                     'distance_m': f'{elastic_bl}e-6',
                     'comment': ''}
        self.config_vals['post_physical_processes'][0]['boundary_layer'].append(fix_dict)

    def build_diffusion(self, substrates: Dict[str,Substrate]) -> None:
        """Build the dictionary used by the jinja template to add substrate diffusion to the inputscript

            Calling the function clears anything already in the relevant data structure.
            i.e., it clobbers
                self.config_vals['chemical_processes'][0]['diffusion_coefficients']
                self.config_vals['chemical_processes'][0]['diffusion_biofilm_ratios']

            Args:
                substrates: Dictionary of substrates, where keys are substrate names and values are instances of Substrate

            Returns:
                None

            Side Effects:
                Mutates in place:
                    self.config_vals['chemical_processes'][0]['diffusion_coefficients']
                    self.config_vals['chemical_processes'][0]['diffusion_biofilm_ratios']
            """
        self.config_vals['chemical_processes'][0]['diffusion_biofilm_ratios'] = []
        self.config_vals['chemical_processes'][0]['diffusion_biofilm_ratios'].append(
            {'title': 'Ratio of diffusion in biofilm as compared to water'})
        for substrate in substrates.values():
            diff_dict = {'name': 'fix',
                         'fix_name': f'coeff_{substrate.name}',
                         'group': 'all',
                         'fix_loc': 'nufeb/diffusion_coeff',
                         'sub1': substrate.name,
                         'coeff1': substrate.biofilm_diffusion_ratio,
                         'comment': ''}
            self.config_vals['chemical_processes'][0]['diffusion_biofilm_ratios'].append(diff_dict)

        self.config_vals['chemical_processes'][0]['diffusion_coefficients'] = []
        self.config_vals['chemical_processes'][0]['diffusion_coefficients'].append({'title': 'Diffusion in water'})
        for substrate in substrates.values():
            diff_dict = {'name': 'fix',
                         'fix_name': f'diff_{substrate.name}',
                         'group': 'all',
                         'fix_loc': 'nufeb/diffusion_reaction',
                         'sub1': substrate.name,
                         'coeff1': substrate.diffusion_coefficient,
                         'comment': ''}
            self.config_vals['chemical_processes'][0]['diffusion_coefficients'].append(diff_dict)



    def build_t6ss(self, t6ss_attackers,t6ss_vulns,seed):
        if not t6ss_attackers:
            return
        if not t6ss_vulns:
            return

        self.config_vals['biological_processes'][0]['t6ss'] = []
        self.config_vals['biological_processes'][0]['t6ss'].append({'name': 'T6SS, requires make yes-T6ss build'})
        entry = {}
        entry["name"] = 'fix'
        entry["ID"] = 't6ss'
        entry['group-ID'] = 'all'
        entry['fix_loc'] = "nufeb/T6SS/contact"
        entry['seed'] = seed
        entry['n_attackers'] =len(t6ss_attackers)
        effector_ids = {}
        for k,v in t6ss_attackers.items():
            entry[f'attacker_{k}'] = self.group_assignments[k]
            effector_ids[t6ss_attackers[k]['effector']] = len(effector_ids)+1
            entry[f'effector_{k}'] = effector_ids[t6ss_attackers[k]['effector']]
            entry[f'harp_len_{k}'] = t6ss_attackers[k]['harpoon_len']
            entry[f'cooldown{k}'] = t6ss_attackers[k]['cooldown']
        entry['n_vulns'] = len(t6ss_vulns)
        for k,v in t6ss_vulns.items():
            entry[f'vuln_{k}'] = self.group_assignments[k]
            entry[f'vuln_to_effector_{k}'] = effector_ids[t6ss_vulns[k]['effector']]
            entry[f'vuln_{k}_intox_prob'] = t6ss_vulns[k]['prob']
            entry[f'vuln_{k}_to_group'] = t6ss_vulns[k]['to_group']
            entry[f'vuln_{k}_to_group_ID'] = self.group_assignments[t6ss_vulns[k]['to_group']]
        self.config_vals['biological_processes'][0]['t6ss'].append(entry)

    def generate(self):
        template = Template(TEMPLATE_STR)
        config_content = template.render(self.config_vals)
        return config_content