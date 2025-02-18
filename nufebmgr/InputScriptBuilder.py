from jinja2 import Template
from datetime import datetime
import math
from functools import reduce
from .SimulationBox import SimulationBox

class InputScriptBuilder:
    # Default configuration data
    DEFAULT_INPUTSCRIPT = {
        "header": ["#----------------------------------------------------------------------#",
                   "#                    NUFEB Simulation                                  #",
                   f"#               Generated on: {datetime.now().strftime('%Y-%m-%d')}   #",
                   "#----------------------------------------------------------------------#"],

        "system_settings": [
            {"title": "#----System Settings----#",
             "content": [{"name": "units", "unit_type": "si", 'comment': "Using si units (m, s, kg)"},
                         {"name": "atom_style", "style": "coccus", 'comment': "Using NUFEB coccus"},
                         {"name": "atom_modify", "method": "map", "struct": "array", 'method2': 'sort', "time": "10",
                          "param": "0", 'comment': 'find atoms using indices, sort every 10 steps'},
                         {"name": "boundary", "xbound": "pp", "ybound": "pp", "zbound": "ff",
                          'comment': 'X,Y periodic boundaries. Z fixed'},
                         {"name": "newton", "toggle": "off",
                          'comment': 'Compute forces between local and ghost atoms in each processor w/out communication'},
                         {"name": "processors", "x": "*", "y": "*", "z": "1", 'comment': 'Processor grid'},
                         {"name": "comm_modify", "param": "vel", "toggle": "yes",
                          'comment': 'Communicate velocities for ghost atoms'},
                         {"name": "read_data", "filename": "atom.in",
                          'comment': 'Read atom.in file defining domain size and initial atoms'},
                         {"name": "region", "num": "1", 'type': 'block', 'params': 'INF INF INF INF INF INF',
                          'comment': 'Read atom.in file defining domain size and initial atoms'},
                         ]
             }
        ],

        "microbes_and_groups": [
            {"title": "#----Microbes and functional groups----#",
             'bug_groups': [
                 {"name": "group", "group_name": "het_1", 'param': "type", 'group_num': '1',
                  'comment': '# Regular heterotrophs'},
                 {"name": "group", "group_name": "het_slow", 'param': "type", 'group_num': '2',
                  'comment': "# Slow growers heterotrophs"},
                 {"name": "group", "group_name": "dead", 'param': "empty"}
             ],

             'neighbors': [
                 {'name': '#Controlling nearest neighbor recalculation'},
                 {"name": "neighbor", "distance": "2e-6", "method": "bin",
                  'comment': '# neighbor skin distance and style'},
                 {"name": "neigh_modify", "when": "every", "interval": "1", 'comment': '# Rebuild neighbor distance'}]

             }
        ],

        "mesh_grid_and_substrates": [
            {"title": "#----Mesh Grid and Substrates----#",
             "content": [{"name": "# define grid style, substrate names, and grid size"},
                         {"name": "grid_style", "loc": "nufeb/chemostat", 'p1': "4", 's1': 'sub', 's2': 'o2',
                          's3': 'no2', 's4': 'no3', 'grid_cell': '10e-6',
                          'comment': ''},
                         {"name": '# set diffusion boundary conditions and initial concentrations (liquid:kg/m3'},
                         {"name": "grid_modify", "action": "set", 'substrate': "sub", 'bound1': 'pp', 'bound2': 'pp',
                          'bound3': 'nd', 'init_conc': '1e-3',
                          'comment': '1 mg/L'},
                         {"name": "grid_modify", "action": "set", 'substrate': "o2", 'bound1': 'pp', 'bound2': 'pp',
                          'bound3': 'nd', 'init_conc': '1e-4',
                          'comment': '0.1 mg/L'},
                         {"name": "grid_modify", "action": "set", 'substrate': "no2", 'bound1': 'pp', 'bound2': 'pp',
                          'bound3': 'nd', 'init_conc': '1e-4',
                          'comment': '0.1 mg/L'},
                         {"name": "grid_modify", "action": "set", 'substrate': "no3", 'bound1': 'pp', 'bound2': 'pp',
                          'bound3': 'nd', 'init_conc': '1e-4',
                          'comment': '0.1 mg/L'}
                         ]
             }
        ],

        "biological_processes": [
            {"title": "#----Biological Processes----#",
             "growth": [{'name': 'Heterotrophic growth'},
                         {"name": "fix", "fix_name": "growth_het_1", 'group': 'het_1', 'fix_loc': "nufeb/growth/het",
                          'l1': 'sub', 'v1': '1e-3', 'l2': 'o2', 'v2': '0', 'l3': 'no2', 'v3': '0',
                          'l4': 'no3', 'v4': '0', 'l5': 'growth', 'v5': '0.000278',
                          'l6': 'yield', 'v6': '0.61', 'l7': 'decay', 'v7': '4e-6',
                          'comment': "#Regular heterotrophs"},
                         {"name": "fix", "fix_name": "growth_het_slow", 'group': 'het_slow',
                          'fix_loc': "nufeb/growth/het",
                          'l1': 'sub', 'v1': '1e-3', 'l2': 'o2', 'v2': '0', 'l3': 'no2', 'v3': '0',
                          'l4': 'no3', 'v4': '0', 'l5': 'growth', 'v5': '0.000144',
                          'l6': 'yield', 'v6': '0.61', 'l7': 'decay', 'v7': '4e-6',
                          'comment': "#Slow heterotrophs"}
                    ],
             'division':[
                         {'name': 'Heterotroph division'},
                         {"name": "fix", "fix_name": "diva", 'group': 'het_1', 'fix_loc': "nufeb/division/coccus",
                          'div_diam': "1.1e-6", 'seed': '1234', 'comment': ""},
                         {"name": "fix", "fix_name": "divv", 'group': 'het_slow', 'fix_loc': "nufeb/division/coccus",
                          'div_diam': "1.3e-6", 'seed': '1234', 'comment': ""},
                    ],
             'lysis':[],
             'death':[
                         {'name': 'Remove dead/dying, indicated by size'},
                         {"name": "fix", "fix_name": "death", 'fix_group': 'all', "fix_loc": "nufeb/death/diameter",
                          'to-group': 'dead', 'diameter': '2.5e-7',
                          'comment': '#Remove all cells from dead group if diameter below 2.5e-7'},
                         {"name": "fix", "fix_name": "rem_dead", "fix_group": "dead", 'fix_loc': 'evaporate', 'p1': '1',
                          'p2': '1000000',
                          'p3': '1', 'seed': '1701', 'comment': ''}
                         ],
             't6ss':[],
             }
        ],

        "physical_processes": [
            {"title": "#----Physical Processes----#",
             "content": [{'name': 'Pairwise interaction between atoms'},
                         {"name": "pair_style", "pair_loc": "gran/hooke/history", 'p1': '1e-4', 'p2': "NULL",
                          'p3': '1e-5', 'p4': 'NULL', 'p5': '0.0', 'p6': '0',
                          'comment': ""},
                         {"name": "pair_coeff", "p1": "*", 'p2': '*', 'comment': ""},
                         {'name': 'Pairwise interaction between z-wall and atoms'},
                         {"name": "fix", "fix_name": "wall", 'group': 'all', 'fix_loc': "wall/gran",
                          'loc2': 'hooke/history',
                          'p1': "1e-3", 'p2': 'NULL', 'p3': '1e-4', 'p5': 'NULL 0 0', 'plane': 'zplane',
                          'p6': '0.0', 'p7': '1e-04', 'comment': ''},
                         {"name": "fix", "fix_name": "vis", 'group': 'all', 'fix_loc': "viscous",
                          'p1': "1e-5", 'comment': 'Viscous damping force'},
                         {"name": "fix", "fix_name": "nve", 'fix_group': 'all', "fix_loc": "nve/limit", 'p1': '1e-7',
                          'comment': 'NVE integration w/max dist. limit'}
                         ]
             }
        ],
        "post_physical_processes": [
            {"title": "#----Post-Physical Processes----#",
             "content": [{'name': 'define diffusion coeff in biofilm'},
                         {"name": "fix", "fix_name": "coeff_sub", 'group': 'all', 'fix_loc': "nufeb/diffusion_coeff",
                          'sub1': 'sub', 'p1': 'ratio', 'coeff1': '0.8', 'comment': 'within biofilm, 80% of bulk coeff'}
                         ]
             }
        ],
        "chemical_processes": [
            {"title": "#----Chemical Processes----#",
             "content": [{'name': 'diffusion reaction for updating substrate concentration distributions'},
                         {"name": "fix", "fix_name": "diff_sub", 'group': 'all', 'fix_loc': "nufeb/diffusion_reaction",
                          'sub1': 'sub', 'coeff1': '1.6e-9', 'comment': ''}
                         ]
             }
        ],
        "computation_output": [
            {"title": "#----Computations and Outputs----#",
             "ab_track": [{'name': 'Abundance Tracking'},
                         {"name": "variable", "varname": "nhet1", 'op': 'equal', 'expression': '"count(het_1)"',
                          'comment': 'total number of regular heterotrophs'},
                         {"name": "variable", "varname": "nhetslow", 'op': 'equal', 'expression': '"count(het_slow)"',
                          'comment': 'total number of slow growing heterotrophs'},
                         ],
             "thermo_output":[{'name': 'screen outputs'},
                         {'name': 'thermo_style', 'p1': 'custom',
                          'dumpvars': 'step cpu atoms v_nhet1 v_nhetslow v_mass', 'comment': ''},
                         {'name': 'thermo', 'p1': '1', 'comment': ''},
                        ],
             "hdf5_output":[],
             "vtk_output":[],
             "content": [{'name': 'Variables used for later output'},
                         {"name": "variable", "varname": "mass", 'op': 'equal', 'expression': '"mass(all)"',
                          'comment': '# total biomass'},
                         ]
             }
        ],
        "run": [
            {"title": "#----Run----#",
             "content": [{'name': 'run with timestep for physical (pairdt) and chemical (diffdt) processes'},
                         {"name": "run_style", "run_type": "nufeb", 'p1': 'diffdt', 'v1': '1e-4',
                          'p2': 'difftol', 'v2': '1e-6', 'p3': 'pairdt', 'v3': '1e-2',
                          'p4': 'pairtol', 'v4': '1', 'p5': 'screen', 'v5': 'no',
                          'p6': 'pairmax', 'v6': '1000', 'p7': 'diffmax', 'v7': '5000',
                          'comment': ''},
                         {"name": "timestep", "val": "900",
                          'comment': '# define biological timestep (900s = 15 min)'},
                         {"name": "run", "val": "86400",
                          'comment': '# run duration (24h)'},
                         {"name": "shell", "val": "touch done.tkn",
                          'comment': '# indicate run completed successfully'}
                         ]
             }
        ]
    }
    TEMPLATE_STR = \
"""{% for line in header -%}
{{ line }}
{% endfor %}

{% for section in system_settings -%}
{{ section.title }}
{% for entry in section.content -%}
{%- set values = entry.values() | list -%}
{% if values | length < 2 %}{{ "\n" }}# {{ values[0] }} {% else %}{{ values[0] }}\t\t{{ values[1:-1] | join(' ') }}{% if values[-1] != '' %}\t\t#{{ values[-1] }}{% endif -%}{% endif %}
{% endfor %}
{% endfor %}


{%- for section in microbes_and_groups -%}
{{ section.title }}
    {% for entry in section.bug_groups %}
        {%- set values = entry.values() | list -%}
{{ values | join(' ') }}
{% endfor %}
    {%- for entry in section.neighbors %}
        {%- set values = entry.values() | list -%}
{{ values | join(' ') }}
{% endfor %}
{% endfor %}

{% for section in mesh_grid_and_substrates -%}
{{ section.title }}
    {%- for entry in section.content -%}
        {%- set values = entry.values() | list -%}
{{ values | join(' ') }}
{% endfor -%}
{% endfor -%}

{%- for section in biological_processes %}
  {%- for subsection, items in section.items() %}
    {%- if subsection == 'title' %}
{{ items }}
    {%- else %}
      {%- for entry in items %}
        {%- set values = entry.values() | list %}
        {%- if values | length < 2 %}

# {{ values[0] }}
        {%- else %}
{{ values | join(' ') }}
        {%- endif %}
      {%- endfor %}
    {%- endif %}
  {%- endfor %}
{%- endfor %}



{% for section in physical_processes -%}
{{ section.title }}
{% for entry in section.content -%}
{%- set values = entry.values() | list -%}
{%- if values | length < 2 -%}
{{ "\n" }}# {{ values[0] }} 
{%- else -%}
{{ values[0] }}\t\t{{ values[1:-1] | join(' ') }}{% if values[-1] != '' %}\t\t#{{ values[-1] }}{% endif -%}
{% endif %}
{% endfor %}
{% endfor %}

{% for section in post_physical_processes -%}
{{ section.title }}
{% for entry in section.content -%}
{%- set values = entry.values() | list -%}
{%- if values | length < 2 -%}
{{ "\n" }}# {{ values[0] }} 
{%- else -%}
{{ values[0] }}\t\t{{ values[1:-1] | join(' ') }}{% if values[-1] != '' %}\t\t#{{ values[-1] }}{% endif -%}
{% endif %}
{% endfor %}
{% endfor %}

{% for section in chemical_processes -%}
{{ section.title }}
{% for entry in section.content -%}
{%- set values = entry.values() | list -%}
{%- if values | length < 2 -%}
{{ "\n" }}# {{ values[0] }} 
{%- else -%}
{{ values[0] }}\t\t{{ values[1:-1] | join(' ') }}{% if values[-1] != '' %}\t\t#{{ values[-1] }}{% endif -%}
{% endif %}
{% endfor %}
{% endfor %}

{%- for section in computation_output %}
  {%- for subsection, items in section.items() %}
    {%- if subsection == 'title' %}
{{ items }}
    {%- else %}
      {%- for entry in items %}
        {%- set values = entry.values() | list %}
        {%- if values | length < 2 %}

# {{ values[0] }}
        {%- else %}
{{ values | join(' ') }}
        {%- endif %}
      {%- endfor %}
    {%- endif %}
  {%- endfor %}
{%- endfor %}

{% for section in run -%}
{{ section.title }}
{% for entry in section.content -%}
{%- set values = entry.values() | list -%}
{%- if values | length < 2 -%}
{{ "\n" }}# {{ values[0] }} 
{%- else -%}
{{ values[0] }}\t\t{{ values[1:-1] | join(' ') }}{% if values[-1] != '' %}\t\t#{{ values[-1] }}{% endif -%}
{% endif %}
{% endfor %}
{% endfor %}

"""

    def __init__(self):
        self.config_vals = self.DEFAULT_INPUTSCRIPT
        self.group_assignments = {}

    def clear_bug_groups(self, keep_dead=True):
        if(keep_dead):
            self.config_vals['microbes_and_groups'][0]['bug_groups'] = [{"name": "group", "group_name": "dead", 'param': "empty", 'comment': "# Dead cells"}]
        else:
            self.config_vals['microbes_and_groups'][0]['bug_groups'] = []

    def build_bug_groups(self, active_taxa, lysis_groups, keep_dead=True):
        self.clear_bug_groups(keep_dead)
        all_groups = {**active_taxa, **lysis_groups}
        for i, k in enumerate(all_groups):
            entry = {'name': 'group', 'group_name': k, 'param': 'type', 'group_num': i+1}
            self.group_assignments[k]=i+1
            if 'description' in all_groups[k]:
                entry['comment'] = f'# {all_groups[k]["description"]}'
            self.config_vals['microbes_and_groups'][0]['bug_groups'].append(entry)

    def build_substrate_grid(self, substrates,simbox):
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
                    grid_style_dict['grid_cell'] = f'{self._pick_grid_size(simbox)}'
                    new_contents.append(grid_style_dict)

        for substrate in substrates:
            #print(substrate)
            new_contents.append(substrates[substrate].as_grid_modify_dict())
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
        possible_sizes = [2.5,2.0,1.5]
        for size in possible_sizes:
            if s.xlen%size==s.ylen%size==s.zlen%size==0:
                return size*1e-6

        # we'd prefer slightly larger, but return a grid size of 1 if needed
        size = 1.0
        if s.xlen % size == s.ylen % size == s.zlen % size == 0:
            return size * 1e-6
        raise ValueError(f'No valid grid size between 5 and 15 microns for a simulation of dimensions {s.dim_string()}')


    def clear_growth_strategy(self):
        self.config_vals['biological_processes'][0]['growth'] = []


    def build_growth_strategy(self, active_taxa):
        self.clear_growth_strategy()
        self.config_vals['biological_processes'][0]['growth'].append({'name':'Growth Strategies'})
        for k,v in self.group_assignments.items():
            if k not in active_taxa:
                break
            if 'growth_strategy' not in active_taxa[k]:
                raise KeyError(f"Growth strategy is not defined for taxon: {k}")
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
            else:
                raise KeyError(f"Taxon {k} has unrecognized growth strategy: {active_taxa[k]['growth_strategy']['name'] }")


    def clear_division(self):
        self.config_vals['biological_processes'][0]['division'] = []

    def build_division(self, active_taxa, seed):
        self.clear_division()
        self.config_vals['biological_processes'][0]['division'].append({'name':'Division'})

        for k,v in self.group_assignments.items():
            if k not in active_taxa:
                break
            if 'division_strategy' not in active_taxa[k]:
                raise KeyError(f"Division strategy is not defined for taxon: {k}")
            if active_taxa[k]['division_strategy']['name'] == 'divide_coccus':
                entry ={}
                entry["name"] = 'fix'
                entry["ID"] = f'division_{k}'
                entry['group-ID'] = k
                entry['fix_loc'] = "nufeb/division/coccus"
                for kj,vj in active_taxa[k]['division_strategy'].items():
                    if kj != 'name':
                        if kj != 'eps_dens':
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

    def add_hdf5_output(self):
        self.config_vals['computation_output'][0]['hdf5_output'] = [
            {'name': 'HDF5 output, efficient binary format for storing many atom properties'},
            {'name': 'requires NUFEB built with HDF5 option'},
            {'name': 'shell', 'command': 'mkdir hdf5', 'comment': '#Create directory for dump'},
            {'name': 'dump', 'dumpname': 'du3', 'group': 'all', 'format': 'nufeb/hdf5', 'p1': '1',
                'loc': 'hdf5/dump.h5',
                'dumpvars': 'id type x y z radius', 'comment': ''},
         ]

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

    def build_run(self, runtime):
        content = self.config_vals['run'][0]['content']
        bio_timestep=900
        for i,item in enumerate(content):
            if item['name']=='run':
                run_hours = runtime/60/60/bio_timestep
                new_item = {"name": "run", "val": f'{runtime}', 'comment': f'# run duration ({run_hours} H)'}
                content[i]=new_item

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
        # Create a Template object
        template = Template(self.TEMPLATE_STR)
        config_content = template.render(self.config_vals)
        return config_content