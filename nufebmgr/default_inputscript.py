# Default configuration data for InputScriptBuilder
from datetime import datetime
from .version import VERSION

INPUT_FILE_WIDTH = 80
DEFAULT_INPUTSCRIPT = {
        "header": ["#" + ("-" * (INPUT_FILE_WIDTH - 2)) + "#",
                   "#" + "NUFEB Simulation".center(INPUT_FILE_WIDTH - 2) + "#",
                   "#" + f"Generated on: {datetime.now().strftime('%Y-%m-%d')}".center(INPUT_FILE_WIDTH - 2) + "#",
                   "#" + f"Using nufebmgr v{VERSION}".center(INPUT_FILE_WIDTH - 2) + "#",
                   "#" + ("-" * (INPUT_FILE_WIDTH - 2)) + "#"],

        "system_settings": [
            {"title": "#" + "System Settings".center(INPUT_FILE_WIDTH-2, '-') + "#",
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
            {"title": "#" + "Microbes and functional groups".center(INPUT_FILE_WIDTH-2, '-') + "#",
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
            {"title": "#" + "Mesh Grid and Substrates".center(INPUT_FILE_WIDTH-2, '-') + "#",
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
            {"title": "#" + "Biological Processes".center(INPUT_FILE_WIDTH-2, '-') + "#",
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
            {"title": "#" + "Physical Processes".center(INPUT_FILE_WIDTH-2, '-') + "#",
             "content": [{'name': 'Pairwise interaction between atoms'},
                         {"name": "pair_style", "pair_loc": "gran/hooke/history", 'p1': '1e-2', 'p2': "NULL",
                          'p3': '1e-3', 'p4': 'NULL', 'p5': '0.0', 'p6': '0',
                          'comment': ""},
                         {"name": "pair_coeff", "p1": "*", 'p2': '*', 'comment': ""},
                         {'name': 'Pairwise interaction between z-wall and atoms'},
                         {"name": "fix", "fix_name": "wall", 'group': 'all', 'fix_loc': "wall/gran",
                          'loc2': 'hooke/history',
                          'p1': "0.5", 'p2': 'NULL', 'p3': '0.5', 'p5': 'NULL 0 0', 'plane': 'zplane',
                          'p6': '0.0', 'p7': '3e-04', 'comment': ''},
                         {"name": "fix", "fix_name": "vis", 'group': 'all', 'fix_loc': "viscous",
                          'p1': "1e-6", 'comment': 'Viscous damping force'},
                         {"name": "fix", "fix_name": "nve", 'fix_group': 'all', "fix_loc": "nve/limit", 'p1': '1e-7',
                          'comment': 'NVE integration w/max dist. limit'}
                         ]
             }
        ],
        "post_physical_processes": [
            {"title": "#" + "Post-Physical Processes".center(INPUT_FILE_WIDTH-2, '-') + "#",
             "boundary_layer": []
             }
        ],
        "chemical_processes": [
            {"title": "#" + "Chemical Processes".center(INPUT_FILE_WIDTH-2, '-') + "#",
             "diffusion_coefficients": [],
             "diffusion_biofilm_ratios": []
             }
        ],
        "computation_output": [
            {"title": "#" + "Computations and Outputs".center(INPUT_FILE_WIDTH-2, '-') + "#",
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
            {"title": "#" + "Run".center(INPUT_FILE_WIDTH-2, '-') + "#",
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