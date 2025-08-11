import numpy as np
import pandas as pd
from jinja2 import Template
import cv2
import csv
import json
from typing import Literal
from .SimulationBox import SimulationBox
from .Substrate import Substrate
from .InputScriptBuilder import InputScriptBuilder
from datetime import datetime
from .poisson import PoissonDisc
from .TaxaAssigmentManager import TaxaAssignmentManager
from .BugPos import BugPos#

class NufebProject:
    taxa_templates = {"basic_heterotroph": {'growth_strategy':
                                                {'name': 'growth_het',
                                                 'sub-ID':'sub',
                                                 'sub-Ks' : 3.5e-5,
                                                 'o2-ID': 'o2',
                                                 'o2-Ks': 0,
                                                 'no2-ID':'no2',
                                                 'no2-Ks' : 0,
                                                 'no3-ID': 'no3',
                                                 'no3-Ks': 0,
                                                 'mu_max-ID': 'growth',
                                                 'mu_max': 0.00028,
                                                 'yield-ID': 'yield',
                                                 'yield': 0.61,
                                                 'decay-ID': 'decay',
                                                 'decay': 0,
                                                 'anoxic-ID': 'anoxic',
                                                 'anoxic': 0,
                                                 'maintain-ID': 'maintain',
                                                 'maintain': 0,
                                                 'epsyield-ID': 'epsyield',
                                                 'epsyield': 0.18,
                                                 'epsdens-ID': 'epsdens',
                                                 'epsdens': 30,
                                                 },
                                            'division_strategy':
                                                {'name':'divide_coccus',
                                                 'diameter':'1.36e-6',
                                                },
                                            'diameter': 1e-6,
                                                'density': 150,
                                                'outer_diameter': 1e-6,
                                                'morphology': 'coccus',

                                            },
                      "slow_heterotroph": {'growth_strategy':
                                               {'name': 'growth_het',
                                                'sub-ID': 'sub',
                                                'sub-Ks': 3.5e-5,
                                                'o2-ID': 'o2',
                                                'o2-Ks': 0,
                                                'no2-ID': 'no2',
                                                'no2-Ks': 0,
                                                'no3-ID': 'no3',
                                                'no3-Ks': 0,
                                                'mu_max-ID': 'growth',
                                                'mu_max': 0.00014,
                                                'yield-ID': 'yield',
                                                'yield': 0.61,
                                                'decay-ID': 'decay',
                                                'decay': 0,
                                                'anoxic-ID': 'anoxic',
                                                'anoxic': 0,
                                                'maintain-ID': 'maintain',
                                                'maintain': 0,
                                                'epsyield-ID': 'epsyield',
                                                'epsyield': 0.18,
                                                'epsdens-ID': 'epsdens',
                                                'epsdens': 30,
                                                },
                                           'division_strategy':
                                               {'name': 'divide_coccus',
                                                'diameter': '1.36e-6',
                                                },
                                           'diameter': 1e-6,
                                                'density': 150,
                                                'outer_diameter': 1e-6,
                                                'morphology': 'coccus',
                                                'description': '# Slow growing heterotroph'

                                        },
                      "small_heterotroph": {'growth_strategy':
                                                {'name': 'growth_het',
                                                 'sub-ID': 'sub',
                                                 'sub-Ks': 3.5e-5,
                                                 'o2-ID': 'o2',
                                                 'o2-Ks': 0,
                                                 'no2-ID': 'no2',
                                                 'no2-Ks': 0,
                                                 'no3-ID': 'no3',
                                                 'no3-Ks': 0,
                                                 'mu_max-ID': 'growth',
                                                 'mu_max': 0.00028,
                                                 'yield-ID': 'yield',
                                                 'yield': 0.61,
                                                 'decay-ID': 'decay',
                                                 'decay': 0,
                                                 'anoxic-ID': 'anoxic',
                                                 'anoxic': 0,
                                                 'maintain-ID': 'maintain',
                                                 'maintain': 0,
                                                 'epsyield-ID': 'epsyield',
                                                 'epsyield': 0.18,
                                                 'epsdens-ID': 'epsdens',
                                                 'epsdens': 30,
                                                 },
                                            'division_strategy':
                                                {'name': 'divide_coccus',
                                                 'diameter': '1.0e-6',
                                                 },
                                            'diameter': 0.8e-6,
                                            'density': 150,
                                            'outer_diameter': 0.8e-6,
                                            'morphology': 'coccus',
                                            },

                      }

    # Default configuration data
    DEFAULT_INPUTSCRIPT = {
        "header": ["#----------------------------------------------------------------------#"
                   "#                    NUFEB Simulation                                  #",
                   f"#               Generated on: {datetime.now().strftime('%Y-%m-%d')}   #",
                    "#----------------------------------------------------------------------#"],

        "system_settings": [
            {"title": "#----System Settings----#",
             "content": {"boundary": {"l": "pp", "b":"pp", "upper":"ff"},
                        "timeout": {"time": "30", "unit": "s"},
                         "cache_size": {"time": "30", "unit": "s"}
                         }
             }
        ]
    }

    def __init__(self,seed=1701):
        self.seed = seed
        np.random.seed(seed)
        self.sim_box = SimulationBox()
        dtype = [('col1', 'int32'), ('col2', 'int32'), ('col3', 'U10')]
        self.bug_locs = []
        self.active_taxa = {}
        self.lysis_groups = {}
        self.track_abs = False
        self.thermo_output = True
        self.thermo_timestep= 1
        self.csv_output = False
        self.csv_timestep = 0
        self.t6ss_attackers ={}
        self.t6ss_vulns = {}
        self.taxa_pre_assigned = False
        self.biostep = 900
        self.runtime = int(24*60*60/self.biostep) # default 24 hours assuming 900 s biological timestep
        self.spatial_distribution_params = {}
        self.substrates = {}  # TODO might be better to refactor as a list, since now every Substrate has its own name
        self.boundary_scenario = "bioreactor"
        self.stop_condition = "runtime"
        self.biomass_percent = None
        self.write_csv = False
        self.max_biofilm_height = None
        self.write_hdf5 = True
        self.write_vtk = True
        self.forced_substrate_grid_size=None
        self.infer_substrates=False


    # __enter__ and __exit__ for handling using project as context
    def __enter__(self):
        # Code to execute when entering the context
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type is not None:
            # If we see this, add exception to issues
            print(f"Exception occurred: {exc_type}, {exc_val}")

    def force_substrate_grid_size(self,size):
        self.forced_substrate_grid_size=size
    def disable_hdf5_output(self):
        self.write_hdf5 = False

    def disable_vtk_output(self):
        self.write_vtk = False

    def enable_inferring_substrates(self) -> None:
        """
        Enable substrate inference based on metabolisms.

        Allows guessing substrates required and setting  correct
        initial concentrations, boundary conditions and diffusion data
        """
        self.infer_substrates=True

    def disable_inferring_substrates(self) -> None:
        """
        Disable substrate inference based on metabolisms.

        Prevents guessing substrates required and setting correct
        initial concentrations, boundary conditions and diffusion data
        """
        self.infer_substrates=False

    def set_substrate(self, name: str, initial: float, bulk: float, diffusion_coefficient: float, biofilm_diffusion_ratio:float) -> None:
        new_sub = Substrate(name=name, init_concentration=initial, bulk_concentration=bulk,
                            diffusion_coefficient=diffusion_coefficient, biofilm_diffusion_ratio=biofilm_diffusion_ratio)
        self.substrates[name]=new_sub

    def stop_at_biomass_percent(self,percent: int):
        self.stop_condition ="percent biomass"
        self.biomass_percent = percent

    def set_boundary_scenario(self,scenario):
        self.boundary_scenario = scenario

    def _infer_substrates(self):
        if not self.infer_substrates:
            raise RuntimeError("Inferring substrates is not enabled. To enable call enable_inferring_substrates before calling _infer_substrates")

        subs_names = []
        # get a list of all substrates associated with growth strategies of taxa
        for taxon_name in self.active_taxa:
            # TODO for 0.1.X redo json for taxa so that we don't need to know specific growth strategy implementations to get substrate names
            growth_strat = self.active_taxa[taxon_name]['growth_strategy']['name']
            if growth_strat == 'growth_het':
                subs_names.append(self.active_taxa[taxon_name]['growth_strategy']['sub-ID'])
                subs_names.append(self.active_taxa[taxon_name]['growth_strategy']['o2-ID'])
                subs_names.append(self.active_taxa[taxon_name]['growth_strategy']['no2-ID'])
                subs_names.append(self.active_taxa[taxon_name]['growth_strategy']['no3-ID'])

        subs_names = set(subs_names)
        for sub_name in subs_names:
            if sub_name not in self.substrates:
                self.set_substrate(sub_name,1e-4,1e-4,2e-9,0.8)


    def use_seed(self,seed=1701):
        self.seed = seed
        np.random.seed(seed)

    def set_box(self,x=100,y=100,z=100,periodic="plane",custom={'x':'','y':'','z':''}):
        self.sim_box = SimulationBox(xlen=x, ylen=y, zlen=z, periodic=periodic, custom=custom)

    def distribute_spatially_even(self):
        self.spatial_distribution = "even"

    def distribute_even_strips(self, direction: Literal["horizontal", "vertical"], noise=0):
        allowed_dirs = {"horizontal", "vertical"}
        if direction not in allowed_dirs:
            raise ValueError(f"Invalid strip direction: {direction}. Must be one of {allowed_dirs}.")
        self.spatial_distribution = "strips"
        self.spatial_distribution_params["direction"] = direction
        self.spatial_distribution_params["strip_proportion"] = "even"
        self.spatial_distribution_params["noise"] = noise

    def distribute_proportional_strips(self, direction: Literal["horizontal", "vertical"], noise=0):
        allowed_dirs = {"horizontal", "vertical"}
        if direction not in allowed_dirs:
            raise ValueError(f"Invalid strip direction: {direction}. Must be one of {allowed_dirs}.")
        self.spatial_distribution = "strips"
        self.spatial_distribution_params["direction"] = direction
        self.spatial_distribution_params["strip_proportion"] = "proportional"
        self.spatial_distribution_params["noise"] = noise

    def layout_poisson(self, radius):
        poisson_disc = PoissonDisc(self.sim_box.xlen*1e-6, self.sim_box.ylen*1e-6, radius*1e-6)
        s = poisson_disc.sample()
        self.bug_locs = [BugPos(x,y,taxon_name="Unassigned") for x,y in s]

    def layout_uniform(self, nbugs):
        base = np.random.rand(nbugs, 2)
        bugs_xy = base * np.array([self.sim_box.xlen, self.sim_box.ylen]).tolist() * 1e-6
        bugs_xy = np.round(bugs_xy, decimals=8)
        self.bug_locs = [BugPos(*row, taxon_name="Unassigned") for row in bugs_xy]


    def add_taxon_by_template(self, name, template):
        self.active_taxa[name] = NufebProject.taxa_templates[template]

    def set_taxa(self,taxa):
        self.active_taxa=taxa

    def add_taxon_by_jsonfile(self, jsonfile):
        with open(jsonfile) as f:
            d = json.load(f)
        for item in d:
            a = 1
            self.active_taxa[item] = d[item]

    def set_composition(self, composition):
        self.composition = composition

    def _n_members(self):
        return len(self.bug_locs)

    def _n_taxa(self):
        return len(self.active_taxa)

    def _n_types(self):
        return len(self.group_assignments.keys())

    def _assign_taxa(self):
        if self.taxa_pre_assigned:
            return

        if self.spatial_distribution == "even":
            # TODO error check that all taxa have entries
            a1 =list(self.active_taxa.keys())
            s1 =self._n_members()
            all_compositions = [float(value) for value in self.composition.values()]
            total = sum(all_compositions)
            p1 = [value / total for value in all_compositions]
            assignments = np.random.choice(a1, size=s1, p= p1, replace=True)
            for bug, taxon in zip(self.bug_locs, assignments):
                bug.taxon_name = taxon
        elif self.spatial_distribution == "strips":
            tam = TaxaAssignmentManager(self.bug_locs)
            if self.spatial_distribution_params["direction"] == "horizontal":
                cutdir = "y"
                cutdim = self.sim_box.ylen*1e-6
            elif self.spatial_distribution_params["direction"] == "vertical":
                cutdir = "x"
                cutdim = self.sim_box.xlen*1e-6
            else:
                raise ValueError(f'Invalid spatial distribution direction parameter for strip layout')
            if self.spatial_distribution_params["strip_proportion"] == "even":
                self.bug_locs = tam.even_strips(list(self.active_taxa.keys()), cutdir,
                                                self.spatial_distribution_params['noise'])
            elif self.spatial_distribution_params["strip_proportion"] == "proportional":
                self.bug_locs = tam.proportional_strips(list(self.active_taxa.keys()),self.composition, cutdim, cutdir, self.spatial_distribution_params['noise'] )
            else:
                raise ValueError(f'Invalid spatial distribution "strip_proportion" parameter for strip layout')
        else:
            raise ValueError(f"Invalid spatial distribution: {self.spatial_distribution}. Must be `strips` or `even`.")

    def limit_biofilm_height(self,max_height):
        self.max_biofilm_height = max_height

    def generate_case(self):
        # because bits of these depend on each other, we enforce order of calling
        inputscript = self._generate_inputscript()
        atom_in = self._generate_atom_in()
        return atom_in, inputscript

    def _generate_atom_in(self):
         self._assign_taxa()

         df = pd.DataFrame([dc.__dict__ for dc in self.bug_locs])
         df["taxon_id"] = df["taxon_name"].map(self.group_assignments)


         df2 = df.assign(**df['taxon_name'].map(self.active_taxa).apply(pd.Series))

         # Define the values for the placeholders
         config_values = {
             'n_atoms': f'{self._n_members()}',
             'n_types': f'{self._n_types()}',
             'xlen_m': f'{self.sim_box.xlen}e-6',
             'ylen_m': f'{self.sim_box.ylen}e-6',
             'zlen_m': f'{self.sim_box.zlen}e-6',
         }

         # Add the atoms list to the config values
         config_values['atoms'] = df2.reset_index().to_dict(orient='records')


         template_str =\
 """NUFEB Simulation

\t\t{{ n_atoms }} atoms
\t\t{{ n_types }} atom types
\t\t0 {{ xlen_m }} xlo xhi
\t\t0 {{ ylen_m }} ylo yhi
\t\t0 {{ zlen_m }} zlo zhi

\tAtoms

 {% for row in atoms -%}
  {{'\t'}}{{ row.index +1}} {{ row.taxon_id }} {{ "%.2e" | format(row.diameter) }} {{row.density}}  {{ row.x }}  {{ "%.2e" | format(row.y )}} {{"%2e" | format(row.diameter)}} {{"%.2e" | format(row.outer_diameter)}}
 {% endfor %}
 """

         # Create a Template object
         template = Template(template_str)
         config_content = template.render(config_values)
         return(config_content)

    def set_runtime(self,time):
        self.runtime = time

    def _generate_inputscript(self):
        isb = InputScriptBuilder()

        if self.infer_substrates:
            self._infer_substrates()

        isb.build_substrate_grid(self.substrates, self.sim_box, self.boundary_scenario, self.forced_substrate_grid_size)
        isb.build_diffusion(self.substrates)
        isb.build_bug_groups(self.active_taxa,self.lysis_groups)
        isb.clear_growth_strategy()
        isb.build_growth_strategy(self.active_taxa)
        isb.clear_division()
        isb.build_division(self.active_taxa, self.seed)

        isb.build_lysis(self.lysis_groups)
        isb.build_t6ss(self.t6ss_attackers, self.t6ss_vulns, self.seed)

        if(self.track_abs):
            isb.add_abs_vars()

        if(self.thermo_output):
            isb.add_thermo_output(self.track_abs, self.thermo_timestep)

        if(self.write_hdf5):
            isb.add_hdf5_output()

        if (self.write_vtk):
            isb.add_vtk_output()

        if(self.write_csv):
            isb.enable_csv_output(self.track_abs, self.stop_condition=="percent biomass")

        if self.stop_condition=="percent biomass":
            isb.build_run(365*24*60*60)
            isb.track_percent_biomass(self.sim_box)
            isb.end_on_biomass(self.biomass_percent)
        elif self.stop_condition=="runtime":
            isb.build_run(self.runtime)
        else:
            raise ValueError(f'Unknown stop condition: {self.stop_condition}')

        if self.max_biofilm_height is not None:
            isb.limit_biofilm_height(self.max_biofilm_height)

        self.group_assignments = isb.group_assignments
        return isb.generate()


    # def match_color(self, bgr_color, predefined_colors):
    #     min_distance = float('inf')
    #     closest_color = None
    #     for color in predefined_colors:
    #         distance = np.linalg.norm(np.array(bgr_color) - np.array(color['bgr']))
    #         if distance < min_distance:
    #             min_distance = distance
    #             closest_color = color
    #     return closest_color

    def simple_image_layout(self, imagefile, mappings):
        self.taxa_pre_assigned = True
        image = cv2.imread(imagefile)
        grey_image = cv2.imread(imagefile, cv2.IMREAD_GRAYSCALE)

        # Count non-white (non-blank) pixels
        nbugs = np.count_nonzero(grey_image != 255)

        height,width,color_components = image.shape
        bug_num = 0
        for x in range(width):
            for y in range(height):
                if np.all(image[y][x] != [255,255,255]):
                    b,g,r = image[y][x]
                    color_code = f'FF{r:02X}{g:02X}{b:02X}'
                    self.bug_locs.append(BugPos(x*1e-6,(height-y-1)*1e-6,mappings[color_code]))
                    bug_num += 1


    # def layout_and_distribute_image(self, imagefile, mappings):
    #     predefined_colors = [
    #         {'name': 'Red', 'bgr': (0, 0, 255)},
    #         {'name': 'Green', 'bgr': (0, 255, 0)},
    #         {'name': 'Blue', 'bgr': (255, 0, 0)},
    #         {'name': 'Yellow', 'bgr': (0, 255, 255)},
    #         {'name': 'Cyan', 'bgr': (255, 255, 0)}
    #     ]
    #     image = cv2.imread(imagefile)
    #
    #     # Convert to grayscale
    #     gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    #
    #     # Threshold the image to create a binary mask
    #     _, binary = cv2.threshold(gray, 1, 255, cv2.THRESH_BINARY)
    #
    #     # Find contours of the blobs
    #     contours, _ = cv2.findContours(binary, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    #
    #     dots = []
    #     output_image = np.zeros_like(image)
    #
    #     for contour in contours:
    #         # Compute the center of the contour using moments
    #         M = cv2.moments(contour)
    #         if M["m00"] > 0:
    #             cx = int(M["m10"] / M["m00"])  # x-coordinate of the center
    #             cy = int(M["m01"] / M["m00"])  # y-coordinate of the center
    #
    #             # Get the color at the center
    #             bgr_color = image[cy, cx].tolist()
    #
    #             # Match the color to predefined colors
    #             matched_color = self.match_color(bgr_color, predefined_colors)
    #
    #             # Save the dot information
    #             dots.append({'x': cx, 'y': cy, 'color_name': matched_color['name'], 'bgr': matched_color['bgr']})
    #
    #             # Draw the matched dot on the output image
    #             cv2.circle(output_image, (cx, cy), radius=0, thickness=-1, color=tuple(matched_color['bgr']) )
    #
    #     # Save the coordinates and colors to a CSV file
    #     output_csv_path = 'dots_output.csv'
    #     with open(output_csv_path, mode='w', newline='') as file:
    #         writer = csv.DictWriter(file, fieldnames=['x', 'y', 'color_name', 'bgr'])
    #         writer.writeheader()
    #         writer.writerows(dots)
    #
    #     # Save the output image
    #     output_image_path = 'output_image.png'
    #     cv2.imwrite(output_image_path, output_image)


    def set_track_abs(self, do_track=True):
        self.track_abs = do_track

    def enable_thermo_output(self, timestep=1):
        self.thermo_output=True
        self.themo_timestep=timestep

    def enable_csv(self):
        self.write_csv = True

    def add_lysis_group_by_json(self,name,definition):
        self.lysis_groups[name]=definition

    def arm_t6ss(self, taxon, effector, harpoon_len, cooldown):
        self.t6ss_attackers[taxon] ={'effector':effector, 'harpoon_len':harpoon_len,'cooldown':cooldown}

    def vuln_t6ss(self, taxon, effector, prob, to_group):
        self.t6ss_vulns[taxon] = {'effector': effector, 'prob': prob, 'to_group': to_group}