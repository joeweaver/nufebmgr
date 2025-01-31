# About nufebmgr
Python package for programmatically defining NUFEB cases in the language of biologists

NUFEB (Newcastle University Frontier in Engineering Biology) is an open source tool for 3D agent-based simulation of microbial systems. The tool is built on top of the molecular dynamic simulator LAMMPS, and extended with features for microbial modelling. As part of this legacy, the inputs to NUFEB are rooted in terms more familiar to molecular dynamics simulations.
https://github.com/joeweaver/nufeb-designer
nufebmgr is a python package which aims to streamline creating NUFEB inputs by 'speaking in the language of biologists'. It allows programmatic genreation of configuration files at a high level of abstraction. It is also used to enable interactive generation via the nufeb-designer GUI application.

Development was funded by the USF NSF Postdoctoral Research Fellowship in Biology Award #2007151.

Because nufebmgr is under active development and experimental, we would like you contact us if you use it. We'd enjoy hearing about how we could make it better AND we'd like to have a list of people to contact if major bugs are found. The best person to contact at this time would be Joseph E. Weaver: joe.weaver@newcastle.ac.uk

# Getting started

nufebmgr was written with Python 3.11 and has the following depdencies:

* numpy 2.2.0
* pandas 2.2.0
* jinja2 3.1.0
* opencv-python 4.9.0

A working installation of [NUFEB](https://github.com/nufeb/NUFEB-2) is not required to use nufebmgr, but generating input files for NUFEB would be rather pointless without it.

Installation is currently via setuptools with a local version. We suggest using git to clone the repository and then, from the top level directory, running ``pip install -e .``  It is also suggested that you create an isolated python environment using the tool of your choice.

# Package overview

## Use case

The motivating reasons for creating nufebmger were to create a formal programmatic way to generate NUFEB input files and to provide a way to specify scenarios at a biological level of abstraction. Many of the 'behind the scenes' and boilerplate parts of a NUFEB setup are abstracted away. You need not know about, for example, the necessary but rarely-modified ``newton`` command. You can also specify things like 'Bacteria are evenly distributed, there are two heterotrophs in the population, with fast growers having a relative abundance of 75%' rather than manually editing an ``atom.in`` file and specifying  multiple ``growth``, ``divsion``, and other fixes.  As part of this flexilbity, nufebmgr has been used to power a graphical user interface for generating NUFEB inputs, the [nufeb-designer](https://github.com/joeweaver/nufeb-designer) sister project.

We expect that some users might want to use either nufebmgr or nufeb-designer as first steps to create working examples of moderately complex scenarios and then use those examples to first understand the NUFEB input file format and then add further custom complexity.

## General operation

nufebmgr works on the basis of modifying a 'working, but boring' base case which is populated with reasonable defaults. This case is encapsulated within a ``NufebProject`` class which acts as a context manager (think of the ``with open("file.txt") as f:`` construct).

Within that context various calls are used to alter the base case, such as to set the simulation box dimensions or specify how bacteria should be physically arranged.  There is also a ``generate_case`` method which returns text suitable for saving as an ``atom.in`` and ``inputscript.nufueb`` (the two files defining a NUFEB case).

We have tried to make the order of calls within the context not matter, and we also try to generate valid case files whenever possible (*i.e.* not requiring any specific method call within the context). However, this is not always possible and given that the package was driven by the specific cases we wanted to generate for our research, there are possibly blind spots. If you run across a problem, we'd appreicate it if you [file an issue](https://github.com/joeweaver/nufebmgr/issues).

# Examples

There are four examples located in the ``examples`` directory. Each contains at least a ``main.py`` script as well as potential supporting files. These examples are intended to be plain showcases of the main nufebmgr functions.

## Hello, world: simple_starter

``simple_starter`` is intended as the equivalent of a `Hello world` example.  It shows you how to:

* specify the random number generator seed
* set the dimensions of the simulation box
* add individual taxa using some built-int templates
* specify an initial uniform random layout of *n* bacteria
* set and distribution the population composition among those individuals
* enable useful tabular output of population structure at each timestep
* generate the actual input files

## Defining your own taxa: taxa_libraries

While the built-in templates are nice for getting started, many users will want to define custom taxa. To do this with nufebmgr, users can create/edit a ``json`` formatted text file which defines the relevant properties of one or more taxa. taxa_libraries shows how to import the json file and includes the file ``example.json`` which illustrates the taxa defintions.

The relevant method call is ``add_taxon_by_jsonfile(filename)``

The json file is very similar to a python dictionary of dictionaries.  Each taxon entry begins with a name which contains a ``growth_strategy`` and ``division_stragey`` compound entry. It must also contain ``diameter``, ``density``, ``outer_diameter``, and ``morphology`` entries. It may optionally contain a ``description`` entry.  The ``growth_strategy`` and ``division_strategy`` entries must contain a ``name`` which maps to the relevant NUFEB biological fix and additional entries used to fill in the fix-specific parameters. 

An easy way to create custom taxa libraries is to use the ``Taxa Manager`` UI within nufeb-designer.

Here is an example of a library containing a single taxon, named 'foo' with a heterotrophic growth strategy and coccus division strategy.

```json
{     
  "foo": {     
    "growth_strategy": {     
      "name": "growth_het",    
      "sub-ID":"sub",    
      "sub-Ks" : 3.5e-5,    
      "o2-ID": "o2",    
      "o2-Ks": 0,    
      "no2-ID":"no2",    
      "no2-Ks" : 0,    
      "no3-ID": "no3",    
      "no3-Ks": 0,    
      "mu_max-ID": "growth",    
      "mu_max": 0.00028,    
      "yield-ID": "yield",    
      "yield": 0.61,    
      "decay-ID": "decay",    
      "decay": 0,    
      "anoxic-ID": "anoxic",    
      "anoxic": 0     
    },    
    "division_strategy": {     
      "name":"divide_coccus",    
      "diameter":1.36e-6     
    },    
    "diameter": 1e-6,    
    "density": 150,    
    "outer_diameter": 1e-6,    
    "morphology": "coccus"    
  }
}
```
## Drawing initial layouts: layout_with_image

Placing initial bacteria positions has long been an issue with NUFEB ease of use. While nufebmgr has some basic random layout functionality, it is still programmatic and hard to visualize. To aid visualization, increase flexiblity, and allow interaction with ``nufeb-designer``, nufebmgr can read images to specify the intial posititions and identities of bacteria used in the ``atom.in`` file.

Each pixel in the image represents 1 micron (so a 100x50 micron simulation plane would be represented by a 100 pixel wide and 50 pixel tall image) and the colour of the pixels represents taxa identies. Fully white pixels indicate no taxa present at that location.

<p align="center">
  <img src="https://github.com/joeweaver/nufebmgr/blob/main/doc/images/zoomed_simple.png" | width=500/>
</p>
The colours are mapped to specific identies by a user-provide python dictionary where each colour code is associated with a taxon name (example below). The taxon names are assumed to be provided either by ingesting a taxa library file or by using the built-in templates.

```python
          mappings = {"FF1B9E77": 'basic_heterotroph',    
                      "FFD95F02": 'slow_heterotroph',    
                      "FF7570B3": 'small_heterotroph'}
```

## Adding a Type VI Secretion System interaction: T6SS

The interactions and behaviours required to model the T6SS in NUFEB require populating two fixes which, while powerful, are also very low level. nufebmgr allows specifying the interactions at a higher level in Python and this gets 'compiled' into the relevant fix calls in the input script.

There are three essential steps:

1. specify the name(s) to represent intoxicated taxa and how each one behaves when intoxicated via one or more ``add_lysis_group_by_json`` calls
2. specify the name(s) of any taxa which will have a T6SS system and the associated paramters defining that system with one or more ``arm_t6ss`` calls
3. specify the name(s) of any taxa which will be vulnerable to a T6SS attack one or more ``vuln_t6SS`` calls

Each one of those calls has named paramters which are filled out in the example code and which are fully defined in the NUFEB documentation for the [T6SS](https://nufeb.readthedocs.io/en/master/fix_T6SS_contact.html) and [lysis](https://nufeb.readthedocs.io/en/master/fix_T6SS_lysis.html) fixes.

The calling code would look like this:

```python
with NufebProject() as prj:
    # snip setup not relevant to T6SS

    # define lysis behaviour
    prj.add_lysis_group_by_json('vuln_intoxicated',
        {'name':'vuln_intoxicated', 'releases':'sub', 'rate':'2e-3', 'percent':'0.2'})

    # define T6SS attacker
    prj.arm_t6ss(taxon="attacker", effector="toxin_a", harpoon_len=1.3e-6, cooldown=100)

    # define T6SS vulnerable taxon
    prj.vuln_t6ss(taxon="immune", effector="toxin_a", prob=1, to_group="vuln_intoxicated")
```

# Limitations

To constrain context nufebmgr was largely developed by 'dogfooding' based on the needs of the cases specific to the developer's research. As such, it is particularly good at setting up Type VI Secretion System simulations.  nufebmgr is currently suitable for generating simulations involving cocci but not the more recent bacillus morphology, and has been generally tested with classic Monod-style heterotrophic growth. Other growth strategies can be represented in the taxa json files, but these have not been as exhaustively tested.

nufebmgr also currently assumes sensible parameter values (*e.g.* don't request an 'apple' x-dimension of the simulation box).
