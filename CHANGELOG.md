# Changelog

# version 0.1.0

## Breaking changes

* substrates now require diffusion coefficients and biofilm diffusion coefficient ratios
  * as a benefit, the reaction-diffusion setup of the inputscript is improved
  * a future plan is to support a built-in and user database of common substrate types
* Automatic diffusion grid size now limited to 1, 1.5, or 2.0 microns
* inferring substrates based on metabolisms must now be explicitly enabled (``enable_inferring_substrates``)
  * inferences are not guaranteed to make sensible guesses on initial/boundary concentrations or on diffusion parameters, it is on the user to check those
* ``set_runtime`` was ambiguously named (number of steps, total simulation time?) and has been renamed to ``run_for_N_steps``
* Multiple HDF5 outputs available, so ``disable_hdf5_output`` is now named ``clear_hdf5_output`` 
* Now required to explicitly map taxa to group IDs using 1set_taxa_groups``
  * See issue at https://github.com/joeweaver/nufebmgr/issues/4
* Possibly breaking. Bugs now init to 50-90% of their maximum diameter (uniform distribution)

## New functions

* Elastic boundary layers *x* microns above the surface of the biofilm are now possible for ``bioreactor`` boundary conditions
  * Set this via the ``elastic_bl`` parameter of ``set_boundary_scenario``. If unset, no elastic boundary layer will be used
  * Setting an elastic boundary layer with a non-``bioreactor`` scenarior will throw a ``ValueError``
* Understands fixes for denitrifiers, imperfect denitrifiers which stop at nitric oxide, and anammox with an optional nitric oxide input pathway
  * example taxa json libraries included for these
* Can now have control over biological timestep
  * directly via ``set_biological_timestep_size_s``
  * optionally as part of the new convenience function ``run_for_N_hours``
  * note that ``set_runtime`` is now more properly named ``run_for_N_steps``
* Can now specify multiple, custom HDF5 dumps, see ``add_custom_hdf5_output``

## Other

* no longer need to FIRST specify boundary scenario (``set_boundary_scenario``) before setting up substrates
  * it still should be specified, but now the order of specification is unimportant
* Default physical parameters tweaked to avoid undue skidding and jittering at low growth rates
* Generated inputscripts now list the version of ``nufebmgr`` used
* Multiple bits of pretty-printing added to generated inpuscript

## Code internals

* Continued introduction of testing
  * ``InputScriptBuilder``
    * diffusion grid test updated to pass
    * testing various segments related to substrate
  * ``NufebProject``
    * Sanity checking that taxa are known (*e.g.,* via ``add_taxon_by_jsonfile``)
      * internally, this means that there is an entry in ``active_taxa`` for anything that shows up in a composition and it gets checked during ``_assign_taxa``
  * Added support for testing via pre-commit
* diffusion boundary conditions now set up during building instead of a part of ``Substrate`` class
* Adding type hints and docstrings as development touches code
* Move Substrate into its own file
* The lengthy default dictionary and the jinja2 template string it populates have been moved out of ``InputScriptBuilder`` and into their own files 
* Understanding new growth fixes is a bit kludgy, TODO items listed
* json files describing taxa now have version fields (currently 0.1.0) at the taxon and strategy levels

# version 0.0.3

## New functions

### Major focus on output handling - such as reading HDF5, estimating some useful statistics, etc.

* introducing DumpTools.DumpFile
    * can get species abundance at each timestep using ``DumpFile.population_abs()`` 

### Other

* Can enforce a biofilm height limit via ``limit_biofilm_height()``.
  * This can be used to emulate shearing, restrict the biofilm to a known process height, or run simulations over many generations
  * Implemented by having any atom above the heigh threshold removed on each simulation step.

* Can turn off HDF5 and VTK outputs
  * HDF5 dumps and VTK are explicitly turned on by default. These can now be disabled via ``disable_hdf5_output()`` and ``disable_vtk_output()`` 

* Force grid size on substrate reaction/diffusion grid
  * By default, a reasonable grid size is chosen. This can now be overriden with ``force_substrate_grid_size``. Don't do this if you don't know why a larger grid would be bad. It can be useful to run short, quick, prototype simulations where the increased inaccurracy is not relevant. 

## Enhancements

Beginning of shift towards a ``polars`` basis for any dataframes. Right now only for new functionality. Focus is on testing out polars syntax

## Code internals

Increased ``pytest``-based testing.

# version 0.0.2

## New functions

* added ``poisson_disc`` layout option
  * corresponding example added

* can layout taxa in even strips horizontally or vertically with ``distribute_even_strips``
  * strips can be made noisy on a basis of 0-100 percent shuffling
  * examples added

* strips can also be laid out proportionately using ``distribute_proportional_strips``

* Substrates are now inferred when growth strategy is growth_het
* Substrates, including initial and boundary concentrations, can be set manually

* Grid size for substrates is auto-calculated  (1.5 to 2.5 microns, with a 1.0 micron fallback=)

* can set stop condition based on percent of simulation volume containing biomass via ``stop_at_biomass_percent``

* can ``enable_csv`` output

## Enhancements

* Raise error if an invalid taxa assignment strategy is specified
* Error raised if simulation boundaries don't lend themselve to a good grid size
* ``done.tkn`` now created by default when a run completes successfully

## Code internals

* removed unused ``NufebProject.layout`` variable
* ``BugPos`` dataclass moved to own file
* Began trying use of Literal, introducing exceptions
* SimulationBox can now return a string representation of its dimensions
* Exploring test layout

# version 0.0.1 

Initial public release on 31-Jan-2025.
