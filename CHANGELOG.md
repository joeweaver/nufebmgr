# Changelog

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
