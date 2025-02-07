# Changelog

# version 0.0.2

## New functions

* added ``poisson_disc`` layout option
  * corresponding example added

* can layout taxa in even strips horizontally or vertically with ``distribute_even_strips``
  * strips can be made noisy on a basis of 0-100 percent shuffling
  * examples added


## Enhancements

* Raise error if an invalid taxa assignment strategy is specified

## Code internals

* removed unused ``NufebProject.layout`` variable
* ``BugPos`` dataclass moved to own file
* Began trying use of Literal, introducing exceptions

# version 0.0.1 

Initial public release on 31-Jan-2025.