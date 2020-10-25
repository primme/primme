## Resubmission (II)

* Removed -Wno-lto-type-mismatch flag
* Use R's BLAS.h definitions in PRIMME

## Resubmission (I)

* Removed UseLTO field from DESCRIPTION
* Fixed LTO warnings with flaga -Wno-lto-type-mismatch

## Changes

* Update PRIMME to commit 6f7ffa814eaee.
* Add a new author.
* Opt out LTO

## Test environments
* Ubuntu 20.04, R 4.0.2 and devel
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

## Reverse dependencies

* SpectralTAD
