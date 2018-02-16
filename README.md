# MedFlyCRISPR-CAS9

This repository contains the data for the MedFlies population suppression model described in the paper: *Consequences of instant induction of resistance evolution on a sex conversion-based suppression gene drive for insect pest management*.

## Datasets

The datasets are compressed in _7z_ format, which can be unarchived using free software such as: [Keka](http://www.kekaosx.com/en/) or [7zip](http://www.7-zip.org/download.html).

### Fertile_WDXX

Zip files naming convention

* N001: Population size of 1000
* N010: Population size of 10000
* N100: Population size of 100000
* Visual: Long term population dynamics for visual inspection

### Infertile_WDXX

* OneReleaseSixValues: One release at various values of rhoR
* ThreeReleaseSixValues: Three releases at various values of rhoR

## File names nomenclature

Files named *CT_N___.7z* store the data for the simulations con the calculation of suppression thresholds, with population size dictated by the digits following the *N* divided by 1000.

Each one of the files contained within these folders is named according to the following nomenclature:

* _ADM_Run[X|Y|Z].csv_: Adult males with experiment identifier dictated by:

    - X: Population size ID
    - Y: Resistant allele regeneration rate ID
    - Z: Experiment's repetition

* _AF1_Aggregate_Run[X|Y|Z].csv_:

    - X: Population size ID
    - Y: Resistant allele regeneration rate ID
    - Z: Experiment's repetition

* _AF1_Run[X|Y|Z].csv_:

    - X: Population size ID
    - Y: Resistant allele regeneration rate ID
    - Z: Experiment's repetition


Where the IDs are specified as follows:

    * X: Population size
        1: Population size of 1000
        2: Population size of 10000
        3: Population size of 100000

    * Y: Resistant allele regeneration rate
        - 01: .00001e-3
        - 02: .00003e-3
        - 03: .00005e-3
        - 04: .00007e-3
        - 05: .00009e-3
        - 06: .0001e-3
        - 07: .0003e-3
        - 08: .0005e-3
        - 09: .0007e-3
        - 10: .0009e-3
        - 11: .001e-3
        - 12: .003e-3
        - 13: .005e-3
        - 14: .007e-3
        - 15: .009e-3
        - 16: .01e-3
        - 17: .03e-3
        - 18: .05e-3
        - 19: .07e-3
        - 20: .09e-3
        - 21: .1e-3
        - 22: .3e-3
        - 23: .5e-3
        - 24: .7e-3
        - 25: .9e-3
        - 26: 1e-3
        - 27: 3e-3
        - 28: 5e-3
        - 29: 7e-3
        - 30: 9e-3
        - 31: 10e-3
        - 32: 30e-3
        - 33: 50e-3
        - 34: 70e-3
        - 35: 90e-3
        - 36: 100e-3
        - 37: 300e-3
        - 38: 500e-3
        - 39: 700e-3
        - 40: 900e-3

    * Z: Repetition number
        - Each experiment was repeated 40 times. The specific number of the iteration is denoted by these digits.

