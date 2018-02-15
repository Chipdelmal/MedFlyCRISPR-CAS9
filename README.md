# MedFlyCRISPR-CAS9

This repository contains the data for the MedFlies population suppression model described in the paper: PAPER TITLE HERE.

## Experiments Description

Two sets of experiments took place as part of this modeling endeavor, which will be described in the following subsections of this document.

## Files Descriptions

As part of the publication, we are releasing the dataset files used to generate our figures for reproducibility purposes.

<!--
as well as the code that describes our model.

### Simulation Files

Lorem ipsum
-->

### Datasets

The datasets are compressed in _7z_ format, which can be unarchived using free software such as: [Keka](http://www.kekaosx.com/en/) or [7zip](http://www.7-zip.org/download.html).

#### Population Dynamics Visual Inspection

Files named *VI_N___.7z* store the data for the visual inspection simulations with population size dictated by the digits following the *N* divided by 1000.

For the nomenclature of these filenames please refer to the following section of this document (as it describes in detail the key for the naming conventions).

#### Population Suppression Thresholds as Function of Resistant Allele Regeneration Rate

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

Where the IDs are specified as follows:

    * X: Population size
        1: Population size of 1000
        2: Population size of 10000
        3: Population size of 100000

    * Y: Resistant allele regeneration rate
        - 1:  .01e-3
        - 2:  .02e-3
        - 3:  .03e-3
        - 4:  .04e-3
        - 5:  .05e-3
        - 6:  .06e-3
        - 7:  .07e-3
        - 8:  .08e-3
        - 9:  .09e-3
        - 10: .1e-3
        - 11: .2e-3
        - 12: .3e-3
        - 13: .4e-3
        - 14: .5e-3
        - 15: .6e-3
        - 16: .7e-3
        - 17: .8e-3
        - 18: .9e-3
        - 19: 1e-3
        - 20: 2e-3
        - 21: 3e-3
        - 22: 4e-3
        - 23: 5e-3
        - 24: 6e-3
        - 25: 7e-3
        - 26: 8e-3
        - 27: 9e-3
        - 28: 10e-3
        - 29: 20e-3
        - 30: 30e-3
        - 31: 40e-3
        - 32: 50e-3

      * Z: Repetition number
        - Each experiment was repeated 52 times. The specific number of the iteration is denoted by these digits.

## Authors

 <!--* [Héctor M. Sánchez C.](chipdelmal.github.io), John M. Marshall-->
