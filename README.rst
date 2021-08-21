==================================
Bayesian Optimisation via CosmoMC
==================================
:Author: Julius Wons
:Contact: j.wons@unsw.edu.au 

Description
============

This code is a fork to CosmoMC and replaces the MCMC sampling alogorithm with a Bayesian Optimisation.
The code finds the global best fit for modulations of the primordial power spectrum efficiently.

Installation and tests
================

Installation is similar to CosmoMC see the `ReadMe <http://cosmologist.info/cosmomc/readme.html>`_.

Requires installation of the  `unbinned Planck likelihood <https://pla.esac.esa.int/#cosmology>`_ (bin1). See `Planck Readme <https://cosmologist.info/cosmomc/readme_planck.html>`_.
Requirements for the compiler are the same as for CosmoMC.
The code was tested for intel 18.0.1 or newer compiler. Technically should work for gcc compiler

After downloading the code, go to ./source and and compile the code using
::
  make
After that go back to base dictonary. There are two test files.
First, test.ini uses likelihood native to CosmoMC. You do not need the Planck likelihoods for this one.
Running ::
  ./cosmomc test.ini
will run a quick example with logarithmic oscillations. This should take around 1-2 hours and requires around 2GB memory. 
If it does not crash immediately, then the installation should be fine.
For a quick check, change in test.ini action=3 to action=4. This will quickly check the likelihood.

When the unbinned planck likelihood has been installed, you can test it by running
::
  ./cosmomc test_planck_lin.ini
This will run an example with linear oscillations using the Planck likelihood code. This will need around 15-20GB
of memory and takes around 10-12 hours. The value of the best fit should be around 6.25. The data will be saved in
./Output/planck_test.txt and can be plotted with plot_test.py. The results should look similar to test.pdf

Implementation of other modulations of the primordial power spectrum
=====================================================================
Four different modulations to the primordial power spectrum are already implemented.
Linear oscillations, logarithmic oscillations, logarithmic oscillations with running frequency,
and three different versions of the primordial standard clock.

To implement your own model, simply overwrite one of the models in /camb/fortran/InitialPower.f90
Five parameters are already implemented: AmpOsc, linfreq, phase, newP4, newP5
These can be used in your own models. Afterwards you need to include prior ranges for each parameter in
the .ini file.


Algorithm details
==================

See the latest `paper <http://arxiv.org/abs/1304.4473>`_.
