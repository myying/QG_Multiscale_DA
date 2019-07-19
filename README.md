Quasi-geostrophic model and Test code for multiscale ensemble data assimilation methods

Authors: Yue (Michael) Ying, NCAR/ASP
Toy EnKF system built by Yue Ying since 2016
QG model code adopted from Dr. Shafer Smith's webpage: https://cims.nyu.edu/~shafer/tools/index.html

References:
Ying, Y., F. Zhang, and J. L. Anderson, 2018: On the selection of localization radius in ensemble filtering for multiscale quasi-geostrophic dynamics. Mon. Wea. Rev., 146, 543-560.
Ying, Y., 2019: A multiscale alignment method for ensemble filtering with displacement errors. Mon. Wea. Rev., accepted.

Code description:

**config/localization/**

    Configuration files for the experiments conducted in Ying et al. 2018.

**config/alignment/**

    Configuration files for the experiments conducted in Ying 2019.

**model/**

    QG model code in Fortran

**namelist_input.sh**

    Generates parameter input file for QG model

**submit_job.sh**

    Batch job submission script

**generate_initial_condition.py**

    Generates initial condition for QG model

**run_truth.sh**

    Run script for generating truth (natrue run)

**make_obs.sh, generate_obs.py**

    Generates synthetic observations for OSSE

**plot_obs.py**

    Plots observation (Fig. 1 in Ying 2019)

**make_initial_ensemble.sh, generate_initial_ensemble**

    Generate initial ensemble with perturbations

**filter.py**

    Runs ensemble filter at a certain cycle

**run_cycle.sh**

    Perform the cycling data assimilation (top level run script)

**data_assimilation.py**

    Contains data assimilation algorithms

**util.py**

    Contains utility functions

**plot_ensemble.py**

    Plots the truth and ensemble state (spaghetti plots of contours, Fig. 3 in Ying 2019)

**calc_spectrum_ens/obs/truth.py, plot_spectrum.py**

    Calculate the error spectrum and plot figures (such as Fig. 4 in Ying et al. 2018)

**calc_errors.py, plot_errors.py**

    Calculate and plot domain-averaged errors (Fig. 4 in Ying 2019)

**demo_error_spectrum.py**

    Plots error spectrum from ensemble forecast, showing error growth behavior

**demo_scale_decomposition.py**

    Decompose a model state into scale components and display them

**demo_filter_multiscale.py**

    Performs Multiscale Alignment method (MSA; Ying 2019) for one cycle and plot intermediate stages
    of the model state at each scales (Fig. 2 in Ying 2019)

**demo_optical_flow.py**

    Demonstrate the calculated optical flows

**demo_random_noise.py**

    Demonstrate random noise field with different spectral distribution and their correlation function
