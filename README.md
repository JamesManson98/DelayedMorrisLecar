# Delayed_Morris_Lecar
Numerical bifurcation analysis of oscillations in a self-coupled neuron with delay

This code runs in MATLAB R2016b and requires DDEBiftoolv3.1.1 and the extensions DDEBiftool_nmfm and DDEBiftool_extra_psol.

Add 'Morris-Lecar' folder to demos folder in DDEBiftool.

mainscript: One and two parameter bifurcation diagrams, followed by continuation of the curve of Hopf bifurcations and computation of branches of fold bifurcations of periodic orbits and torus bifurcations. Then the full bifurcation diagram is plotted.

behaviour_of_system: Plots of the  behaviour of the delayed coupling.

br_plot3, br_splot3, br_splot, br_splot_l1, p_dststb_l1: Functions to assist with plotting.

complete_bifurcation_plotting: The complete bifurcation diagram in the (tau,kappa) plane.

dde23_POs: Time series plots of periodic orbits computed using dde23.

floquet_multipliers: Records a GIF of the floquet multipliers of a branch of periodic orbits and saves in the current directory.

Folds_of_POs: Computation of the fold bifurcations of periodic orbits.

Initialise_PO_Fold: Finds branches of periodic orbits to start the continuation of the fold bifurcations of periodic orbits.

Initialise_Torus: Finds branches of periodic orbits to start the continuation of torus bifurcations.

Kappa_15_Bistable_Taus, Kappa_45_Bistable_Taus, Kappa_85_Bistable_Taus: Plotting phase portraits and time series at multistable values of tau for kappa = 15, 45, and 85 respectively.

Long_contination_kappa_tau: The continuation of the curve of Hopf bifurcations up to Kappa=600.

One_Kappa_POs_15, One_Kappa_POs_45, One_Kappa_POs_85: Plotting of branches of periodic orbits emerging from the Hopf points close to kappa = 15, 45, and 85 respectively.

One_par_I: One parameter bifurcation diagram in the (I,V) plane.

POs_check: Checking the criticality of the Hopf bifurcations along the curve of Hopf bifurcations.

POs_near_swallowtail: Checking that the swallowtail exists and is not a numerical error.

Stability_large_kappa: Computation and plotting of the stability of the curve of Hopf bifurcations up to Kappa=600.

Torus_bifurcations: Computation of torus bifurcations.

Two_par_I_kappa: Computation of the two parameter bifurcation diagrams and corresponding 3-dimensional plots in the (I,kappa,V) plane. The curve of Hopf bifurcations are also computed. 

Two_par_I_tau: Computation of the two parameter bifurcation diagrams and corresponding 3-dimensional plots in the (I,tau,V) plane. The curve of Hopf bifurcations are also computed. 

Two_par_kappa_tau: Computation of the two parameter bifurcation diagrams and corresponding 3-dimensional plots in the (kappa,tau,V) plane. The curve of Hopf bifurcations are also computed. 


