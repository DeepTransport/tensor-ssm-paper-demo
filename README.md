# TT-SSM
This git repo includes the basic and preconditioned versions of the tensor-train-based (TT) estimation approach in state-space models (SSM).

The tensor-train-based approach solves the filtering, smoothing and parameter estimation problem in state-space models. The approach relies on continuous tensor-train approximation and aims at the exact Bayesian solutions. Consequently, it does not suffer from particle dengeneracy which is the common problem in the popular particle-based methods. The sampling methods and preconditioning technique accompanying TT are also presented.

In our approach, we use the continuous TT package ftt.m-master developed by Cui and Dolgov (https://github.com/fastfins/ftt.m). The latest version is also published (https://github.com/fastfins/fastfins.m).

# Usage

* Run ftt.m_master/load_dir.m and ref/load_ref.m to add the path of TT and models.
* The filtering approach has an interface in TTinSSM/filter_main. The first section contains the parameters of the approach. The section section generate the synthetic data. The third section contains the parameterization of TT package. The origin setup in the file is for the basic version of TT in Kalman filter.
* The filtering results are stored in the structure 'model'.
* The smoothing approach has an interface in TTinSSM/smooth_main. The origin setup in the file is for the basic version of TT in Kalman filter.
* The smoothing results are stored in the structure 'stats'.

# Comments
* The synthetic data used in the paper for Kalman filter is stored in data50kalman.mat.
* The codes for the SMC$^2$ used as a counterpart of TT-based method in our work are presented in SMC2.
* Currently the predator and prey model simply applies Gaussian perturbation. However this can lead to negative state values. A predator and prey model on the log scale is under development, which solves the above problem.
