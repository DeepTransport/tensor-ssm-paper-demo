# About The Project
The tensor-ssm git repo includes the basic and preconditioned versions of the tensor-train-based (TT) estimation approach in state-space models (SSM).

The tensor-train-based approach solves the filtering, smoothing and parameter estimation problem in state-space models. The approach relies on continuous tensor-train approximation and aims at the exact Bayesian solutions. Consequently, it does not suffer from particle dengeneracy which is the common problem in the popular particle-based methods. The sampling methods and preconditioning technique accompanying TT are also presented.

In our approach, we use the continuous TT package ftt.m-master developed by Cui and Dolgov (https://github.com/fastfins/ftt.m) with the latest version published on (https://github.com/fastfins/fastfins.m).

The project is built on matlab 2021a, and proves to work on 2023a. 

# Usage
The folders [ftt.m-master](https://github.com/DeepTransport/tensor-ssm-paper-demo/tree/main/ftt.m-master) and [ref](https://github.com/DeepTransport/tensor-ssm-paper-demo/tree/main/ref) are the toolboxes for continuous TT package and model configs. To implement the project, first run the file [ftt.m-master/load_dir.m](https://github.com/DeepTransport/tensor-ssm-paper-demo/blob/main/ftt.m-master/load_dir.m) and [ref/load_ref.m](https://github.com/DeepTransport/tensor-ssm-paper-demo/blob/main/ref/load_ref.m) to add the paths.

The folder [TTinSSM](https://github.com/DeepTransport/tensor-ssm-paper-demo/tree/main/TTinSSM) is the main project of the TT estimation approach in SSM.
Run the script [TTinSSM/filter_main.m](https://github.com/DeepTransport/tensor-ssm-paper-demo/blob/main/TTinSSM/filter_main.m) to obtain the filtering result, stored in model_name.mat file in the folder ``data``. 
After storing the filtering result in model_name.mat file, input the model_name to the variable ``file`` in [TTinSSM/smooth_main.m](https://github.com/DeepTransport/tensor-ssm-paper-demo/blob/main/TTinSSM/smooth_main.m):
```
file = 'test';
```
run the script to obtain the smoothing result, stored in the structure 'stats' in the same model_name.mat file.

The output ``stats`` contains:
* ``stats.samples``: the simulated samples obtained from a series of tensor trains.
* ``stats.ess``: the effective sample size.
* ``stats.time_sample``: the computation time to generate samples.

# Config
The project includes the three models presented in the paper "Tensor-train methods for sequential state and parameter learning in state-space models".
The user can switch the model and change the parameters using the interface in the filter_main script:
```
model.name = 'kalman';   % model names: 'kalman','pp', 'sv'
model.pre = 'NL'; % two different options: NL and DIRT. DIRT stands for preconditioned approach.
model.steps = 20;   % the steps for sequential inference
                   % for cts model, the total time is model.dt * steps 
model.sam_size = 1e4; % the sample size to estimate precondition in filter
```

The config for ftt.m-master package can also be modified in the filter_main script:
```
poly1 = Legendre(40, [-4, 4]);
poly2 = Lagrange1(50, [-4, 4]);
poly3 = Lagrangep(4, 8, [0, 1], 'ghost_size', 1E-2, 'bc', 'Neumann');
poly4 = Fourier(40, [-4,4]);
poly5 = Lagrange1(60, [-4, 4], 'ghost_size', 1E-3, 'bc', 'Neumann');
poly6 = Lagrangep(2, 20, [-4, 4], 'ghost_size', 1E-3, 'bc', 'Neumann');
poly7 = Lagrange1(80, [-4, 4], 'ghost_size', 1E-3, 'bc', 'Neumann');
poly8 = Lagrangep(4, 8, [-4, 4], 'ghost_size', 1E-3, 'bc', 'Neumann');

opt1 = FTToption('tt_method', 'random', 'max_als', 6, 'als_tol', 1E-5, 'local_tol', 1E-5, 'kick_rank', 5, 'init_rank', 5, 'max_rank', 15, 'sqrt_flag',true);
opt2 = FTToption('tt_method', 'amen', 'max_als', 6, 'als_tol', 1E-5, 'local_tol', 1E-5, 'kick_rank', 5, 'init_rank', 5, 'max_rank', 15, 'sqrt_flag',true);
opt3 = FTToption('tt_method', 'amen', 'max_als', 6, 'als_tol', 1E-5, 'local_tol', 1E-5, 'kick_rank', 5, 'init_rank', 5, 'max_rank', 15, 'sqrt_flag',true);
opt4 = FTToption('tt_method', 'amen', 'max_als', 6, 'als_tol', 1E-5, 'local_tol', 1E-5, 'kick_rank', 5, 'init_rank', 5, 'max_rank', 15, 'sqrt_flag',true);
```

The config for the built-in models can be changed in the script ``truetheta`` in each model folder.
For example, the parameters for Kalman filter example is stored in [ref/kalman/truetheta](https://github.com/DeepTransport/tensor-ssm-paper-demo/blob/main/ref/kalman/truetheta.mlx)

```
model.dimension = 3; % dimension of the state
model.C = eye(model.dimension);
model.td = 2;  % dimension of the parameter
model.theta1 = 0.8;  % process error
model.theta2 = 0.5;   % observation error
model.theta = [model.theta1;model.theta2];
model.theta_transformed = norminv((model.theta -0.4)/0.6);
```
The default parameters are used to generate the results in the paper.


# Comments
* The codes for the SMC$^2$ used as a counterpart of TT-based method in our work are presented in the folder ``smc2``.
* Currently the predator and prey model simply applies Gaussian perturbation. However this can lead to negative state values. A predator and prey model on the log scale is under development, which solves the above problem.
