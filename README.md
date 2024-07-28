# About The Project
The tensor-ssm git repo includes the basic and preconditioned versions of the tensor-train-based (TT) estimation approach in state-space models (SSM).

The tensor-train-based approach solves the filtering, smoothing and parameter estimation problem in state-space models. The approach relies on continuous tensor-train approximation and aims at the exact Bayesian solutions. Consequently, it does not suffer from particle degeneracy which is the common problem in the popular particle-based methods. The sampling methods and the preconditioning technique accompanying TT are also presented.

In our approach, we use the continuous TT package deep-tensor developed by Cui and Dolgov (https://github.com/DeepTransport/deep-tensor) with the latest version published on (https://github.com/DeepTransport/deep-tensor).

The project is built on Matlab 2021a, and proves to work on Matlab 2023a. 

# Usage
The folders [deep-tensor.dev](https://github.com/DeepTransport/tensor-ssm-paper-demo/tree/main/deep-tensor.dev) and [models](https://github.com/DeepTransport/tensor-ssm-paper-demo/tree/main/models) are the toolboxes for continuous TT package and model configs. To implement the project, first run the file [deep-tensor.dev/load_dir.m](https://github.com/DeepTransport/tensor-ssm-paper-demo/blob/main/deep-tensor.dev/load_dir.m) to add the paths.

Then, the folders with prefix ``eg'' contain the files for four demonstrating examples in the paper "Tensor-train methods for sequential state and parameter learning in state-space models".: Kalman filter, stochastic volatility model, SIR model and predator-prey model. Simply run the script named ``main_script'' in each folder to check the outputs.

The output ``stats`` contains:
* ``stats.samples``: the simulated samples obtained from a series of tensor trains.
* ``stats.ess``: the effective sample size.
* ``stats.time_sample``: the computation time to generate samples.

# Config
The user can switch the model and change the parameters using the interface in the main script:
```
name = 'kalman';   % model names: 'kalman','pp', 'sv', 'sir'
T = 50;   % the steps for sequential inference
N = 1e4; % the sample size to estimate precondition in filter
```

The config for tensor-train-specific commands can also be modified in the main script:
```
dom = BoundedDomain([-1, 1]);
dom2 = AlgebraicMapping(1);

poly = ApproxBases(Lagrange1(32), dom, d + 2*m); 
poly2 = ApproxBases(Lagrangep(4,8), dom2, d + 2*m); 

opt = TTOption('tt_method', 'random', ...
    'als_tol', 1E-10, 'local_tol', 1E-4, 'max_rank', rank, 'max_als', 6, 'init_rank', 5, 'kick_rank', 5);
lowopt = TTOption('tt_method', 'random', ...
    'als_tol', 1E-10, 'local_tol', 1E-4, 'max_rank', rank, 'max_als', 2, 'init_rank', 5, 'kick_rank', 5);
```
We set the default domain for the codes to be $[-1,1]$, whereas the user can switch between the BoundedDomain or the AlgebraicMapping as above. The freedom of the polynomials can be modified in ``ApproxBases'' command. The rank of the tensor train, which is a crucial factor for accuracy and efficiency, can be modified in ``opt'' command.

The config for the built-in models can be changed in the script ``setup`` in each model folder in [models].
For example, the parameters for the Kalman filter example is stored in [models/kalman/setup](https://github.com/DeepTransport/tensor-ssm-paper-demo/blob/main/models/kalman/setup.m)

```
model.theta = norminv(([.8; .5]-0.4)/0.6);
model.pre.C = rand(model.n, model.m);
```
The default parameters are used to generate the results in the paper.


# Comments
* Currently, the predator-prey model simply applies Gaussian perturbation. However, this can lead to negative state values. A predator-prey model on the log scale is under development, which solves the above problem.
