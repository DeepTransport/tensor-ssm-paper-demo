% the main file for FTT estimation in SSM
clear
clc

%%  model configuration
root = pwd();
cd ..
proot = pwd();
cd(root)
addpath([proot, '/models']) % all models and classes are stored here
addpath([proot, '/models/tensordot'])
% name = 'svmodels';
name = 'svmodels';


m = 1;
n = 1;

sp = load("SP500").SP500;

rsp = diff(log(sp.Close(end:-1:251)))';
rsp = (rsp - mean(rsp))/std(rsp);

T = length(rsp);


%% complete model
addpath([proot, '/models/', name])
% For SV, model.pre.para contains all the parameters
% They are (gamma, mu, tau, phi, a, delta, nu1, nu2) in order
% NaN indicates this parameter is to estimate by FTT

% set up parameter space
gamma = 0.95; % 0.95
tau = sqrt(3/64); % sqrt(3/64)
mu = 0; % -0.95
phi = 0; % -1/8
a = 0; % 0.03
delta = 0; % .2
nu1 = inf; % 25
nu2 = inf; % 25

pre.para = [gamma, tau, mu, phi, a, delta, nu1, nu2];

ind = [1,2,3];
d = length(ind);

myModel = ssmodel(name, d, m, n, T, pre);
myModel.pre.ind = ind;
d = length(myModel.pre.ind);
myModel = setup(myModel);

% complete model
rng(1);
myModel = complete(myModel);
rng('shuffle')
plot(myModel)

myModel.Y = rsp;

%% FTT configuration

N = 5e3; % default 5e3
tau = 10; % only used when process error is absent
sqr = 1;
% n_nodes = 50;
dom = BoundedDomain([-1, 1]);
dom2 = AlgebraicMapping(1);
poly = ApproxBases(Lagrangep(4, 8), dom, d + 2*m);
epd_Lag = 6;


opt = TTOption('tt_method', 'random', ...
    'als_tol', 1E-10, 'local_tol', 1E-4, 'max_rank', 20, 'max_als', 6, 'init_rank', 5, 'kick_rank', 5);
lowopt = TTOption('tt_method', 'random', ...
    'als_tol', 1E-10, 'local_tol', 1E-4, 'max_rank', 20, 'max_als', 2, 'init_rank', 5, 'kick_rank', 5);

rng(2)
mySol = full_sol(myModel, sqr, poly, ...
                opt, lowopt, N, epd_Lag);

mySol = solve(mySol);
rng('shuffle')
plot_sirt(mySol, T)

%%

[thetas, sams, w] = smooth(mySol, N, T);
plot_stats(mySol, thetas, sams, w);

rmpath([proot, '/models/', name])

