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
name = 'sv';

d = 2;
m = 1;
n = 1;
T = 1000;


%% complete model

addpath([proot, '/models/', name])
myModel = ssmodel(name, d, m, n, T);
myModel = setup(myModel);

rng(1); 
myModel = complete(myModel);
rng('shuffle')

%% FTT configuration

N = 5e3; % default 5e3
tau = 10; % only used when process error is absent
sqr = 1;
% n_nodes = 50;
dom = BoundedDomain([-1, 1]);
dom2 = AlgebraicMapping(1);
poly = ApproxBases(Lagrangep(4, 8), dom, d + 2*m);
poly2 = ApproxBases(Lagrangep(4, 8), dom2, d + 2*m);
epd_Lag = 6;


opt = TTOption('tt_method', 'random', ...
    'als_tol', 1E-10, 'local_tol', 1E-4, 'max_rank', 10, 'max_als', 6, 'init_rank', 5, 'kick_rank', 5);
lowopt = TTOption('tt_method', 'random', ...
    'als_tol', 1E-10, 'local_tol', 1E-4, 'max_rank', 10, 'max_als', 2, 'init_rank', 5, 'kick_rank', 5);

rng(2)
mySol = full_sol(myModel, sqr, poly2, ...
                opt, lowopt, N, epd_Lag);

mySol = solve(mySol);
rng('shuffle')
plot_sirt(mySol, T)

%%
NN = 1e4;
[thetas, sams, w] = smooth(mySol, N, T);
plot_stats(mySol, thetas, sams, w);

rmpath([proot, '/models/', name])


