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
name = 'sir_austria';

d = 0;
m = 18;
n = m/2;
T = 20;

%% complete model
addpath([proot, '/models/', name])

myModel = ssmodel(name, d, m, n, T);
myModel = setup(myModel);

% complete model
rng(1);
myModel = complete(myModel);
rng('shuffle')
plot(myModel)

flag_y = any(isnan(myModel.Y));
if flag_y == 1
    error("Sythetic data contain NaN.")
end


%% FTT configuration

N = 5e3; % default 5e3
tau = 10; % only used when process error is absent
sqr = 1;

dom = BoundedDomain([-1, 1]);
poly1 = ApproxBases(Lagrangep(4, 8), dom, d + 2*m);
poly2 = ApproxBases(Lagrangep(4, 8),  AlgebraicMapping(1), d + 2*m);
epd_Lag = 6;

opt = TTOption('tt_method', 'random', ...
    'als_tol', 1E-10, 'local_tol', 1E-4, 'max_rank', 40, 'max_als', 8, 'init_rank', 20, 'kick_rank', 5);
lowopt = TTOption('tt_method', 'random', ...
    'als_tol', 1E-10, 'local_tol', 1E-4, 'max_rank', 40, 'max_als', 2, 'init_rank', 20, 'kick_rank', 5);

rng(2)
mySol_unbounded = full_sol(myModel, sqr, poly2, ...
                opt, lowopt, N, 4);
mySol_unbounded = solve(mySol_unbounded);

rng('shuffle')

% plot_sirt(mySol_unbounded, T)

%%
N_sample = 1e4;
[thetas, sams, w, logpdf_eall] = smooth(mySol_unbounded, N_sample, T);
stats = plot_stats(mySol_unbounded, thetas, sams, w);


rmpath([proot, '/models/', name])


