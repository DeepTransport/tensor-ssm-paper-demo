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
name = 'pp';

d = 6;
m = 2;
n = 2;
T = 20;
% theta = 1;

%% complete model
addpath([proot, '/models/', name])

myModel = ssmodel(name, d, m, n, T);
myModel = setup(myModel);


% complete model
rng(1); 
myModel = complete(myModel);
rng('shuffle')
plot(myModel)

opt = odeset('RelTol',1e-2);
[t_temp, x_temp] = ode45(@(t, x) odefun(t, x, myModel.pre.theta), [0, T*myModel.pre.dt], myModel.pre.init, opt);
figure
plot(t_temp, x_temp(:, 1))
title("ODE for prey")
figure
plot(t_temp, x_temp(:, 2))
title("ODE for predator")


%% FTT configuration

N = 1e4; % default 5e3
tau = 10; % only used when process error is absent
sqr = 1;
dom = BoundedDomain([-1, 1]);
poly1 = ApproxBases(Lagrangep(4, 8), dom, d + 2*m);
poly2 = ApproxBases(Lagrangep(4, 8),  AlgebraicMapping(1), d + 2*m);
epd_Lag = 5;


%%
precond.c = 0.4;
epd2 = 1;
precond.opt = "pifg"; % three options "pifg", "fg", "g" 
tail = 3;
precond.R = @(x) tg_cdf(tail*x, epd2 * tail);
precond.Rinv = @(x) tail\tg_inv(x, epd2 * tail);
precond.r = @(x) tail*tg_pdf(tail*x, epd2 * tail);
precond.logr = @(x) log(tail) + tg_logpdf(tail*x, epd2 * tail);

opt = TTOption('tt_method', 'random', ...
    'als_tol', 1E-10, 'local_tol', 1E-4, 'max_rank', 20, 'max_als', 4, 'init_rank', 10, 'kick_rank', 10);
lowopt = TTOption('tt_method', 'random', ...
    'als_tol', 1E-10, 'local_tol', 1E-4, 'max_rank', 20, 'max_als', 2, 'init_rank', 10, 'kick_rank', 10);

rng(2)
mySol_bounded = full_sol(myModel, sqr, poly1, ...
                opt, lowopt, N, epd_Lag);
mySol_bounded = solve(mySol_bounded);



rng(2)
mySol_pre = pre_sol(myModel, poly1, ...
    opt, lowopt, N, epd_Lag, precond);
mySol_pre = solve(mySol_pre);

rng('shuffle')



%%
N_sample = 1e4;

[thetas, sams, w] = smooth(mySol_bounded, N_sample, T);
plot_sirt(mySol_bounded, T)
stats = plot_stats(mySol_bounded, thetas, sams, w);

[thetas, sams, w] = smooth(mySol_pre, N_sample, T);
plot_sirt(mySol_pre, T)
stats = plot_stats(mySol_pre, thetas, sams, w);

rmpath([proot, '/models/', name])

