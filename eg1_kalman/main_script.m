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
name = 'kalman';

d = 2;
m = 3;
n = 3;
T = 50;

%% complete model
addpath([proot, '/models/', name])
myModel = ssmodel(name, d, m, n, T);
myModel = setup(myModel);

myModel = complete(myModel);


%% FTT configuration

N = 5e3; % default 5e3
rank = 30;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sqr = 1;
dom2 = AlgebraicMapping(1);
dom = BoundedDomain([-1, 1]);
poly = ApproxBases(Lagrange1(32), dom, d + 2*m); 
poly2 = ApproxBases(Lagrangep(4,8), dom2, d + 2*m);  % use poly2 for best results in squared version (sqr = 1)
epd_Lag = 6;

opt = TTOption('tt_method', 'random', ...
    'als_tol', 1E-10, 'local_tol', 1E-4, 'max_rank', rank, 'max_als', 6, 'init_rank', 5, 'kick_rank', 5);
lowopt = TTOption('tt_method', 'random', ...
    'als_tol', 1E-10, 'local_tol', 1E-4, 'max_rank', rank, 'max_als', 2, 'init_rank', 5, 'kick_rank', 5);


mySolsqr = full_sol(myModel, sqr, poly, ...
                opt, lowopt, N, epd_Lag);
mySolsqr = solve(mySolsqr);
% plot_sirt(mySolsqr, T)

mySolsqr2 = full_sol(myModel, sqr, poly, ...
                opt, lowopt, N, epd_Lag);
mySolsqr2 = solve(mySolsqr2);
% plot_sirt(mySolsqr2, T)

sqr = 0;
mySolroot = full_sol(myModel, sqr, poly, ...
                opt, lowopt, N, epd_Lag);
mySolroot = solve(mySolroot);

% plot_sirt(mySolroot, T)

%% Compute the L1 distance
L1sqr = zeros(1, T);
L1root = zeros(1, T);

for t = 1:T
    f1 = @(x) theta_pdf(myModel, x, t);

    sirt = mySolsqr.SIRTs{t};
    f2 = @(x) eval_pdf(sirt, mySolsqr.L(1:2,1:2,t)\( norminv((x-0.4)/0.6) - mySolsqr.mu(1:2,t)))...
        ./ normpdf(norminv((x(1,:)-0.4)/0.6))./ normpdf(norminv((x(2,:)-0.4)/0.6));

    ftt = mySolroot.SIRTs{t};
    f3 = @(x) eval_pdf(ftt, mySolroot.L(1:2,1:2,t)\( norminv((x-0.4)/0.6) - mySolroot.mu(1:2,t)))...
        ./ normpdf(norminv((x(1,:)-0.4)/0.6))./ normpdf(norminv((x(2,:)-0.4)/0.6));

    L1sqr(t) = computeL1(f1, f2);
    L1root(t) = computeL1(f1, f3);
end

figure
plot(1:T, L1sqr, 1:T, L1root)
title("L1 error")
legend("Squared Version", "Non-squared Version")

[thetas, sams, w] = smooth(mySolsqr2, N, T);
plot_stats(mySolsqr2, thetas, sams, w);


rmpath([proot, '/models/', name])


