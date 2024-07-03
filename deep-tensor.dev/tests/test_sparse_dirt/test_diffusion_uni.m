%%
root = pwd();
%
addpath(genpath([root '/fem_hp']));
addpath(genpath([root '/poisson']));

%%
test_case = 'scalar';
sq_flag = false;

% boundary conditions
tol = 1E-5;
%bnd_funcs = {@(p) find(abs(p(:,1))<tol), @(p) find(abs(p(:,1)-1)<tol)};
%bc_types  = {'essential', 'essential'};
%bc_funcs  = {@(x1,x2) 1, @(x1,x2) -1};

bnd_funcs = {@(p) find(abs(p(:,1))<tol), @(p) find(abs(p(:,1)-1)<tol)};
bc_types  = {'essential', 'essential'};
bc_funcs  = {@(x1,x2)-(cos(x2.*pi.*2)-1)/sqrt(2), @(x1,x2)cos(x2.*pi.*2)-1};

% flux function for qoi
flux_fx = @(x) uniform_fun (x, 0, 1);
flux_fy = @(y) triangle_fun(y, 0, 1, 'right');
flux = @(x1, x2) flux_fx(x1).*flux_fy(x2);

% obs locations
gx  = linspace(0, 1, 9);
gy  = linspace(0, 1, 9);
[xx0,yy0] = meshgrid(gx(2:end-1), gy(2:end-1));
xx0 = xx0(:);
yy0 = yy0(:);

obs_locs  = unique([xx0, yy0], 'rows');

%
model_opts = hp_options('h', 1/16, 'poly_order', 2, 'quad_order', 15, ...
    'xyratio', 1, 'qoi_func', flux, 'obs_locs', obs_locs, 'sq_param', sq_flag, ...
    'bnd_funcs', bnd_funcs,'bc_types', bc_types, 'bc_funcs', bc_funcs);

%inv_opts = inverse_options('s2n', 10, 'cov_type', 'MRF', 'mean', 0, ...
%    'gamma', 20, 'cond', [1, 1, 0], 'sigma', 1);

%inv_opts = inverse_options('s2n', 10, 'cov_type', 'GP', 'mean', 0, ...
%    'scale', 50, 'power', 2, 'sigma', 2);

inv_opts = inverse_options('s2n', 10, 'cov_type', 'GP', 'mean', 0, ...
    'scale', 1, 'power', -2, 'sigma', 1);

[model, obs, prior] = setup_poisson_uni(model_opts, inv_opts, 6, []);

solt = forward_solve(model, obs.true_u);

%%
figure
subplot(1,2,1)
trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), solt.state, 'edgecolor', 'none')
subplot(1,2,2)
u = reshape(obs.true_u, [], 1);
trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), u, 'edgecolor', 'none')

mlpost  = @(v) minus_log_post_uni(model, obs, prior, v);

%%

% temp1 = Tempering1(1);
%temp1 = Tempering1('min_beta', min(beta), 'ess_tol', 0.5, 'ess_tol_init', 0.5);  
temp1 = Tempering1('min_beta', 0.01, 'ess_tol', 0.5, 'ess_tol_init', 0.5);  
mask = [];
% temp1 = DataBatchAdapt(1:numel(Q_obs), 'ess_tol', 0.5, 'ess_tol_init', 0.5);

d = prior.dof;
p = 20;
%dom = BoundedDomain([-sqrt(3),sqrt(3)]);
dom = BoundedDomain([-1, 1]);
% dom = BoundedDomain([-5,5]);
% dom = LogarithmicMapping(2);
%dom = AlgebraicMapping(4);
base = ApproxBases(Legendre(p), dom, d);
% base = ApproxBases(Chebyshev2nd(p), dom, d);
diag = UniformReference(dom);

dirt_opt = DIRTOption('method', 'Aratio', 'defensive', 1E-8);

%%
tt_opt = TTOption('tt_method', 'amen', ...
    'als_tol', 0, 'local_tol', 1e-12, 'kick_rank', 1, 'max_rank', 10, 'max_als', 2, 'init_rank', 6);

tic;
tt_dirt = TTDIRT(mlpost, base, temp1, diag, tt_opt, dirt_opt); % , 'init_samples', init_samples1);
toc
