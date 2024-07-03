function [model, obs, prior] = setup_poisson_uni(model_opts, inv_opts, nKL, output)

% Tiangang Cui, 10/Sep/2020

[p,t]   = simplemesh2d(model_opts.h, 1, model_opts.xyratio);
model   = setup_2nd_order(p, t, model_opts);

%
model.gmres_flag  = model_opts.gmres;
model.res_tol = model_opts.res_tol;
model.explicit_ja = false;

% assemble forcing term
force   = model_opts.force_func(model.mesh.nodes(:,1), model.mesh.nodes(:,2));
% weak form of the force, on nodes
if length(force) == 1
    model.b = full(sum(model.mass, 2))*force;
else
    model.b = model.mass*force(:);
end

% assemble boundary condition
model.bnd_b = zeros(model.mesh.dof, 1);
model.bmass = spalloc(model.mesh.dof,model.mesh.dof,0);
for bi = 1:length(model_opts.bnd_funcs)
    b_ind   = model_opts.bnd_funcs{bi}(model.mesh.nodes);
    b_force = model_opts.bc_funcs{bi}(model.mesh.nodes(:,1),model.mesh.nodes(:,2));
    if length(b_force) == 1
        b_force = ones(size(b_ind(:)))*b_force;
    else
        b_force = b_force(b_ind);
    end
    bf_weak = model.mass_bnds{bi}*sparse(b_ind(:),ones(size(b_ind(:))),b_force(:),model.mesh.dof, 1);
    switch model_opts.bc_types{bi}
        case {'flux'}
            disp('double check non-zero flux b.c.')
            model.bnd_b = model.bnd_b + bf_weak;
        case {'essential'}
            model.bmass = model.bmass + model.mass_bnds{bi}*model.penalty;
            model.bnd_b = model.bnd_b + bf_weak*model.penalty;
        case {'mixed'}
            disp('Mixed b.c. is not implemented')
    end
end

model.b = model.b + model.bnd_b;

% apply observation operator
model.obs_operator  = mesh_interpolate(model.mesh, model.local_elem, model_opts.obs_locs, model.h, false);
model.n_sensors     = size(model_opts.obs_locs,1);
model.n_datasets    = 1;

% apply qoi function
tol = 1E-10;
disp('flux QoI')
model.qoi_flag = false;
if ~isempty(model_opts.qoi_func)
    phi = reshape(model_opts.qoi_func(model.mesh.nodes(:,1),model.mesh.nodes(:,2)), 1, []);
    phi(abs(phi)<tol) = 0;
    model.phi = sparse(phi);
    model.qoi_flag  = true;
end

model.dof = model.mesh.dof; % for consistency
% transformation
model.exp_param = model_opts.exp_param;
model.sq_param  = model_opts.sq_param;


%%%%%%%%%%%%%%%%%%% build prior %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prior = make_prior_gp  (model.mesh, inv_opts.scale, inv_opts.power, inv_opts.sigma);
prior.mesh_dof = model.mesh.dof;
%
disp('set zero mean')
if length(inv_opts.mean) == 1
    prior.mean_u = 0*ones(prior.mesh_dof,1);
else
    prior.mean_u = zeros(size(inv_opts.mean));
end
%
[V, d]          = prior_cov_eig(prior);
prior.dof       = ceil(nKL);
jnd             = 1:prior.dof;
prior.P         = V(:, jnd);
S               = d(jnd).^(0.5);
prior.basis     = prior.P.*S(:)';
prior.basis_w   = prior.P.*(S(:)'.^(-1));
prior.mean_v    = zeros(prior.dof, 1);
prior.d         = d;
prior.chol2w    = prior_cov_lt(prior, prior.basis_w);
prior.w2chol    = prior_cov_il(prior, prior.basis);

prior.type      = 'Basis';
prior.note      = 'KL';

mm = 0;
for i = 1:nKL
    mm = mm + abs(prior.basis(:,i));
end
mm = ceil(max(mm));

model.exp_thres = mm;

%%%%%%%%%%%%%%%%%%% load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obs.jacobian = inv_opts.jacobian;
obs.like = inv_opts.like;
%
if  ~isempty(output) && exist(output,'file')
    tmp         = load(output);
    obs.data    = tmp.data;
    obs.std     = tmp.std;
    obs.true_u  = tmp.true_u;
else
    r = rand(prior.dof, 1)*2-1;
    true_u  = matvec_prior_L(prior, r) + prior.mean_u;
    soln    = forward_solve(model, true_u);  % reference solution
    
    % the s.t.d. is calculated from the signal to noise ratio
    if inv_opts.s2n > 0
        std = mean(abs(soln.d(:)))/inv_opts.s2n;
    else
        std = inv_opts.std;
    end
    
    % generate data
    data    = soln.d + randn(model.n_sensors, model.n_datasets)*std;
    
    obs.data    = data;
    obs.std     = std;
    obs.true_u  = true_u;
end

if  ~isempty(output) && ~exist(output,'file')
    save(output, 'data', 'std', 'true_u');
end

obs.n_sensors    = model.n_sensors;
obs.n_datasets   = model.n_datasets;
obs.n_data       = obs.n_sensors*obs.n_datasets;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

