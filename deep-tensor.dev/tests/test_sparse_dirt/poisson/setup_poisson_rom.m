function rom = setup_poisson_rom(model, prior_redu, sub_vs, states, weights, kappa_type, rom_opts)
%set up the reduced order model for the 2nd order PDE
%
%inputs:
%  U:  nxr, the projection basis for the states
%  Vs: mx1, the projection basis for the tensor field
%      scalar: m = 1
%      vector: m = d
%      tensor: m = d+1

rom.sq_param = model.sq_param;

% build reduced basis for the states
% U needs to be orthogonal w.r.t. the mass matrix
% U'*M*U = I, M = LL' => U'*L*L'*U = I => L'*U is orthogonal
%
[U,S] = svd((model.mass_L'*states).*sqrt(weights(:)'), 'econ');
s  = diag(S);
nU = truncate_energy(s(2:end).^2, rom_opts.pod_state_tol) + 1; % take out the first basis from the energy
U  = model.mass_L'\U(:,1:nU);
rom.dof = nU;
%
nd = model.local_elem.num_dim;

% build reduced basis for the params
n = length(prior_redu.bases);
switch kappa_type
    case {'scalar'}
        rom.es = cell(1);
        rom.Ks = cell(1);
        rom.redu_us = cell(1);
        rom.redu_means = cell(1);
        Vs = cell(1);
        xs = exp(prior_redu.bases{1}*sub_vs{1} + prior_redu.mean_us{1}(:));
        [V,S,~] = svd((model.mass_L'*xs).*sqrt(weights(:)'), 'econ');
        V  = model.mass_L'\V;
        d  = diag(S).^2;
        ii = d >= d(2)*rom_opts.prior_deim_tol;
        nP = sum(ii);
        [rom.es{1},~,rom.Ks{1}] = deim(V, nP, rom_opts.deim_reg_factor);
        rom.redu_us{1} = prior_redu.bases{1}(rom.es{1},:);
        rom.redu_means{1} = prior_redu.mean_us{1}(rom.es{1});
        Vs{1} = V(:,1:nP);
    case {'vector'}
        rom.es = cell(nd, 1);
        rom.Ks = cell(nd, 1);
        rom.redu_us = cell(nd, 1);
        rom.redu_means = cell(nd, 1);
        Vs = cell(nd, 1);
        %
        k = 0;
        for di = 1:nd
            %
            % build indices to the input reduced dimensional parameter vector
            rv = size(prior_redu.bases{di}, 2);
            rom.redu_param_ind{di} = (1:rv) + k;
            k = k+rv;
            %
            xs = exp(prior_redu.bases{di}*sub_vs{di} + prior_redu.mean_us{di}(:));
            [V,S,~] = svd((model.mass_L'*xs).*sqrt(weights(:)'), 'econ');
            V  = model.mass_L'\V;
            d  = diag(S).^2;
            ii = d >= d(2)*rom_opts.prior_deim_tol;
            nP = sum(ii);
            [rom.es{di},~,rom.Ks{di}] = deim(V, nP, rom_opts.deim_reg_factor);
            rom.redu_us{di} = prior_redu.bases{di}(rom.es{di},:);
            rom.redu_means{di} = prior_redu.mean_us{di}(rom.es{di});
            Vs{di} = V(:,1:nP);
        end
    case {'tensor'}
        rom.es = cell(1);
        rom.Ks = cell(1);
        rom.redu_us = cell(1);
        rom.redu_means = cell(1);
        Vs = cell(nd+1, 1);
        %
        rom.redu_param_ind{1} = 1:size(prior_redu.bases{1}, 2);
        rom.redu_param_ind{2} = (1:size(prior_redu.bases{2}, 2)) + size(prior_redu.bases{1}, 2);
        %
        xs = exp(prior_redu.bases{1}*sub_vs{1} + prior_redu.mean_us{1}(:));
        [V,S,~] = svd((model.mass_L'*xs).*sqrt(weights(:)'), 'econ');
        V  = model.mass_L'\V;
        d  = diag(S).^2;
        ii = d >= d(2)*rom_opts.prior_deim_tol;
        nP = sum(ii);
        [rom.es{1},~,rom.Ks{1}] = deim(V, nP, rom_opts.deim_reg_factor);
        rom.redu_us{1} = prior_redu.bases{1}(rom.es{1},:);
        rom.redu_means{1} = prior_redu.mean_us{1}(rom.es{1});
        Vs{1} = V(:,1:nP);
        %
        dof = model.mesh.dof;
        for di = 1:nd
            ind = (1:dof) + (di-1)*dof;
            Vs{di+1} = prior_redu.bases{2}(ind,:);
        end
    otherwise
        error('kappa not implemented')
end

rom.kappa_type = kappa_type;

if rom.sq_param
    [rom.As, rom.inds] = build_poisson_rom_matrices_sq(model, U(:,1:nU), Vs);
elseif (n == nd) || (n == 1)
    rom.As = build_poisson_rom_matrices(model, U, Vs);
else
    disp('cross term needs a quadratic rom')
end

% other matrices, rhs, and observation operators
rom.A = U'*( (model.bmass + model.exp_thres*model.stiff)*U );
rom.b = U'*model.b;
rom.obs_operator = model.obs_operator*U;
rom.qoi_flag = model.qoi_flag;
%
% one should build phi as a reduced basis
if rom.qoi_flag
    rom.phi = model.phi*U;
end

rom.states = U;

end