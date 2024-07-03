function prior = make_prior_mrf(model, gamma, cond, sigma)
%UPDATE_MRF
%
% Q u = w, Q = inv(sigma) * (K(cond) + k I)
%
% Tiangang Cui, 09/May/2014

if length(cond) == 3 || length(cond) == 2 % for homogeneous case
    cond = cond(:)';
end

prior.stiff = assemble_stiff(model.mesh, model.local_elem, model.grad, model.fill_i, cond, model.fill_d);
prior.M = model.mass*sigma^2;
prior.K = prior.stiff + model.mass*gamma;

%Q = K*inv(M)*K

[prior.Lk,~,prior.pk] = chol(prior.K, 'lower', 'vector');
[prior.Lm,~,prior.pm] = chol(prior.M, 'lower', 'vector');

prior.dof       = model.mesh.dof;
prior.cov_type  = 'MRF';
prior.type      = 'Field';

end