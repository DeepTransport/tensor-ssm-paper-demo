function J = explicit_jacobian(model, sol)
%
% compute the Jacobian
% Tiangang Cui, 10/Mar/2013

lambda = zeros(size(sol.state,1), model.n_sensors);
lambda(sol.p,:) = -sol.L'\(sol.L\model.obs_operator(:,sol.p)');

J = zeros(model.n_sensors*model.n_datasets, length(sol.dxdu(:)));

for k = 1:model.n_datasets
    bi = (k-1)*model.n_sensors;
    for j = 1:model.n_sensors
        if model.sq_param
            J(bi+j,:) = deri_adjoint_stiff_sol_sq(model.mesh, model.local_elem, sol.kappa_type, sol.kappa, lambda(:,j), sol.state(:,k));
        else
            J(bi+j,:) = deri_adjoint_stiff_sol(model.mesh, model.local_elem, sol.kappa_type, sol.kappa, lambda(:,j), sol.state(:,k));
        end
    end
end

% transform to physical space x
J = J.*sol.dxdu(:)';

end
