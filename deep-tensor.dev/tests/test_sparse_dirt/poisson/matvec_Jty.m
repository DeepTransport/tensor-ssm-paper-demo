function Jty = matvec_Jty(model, sol, dy)

dy = reshape(dy, model.n_sensors, []);
t  = model.obs_operator'*dy;  % for multiple RHS, size(dy, 2) = N * N_RHS
lambda = zeros(size(sol.state,1), size(dy,2));
lambda(sol.p,:) = -sol.L'\(sol.L\t(sol.p,:));

N   = size(dy,2)/model.n_datasets;
Jty = zeros(length(sol.dxdu(:)), N);

for k = 1:N
    bi = (k-1)*model.n_datasets;
    for j = 1:model.n_datasets
        if model.sq_param
            g = deri_adjoint_stiff_sol_sq(model.mesh, model.local_elem, sol.kappa_type, sol.kappa, lambda(:,bi+j), sol.state(:,j));
        else
            g = deri_adjoint_stiff_sol(model.mesh, model.local_elem, sol.kappa_type, sol.kappa, lambda(:,bi+j), sol.state(:,j));
        end
        Jty(:,k) = Jty(:,k) + g;
    end
end

% transform to physical space x
Jty = sol.dxdu.*Jty;

end

