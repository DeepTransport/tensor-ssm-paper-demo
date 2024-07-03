function Qty = matvec_Qty(model, sol)

lambda = zeros(size(sol.state,1), 1);
lambda(sol.p,:) = -sol.L'\(sol.L\sol.Q(:,sol.p)');
lambda = lambda - model.phi';

if size(model.phi, 1) > 1 && size(sol.state, 2) > 1
    disp('not implemented');
end

if model.sq_param
    g = deri_adjoint_stiff_sol_sq(model.mesh, model.local_elem, sol.kappa_type, sol.kappa, lambda(:), sol.state(:));
else
    g = deri_adjoint_stiff_sol(model.mesh, model.local_elem, sol.kappa_type, sol.kappa, lambda(:), sol.state(:));
end

% transform to physical space x
Qty = sol.dxdu.*g;

end

