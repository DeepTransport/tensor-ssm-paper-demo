function Ju = matvec_Ju(model, sol, du)

% transform to physical space x
dx = sol.dxdu.*du;
Ju = zeros(model.n_sensors, model.n_datasets*size(dx,2));
dstate = zeros(size(sol.state));

for i = 1:size(dx,2)
    dk = reshape(dx(:,i), model.mesh.dof, []);
    
    if model.sq_param
        t = matvec_dstiff_sol_sq(model.mesh, model.local_elem, sol.kappa_type, sol.kappa, dk, sol.state);
    else
        t = matvec_dstiff_sol(model.mesh, model.local_elem, sol.kappa_type, sol.kappa, dk, sol.state);
    end
    dstate(sol.p,:) = -sol.L'\(sol.L\t(sol.p,:));
    ind = (1:model.n_datasets) + (i-1)*model.n_datasets;
    Ju(:,ind) = model.obs_operator*dstate;
end

Ju = reshape(Ju, [], size(dx,2));

end