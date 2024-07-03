function dyJ2du = matvec_dyJ2du(model, sol, dy, du)
%
% compute the Jacobian
% Tiangang Cui, 10/Jul/2023

if model.sq_param
    warning('there is a known bug with squared parameters');
end

dx = sol.dxdu.*du;
dyJ2du = zeros(size(du));
% derivative of the linearization, acting on pertubations (on the right) and data misfit (on the left)
dy = reshape(dy, model.n_sensors, []);
ndy = size(dy,2)/model.n_datasets;
if ndy ~= 1
    error('number of misfit should be one for the action of second derivative');
end
t = model.obs_operator'*dy;  % for multiple RHS, size(dy, 2) = N_RHS
lambda = zeros(size(sol.state,1), ndy);
lambda(sol.p,:) = -sol.L'\(sol.L\t(sol.p,:));

dstate  = zeros(size(sol.state));
dlambda = zeros(size(lambda));
for i = 1:size(dx,2)
    % dstate = action of A(x)^(-1)*dA*u
    dk = reshape(dx(:,i), model.mesh.dof, []);
    %if model.sq_param
    %    ts = matvec_dstiff_sol_sq(model.mesh, model.local_elem, sol.kappa_type, sol.kappa, dk, sol.state);
    %else
    ts = matvec_dstiff_sol(model.mesh, model.local_elem, sol.kappa_type, sol.kappa, dk, sol.state);
    %end
    dstate(sol.p,:) = -sol.L'\(sol.L\ts(sol.p,:));

    % dlambda^T = action of lambda^T*dA*A(x)^(-1)
    %if model.sq_param
    %    tl = matvec_dstiff_sol_sq(model.mesh, model.local_elem, sol.kappa_type, sol.kappa, dk, lambda);
    %else
    tl = matvec_dstiff_sol(model.mesh, model.local_elem, sol.kappa_type, sol.kappa, dk, lambda);
    %tl = tl + model.obs_operator'*(sol.I*(model.obs_operator*dstate));
    %end
    dlambda(sol.p,:) = -sol.L'\(sol.L\tl(sol.p,:));

    % action of lamba^T*dAdk*dstate + dlambda^T*dAdk*state
    for j = 1:model.n_datasets
        %if model.sq_param
        %    g1 = deri_adjoint_stiff_sol_sq(model.mesh, model.local_elem, sol.kappa_type, sol.kappa, lambda(:,j), dstate(:,j));
        %    g2 = deri_adjoint_stiff_sol_sq(model.mesh, model.local_elem, sol.kappa_type, sol.kappa, dlambda(:,j), sol.state(:,j));
        %else
        g1 = deri_adjoint_stiff_sol(model.mesh, model.local_elem, sol.kappa_type, sol.kappa, lambda(:,j), dstate(:,j));
        g2 = deri_adjoint_stiff_sol(model.mesh, model.local_elem, sol.kappa_type, sol.kappa, dlambda(:,j), sol.state(:,j));
        %end
        dyJ2du(:,i) = g1 + g2;
    end

end
if model.exp_param
    dyJ2du = sol.dxdu.*dyJ2du + sol.gu.*du;
end
end
