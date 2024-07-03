function J2du = matvec_J2du(model, sol, du)
%
% compute the Jacobian
% Tiangang Cui, 10/Jul/2023

if model.sq_param
    warning('there is a known bug with squared parameters');
end

dx = sol.dxdu.*du;
ndx = size(dx,2);
J2du = cell(ndx,1);

% adjoint states
t  = model.obs_operator';
lambda = zeros(size(sol.state,1), size(t,2));
lambda(sol.p,:) = -sol.L'\(sol.L\t(sol.p,:));

% Jacobian
J = zeros(model.n_sensors*model.n_datasets, length(sol.dxdu(:)));
for j = 1:model.n_datasets
    for k = 1:model.n_sensors
        J(k+(j-1)*model.n_sensors,:) = deri_adjoint_stiff_sol(model.mesh, model.local_elem, sol.kappa_type, sol.kappa, lambda(:,k), sol.state(:,j));
    end
end

dstate  = zeros(size(sol.state));
dlambda = zeros(size(lambda));
for i = 1:size(dx,2)
    J2du{i} = zeros(model.n_datasets*model.n_sensors, model.mesh.dof);

    % dstate = action of A(x)^(-1)*dA*u
    dk = reshape(dx(:,i), model.mesh.dof, []);
    %if model.sq_param
    %    ts = matvec_dstiff_sol_sq(model.mesh, model.local_elem, sol.kappa_type, sol.kappa, dk, sol.state);
    %else
    ts = matvec_dstiff_sol(model.mesh, model.local_elem, sol.kappa_type, sol.kappa, dk, sol.state);
    %end
    dstate(sol.p,:) = -sol.L'\(sol.L\ts(sol.p,:));

    % dlambda = action of lambda^T*dA*A(x)^(-1)
    for k = 1:model.n_sensors
        %if model.sq_param
        %    tl = matvec_dstiff_sol_sq(model.mesh, model.local_elem, sol.kappa_type, sol.kappa, dk, lambda);
        %else
        tl = matvec_dstiff_sol(model.mesh, model.local_elem, sol.kappa_type, sol.kappa, dk, lambda(:,k));
        %end
        dlambda(sol.p,k) = -sol.L'\(sol.L\tl(sol.p,:));
    end

    % action of lambda^T*dAdk*dstate + dlambda^T*dAdk*state
    for j = 1:model.n_datasets
        for k = 1:model.n_sensors
            %if model.sq_param
            %    g1 = deri_adjoint_stiff_sol_sq(model.mesh, model.local_elem, sol.kappa_type, sol.kappa, lambda(:,k), dstate(:,j));
            %    g2 = deri_adjoint_stiff_sol_sq(model.mesh, model.local_elem, sol.kappa_type, sol.kappa, dlambda(:,k), sol.state(:,j));
            %else
            g1 = deri_adjoint_stiff_sol(model.mesh, model.local_elem, sol.kappa_type, sol.kappa, lambda(:,k), dstate(:,j));
            g2 = deri_adjoint_stiff_sol(model.mesh, model.local_elem, sol.kappa_type, sol.kappa, dlambda(:,k), sol.state(:,j));
            %end
            J2du{i}(k+(j-1)*model.n_sensors,:) = (g1+g2)';
        end
    end

    J2du{i} = J2du{i}.*sol.dxdu(:)' + J.*sol.dx2du2(:)'.*du(:,i)';
end

end
