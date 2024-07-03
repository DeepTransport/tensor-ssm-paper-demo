function stiff = assemble_stiff_sq(mesh, local_elem, b_grad_x, fill_i, kappa, fill_d)
%assemble_stiff
%
% assemble the stiffness matrix for a tensor field viscosity/permeability
% three forms for permeabilty: 
%   scalar: kappa is nx1, the square root of the scalar field
%   vector: kappa is nxd, the square root of the diagonal of the field
%   tensor: kappa is nx(d+1), each row forms K = gamma I + b x b
%           the first col is the sqaure root of the scalar field gamma, 
%           the rest of cols are for the vector field b
%           
% Tiangang Cui, 02/Sep/2020

% 
ne = mesh.num_elem;
nd = local_elem.num_dim;
nb = local_elem.num_node;
nq = local_elem.num_quad;


% bring the kappa in
[m,n] = size(kappa);
if m ~= mesh.dof && m ~= 1
    error('kappa with wrong dimension')
end

if n == 1
    % for any d
    % evaluate sqrt(K) grad(phi), without quadrature and detJ
    grad = cell(nd, 1);
    for di = 1:nd
        if m == 1
            grad{di} = b_grad_x{di}*kappa;
        else
            loc_ka   = local_elem.f_quad_pts*kappa(mesh.node_map');
            grad{di} = b_grad_x{di}.*repmat(loc_ka, nb, 1);
        end
    end
    n_grad = nd;
elseif n == nd && nd > 1
    % for any d > 1
    % evaluate sqrt(K) grad(phi), without quadrature and detJ
    grad = cell(nd, 1);
    for di = 1:nd
        if m == 1
            grad{di} = b_grad_x{di}*kappa(di);
        else
            tmp = kappa(:,di);
            loc_ka   = local_elem.f_quad_pts*tmp(mesh.node_map');
            grad{di} = b_grad_x{di}.*repmat(loc_ka, nb, 1);
        end
    end
    n_grad = nd;
elseif n == (nd+1)
    % for 2d
    % the diagonal part, using sqrt
    grad = cell(nd+1, 1);
    if m == 1
        for di = 1:nd
            grad{di} = b_grad_x{di}*kappa(1);
        end
    else
        tmp = kappa(:,1);
        loc_ka = local_elem.f_quad_pts*tmp(mesh.node_map');
        for di = 1:nd
            grad{di} = b_grad_x{di}.*repmat(loc_ka, nb, 1);
        end
    end
    % the cross part
    grad{nd+1} = zeros(nq*nb, ne);
    for di = 1:nd
        if m == 1
            grad{nd+1} = grad{nd+1} + b_grad_x{di}*kappa(1+di);
        else
            tmp = kappa(:,1+di);
            loc_ka  = local_elem.f_quad_pts*tmp(mesh.node_map');
            grad{nd+1} = grad{nd+1} + b_grad_x{di}.*repmat(loc_ka, nb, 1);
        end
    end
    n_grad = nd+1;
else
    error('kappa not implemented')
end

% evaluate the dot product
fill = zeros(nb^2, ne);
for ei = 1:ne
    loc = zeros(nb, nb);
    for di = 1:n_grad
        tmp = reshape(grad{di}(:,ei), nq, nb);
        loc = loc + tmp'*(local_elem.quad_weights(:).*tmp);
    end
    fill(:,ei) = loc(:)*mesh.detJ(ei);
end

stiff = accumarray(fill_i, fill(:), fill_d, [], [], true);

end
