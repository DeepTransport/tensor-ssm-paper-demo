function As = build_poisson_rom_matrices(model, U, Vs)
%set up the reduced order model for the 2nd order PDE
%
%inputs:
%  U:  nxr, the projection basis for the states
%  Vs: mx1, the projection basis for the tensor field
%      scalar: m = 1
%      vector: m = d
%      tensor: m = d+1

ne = model.mesh.num_elem;
nd = model.local_elem.num_dim;
nb = model.local_elem.num_node;
nq = model.local_elem.num_quad;
nU = size(U,2);

% evaluate grad of u for each i
gradu = cell(nd,1);
for di = 1:nd
    gradu{di} = zeros(nq*ne, nU);
end
%
for i = 1:nU
    tmp = reshape(U(model.mesh.node_map',i), nb, ne);
    for di = 1:nd
        for dj = 1:nd
            invJtij_duj = model.local_elem.grad_f_quad_pts{dj}*(tmp.*reshape(model.mesh.inv_Jt{di,dj},1,[]));
            gradu{di}(:,i) = gradu{di}(:,i) + reshape(invJtij_duj,[],1);
        end
    end
end
%
WdetJ = model.local_elem.quad_weights(:)*model.mesh.detJ(:)';
%
n = length(Vs);
if n == 1
    nK = size(Vs{1},2);
    As = cell(1);
    As{1} = zeros(nU^2, nK);
    for i = 1:nK
        k_atq = model.local_elem.f_quad_pts*reshape(Vs{1}(model.mesh.node_map',i), nb, ne);
        tmp = zeros(nU);
        for di = 1:nd
            tmp = tmp + (gradu{di}.*k_atq(:).*WdetJ(:))'*(gradu{di});
        end
        As{1}(:,i) = tmp(:);
    end
elseif n == nd && nd > 1
    As = cell(nd,1);
    for di = 1:nd
        nK = size(Vs{di},2);
        As{di} = zeros(nU^2, nK);
        for i = 1:nK
            k_atq = model.local_elem.f_quad_pts*reshape(Vs{di}(model.mesh.node_map',i), nb, ne);
            tmp = (gradu{di}.*k_atq(:).*WdetJ(:))'*(gradu{di});
            As{di}(:,i) = tmp(:);
        end
    end
elseif n == (nd+1)
    error('the cross term needs a quadratic assembly')
else
    error('kappa not implemented')
end


end