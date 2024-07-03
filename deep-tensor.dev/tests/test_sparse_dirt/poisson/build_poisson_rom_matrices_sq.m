function [As, inds] = build_poisson_rom_matrices_sq(model, U, Vs)
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
    As = cell(1);
    inds = cell(1);
    % for any d
    nK = size(Vs{1},2);
    K  = zeros(nq*ne, nK);
    for i = 1:nK
        k_atq = model.local_elem.f_quad_pts*reshape(Vs{1}(model.mesh.node_map',i), nb, ne);
        K(:,i) = reshape(k_atq,[],1);
    end
    % assemble squared K basis
    nKs = (nK+1)*nK/2;
    As{1} = zeros(nU^2, nKs);
    inds{1} = zeros(nKs, 1);
    k = 0;
    for i = 1:nK
        for j = i:nK
            tmp = zeros(nU);
            for di = 1:nd
                if i == j
                    tmp = tmp + (gradu{di}.*K(:,i).*WdetJ(:))'*(gradu{di}.*K(:,j));
                else
                    tmp = tmp + (gradu{di}.*K(:,i).*WdetJ(:))'*(gradu{di}.*K(:,j));
                    tmp = tmp + (gradu{di}.*K(:,j).*WdetJ(:))'*(gradu{di}.*K(:,i));
                end
            end
            k = k + 1;
            As{1}(:,k) = tmp(:);
            inds{1}(k) = (j-1)*nK + i;
        end
    end
elseif n == nd && nd > 1
    As = cell(nd,1);
    inds = cell(nd,1);
    for di = 1:nd
        nK = size(Vs{di},2);
        K  = zeros(nq*ne, nK);
        for i = 1:nK
            k_atq = model.local_elem.f_quad_pts*reshape(Vs{di}(model.mesh.node_map',i), nb, ne);
            K(:,i) = reshape(k_atq,[],1);
        end
        % assemble squared K basis
        nKs = (nK+1)*nK/2;
        As{di} = zeros(nU^2, nKs);
        inds{di} = zeros(nKs, 1);
        k = 0;
        for i = 1:nK
            for j = i:nK
                if i == j
                    tmp = (gradu{di}.*K(:,i).*WdetJ(:))'*(gradu{di}.*K(:,j));
                else
                    tmp = (gradu{di}.*K(:,i).*WdetJ(:))'*(gradu{di}.*K(:,j)) + (gradu{di}.*K(:,j).*WdetJ(:))'*(gradu{di}.*K(:,i));
                end
                k = k + 1;
                As{di}(:,k) = tmp(:);
                inds{di}(k) = (j-1)*nK + i;
            end
        end
    end
elseif n == (nd+1)
    As = cell(2,1);
    inds = cell(2,1);
    % for the diagonal terms
    nK = size(Vs{1},2);
    K  = zeros(nq*ne, nK);
    for i = 1:nK
        k_atq = model.local_elem.f_quad_pts*reshape(Vs{1}(model.mesh.node_map',i), nb, ne);
        K(:,i) = reshape(k_atq,[],1);
    end
    % assemble squared K basis
    nKs = (nK+1)*nK/2;
    As{1} = zeros(nU^2, nKs);
    inds{1} = zeros(nKs, 1);
    k = 0;
    for i = 1:nK
        for j = i:nK
            tmp = zeros(nU);
            for di = 1:nd
                if i == j
                    tmp = tmp + (gradu{di}.*K(:,i).*WdetJ(:))'*(gradu{di}.*K(:,j));
                else
                    tmp = tmp + (gradu{di}.*K(:,i).*WdetJ(:))'*(gradu{di}.*K(:,j));
                    tmp = tmp + (gradu{di}.*K(:,j).*WdetJ(:))'*(gradu{di}.*K(:,i));
                end
            end
            k = k + 1;
            As{1}(:,k) = tmp(:);
            inds{1}(k) = (j-1)*nK + i;
        end
    end
    % for the cross terms, Vs{2} to Vs{nd+1} are columns of the basis matrices
    disp('the current parametrisation of the cross term is high-dimensional')
    nB = size(Vs{2},2);
    Bs = cell(nd, 1);
    for di = 1:nd
        Bs{di} = zeros(nq*ne, nB);
        for i = 1:nB
            b_atq = model.local_elem.f_quad_pts*reshape(Vs{1+di}(model.mesh.node_map',i), nb, ne);
            Bs{di}(:,i) = reshape(b_atq,[],1);
        end
    end
    %
    nBs = (nB+1)*nB/2;
    As{2} = zeros(nU^2, nBs);
    inds{2} = zeros(nBs, 1);
    k = 0;
    for i = 1:nB
        for j = i:nB
            tmp = zeros(nU);
            left  = zeros(size(gradu{di}));
            right = zeros(size(gradu{di}));
            for di = 1:nd
                left  = left  + gradu{di}.*Bs{di}(:,i);
                right = right + gradu{di}.*Bs{di}(:,j);
            end
            if i == j
                tmp = tmp + (left.*WdetJ(:))'*right;
            else
                tmp = tmp + (left.*WdetJ(:))'*right + (right.*WdetJ(:))'*left;
            end      
            k = k + 1;
            As{2}(:,k) = tmp(:);
            inds{2}(k) = (j-1)*nB + i;
        end
    end
else
    error('kappa not implemented')
end


end