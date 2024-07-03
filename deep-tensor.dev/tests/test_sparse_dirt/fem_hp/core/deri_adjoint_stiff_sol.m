function g = deri_adjoint_stiff_sol(mesh, local_elem, kappa_type, kappa, lambda, state)

fill_i  = mesh.node_map';
grad_l  = calc_grad(mesh, local_elem, lambda);
grad_s  = calc_grad(mesh, local_elem, state);
%
fill    = zeros(local_elem.num_node*mesh.num_elem, local_elem.num_dim);
%
for di = 1:local_elem.num_dim
    %fill_ka{di} = ( local_elem.f_quad_pts.*local_elem.quad_weights(:) )'...
    %             *( (grad_l{di}.*grad_s{di}).*mesh.detJ(:)' );
    tmp = ( local_elem.f_quad_pts.*local_elem.quad_weights(:) )'...
         *( (grad_l{di}.*grad_s{di}).*mesh.detJ(:)' );
    fill(:,di) = tmp(:);
end

switch kappa_type
    case {'scalar'}
        g = accumarray(fill_i(:), sum(fill,2), [mesh.dof,1], [], [], false);
    case {'vector'}
        g = zeros(mesh.dof, local_elem.num_dim);
        for di = 1:local_elem.num_dim
            g(:,di) = accumarray(fill_i(:), fill(:,di), [mesh.dof,1], [], [], false);
        end
    case {'tensor'}
        g = zeros(mesh.dof, local_elem.num_dim+1);
        % the diagonal
        g(:,1) = accumarray(fill_i(:), sum(fill,2), [mesh.dof,1], [], [], false);
        %{
        fill = zeros(size(fill_ka{1}));
        for di = 1:local_elem.num_dim
            fill = fill + fill_ka{di};
        end
        Jx(:,1) = accumarray(mesh.node_map(:), fill(:), [mesh.dof,1], [], [], false);
        %}
        % the cross
        grad_sb = zeros(size(grad_s{1}));
        grad_lb = zeros(size(grad_l{1}));
        for di = 1:local_elem.num_dim
            tmp = kappa(:,di+1);
            loc_bi = local_elem.f_quad_pts*tmp(fill_i);
            grad_sb = grad_sb + grad_s{di}.*loc_bi;
            grad_lb = grad_lb + grad_l{di}.*loc_bi;
        end
        for di = 1:local_elem.num_dim
            tmp = ( local_elem.f_quad_pts.*local_elem.quad_weights(:) )'...
                  *( (grad_l{di}.*grad_sb + grad_lb.*grad_s{di}).*mesh.detJ(:)' );
            g(:,1+di) = accumarray(fill_i(:), tmp(:), [mesh.dof,1], [], [], false);
        end
end

g = g(:);

end
