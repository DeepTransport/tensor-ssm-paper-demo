function g = deri_adjoint_stiff_sol_sq(mesh, local_elem, kappa_type, kappa, lambda, state)

fill_i  = mesh.node_map';
grad_l  = calc_grad(mesh, local_elem, lambda);
grad_s  = calc_grad(mesh, local_elem, state);
%

switch kappa_type
    case {'scalar'}
        fill = zeros(local_elem.num_node*mesh.num_elem, 1);
        loc_ki = local_elem.f_quad_pts*kappa(fill_i);
        for di = 1:local_elem.num_dim
            tmp = ( local_elem.f_quad_pts.*local_elem.quad_weights(:) )'...
                *( 2*(grad_l{di}.*loc_ki.*grad_s{di}).*mesh.detJ(:)' );
            fill = fill + tmp(:);
        end
        g = accumarray(fill_i(:), fill(:), [mesh.dof,1], [], [], false);
    case {'vector'}
        g = zeros(mesh.dof, local_elem.num_dim);
        for di = 1:local_elem.num_dim
            loc_ki = local_elem.f_quad_pts*reshape(kappa(fill_i, di), local_elem.num_node, mesh.num_elem);
            tmp = ( local_elem.f_quad_pts.*local_elem.quad_weights(:) )'...
                *( 2*(grad_l{di}.*loc_ki.*grad_s{di}).*mesh.detJ(:)' );
            g(:,di) = accumarray(fill_i(:), tmp(:), [mesh.dof,1], [], [], false);
        end
    case {'tensor'}
        g = zeros(mesh.dof, local_elem.num_dim+1);
        % the diagonal
        fill = zeros(local_elem.num_node*mesh.num_elem, 1);
        loc_ki = local_elem.f_quad_pts*reshape(kappa(fill_i, 1), local_elem.num_node, mesh.num_elem);
        for di = 1:local_elem.num_dim
            tmp = ( local_elem.f_quad_pts.*local_elem.quad_weights(:) )'...
                *( 2*(grad_l{di}.*loc_ki.*grad_s{di}).*mesh.detJ(:)' );
            fill = fill + tmp(:);
        end
        g(:,1) = accumarray(fill_i(:), fill(:), [mesh.dof,1], [], [], false);
        % the cross
        grad_sb = zeros(size(grad_s{1}));
        grad_lb = zeros(size(grad_l{1}));
        for di = 1:local_elem.num_dim
            loc_bi = local_elem.f_quad_pts*reshape(kappa(fill_i, di+1), local_elem.num_node, mesh.num_elem);
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
