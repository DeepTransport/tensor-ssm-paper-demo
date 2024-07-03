function g = matvec_dstiff_sol(mesh, local_elem, kappa_type, kappa, d_kappa, state)

%
ne = mesh.num_elem;
nd = local_elem.num_dim;
nb = local_elem.num_node;
nq = local_elem.num_quad;
%
fill_i  = mesh.node_map';
grad_s  = calc_grad(mesh, local_elem, state);
fill    = zeros(nb, ne);
%
switch kappa_type
    case {'scalar'}
        grad_sdk = cell(nd,1);
        loc_ka = local_elem.f_quad_pts*d_kappa(fill_i);
        for di = 1:nd
            grad_sdk{di} = (loc_ka.*grad_s{di});
        end
        for di = 1:nd
            for dj = 1:nd
                geo = reshape(mesh.detJ,1,[]).*reshape(mesh.inv_Jt{di,dj},1,[]);
                tmp = (local_elem.grad_f_quad_pts{dj}.*local_elem.quad_weights(:))'*(grad_sdk{di}.*geo);
                fill = fill + tmp;
            end
        end
    case {'vector'}        
        grad_sdk = cell(nd,1);
        for di = 1:nd
            dk = d_kappa(:,di);
            loc_ka = local_elem.f_quad_pts*dk(fill_i);
            grad_sdk{di} = (loc_ka.*grad_s{di});
        end
        for di = 1:nd
            for dj = 1:nd
                geo = reshape(mesh.detJ,1,[]).*reshape(mesh.inv_Jt{di,dj},1,[]);
                tmp = (local_elem.grad_f_quad_pts{dj}.*local_elem.quad_weights(:))'*(grad_sdk{di}.*geo);
                fill = fill + tmp;
            end
        end
    case {'tensor'}
        grad_sdk = cell(nd,1);
        dk = d_kappa(:,1);
        loc_ka = local_elem.f_quad_pts*dk(fill_i);
        for di = 1:nd
            grad_sdk{di} = (loc_ka.*grad_s{di});
        end
        for di = 1:nd
            for dj = 1:nd
                geo = reshape(mesh.detJ,1,[]).*reshape(mesh.inv_Jt{di,dj},1,[]);
                tmp = (local_elem.grad_f_quad_pts{dj}.*local_elem.quad_weights(:))'*(grad_sdk{di}.*geo);
                fill = fill + tmp;
            end
        end
        % the cross
        %
        loc_db = cell(nd,1);
        loc_b  = cell(nd,1);
        for di = 1:nd
            db = d_kappa(:,di+1);
            b  = kappa(:,di+1);
            loc_db{di} = local_elem.f_quad_pts*db(fill_i);
            loc_b{di}  = local_elem.f_quad_pts*b(fill_i);
        end
        %
        grad_sdb = zeros(nq, ne);
        grad_sb  = zeros(nq, ne);
        for di = 1:nd
            grad_sdb = grad_sdb + grad_s{di}.*loc_db{di};
            grad_sb  = grad_sb  + grad_s{di}.*loc_b{di};
        end
        %
        loc_dbg = cell(nd,1);
        loc_bg  = cell(nd,1);
        for dj = 1:nd
            loc_dbg{dj} = zeros(nq, ne);
            loc_bg{dj}  = zeros(nq, ne);
            for di = 1:nd
                geo = reshape(mesh.detJ,1,[]).*reshape(mesh.inv_Jt{di,dj},1,[]);
                loc_dbg{dj} = loc_dbg{dj} + loc_db{di}.*geo;
                loc_bg{dj}  = loc_bg{dj}  + loc_b{di}.*geo;
            end
        end
        %
        for di = 1:nd
            %
            tmp = (local_elem.grad_f_quad_pts{di}.*local_elem.quad_weights(:))'*((grad_sdb.*loc_bg{di}));
            fill = fill + tmp;
            %
            tmp = (local_elem.grad_f_quad_pts{di}.*local_elem.quad_weights(:))'*((grad_sb.*loc_dbg{di}));
            fill = fill + tmp;
        end
end

g = accumarray(fill_i(:), fill(:), [mesh.dof,1], [], [], false);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

