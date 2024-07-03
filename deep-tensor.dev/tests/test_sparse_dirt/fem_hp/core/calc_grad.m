function grad_x = calc_grad(mesh, local_elem, u)

grad_b = cell(local_elem.num_dim, 1);
grad_x = cell(local_elem.num_dim, 1);

for di = 1:local_elem.num_dim
    grad_b{di} = local_elem.grad_f_quad_pts{di}*u(mesh.node_map');
end

for di = 1:local_elem.num_dim
    grad_x{di} = zeros(size(grad_b{di}));
    for dj = 1:local_elem.num_dim
        grad_x{di} = grad_x{di} + grad_b{dj}.*reshape(mesh.inv_Jt{di,dj},1,[]);
    end
end

end