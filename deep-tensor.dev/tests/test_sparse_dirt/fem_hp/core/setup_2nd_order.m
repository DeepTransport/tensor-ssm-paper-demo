function FEM = setup_2nd_order(p, t, opts)
%setup_FEM
%
% setup data structures used for FEM computation
%
% Tiangang Cui, 10/Sep/2020

% setup the reference element
[FEM.local_elem, FEM.local_elem_bnd] = local_elem_simplex(size(p,2), opts.poly_order, opts.quad_order);

% setup the mesh
FEM.h       = opts.h;
h_tol       = opts.h/opts.poly_order/20;
FEM.mesh    = calc_mesh_geometry(p, t, h_tol, FEM.local_elem);

% fill indicies 
[FEM.fill_i, FEM.fill_vec]  = calc_fill_index(FEM.mesh.node_map);
FEM.fill_d  = [FEM.mesh.dof, FEM.mesh.dof];

% mass matrix
f_ones      = ones(FEM.mesh.dof,1);
FEM.mass    = assemble_mass(FEM.mesh, FEM.local_elem, FEM.fill_i, f_ones, FEM.fill_d);
FEM.mass_L  = chol(FEM.mass, 'lower');
% setup boundaries
FEM.n_bnd       = length(opts.bnd_funcs);
FEM.bnd_mesh    = cell(FEM.n_bnd, 1);
FEM.bnd2domain  = cell(FEM.n_bnd, 1);
FEM.bnd_quadWJ  = cell(FEM.n_bnd, 1);
FEM.bnd_fill_i  = cell(FEM.n_bnd, 1);
FEM.mass_bnds   = cell(FEM.n_bnd, 1);

for i = 1:FEM.n_bnd
    % mesh
    FEM.bnd_mesh{i}   = calc_mesh_geometry_bnd(FEM.mesh, FEM.local_elem_bnd, opts.bnd_funcs{i});
    % map from boundary note to interior
    FEM.bnd2domain{i} = sparse(FEM.bnd_mesh{i}.node_indices, 1:FEM.bnd_mesh{i}.dof, ...
        ones(1,FEM.bnd_mesh{i}.dof), FEM.mesh.dof, FEM.bnd_mesh{i}.dof);
    % fill index
    [FEM.bnd_fill_i{i}, FEM.bnd_fill_vec{i}] = calc_fill_index(FEM.bnd_mesh{i}.node_map);
    % mass
    f_bnd = FEM.bnd2domain{i}*ones(FEM.bnd_mesh{i}.dof,1);
    FEM.mass_bnds{i}  = assemble_mass(FEM.bnd_mesh{i}, FEM.local_elem_bnd, FEM.bnd_fill_i{i}, f_bnd, FEM.fill_d);
end

% compute inv(J')*grad for all quadrature points
FEM.grad = cell(FEM.local_elem.num_dim, 1);
for di = 1:FEM.local_elem.num_dim
    FEM.grad{di} = zeros(FEM.local_elem.num_quad*FEM.local_elem.num_node, FEM.mesh.num_elem);
    for dj = 1:FEM.local_elem.num_dim
        FEM.grad{di} = FEM.grad{di} + reshape(FEM.local_elem.grad_f_quad_pts{dj},[],1)*reshape(FEM.mesh.inv_Jt{di,dj},1,[]);
    end
end

k_ones      = ones(FEM.mesh.dof,1);
FEM.stiff   = assemble_stiff(FEM.mesh, FEM.local_elem, FEM.grad, FEM.fill_i, k_ones, FEM.fill_d);

FEM.penalty = 1E2*opts.h^(-opts.poly_order-1);

end

% compute the mass and stiffness matrices on the boundary 
% local mass on the boundary
%fill = FEM.local_elem_bnd.mass(:) * FEM.mesh_botbnd.detJ';
%FEM.mass_botbnd = global_fill(FEM.quad_botbnd, FEM.mesh.dof, fill);
%
%fill = FEM.local_elem_bnd.mass(:) * FEM.mesh_topbnd.detJ';
%FEM.mass_topbnd = global_fill(FEM.quad_topbnd, FEM.mesh.dof, fill);

% local stiff on the boundary
%fill = FEM.local_elem_bnd.stiff_iso(:) * FEM.mesh_botbnd.detJ';
% this is on the glocal scale
%FEM.stiff_botbnd = global_fill(FEM.quad_botbnd, FEM.mesh.dof, fill);