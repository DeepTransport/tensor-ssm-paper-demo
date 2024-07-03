function mesh_bnd = calc_mesh_geometry_bnd(mesh, local_elem_bnd, func_bnd)

ind_p_bnd = func_bnd(mesh.grid_p);

% identify elements that are on the boundary 
mask = false(size(mesh.grid_t)); % false for not on the boundary
for i = 1:length(ind_p_bnd)
    mask_i = mesh.grid_t == ind_p_bnd(i);
    mask = or(mask, mask_i);
end

int_t_bnd = find(sum(mask, 2) >= mesh.num_dim); % avoiding boundary lines and points

if sum(find(sum(mask, 2) > mesh.num_dim))
    disp('Not implemented for any elements that have more than one boduanry face')
end

% boundary mesh
t = zeros(length(int_t_bnd), mesh.num_dim);
for i = 1:size(int_t_bnd)
    j = int_t_bnd(i);
    t(i,:) = mesh.grid_t(j, mask(j,:));
end

% nd_bnd = nd-1;
mesh_bnd.num_elem = size(t, 1);
mesh_bnd.num_dim  = mesh.num_dim-1;

mesh_bnd.node_indices = func_bnd(mesh.nodes);
mesh_bnd.dof = length(mesh_bnd.node_indices);

mesh_bnd.root = cell(mesh_bnd.num_elem);
mesh_bnd.axis = cell(mesh_bnd.num_elem);
mesh_bnd.detJ = zeros(mesh_bnd.num_elem, 1);
% set the transformed coordinate axis on the surface 
for i = 1:mesh_bnd.num_elem
    tmp = mesh.grid_p(t(i,:), :);
    mesh_bnd.root{i} = tmp(1,:);
    mesh_bnd.axis{i} = tmp(2:end,:) - repmat(mesh_bnd.root{i}, mesh_bnd.num_dim, 1);
    mesh_bnd.detJ(i) = sqrt( det(mesh_bnd.axis{i}*mesh_bnd.axis{i}') );
end

% this way, the transformed coordinates coincide with barycentric 
% mesh_bnd.root{i} + s * mesh_bnd.axis{i}, each row of is the plane coordinate
% s = [zeros(1,nd-1); eye(nd-1)] for the vertices in the new axis

% do we need the catesian to barycentric mapping?
% mesh_bnd.vertices_invT = eye(nd_bnd+1)/[zeros(1,nd_bnd); eye(nd_bnd)];

% map points with global interpolation points
mesh_bnd.node_map = zeros(mesh_bnd.num_elem, local_elem_bnd.num_node);
nodes_bnd = mesh.nodes(mesh_bnd.node_indices,:);
for i = 1:mesh_bnd.num_elem
    local_x_pts = local_elem_bnd.lambda_pts*mesh_bnd.axis{i} ...
        + repmat(mesh_bnd.root{i}, local_elem_bnd.num_node, 1);
    index1 = index_matching(nodes_bnd, local_x_pts, mesh.h_tol);
    %[t_index1, index2] = robust_intersect_old(nodes_bnd, local_x_pts, 'rows', 4);
    %if sum(abs(sort(index1) - sort(t_index1))) > 1E-8
    %    disp('Error in matching')
    %end
    if length(index1) ~= local_elem_bnd.num_node
        disp('>>>> Error: cannot match points')
        return
    end
    
    mesh_bnd.node_map(i,:) = mesh_bnd.node_indices(index1);
end

if sum( abs(unique(mesh_bnd.node_map(:)) - mesh_bnd.node_indices) ) > 1E-10
    disp('>>>> Error: wrong boundary nodes mapping')
end

mesh_bnd.grid_t = t;

mesh_bnd.quad_pts = cell(mesh_bnd.num_dim+1, 1);
for d = 1:(mesh.num_dim+1)
    mesh_bnd.quad_pts{d} = zeros(local_elem_bnd.num_quad, mesh_bnd.num_elem);
end

for i = 1:mesh_bnd.num_elem
    quad_t = local_elem_bnd.quad_lambdas*mesh_bnd.axis{i} ...
           + repmat(mesh_bnd.root{i}, local_elem_bnd.num_quad, 1);
    for d = 1:(mesh_bnd.num_dim+1)
        mesh_bnd.quad_pts{d}(:,i) = quad_t(:,d);
    end
end

end