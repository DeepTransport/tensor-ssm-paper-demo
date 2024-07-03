function mesh = calc_mesh_geometry(p, t, h_tol, local_elem)

mesh.num_dim  = size(p, 2);
mesh.num_elem = size(t, 1);
mesh.h_tol    = h_tol;
% set the local transformation matrix and its inversion for coordinate
% transformation between baryentric and cartesian coordinates
mesh.vertices_invT = cell(mesh.num_elem, 1);
mesh.detJ   = zeros(mesh.num_elem, 1);
mesh.inv_J  = cell(mesh.num_dim);
mesh.inv_Jt = cell(mesh.num_dim);
mesh.inv_JtJ  = cell(mesh.num_dim);
tmp_inv_J   = zeros(mesh.num_elem, mesh.num_dim^2);
tmp_inv_Jt  = zeros(mesh.num_elem, mesh.num_dim^2);
tmp_inv_JtJ = zeros(mesh.num_elem, mesh.num_dim^2);
for i = 1:mesh.num_elem
    tmp0 = align_cartesian_pts(p(t(i,:),:));
    mesh.vertices_invT{i} = eye(mesh.num_dim+1)/tmp0;
    
    tmp1 = p(t(i,:),:);
    % dx_dlambda
    J = [tmp1(1,:)' - tmp1(3,:)', tmp1(2,:)' - tmp1(3,:)'];
    mesh.detJ(i) = det(J);
    % 
    inv_J   = inv(J);
    inv_Jt  = inv_J';
    tmp2    = inv_J*inv_Jt;
    
    tmp_inv_J(i,:)   = inv_J(:);
    tmp_inv_Jt(i,:)  = inv_Jt(:);
    tmp_inv_JtJ(i,:) = tmp2(:);
end

for di = 1:mesh.num_dim
    for dj = 1:mesh.num_dim
        k = (dj-1)*mesh.num_dim + di;
        mesh.inv_JtJ{di, dj} = tmp_inv_JtJ(:, k);
        mesh.inv_J  {di, dj} = tmp_inv_J (:, k);
        mesh.inv_Jt {di, dj} = tmp_inv_Jt(:, k);
    end
end
% add points
node_map   = zeros(mesh.num_elem, local_elem.num_node);
buffer_pts = zeros(local_elem.num_node*mesh.num_elem, mesh.num_dim);
for i = 1:mesh.num_elem 
    ind = (i-1)*local_elem.num_node + (1:local_elem.num_node);
    local_x_pts = barycentric2cartesian(p(t(i,:), :), local_elem.lambda_pts);
    buffer_pts(ind,:) = local_x_pts;
    node_map(i,:) = ind;
end
% delete repetitive points
[mesh.nodes, ~, ind_nodes] = robust_unique(buffer_pts, h_tol);

%[t_ind_nodes, ind_nodes]
%sum(abs(t_ind_nodes - ind_nodes))

% rebuild node_map
mesh.node_map = ind_nodes(node_map);

% triangulation for ploting
if mesh.num_dim == 2
    mesh.node_tri = delaunay(mesh.nodes(:,1), mesh.nodes(:,2));
end

mesh.dof = size(mesh.nodes, 1);
mesh.grid_p = p;
mesh.grid_t = t;

mesh.quad_pts = cell(mesh.num_dim, 1);
for d = 1:mesh.num_dim
    mesh.quad_pts{d} = zeros(local_elem.num_quad, mesh.num_elem);
end

quad_t = [local_elem.quad_lambdas, 1-sum(local_elem.quad_lambdas,2)];
for i = 1:mesh.num_elem
    tmp_pts = quad_t*mesh.grid_p(mesh.grid_t(i,:),:);
    for d = 1:mesh.num_dim
        mesh.quad_pts{d}(:,i) = tmp_pts(:,d);
    end
end

end